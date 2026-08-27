# Library QC metrics: concatemers, jackpotting, barcode hopping

Per-library QC for a run that has been through `dorado demux`. Every check is self-contained —
paste and run.

## Setup

```bash
export SCRIPTS=/path/to/scripts      # pipeline scripts
export DATA=/data                    # parent for per-experiment directories
export EXP=EXP26001055               # the run being checked
export R=$DATA/$EXP
```

Python 3.11 specifically — 3.13 resolves to an ancient `umi_tools` that fails at import:

```bash
conda create -n hbv_lr -c conda-forge -c bioconda \
  python=3.11 minimap2 'samtools>=1.16' pychopper \
  'umi_tools=1.1.6' pysam 'pandas<3' numpy scipy statsmodels matplotlib pyarrow setuptools
conda activate hbv_lr
pip install edlib                    # no conda build for py3.11
```

`setuptools` explicitly or `umi_tools` fails on `pkg_resources`; `pyarrow` for parquet; `edlib` for
the concatemer filter's fuzzy matching; `samtools ≥ 1.16` because phase 2 exits below that. If `pip`
hits PEP 668, the environment has no pip of its own — install pip with conda first. `dorado` is a
separate ONT binary. Verify once that `pychopper` really solved alongside `umi_tools=1.1.6` under
3.11.

Every session — `umi_tools` and `pyarrow` fail late otherwise:

```bash
conda activate hbv_lr
which samtools minimap2 pychopper umi_tools \
  && python3 -c "import pysam, numpy, pandas, pyarrow, scipy, edlib; print('ok')"
```

## Where each metric comes from

```
dorado demux --barcode-both-ends
  │  heterologous fusions (two libraries) → end barcodes disagree → unclassified
  │  homologous fusions  (one library)    → end barcodes agree    → pass through
  │
  ├── 1  assignment rate
  └── 2  barcode identity + sequence-level bleed          [GATE]
        ↓
phase1_align → phase1_5_readqc                            (read length N50, median, mean Q)
        ↓
phase2_hbv_umi
  ├─ extract HBV-aligned reads
  ├─ filter_concatemers.py ─── 3  concatemer rate         [GATE on the summary schema]
  ├─ pychopper                     (orient, extract UMI, trim primers)
  └─ umi_tools dedup
        ↓
  ├── 4  jackpotting
  ├── 5  duplication, merge ratio, coordinate complexity
  └── 6  molecule-level sharing
```

`--barcode-both-ends` removes only heterologous concatemers: a fusion across libraries carries
disagreeing end barcodes and lands in `unclassified`, while a fusion within one library agrees at
both ends and survives — which is why `filter_concatemers.py` exists.

**The concatemer filter must run before pychopper.** Pychopper trims the primers, so downstream of
it the adapter sequence is gone whether or not the read was a fusion. Detection is only possible
between HBV extraction and pychopper.

---

## 1. Demux assignment rate

```bash
find "$R/demux" -path '*bam_pass*' -name '*.bam' | while read -r f; do
  printf "%s\t%s\n" "$(basename "$(dirname "$f")")" "$(samtools view -c -F 0x900 "$f")"
done | sort | awk -F'\t' '
{ if (!($1 in n)) ord[++k] = $1
  n[$1] += $2; tot += $2
  if ($1 != "unclassified") asg += $2
  if ($1 !~ /^(barcode[0-9]+|unclassified)$/) badlabel = $1 }
END { for (i = 1; i <= k; i++) printf "%-16s %12d\n", ord[i], n[ord[i]]
      if (!tot) { print "no BAMs found"; exit 1 }
      printf "\nassigned: %d / %d = %.1f%%\n", asg, tot, 100*asg/tot
      if (badlabel) {
        printf "\nWARNING: label \"%s\" is not a barcode directory.\n", badlabel
        print  "         The BAM parent directory is not the barcode level in this layout,"
        print  "         so counts are being summed across barcodes and the rate above is wrong." }
      else if (k < 2) print "\nWARNING: only one label found - check the demux output layout" }'
```

Counts each BAM separately — `samtools view` takes one file plus *regions*, so a glob would parse
the second filename as a region. The script warns if the BAM's parent directory isn't the barcode
level, which would otherwise sum across barcodes and report 100%.

**Expect 60–80%**; near 50% is consistent with a very short-fragment run. To turn that into a test,
run demux once without `--barcode-both-ends` and check whether the both-ends rate ≈ the square of
the single-end rate — an arrangement-mask problem depresses both, a library problem doesn't.

Read the spread as well as the total: one order of magnitude across the pool is workable, two or
more means the small libraries are being measured under conditions the large ones are not, and it
drives both hopping measures.

| observation | reading |
|---|---|
| 60–80% | expected |
| low on **one** barcode | that library |
| low **across every barcode at once** | arrangement masks don't match the construct, or adapters didn't survive basecalling |

---

## 2. Barcode identity and sequence-level bleed — GATE

```bash
python3 - "$R" "$SCRIPTS/barcode_sequences_$EXP.fasta" <<'PY' || echo "GATE FAILED - stop here"
import pysam, glob, os, sys, random
root, fasta = sys.argv[1], sys.argv[2]
N = int(sys.argv[3]) if len(sys.argv) > 3 else 50000
IUPAC = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")

def load_barcodes(path):                       # multi-line safe, case safe, validated
    seqs, name, buf = {}, None, []
    def flush():
        if name and buf:
            s = "".join(buf).upper()
            if name in seqs: sys.exit(f"duplicate barcode name: {name}")
            bad = set(s) - set("ACGTRYSWKMBDHVN")
            if bad: sys.exit(f"{name}: unexpected characters {sorted(bad)}")
            seqs[name] = s
    for line in open(path):
        line = line.strip()
        if line.startswith(">"): flush(); name, buf = line[1:].split()[0], []
        elif line: buf.append(line)
    flush()
    if not seqs: sys.exit(f"no sequences parsed from {path}")
    L = {len(s) for s in seqs.values()}
    if len(L) != 1: sys.exit(f"barcodes differ in length: {sorted(L)}")
    return [(n, s, s.translate(IUPAC)[::-1]) for n, s in seqs.items()], L.pop()

def load_expected(root):                       # sample_key -> barcode, from the sample sheet
    m, hdr = {}, None
    p = os.path.join(root, "config/samples.tsv")
    if not os.path.exists(p): return m
    for line in open(p):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if hdr is None: hdr = {v: i for i, v in enumerate(f)}; continue
        if len(f) > max(hdr.get("barcode", 0), hdr.get("sample_key", 0)):
            m[f[hdr["sample_key"]]] = f[hdr["barcode"]]
    return m

bcs, BL = load_barcodes(fasta)
expected = load_expected(root)
print(f"{len(bcs)} barcodes x {BL} nt; sample sheet maps {len(expected)} libraries\n")
hdr = (f"{'library':<14}{'expected':>11}{'observed':>11}{'diag':>7}{'Q-pred':>8}{'ratio':>6}"
       f"{'offdiag':>9}{'sampled':>9}  flags")
print(hdr); print("-"*len(hdr))
gate = []
for d in sorted(glob.glob(os.path.join(root, 'analysis/samples/*/'))):
    bam = os.path.join(d, 'aligned_sorted.bam')
    if not os.path.exists(bam): continue
    name = os.path.basename(d.rstrip('/'))
    rng, res, total = random.Random(0), [], 0
    with pysam.AlignmentFile(bam) as af:
        for r in af.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary: continue   # one row per read
            s = r.query_sequence
            if not s: continue
            total += 1
            if len(res) < N: res.append((s, r.query_qualities))  # reservoir sample:
            else:                                               # representative, reproducible
                j = rng.randrange(total)
                if j < N: res[j] = (s, r.query_qualities)
    if not res: print(f"{name:<14} no reads"); continue
    counts, exp_sum, qn = [0]*len(bcs), 0.0, 0
    for s, q in res:
        for i, (_, fwd, rev) in enumerate(bcs):
            if fwd in s or rev in s: counts[i] += 1
        if q is not None and len(q) >= 2*BL:   # this read's own P(exact match at either end)
            p5 = p3 = 1.0
            for x in q[:BL]:  p5 *= 1 - 10**(-x/10)
            for x in q[-BL:]: p3 *= 1 - 10**(-x/10)
            exp_sum += 1 - (1-p5)*(1-p3); qn += 1
    n = len(res); qpred = 100*exp_sum/qn if qn else 0
    obs_i = max(range(len(bcs)), key=lambda i: counts[i])
    obs = bcs[obs_i][0]; exp_bc = expected.get(name)
    diag_i = next((i for i, b in enumerate(bcs) if b[0] == exp_bc), obs_i)
    diag = 100*counts[diag_i]/n
    off = [(100*c/n, bcs[i][0]) for i, c in enumerate(counts) if i != diag_i]
    offmax, offname = max(off) if off else (0.0, "-")
    flags = []
    if exp_bc and obs != exp_bc:
        flags.append(f"*** IDENTITY MISMATCH: reads carry {obs} ***"); gate.append(name)
    elif offmax > 0.5*diag and diag > 0:
        flags.append(f"*** off-diagonal {offname} rivals diagonal ***"); gate.append(name)
    if qpred and diag < 0.75*qpred: flags.append("diagonal below Q prediction")
    if n < 1000: flags.append(f"off-diag resolution only {100/n:.2f}%")
    print(f"{name:<14}{str(exp_bc):>11}{obs:>11}{diag:>6.1f}%{qpred:>7.1f}%"
          f"{diag/qpred if qpred else 0:>6.2f}{offmax:>8.2f}%{n:>9,}  {' '.join(flags)}")
print()
if gate: sys.exit(f"STOP: identity unconfirmed for {', '.join(gate)} — "
                  "fix the sample sheet before interpreting anything downstream")
print("identity confirmed for all libraries")
PY
```

Counts which barcode *sequence* appears in each library's reads, and checks it against
`config/samples.tsv`. Identity needs confirming because dorado names output directories
positionally (`barcode01`, `barcode02`) and the demux summary uses the same names, so nothing else
in the pipeline catches a transposed sample sheet. Exits non-zero if any library's reads carry a
barcode the sheet doesn't claim.

Runs against `aligned_sorted.bam` — `--no-trim` retains adapters and minimap2 soft-clips rather
than removing sequence. Before phase 1, point it at the demux BAMs. Reads are reservoir-sampled,
not taken as the first N: a coordinate-sorted BAM's first N records are the leftmost N alignments,
and coverage concentrates at hotspots.

**The diagonal will not be 100%** — it is an exact 24-mer search, so one error defeats it. The
script predicts the diagonal from each sampled read's own base qualities and prints the ratio;
compare that, not the raw figure. For orientation, uniform-quality predictions are 0.71 at Q15,
0.90 at Q18, 0.95 at Q20, 0.99 at Q23 — so 0.90 is clean on Q18 reads and means 10% of reads lack
a detectable barcode on Q23 reads.

**The off-diagonal is the bleed rate.** ONT has no ExAmp index hopping, but libraries are barcoded
individually by PCR and pooled *before* hybridisation capture, so any post-pooling amplification
can form cross-library chimeras. `--barcode-both-ends` already removed fusions that put a
disagreeing barcode at a read *end*, so what remains is internal. The rate scales with the share of
the pool held by the dominant library; a balanced pool is the control you have.

| observation | reading |
|---|---|
| ratio (diag / Q-prediction) ≈ 1 | as good as read quality allows |
| ratio < 0.75 | reads genuinely lack a detectable barcode — check trimming and adapter survival |
| off-diagonal < 0.1% | background (needs ≥ 1,000 sampled reads to be measurable; the script flags this) |
| off-diagonal 0.1–1% | bleed from a dominant library; rebalance next time |
| off-diagonal > 1% | investigate before trusting per-library numbers |
| **identity mismatch, or one off-diagonal entry rivalling the diagonal** | **transposed sample sheet. Stop — every downstream result is attached to the wrong sample.** |

---

## 3. Concatemer rate — GATE on the summary schema

Two or more cDNA molecules ligated into one read. Homologous fusions survive demultiplexing and do
two kinds of damage: a fused pair carries **two UMIs**, so dedup sees whichever pychopper extracted
and the second molecule is lost; and on the 2× reference a fused pair can align contiguously across
the copy junction and satisfy the ≥ 3,982 bp `pgRNA_RT` span gate, scoring a PCR artifact as
tandem/concatemeric readthrough. The filter prevents the second.

```bash
PROJECT_ROOT=$R "$SCRIPTS/phase2_hbv_umi_v2.sh"
awk '
BEGIN { FS="\t"
        split("sample_key hbv_primary_reads hbv_concatemers_removed hbv_concatemer_pct", req, " ") }
/^#/ { next }                                  # these summaries open with a comment line
!hdr {
  for (i = 1; i <= NF; i++) h[$i] = i
  hdr = 1
  for (j in req) if (!(req[j] in h)) missing = missing (missing ? ", " : "") req[j]
  if (missing) { printf "WARNING: column(s) not present: %s\n", missing
                 print  "         this summary predates the concatemer filter, or the schema changed."
                 print  "         Molecule counts downstream INCLUDE fused reads. Reprocess, or"
                 print  "         measure retrospectively, before comparing with a filtered run."
                 exit 3 }
  printf "%-14s %10s %10s %8s   %s\n", "library", "reads", "concat", "pct", "reading"
  next
}
NF < 2 { next }
{ pct = $h["hbv_concatemer_pct"] + 0
  reading = (pct < 1)  ? "background" \
          : (pct <= 5) ? "present - confirm the filter ran before dedup" \
                       : "material loss - check ligation, blockers, pool composition"
  printf "%-14s %10s %10s %7.2f%%   %s\n", $h["sample_key"], $h["hbv_primary_reads"],
         $h["hbv_concatemers_removed"], pct, reading }' \
    "$R/analysis/comparison/phase2_hbv_umi_summary.tsv" || echo "schema gate failed"
```

**An absent `hbv_concatemer_pct` is not cosmetic** — the library was processed by an older phase 2
and **its molecule counts include fused reads**, so they aren't comparable with a filtered run. The
script exits 3 rather than printing a partial table. Reprocess, or measure retrospectively:

```bash
cd "$SCRIPTS"
HBV_CONTIG=U95551.1_2x        # check yours: samtools idxstats <a bam> | head

for d in "$R"/analysis/samples/*/; do
  bam="$d/aligned_sorted.bam"
  [ -f "$bam" ] || continue
  [ -f "$bam.bai" ] || samtools index "$bam"        # region fetch needs an index
  samtools view -b -F 0x900 "$bam" "$HBV_CONTIG" \
    | samtools fastq -n - 2>/dev/null \
    | python3 filter_concatemers.py --barcodes-fasta "barcode_sequences_$EXP.fasta" \
        --label "$(basename "$d")" > /dev/null
done
```

Counts print to stderr. Use `samtools fastq -n`, not an `awk` on the SAM — reverse-strand reads need
reverse-complementing back to original orientation before the motif scan. A wrong `HBV_CONTIG`
returns nothing and the rate reads 0%.

### Reading it

Detection needs an adapter or barcode motif more than 200 bp from **both** ends, so a read must
exceed ~427 bp to be scannable at all. Report the rate conditioned on scannable reads, with the
fraction eligible, or the denominator is doing the work. Note the limit is set by the length of the
*fused* read, roughly twice the insert — a 250 bp-insert library still makes ~500 bp fusions that
are detectable, so "short library, therefore low rate" needs the eligible fraction to back it up.
**The reported figure is a floor**: junctions within 200 bp of either end and degraded internal
adapters are both invisible.

**Expect under 2% with blocking oligos.** Set by ligation chemistry, not by the sample — without
blockers, or with a large molar excess of adapter over insert (which is what low input produces),
it rises. **Measure it per run and never carry a rate over:** the same physical library measured
0.97% and 5.43% in two runs, 5.6× apart on identical code, the second in a pool dominated by
libraries ~1,600× larger.

| rate | reading |
|---|---|
| < 1% | background |
| 1–5% | present; confirm the filter ran before dedup |
| > 5% | material molecule loss; check ligation conditions and blocking oligos |

---

## 4. Jackpotting

PCR amplifies molecules unevenly. A molecule that amplifies early is amplified again every
subsequent cycle, so final abundance is exponential in how early that happened; with few starting
molecules and many cycles one molecule can dominate. The consequence is that **read depth stops
corresponding to molecular abundance**.

```bash
python3 - "$R" 50 <<'PY'      # 50 = position clustering tolerance, bp
import pysam, numpy as np, glob, os, sys, math
from collections import defaultdict
root = sys.argv[1]; TOL = int(sys.argv[2]) if len(sys.argv) > 2 else 50

def ztp_lambda(kbar):                      # solve kbar = lam/(1-exp(-lam))
    if kbar <= 1.0000001: return 1e-6
    lo, hi = 1e-9, 1e4
    for _ in range(200):
        mid = (lo+hi)/2
        if mid/(1-math.exp(-mid)) < kbar: lo = mid
        else: hi = mid
    return (lo+hi)/2

def ztp_top1(lam):                         # expected top-1% read share under ZTP(lam)
    if lam < 1e-6: return 1.0
    kmax = int(lam + 40*math.sqrt(lam) + 60)
    den = 1-math.exp(-lam); ll = math.log(lam)
    p = {k: math.exp(-lam + k*ll - math.lgamma(k+1))/den for k in range(1, kmax+1)}
    mean = sum(k*v for k, v in p.items())
    rem, reads = 0.01, 0.0
    for k in sorted(p, reverse=True):
        take = min(p[k], rem); reads += take*k; rem -= take
        if rem <= 1e-15: break
    return 100*reads/mean

hdr = (f"{'library':<14}{'mol':>9}{'reads':>10}{'med':>5}{'mean':>7}{'max':>8}"
       f"{'top1%':>7}{'null':>6}{'ratio':>6}{'jack>=1k':>9}  flags")
print(hdr); print("-"*len(hdr))
for d in sorted(glob.glob(os.path.join(root, 'analysis/samples/*/'))):
    b = os.path.join(d, 'hbv.umi.bam')
    if not os.path.exists(b): continue
    name = os.path.basename(d.rstrip('/'))
    pos_by_umi, malformed, total = defaultdict(list), 0, 0
    with pysam.AlignmentFile(b) as af:
        for r in af.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
            total += 1
            if r.has_tag('RX') and len(r.get_tag('RX')) == 28:
                pos_by_umi[r.get_tag('RX')].append(r.reference_start)
            else: malformed += 1
    if not pos_by_umi:
        print(f"{name:<14} no reads with a well-formed 28 nt UMI"); continue
    counts = []                            # single-linkage cluster positions within TOL
    for ps in pos_by_umi.values():
        ps.sort(); n, prev = 1, ps[0]
        for p in ps[1:]:
            if p - prev > TOL: counts.append(n); n = 1
            else: n += 1
            prev = p
        counts.append(n)
    v = np.array(sorted(counts, reverse=True)); t = int(v.sum())
    kbar = t/len(v); lam = ztp_lambda(kbar); null = ztp_top1(lam)
    top1 = 100*v[:max(1, len(v)//100)].sum()/t
    flags = []
    if len(v) < 300: flags.append("N<300:top1%-unreliable")
    if v.max() == 1: flags.append("INPUT-LOOKS-DEDUPED")
    if malformed/max(1, total) > 0.5: flags.append(f"{100*malformed/total:.0f}%-UMIs-dropped")
    print(f"{name:<14}{len(v):>9,}{t:>10,}{np.median(v):>5.0f}{kbar:>7.2f}{v.max():>8,}"
          f"{top1:>6.1f}%{null:>5.1f}%{top1/null:>6.1f}{int((v>=1000).sum()):>9}  {' '.join(flags)}")
PY
```

A molecule is a UMI plus a position, and positions are **clustered within a tolerance**, not binned.
A hard `reference_start // 50` splits a molecule whose reads straddle a boundary — at 6 bp of
alignment-start jitter that inflates molecule counts by up to ~1.5× and deflates reads per molecule
by up to a third, biasing everything here toward "less jackpotted than it is". Use the same
tolerance as section 6.

Three flags cover the ways this reads clean while being wrong: `INPUT-LOOKS-DEDUPED` means
`hbv.umi.bam` is post-dedup, so every count is 1 and jackpotting is invisible; `N<300` means the
top-1% slice is taking more than 1% of molecules, so the bands below don't apply; `%-UMIs-dropped`
means the 28 nt length filter is discarding most reads.

### Reading it

**Expected top-1% share is 1.2–3.0%** under Poisson sampling and nearly flat in λ, which is what
makes the metric usable without knowing depth in advance. The script solves for your run's λ,
prints that null and the ratio. It is a floor — well-behaved libraries sit somewhat above ratio 1.

| ratio (observed top 1% / null) | reading |
|---|---|
| ≈ 1 | duplication is Poisson — near-uniform |
| 1–10 | normal for an enriched library |
| 10–30 | jackpotting present; use molecules, never reads |
| > 30 | dominated by a few molecules; counts unreliable |

**Median versus mean is the second diagnostic.** In a jackpotted library the median stays at 1–2
while the mean climbs into the tens or hundreds; a mean or average-coverage figure alone hides this.
A median of 1 alongside a mean in the tens cannot arise from any single Poisson.

**Jackpot count can fall while severity rises** — a library with very few starting molecules may
show fewer molecules above 1,000 reads because amplification concentrated on fewer templates and
made each jackpot larger. Read the count next to `max` and the ratio.

**Watch for wasted sequencing.** A library that received the most reads and returned the fewest
molecules has been sequenced past the point of return — observed as bad as 16× worse molecular yield
per read at the lowest input of a titration.

---

## 5. Duplication, merge ratio, coordinate complexity

Two questions: has sequencing exhausted the library, and is the molecule count even well defined?

The second is easy to miss. `umi_tools dedup` in `directional` mode collapses UMIs within one edit
distance when read counts are consistent with one being a sequencing error of the other;
`unique`/`exact` does not. **`merge_ratio` = `unique_exact` / `unique_directional`** measures how
much collapsing happens. At high duplication, amplification generates enough reads per molecule
that UMI errors produce many near-neighbour UMIs: `exact` counts each as new and inflates,
`directional` collapses them and may over-collapse genuinely distinct ones. When they disagree
substantially the molecule count is **method-dependent, not merely uncertain**.

```bash
cd "$SCRIPTS"
for sk in $(awk -F'\t' '/^#/ || $1=="barcode" {next} NF>3 {print $4}' \
              "$R/config/samples.tsv"); do
  PROJECT_ROOT="$R" ./umi_saturation_curve.sh "$sk" 2>&1 | tee "/tmp/sat_$sk.log"
done
```

`$1=="barcode"` skips the header; `NR>1` would fail because these sheets open with comment lines,
and the header row itself would be passed to the script as a sample key.

> The script **overwrites** `analysis/comparison/umi_saturation.tsv` on each call rather than
> appending, so the file holds only the last library — hence the per-library logs.

Per library, at five subsampled depths: `dup_rate`, `molecules_per_1k_reads`, `merge_ratio`.

### Reading it

The bands below are **operating limits, not derived values** — no closed form fits, because error
variants are themselves amplified and `directional` collapses a fraction of them. Establish the
crossing point for your own chemistry from a titration; the duplication level at which the ratio
passes 1.3 is the number worth recording (observed once just above ~96% duplication, with a 2.2×
disagreement at 99.7%). The ratio also rises with depth *within* a library, so compare at matched
depth.

| merge ratio | reading |
|---|---|
| < 1.15 | negligible collapsing; counts robust to method |
| 1.15–1.30 | mild; note it |
| 1.30–1.60 | method-dependent; report which method was used |
| > 1.60 | molecule count is not well defined; report it as a range |

**Falling `mol/1k` means saturating; flat means undersampled** — a flat curve says the molecule
count is a floor. Whole runs have come back with nothing saturated, in which case every molecule
count in them is a lower bound.

### Coordinate complexity

At low input, UMI counts flatten while alignment coordinates keep resolving — UMI counts have
differed by only 1.2× across a 10× input change at the bottom of a titration. Distinct coordinates
cross-check UMI over-splitting:

```bash
python3 - "$R" 25 <<'PY'
import pandas as pd, glob, os, sys
root = sys.argv[1]; BIN = int(sys.argv[2]) if len(sys.argv) > 2 else 25
paths = sorted(glob.glob(os.path.join(root, 'analysis/samples/*/hbv_classified_molecules.parquet')))
if not paths: sys.exit("no hbv_classified_molecules.parquet found — phase3b has not run yet")
print(f"coordinate bin = {BIN} bp (match this to the jackpotting and sharing tolerances)\n")
print(f"{'library':<14}{'molecules':>11}{'coords_exact':>13}{'coords_binned':>14}"
      f"{'mol/coord':>10}{'top5':>7}  reading")
print("-"*90)
for p in paths:
    d = pd.read_parquet(p); name = os.path.basename(os.path.dirname(p))
    if len(d) == 0: print(f"{name:<14} empty"); continue
    ex = len(d.groupby(['ref_start','ref_end']))
    g  = d.assign(_s=d.ref_start//BIN, _e=d.ref_end//BIN).groupby(['_s','_e']).size()
    r  = len(d)/len(g)
    reading = ("counts agree" if r < 1.15 else
               "mild UMI over-splitting" if r < 1.6 else
               "molecule count exceeds distinct positions substantially - UMI count inflated")
    print(f"{name:<14}{len(d):>11,}{ex:>13,}{len(g):>14,}{r:>10.2f}"
          f"{100*g.nlargest(5).sum()/len(d):>6.1f}%  {reading}")
    if ex/len(g) > 1.5:
        print(f"{'':<14}exact grouping gives {ex:,} vs {len(g):,} binned - exact coordinates are "
              f"resolving alignment jitter, not distinct molecules")
PY
```

Needs `phase3b_classify_molecules.py` to have run. Two bounds have to be stated:

- **Coordinates must be binned, not matched exactly.** With exact `(ref_start, ref_end)` grouping,
  1 bp of jitter drives molecules/coords from a true 1.40 to 1.02 — the metric becomes a near-copy
  of the UMI count and detects nothing. It is then also *not* amplification-independent: with true
  complexity held fixed, exact `distinct_coords` moved 23× across a 400-fold amplification range,
  because over-splitting rises with reads per molecule. Binning at 25 bp cuts that to 1.6×.
- **The binned count is a floor on complexity and `mol/coord` is a ceiling on over-splitting**,
  since distinct fragments can share a bin — in a 6.4 kb reference a clean library reads ~1.4
  rather than 1.0. Monotonic in true over-splitting, so use it to compare libraries at matched
  depth, not as an absolute.

---

## 6. Molecule-level sharing

A molecule is a (UMI, position) pair. If the same one appears in two libraries, either the UMI
recurred by chance *and* landed at the same coordinate, or material physically moved between
libraries. **This is the measure that matters for quantification.** It is a different phenomenon
from section 2, not the same one at finer resolution: section 2 finds a foreign barcode inside a
read, this finds a molecule counted in the wrong library. A clean off-diagonal is not reassurance
about this.

```bash
cd "$SCRIPTS"
python3 cross_library_umi_check_v2.py --root "$R" --tol 50 \
  --out-tsv "$R/analysis/comparison/umi_sharing_v2_tol50.tsv"

python3 cross_library_umi_check_v2.py --root "$R" --tol 10 \
  --out-tsv "$R/analysis/comparison/umi_sharing_v2_tol10.tsv"
```

Pass `--tol` explicitly in both so the log records it. Reads `config/samples.tsv` and each library's
`hbv.umi.bam`, so phase 2 must have completed; runtime scales with the product of the two libraries'
molecule counts. Columns: `shared_umi_and_position` (raw count), `permuted_expectation`,
`excess_molecules`, `pct_of_smaller_library`.

Chance sharing is not negligible and grows with library size, so **a raw shared-UMI count is not
interpretable on its own** — it will mislead in the direction of seeing contamination that isn't
there. `permuted_expectation` shuffles UMIs within each library while every molecule keeps its
position and read count. Each library's UMI set is unchanged, so the shared-UMI set is unchanged
and the test measures only whether shared UMIs sit at the same coordinate. That is what makes it
robust to non-uniform UMI synthesis: synthesis bias can put a UMI in two libraries, it cannot put
it at the same coordinate.

> Residual bias: the shuffle breaks the UMI↔read-count link, so the null is mildly optimistic if
> concordance correlates with read count. Stratifying the shuffle within read-count bins closes it.

### Reading it

**Position concordance is the discriminating statistic**, not the shared-UMI count — libraries that
shared a capture put shared molecules at the same coordinate (observed 88–100%) while libraries from
different runs sit at 5–16%.

**Statistical significance is not importance.** A z-score above 20 on an excess of 36 molecules is
real, tightly localised and negligible. Read `pct_of_smaller_library`, not the z.

**Check the tolerance.** Identical counts at ±10 and ±50 bp means the tolerance contributed nothing
and the matches are genuinely tight; if they collapse as the window narrows, the tolerance was
generating the signal.

**Magnitude scales with pool imbalance, not with technique.** A pool spanning 11.5× in molecule
count gave worst-pair sharing of 0.75%; a pool spanning 244× gave 76.6%.

| excess as % of smaller library | reading |
|---|---|
| < 1% | background |
| 1–10% | present; note it when comparing libraries of different size |
| > 10% | the smaller library's counts are substantially foreign material |

**Use a negative control when the answer matters** — libraries from different runs were never in the
same capture, so any excess between them is method artifact:

```bash
OTHER=$DATA/EXP_FROM_A_DIFFERENT_RUN     # a real experiment directory
SHARED=SeqLibXXXX_something              # any library present in BOTH runs; omit if none

python3 cross_library_umi_check_v2.py --root "$R" --root-b "$OTHER" --tol 50 \
    --exclude-b "$SHARED" \
    --out-tsv "$R/analysis/comparison/umi_negative_control.tsv"
```

`--exclude-b` must name any library physically present in both runs, since it legitimately shares
molecules with its own prior run. Repeatable. **One weakness:** a negative control only tests the
size regime it covers, so if its largest partner is far smaller than the within-run partner it never
probes the regime the finding lives in. Normalise to excess per partner molecule to close that gap.

---

## Expectations

| metric | expectation | strength |
|---|---|---|
| assignment rate | 60–80%; ≈ *p*², *p* = single-end detection | derived form; *p* measurable by one demux without `--barcode-both-ends` |
| barcode diagonal | the run's own Q-based prediction | derived, computed per run |
| sequence-level bleed | < 1%, scales with pool imbalance | mechanism only |
| concatemers | < 2% with blockers; a floor; needs reads > ~427 bp | chemistry; per run only, never carried over |
| top 1% share | 1.2–3.0%, flat in λ | derived, computed per run |
| median reads/molecule | ≈ λ | derived |
| merge ratio | rises with reads/molecule and with depth | direction only; crossing point from a titration |
| coordinate complexity | relative, at matched depth | bounded, not absolute |
| molecule sharing | must be computed | permutation null required |

Record what each run gives. Once several runs share a chemistry and a depth, observed ranges become
the better baseline for the metrics where theory gives only a direction.

## When a check fails

| failure | what to change |
|---|---|
| identity mismatch, or off-diagonal rivalling the diagonal | **transposed sample sheet — fix before interpreting anything** |
| assignment low across all barcodes at once | arrangement masks; verify adapters survived basecalling |
| assignment low on one barcode | that library — check input and adapter ligation |
| diagonal well below the Q prediction | reads lack a detectable barcode; check trimming |
| off-diagonal > 1% | rebalance the pool; don't co-pool dominant and minor libraries |
| `hbv_concatemer_pct` absent | filter did not run; reprocess or measure retrospectively |
| concatemers > 5% | check ligation conditions and blocking oligos |
| top-1% ratio > 30 | more input, or fewer PCR cycles; report molecules, never reads |
| `INPUT-LOOKS-DEDUPED` | pointing at a post-dedup BAM — jackpotting can't be measured from it |
| merge ratio > 1.6 | report which dedup method was used; treat the count as a range |
| `mol/1k` flat across depths | undersampled; counts are floors, sequence deeper |
| `mol/coord` > 1.6 at matched depth | UMI over-splitting; prefer the binned coordinate count |
| molecule sharing > 10% of the smaller library | rebalance; do not co-pool high and low viral load |

Two root causes recur: **too few starting molecules relative to PCR cycles** (jackpotting, merge
ratio, saturation) and **pool imbalance** (both hopping measures, plus capture competition). Neither
is fixable downstream.

## Record per run

- assignment rate and per-barcode read counts; the single-end rate if *p* was measured
- barcode diagonal, its Q-based prediction and the ratio; largest off-diagonal entry
- concatemer rate, the fraction of reads long enough to be scannable, and the N50 and median
- reads per molecule: median, mean, max, top-1% share, the null, and the ratio
- `merge_ratio` and duplication at full depth, plus the saturation curve shape
- molecules, binned distinct coordinates and `mol/coord`, with the bin width used
- molecule-level sharing at both tolerances, with the permutation null
- PCR cycle count and input ng per library
- position tolerance used in sections 4, 5 and 6 — they should match
- confirmation that the demux `[scoring]` block was unchanged

The 200 bp barcode windows in the demux `[scoring]` block are shared with concatemer detection, so
changing them makes concatemer and duplication rates non-comparable between runs. Record any change
loudly.
