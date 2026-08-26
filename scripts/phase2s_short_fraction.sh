#!/usr/bin/env bash
# Phase 2s — HBV in the SHORT fraction (the reads Phases 1-4 silently discard).
#
#   PROJECT_ROOT=/data/EXP26000993 bash phase2s_short_fraction.sh
#   PROJECT_ROOT=/data/EXP26000993 bash phase2s_short_fraction.sh SeqLib5576_yecuris1_1ng_20
#
# WHY THIS EXISTS
# phase1_align.sh aligns with `-x splice`, which cannot place a ~66 bp insert. In the
# serum/plasma libraries that discards 56-88% of all reads. Phase 2 then draws only from
# HBV-aligning PRIMARY reads, so every molecule count in Phases 2-4 is the LONG fraction
# and nothing else.
#
# That matters for two open questions:
#   1. TREATMENT COMPARISON. Total HBV per library = long + short. If the treated and
#      control animals differ in fragment length (they do: N50 550 vs 450), comparing
#      only the long fraction is biased by an unknown amount and in an unknown direction.
#   2. RNase H HYPOTHESIS. HBV polymerase degrades pgRNA during minus-strand synthesis.
#      If the short HBV fragments are RNase H products they should cluster in a pattern
#      tracking RT progression from DR1, not tile the genome uniformly. Positional
#      coverage of short vs long fragments tests that.
#
# WHAT IT DOES, per library
#   unmapped reads (-f 4 from aligned_sorted.bam)
#     -> pychopper -U    orient + extract the 28-nt UMI (still present: the UMI lives in
#                        the SSPII primer, which short fragments retain)
#     -> minimap2 -ax sr  short-read preset against the HBV-only contig. NOT splice —
#                        that is the whole point.
#     -> umi_tools dedup  molecules, directly comparable to the long-fraction counts
#
# Outputs, per library under analysis/samples/<key>/:
#   hbv_short.umi.bam(.bai)         RX-tagged short HBV reads (RETAINED)
#   hbv_short_pychopper_stats.tsv
#   hbv_short_coverage.txt          per-base depth, for the positional comparison
# Shared: analysis/comparison/phase2s_short_fraction_summary.tsv
#         analysis/comparison/phase2s_coverage_long_vs_short.tsv
#
# CAVEATS
#   - pychopper discards fragments whose primers it cannot resolve, which is plausibly
#     biased against the shortest and most degraded reads. Counts are a floor.
#   - `-x sr` on ~66 bp inserts against a 6.4 kb target will produce some spurious
#     alignments. MIN_MAPQ filters the worst of it; treat low-complexity hits with care.
#   - The short fraction was never capture-enriched, so its duplication should be much
#     lower than the long fraction's. A modest read count can therefore represent a
#     large molecule count — that is the point of measuring it.
set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT:-/data/EXP26000993}
ANALYSIS=$PROJECT_ROOT/analysis
SAMPLES=$PROJECT_ROOT/config/samples.tsv
REF_2X=$ANALYSIS/refs/hg38_hbv_2x.fa
HBV_2X=U95551.1_2x
HBV_LEN=3182                 # single-copy genome length
HBV=U95551.1                 # contig name in the 1x reference we build below
HBV_FA=${HBV_FA:-/tmp/hbv_1x_$$.fa}

# ---------------------------------------------------------------------------
# 1x REFERENCE, deliberately — not the 2x contig the rest of the pipeline uses.
#
# The 2x contig exists so genome-WRAPPING transcripts align as one continuous record.
# A ~66 bp fragment cannot wrap anything, so doubling buys nothing here and costs two
# things: every short read matches both copies equally and receives MAPQ 0 by
# construction (destroying MAPQ as a signal), and coverage splits across two copies,
# which makes the positional analysis harder to interpret.
#
# On a single copy short reads map uniquely, MAPQ is informative again, and coordinates
# are directly genomic with no `mod 3182` needed.
#
# COST: fragments spanning the circular junction (~66/3182 = 2%) cannot align to a
# linear single copy and are lost. Accepted.
# ---------------------------------------------------------------------------
SUMMARY=$ANALYSIS/comparison/phase2s_short_fraction_summary.tsv
COVOUT=$ANALYSIS/comparison/phase2s_coverage_long_vs_short.tsv
VER=phase2s-v1

Q=${Q:-0.4}
# On the 1x reference MAPQ is informative again — this is one of the reasons for using it.
# (On the 2x contig every non-junction-spanning read is MAPQ 0 by construction, which is
# why only 4-19% of long HBV reads clear MAPQ 20.) Default is still 0 so nothing is
# silently dropped; short reads in low-complexity or repeat regions (DR1/DR2) can
# legitimately be MAPQ 0. Raise it to filter spurious sr alignments once you have looked
# at the MAPQ distribution.
MIN_MAPQ=${MIN_MAPQ:-0}
MAX_INSERT=${MAX_INSERT:-400}    # only consider reads shorter than this as "short"
UMI_METHOD=${UMI_METHOD:-directional}

NPROC=$(nproc)
PY_THREADS=${PY_THREADS:-$(( NPROC > 8 ? NPROC - 6 : NPROC ))}
MM_THREADS=${MM_THREADS:-6}
SORT_THREADS=${SORT_THREADS:-4}

for t in pychopper minimap2 samtools umi_tools; do
    command -v "$t" >/dev/null || { echo "ERROR: $t not on PATH" >&2; exit 1; }
done
[ -s "$REF_2X" ] || { echo "ERROR: reference not found: $REF_2X" >&2; exit 1; }

# Build the 1x HBV reference: first copy only, renamed so downstream coordinates are
# unambiguously single-copy genomic. Built fresh — the existing .mmi files were made with
# splice preset -k15 and must not be reused for -x sr.
samtools faidx "$REF_2X" "${HBV_2X}:1-${HBV_LEN}" \
  | sed "1s/.*/>${HBV}/" > "$HBV_FA"
samtools faidx "$HBV_FA"
ACTUAL_LEN=$(awk -v c="$HBV" '$1==c{print $2}' "$HBV_FA.fai")
[ "$ACTUAL_LEN" = "$HBV_LEN" ] || { echo "ERROR: 1x reference is ${ACTUAL_LEN}bp, expected ${HBV_LEN}" >&2; exit 1; }
trap 'rm -f "$HBV_FA" "$HBV_FA".*' EXIT

echo "=== Phase 2s — short-fraction HBV ==="
echo "  project:   $PROJECT_ROOT"
echo "  reference: 1x $HBV (${ACTUAL_LEN} bp) — single copy, NOT the 2x contig"
echo "  preset:    minimap2 -ax sr   (NOT splice)"
echo "  max insert: $MAX_INSERT bp   min MAPQ: $MIN_MAPQ"
echo "  threads:   pychopper=$PY_THREADS minimap2=$MM_THREADS sort=$SORT_THREADS"
echo

WANT="$*"
want(){ [ -z "$WANT" ] && return 0; for w in $WANT; do [ "$w" = "$1" ]||[ "$w" = "$2" ] && return 0; done; return 1; }

if [ -z "$WANT" ] || [ ! -f "$SUMMARY" ]; then
    printf '# experiment=%s; pipeline_version=%s; generated=%s\n' \
        "$(basename "$PROJECT_ROOT")" "$VER" "$(date +%F)" > "$SUMMARY"
    printf 'sample_key\tbarcode\tunmapped_reads\tshort_reads_in\tpychopper_out\thbv_short_aligned\thbv_short_umi28\thbv_short_molecules\tshort_dup_rate\n' >> "$SUMMARY"
fi

hms(){ printf '%dm%02ds' $(( $1 / 60 )) $(( $1 % 60 )); }

awk -F'\t' 'NF>=4 && $1!~/^#/ && $1!="barcode"' "$SAMPLES" | while IFS=$'\t' read -r barcode lib_id sample_name sample_key rest; do
    want "$sample_key" "$barcode" || continue
    D=$ANALYSIS/samples/$sample_key
    BAM=$D/aligned_sorted.bam
    [ -s "$BAM" ] || { echo "WARN $sample_key: no aligned_sorted.bam"; continue; }

    T0=$SECONDS
    echo "[$(date +%T)] === $sample_key ($barcode) ==="

    NUNMAP=$(samtools view -c -f 4 "$BAM")
    echo "[$(date +%T)]   unmapped reads: $NUNMAP"
    if [ "$NUNMAP" -eq 0 ]; then
        printf '%s\t%s\t0\t0\t0\t0\t0\t0\tNA\n' "$sample_key" "$barcode" >> "$SUMMARY"
        continue
    fi

    # ---- extract short unmapped reads -> pychopper (UMI) -> sr align -> sort ----
    echo "[$(date +%T)]   stage 1/3 pychopper + sr alignment"
    S1=$SECONDS
    set -o pipefail
    samtools view -f 4 "$BAM" \
      | awk -v MAX="$MAX_INSERT" 'length($10) <= MAX {print "@"$1"\n"$10"\n+\n"$11}' \
      > "$D/.short_in.fastq"
    NSHORT=$(( $(wc -l < "$D/.short_in.fastq") / 4 ))
    echo "[$(date +%T)]   short reads (<=${MAX_INSERT}bp): $NSHORT"

    if [ "$NSHORT" -eq 0 ]; then
        rm -f "$D/.short_in.fastq"
        printf '%s\t%s\t%s\t0\t0\t0\t0\t0\tNA\n' "$sample_key" "$barcode" "$NUNMAP" >> "$SUMMARY"
        continue
    fi

    pychopper -k PCB114 -m edlib -U -y -q "$Q" -t "$PY_THREADS" \
        -S "$D/hbv_short_pychopper_stats.tsv" \
        "$D/.short_in.fastq" "$D/.short_trimmed.fastq" 2>"$D/hbv_short.log" || true
    NPY=$(( $(wc -l < "$D/.short_trimmed.fastq" 2>/dev/null || echo 0) / 4 ))
    echo "[$(date +%T)]   pychopper output: $NPY"

    minimap2 -ax sr --secondary=no -y -t "$MM_THREADS" "$HBV_FA" \
        "$D/.short_trimmed.fastq" 2>>"$D/hbv_short.log" \
      | samtools sort -@ "$SORT_THREADS" -o "$D/hbv_short.umi.bam" -
    samtools index "$D/hbv_short.umi.bam"
    rm -f "$D/.short_in.fastq" "$D/.short_trimmed.fastq"

    NALN=$(samtools view -c -F 0x904 -q "$MIN_MAPQ" "$D/hbv_short.umi.bam")
    echo "[$(date +%T)]   stage 1 done in $(hms $(( SECONDS - S1 )))  HBV-aligned: $NALN"

    # ---- 28-nt UMIs ----
    echo "[$(date +%T)]   stage 2/3 filtering to 28-nt UMIs"
    samtools view -b -F 0x904 -q "$MIN_MAPQ" -e 'length([RX])==28' "$D/hbv_short.umi.bam" \
      | samtools sort -@ "$SORT_THREADS" -o "$D/.short_umi28.bam" -
    samtools index "$D/.short_umi28.bam"
    WF=$(samtools view -c "$D/.short_umi28.bam")
    echo "[$(date +%T)]   well-formed 28-nt UMIs: $WF"

    # ---- dedup ----
    UNIQ=0
    if [ "$WF" -gt 0 ]; then
        echo "[$(date +%T)]   stage 3/3 umi_tools dedup (--method=$UMI_METHOD)"
        umi_tools dedup -I "$D/.short_umi28.bam" -S "$D/.short_dedup.bam" \
            --extract-umi-method=tag --umi-tag=RX --method="$UMI_METHOD" \
            -L "$D/hbv_short_umi_dedup.log" >/dev/null 2>&1 || \
            echo "           WARN umi_tools returned non-zero"
        UNIQ=$(samtools view -c -F 0x904 "$D/.short_dedup.bam" 2>/dev/null || echo 0)
    fi
    DUP=$(awk -v u="$UNIQ" -v w="$WF" 'BEGIN{printf (w>0)?"%.4f":"NA",(w>0)?1-u/w:0}')

    # ---- coverage, for the positional comparison against the long fraction ----
    samtools depth -a -r "$HBV" "$D/hbv_short.umi.bam" > "$D/hbv_short_coverage.txt" 2>/dev/null || true

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample_key" "$barcode" "$NUNMAP" "$NSHORT" "$NPY" "$NALN" "$WF" "$UNIQ" "$DUP" >> "$SUMMARY"
    echo "[$(date +%T)] DONE $sample_key in $(hms $(( SECONDS - T0 )))"
    echo "           unmapped $NUNMAP -> short $NSHORT -> pychopper $NPY -> HBV $NALN -> UMI28 $WF -> molecules $UNIQ (dup $DUP)"
    echo
    rm -f "$D/.short_umi28.bam"* "$D/.short_dedup.bam"*
done

echo "=== Phase 2s complete ==="
column -t -s$'\t' "$SUMMARY" | grep -v '^#'

# ---------------------------------------------------------------------------
# positional coverage: short vs long, per library, binned to 50 bp
# ---------------------------------------------------------------------------
echo
echo "--- building long-vs-short coverage comparison"
python3 - "$ANALYSIS" "$SAMPLES" "$COVOUT" "$HBV" "$HBV_2X" "$HBV_LEN" <<'PYEOF'
import sys, subprocess, csv
from pathlib import Path
analysis, samples_tsv, out, hbv_1x, hbv_2x, hbv_len = sys.argv[1:7]
analysis = Path(analysis)
hbv_len = int(hbv_len)
BIN = 50

# The long fraction is aligned to the 2x contig (0-6363) and the short fraction to the 1x
# reference (0-3181). Coordinates are NOT comparable as-is — the long fraction's positions
# must be folded to single-copy space before the two profiles can be overlaid.

keys = []
for line in Path(samples_tsv).read_text().splitlines():
    if line.startswith("#") or line.startswith("barcode\t") or not line.strip():
        continue
    f = line.split("\t")
    if len(f) >= 4:
        keys.append(f[3])

def binned_depth(bam, contig, fold):
    """-> {bin_index: summed depth} in SINGLE-COPY space, or None.

    fold=True folds 2x coordinates into single-copy space (positions p and p+3182 are the
    same genomic base). Depth is summed rather than averaged across the two copies, since
    a wrapping read contributes to both."""
    if not Path(bam).exists():
        return None
    try:
        r = subprocess.run(["samtools", "depth", "-a", "-r", contig, str(bam)],
                           capture_output=True, text=True, check=True)
    except Exception:
        return None
    acc, cnt = {}, {}
    for line in r.stdout.splitlines():
        p = line.split("\t")
        if len(p) < 3:
            continue
        pos = int(p[1])
        if fold:
            pos = ((pos - 1) % hbv_len) + 1
        b = pos // BIN
        acc[b] = acc.get(b, 0) + int(p[2]); cnt[b] = cnt.get(b, 0) + 1
    return {b: acc[b] / max(1, cnt[b]) for b in acc}

rows = []
for k in keys:
    d = analysis / "samples" / k
    long_d = binned_depth(d / "hbv.umi.bam", hbv_2x, True)      # 2x -> folded
    short_d = binned_depth(d / "hbv_short.umi.bam", hbv_1x, False)  # already 1x
    if long_d is None and short_d is None:
        continue
    bins = sorted(set(long_d or {}) | set(short_d or {}))
    ltot = sum((long_d or {}).values()) or 1
    stot = sum((short_d or {}).values()) or 1
    for b in bins:
        lv = (long_d or {}).get(b, 0.0)
        sv = (short_d or {}).get(b, 0.0)
        rows.append({
            "sample_key": k, "bin_start": b * BIN,
            "long_depth": round(lv, 2), "short_depth": round(sv, 2),
            # normalised so the two fractions are comparable in shape, not magnitude
            "long_frac_of_total": round(lv / ltot, 6),
            "short_frac_of_total": round(sv / stot, 6),
        })

if rows:
    with open(out, "w", newline="") as fh:
        fh.write(f"# per-{BIN}bp mean depth in SINGLE-COPY coordinates (0-{hbv_len-1}); "
                 f"long = hbv.umi.bam on {hbv_2x} folded mod {hbv_len}; "
                 f"short = hbv_short.umi.bam on {hbv_1x}\n")
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print(f"    written: {out}  ({len(rows)} rows)")
    print()
    print("    Interpretation: compare long_frac_of_total against short_frac_of_total")
    print("    per bin. If the short fragments are RNase H products they should be")
    print("    CONCENTRATED in particular regions (tracking minus-strand synthesis from")
    print("    DR1 near 1820) rather than following the long-fragment profile. If the two")
    print("    normalised profiles overlay, the short fraction is just degraded versions")
    print("    of the same transcripts and carries no extra positional information.")
else:
    print("    no coverage data produced")
PYEOF
