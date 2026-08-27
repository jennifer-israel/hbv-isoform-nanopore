# From sequencing run to phase 1: demultiplexing and QC

How to take a nanopore run off S3 and get it ready for the phase 1 analysis scripts.

These libraries use **custom PCR barcodes** carried inside the cDNA construct, not ONT
native barcoding. MinKNOW therefore does no demultiplexing — every read lands in
`bam_pass/` unclassified — and we demultiplex computationally with `dorado demux` using a
custom barcode arrangement.

Everything below stops short of `phase1_align.sh`. Runs covered so far: cDNA003
(EXP26000896), cDNA005 (EXP26000993), cDNA006 (EXP26001054), cDNA007 (EXP26001055).

---

## Before you start

### Set your paths

Examples below use the paths from the original analysis instance. Set these three to
match wherever you're working and the commands transfer unchanged:

```bash
export EXP=EXP26001055                      # experiment ID, change per run
export SCRIPTS=/path/to/scripts             # where the pipeline scripts live
export DATA=/data                           # parent for per-experiment directories
```

Everything else derives from those. `$DATA` needs to be on a volume with room for roughly
2.5× the raw run size — a 70 GB run wants about 175 GB free, since dorado writes new BAMs
rather than filtering in place and the reshaped FASTQs are another copy.

### Build the analysis environment

The scripts expect a conda environment with these tools. The Python version matters:
3.13 pulls an ancient `umi_tools` build that fails at import, so pin 3.11.

```bash
conda create -n hbv_lr -c conda-forge -c bioconda \
  python=3.11 \
  minimap2 'samtools>=1.16' pychopper \
  'umi_tools=1.1.6' pysam \
  'pandas<3' numpy scipy statsmodels matplotlib pyarrow \
  setuptools
conda activate hbv_lr
pip install edlib          # no conda build for py3.11; pip works
```

Notes from setting this up the first time:

- `samtools` must be ≥1.16 — phase 2 checks and exits if not.
- `setuptools` is needed explicitly, otherwise `umi_tools` fails with
  `ModuleNotFoundError: No module named 'pkg_resources'`.
- `pyarrow` is what lets the phase 3/4 scripts write `.parquet`.
- If `pip install` hits PEP 668 (`externally-managed-environment`), the environment has no
  pip of its own — install it with conda first.

### Install dorado separately

`dorado` is **not** in the conda environment. Download the ONT binary and put it on PATH:

```bash
# check what you have
which dorado && dorado --version
```

If it's missing, get the current release from Oxford Nanopore's dorado GitHub releases,
unpack it, and add its `bin/` to PATH. **Version matters:** dorado 2.x requires
`--kit-name` even with a fully custom barcode arrangement, whereas some 0.x builds did not.
The instructions here assume 2.x.

### Reference files

Phase 1 needs a composite hg38 + 2×HBV minimap2 index at
`$DATA/$EXP/analysis/refs/hg38_hbv_2x_splice.mmi`. In the original setup `analysis/refs`
is a symlink to a shared `/data/refs` directory holding:

```
hg38_hbv_2x.fa  hg38_hbv_2x.fa.fai  hg38_hbv_2x_splice.mmi
hbv_2x_only.fa  hbv_2x_only.mmi
```

The HBV component is U95551.1 doubled end-to-end to 6,364 bp, so transcripts that wrap the
circular genome align contiguously. `phase0_build_ref.sh` builds these; budget time and
about 20 GB. Nothing in *this* document needs them — they're required from phase 1 onward.

### Verify before starting

```bash
conda activate hbv_lr
which dorado samtools minimap2 && samtools --version | head -1
df -h $DATA | tail -1
```

### What you need from Benchling

| what | why |
|---|---|
| library ↔ barcode sequence table | builds the barcode FASTA and sample sheet |
| input amount per library (ng) | required for any recovery-per-input comparison |
| PCR cycle count | required to interpret duplication across runs |
| poly(A) route (enzymatic vs oligo-d(T)) | confounds cross-library comparison |
| what any spike-in or carrier actually is | determines whether it maps to the reference |

---

## Directory layout

Set this up first and the rest follows. The single most important rule: **`rundata/`
contains only demultiplexed reads.** Phase 1 globs `$PROJECT_ROOT/rundata/fastq_pass/<barcode>/*.fastq.gz`,
so if the S3 sync also lands there you will end up bridging paths with symlinks.

```
$DATA/$EXP/
  raw/                        S3 sync: bam_pass/, other_reports/, sample_sheet CSV
  demux/                      dorado demux output
    strict/sequencing_summary.txt   where phase1_5 looks for the summary
  rundata/fastq_pass/<bc>/    demultiplexed reads — NOTHING ELSE IN rundata/
  config/samples.tsv          barcode → sample map
  analysis/
    samples/  comparison/  reports/  logs/
    refs -> shared reference directory
```

```bash
mkdir -p $DATA/$EXP/{raw,config,rundata,demux/strict}
mkdir -p $DATA/$EXP/analysis/{samples,comparison,reports,logs}
ln -sfn $DATA/refs $DATA/$EXP/analysis/refs     # point at your reference directory
```

---

## Step 1 — Get the data off S3

Inventory before copying. `bam_pass` and `fastq_pass` are the same reads in two formats;
we use BAM, so skip `fastq_pass` and both fail folders. That typically halves the transfer.

```bash
S3=s3://tunetx-raw-ngs-.../nanopore/<run_folder>/.../<run_dir>

aws s3 ls --recursive --summarize "$S3/" | tail -3
aws s3 ls --recursive "$S3/" | awk '{n=split($4,p,"/"); d=p[n-1]; s[d]+=$3; c[d]++}
  END{for(k in s) printf "%9.1f GB  %5d files  %s\n", s[k]/1073741824, c[k], k}' | sort -rn
df -h $DATA | tail -1
```

```bash
aws s3 sync "$S3/" $DATA/$EXP/raw/ \
  --exclude "*" --include "bam_pass/*" --include "other_reports/*" --include "*.csv" \
  --no-progress
```

Raise concurrency if the transfer is slow (default is 10):

```bash
aws configure set default.s3.max_concurrent_requests 32
aws configure set default.s3.max_queue_size 10000
```

### Verify the transfer

```bash
# nothing to copy means everything is present at matching size
aws s3 sync "$S3/" $DATA/$EXP/raw/ --exclude "*" --include "bam_pass/*" --dryrun

# file count and bytes
aws s3 ls --recursive "$S3/bam_pass/" | awk '{n++; b+=$3} END{printf "S3:    %d files %d bytes\n",n,b}'
find $DATA/$EXP/raw/bam_pass -name '*.bam' -printf '%s\n' | awk '{n++; b+=$1} END{printf "local: %d files %d bytes\n",n,b}'

# BAMs intact — catches truncation, which a size match will not
samtools quickcheck -v $DATA/$EXP/raw/bam_pass/*.bam && echo "all BAMs OK"
```

---

## Step 2 — Build the barcode configuration

Two files, both in `$SCRIPTS`.

### The barcode FASTA

One entry per library, in the order they appear in the Benchling table. Barcodes are 24 nt.

```
>TY01
AGAACGACTTCCATACTCGTGTGA
>TY02
AACGAGTCTCTTGGGACCCATAGA
...
```

**Barcode names must not collide with dorado's built-ins.** Dorado ships the standard ONT
barcode sets and rejects the arrangement outright with `Custom barcode names already exist`
if you reuse a name. Our barcodes frequently *are* ONT catalogue sequences. Use a fresh
prefix each run:

| run | prefix | note |
|---|---|---|
| cDNA003 / cDNA005 | `MW01`–`MW04` | sequences are ONT `BC01`–`BC04` |
| cDNA006 | `TX01`–`TX06` | reverse complements of `BC05`–`BC10` |
| cDNA007 | `TY01`–`TY06` | new sequences |

Never name them `BC..`, `NB..`, `RB..`, `BP..` or `RLB..`.

**Orientation.** Benchling has recorded these inconsistently — cDNA003/005 in ONT-canonical
orientation, cDNA006 as reverse complements. Dorado searches reads in both orientations so
either resolves, but enter them exactly as recorded and confirm empirically (Step 5). A
systematic orientation error produces near-zero assignment with no error message.

### The arrangement TOML

This tells dorado where in the read to look and what flanks the barcode. Complete working
example — `barcode_arrangement_EXP26001055.toml`, as used for cDNA007:

```toml
[arrangement]
name = "EXP26001055"
kit  = "EXP26001055"

mask1_front = "AATGATACGGCGACCACCGAGGTGCTG"
mask1_rear  = "TTAACCTTTTCTGTTGGTGCTGATATTGC"
mask2_front = "CAAGCAGAAGACGGCATACGAGATGGTGCTG"
mask2_rear  = "TTAACCTCTTGCCTGTCGCTCTATC"

barcode1_pattern = "TY%02i"
barcode2_pattern = "TY%02i"
first_index = 1
last_index  = 6

[scoring]
max_barcode_penalty      = 6
min_barcode_penalty_dist = 4
min_separation_only_dist = 6
flank_left_pad           = 5
flank_right_pad          = 5
front_barcode_window     = 200
rear_barcode_window      = 200
barcode_end_proximity    = 200
min_flank_score          = 0.6
midstrand_flank_score    = 0.95
```

Only four things change between runs:

| field | what to set |
|---|---|
| `name`, `kit` | the experiment ID. `kit` is what you pass to `--kit-name`, and it must not match a built-in ONT kit name |
| `barcode1_pattern`, `barcode2_pattern` | the FASTA name prefix plus `%02i`, e.g. `TY%02i` for `TY01`…`TY06`. Both are the same because the same barcode appears at both ends |
| `last_index` | number of libraries |
| `first_index` | leave at 1 |

**The masks are the library construct** — the P5/P7, SSPII and CRTA sequences flanking the
barcode, front and rear, on each of the two adapter ends. They were verified base-for-base
against the cDNA003 oligos and only change if the construct itself changes. Blocking-oligo
changes (cDNA007) do not affect them. If demux assignment collapses across *all* barcodes
at once, suspect the masks before touching anything else.

For a run that follows a previous one, deriving the file is faster than writing it:

```bash
cd $SCRIPTS
sed -e 's/EXP26001054/EXP26001055/g' -e 's/TX%02i/TY%02i/g' \
    barcode_arrangement_EXP26001054.toml > barcode_arrangement_EXP26001055.toml
grep -E 'name|kit|pattern|last_index' barcode_arrangement_EXP26001055.toml
```

Always grep the result — a missed substitution leaves the previous run's kit name in place,
and dorado will run happily with the wrong barcode pattern and assign nothing.

> **Do not tune the `[scoring]` block between runs.** Those thresholds — particularly the
> 200 bp barcode windows — are shared with the concatemer analysis, which detects
> adapter/barcode motifs more than 200 bp from both read ends. Changing them makes
> concatemer rates non-comparable across experiments, and that comparison is the primary
> endpoint of some runs.

---

## Step 3 — Pre-flight: are the adapters still there?

**Do this before every demux.** If MinKNOW trimmed adapters during basecalling, the custom
barcodes are gone from the read ends and no configuration can recover them — the run has to
be re-basecalled from POD5 with `--no-trim`.

```bash
B=$(ls $DATA/$EXP/raw/bam_pass/*.bam | head -1)
samtools view "$B" | head -50000 | awk '{print $10}' \
  | grep -c "AATGATACGGCGACCACCGAG"
```

Interpretation: an exact 21-mer grep on ONT reads succeeds roughly 0.97²¹ ≈ 53% of the time
by chance, so **~25,000 out of 50,000 means essentially every read retains its adapter.**
Near zero means trimming happened. Stop and talk to the sequencing team.

---

## Step 4 — Run dorado demux

```bash
cd $SCRIPTS
screen -dmS demux bash -c "
dorado demux --kit-name $EXP \
  --barcode-arrangement $SCRIPTS/barcode_arrangement_$EXP.toml \
  --barcode-sequences  $SCRIPTS/barcode_sequences_$EXP.fasta \
  --barcode-both-ends \
  --no-trim \
  --emit-summary \
  --output-dir $DATA/$EXP/demux \
  $DATA/$EXP/raw/bam_pass/ \
  2>&1 | tee $DATA/$EXP/demux.log; exec bash"
screen -ls
```

Reattach with `screen -r demux`, detach with `Ctrl-A` then `d`. The job survives a dropped
SSH session, which matters — these runs take hours.

Flags that are not optional:

- **`--kit-name`** — dorado 2.x requires it even with a custom arrangement. Pass the
  `kit` value from the TOML.
- **`--barcode-both-ends`** — requires the same barcode at both ends. This rejects fusions
  of molecules from *different* libraries, which is a large part of why we trust the
  assignments. It costs assignment rate (see below) and is worth it.
- **`--no-trim`** — keeps adapters in the output so downstream concatemer detection can
  still see them.
- **`--emit-summary`** — writes `sequencing_summary.txt`, which phase 1.5 needs for read
  length and quality.

### Watch disk

Dorado writes new BAMs rather than filtering in place, so the demux output is roughly the
size of the input. With the reshaped FASTQs on top, budget about 2.5× the input size.

---

## Step 5 — QC the demux

### 5a. Assignment rate

```bash
D=$(find $DATA/$EXP/demux -type d -name 'bam_pass' | head -1)
TOT=0; ASG=0
for b in $(ls $D); do
  n=$(samtools view -c $D/$b/*.bam 2>/dev/null || echo 0)
  printf "%-14s %12d\n" "$b" "$n"
  TOT=$((TOT+n)); [ "$b" != "unclassified" ] && ASG=$((ASG+n))
done
echo "assigned: $ASG / $TOT = $(echo "scale=1; 100*$ASG/$TOT" | bc)%"
```

Observed so far, for reference:

| run | assignment rate | note |
|---|---|---|
| cDNA006 | 48.1% | short mouse serum fragments |
| cDNA007 | 75.4% | intracellular PHH, new blocking oligos |

Both-ends demultiplexing costs assignment rate by design, so 50–75% is normal. A low rate
is not necessarily an error, but it is informative: barcodes are harder to detect on short
reads, and mid-read barcode/chimera events cause the two ends to disagree. A rate far below
the range above, across *all* barcodes, suggests the masks don't match the construct.

Also check the distribution. A pool where one library holds most of the reads will suppress
the others during capture — cDNA005 showed a library dropping from 2.58% to 0.27% HBV
purely from pool composition, and cDNA006 had a 167× spread across six libraries.

### 5b. Confirm which library landed where — do not skip this

Dorado names output directories `barcode01`…`barcodeNN` positionally. That the first FASTA
entry maps to `barcode01` is an assumption, and the demux summary cannot confirm it because
it uses the same positional names. A transposition silently attaches every result to the
wrong sample.

The check is to count which barcode *sequence* appears in each directory's reads:

```bash
cd $D
declare -A F=( [TY01]=AGAACGACTTCCATACTCGTGTGA [TY02]=AACGAGTCTCTTGGGACCCATAGA \
               [TY03]=AGGTCTACCTCGCTAACACCACTG [TY04]=CGTCAACTGACAGTGGTTCGTACT \
               [TY05]=ACCCTCCAGGAAAGTACCTCTGAT [TY06]=CCAAACCCAACAACCTAGATAGGC )
declare -A R=( [TY01]=TCACACGAGTATGGAAGTCGTTCT [TY02]=TCTATGGGTCCCAAGAGACTCGTT \
               [TY03]=CAGTGGTGTTAGCGAGGTAGACCT [TY04]=AGTACGAACCACTGTCAGTTGACG \
               [TY05]=ATCAGAGGTACTTTCCTGGAGGGT [TY06]=GCCTATCTAGGTTGTTGGGTTTGG )

printf "%-12s" "dir"; for t in TY0{1..6}; do printf "%8s" $t; done; echo
for i in 1 2 3 4 5 6; do
  samtools view barcode0$i/*.bam | head -50000 | awk '{print $10}' > /tmp/s.txt
  printf "%-12s" "barcode0$i"
  for t in TY0{1..6}; do printf "%8d" "$(grep -c -e "${F[$t]}" -e "${R[$t]}" /tmp/s.txt)"; done
  echo
done; rm -f /tmp/s.txt
```

`R` holds the reverse complements — reads come off in both orientations, so both must be
counted. Expect a clean diagonal near 48,000–49,500 out of 50,000, with off-diagonal
entries in the tens or low hundreds. Anything off-diagonal and large means the sample sheet
needs transposing before you go further.

A dominant library's barcode showing up in all the others at a few tenths of a percent is
normal low-level bleed and not a transposition.

---

## Step 6 — Reshape and place the outputs

Dorado writes a nested structure: `<sample>/<run>/bam_pass/barcodeNN/`.
`reshape_demux_output.sh` flattens it and renames barcodes to `custom_bcNN`, which is what
`config/samples.tsv` and the phase scripts use.

```bash
cd $SCRIPTS
PROJECT_ROOT=$DATA/$EXP ./reshape_demux_output.sh $DATA/$EXP/demux /tmp/rs
mv /tmp/rs/fastq_pass $DATA/$EXP/rundata/fastq_pass
cp $DATA/$EXP/demux/sequencing_summary.txt $DATA/$EXP/demux/strict/
ls $DATA/$EXP/rundata/fastq_pass/
```

The script creates a `fastq_pass/` subdirectory *inside* whatever output directory you give
it, so pass a scratch path and move the result — don't pass `$DATA/$EXP/rundata/fastq_pass`
or you get `fastq_pass/fastq_pass`.

> Both source and destination must be on `/data`. Moving tens of GB via `/tmp` fails with
> `Disk quota exceeded` — `/tmp` is a small tmpfs, not the data volume. Moves within `/data`
> are instant renames.

---

## Step 7 — Write the sample sheet

`$DATA/$EXP/config/samples.tsv`, tab-separated, eight columns, `barcode` as `custom_bcNN`:

```
barcode	lib_id	sample_name	sample_key	input_ng	polya	pcr	barcode_seq
custom_bc01	SeqLib5666	Blocker_v1_10ng	SeqLib5666_blockerV1_10ng	10	oligo_dT	NA	AGAACGACTTCCATACTCGTGTGA
```

`sample_key` names the output directory under `analysis/samples/`, so keep it stable —
renaming it later orphans everything already computed. Comment lines starting with `#` are
ignored, and are the right place to record what Benchling didn't tell you.

If you ever need to rebuild this sheet, derive the pairing from what was actually processed
rather than from Benchling:

```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next} {print $h["barcode"]"\t"$h["sample_key"]}' \
  $DATA/$EXP/analysis/comparison/phase1_align_summary.tsv
```

---

## Step 8 — Hand off to phase 1

```bash
cd $SCRIPTS && conda activate hbv_lr
screen -S p1
PROJECT_ROOT=$DATA/$EXP ./phase1_align.sh
PROJECT_ROOT=$DATA/$EXP ./phase1_5_readqc_demux.sh
```

**Check `PROJECT_ROOT` is actually overridable first.** Several phase scripts hardcode it
and silently ignore the environment variable, which means they run against the *previous*
experiment and produce plausible-looking output for the wrong samples:

```bash
grep -n '^PROJECT_ROOT' $SCRIPTS/*.sh
```

Every line should read `PROJECT_ROOT="${PROJECT_ROOT:-/data/EXP...}"`. To fix one:

```bash
cp -n script.sh script.sh.bak
sed -i 's|^PROJECT_ROOT=.*|PROJECT_ROOT="${PROJECT_ROOT:-/data/EXP26000993}"|' script.sh
```

Then sanity-check the first output table: **the sample keys must belong to the experiment
you meant to run.** This is how a phase 1.5 run silently reported cDNA005 results during
cDNA006 processing.

---

## Troubleshooting

| symptom | cause | fix |
|---|---|---|
| `Please specify either --no-classify or --kit-name` | dorado 2.x requires `--kit-name` even for custom arrangements | pass the TOML's `kit` value |
| `Custom barcode names already exist: BC01 …` | barcode names collide with dorado built-ins | rename with an unused prefix |
| Near-zero assignment across all barcodes | masks don't match the construct, or adapters were trimmed at basecalling | run the Step 3 adapter check |
| `ERROR: summary contains no <NAME> classifications` | `run_dorado_demux.sh`'s `verify_pass()` greps barcode *names* while dorado writes `barcodeNN` directories | known false alarm; ignore |
| `no fastq_pass files — skipping` for every library | `rundata/fastq_pass/<barcode>/` path wrong, or files are `.bam` not `.fastq.gz` | check the reshape output landed where phase 1 globs |
| `phase1.log: No such file or directory` | `tee` target's directory doesn't exist yet | `mkdir -p $DATA/$EXP/analysis/logs` |
| `demux summary not found: .../demux/strict/sequencing_summary.txt` | ran dorado directly rather than via `run_dorado_demux.sh` | copy `demux/sequencing_summary.txt` into `demux/strict/` |
| `cp: cannot stat ''` | a `find` inside `$( )` returned nothing | locate the file first: `find $DATA/$EXP/demux -type f ! -name '*.bam'` |
| `mv: Disk quota exceeded` | moving large files through `/tmp` | keep moves within `/data` |
| Results attributed to the wrong experiment | script hardcodes `PROJECT_ROOT` | see Step 8 |
| `ExpiredToken` / `NoCredentials` on S3 | temporary AWS credentials lapsed | re-login; if SSO federates through Microsoft, a password problem breaks both |

---

## What to record for each run

Keep these with the run — they are what make cross-experiment comparison possible:

- assignment rate and per-barcode read counts
- the barcode diagonal, confirming library identity
- adapter retention count from Step 3
- read length N50 and median per library (from phase 1.5)
- barcode prefix used, so the next run picks a fresh one
- confirmation that the `[scoring]` block was unchanged

Anything that changes detection thresholds should be recorded loudly, because concatemer
and duplication rates are compared across runs and those comparisons assume the thresholds
held constant.
