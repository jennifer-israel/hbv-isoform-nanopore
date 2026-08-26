#!/usr/bin/env bash
# Demultiplex EXP26000892 custom PCR barcodes with dorado demux.
#
#   bash run_dorado_demux.sh /path/to/run /path/to/demux_out
#
# Runs dorado demux TWICE on purpose:
#   pass A  --barcode-both-ends   -> STRICT. Both ends must carry the same barcode.
#   pass B  (single-end allowed)  -> PERMISSIVE. One confident end is enough.
#
# Why two passes: the 4 libraries were pooled BEFORE Twist hyb capture and then
# co-amplified in the post-capture PCR, and they share identical P5/P7/flank
# sequences. That is textbook conditions for template switching / index hopping.
# (B minus A) is your empirical hopping+chimera rate. It matters enormously here:
# the 0.1 ng library got 24 PCR cycles, so even ~1% leakage from the 100 ng library
# could masquerade as real HBV signal at the lowest input and produce a false
# "we can go down to 0.1 ng" conclusion. STRICT (pass A) is what you should carry
# into the phase1-5 pipeline; pass B exists to quantify what STRICT threw away.
#
# Also reshapes dorado's flat per-barcode output into the nested directory layout
# phase1_align.sh already globs ($RUNDATA/fastq_pass/$barcode/*.fastq.gz).
set -euo pipefail

RUN_DIR=${1:?usage: run_dorado_demux.sh <run_dir> <out_dir>}
OUT_DIR=${2:?usage: run_dorado_demux.sh <run_dir> <out_dir>}
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

ARRANGEMENT=$SCRIPT_DIR/barcode_arrangement_EXP26000892.toml
SEQUENCES=$SCRIPT_DIR/barcode_sequences_EXP26000892.fasta
ARR_NAME=EXP26000892           # must match `name` and `kit` in the TOML
THREADS=${THREADS:-16}

# Input preference: bam_pass over fastq_pass. The BAMs carry the basecaller's tags
# (and dorado demux writes the BC tag back into BAM), so staying in BAM keeps the
# run auditable. Set INPUT_DIR explicitly to override.
INPUT_DIR=${INPUT_DIR:-}
if [ -z "$INPUT_DIR" ]; then
    if compgen -G "$RUN_DIR/bam_pass/*.bam" >/dev/null; then INPUT_DIR=$RUN_DIR/bam_pass
    elif compgen -G "$RUN_DIR/fastq_pass/*.fastq*" >/dev/null; then INPUT_DIR=$RUN_DIR/fastq_pass
    else echo "ERROR: no bam_pass/*.bam or fastq_pass/*.fastq* under $RUN_DIR" >&2; exit 1
    fi
fi

command -v dorado >/dev/null || { echo "ERROR: dorado not on PATH" >&2; exit 1; }
for f in "$ARRANGEMENT" "$SEQUENCES"; do
    [ -s "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done

echo "=== dorado demux — EXP26000892 ==="
dorado --version 2>&1 | sed 's/^/  dorado version: /'
echo "  input:       $INPUT_DIR"
echo "  arrangement: $ARRANGEMENT"
echo "  output:      $OUT_DIR"
echo

# ---------------------------------------------------------------------------
# NOTE ON --kit-name: it is MANDATORY. dorado demux refuses to run without either
# --kit-name or --no-classify:
#   "Please specify either --no-classify or --kit-name to use the demux tool."
# We pass our OWN arrangement name, which must match the `kit` field in the TOML.
#
# Re dorado#1548 (custom barcodes silently ignored when --kit-name is set): there,
# --kit-name named a real BUILT-IN ONT kit, so dorado preferred the built-in
# barcodes. Using a name that cannot collide with any built-in kit avoids that.
# verify_pass() below is the backstop that catches it regardless.
#
# NOTE ON BARCODE NAMES: they are MW01..MW04, not BC01..BC04. dorado ships built-in
# barcodes named BC01..BC96 and rejects a custom file that redefines those names
# ("Custom barcode names already exist"). Sequences are unchanged; only labels differ.
# ---------------------------------------------------------------------------

run_pass () {
    local label=$1 outsub=$2; shift 2
    local dest=$OUT_DIR/$outsub
    mkdir -p "$dest"
    echo "--- pass $label -> $dest"
    dorado demux \
        --kit-name "$ARR_NAME" \
        --barcode-arrangement "$ARRANGEMENT" \
        --barcode-sequences   "$SEQUENCES" \
        --output-dir "$dest" \
        --threads "$THREADS" \
        --emit-summary \
        --no-trim \
        "$@" \
        "$INPUT_DIR"
    echo
}

# Fail loudly if the custom arrangement did not actually take effect.
verify_pass () {
    local dest=$OUT_DIR/$1
    local summary
    summary=$(find "$dest" -maxdepth 1 -name '*summary*' -type f | head -1)
    if [ -z "$summary" ]; then
        echo "ERROR: no barcoding summary emitted in $dest — cannot verify." >&2; exit 1
    fi
    # The classified barcode names must be OUR names (MW01..MW04 / EXP26000892_*).
    # If dorado fell back to a built-in kit we would see e.g. barcode01 / SQK-* here.
    if ! grep -qE "(${ARR_NAME}_)?barcode0[1-9]" "$summary"; then
        echo "ERROR: summary $summary contains no barcode01-04 classifications." >&2
        echo "       The custom arrangement was probably ignored (see dorado#1548)." >&2
        echo "       Observed barcode values:" >&2
        awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="barcode") c=i; next} c{print $c}' "$summary" \
            | sort | uniq -c | sort -rn | head >&2
        exit 1
    fi
    echo "  verify OK: custom barcode names present in $(basename "$summary")"
}

mkdir -p "$OUT_DIR"

run_pass "A (STRICT, both ends)" strict --barcode-both-ends
verify_pass strict

run_pass "B (PERMISSIVE, single end ok)" permissive
verify_pass permissive

# ---------------------------------------------------------------------------
# Reshape STRICT output into the layout phase1_align.sh expects.
# dorado writes flat files (one per barcode) named after the arrangement, e.g.
# EXP26000892_MW01.bam / .fastq. Discover them rather than hardcoding, since the
# exact naming has shifted between dorado releases.
# ---------------------------------------------------------------------------
echo "--- reshaping STRICT output into fastq_pass/<barcode>/ layout"
PIPE_DIR=$OUT_DIR/pipeline_input/fastq_pass
mkdir -p "$PIPE_DIR"

shopt -s nullglob
for f in "$OUT_DIR"/strict/*MW0[1-4]* "$OUT_DIR"/strict/unclassified.*; do
    base=$(basename "$f")
    if [[ "$base" == unclassified.* ]]; then
        bc=unclassified
    else
        mw=$(echo "$base" | grep -oE 'MW0[1-4]')
        # MW01 -> custom_bc01, so samples_EXP26000892.tsv's barcode column is unchanged
        bc="custom_bc${mw#MW}"
    fi
    mkdir -p "$PIPE_DIR/$bc"
    case "$f" in
        *.bam)
            samtools fastq -@ 4 -n "$f" 2>/dev/null | gzip -1 > "$PIPE_DIR/$bc/${bc}.fastq.gz" ;;
        *.fastq.gz)
            cp "$f" "$PIPE_DIR/$bc/${bc}.fastq.gz" ;;
        *.fastq)
            gzip -1 -c "$f" > "$PIPE_DIR/$bc/${bc}.fastq.gz" ;;
    esac
    echo "  $base -> fastq_pass/$bc/"
done
shopt -u nullglob

# ---------------------------------------------------------------------------
# Hopping / ambiguity report: STRICT vs PERMISSIVE.
# ---------------------------------------------------------------------------
echo
echo "--- STRICT vs PERMISSIVE (index-hopping / chimera estimate)"
python3 - "$OUT_DIR" <<'PYEOF'
import sys, glob, os, csv
out = sys.argv[1]

def load(sub):
    """read_id -> barcode, from whichever summary file dorado emitted."""
    hits = glob.glob(os.path.join(out, sub, "*summary*"))
    if not hits:
        return None
    d = {}
    with open(hits[0]) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        idc = next((c for c in (r.fieldnames or []) if c.lower() in ("read_id","readid")), None)
        bcc = next((c for c in (r.fieldnames or []) if c.lower() == "barcode"), None)
        if not idc or not bcc:
            print(f"  (could not find read_id/barcode columns in {hits[0]}; "
                  f"saw {r.fieldnames})")
            return None
        for row in r:
            d[row[idc]] = row[bcc]
    return d

S, P = load("strict"), load("permissive")
if not S or not P:
    print("  skipped (summaries unavailable)")
    raise SystemExit

def norm(b):
    return "unclassified" if not b or "unclassified" in b.lower() else b

tot = len(P)
s_class = sum(1 for v in S.values() if norm(v) != "unclassified")
p_class = sum(1 for v in P.values() if norm(v) != "unclassified")

only_p = [k for k in P if norm(P[k]) != "unclassified" and norm(S.get(k, "")) == "unclassified"]
disagree = [k for k in P if norm(P[k]) != "unclassified"
            and norm(S.get(k, "")) not in ("unclassified", norm(P[k]))]

print(f"  total reads:                      {tot:>12,}")
print(f"  classified, STRICT (both ends):   {s_class:>12,}  ({100*s_class/tot:5.2f}%)")
print(f"  classified, PERMISSIVE (1 end):   {p_class:>12,}  ({100*p_class/tot:5.2f}%)")
print(f"  single-end-only (B not A):        {len(only_p):>12,}  ({100*len(only_p)/tot:5.2f}%)")
print(f"  STRICT/PERMISSIVE disagreement:   {len(disagree):>12,}  ({100*len(disagree)/tot:5.2f}%)")
print()
print("  Per-barcode read counts:")
print(f"    {'barcode':<22}{'STRICT':>12}{'PERMISSIVE':>14}{'ratio':>9}")
for bc in sorted({norm(v) for v in P.values()}):
    s = sum(1 for v in S.values() if norm(v) == bc)
    p = sum(1 for v in P.values() if norm(v) == bc)
    ratio = f"{s/p:.3f}" if p else "-"
    print(f"    {bc:<22}{s:>12,}{p:>14,}{ratio:>9}")
print()
print("  Interpretation: 'single-end-only' reads are ones where only one end gave a")
print("  confident barcode. Some are ordinary read truncation; some are chimeras that")
print("  genuinely carry two different barcodes. A LOW STRICT/PERMISSIVE ratio that is")
print("  concentrated in the low-input barcodes (MW03 1ng / MW04 0.1ng) is the warning")
print("  sign for leakage from the high-input libraries -- exactly the artifact that")
print("  would fake a positive result at 0.1 ng. Carry STRICT forward.")
PYEOF

echo
echo "=== done ==="
echo "Feed the pipeline:  RUNDATA=$OUT_DIR/pipeline_input"
echo "Sample sheet:       samples_EXP26000892.tsv  (barcode column = custom_bc01..04)"
