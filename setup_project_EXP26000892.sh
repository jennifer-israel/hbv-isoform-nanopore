#!/usr/bin/env bash
# Wire the EXP26000892 demux output into the phase1-5 pipeline layout.
#
#   bash setup_project_EXP26000892.sh
#
# Creates the directory structure the phase scripts expect, and repoints their
# hardcoded PROJECT_ROOT (currently EXP26000559) at this experiment.
#
# Layout produced:
#   /data/EXP26000896/                     <- PROJECT_ROOT
#     rundata/fastq_pass/custom_bc0N/      -> symlink to demux/pipeline_input
#     analysis/refs/                       -> symlink to /data/refs (phase 0 output)
#     analysis/{samples,comparison,logs,reports}/
#     config/samples.tsv                   <- from samples_EXP26000892.tsv
#   /data/shared/hbv_transcript_classify.py   <- phase3 imports from PROJECT_ROOT.parent/shared
set -euo pipefail

PROJECT_ROOT=/data/EXP26000896
SCRIPTS=/data/scripts
REFS=/data/refs
DEMUX_INPUT=$PROJECT_ROOT/demux/pipeline_input
OLD_ROOT='/home/ubuntu/matt_wolpert_claude_code_analysis/2026_07_02_EXP26000559_cDNA001'

echo "=== EXP26000892 project setup ==="

# ---- 1. directories ----
mkdir -p "$PROJECT_ROOT"/analysis/{samples,comparison,logs,reports} "$PROJECT_ROOT"/config /data/shared

# ---- 2. rundata -> demuxed fastqs ----
# phase1_align.sh globs $RUNDATA/fastq_pass/$barcode/*.fastq.gz
if [ ! -e "$PROJECT_ROOT/rundata" ]; then
    ln -s "$DEMUX_INPUT" "$PROJECT_ROOT/rundata"
    echo "  rundata -> $DEMUX_INPUT"
else
    echo "  rundata already exists, leaving as-is"
fi

# ---- 3. analysis/refs -> phase 0 output ----
if [ ! -e "$PROJECT_ROOT/analysis/refs" ]; then
    ln -s "$REFS" "$PROJECT_ROOT/analysis/refs"
    echo "  analysis/refs -> $REFS"
fi

# ---- 4. sample sheet ----
if [ -s "$SCRIPTS/samples_EXP26000892.tsv" ]; then
    cp "$SCRIPTS/samples_EXP26000892.tsv" "$PROJECT_ROOT/config/samples.tsv"
    echo "  config/samples.tsv installed"
else
    echo "  WARN: $SCRIPTS/samples_EXP26000892.tsv not found — copy it there first" >&2
fi

# ---- 5. repoint PROJECT_ROOT in every phase script ----
echo
echo "--- repointing PROJECT_ROOT in phase scripts (backups: *.bak) ---"
for f in "$SCRIPTS"/phase*.sh "$SCRIPTS"/phase*.py; do
    [ -e "$f" ] || continue
    if grep -q "$OLD_ROOT" "$f"; then
        cp -n "$f" "$f.bak"
        sed -i "s|$OLD_ROOT|$PROJECT_ROOT|g" "$f"
        sed -i "s|EXP26000559_cDNA001|EXP26000892_cDNA003|g" "$f"
        echo "  patched $(basename "$f")"
    fi
done

# ---- 6. verification ----
echo
echo "=== verification ==="
MMI=$PROJECT_ROOT/analysis/refs/hg38_hbv_2x_splice.mmi
FA=$PROJECT_ROOT/analysis/refs/hg38_hbv_2x.fa
[ -s "$MMI" ] && echo "  splice index      OK  ($(du -h "$MMI" | cut -f1))" || echo "  splice index      MISSING"
[ -s "$FA" ]  && echo "  composite fasta   OK  ($(du -h "$FA" | cut -f1))"  || echo "  composite fasta   MISSING"

if [ -s "$FA.fai" ]; then
    if grep -q '^U95551.1_2x' "$FA.fai"; then
        LEN=$(awk '$1=="U95551.1_2x"{print $2}' "$FA.fai")
        echo "  HBV contig        OK  U95551.1_2x, ${LEN} bp (expect 6364)"
    else
        echo "  HBV contig        MISSING from .fai — check the reference build" >&2
    fi
fi

echo "  sample sheet:"
awk -F'\t' '$1!~/^#/ && $1!="barcode" && NF>=4 {printf "    %-14s %-12s %-28s input=%-6s pcr=%s\n",$1,$2,$4,$5,$7}' \
    "$PROJECT_ROOT/config/samples.tsv" 2>/dev/null

echo "  input fastqs:"
for d in "$PROJECT_ROOT"/rundata/fastq_pass/custom_bc0*; do
    [ -d "$d" ] || continue
    n=$(ls "$d"/*.fastq.gz 2>/dev/null | wc -l)
    sz=$(du -sh "$d" 2>/dev/null | cut -f1)
    printf "    %-14s %s files, %s\n" "$(basename "$d")" "$n" "$sz"
done

echo
echo "  NOTE: phase3_classify.py imports hbv_transcript_classify.py from"
echo "        /data/shared/ — that file is NOT present. Phase 3 will fail until"
echo "        it is copied from the EXP26000559 project. Phases 1-2 are unaffected."
echo
echo "=== next ==="
echo "  conda activate hbv_lr_analysis"
echo "  bash $SCRIPTS/phase1_align.sh SeqLib5553_0.1ng_polyA_24    # smallest library first"
