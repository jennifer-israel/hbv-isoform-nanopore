#!/usr/bin/env bash
# Phase 0 — build the hg38 + 2×HBV composite reference and its minimap2 splice index.
#
#   conda activate hbv_lr_analysis
#   bash scripts/phase0_build_ref.sh
#
# Idempotent: each step is skipped if its output already exists and looks valid.
# Run inside tmux — index build reads/writes ~3 GB and takes several minutes.
#
# Outputs (analysis/refs/):
#   hg38_hbv_2x.fa            hg38 contigs + U95551.1_2x (6,364 bp)
#   hg38_hbv_2x.fa.fai        samtools faidx
#   hg38_hbv_2x_splice.mmi    minimap2 splice index (k15/w5, matches `-ax splice`)
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/data/EXP26000896}"
ANALYSIS_ROOT=$PROJECT_ROOT/analysis
REF_SRC=/home/ubuntu/matt_wolpert_claude_code_analysis/2026_05_14_HBV_LR_Transcript_Detect_EXP26000465/refs/hg38_hbv.fa
REF_2X=$ANALYSIS_ROOT/refs/hg38_hbv_2x.fa
MMI_INDEX=$ANALYSIS_ROOT/refs/hg38_hbv_2x_splice.mmi
LOG=$ANALYSIS_ROOT/logs/phase0_build_ref.log

HBV_2X_CONTIG=U95551.1_2x
HBV_2X_LEN=6364

mkdir -p "$ANALYSIS_ROOT/refs" "$ANALYSIS_ROOT/logs"
exec > >(tee -a "$LOG") 2>&1
echo "[$(date)] Phase 0 start. nproc=$(nproc)"

command -v minimap2 >/dev/null || { echo "ERROR: minimap2 not on PATH (conda activate hbv_lr_analysis)"; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools not on PATH"; exit 1; }
[ -f "$REF_SRC" ] || { echo "ERROR: source composite not found: $REF_SRC"; exit 1; }

# 1. Build the 2× FASTA -------------------------------------------------------
if [ -s "$REF_2X" ] && grep -q "^>$HBV_2X_CONTIG" "$REF_2X"; then
    echo "[$(date)] 2× FASTA exists, skipping build: $REF_2X"
else
    echo "[$(date)] Building 2× FASTA ..."
    python "$PROJECT_ROOT/scripts/make_hbv_2x_ref.py" --in-fasta "$REF_SRC" --out-fasta "$REF_2X"
fi

# 2. faidx + verify the doubled HBV contig length ----------------------------
[ -f "$REF_2X.fai" ] || { echo "[$(date)] samtools faidx ..."; samtools faidx "$REF_2X"; }
GOT_LEN=$(awk -v c="$HBV_2X_CONTIG" '$1==c{print $2}' "$REF_2X.fai")
if [ "$GOT_LEN" != "$HBV_2X_LEN" ]; then
    echo "ERROR: $HBV_2X_CONTIG length is '$GOT_LEN', expected $HBV_2X_LEN"; exit 1
fi
echo "[$(date)] Verified $HBV_2X_CONTIG = $GOT_LEN bp."

# 3. minimap2 splice index (k15/w5 to match align-time `-ax splice`) ----------
if [ -s "$MMI_INDEX" ]; then
    echo "[$(date)] Splice index exists, skipping: $MMI_INDEX"
else
    echo "[$(date)] Building minimap2 splice index ..."
    minimap2 -x splice -d "$MMI_INDEX" "$REF_2X"
fi

echo "[$(date)] Phase 0 complete."
echo "  REF_2X=$REF_2X ($(du -h "$REF_2X" | cut -f1))"
echo "  MMI_INDEX=$MMI_INDEX ($(du -h "$MMI_INDEX" | cut -f1))"
