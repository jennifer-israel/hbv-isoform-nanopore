#!/usr/bin/env bash
# Phase 1 — align every barcode's PASS reads to hg38+2×HBV (cDNA splice mode).
#
# Headline alignment = ALL pass reads (max HBV-detection sensitivity). The
# pychopper-tagged / UMI alignment for duplication metrics is Phase 1.5, separate.
#
#   conda activate hbv_lr_analysis
#   bash scripts/phase1_align.sh                 # all libraries in config/samples.tsv
#   bash scripts/phase1_align.sh SeqLib5543_150ng_NOpolyA_20   # one (sample_key or barcode)
#   FORCE=1 bash scripts/phase1_align.sh <key>   # re-run even if BAM exists
#
# Per library (analysis/samples/<sample_key>/):
#   aligned_sorted.bam(.bai)   coordinate-sorted, all mapped reads (RETAINED)
#   flagstat.txt, idxstats.txt
# Plus one shared table: analysis/comparison/phase1_align_summary.tsv
#
# One library at a time (RAM + sort-temp disk). Run the big ones (bc03, bc07)
# inside tmux. Target box: c7i.8xlarge (32 vCPU / 64 GB) + 1 TB gp3 for analysis/.
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/data/EXP26000993}"
RUNDATA=$PROJECT_ROOT/rundata
ANALYSIS_ROOT=$PROJECT_ROOT/analysis
MMI_INDEX=$ANALYSIS_ROOT/refs/hg38_hbv_2x_splice.mmi
SAMPLES=$PROJECT_ROOT/config/samples.tsv
SUMMARY=$ANALYSIS_ROOT/comparison/phase1_align_summary.tsv
HBV_CONTIG=U95551.1_2x
PIPELINE_VERSION=phase1-v1
MAPQ=20

# Thread split: reserve a few cores for the sort in the pipe.
NPROC=$(nproc)
MM_THREADS=$(( NPROC > 6 ? NPROC - 4 : NPROC ))
SORT_THREADS=$(( NPROC > 6 ? 4 : 1 ))

mkdir -p "$ANALYSIS_ROOT/comparison" "$ANALYSIS_ROOT/logs"
command -v minimap2 >/dev/null || { echo "ERROR: minimap2 not on PATH (conda activate hbv_lr_analysis)"; exit 1; }
[ -s "$MMI_INDEX" ] || { echo "ERROR: splice index missing — run scripts/phase0_build_ref.sh first"; exit 1; }
[ -f "$SAMPLES" ]   || { echo "ERROR: sample map missing: $SAMPLES"; exit 1; }

# Optional filter list (sample_key or barcode) from CLI args.
WANT="$*"
want() { [ -z "$WANT" ] && return 0; for w in $WANT; do [ "$w" = "$1" ] || [ "$w" = "$2" ] && return 0; done; return 1; }

# (Re)create the summary header if starting a full run.
if [ -z "$WANT" ] || [ ! -f "$SUMMARY" ]; then
    printf '# experiment=EXP26000892_cDNA003; pipeline_version=%s; generated=%s\n' \
        "$PIPELINE_VERSION" "$(date +%F)" > "$SUMMARY"
    printf 'sample_key\tbarcode\tlib_id\tprimary_reads\tprimary_mapped\tmapping_rate\tsupplementary\thbv_alignments\thbv_reads_mapq%s\n' "$MAPQ" >> "$SUMMARY"
fi

# Iterate the sample map (skip comments / header).
awk -F'\t' 'NF>=4 && $1!~/^#/ && $1!="barcode"' "$SAMPLES" | while IFS=$'\t' read -r barcode lib_id sample_name sample_key rest; do
    want "$sample_key" "$barcode" || continue

    OUT_DIR=$ANALYSIS_ROOT/samples/$sample_key
    BAM=$OUT_DIR/aligned_sorted.bam
    LOG=$ANALYSIS_ROOT/logs/phase1_${sample_key}.log
    mkdir -p "$OUT_DIR"

    if [ -s "$BAM" ] && [ -s "$BAM.bai" ] && [ "${FORCE:-0}" != "1" ]; then
        echo "[$(date)] SKIP $sample_key (BAM exists; FORCE=1 to redo)"; continue
    fi

    shopt -s nullglob
    FQ=( "$RUNDATA/fastq_pass/$barcode"/*.fastq.gz )
    shopt -u nullglob
    if [ ${#FQ[@]} -eq 0 ]; then
        echo "[$(date)] WARN $sample_key ($barcode): no fastq_pass files — skipping"; continue
    fi

    echo "[$(date)] ALIGN $sample_key ($barcode): ${#FQ[@]} fastq files, mm=$MM_THREADS sort=$SORT_THREADS"
    minimap2 -ax splice --secondary=no -t "$MM_THREADS" "$MMI_INDEX" "${FQ[@]}" 2>"$LOG" \
        | samtools sort -@ "$SORT_THREADS" -m 2G -T "$OUT_DIR/sorttmp.$$" -o "$BAM" -
    samtools index "$BAM"
    samtools flagstat "$BAM" > "$OUT_DIR/flagstat.txt"
    samtools idxstats "$BAM" > "$OUT_DIR/idxstats.txt"

    PRIM=$(awk '/primary$/{print $1; exit}' "$OUT_DIR/flagstat.txt")
    PMAP=$(awk '/primary mapped/{print $1; exit}' "$OUT_DIR/flagstat.txt")
    SUPP=$(awk '/supplementary/{print $1; exit}' "$OUT_DIR/flagstat.txt")
    RATE=$(awk -v m="$PMAP" -v t="$PRIM" 'BEGIN{printf (t>0)?"%.4f":"NA", (t>0)?m/t:0}')
    HBV_ALN=$(awk -v c="$HBV_CONTIG" '$1==c{print $3}' "$OUT_DIR/idxstats.txt")
    HBV_Q=$(samtools view -c -q "$MAPQ" "$BAM" "$HBV_CONTIG")

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample_key" "$barcode" "$lib_id" "$PRIM" "$PMAP" "$RATE" "$SUPP" "${HBV_ALN:-0}" "$HBV_Q" >> "$SUMMARY"
    echo "[$(date)] DONE $sample_key — mapped $PMAP/$PRIM (rate $RATE); HBV aln ${HBV_ALN:-0}, HBV MAPQ≥$MAPQ $HBV_Q"
done

echo "[$(date)] Phase 1 complete. Summary: $SUMMARY"
column -t -s$'\t' "$SUMMARY" | grep -v '^#' || true
