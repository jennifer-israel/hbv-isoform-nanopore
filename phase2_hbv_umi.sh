#!/usr/bin/env bash
# Phase 2 (HBV UMI) — unique HBV molecules per library, the fast targeted way.
#
# HBV reads are rare (tens–hundreds/library), so instead of genome-wide pychopper we:
#   1. pull the HBV-aligning PRIMARY reads from the Phase 1 BAM, with FULL sequence
#      (soft-clipped ends retained — the UMI lives in the strand-switching primer that
#      minimap2 soft-clips off the genomic alignment);
#   2. run pychopper UMI extraction on just those reads (seconds; fixed -q avoids the
#      autotune-on-tiny-input crash);
#   3. re-align to a tiny HBV-only index (instant load) so umi_tools has coordinates;
#   4. umi_tools dedup the well-formed 28-nt UMIs → unique HBV molecules.
#
#   conda activate hbv_lr_analysis        # needs pychopper + umi_tools
#   bash scripts/phase2_hbv_umi.sh                     # all libraries
#   bash scripts/phase2_hbv_umi.sh <sample_key> ...    # subset
#
# Per library (analysis/samples/<sample_key>/): hbv.umi.bam(.bai) (RX-tagged, retained),
#   hbv_pychopper_stats.tsv, hbv_umi.log
# Shared: analysis/comparison/phase2_hbv_umi_summary.tsv
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/data/EXP26000896}"
ANALYSIS_ROOT=$PROJECT_ROOT/analysis
REF_2X=$ANALYSIS_ROOT/refs/hg38_hbv_2x.fa
HBV=U95551.1_2x
HBV_REF=$ANALYSIS_ROOT/refs/hbv_2x_only.fa
HBV_MMI=$ANALYSIS_ROOT/refs/hbv_2x_only.mmi
SAMPLES=$PROJECT_ROOT/config/samples.tsv
SUMMARY=$ANALYSIS_ROOT/comparison/phase2_hbv_umi_summary.tsv
Q=0.4
VER=phase2hbvumi-v1

for t in pychopper minimap2 samtools umi_tools; do command -v "$t" >/dev/null || { echo "ERROR: $t not on PATH"; exit 1; }; done

# Tiny HBV-only index for instant re-alignment (built once).
if [ ! -s "$HBV_MMI" ]; then
    echo "[$(date)] building HBV-only index"
    samtools faidx "$REF_2X" "$HBV" > "$HBV_REF"
    minimap2 -x splice -k15 -d "$HBV_MMI" "$HBV_REF" 2>/dev/null
fi

WANT="$*"
want(){ [ -z "$WANT" ] && return 0; for w in $WANT; do [ "$w" = "$1" ]||[ "$w" = "$2" ] && return 0; done; return 1; }

if [ -z "$WANT" ] || [ ! -f "$SUMMARY" ]; then
    printf '# experiment=EXP26000892_cDNA003; pipeline_version=%s; generated=%s\n' "$VER" "$(date +%F)" > "$SUMMARY"
    printf 'sample_key\tbarcode\thbv_primary_reads\thbv_umi_wellformed28\thbv_unique_molecules\thbv_dup_rate\n' >> "$SUMMARY"
fi

awk -F'\t' 'NF>=4 && $1!~/^#/ && $1!="barcode"' "$SAMPLES" | while IFS=$'\t' read -r barcode lib_id sample_name sample_key rest; do
    want "$sample_key" "$barcode" || continue
    D=$ANALYSIS_ROOT/samples/$sample_key
    BAM=$D/aligned_sorted.bam
    [ -s "$BAM" ] || { echo "WARN $sample_key: no aligned_sorted.bam"; continue; }

    NHBV=$(samtools view -c -F 0x900 "$BAM" "$HBV")
    if [ "$NHBV" -eq 0 ]; then
        printf '%s\t%s\t0\t0\t0\tNA\n' "$sample_key" "$barcode" >> "$SUMMARY"
        echo "[$(date)] $sample_key: 0 HBV reads"; continue
    fi

    set -o pipefail
    samtools view -b -F 0x900 "$BAM" "$HBV" | samtools fastq -n - 2>/dev/null \
      | pychopper -k PCB114 -m edlib -U -y -q "$Q" -t 4 -S "$D/hbv_pychopper_stats.tsv" - - 2>"$D/hbv_umi.log" \
      | minimap2 -ax splice --secondary=no -y -t 4 "$HBV_MMI" - 2>>"$D/hbv_umi.log" \
      | samtools sort -o "$D/hbv.umi.bam" -
    samtools index "$D/hbv.umi.bam"

    samtools view -b -e 'length([RX])==28' "$D/hbv.umi.bam" | samtools sort -o "$D/hbv.umi28.bam" - && samtools index "$D/hbv.umi28.bam"
    WF=$(samtools view -c -F 0x904 "$D/hbv.umi28.bam")
    umi_tools dedup -I "$D/hbv.umi28.bam" -S "$D/hbv.umi.dedup.bam" \
        --extract-umi-method=tag --umi-tag=RX --method=directional -L "$D/hbv_umi_dedup.log" >/dev/null 2>&1 || true
    UNIQ=$(samtools view -c -F 0x904 "$D/hbv.umi.dedup.bam" 2>/dev/null || echo 0)
    DUP=$(awk -v u="$UNIQ" -v w="$WF" 'BEGIN{printf (w>0)?"%.4f":"NA",(w>0)?1-u/w:0}')

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$sample_key" "$barcode" "$NHBV" "$WF" "$UNIQ" "$DUP" >> "$SUMMARY"
    echo "[$(date)] $sample_key: HBV reads $NHBV → wellformed-UMI $WF → unique $UNIQ (dup $DUP)"
    # persist the surviving (one-per-molecule) read IDs before discarding the dedup BAM,
    # so downstream transcript figures can restrict to de-duplicated molecules (see phase2b).
    samtools view -F 0x904 "$D/hbv.umi.dedup.bam" 2>/dev/null | cut -f1 | sort -u > "$D/hbv.dedup_read_ids.txt" || : > "$D/hbv.dedup_read_ids.txt"
    rm -f "$D/hbv.umi28.bam"* "$D/hbv.umi.dedup.bam"*
done

echo "[$(date)] Phase 2 HBV-UMI complete."
column -t -s$'\t' "$SUMMARY" | grep -v '^#'