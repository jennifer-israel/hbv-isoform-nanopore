#!/usr/bin/env bash
# Phase 2b — persist the UMI-dedup survivor read IDs per library.
#
# phase2_hbv_umi.sh counts unique HBV molecules by umi_tools dedup, but DELETES the
# deduplicated BAM afterwards, keeping only the count. Downstream transcript figures
# (composition / TSS / read-span, report figs 4–6) therefore had no way to restrict to
# the de-duplicated molecules and were drawn over ALL HBV reads (raw, with PCR repeats).
#
# This step reproduces that exact dedup — well-formed 28-nt UMIs, umi_tools
# --method=directional on the RX tag — from the RETAINED hbv.umi.bam (untouched), and
# writes the surviving (one-per-molecule) read IDs to:
#     analysis/samples/<sample_key>/hbv.dedup_read_ids.txt
# The count in this file equals hbv_unique_molecules in phase2_hbv_umi_summary.tsv.
#
#   conda activate hbv_lr_analysis        # needs umi_tools + samtools
#   bash scripts/phase2b_dedup_read_ids.sh              # all libraries
#   bash scripts/phase2b_dedup_read_ids.sh <sample_key> ...   # subset
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/data/EXP26000896}"
ANALYSIS_ROOT=$PROJECT_ROOT/analysis

for t in samtools umi_tools; do command -v "$t" >/dev/null || { echo "ERROR: $t not on PATH"; exit 1; }; done

WANT="$*"
want(){ [ -z "$WANT" ] && return 0; for w in $WANT; do [ "$w" = "$1" ] && return 0; done; return 1; }

for D in "$ANALYSIS_ROOT"/samples/*/; do
    sk=$(basename "$D")
    want "$sk" || continue
    BAM=$D/hbv.umi.bam
    OUT=$D/hbv.dedup_read_ids.txt
    if [ ! -s "$BAM" ]; then
        : > "$OUT"; echo "$sk: no hbv.umi.bam → 0 dedup ids"; continue
    fi
    tmp=$(mktemp -d)
    # same filter + dedup as phase2_hbv_umi.sh (well-formed 28-nt UMIs, directional)
    samtools view -b -e 'length([RX])==28' "$BAM" | samtools sort -o "$tmp/umi28.bam" - && samtools index "$tmp/umi28.bam"
    if [ "$(samtools view -c -F 0x904 "$tmp/umi28.bam")" -eq 0 ]; then
        : > "$OUT"; rm -rf "$tmp"; echo "$sk: 0 well-formed UMIs → 0 dedup ids"; continue
    fi
    umi_tools dedup -I "$tmp/umi28.bam" -S "$tmp/dedup.bam" \
        --extract-umi-method=tag --umi-tag=RX --method=directional -L /dev/null >/dev/null 2>&1
    samtools view -F 0x904 "$tmp/dedup.bam" | cut -f1 | sort -u > "$OUT"
    rm -rf "$tmp"
    echo "$sk: $(wc -l < "$OUT") unique-molecule read ids → $OUT"
done
