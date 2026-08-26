#!/usr/bin/env bash
# UMI saturation curve — are the unique-molecule counts complete, and is the reported
# duplication rate real or an artefact of over-sequencing?  (EXP26000892)
#
#   bash umi_saturation_curve.sh                          # all libraries
#   bash umi_saturation_curve.sh SeqLib5551_10ng_polyA_17  # one
#
# WHY
# Phase 2 reported 91.7-99.7% duplication. Those are depth/complexity ratios, so they
# are only meaningful if the molecule counts are COMPLETE. Two failure modes:
#
#   (a) undersampling — if more reads would still reveal more molecules, the unique
#       count is a floor and duplication is overstated.
#   (b) umi_tools over-merging — --method=directional collapses UMIs within 1 edit,
#       and in pileups of thousands that can merge genuinely distinct molecules.
#
# WHAT IT DOES
# Subsamples each library's UMI-tagged HBV reads to 10/25/50/75/100% and re-runs dedup
# at each depth, with BOTH methods:
#   directional  the production setting (error-tolerant, can over-merge)
#   unique       exact UMI matching, no merging — an upper bound
#
# READING THE OUTPUT
#   curve PLATEAUS      -> molecules are essentially all recovered; the counts and the
#                          duplication rates stand, and you are over-sequencing.
#   curve STILL CLIMBING -> undercounting; more depth would find more molecules and the
#                          quoted duplication is too high.
#   directional ~ unique -> merging is modest; directional is safe.
#   unique >> directional -> directional is collapsing distinct molecules; revisit
#                          --edit-distance-threshold before trusting the counts.
#
# Also reports molecules-per-1000-reads at each depth, which is the number to use when
# comparing libraries sequenced to different depths (raw duplication rate is NOT
# comparable across libraries — that is what this normalises).
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/data/EXP26000993}"
ANALYSIS=$PROJECT_ROOT/analysis
SAMPLES=$PROJECT_ROOT/config/samples.tsv
OUT=$ANALYSIS/comparison/umi_saturation.tsv
FRACTIONS=${FRACTIONS:-"0.10 0.25 0.50 0.75 1.00"}
SEED=${SEED:-42}
THREADS=${THREADS:-4}
WORK=${WORK:-/tmp/umisat.$$}

for t in samtools umi_tools; do command -v $t >/dev/null || { echo "ERROR: $t not on PATH" >&2; exit 1; }; done
mkdir -p "$WORK" "$ANALYSIS/comparison"
trap 'rm -rf "$WORK"' EXIT

printf '# experiment=%s; generated=%s; seed=%s\n' "$(basename "$PROJECT_ROOT")" "$(date +%F)" "$SEED" > "$OUT"
printf 'sample_key\tfraction\treads\tunique_directional\tunique_exact\tdup_rate_directional\tmolecules_per_1k_reads\tmerge_ratio\n' >> "$OUT"

WANT="$*"
want(){ [ -z "$WANT" ] && return 0; for w in $WANT; do [ "$w" = "$1" ] && return 0; done; return 1; }

echo "=== UMI saturation curve ==="
echo "  fractions: $FRACTIONS"
echo

awk -F'\t' 'NF>=4 && $1!~/^#/ && $1!="barcode"' "$SAMPLES" | while IFS=$'\t' read -r barcode lib_id sample_name sample_key rest; do
    want "$sample_key" || continue
    D=$ANALYSIS/samples/$sample_key
    [ -s "$D/hbv.umi.bam" ] || { echo "SKIP $sample_key: no hbv.umi.bam"; continue; }

    # Well-formed 28-nt UMIs only — same population Phase 2 deduplicated.
    BASE=$WORK/$sample_key.28.bam
    samtools view -b -e 'length([RX])==28' "$D/hbv.umi.bam" \
      | samtools sort -@ "$THREADS" -o "$BASE" -
    samtools index "$BASE"
    TOTAL=$(samtools view -c -F 0x904 "$BASE")
    echo "--- $sample_key: $TOTAL well-formed reads"
    printf '    %-8s %-10s %-12s %-12s %-10s %-12s %s\n' \
           frac reads directional exact dup_rate mol/1k merge
    echo "    ---------------------------------------------------------------------------"

    for FRAC in $FRACTIONS; do
        if [ "$FRAC" = "1.00" ]; then
            SUB=$BASE
        else
            SUB=$WORK/sub.bam
            samtools view -b -s "${SEED}${FRAC#0}" "$BASE" > "$SUB"
            samtools index "$SUB"
        fi
        N=$(samtools view -c -F 0x904 "$SUB")
        [ "$N" -eq 0 ] && continue

        umi_tools dedup -I "$SUB" -S "$WORK/dd.bam" --extract-umi-method=tag \
            --umi-tag=RX --method=directional -L "$WORK/dd.log" >/dev/null 2>&1 || true
        UD=$(samtools view -c -F 0x904 "$WORK/dd.bam" 2>/dev/null || echo 0)

        umi_tools dedup -I "$SUB" -S "$WORK/du.bam" --extract-umi-method=tag \
            --umi-tag=RX --method=unique -L "$WORK/du.log" >/dev/null 2>&1 || true
        UU=$(samtools view -c -F 0x904 "$WORK/du.bam" 2>/dev/null || echo 0)

        DUP=$(awk -v u="$UD" -v n="$N" 'BEGIN{printf (n>0)?"%.4f":"NA",(n>0)?1-u/n:0}')
        PERK=$(awk -v u="$UD" -v n="$N" 'BEGIN{printf (n>0)?"%.1f":"NA",(n>0)?1000*u/n:0}')
        MERGE=$(awk -v d="$UD" -v e="$UU" 'BEGIN{printf (d>0)?"%.2f":"NA",(d>0)?e/d:0}')

        printf '    %-8s %-10s %-12s %-12s %-10s %-12s %s\n' "$FRAC" "$N" "$UD" "$UU" "$DUP" "$PERK" "$MERGE"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$sample_key" "$FRAC" "$N" "$UD" "$UU" "$DUP" "$PERK" "$MERGE" >> "$OUT"
        rm -f "$WORK/dd.bam"* "$WORK/du.bam"* "$WORK/sub.bam"*
    done
    echo
    rm -f "$BASE"*
done

echo "=== written: $OUT ==="
echo
echo "Interpretation:"
echo "  molecules_per_1k_reads FALLING steeply as fraction rises  -> saturating (good)"
echo "  molecules_per_1k_reads roughly FLAT                       -> still undersampled;"
echo "                                                              unique counts are a floor"
echo "  merge_ratio (exact/directional) near 1.0                  -> little merging"
echo "  merge_ratio > ~1.3                                        -> directional collapsing"
echo "                                                              distinct molecules"
