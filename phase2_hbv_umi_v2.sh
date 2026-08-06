#!/usr/bin/env bash
# Phase 2 (HBV UMI) v2 — unique HBV molecules per library.
#
# Same logic, inputs and OUTPUT FORMAT as phase2_hbv_umi.sh (phase4_quantify.py reads
# phase2_hbv_umi_summary.tsv unchanged), with three changes forced by this experiment:
#
#   1. THREADS. The original hardcodes -t 4 throughout because EXP26000559 had
#      tens-to-hundreds of HBV reads per library. EXP26000892 has 439k-1.29M after
#      hybridisation capture. Threads now come from nproc.
#
#   2. TIMING + PROGRESS. Each stage is timed and announced, so a stall is visible
#      rather than being indistinguishable from slow progress.
#
#   3. DEDUP GUARDRAILS. umi_tools dedup --method=directional builds a similarity
#      network per alignment position; on a 6,364 bp contig carrying >1M reads the
#      pileups are enormous and this is roughly quadratic within each. Set
#      UMI_METHOD=unique for a fast (error-intolerant) lower bound first if
#      directional is not converging. --output-stats is deliberately NOT used; it is
#      very expensive at this depth.
#
#   conda activate hbv_lr_analysis
#   bash phase2_hbv_umi_v2.sh                                  # all libraries
#   bash phase2_hbv_umi_v2.sh SeqLib5551_10ng_polyA_17         # one (start here)
#   UMI_METHOD=unique bash phase2_hbv_umi_v2.sh <key>          # fast lower bound
#   FORCE=1 bash phase2_hbv_umi_v2.sh <key>                    # redo despite existing bam
set -euo pipefail

PROJECT_ROOT=${PROJECT_ROOT:-/data/EXP26000896}
ANALYSIS_ROOT=$PROJECT_ROOT/analysis
REF_2X=$ANALYSIS_ROOT/refs/hg38_hbv_2x.fa
HBV=U95551.1_2x
HBV_REF=$ANALYSIS_ROOT/refs/hbv_2x_only.fa
HBV_MMI=$ANALYSIS_ROOT/refs/hbv_2x_only.mmi
SAMPLES=$PROJECT_ROOT/config/samples.tsv
SUMMARY=$ANALYSIS_ROOT/comparison/phase2_hbv_umi_summary.tsv
Q=0.4
VER=phase2hbvumi-v3
UMI_METHOD=${UMI_METHOD:-directional}

# --- homologous concatemer filter -------------------------------------------
# dorado --barcode-both-ends rejects fusions of two DIFFERENT libraries (the ends
# disagree). A fusion of two molecules from the SAME library agrees at both ends and is
# assigned normally — nothing downstream removes it. On the 2x reference such a read can
# align contiguously across the copy junction and be classified pgRNA_RT (span >= 3,982 bp,
# "tandem/concatemeric readthrough"), i.e. a PCR artifact scored as biology.
#
# The filter must run BEFORE pychopper, which trims the primers away: after trimming there
# is no adapter left to detect. Set CONCAT_FILTER=0 to disable.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
FILTER_PY=${FILTER_PY:-$SCRIPT_DIR/filter_concatemers.py}
BARCODES_FA=${BARCODES_FA:-$SCRIPT_DIR/barcode_sequences_EXP26000892.fasta}
END_MARGIN=${END_MARGIN:-200}
CONCAT_FILTER=${CONCAT_FILTER:-1}

NPROC=$(nproc)
PY_THREADS=${PY_THREADS:-$(( NPROC > 8 ? NPROC - 6 : NPROC ))}   # pychopper is the bottleneck
MM_THREADS=${MM_THREADS:-$(( NPROC > 8 ? 6 : NPROC ))}
SORT_THREADS=${SORT_THREADS:-4}

for t in pychopper minimap2 samtools umi_tools; do
    command -v "$t" >/dev/null || { echo "ERROR: $t not on PATH" >&2; exit 1; }
done

# samtools -e filter expressions need >= 1.16
SV=$(samtools --version | head -1 | awk '{print $2}')
awk -v v="$SV" 'BEGIN{split(v,a,"."); if (a[1]<1 || (a[1]==1 && a[2]<16)) exit 1}' \
    || { echo "ERROR: samtools $SV too old; 'view -e' needs >= 1.16" >&2; exit 1; }

if [ "$CONCAT_FILTER" = "1" ]; then
    [ -s "$FILTER_PY" ]   || { echo "ERROR: filter not found: $FILTER_PY (CONCAT_FILTER=0 to skip)" >&2; exit 1; }
    [ -s "$BARCODES_FA" ] || { echo "ERROR: barcode FASTA not found: $BARCODES_FA" >&2; exit 1; }
    python3 -c "import edlib" 2>/dev/null || { echo "ERROR: filter needs edlib (pip install edlib)" >&2; exit 1; }
fi

echo "=== Phase 2 HBV-UMI v3 ==="
echo "  nproc=$NPROC  pychopper=$PY_THREADS  minimap2=$MM_THREADS  sort=$SORT_THREADS"
echo "  samtools=$SV  umi_tools method=$UMI_METHOD"
if [ "$CONCAT_FILTER" = "1" ]; then
    echo "  concatemer filter: ON  (end_margin=$END_MARGIN, barcodes=$(basename "$BARCODES_FA"))"
else
    echo "  concatemer filter: OFF"
fi
echo "  free disk: $(df -h "$ANALYSIS_ROOT" | awk 'NR==2{print $4}')"
echo

# Tiny HBV-only index (built once).
if [ ! -s "$HBV_MMI" ]; then
    echo "[$(date +%T)] building HBV-only index"
    samtools faidx "$REF_2X" "$HBV" > "$HBV_REF"
    minimap2 -x splice -k15 -d "$HBV_MMI" "$HBV_REF" 2>/dev/null
fi

WANT="$*"
want(){ [ -z "$WANT" ] && return 0; for w in $WANT; do [ "$w" = "$1" ]||[ "$w" = "$2" ] && return 0; done; return 1; }

if [ -z "$WANT" ] || [ ! -f "$SUMMARY" ]; then
    printf '# experiment=EXP26000892_cDNA003; pipeline_version=%s; generated=%s\n' "$VER" "$(date +%F)" > "$SUMMARY"
    printf 'sample_key\tbarcode\thbv_primary_reads\thbv_concatemers_removed\thbv_concatemer_pct\thbv_umi_wellformed28\thbv_unique_molecules\thbv_dup_rate\n' >> "$SUMMARY"
fi

hms(){ printf '%dm%02ds' $(( $1 / 60 )) $(( $1 % 60 )); }

awk -F'\t' 'NF>=4 && $1!~/^#/ && $1!="barcode"' "$SAMPLES" | while IFS=$'\t' read -r barcode lib_id sample_name sample_key rest; do
    want "$sample_key" "$barcode" || continue
    D=$ANALYSIS_ROOT/samples/$sample_key
    BAM=$D/aligned_sorted.bam
    [ -s "$BAM" ] || { echo "WARN $sample_key: no aligned_sorted.bam — run Phase 1 first"; continue; }

    if [ -s "$D/hbv.umi.bam" ] && [ "${FORCE:-0}" != "1" ]; then
        echo "[$(date +%T)] SKIP $sample_key (hbv.umi.bam exists; FORCE=1 to redo)"; continue
    fi

    T0=$SECONDS
    NHBV=$(samtools view -c -F 0x900 "$BAM" "$HBV")
    echo "[$(date +%T)] === $sample_key ($barcode): $NHBV HBV primary reads ==="
    if [ "$NHBV" -eq 0 ]; then
        printf '%s\t%s\t0\t0\tNA\t0\t0\tNA\n' "$sample_key" "$barcode" >> "$SUMMARY"
        echo "   0 HBV reads — skipping"; continue
    fi

    # ---- stage 1: extract -> [concatemer filter] -> pychopper (UMI) -> realign -> sort ----
    echo "[$(date +%T)]   stage 1/3 filter + pychopper + realign  (this is the long one)"
    S1=$SECONDS
    FSTATS=$D/hbv_concatemer_filter.tsv
    rm -f "$FSTATS"
    if [ "$CONCAT_FILTER" = "1" ]; then
        FILT=(python3 "$FILTER_PY" --barcodes-fasta "$BARCODES_FA"
              --end-margin "$END_MARGIN" --stats "$FSTATS" --label "$sample_key")
    else
        FILT=(cat)
    fi
    set -o pipefail
    samtools view -b -F 0x900 "$BAM" "$HBV" \
      | samtools fastq -@ "$SORT_THREADS" -n - 2>/dev/null \
      | "${FILT[@]}" 2>>"$D/hbv_umi.log" \
      | pychopper -k PCB114 -m edlib -U -y -q "$Q" -t "$PY_THREADS" \
            -S "$D/hbv_pychopper_stats.tsv" - - 2>>"$D/hbv_umi.log" \
      | minimap2 -ax splice --secondary=no -y -t "$MM_THREADS" "$HBV_MMI" - 2>>"$D/hbv_umi.log" \
      | samtools sort -@ "$SORT_THREADS" -o "$D/hbv.umi.bam" -
    samtools index "$D/hbv.umi.bam"

    CONCAT_N=NA; CONCAT_PCT=NA
    if [ -s "$FSTATS" ]; then
        CONCAT_N=$(awk -F'\t' 'NR==2{print $4}' "$FSTATS")
        CONCAT_PCT=$(awk -F'\t' 'NR==2{printf "%.4f",$5}' "$FSTATS")
        echo "[$(date +%T)]   concatemers removed: $CONCAT_N of $NHBV (${CONCAT_PCT}%)"
    fi
    echo "[$(date +%T)]   stage 1 done in $(hms $(( SECONDS - S1 )))"

    # ---- stage 2: keep well-formed 28-nt UMIs ----
    echo "[$(date +%T)]   stage 2/3 filtering to 28-nt UMIs"
    S2=$SECONDS
    samtools view -b -e 'length([RX])==28' "$D/hbv.umi.bam" \
      | samtools sort -@ "$SORT_THREADS" -o "$D/hbv.umi28.bam" -
    samtools index "$D/hbv.umi28.bam"
    WF=$(samtools view -c -F 0x904 "$D/hbv.umi28.bam")
    TAGGED=$(samtools view -c -F 0x904 "$D/hbv.umi.bam")
    echo "[$(date +%T)]   stage 2 done in $(hms $(( SECONDS - S2 )))  well-formed $WF / $TAGGED mapped"

    # ---- stage 3: dedup ----
    echo "[$(date +%T)]   stage 3/3 umi_tools dedup (--method=$UMI_METHOD) on $WF reads"
    if [ "$WF" -gt 500000 ] && [ "$UMI_METHOD" = "directional" ]; then
        echo "           NOTE: >500k reads on a 6.4 kb contig. Directional clustering is"
        echo "           ~quadratic within each position pileup and may take hours. If it"
        echo "           stalls, rerun this library with UMI_METHOD=unique for a lower bound."
    fi
    S3=$SECONDS
    umi_tools dedup -I "$D/hbv.umi28.bam" -S "$D/hbv.umi.dedup.bam" \
        --extract-umi-method=tag --umi-tag=RX --method="$UMI_METHOD" \
        -L "$D/hbv_umi_dedup.log" >/dev/null 2>&1 || \
        echo "           WARN: umi_tools returned non-zero; see $D/hbv_umi_dedup.log"
    UNIQ=$(samtools view -c -F 0x904 "$D/hbv.umi.dedup.bam" 2>/dev/null || echo 0)
    echo "[$(date +%T)]   stage 3 done in $(hms $(( SECONDS - S3 )))"

    DUP=$(awk -v u="$UNIQ" -v w="$WF" 'BEGIN{printf (w>0)?"%.4f":"NA",(w>0)?1-u/w:0}')
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample_key" "$barcode" "$NHBV" "$CONCAT_N" "$CONCAT_PCT" "$WF" "$UNIQ" "$DUP" >> "$SUMMARY"
    echo "[$(date +%T)] DONE $sample_key in $(hms $(( SECONDS - T0 )))"
    echo "           HBV reads $NHBV -> concatemers removed $CONCAT_N -> 28nt-UMI $WF -> unique $UNIQ (dup $DUP)"
    echo

    rm -f "$D/hbv.umi28.bam"* "$D/hbv.umi.dedup.bam"*
done

echo "=== Phase 2 complete ==="
column -t -s$'\t' "$SUMMARY" | grep -v '^#'
