#!/usr/bin/env bash
# Phase 1.5B — UMI-based PCR-duplication & library-complexity QC (per library).
#
# Flow (UMIs must be read while primers are on the read → pychopper BEFORE minimap2):
#   pychopper -k PCB114 -U -y  → oriented/trimmed full-length reads, UMI in comment
#   minimap2 -y                → UMI carried into BAM as RX:Z:
#   umi_tools dedup            → unique molecules (error-tolerant, directional)
#
#   conda activate hbv_lr_analysis      # needs pychopper (>=2.7.10) + umi_tools
#   bash scripts/phase1_5_umi.sh                 # all libraries
#   bash scripts/phase1_5_umi.sh <sample_key|barcode> ...   # subset
#   KEEP_UMI_BAM=1 bash scripts/phase1_5_umi.sh <key>       # keep full UMI BAM
#
# Per library (analysis/samples/<sample_key>/):
#   pychopper_stats.tsv, pychopper_report.pdf   full-length / rescue / UMI rates
#   umi_tagged.bam(.bai)                        RX:Z: tagged (deleted unless KEEP_UMI_BAM=1)
#   hbv.umi.bam(.bai)                           HBV-region UMI-tagged reads → Phase 2 (RETAINED)
# Shared: analysis/comparison/phase1_5_umi_summary.tsv
#   dup rate = 1 - deduped/mapped ; unique molecules = deduped mapped reads
set -euo pipefail

PROJECT_ROOT=/home/ubuntu/matt_wolpert_claude_code_analysis/2026_07_02_EXP26000559_cDNA001
RUNDATA=$PROJECT_ROOT/rundata
ANALYSIS_ROOT=$PROJECT_ROOT/analysis
MMI_INDEX=$ANALYSIS_ROOT/refs/hg38_hbv_2x_splice.mmi
SAMPLES=$PROJECT_ROOT/config/samples.tsv
SUMMARY=$ANALYSIS_ROOT/comparison/phase1_5_umi_summary.tsv
HBV_CONTIG=U95551.1_2x
PIPELINE_VERSION=phase1_5umi-v1
NPROC=$(nproc)
# env-overridable so the run can coexist with other jobs on a shared box
PY_THREADS=${PY_THREADS:-$(( NPROC>8 ? NPROC-8 : NPROC ))}   # pychopper is the pipe bottleneck → most threads
MM_THREADS=${MM_THREADS:-$(( NPROC>8 ? 6 : NPROC ))}
SORT_THREADS=${SORT_THREADS:-2}

mkdir -p "$ANALYSIS_ROOT/comparison" "$ANALYSIS_ROOT/logs"
for t in pychopper minimap2 samtools umi_tools; do
  command -v "$t" >/dev/null || { echo "ERROR: '$t' not on PATH. conda install -c bioconda pychopper umi_tools"; exit 1; }
done
[ -s "$MMI_INDEX" ] || { echo "ERROR: run scripts/phase0_build_ref.sh first"; exit 1; }

WANT="$*"
want(){ [ -z "$WANT" ] && return 0; for w in $WANT; do [ "$w" = "$1" ]||[ "$w" = "$2" ] && return 0; done; return 1; }

if [ -z "$WANT" ] || [ ! -f "$SUMMARY" ]; then
  printf '# experiment=EXP26000559_cDNA001; pipeline_version=%s; generated=%s\n' "$PIPELINE_VERSION" "$(date +%F)" > "$SUMMARY"
  printf 'sample_key\tbarcode\tinput_reads\tfull_length_reads\tumi_wellformed_28nt\tunique_molecules\tdup_rate\n' >> "$SUMMARY"
fi

awk -F'\t' 'NF>=4 && $1!~/^#/ && $1!="barcode"' "$SAMPLES" | while IFS=$'\t' read -r barcode lib_id sample_name sample_key rest; do
  want "$sample_key" "$barcode" || continue
  OUT_DIR=$ANALYSIS_ROOT/samples/$sample_key; mkdir -p "$OUT_DIR"
  LOG=$ANALYSIS_ROOT/logs/phase1_5umi_${sample_key}.log
  UBAM=$OUT_DIR/umi_tagged.bam
  HBVUMI=$OUT_DIR/hbv.umi.bam
  COMBINED=$OUT_DIR/combined.fastq.gz

  shopt -s nullglob; FQ=( "$RUNDATA/fastq_pass/$barcode"/*.fastq.gz ); shopt -u nullglob
  [ ${#FQ[@]} -eq 0 ] && { echo "[$(date)] WARN $sample_key: no fastq — skip"; continue; }

  # pychopper's autotune reads its input TWICE (count + sample), so it needs a SEEKABLE
  # file, not a pipe. Concatenate the per-barcode gzips into one (valid multi-member
  # gzip), let pychopper autotune + classify from it, and stream its OUTPUT to minimap2
  # (no uncompressed FASTQ on disk). combined.fastq.gz is transient, deleted right after.
  # SUBSAMPLE=N (env) → random ~N-read subset for duplication-RATE QC on huge libraries
  # (dup rate is depth-dependent, so we normalize the big NOpolyA libs to a common depth).
  # Probabilistic keep at p=N/total (total from Phase 1 summary); no seqtk/seqkit needed.
  if [ "${SUBSAMPLE:-0}" -gt 0 ]; then
    TOT=$(awk -F'\t' -v k="$sample_key" '$1==k{print $4}' "$ANALYSIS_ROOT/comparison/phase1_align_summary.tsv" 2>/dev/null)
    if [ -n "$TOT" ] && [ "$TOT" -gt "$SUBSAMPLE" ]; then
      P=$(awk -v n="$SUBSAMPLE" -v t="$TOT" 'BEGIN{printf "%.6f", n/t}')
      echo "[$(date)] CONCAT+SUBSAMPLE $sample_key: ~$SUBSAMPLE of $TOT reads (p=$P)"
      zcat "${FQ[@]}" | awk -v p="$P" 'BEGIN{srand(42)} (NR-1)%4==0{keep=(rand()<=p)} keep' | gzip -1 > "$COMBINED"
    else
      echo "[$(date)] CONCAT $sample_key: full (total ${TOT:-?} <= target $SUBSAMPLE)"
      cat "${FQ[@]}" > "$COMBINED"
    fi
  else
    echo "[$(date)] CONCAT $sample_key ($barcode): ${#FQ[@]} files → combined.fastq.gz"
    cat "${FQ[@]}" > "$COMBINED"
  fi

  echo "[$(date)] PYCHOPPER+ALIGN+TAG $sample_key; py=$PY_THREADS mm=$MM_THREADS"
  set -o pipefail
  # NB: pychopper 2.7.10's _plot_stats crashes under pandas 3.x (float(Series)); we
  # patched the installed pychopper.py line 184 (backup at pychopper.py.orig).
  pychopper -k PCB114 -m edlib -U -y -t "$PY_THREADS" \
        -S "$OUT_DIR/pychopper_stats.tsv" -r "$OUT_DIR/pychopper_report.pdf" "$COMBINED" - 2>"$LOG" \
    | minimap2 -ax splice --secondary=no -y -t "$MM_THREADS" "$MMI_INDEX" - 2>>"$LOG" \
    | samtools sort -@ "$SORT_THREADS" -m 2G -T "$OUT_DIR/umisort.$$" -o "$UBAM" -
  samtools index "$UBAM"
  rm -f "$COMBINED"

  # HBV-region UMI reads → tiny BAM handed to Phase 2 (retained even if full UMI BAM deleted)
  samtools view -b "$UBAM" "$HBV_CONTIG" | samtools sort -o "$HBVUMI" - && samtools index "$HBVUMI"

  # umi_tools requires equal-length UMIs; ONT indels make pychopper UMIs vary (mostly
  # 28 nt = full structured probe, plus 27/29/… and RX:Z:None for undetected). Dedup the
  # well-formed 28-nt UMIs only — padding indel-shifted UMIs would frameshift + corrupt
  # clustering. ~93% of tagged reads are the clean 28-nt length.
  echo "[$(date)] DEDUP $sample_key (umi_tools directional; 28-nt UMIs only)"
  UMI28=$OUT_DIR/umi28.bam
  samtools view -b -e 'length([RX])==28' "$UBAM" | samtools sort -@ "$SORT_THREADS" -o "$UMI28" - && samtools index "$UMI28"
  umi_tools dedup -I "$UMI28" -S "$OUT_DIR/umi_dedup.bam" \
      --extract-umi-method=tag --umi-tag=RX --method=directional \
      -L "$OUT_DIR/umi_tools_dedup.log" >/dev/null 2>&1 || true

  INREADS=$(awk -F'\t' -v k="$sample_key" '$1==k{print $4}' "$ANALYSIS_ROOT/comparison/phase1_align_summary.tsv" 2>/dev/null); [ -z "$INREADS" ] && INREADS=NA
  FLREADS=$(samtools view -c -F 0x900 "$UBAM")           # primary reads = pychopper full-length output
  WELLFORMED=$(samtools view -c -F 0x904 "$UMI28")       # primary mapped w/ 28-nt UMI (dedup denominator)
  UNIQ=$(samtools view -c -F 0x904 "$OUT_DIR/umi_dedup.bam" 2>/dev/null || echo 0)
  DUP=$(awk -v u="$UNIQ" -v w="$WELLFORMED" 'BEGIN{printf (w>0)?"%.4f":"NA",(w>0)?1-u/w:0}')

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$sample_key" "$barcode" "$INREADS" "$FLREADS" "$WELLFORMED" "$UNIQ" "$DUP" >> "$SUMMARY"
  echo "[$(date)] DONE $sample_key — full-length $FLREADS/$INREADS; unique $UNIQ/$WELLFORMED (dup $DUP)"

  rm -f "$OUT_DIR/umi_dedup.bam" "$UMI28" "$UMI28.bai"
  [ "${KEEP_UMI_BAM:-0}" = "1" ] || { rm -f "$UBAM" "$UBAM.bai"; echo "  (removed full umi_tagged.bam; hbv.umi.bam kept)"; }
done

echo "[$(date)] Phase 1.5B complete. Summary: $SUMMARY"
column -t -s$'\t' "$SUMMARY" | grep -v '^#' || true
