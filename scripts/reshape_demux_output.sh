#!/usr/bin/env bash
# Standalone reshape: dorado 2.x demux output -> the layout phase1_align.sh globs.
# Run this AFTER a successful `dorado demux`; it does NOT re-run the demux.
#
#   bash reshape_demux_output.sh /data/EXP26000896/demux/strict /data/EXP26000896/demux/pipeline_input
#
# dorado 2.x mirrors MinKNOW's output structure and names barcode dirs generically:
#   <demux_out>/no_sample/<run_id>/bam_pass/barcode01/*.bam
#                                           barcode02/*.bam ...
#                                           unclassified/*.bam
# The barcodeNN index corresponds to the arrangement's first_index..last_index, i.e.
# barcode01 -> MW01 -> 100 ng (SeqLib5550), barcode02 -> MW02 -> 10 ng, etc.
#
# Emits:  <out>/fastq_pass/custom_bc01/custom_bc01.fastq.gz   (etc.)
# so config/samples_EXP26000892.tsv needs no change.
set -euo pipefail

SRC=${1:?usage: reshape_demux_output.sh <demux_pass_dir> <pipeline_input_dir>}
DEST=${2:?usage: reshape_demux_output.sh <demux_pass_dir> <pipeline_input_dir>}
GZIP_LEVEL=${GZIP_LEVEL:-1}
THREADS=${THREADS:-4}

command -v samtools >/dev/null || { echo "ERROR: samtools not on PATH" >&2; exit 1; }

PIPE_DIR=$DEST/fastq_pass
mkdir -p "$PIPE_DIR"

# Find every barcode directory anywhere under SRC (handles the nested MinKNOW layout).
mapfile -t BCDIRS < <(find "$SRC" -type d \( -name 'barcode*' -o -name 'unclassified' \) | sort)
if [ ${#BCDIRS[@]} -eq 0 ]; then
    echo "ERROR: no barcode*/unclassified directories found under $SRC" >&2; exit 1
fi

echo "=== reshaping $SRC -> $PIPE_DIR ==="
printf '  found %d barcode directories\n\n' "${#BCDIRS[@]}"

TOTAL=0
for d in "${BCDIRS[@]}"; do
    name=$(basename "$d")
    case "$name" in
        unclassified) out=unclassified ;;
        barcode*)     out="custom_bc${name#barcode}" ;;   # barcode01 -> custom_bc01
        *)            echo "  skip $name"; continue ;;
    esac

    shopt -s nullglob
    bams=( "$d"/*.bam )
    fqs=(  "$d"/*.fastq.gz "$d"/*.fastq )
    shopt -u nullglob

    mkdir -p "$PIPE_DIR/$out"
    target=$PIPE_DIR/$out/$out.fastq.gz

    if [ ${#bams[@]} -gt 0 ]; then
        echo "  $name -> $out  (${#bams[@]} bam files)"
        # samtools cat concatenates BAMs cheaply, then one fastq conversion.
        samtools cat "${bams[@]}" \
          | samtools fastq -@ "$THREADS" -n - 2>/dev/null \
          | gzip -"$GZIP_LEVEL" > "$target"
    elif [ ${#fqs[@]} -gt 0 ]; then
        echo "  $name -> $out  (${#fqs[@]} fastq files)"
        : > "$target"
        for f in "${fqs[@]}"; do
            case "$f" in
                *.gz) cat "$f" >> "$target" ;;                      # valid multi-member gzip
                *)    gzip -"$GZIP_LEVEL" -c "$f" >> "$target" ;;
            esac
        done
    else
        echo "  $name -> $out  (EMPTY - no bam/fastq found)"
        : | gzip -"$GZIP_LEVEL" > "$target"
    fi

    n=$(( $(zcat "$target" | wc -l) / 4 ))
    TOTAL=$(( TOTAL + n ))
    printf '      %-16s %12d reads\n' "$out" "$n"
done

echo
printf '  TOTAL %d reads written\n' "$TOTAL"
echo
echo "=== done ==="
echo "Feed the pipeline:  RUNDATA=$DEST"
echo "Sample sheet:       samples_EXP26000892.tsv  (barcode column = custom_bc01..04)"
