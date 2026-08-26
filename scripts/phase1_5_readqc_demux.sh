#!/usr/bin/env bash
# Phase 1.5A — per-library read QC (one streaming pass).  EXP26000892 version.
#
# CHANGES vs the EXP26000559 original, all forced by this run:
#
#   INPUT. The original reads $RUNDATA/sequencing_summary_*.txt. This run was
#   demultiplexed computationally, so RUNDATA holds only the reshaped fastqs and has
#   no summary. The per-library barcode assignment lives in dorado demux's summary,
#   which is what we read instead.
#
#   COLUMN POSITIONS. The demux summary orders columns differently from the MinKNOW
#   summary the original was written against:
#       passes_filtering          13 -> 11
#       sequence_length_template  17 -> 15
#       mean_qscore_template      18 -> 16
#       barcode_arrangement       30 -> 23
#
#   POLY(A) REMOVED. There is no poly_tail_length column in either the demux summary
#   or the basecaller summary for this run — poly(A) estimation was not enabled at
#   basecall time, so the data does not exist. frac_with_polyA / median_polyA_len are
#   therefore dropped. phase5_report.py's Figure 1 must drop that panel to match.
#
#   BARCODE NAMES. dorado reports EXP26000892_barcode01..04 (or barcode01..04). These
#   are mapped to the custom_bc01..04 keys used in config/samples.tsv.
#
#   conda activate hbv_lr
#   bash phase1_5_readqc.sh
#
# Memory-safe: one awk pass holding only per-barcode length histograms, never
# per-read arrays. The summary is ~46 GB; expect 10-25 min.
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/data/EXP26000993}"
ANALYSIS_ROOT=$PROJECT_ROOT/analysis
SAMPLES=$PROJECT_ROOT/config/samples.tsv
SS=${SS:-$PROJECT_ROOT/demux/strict/sequencing_summary.txt}
OUT=$ANALYSIS_ROOT/comparison/phase1_5_read_qc.tsv
PIPELINE_VERSION=phase1_5-v2-demux

[ -f "$SS" ] || { echo "ERROR: demux summary not found: $SS" >&2; exit 1; }
mkdir -p "$ANALYSIS_ROOT/comparison"

echo "[$(date)] Phase 1.5A read QC — streaming $(du -h "$SS" | cut -f1): $SS"

awk -F'\t' -v OUT="$OUT" -v SAMPLES="$SAMPLES" -v VER="$PIPELINE_VERSION" '
BEGIN {
  # barcode -> lib_id, sample_name, sample_key  from config/samples.tsv
  while ((getline line < SAMPLES) > 0) {
    if (line ~ /^#/ || line ~ /^barcode\t/) continue
    n = split(line, f, "\t"); if (n < 4) continue
    lib[f[1]] = f[2]; sname[f[1]] = f[3]; skey[f[1]] = f[4]
  }
  LB = 100      # length histogram bin (N50/median good to +/-100 bp)
}
NR == 1 {
  # locate columns BY NAME so a future dorado release reordering them is not silent
  for (i = 1; i <= NF; i++) col[$i] = i
  c_pass = col["passes_filtering"]
  c_len  = col["sequence_length_template"]
  c_q    = col["mean_qscore_template"]
  c_bc   = col["barcode_arrangement"]
  if (!c_len || !c_q || !c_bc) {
    print "ERROR: required columns not found in " FILENAME > "/dev/stderr"
    print "  need sequence_length_template, mean_qscore_template, barcode_arrangement" > "/dev/stderr"
    exit 1
  }
  next
}
{
  # passes_filtering may be absent or already filtered upstream; only test if present
  if (c_pass && $c_pass != "TRUE" && $c_pass != "true" && $c_pass != "1") next

  raw = $c_bc
  if (raw == "" || raw ~ /unclassified/) next

  # EXP26000892_barcode03 / barcode03 -> custom_bc03
  bc = raw
  if (match(bc, /barcode0?[0-9]+$/)) {
    num = substr(bc, RSTART); gsub(/barcode/, "", num)
    bc = sprintf("custom_bc%02d", num + 0)
  }

  reads[bc]++; bases[bc] += $c_len; qsum[bc] += $c_q
  lh[bc, int($c_len / LB)]++
  if ($c_len + 0 > maxlb[bc]) maxlb[bc] = $c_len
}
function pct_from_hist(bcx, targetfrac, total,   i, cum, goal, mx, k) {
  goal = targetfrac * total; cum = 0; mx = int(maxlb[bcx] / LB)
  for (i = 0; i <= mx; i++) { k = bcx SUBSEP i; cum += lh[k]; if (cum >= goal) return i * LB + LB / 2 }
  return 0
}
function n50(bcx,   i, cum, goal, mx, mid) {     # N50 by BASES over the length histogram
  goal = bases[bcx] / 2; cum = 0; mx = int(maxlb[bcx] / LB)
  for (i = mx; i >= 0; i--) { mid = i * LB + LB / 2; cum += lh[bcx SUBSEP i] * mid; if (cum >= goal) return mid }
  return 0
}
END {
  printf("# experiment=EXP26000892_cDNA003; pipeline_version=%s; source=demux_summary; generated=%s\n",
         VER, strftime("%Y-%m-%d")) > OUT
  printf("sample_key\tbarcode\tlib_id\treads\ttotal_bases\tlen_N50\tlen_median\tlen_p10\tlen_p90\tmean_qscore\n") >> OUT
  for (bc in reads) {
    sk = (bc in skey) ? skey[bc] : bc
    lb = (bc in lib)  ? lib[bc]  : "NA"
    mq = (reads[bc] > 0) ? qsum[bc] / reads[bc] : 0
    printf("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n",
           sk, bc, lb, reads[bc], bases[bc], n50(bc),
           pct_from_hist(bc, 0.5, reads[bc]),
           pct_from_hist(bc, 0.1, reads[bc]),
           pct_from_hist(bc, 0.9, reads[bc]), mq) >> OUT
  }
}
' "$SS"

echo "[$(date)] Phase 1.5A complete: $OUT"
echo
column -t -s$'\t' "$OUT" | grep -v '^#'
echo
echo "NOTE: no poly(A) columns — poly_tail_length is absent from this run's summaries"
echo "      (poly(A) estimation was not enabled at basecalling). phase5_report.py's"
echo "      fig_qc must drop its 'Poly(A) detection' panel accordingly."
