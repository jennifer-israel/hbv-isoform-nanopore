#!/usr/bin/env bash
# Phase 1.5A — per-library read QC from the sequencing_summary (one streaming pass).
#
# Produces, per barcode (passes_filtering==TRUE only):
#   reads, total_bases, len_N50, len_median, mean_qscore,
#   polyA tail: frac_with_tail, median_tail_len (tail>0)
# grouped by barcode_arrangement, joined to config/samples.tsv for sample_key.
#
#   conda activate hbv_lr_analysis
#   bash scripts/phase1_5_readqc.sh
#
# The summary is ~36 GB; this makes ONE awk pass and holds only per-barcode
# length/tail HISTOGRAMS (memory-safe — never stores per-read arrays). ~10-20 min
# on the big box. Output: analysis/comparison/phase1_5_read_qc.tsv
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/data/EXP26000896}"
RUNDATA=$PROJECT_ROOT/rundata
ANALYSIS_ROOT=$PROJECT_ROOT/analysis
SAMPLES=$PROJECT_ROOT/config/samples.tsv
OUT=$ANALYSIS_ROOT/comparison/phase1_5_read_qc.tsv
PIPELINE_VERSION=phase1_5-v1

SS=$(ls "$RUNDATA"/sequencing_summary_*.txt 2>/dev/null | head -1)
[ -f "$SS" ] || { echo "ERROR: sequencing_summary not found under $RUNDATA"; exit 1; }
mkdir -p "$ANALYSIS_ROOT/comparison"

echo "[$(date)] Phase 1.5A read QC — streaming $(du -h "$SS" | cut -f1): $SS"

# Single awk pass. Columns (from header, 1-based):
#   13 passes_filtering  17 sequence_length_template  18 mean_qscore_template
#   19 poly_tail_length  30 barcode_arrangement
# Length histogram bin = 100 bp (N50/median good to ±100 bp). Tail histogram bin = 5 bp.
awk -F'\t' -v OUT="$OUT" -v SAMPLES="$SAMPLES" -v VER="$PIPELINE_VERSION" '
BEGIN {
  # load sample_key map: barcode -> lib_id, sample_name, sample_key
  while ((getline line < SAMPLES) > 0) {
    if (line ~ /^#/ || line ~ /^barcode\t/) continue
    n=split(line, f, "\t"); if (n<4) continue
    lib[f[1]]=f[2]; sname[f[1]]=f[3]; skey[f[1]]=f[4]
  }
  LB=100; TB=5
}
NR==1 {   # verify expected header positions
  if ($13!="passes_filtering"||$17!="sequence_length_template"||$18!="mean_qscore_template"||$19!="poly_tail_length"||$30!="barcode_arrangement") {
    print "ERROR: sequencing_summary columns not in expected positions" > "/dev/stderr"; exit 1
  }
  next
}
{
  if ($13!="TRUE") next
  bc=$30
  reads[bc]++; bases[bc]+=$17; qsum[bc]+=$18
  lh[bc,int($17/LB)]++            # length histogram
  if ($19+0>0) { twith[bc]++; th[bc,int($19/TB)]++ }   # polyA tail histogram (tail>0)
  if ($17+0>maxlb[bc]) maxlb[bc]=$17
  if ($19+0>maxtb[bc]) maxtb[bc]=$19
}
function pct_from_hist(bcx, arrname, binw, targetfrac, total,    i,cum,goal,mx) {
  # generic: walk histogram (ascending bins) to the bin covering targetfrac of `total`
  goal=targetfrac*total; cum=0
  mx = (arrname=="L") ? int(maxlb[bcx]/binw) : int(maxtb[bcx]/binw)
  for (i=0;i<=mx;i++){ k=bcx SUBSEP i; c=(arrname=="L")?lh[k]:th[k]; cum+=c; if(cum>=goal) return (i*binw + binw/2) }
  return 0
}
function n50(bcx,    i,cum,goal,mx) {   # N50 by BASES over length histogram
  goal=bases[bcx]/2; cum=0; mx=int(maxlb[bcx]/LB)
  for (i=mx;i>=0;i--){ mid=i*LB+LB/2; cum+=lh[bcx SUBSEP i]*mid; if(cum>=goal) return mid }
  return 0
}
END {
  printf("# experiment=EXP26000892_cDNA003; pipeline_version=%s; generated=%s\n", VER, strftime("%Y-%m-%d")) > OUT
  printf("sample_key\tbarcode\tlib_id\treads\ttotal_bases\tlen_N50\tlen_median\tmean_qscore\tfrac_with_polyA\tmedian_polyA_len\n") >> OUT
  for (bc in reads) {
    sk=(bc in skey)?skey[bc]:bc; lb=(bc in lib)?lib[bc]:"NA"
    mq=(reads[bc]>0)?qsum[bc]/reads[bc]:0
    med=pct_from_hist(bc,"L",LB,0.5,reads[bc])
    fwt=(reads[bc]>0)?twith[bc]/reads[bc]:0
    tmed=(twith[bc]>0)?pct_from_hist(bc,"T",TB,0.5,twith[bc]):0
    printf("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.2f\t%.4f\t%d\n", sk,bc,lb,reads[bc],bases[bc],n50(bc),med,mq,fwt,tmed) >> OUT
  }
}
' "$SS"

echo "[$(date)] Phase 1.5A complete: $OUT"
# show only the real libraries (present in samples.tsv), sorted by yield
( head -2 "$OUT"; grep -Ff <(awk -F'\t' 'NR>2{print $1}' "$SAMPLES" 2>/dev/null; awk -F'\t' '$1!~/^#/&&$1!="barcode"{print $4}' "$SAMPLES") "$OUT" 2>/dev/null | sort -t$'\t' -k5,5 -nr ) | column -t -s$'\t' || cat "$OUT"
