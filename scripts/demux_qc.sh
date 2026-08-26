#!/usr/bin/env bash
# Demux QC — EXP26000892. One streaming pass over dorado's sequencing_summary.txt.
#
#   bash demux_qc.sh /data/EXP26000896/demux/strict/sequencing_summary.txt [out.tsv]
#
# The summary is ~46 GB, so this makes a SINGLE awk pass and holds only per-barcode
# HISTOGRAMS in memory (never per-read arrays) — same approach as your
# phase1_5_readqc.sh. Runtime is disk-bound, roughly 10-25 min.
#
# Columns are located BY NAME from the header, so this survives dorado changing
# column order between releases.
#
# What it answers, in order of how much it should change your mind:
#
#  1. UNCLASSIFIED COMPOSITION. Are unclassified reads shorter than classified?
#     Shorter  -> benign truncation; strict both-ends correctly rejected them.
#     Same len -> your thresholds are rejecting good reads and the loss is
#                 non-random. Consider raising max_barcode_penalty.
#
#  2. SCORE DISTRIBUTION SHAPE. You want BIMODAL: a tight cluster of confident
#     hits, then a gap, then nothing. A smooth continuum piled up against the
#     cutoff means the threshold is arbitrary and read counts will move around
#     with small parameter changes.
#
#  3. FRONT vs REAR AGREEMENT. Both ends are scored even in strict mode. Reads
#     where one end scores well and the other poorly are chimera candidates —
#     the index-hopping signal, without needing the permissive pass.
#
#  4. BARCODE POSITION. barcode_front_begin_index tells you where the barcode
#     actually sits. The TOML assumes the construct ends by ~80 bp (adapter
#     trimmed) or ~140 bp (untrimmed), with windows of 200. If observed
#     positions crowd the window edge, widen the windows.
set -euo pipefail

SUMMARY=${1:?usage: demux_qc.sh <sequencing_summary.txt> [out.tsv]}
OUT=${2:-demux_qc_summary.tsv}
[ -s "$SUMMARY" ] || { echo "ERROR: not found or empty: $SUMMARY" >&2; exit 1; }

echo "=== demux QC ==="
echo "  summary: $SUMMARY ($(du -h "$SUMMARY" | cut -f1))"
echo "  output:  $OUT"
echo "  single streaming pass, histograms only — this takes a while."
echo

awk -F'\t' -v OUT="$OUT" '
function pct(bc, which, binw, frac, total,   i, cum, goal, mx, k, c) {
    # walk a histogram ascending to the bin covering `frac` of `total`
    goal = frac * total; cum = 0; mx = maxbin[bc SUBSEP which]
    for (i = 0; i <= mx; i++) {
        k = bc SUBSEP which SUBSEP i
        c = H[k]; cum += c
        if (cum >= goal) return i * binw + binw/2
    }
    return 0
}
function add(bc, which, binw, val,   b) {
    if (val == "" || val == "-" || val !~ /^-?[0-9.]+$/) return
    b = int(val / binw); if (b < 0) b = 0
    H[bc SUBSEP which SUBSEP b]++
    if (b > maxbin[bc SUBSEP which]) maxbin[bc SUBSEP which] = b
    S[bc SUBSEP which] += val; N[bc SUBSEP which]++
    if (!(bc SUBSEP which in MIN) || val+0 < MIN[bc SUBSEP which]) MIN[bc SUBSEP which] = val+0
    if (val+0 > MAX[bc SUBSEP which]) MAX[bc SUBSEP which] = val+0
}

NR == 1 {
    for (i = 1; i <= NF; i++) col[$i] = i
    # required
    split("sequence_length_template mean_qscore_template barcode_arrangement", req, " ")
    for (j in req) if (!(req[j] in col)) { print "ERROR: missing column " req[j] > "/dev/stderr"; exit 1 }
    c_len=col["sequence_length_template"]; c_q=col["mean_qscore_template"]
    c_bc=col["barcode_arrangement"]
    c_sc=("barcode_score" in col)?col["barcode_score"]:0
    c_fs=("barcode_front_score" in col)?col["barcode_front_score"]:0
    c_rs=("barcode_rear_score" in col)?col["barcode_rear_score"]:0
    c_fb=("barcode_front_begin_index" in col)?col["barcode_front_begin_index"]:0
    c_re=("barcode_rear_end_index" in col)?col["barcode_rear_end_index"]:0
    c_pf=("passes_filtering" in col)?col["passes_filtering"]:0
    next
}
{
    bc = $c_bc; if (bc == "") bc = "(blank)"
    reads[bc]++; bases[bc] += $c_len
    add(bc, "L",  100, $c_len)
    add(bc, "Q",    1, $c_q)
    if (c_sc) add(bc, "S", 1, $c_sc)
    if (c_fs) { add(bc, "F", 1, $c_fs); if ($c_fs+0 > 0) nf[bc]++ }
    if (c_rs) { add(bc, "R", 1, $c_rs); if ($c_rs+0 > 0) nr_[bc]++ }
    if (c_fs && c_rs && $c_fs+0 > 0 && $c_rs+0 > 0) both[bc]++
    if (c_fb) add(bc, "B", 10, $c_fb)
    if (c_re) add(bc, "E", 10, $c_re)
    total++
}
END {
    printf("# experiment=EXP26000892; source=%s; generated=%s\n", "sequencing_summary", strftime("%Y-%m-%d")) > OUT
    printf("barcode\treads\tpct_of_total\ttotal_bases\tlen_median\tlen_p10\tlen_p90\tmean_q\t") >> OUT
    printf("score_median\tscore_p10\tscore_min\tscore_max\tfront_median\trear_median\t") >> OUT
    printf("frac_front_scored\tfrac_rear_scored\tfrac_both_scored\tfront_begin_median\tfront_begin_p90\trear_end_median\n") >> OUT

    for (bc in reads) {
        n = reads[bc]
        lm = pct(bc,"L",100,0.5,n); l10 = pct(bc,"L",100,0.1,n); l90 = pct(bc,"L",100,0.9,n)
        mq = (N[bc SUBSEP "Q"] > 0) ? S[bc SUBSEP "Q"]/N[bc SUBSEP "Q"] : 0
        ns = N[bc SUBSEP "S"]; sm = ns ? pct(bc,"S",1,0.5,ns) : 0; s10 = ns ? pct(bc,"S",1,0.1,ns) : 0
        nfs = N[bc SUBSEP "F"]; fm = nfs ? pct(bc,"F",1,0.5,nfs) : 0
        nrs = N[bc SUBSEP "R"]; rm = nrs ? pct(bc,"R",1,0.5,nrs) : 0
        nb = N[bc SUBSEP "B"]; bm = nb ? pct(bc,"B",10,0.5,nb) : 0; b90 = nb ? pct(bc,"B",10,0.9,nb) : 0
        ne = N[bc SUBSEP "E"]; em = ne ? pct(bc,"E",10,0.5,ne) : 0
        printf("%s\t%d\t%.3f\t%d\t%d\t%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%d\t%d\t%d\n",
               bc, n, 100*n/total, bases[bc], lm, l10, l90, mq,
               sm, s10, MIN[bc SUBSEP "S"], MAX[bc SUBSEP "S"], fm, rm,
               n?nf[bc]/n:0, n?nr_[bc]/n:0, n?both[bc]/n:0, bm, b90, em) >> OUT
    }

    # ---- score histograms, for judging bimodality by eye ----
    print "" >> OUT
    print "# barcode_score histogram (bin=1) per barcode: barcode<TAB>score_bin<TAB>count" >> OUT
    for (k in H) {
        split(k, p, SUBSEP)
        if (p[2] == "S") printf("#HIST_SCORE\t%s\t%d\t%d\n", p[1], p[3], H[k]) >> OUT
    }
    print "" >> OUT
    print "# read-length histogram (bin=100) per barcode" >> OUT
    for (k in H) {
        split(k, p, SUBSEP)
        if (p[2] == "L") printf("#HIST_LEN\t%s\t%d\t%d\n", p[1], p[3]*100, H[k]) >> OUT
    }
    printf("\n# total_records\t%d\n", total) >> OUT
}
' "$SUMMARY"

echo
echo "=== per-barcode summary ==="
grep -v '^#' "$OUT" | grep -v '^$' | column -t -s$'\t'

echo
echo "=== interpretation ==="
python3 - "$OUT" <<'PYEOF'
import sys, csv
rows=[]
with open(sys.argv[1]) as fh:
    hdr=None
    for line in fh:
        if line.startswith("#") or not line.strip(): continue
        f=line.rstrip("\n").split("\t")
        if hdr is None: hdr=f; continue
        rows.append(dict(zip(hdr,f)))
if not rows:
    print("  (no rows parsed)"); raise SystemExit

def num(r,k,d=0.0):
    try: return float(r[k])
    except: return d

cls=[r for r in rows if "unclassified" not in r["barcode"].lower()]
unc=[r for r in rows if "unclassified" in r["barcode"].lower()]

print("\n1. UNCLASSIFIED COMPOSITION")
if unc and cls:
    u=unc[0]; ul=num(u,"len_median")
    cl=sum(num(r,"len_median")*num(r,"reads") for r in cls)/max(1,sum(num(r,"reads") for r in cls))
    print(f"   unclassified median length : {ul:,.0f} bp")
    print(f"   classified   median length : {cl:,.0f} bp")
    if ul < 0.7*cl:
        print("   -> unclassified are markedly SHORTER. Benign: truncated reads lack one")
        print("      end's barcode, which strict both-ends is supposed to reject.")
    elif ul > 0.9*cl:
        print("   -> unclassified are the SAME length as classified. That is a warning:")
        print("      you are rejecting full-length reads, so the loss is not just")
        print("      truncation. Consider raising max_barcode_penalty and re-checking.")
    else:
        print("   -> intermediate. Mixed truncation and threshold loss.")
else:
    print("   (no unclassified row found)")

print("\n2. SCORE DISTRIBUTION")
for r in cls:
    print(f"   {r['barcode']:<26} median={num(r,'score_median'):>6.0f}  p10={num(r,'score_p10'):>6.0f}"
          f"  range=[{num(r,'score_min'):.0f}, {num(r,'score_max'):.0f}]")
print("   -> Check direction: if higher score = better match, you want median high and")
print("      p10 not far below it (tight, confident cluster). A long tail down toward")
print("      the cutoff means the threshold is slicing through a continuum.")
print("      Plot the #HIST_SCORE rows to see whether it is genuinely bimodal.")

print("\n3. FRONT vs REAR AGREEMENT (chimera / index-hopping signal)")
for r in cls:
    fb,rb,bo=num(r,"frac_front_scored"),num(r,"frac_rear_scored"),num(r,"frac_both_scored")
    print(f"   {r['barcode']:<26} front={fb:6.1%}  rear={rb:6.1%}  both={bo:6.1%}")
print("   -> In strict mode both ends should be scored for nearly every read. A barcode")
print("      with a notably lower 'both' fraction is enriched for single-ended or")
print("      chimeric molecules. Watch custom_bc03 (1 ng) and custom_bc04 (0.1 ng):")
print("      leakage from the high-input libraries shows up here first.")

print("\n4. BARCODE POSITION vs TOML WINDOWS")
for r in cls:
    print(f"   {r['barcode']:<26} front_begin median={num(r,'front_begin_median'):>5.0f}  "
          f"p90={num(r,'front_begin_p90'):>5.0f}   rear_end median={num(r,'rear_end_median'):>7.0f}")
print("   -> TOML uses front_barcode_window=200, barcode_end_proximity=200. If the p90")
print("      of front_begin approaches 200 you are clipping real barcodes and should")
print("      widen the windows. Well under 200 means there is comfortable headroom.")
PYEOF

echo
echo "Full table + histograms: $OUT"
echo "To plot the score histogram:  grep '^#HIST_SCORE' $OUT | cut -f2-"
