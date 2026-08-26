#!/usr/bin/env python3
"""
Phase 2 — extract HBV reads from each library's genome BAM, normalize 2× coords,
and emit per-read stats + per-library HBV recovery counts.

For each library (analysis/samples/<sample_key>/):
  aligned_sorted.bam  (Phase 1, all reads)  → HBV reads on U95551.1_2x
Outputs:
  hbv.bam(.bai)            HBV-only, coordinate-sorted (RETAINED, for IGV)
  hbv_per_read.tsv/.parquet  read_id, length, mean_qscore, mapq, strand,
                             ref_start, ref_end (2× coords), start_norm/end_norm
                             (mod 3182), spans_wrap
Shared: analysis/comparison/phase2_hbv_counts.tsv
  raw HBV reads (all + MAPQ≥20), unique HBV molecules (UMI, from hbv.umi.bam if
  present), composite primary-mapped, HBV per-million (raw + unique), Wilson 95% CI
  on the HBV fraction.

  conda activate hbv_lr_analysis
  python scripts/phase2_extract_hbv.py                 # all libraries
  python scripts/phase2_extract_hbv.py --sample-key SeqLib5543_150ng_NOpolyA_20
"""
import argparse, csv, math, subprocess, sys, datetime, shutil
from pathlib import Path

import pysam
import pandas as pd
try:
    from statsmodels.stats.proportion import proportion_confint
except Exception:
    proportion_confint = None

PROJECT_ROOT = Path("/data/EXP26000896")
ANALYSIS = PROJECT_ROOT / "analysis"
SAMPLES  = PROJECT_ROOT / "config" / "samples.tsv"
HBV_CONTIG = "U95551.1_2x"
HBV_LEN    = 3182          # single-copy length; 2× contig = 6364
MAPQ_MIN   = 20
VERSION    = "phase2-v1"


def load_samples():
    rows = []
    for line in SAMPLES.read_text().splitlines():
        if line.startswith("#") or line.startswith("barcode\t") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) >= 4:
            rows.append({"barcode": f[0], "lib_id": f[1], "sample_name": f[2], "sample_key": f[3]})
    return rows


def mean_qscore(quals):
    if not quals:
        return float("nan")
    err = sum(10 ** (-q / 10) for q in quals) / len(quals)
    return -10 * math.log10(err) if err > 0 else float("nan")


def wilson(k, n):
    if not n or proportion_confint is None:
        return (float("nan"), float("nan"))
    return proportion_confint(k, n, alpha=0.05, method="wilson")


def primary_mapped_from_flagstat(path: Path):
    if not path.exists():
        return None
    for line in path.read_text().splitlines():
        if "primary mapped" in line:
            return int(line.split()[0])
    return None


def unique_molecules(hbv_umi_bam: Path, work: Path):
    """umi_tools dedup on the HBV UMI BAM → unique molecule count. NA if unavailable.

    umi_tools requires equal-length UMIs; ONT indels make pychopper's structured UMI
    vary (mostly 28 nt, plus 27/29/RX:Z:None). Restrict to the well-formed 28-nt UMIs
    (same rule as Phase 1.5B) before dedup.
    """
    if not hbv_umi_bam.exists() or shutil.which("umi_tools") is None:
        return None
    filt = work / "hbv.umi28.bam"
    ded  = work / "hbv.umi.dedup.bam"
    try:
        subprocess.run(["samtools", "view", "-b", "-e", "length([RX])==28",
                        "-o", str(filt), str(hbv_umi_bam)], check=True)
        pysam.index(str(filt))
        subprocess.run(["umi_tools", "dedup", "-I", str(filt), "-S", str(ded),
                        "--extract-umi-method=tag", "--umi-tag=RX", "--method=directional"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        n = int(pysam.view("-c", "-F", "0x904", str(ded)).strip())
        for f in (filt, ded, Path(str(filt) + ".bai"), Path(str(ded) + ".bai")):
            f.unlink(missing_ok=True)
        return n
    except Exception as e:
        sys.stderr.write(f"  WARN umi_tools dedup failed on {hbv_umi_bam.name}: {e}\n")
        return None


def extract_one(s, counts_writer):
    sk = s["sample_key"]
    d = ANALYSIS / "samples" / sk
    bam = d / "aligned_sorted.bam"
    if not bam.exists():
        sys.stderr.write(f"SKIP {sk}: no aligned_sorted.bam (run Phase 1 first)\n")
        return
    print(f"[{datetime.datetime.now():%H:%M:%S}] Phase 2 {sk}")

    hbv_bam = d / "hbv.bam"
    rows = []
    raw_all = raw_q20 = 0
    with pysam.AlignmentFile(bam) as af, \
         pysam.AlignmentFile(hbv_bam, "wb", template=af) as out:
        if HBV_CONTIG not in af.references:
            sys.stderr.write(f"  WARN {HBV_CONTIG} not in {bam.name}\n")
        for r in af.fetch(HBV_CONTIG) if HBV_CONTIG in af.references else []:
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            out.write(r)
            raw_all += 1
            if r.mapping_quality >= MAPQ_MIN:
                raw_q20 += 1
            rs, re = r.reference_start, r.reference_end  # 0-based, end exclusive
            rows.append({
                "read_id": r.query_name,
                "length": r.query_length or (r.infer_read_length() or 0),
                "mean_qscore": round(mean_qscore(list(r.query_qualities or [])), 3),
                "mapq": r.mapping_quality,
                "strand": "-" if r.is_reverse else "+",
                "ref_start": rs, "ref_end": re,
                "start_norm": rs % HBV_LEN,
                "end_norm": (re - 1) % HBV_LEN + 1 if re else None,
                "spans_wrap": bool(rs is not None and re and rs < HBV_LEN <= re),
            })
    pysam.sort("-o", str(hbv_bam) + ".tmp", str(hbv_bam)); Path(str(hbv_bam) + ".tmp").replace(hbv_bam)
    pysam.index(str(hbv_bam))

    df = pd.DataFrame(rows)
    hdr = f"# sample={sk}; pipeline_version={VERSION}; generated={datetime.date.today()}"
    tsv = d / "hbv_per_read.tsv"
    with open(tsv, "w") as fh:
        fh.write(hdr + "\n")
        df.to_csv(fh, sep="\t", index=False)
    if not df.empty:
        df.to_parquet(d / "hbv_per_read.parquet", index=False)

    comp_mapped = primary_mapped_from_flagstat(d / "flagstat.txt")
    uniq = unique_molecules(d / "hbv.umi.bam", d)
    def per_m(k):
        return round(k / comp_mapped * 1e6, 3) if comp_mapped else float("nan")
    lo, hi = wilson(raw_all, comp_mapped) if comp_mapped else (float("nan"), float("nan"))

    counts_writer.writerow([
        sk, s["barcode"], s["lib_id"], raw_all, raw_q20,
        uniq if uniq is not None else "NA",
        comp_mapped if comp_mapped is not None else "NA",
        per_m(raw_all), per_m(uniq) if uniq is not None else "NA",
        f"{lo:.3e}", f"{hi:.3e}",
    ])
    print(f"  HBV reads {raw_all} (MAPQ≥{MAPQ_MIN} {raw_q20}); unique {uniq}; "
          f"per-million {per_m(raw_all)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-key", action="append", help="restrict to sample_key(s)")
    args = ap.parse_args()

    samples = load_samples()
    if args.sample_key:
        samples = [s for s in samples if s["sample_key"] in set(args.sample_key)]
    (ANALYSIS / "comparison").mkdir(parents=True, exist_ok=True)
    out_tsv = ANALYSIS / "comparison" / "phase2_hbv_counts.tsv"
    with open(out_tsv, "w", newline="") as fh:
        fh.write(f"# experiment=EXP26000892_cDNA003; pipeline_version={VERSION}; generated={datetime.date.today()}\n")
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_key", "barcode", "lib_id", "hbv_reads_all", f"hbv_reads_mapq{MAPQ_MIN}",
                    "hbv_unique_molecules", "composite_primary_mapped",
                    "hbv_per_million_raw", "hbv_per_million_unique",
                    "hbv_frac_wilson_lo", "hbv_frac_wilson_hi"])
        for s in samples:
            extract_one(s, w)
    print(f"\nPhase 2 complete. Counts: {out_tsv}")


if __name__ == "__main__":
    main()
