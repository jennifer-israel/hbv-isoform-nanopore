#!/usr/bin/env python3
"""
Phase 3 — classify HBV reads into transcript/protein classes (cDNA, 2× reference).

Input per library: analysis/samples/<sample_key>/hbv.umi.bam
  — pychopper-ORIENTED HBV reads (forward = transcript sense), re-aligned to the
    HBV-only 2× contig (U95551.1_2x, 6364 bp), carrying the UMI in the RX tag.

Why this input (not the raw genome BAM): cDNA reads occur in both orientations, so a
raw minus-strand alignment's ref_start is the transcript 3′ end, not its TSS. pychopper
resolves orientation, so on hbv.umi.bam a +-strand alignment's ref_start IS the 5′ TSS.

2× coordinate handling (decided): coordinates are KEPT in 0–6363 space (the wrap stays
visible for plots). `mod 3182` is used ONLY to derive the genomic TSS / 3′ end for
binning — never to rewrite stored coordinates.

Classification (adapted from EXP26000465 phase3_classify.py; same U95551.1 windows).
Each transcript class needs BOTH a 5′ TSS in its window AND a plausible transcript
length — a read too short to be that transcript (a fragment) falls through to unclassified:
  precore / pgRNA  TSS in [1730,1880] AND crosses the 3182 junction AND span ≥ 2600 bp
                   (they wrap the whole genome, ~3.3 kb); precore if TSS ≤ 1815, else pgRNA
  preS1            TSS in [2700,3100] AND span ≥ 1600 bp
  preS2_S          TSS ≥ 3100 or ≤ 150 (circular) AND span ≥ 1200 bp
  X                TSS in [1260,1450] AND span in [300,1000] (X mRNA ~0.7 kb)
  pgRNA_RT         span ≥ ~3982 bp (tandem/concatemeric readthrough; >1 genome copy)
  spliced          carries a real intron (ref gap > MIN_INTRON=200 bp) → grouped here
                   regardless of TSS; the `splice_junction` column flags canonical SP1
                   (donor ~2447 / acceptor ~489, the minor HBSP-encoding variant) vs
                   non-canonical. Major HBV mRNAs above are unspliced.
  unclassified     none of the above (fragments failing a length gate, truncated 5′, ambiguous)
  antisense        aligned to the minus strand after pychopper orientation (flagged, not TSS-binned)

Per-read exon-block coordinates (2× space, small indels merged) are stored in the `blocks`
column so downstream plots can render intron gaps for spliced reads.

Outputs per library: hbv_classified.tsv (+ .parquet) — one row per read.
Cross-library: analysis/comparison/phase3_classification_counts.tsv (reads per class × condition).

Usage:
    conda activate hbv_lr_analysis
    python scripts/phase3_classify.py                       # all libraries
    python scripts/phase3_classify.py --sample-key <key>    # subset
"""
import argparse, datetime, sys
from pathlib import Path
import pysam
import pandas as pd

PROJECT_ROOT = Path("/data/EXP26000993")
ANALYSIS = PROJECT_ROOT / "analysis"
SAMPLES  = PROJECT_ROOT / "config" / "samples.tsv"
VERSION  = "phase3-v1"

# --- transcript classifier: use the CANONICAL shared implementation (single source of
# truth across all HBV projects). Do NOT re-inline the rules here; edit shared/ instead. ---
sys.path.insert(0, str(PROJECT_ROOT.parent / "shared"))
from hbv_transcript_classify import classify, VERSION as CLASSIFY_VERSION  # noqa: E402


def process_library(sample_key, barcode):
    d = ANALYSIS / "samples" / sample_key
    bam = d / "hbv.umi.bam"
    if not bam.exists():
        sys.stderr.write(f"SKIP {sample_key}: no hbv.umi.bam\n")
        return None
    rows = []
    with pysam.AlignmentFile(bam) as af:
        for r in af.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            umi = r.get_tag("RX") if r.has_tag("RX") else ""
            c = classify(r.reference_start, r.reference_end, r.get_blocks(), r.is_reverse)
            rows.append({
                "read_id": r.query_name, "umi": umi,
                "strand": "-" if r.is_reverse else "+",
                "mapq": r.mapping_quality,
                "read_length": r.query_length or 0,
                "ref_start": r.reference_start, "ref_end": r.reference_end,  # 2× coords, kept
                "span_len": r.reference_end - r.reference_start,
                **c,
            })
    if not rows:
        return None
    df = pd.DataFrame(rows)
    df.insert(0, "sample_key", sample_key)
    df.insert(1, "barcode", barcode)
    hdr = f"# sample={sample_key}; pipeline_version={VERSION}; generated={datetime.date.today()}"
    with open(d / "hbv_classified.tsv", "w") as fh:
        fh.write(hdr + "\n"); df.to_csv(fh, sep="\t", index=False)
    df.to_parquet(d / "hbv_classified.parquet", index=False)
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-key", action="append")
    args = ap.parse_args()

    libs = []
    for line in SAMPLES.read_text().splitlines():
        if line.startswith("#") or line.startswith("barcode\t") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) >= 4 and (not args.sample_key or f[3] in set(args.sample_key)):
            libs.append((f[3], f[0]))

    (ANALYSIS / "comparison").mkdir(parents=True, exist_ok=True)
    all_rows = []
    for sk, bc in libs:
        df = process_library(sk, bc)
        if df is None:
            print(f"{sk}: 0 HBV UMI reads"); continue
        vc = df["tx_class"].value_counts()
        print(f"{sk}: {len(df)} reads — " + ", ".join(f"{k}={v}" for k, v in vc.items()))
        for cls, n in vc.items():
            all_rows.append({"sample_key": sk, "barcode": bc, "tx_class": cls, "n_reads": int(n)})

    if all_rows:
        out = ANALYSIS / "comparison" / "phase3_classification_counts.tsv"
        with open(out, "w") as fh:
            fh.write(f"# experiment=EXP26000892_cDNA003; pipeline_version={VERSION}; generated={datetime.date.today()}\n")
            pd.DataFrame(all_rows).to_csv(fh, sep="\t", index=False)
        print(f"\nWritten: {out}")


if __name__ == "__main__":
    main()
