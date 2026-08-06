#!/usr/bin/env python3
"""
Build the hg38 + 2×HBV composite reference for EXP26000559_cDNA001.

Reads a single-copy hg38+HBV composite FASTA and writes a new FASTA where the
HBV contig U95551.1 is replaced by U95551.1_2x — the same 3,182 bp sequence
concatenated with itself (6,364 bp total). All hg38 contigs pass through
unchanged.

Why 2×: HBV pgRNA and other reads that wrap the circular genome are longer than
the 3,182 bp linearised reference. Against a single copy they split into two
alignment records at the linearisation point (supplementary alignments). Against
the doubled reference they align as one continuous record spanning the junction,
which makes the wrap unambiguous in IGV and in coverage plots, and simplifies
downstream per-read handling.

Coordinate mapping (used in Phase 2):
    positions 0–3181   → first copy  (== original coordinates)
    positions 3182–6363 → second copy (wrapped portion)
    normalise back to 0–3181 space with `coord mod 3182`.

Adapted from EXP26000465 scripts/make_hbv_2x_ref.py; made self-contained for
this project (explicit --in-fasta / --out-fasta, defaulting to CLAUDE.md anchors).

Usage
-----
    conda activate hbv_lr_analysis
    python scripts/make_hbv_2x_ref.py            # uses default anchor paths
    python scripts/make_hbv_2x_ref.py --in-fasta <src> --out-fasta <dst>

Default input  ($REF_SRC): the single-copy hg38+HBV composite from EXP26000465.
Default output ($REF_2X):  analysis/refs/hg38_hbv_2x.fa in this project.
"""

import argparse
import datetime
import sys
from pathlib import Path

HBV_CONTIG    = "U95551.1"
HBV_CONTIG_2X = "U95551.1_2x"
LINE_WIDTH    = 60

# Anchor defaults (keep in sync with CLAUDE.md)
PROJECT_ROOT   = Path("/home/ubuntu/matt_wolpert_claude_code_analysis/2026_07_02_EXP26000559_cDNA001")
DEFAULT_IN     = Path("/home/ubuntu/matt_wolpert_claude_code_analysis/"
                      "2026_05_14_HBV_LR_Transcript_Detect_EXP26000465/refs/hg38_hbv.fa")
DEFAULT_OUT    = PROJECT_ROOT / "analysis" / "refs" / "hg38_hbv_2x.fa"


def _write_doubled(fout, seq: str) -> int:
    doubled = seq + seq
    fout.write(f">{HBV_CONTIG_2X}  len={len(doubled)} [{HBV_CONTIG} x2]\n")
    for i in range(0, len(doubled), LINE_WIDTH):
        fout.write(doubled[i:i + LINE_WIDTH] + "\n")
    return len(doubled)


def make_2x_ref(in_path: Path, out_path: Path) -> None:
    print(f"[{datetime.datetime.now():%H:%M:%S}] Reading {in_path}")
    print(f"  Output → {out_path}")
    if not in_path.exists():
        print(f"ERROR: input FASTA not found: {in_path}", file=sys.stderr)
        sys.exit(1)

    hbv_found = False
    in_hbv    = False
    hbv_seq   = []

    with open(in_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                if in_hbv and hbv_seq:                      # flush buffered HBV
                    n = _write_doubled(fout, "".join(hbv_seq))
                    print(f"  Written {HBV_CONTIG_2X}: {n:,} bp "
                          f"(original {n // 2:,} bp × 2)")
                    hbv_seq = []
                contig = line.strip().lstrip(">").split()[0]
                if contig == HBV_CONTIG:
                    in_hbv = hbv_found = True
                else:
                    in_hbv = False
                    fout.write(line)
            else:
                (hbv_seq.append(line.strip()) if in_hbv else fout.write(line))

        if in_hbv and hbv_seq:                              # flush if HBV was last
            n = _write_doubled(fout, "".join(hbv_seq))
            print(f"  Written {HBV_CONTIG_2X}: {n:,} bp "
                  f"(original {n // 2:,} bp × 2)")

    if not hbv_found:
        print(f"ERROR: contig {HBV_CONTIG} not found in {in_path}", file=sys.stderr)
        sys.exit(1)

    size_gb = out_path.stat().st_size / 1e9
    print(f"[{datetime.datetime.now():%H:%M:%S}] Done. Output size: {size_gb:.2f} GB")
    print("\nNext: build the splice index (Phase 0):")
    print(f"  minimap2 -x splice -k14 -d {out_path.parent}/hg38_hbv_2x_splice.mmi {out_path}")


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in-fasta",  default=str(DEFAULT_IN),
                   help="single-copy hg38+HBV composite (default: $REF_SRC)")
    p.add_argument("--out-fasta", default=str(DEFAULT_OUT),
                   help="output hg38+2×HBV FASTA (default: $REF_2X)")
    args = p.parse_args()

    out_path = Path(args.out_fasta)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    make_2x_ref(Path(args.in_fasta), out_path)


if __name__ == "__main__":
    main()
