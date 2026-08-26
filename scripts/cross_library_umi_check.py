#!/usr/bin/env python3
"""
Cross-library UMI sharing — index hopping / carry-over check (EXP26000892).

THE QUESTION
Phase 2 gave 3,800 unique HBV molecules at 1 ng and 3,104 at 0.1 ng — a 10x input
difference producing a 1.2x molecule difference. Either the assay floors out around
~3,000 molecules, or the low-input libraries are partly filled with molecules that
actually originated in the 100 ng library (70,943 molecules) and leaked across during
the pooled hybridisation capture and post-capture PCR.

THE TEST
A molecule is identified by its UMI *and* the genomic position it aligns to. If the
SAME 28-nt UMI appears at the SAME position in two libraries, that is one physical
molecule counted twice — leakage. Two independently-captured molecules would have to
collide in a ~43 M UMI space AND land at the same coordinate.

Crucially this needs no external null. Chance UMI collisions land at RANDOM positions;
real leaked molecules land at the SAME position. So we compare, among UMIs shared
between two libraries, how often positions coincide, against how often positions
coincide for random cross-library molecule pairs. The excess is the leakage estimate.

  estimated leaked molecules = n_shared_umi * (p_position_match_observed
                                               - p_position_match_null)

Reported as a fraction of the SMALLER library, since contamination flows from
abundant to sparse and matters most where the recipient has few molecules of its own.

INPUT
analysis/samples/<sample_key>/hbv.umi.bam  (retained by phase2_hbv_umi_v2.sh)
Reads are collapsed to unique (UMI, position) molecules first — i.e. dedup by exact
match, which is stricter than umi_tools directional and is the right unit here.

  conda activate hbv_lr
  python cross_library_umi_check.py
  python cross_library_umi_check.py --tol 50 --out-tsv umi_sharing.tsv
"""
import argparse, itertools, random, sys
from collections import defaultdict
from pathlib import Path

import pysam

PROJECT_ROOT = Path("/data/EXP26000993")
ANALYSIS = PROJECT_ROOT / "analysis"
SAMPLES = PROJECT_ROOT / "config" / "samples.tsv"
UMI_LEN = 28


def load_samples():
    rows = []
    for line in SAMPLES.read_text().splitlines():
        if line.startswith("#") or line.startswith("barcode\t") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) >= 5:
            rows.append({"barcode": f[0], "lib_id": f[1], "sample_key": f[3], "input_ng": f[4]})
    return rows


def collapse(positions, tol):
    """Merge alignment starts within `tol` into one molecule.

    Reads from a single molecule do NOT share an exact start coordinate — nanopore
    alignment ends jitter by a few bp. Collapsing on exact position would count one
    molecule many times and inflate every per-library denominator.
    """
    out = []
    for p in sorted(positions):
        if not out or p - out[-1] > tol:
            out.append(p)
    return out


def load_molecules(sample_key, tol):
    """-> dict umi -> list of distinct positions (one entry per physical molecule)."""
    bam = ANALYSIS / "samples" / sample_key / "hbv.umi.bam"
    if not bam.exists():
        sys.stderr.write(f"  SKIP {sample_key}: no hbv.umi.bam\n")
        return None
    raw = defaultdict(set)
    n_reads = 0
    with pysam.AlignmentFile(bam) as af:
        for r in af.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            if not r.has_tag("RX"):
                continue
            umi = r.get_tag("RX")
            if len(umi) != UMI_LEN:
                continue
            n_reads += 1
            raw[umi].add(r.reference_start)
    mols = {u: collapse(p, tol) for u, p in raw.items()}
    n_mol = sum(len(v) for v in mols.values())
    print(f"  {sample_key:<30} {n_reads:>10,} reads -> {n_mol:>8,} molecules "
          f"({len(mols):,} distinct UMIs)")
    return mols


def pos_match(pa, pb, tol):
    """any position in pa within tol of any in pb"""
    for a in pa:
        for b in pb:
            if abs(a - b) <= tol:
                return True
    return False


def null_pos_match_rate(A, B, tol, n_draws, rng):
    """P(positions coincide) for RANDOM cross-library molecule pairs — the chance rate."""
    ua, ub = list(A.keys()), list(B.keys())
    if not ua or not ub:
        return 0.0
    hits = 0
    for _ in range(n_draws):
        pa = A[rng.choice(ua)]
        pb = B[rng.choice(ub)]
        if pos_match(pa, pb, tol):
            hits += 1
    return hits / n_draws


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tol", type=int, default=50,
                    help="bp tolerance for calling two alignment starts the same position "
                         "(default 50; nanopore start positions jitter)")
    ap.add_argument("--null-draws", type=int, default=200000,
                    help="random cross-library pairs used to estimate the chance rate")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out-tsv")
    args = ap.parse_args()
    rng = random.Random(args.seed)

    samples = load_samples()
    print(f"=== loading molecules (UMI={UMI_LEN}nt, primary mapped, "
          f"starts within {args.tol}bp merged) ===")
    mols, meta = {}, {}
    for s in samples:
        m = load_molecules(s["sample_key"], args.tol)
        if m:
            mols[s["sample_key"]] = m
            meta[s["sample_key"]] = s
    if len(mols) < 2:
        sys.exit("need at least two libraries with hbv.umi.bam")

    rows = []
    print(f"\n=== pairwise UMI sharing (position tolerance {args.tol} bp) ===\n")
    hdr = (f"  {'library A':<26}{'library B':<26}{'shared UMI':>11}{'+pos':>8}"
           f"{'null':>8}{'excess':>9}{'% of smaller':>14}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    for ka, kb in itertools.combinations(mols, 2):
        A, B = mols[ka], mols[kb]
        nA = sum(len(v) for v in A.values())
        nB = sum(len(v) for v in B.values())
        shared = set(A) & set(B)
        n_shared = len(shared)
        n_pos = sum(1 for u in shared if pos_match(A[u], B[u], args.tol))
        p_obs = n_pos / n_shared if n_shared else 0.0
        p_null = null_pos_match_rate(A, B, args.tol, args.null_draws, rng)
        excess = max(0.0, (p_obs - p_null)) * n_shared
        smaller = min(nA, nB)
        pct = 100 * excess / smaller if smaller else 0.0
        print(f"  {ka:<26}{kb:<26}{n_shared:>11,}{n_pos:>8,}"
              f"{p_null*n_shared:>8.1f}{excess:>9.1f}{pct:>13.2f}%")
        rows.append({
            "library_a": ka, "library_b": kb,
            "input_a": meta[ka]["input_ng"], "input_b": meta[kb]["input_ng"],
            "molecules_a": nA, "molecules_b": nB,
            "shared_umis": n_shared, "shared_umi_and_position": n_pos,
            "expected_by_chance": round(p_null * n_shared, 2),
            "excess_molecules": round(excess, 1),
            "pct_of_smaller_library": round(pct, 3),
        })

    if args.out_tsv:
        import csv
        with open(args.out_tsv, "w", newline="") as fh:
            fh.write("# EXP26000892 cross-library UMI sharing; "
                     f"tol={args.tol}bp; umi_len={UMI_LEN}\n")
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader(); w.writerows(rows)
        print(f"\nWritten: {args.out_tsv}")

    print("""
How to read this
  'shared UMI'   UMIs present in both libraries — mostly chance collisions.
  '+pos'         of those, how many ALSO align to the same position (within tol).
  'null'         how many would coincide by chance, from random cross-library pairs.
  'excess'       +pos minus null = molecules plausibly counted in both libraries.
  '% of smaller' excess as a share of the smaller library's molecules — the number
                 that matters, since leakage from an abundant library can dominate
                 a sparse one.

  <1%   negligible; the 1 ng / 0.1 ng floor is real assay behaviour.
  5-20% material; low-input molecule counts are inflated and need correcting.
  >20%  the 0.1 ng result is largely carry-over from the 100 ng library and should
        not be reported as a detection limit.

  Note this measures only leakage that survived strict both-ends demultiplexing —
  i.e. molecules carrying a CONSISTENT wrong barcode at both ends, which is what an
  early-PCR template switch produces. Reads with two different barcodes were already
  removed at demux (4-10% of reads) and are not counted again here.
""")


if __name__ == "__main__":
    main()
