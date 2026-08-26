#!/usr/bin/env python3
"""
Cross-library UMI sharing, v2 — with a negative control and a corrected null.

  # within-run comparison, all pairs (as v1 did)
  python3 cross_library_umi_check_v2.py --root /data/EXP26000993

  # NEGATIVE CONTROL: libraries that were never in the same capture
  python3 cross_library_umi_check_v2.py \
      --root /data/EXP26000993 --root-b /data/EXP26000896 \
      --exclude-b SeqLib5552_1ng_polyA_21 \
      --out-tsv negative_control.tsv

  # tolerance sensitivity
  python3 cross_library_umi_check_v2.py --root /data/EXP26000993 --tol 10


WHY v2 EXISTS — the flaw in v1's null
v1 estimated the chance position-match rate by drawing RANDOM cross-library molecule
pairs. But shared UMIs are not a random subset of molecules: a UMI can only match across
libraries if it was sequenced accurately in BOTH, which biases the shared set toward
high-read-count molecules. High-read molecules cluster at abundance hotspots — exactly
where position matching is most likely. So v1's null samples uniformly while the observed
set is enriched for the easiest-to-match positions, inflating the apparent excess.

THREE CORRECTIONS
  1. NEGATIVE CONTROL (--root-b). Libraries from two different sequencing runs were never
     in the same hybridisation capture and CANNOT have exchanged molecules. Any excess
     measured between them is method artifact, not leakage. This is the decisive
     specificity test: if the negative control is near zero, within-run findings stand.
     NOTE: exclude any library physically shared between the runs (--exclude-b), e.g.
     SeqLib5552, which appears in both EXP26000892 and EXP26000993 and legitimately
     shares molecules.

  2. PERMUTATION NULL. UMIs are shuffled WITHIN each library while each molecule keeps
     its read count and position. This destroys real cross-library identity but preserves
     the abundance-position structure that biases the random-pair null. Repeated
     --permutations times; the mean is the expectation and the spread gives a z-score.

  3. ABUNDANCE STRATIFICATION. Reports median read count of shared vs unshared molecules.
     If shared molecules are systematically higher-abundance, the v1 bias was operating
     and its magnitude is visible here.

Molecules are (UMI, position) after collapsing starts within --tol bp; distinct molecules
can share coordinates, so counts are conservative.
"""
import argparse, itertools, random, sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pysam

UMI_LEN = 28


def load_samples(root):
    p = Path(root) / "config" / "samples.tsv"
    rows = []
    for line in p.read_text().splitlines():
        if line.startswith("#") or line.startswith("barcode\t") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) >= 5:
            rows.append({"barcode": f[0], "sample_key": f[3], "input_ng": f[4],
                         "root": str(root)})
    return rows


def collapse(positions, tol):
    out = []
    for p in sorted(positions):
        if not out or p - out[-1] > tol:
            out.append(p)
    return out


def load_molecules(root, sample_key, tol):
    """-> {umi: {pos: read_count}} — read counts retained for the permutation null
    and for abundance stratification."""
    bam = Path(root) / "analysis" / "samples" / sample_key / "hbv.umi.bam"
    if not bam.exists():
        sys.stderr.write(f"  SKIP {sample_key}: no hbv.umi.bam under {root}\n")
        return None
    raw = defaultdict(lambda: defaultdict(int))
    n_reads = 0
    with pysam.AlignmentFile(bam) as af:
        for r in af.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            if not r.has_tag("RX"):
                continue
            u = r.get_tag("RX")
            if len(u) != UMI_LEN:
                continue
            n_reads += 1
            raw[u][r.reference_start] += 1
    mols = {}
    for u, posdict in raw.items():
        merged = collapse(posdict.keys(), tol)
        acc = {}
        for m in merged:
            acc[m] = sum(c for p, c in posdict.items() if abs(p - m) <= tol)
        mols[u] = acc
    n_mol = sum(len(v) for v in mols.values())
    print(f"  {Path(root).name}/{sample_key:<30} {n_reads:>10,} reads -> "
          f"{n_mol:>8,} molecules ({len(mols):,} UMIs)")
    return mols


def pos_match(pa, pb, tol):
    for a in pa:
        for b in pb:
            if abs(a - b) <= tol:
                return True
    return False


def count_shared(A, B, tol):
    """-> (n_shared_umis, n_position_matched, shared_read_counts_A)"""
    shared = set(A) & set(B)
    matched, abund = 0, []
    for u in shared:
        if pos_match(A[u].keys(), B[u].keys(), tol):
            matched += 1
            abund.append(max(A[u].values()))
    return len(shared), matched, abund


def permuted_expectation(A, B, tol, n_perm, rng):
    """Shuffle UMIs within each library, preserving each molecule's positions and read
    counts. Destroys real cross-library identity; keeps abundance-position structure."""
    ua, ub = list(A.keys()), list(B.keys())
    va, vb = [A[u] for u in ua], [B[u] for u in ub]
    out = []
    for _ in range(n_perm):
        pa = ua[:]; pb = ub[:]
        rng.shuffle(pa); rng.shuffle(pb)
        Ap = dict(zip(pa, va))
        Bp = dict(zip(pb, vb))
        _, m, _ = count_shared(Ap, Bp, tol)
        out.append(m)
    return float(np.mean(out)), float(np.std(out))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, help="project root (experiment A)")
    ap.add_argument("--root-b", help="second project root — enables NEGATIVE CONTROL mode: "
                                     "every A library vs every B library")
    ap.add_argument("--exclude-b", action="append", default=[],
                    help="sample_key in root-b to skip (libraries physically shared "
                         "between the two runs, e.g. SeqLib5552)")
    ap.add_argument("--tol", type=int, default=50)
    ap.add_argument("--permutations", type=int, default=20)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out-tsv")
    args = ap.parse_args()
    rng = random.Random(args.seed)

    print(f"=== cross-library UMI sharing v2 (tol={args.tol} bp, "
          f"{args.permutations} permutations) ===")
    mode = "NEGATIVE CONTROL" if args.root_b else "within-run"
    print(f"    mode: {mode}\n")

    A_samples = load_samples(args.root)
    A = {}
    for s in A_samples:
        m = load_molecules(args.root, s["sample_key"], args.tol)
        if m:
            A[s["sample_key"]] = m

    if args.root_b:
        B_samples = [s for s in load_samples(args.root_b)
                     if s["sample_key"] not in set(args.exclude_b)]
        B = {}
        for s in B_samples:
            m = load_molecules(args.root_b, s["sample_key"], args.tol)
            if m:
                B[s["sample_key"]] = m
        pairs = [(ka, kb) for ka in A for kb in B]
        get = lambda k, side: (A if side == 0 else B)[k]
        if args.exclude_b:
            print(f"\n    excluded from root-b: {', '.join(args.exclude_b)}")
    else:
        B = A
        pairs = list(itertools.combinations(A, 2))
        get = lambda k, side: A[k]

    print(f"\n=== {len(pairs)} pairwise comparisons ===\n")
    hdr = (f"  {'library A':<30}{'library B':<30}{'sharedUMI':>10}{'+pos':>8}"
           f"{'perm_exp':>10}{'perm_sd':>8}{'excess':>9}{'z':>7}{'%smaller':>10}")
    print(hdr); print("  " + "-" * (len(hdr) - 2))

    rows = []
    for ka, kb in pairs:
        Am, Bm = get(ka, 0), get(kb, 1)
        nA = sum(len(v) for v in Am.values())
        nB = sum(len(v) for v in Bm.values())
        n_shared, n_pos, abund_shared = count_shared(Am, Bm, args.tol)
        exp_m, sd_m = permuted_expectation(Am, Bm, args.tol, args.permutations, rng)
        excess = max(0.0, n_pos - exp_m)
        z = (n_pos - exp_m) / sd_m if sd_m > 0 else float("nan")
        smaller = min(nA, nB)
        pct = 100 * excess / smaller if smaller else 0.0

        all_abund = [max(v.values()) for m in (Am,) for v in m.values()]
        med_shared = float(np.median(abund_shared)) if abund_shared else float("nan")
        med_all = float(np.median(all_abund)) if all_abund else float("nan")

        print(f"  {ka[:29]:<30}{kb[:29]:<30}{n_shared:>10,}{n_pos:>8,}"
              f"{exp_m:>10.1f}{sd_m:>8.1f}{excess:>9.1f}{z:>7.1f}{pct:>9.2f}%")
        rows.append({
            "library_a": ka, "library_b": kb, "mode": mode,
            "molecules_a": nA, "molecules_b": nB,
            "shared_umis": n_shared, "shared_umi_and_position": n_pos,
            "permuted_expectation": round(exp_m, 2), "permuted_sd": round(sd_m, 2),
            "excess_molecules": round(excess, 1), "z_score": round(z, 2),
            "pct_of_smaller_library": round(pct, 3),
            "median_reads_shared_molecules": med_shared,
            "median_reads_all_molecules_a": med_all,
            "tol_bp": args.tol,
        })

    if args.out_tsv and rows:
        import csv
        with open(args.out_tsv, "w", newline="") as fh:
            fh.write(f"# cross-library UMI sharing v2; mode={mode}; tol={args.tol}; "
                     f"permutations={args.permutations}\n")
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader(); w.writerows(rows)
        print(f"\nWritten: {args.out_tsv}")

    print(f"""
How to read this

  perm_exp / perm_sd   position-matched sharing expected when UMIs are shuffled within
                       each library, preserving read count and position. Unlike v1's
                       random-pair null this keeps the abundance-position structure that
                       makes shared (high-abundance) molecules easier to match.
  z                    (observed - expected) / sd. Above ~5 is a real excess; near 0 is
                       fully explained by chance plus structure.
  median_reads_shared  vs median_reads_all_molecules_a. If shared molecules are markedly
                       higher-abundance, v1's null was biased and its excess overstated.

  NEGATIVE CONTROL mode: these libraries were never in the same capture, so any excess is
  method artifact. Near-zero excess and z near 0 validates the within-run measurements.
  A large excess here means the method is detecting something other than leakage and no
  within-run conclusion should be drawn from it.
""")


if __name__ == "__main__":
    main()
