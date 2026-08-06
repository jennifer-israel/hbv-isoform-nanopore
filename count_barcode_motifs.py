#!/usr/bin/env python3
"""
TEST 2 — count barcode motif occurrences across the WHOLE read.

Question this answers: are the long unclassified reads concatemers/chimeras, or
real long molecules whose barcode simply basecalled poorly?

The two hypotheses make opposite predictions:

  concatemer / chimera   -> 3+ barcode motifs per read, and/or 2 motifs that are
                            DIFFERENT barcodes, and/or motifs far from both ends.
                            Multi-motif fraction RISES with read length.
  real long molecule     -> at most 2 motifs, both the SAME barcode, both near the
                            ends. Multi-motif fraction FLAT with read length.

Unlike demux_custom_barcodes.py (which only inspects a window at each end), this
scans the entire read and reports EVERY occurrence, by iteratively finding the best
match and masking it before searching again.

It also looks for internal P5 / P7 adapter sequence, since a genuine fusion of two
library molecules necessarily contains an adapter junction mid-read. That is close
to definitive on its own.

Usage:
    pip install edlib
    python count_barcode_motifs.py --fastq unc_sample.fastq.gz --max-reads 200000
    python count_barcode_motifs.py --fastq unc_sample.fastq.gz --out-tsv per_read.tsv

Compare unclassified against a classified barcode as a control:
    python count_barcode_motifs.py --fastq classified_bc01.fastq.gz --max-reads 200000
"""
import argparse, gzip, sys
from collections import Counter, defaultdict

try:
    import edlib
except ImportError:
    sys.exit("ERROR: needs edlib. pip install edlib --break-system-packages")

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
def rc(s): return s.translate(COMPLEMENT)[::-1]

FLANK_UP, FLANK_DOWN = "GGTGCTG", "TTAACCT"
BARCODES = {
    "bc01": "AAGAAAGTTGTCGGTGTCTTTGTG",
    "bc02": "TCGATTCCGTTTGTAGTCGTCTGT",
    "bc03": "GAGTCTTGTGTCCCAGTTACCAGG",
    "bc04": "TTCGGATTCTATCGTGTTTCCCTA",
}
P5 = "AATGATACGGCGACCACCGA"
P7_RC = "CAAGCAGAAGACGGCATACGAGAT"      # as it appears in the rev1 primer

# read-length strata for the key stratified table
STRATA = [(0, 1000), (1000, 2000), (2000, 3000), (3000, 4000), (4000, 10**9)]
def stratum(n):
    for lo, hi in STRATA:
        if lo <= n < hi:
            return f"{lo}-{hi}" if hi < 10**9 else f"{lo}+"
    return "?"


def build_motifs():
    m = {}
    for bc, seq in BARCODES.items():
        fwd = FLANK_UP + seq + FLANK_DOWN
        m[bc] = {"fwd": fwd, "rev": rc(fwd)}
    return m


def find_all(seq, query, max_ed, max_hits=8):
    """All non-overlapping occurrences of `query` in `seq` within max_ed.
    Iteratively takes the best match, records it, masks it, repeats."""
    hits = []
    work = seq
    for _ in range(max_hits):
        r = edlib.align(query, work, mode="HW", task="locations", k=max_ed)
        if r["editDistance"] < 0 or not r.get("locations"):
            break
        s, e = r["locations"][0]
        if s is None:
            break
        hits.append((s, e, r["editDistance"]))
        work = work[:s] + "N" * (e - s + 1) + work[e + 1:]
    return hits


def scan_read(seq, motifs, max_ed, end_margin):
    """Return (all_hits, n_internal_adapter). all_hits = [(pos, barcode, orient, ed)]"""
    hits = []
    for bc, m in motifs.items():
        for orient in ("fwd", "rev"):
            for s, e, ed in find_all(seq, m[orient], max_ed):
                hits.append((s, bc, orient, ed))
    hits.sort()
    # de-duplicate overlapping calls from different barcodes: keep lowest edit distance
    dedup, used = [], []
    for pos, bc, orient, ed in sorted(hits, key=lambda h: h[3]):
        if any(abs(pos - p) < 30 for p, _, _, _ in dedup):
            continue
        dedup.append((pos, bc, orient, ed))
    dedup.sort()

    # internal adapter: P5 or rc(P7) found away from both ends
    n_internal_adapter = 0
    L = len(seq)
    for ad in (P5, P7_RC):
        for a in (ad, rc(ad)):
            for s, e, ed in find_all(seq, a, max(2, len(a) // 6), max_hits=4):
                if end_margin < s < L - end_margin:
                    n_internal_adapter += 1
    return dedup, n_internal_adapter


def iter_fastq(path, limit):
    op = gzip.open if str(path).endswith(".gz") else open
    n = 0
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().strip(); fh.readline(); fh.readline()
            n += 1
            if limit and n > limit:
                break
            yield h[1:].split()[0], s


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--max-reads", type=int, default=100000)
    ap.add_argument("--max-edit-distance", type=int, default=8,
                    help="per-motif edit distance ceiling (38bp motif; default 8, looser "
                         "than demux on purpose — we want to FIND motifs, not assign them)")
    ap.add_argument("--end-margin", type=int, default=200,
                    help="a motif/adapter beyond this distance from both ends is 'internal'")
    ap.add_argument("--out-tsv", help="optional per-read output")
    args = ap.parse_args()

    motifs = build_motifs()
    cat = Counter()
    by_stratum = defaultdict(Counter)
    strat_n = Counter()
    internal_by_stratum = Counter()
    nhits_hist = Counter()

    out = open(args.out_tsv, "w") if args.out_tsv else None
    if out:
        out.write("read_id\tlength\tn_motifs\tbarcodes\tpositions\tn_internal_adapter\tcategory\n")

    total = 0
    for rid, seq in iter_fastq(args.fastq, args.max_reads):
        total += 1
        L = len(seq)
        hits, n_int_ad = scan_read(seq, motifs, args.max_edit_distance, args.end_margin)
        bcs = [h[1] for h in hits]
        uniq = set(bcs)
        n = len(hits)
        nhits_hist[min(n, 5)] += 1

        if n == 0:
            c = "0_motifs"
        elif n == 1:
            c = "1_motif"
        elif n == 2 and len(uniq) == 1:
            c = "2_same"          # normal, well-formed molecule
        elif n == 2 and len(uniq) == 2:
            c = "2_DIFFERENT"     # chimera
        else:
            c = f"{min(n,5)}plus_motifs"  # concatemer

        st = stratum(L)
        cat[c] += 1
        by_stratum[st][c] += 1
        strat_n[st] += 1
        if n_int_ad > 0:
            internal_by_stratum[st] += 1

        if out:
            out.write(f"{rid}\t{L}\t{n}\t{','.join(bcs)}\t"
                      f"{','.join(str(h[0]) for h in hits)}\t{n_int_ad}\t{c}\n")
    if out:
        out.close()

    print(f"\n=== {args.fastq} — {total:,} reads scanned ===\n")

    print("Motifs found per read:")
    for k in sorted(nhits_hist):
        lbl = f"{k}" if k < 5 else "5+"
        print(f"  {lbl:>3} motifs : {nhits_hist[k]:>9,}  ({100*nhits_hist[k]/total:5.2f}%)")

    print("\nCategory totals:")
    for c, n in cat.most_common():
        print(f"  {c:<18} {n:>9,}  ({100*n/total:5.2f}%)")

    print("\n*** KEY TABLE — category by read length ***")
    order = ["0_motifs", "1_motif", "2_same", "2_DIFFERENT",
             "3plus_motifs", "4plus_motifs", "5plus_motifs"]
    present = [c for c in order if any(by_stratum[s].get(c) for s in by_stratum)]
    hdr = f"  {'length':<10}{'n':>10}" + "".join(f"{c:>14}" for c in present) + f"{'internal_ad':>13}"
    print(hdr)
    for lo, hi in STRATA:
        st = f"{lo}-{hi}" if hi < 10**9 else f"{lo}+"
        n = strat_n.get(st, 0)
        if not n:
            continue
        row = f"  {st:<10}{n:>10,}"
        for c in present:
            row += f"{100*by_stratum[st].get(c,0)/n:>13.1f}%"
        row += f"{100*internal_by_stratum.get(st,0)/n:>12.1f}%"
        print(row)

    print("\nHow to read it:")
    print("  CONCATEMER/CHIMERA  -> '3plus_motifs', '2_DIFFERENT' and 'internal_ad' all")
    print("                         RISE sharply with read length. Long reads are fusions,")
    print("                         and strict demux was right to reject them.")
    print("  REAL LONG READS     -> those columns stay FLAT across strata and '2_same' or")
    print("                         '1_motif' dominates even at 3kb+. Then the long tail is")
    print("                         genuine molecules you are losing, and the demux")
    print("                         thresholds (not the biology) need revisiting.")
    print("  '0_motifs' high     -> reads with no detectable barcode anywhere: neither")
    print("                         hypothesis; likely adapter-less or badly degraded ends.")


if __name__ == "__main__":
    main()
