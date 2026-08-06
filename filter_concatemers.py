#!/usr/bin/env python3
"""
Drop homologous concatemers from a FASTQ stream, before pychopper.

  samtools fastq -n in.bam | python3 filter_concatemers.py --stats s.tsv | pychopper ...

WHY THIS EXISTS
`dorado demux --barcode-both-ends` rejects a fusion of two molecules from DIFFERENT
libraries, because the two ends disagree on the barcode. A fusion of two molecules from
the SAME library carries a consistent barcode at both outer ends, agrees, and is assigned
normally. Nothing downstream removes it. Those reads then reach transcript
classification, where on a 2x reference a fused pair of HBV molecules can align
contiguously across the copy junction and be called `pgRNA_RT` (span >= 3,982 bp,
"tandem/concatemeric readthrough") — a PCR artifact scored as biology.

THE TEST
An adapter or barcode motif found far from BOTH read ends cannot occur within a single
intact molecule; the library construct places them only at the termini. Its presence
mid-read is a fusion junction. This is the same criterion used to characterise the
unclassified bucket, applied here to reads that survived demultiplexing.

WHERE IT MUST GO
Immediately before pychopper. pychopper TRIMS the primers, so its output no longer
contains the adapter sequence — scanning downstream of it would find nothing whether or
not the read was a concatemer. Placed in the phase 2 pipe it only ever sees HBV-aligning
reads (hundreds of thousands, not the full run), so it is cheap.

WHAT IT DOES NOT DO
It cannot detect a fusion whose junction happens to fall within END_MARGIN of a read end,
nor one where the internal adapter was degraded beyond the edit-distance threshold. The
reported rate is therefore a floor.

  --end-margin      a motif beyond this distance from both ends is internal (default 200)
  --max-edit-distance   per-motif tolerance; adapters use len/6 (default 8 for barcodes)
  --stats           TSV of counts written here (also to stderr)
  --keep-removed    write the discarded reads to this FASTQ for inspection
"""
import argparse, sys

try:
    import edlib
except ImportError:
    sys.exit("ERROR: needs edlib (pip install edlib)")

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
def rc(s): return s.translate(COMPLEMENT)[::-1]

FLANK_UP, FLANK_DOWN = "GGTGCTG", "TTAACCT"
P5 = "AATGATACGGCGACCACCGA"
P7_RC = "CAAGCAGAAGACGGCATACGAGAT"

DEFAULT_BARCODES = [
    "AAGAAAGTTGTCGGTGTCTTTGTG",
    "TCGATTCCGTTTGTAGTCGTCTGT",
    "GAGTCTTGTGTCCCAGTTACCAGG",
    "TTCGGATTCTATCGTGTTTCCCTA",
]


def load_barcodes(path):
    if not path:
        return DEFAULT_BARCODES
    seqs, cur = [], []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur)); cur = []
        elif line:
            cur.append(line)
    if cur:
        seqs.append("".join(cur))
    return seqs or DEFAULT_BARCODES


def build_probes(barcodes):
    """(sequence, max_edit_distance) pairs to search for, both orientations."""
    probes = []
    for b in barcodes:
        m = FLANK_UP + b + FLANK_DOWN
        probes += [m, rc(m)]
    probes += [P5, rc(P5), P7_RC, rc(P7_RC)]
    return probes


def has_internal(seq, probes, end_margin, bc_ed):
    """True if any probe matches at a position beyond end_margin from both ends."""
    n = len(seq)
    if n < 2 * end_margin + 40:
        return False
    # only the interior can contain a junction; searching it directly is cheaper
    # and removes any chance of a terminal match being miscounted
    interior = seq[end_margin:n - end_margin]
    if len(interior) < 30:
        return False
    for p in probes:
        k = bc_ed if len(p) > 30 else max(2, len(p) // 6)
        if edlib.align(p, interior, mode="HW", task="distance", k=k)["editDistance"] >= 0:
            return True
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--barcodes-fasta")
    ap.add_argument("--end-margin", type=int, default=200)
    ap.add_argument("--max-edit-distance", type=int, default=8)
    ap.add_argument("--stats")
    ap.add_argument("--label", default="")
    ap.add_argument("--keep-removed")
    args = ap.parse_args()

    probes = build_probes(load_barcodes(args.barcodes_fasta))
    removed_fh = open(args.keep_removed, "w") if args.keep_removed else None

    kept = removed = 0
    out = sys.stdout
    r = sys.stdin
    while True:
        h = r.readline()
        if not h:
            break
        s = r.readline(); p = r.readline(); q = r.readline()
        if not q:
            break
        if has_internal(s.rstrip("\n"), probes, args.end_margin, args.max_edit_distance):
            removed += 1
            if removed_fh:
                removed_fh.write(h + s + p + q)
        else:
            kept += 1
            out.write(h + s + p + q)
    out.flush()
    if removed_fh:
        removed_fh.close()

    total = kept + removed
    pct = (100 * removed / total) if total else 0.0
    msg = (f"[filter_concatemers] {args.label} total={total} kept={kept} "
           f"removed={removed} ({pct:.2f}%)")
    sys.stderr.write(msg + "\n")
    if args.stats:
        with open(args.stats, "w") as fh:
            fh.write("label\ttotal_reads\tkept\tremoved\tremoved_pct\tend_margin\n")
            fh.write(f"{args.label}\t{total}\t{kept}\t{removed}\t{pct:.4f}\t{args.end_margin}\n")


if __name__ == "__main__":
    main()
