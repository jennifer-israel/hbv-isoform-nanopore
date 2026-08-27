#!/usr/bin/env python3
"""
HBV transcript classifier.

COORDINATES
-----------
Alignments are against U95551.1_2x — the 3,182 bp HBV genome doubled to 6,364 bp so
that genome-wrapping transcripts align as one continuous record. Coordinates are KEPT
in 0-6363 space; `mod 3182` is used ONLY to derive the genomic TSS / 3' end for
window tests, never to rewrite stored coordinates.

Input reads are pychopper-oriented, so for a forward-strand alignment `reference_start`
is the transcript 5' end (TSS) and `reference_end` is the 3' end.

RULES (as supplied)
-------------------
Transcripts
  pgRNA_RT  genomic footprint >= 3,982 bp (3182 + 800): tandem/concatemeric
            readthrough spanning ~>= 1.25 HBV genomes
  spliced   at least two exon blocks remain after merging gaps <= 200 bp, i.e. at
            least one alignment gap > 200 bp
  precore   TSS 1,730-1,815 inclusive, crosses the 3,182 junction, span >= 2,600 bp
  pgRNA     TSS 1,816-1,880 inclusive, crosses the 3,182 junction, span >= 2,600 bp
  preS1     TSS >= 2,700 and < 3,100, span >= 1,600 bp
  preS2_S   circular window: TSS >= 3,100 or <= 150, span >= 1,200 bp
  X         TSS 1,260-1,450 inclusive, span 300-1,000 bp inclusive

Splicing
  All reads with an intron > 200 bp are assigned `spliced` EXCEPT those already
  assigned `pgRNA_RT`. The largest detected intron is annotated in `splice_junction`:
     SP1              donor within +/-60 bp of 2,447 (2,387-2,507)
                      AND acceptor within +/-60 bp of 489 (429-549)
     non_canonical    otherwise

Poly-A status (reported separately; does NOT determine transcript class)
  canonical    3' end within 1,936 +/- 80  -> 1,856-2,016
  alt_polya    3' end within 1,808 +/- 25  -> 1,783-1,833
  truncated    ends before the canonical window and not in the alternative window
  readthrough  ends after the canonical window
  NA           antisense read

Antisense reads (minus strand after pychopper orientation) are flagged and not
TSS-binned.

INTERFACE
---------
classify(reference_start, reference_end, blocks, is_reverse) -> dict with keys:
    tx_class, tss, polya_end, polya_status, blocks, splice_junction,
    n_exon_blocks, max_intron
`blocks` is the pysam AlignedSegment.get_blocks() list of (start, end) tuples.
"""

VERSION = "hbv-transcript-classify-REIMPL-v1"

HBV_LEN = 3182          # single-copy genome length; 2x contig = 6364
MIN_INTRON = 200        # gaps <= this are merged (minimap2 splice-mode artifacts)
JUNCTION = 3182         # linearisation point / copy boundary in 2x space

# --- transcript class thresholds ---
PGRNA_RT_MIN_SPAN = 3982            # 3182 + 800
PRECORE_TSS = (1730, 1815)          # inclusive
PGRNA_TSS = (1816, 1880)            # inclusive
WRAP_MIN_SPAN = 2600                # precore + pgRNA
PRES1_TSS_LO, PRES1_TSS_HI = 2700, 3100   # >= lo and < hi
PRES1_MIN_SPAN = 1600
PRES2S_TSS_HI, PRES2S_TSS_LO = 3100, 150  # >= 3100 or <= 150 (circular)
PRES2S_MIN_SPAN = 1200
X_TSS = (1260, 1450)                # inclusive
X_SPAN = (300, 1000)                # inclusive

# --- splice junction (SP1) windows, single-copy coords ---
SP1_DONOR, SP1_ACCEPTOR, SP1_TOL = 2447, 489, 60
SP1_DONOR_WIN = (SP1_DONOR - SP1_TOL, SP1_DONOR + SP1_TOL)        # 2387-2507
SP1_ACCEPTOR_WIN = (SP1_ACCEPTOR - SP1_TOL, SP1_ACCEPTOR + SP1_TOL)  # 429-549

# --- poly-A windows, single-copy coords ---
POLYA_CANONICAL = (1856, 2016)      # 1936 +/- 80
POLYA_ALT = (1783, 1833)            # 1808 +/- 25


def merge_blocks(blocks, max_gap=MIN_INTRON):
    """Merge aligned blocks separated by <= max_gap. Returns list of (start, end).

    minimap2 splice mode emits small spurious gaps (72-141 bp observed); merging at
    200 bp removes those while retaining the ~1.2 kb SP1 intron.
    """
    if not blocks:
        return []
    ordered = sorted(tuple(b) for b in blocks)
    merged = [list(ordered[0])]
    for s, e in ordered[1:]:
        if s - merged[-1][1] <= max_gap:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(int(s), int(e)) for s, e in merged]


def introns_from(merged):
    """[(donor, acceptor, length), ...] for the gaps between merged blocks."""
    out = []
    for i in range(len(merged) - 1):
        donor = merged[i][1]          # first intronic base (2x coords)
        acceptor = merged[i + 1][0]   # first exonic base after the intron
        out.append((donor, acceptor, acceptor - donor))
    return out


def _norm(pos):
    """2x coordinate -> single-copy genomic coordinate."""
    return pos % HBV_LEN


def _in(v, lo, hi):
    return lo <= v <= hi


def classify_splice_junction(merged):
    """Annotate the LARGEST intron as SP1 or non_canonical. '' if unspliced."""
    ints = introns_from(merged)
    if not ints:
        return "", None
    donor, acceptor, length = max(ints, key=lambda x: x[2])
    d, a = _norm(donor), _norm(acceptor)
    if _in(d, *SP1_DONOR_WIN) and _in(a, *SP1_ACCEPTOR_WIN):
        return "SP1", length
    return "non_canonical", length


def polya_status_for(polya_end, is_reverse):
    """Poly-A status. Reported separately; does not affect transcript class."""
    if is_reverse:
        return "NA"
    if _in(polya_end, *POLYA_CANONICAL):
        return "canonical"
    if _in(polya_end, *POLYA_ALT):
        return "alt_polya"
    if polya_end > POLYA_CANONICAL[1]:
        return "readthrough"
    return "truncated"


def classify(reference_start, reference_end, blocks, is_reverse):
    """Classify one alignment. See module docstring for the rule set."""
    rs, re_ = int(reference_start), int(reference_end)
    span = re_ - rs

    merged = merge_blocks(blocks)
    if not merged:
        merged = [(rs, re_)]
    sj, max_intron = classify_splice_junction(merged)
    n_blocks = len(merged)

    tss = _norm(rs)
    polya_end = (re_ - 1) % HBV_LEN + 1 if re_ else 0
    pa_status = polya_status_for(polya_end, is_reverse)

    base = {
        "tss": tss,
        "polya_end": polya_end,
        "polya_status": pa_status,
        "blocks": ";".join(f"{s}-{e}" for s, e in merged),
        "splice_junction": sj,
        "n_exon_blocks": n_blocks,
        "max_intron": max_intron if max_intron is not None else 0,
    }

    # 1. antisense — flagged, not TSS-binned
    if is_reverse:
        return {**base, "tx_class": "antisense"}

    # 2. pgRNA_RT — takes precedence over spliced, per spec
    if span >= PGRNA_RT_MIN_SPAN:
        return {**base, "tx_class": "pgRNA_RT"}

    # 3. spliced — any real intron (>200 bp), regardless of TSS
    if n_blocks >= 2:
        return {**base, "tx_class": "spliced"}

    # 4. TSS + span gated classes
    crosses = rs < JUNCTION <= re_

    if crosses and span >= WRAP_MIN_SPAN:
        if _in(tss, *PRECORE_TSS):
            return {**base, "tx_class": "precore"}
        if _in(tss, *PGRNA_TSS):
            return {**base, "tx_class": "pgRNA"}

    if PRES1_TSS_LO <= tss < PRES1_TSS_HI and span >= PRES1_MIN_SPAN:
        return {**base, "tx_class": "preS1"}

    if (tss >= PRES2S_TSS_HI or tss <= PRES2S_TSS_LO) and span >= PRES2S_MIN_SPAN:
        return {**base, "tx_class": "preS2_S"}

    if _in(tss, *X_TSS) and _in(span, *X_SPAN):
        return {**base, "tx_class": "X"}

    # 5. everything else
    return {**base, "tx_class": "unclassified"}


# ---------------------------------------------------------------------------
# self-test: python hbv_transcript_classify.py
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    def blk(*pairs):
        return [tuple(p) for p in pairs]

    cases = [
        # (label, start, end, blocks, is_reverse, expected_class)
        ("antisense",            1750, 4500, blk((1750, 4500)), True,  "antisense"),
        ("pgRNA_RT span 4000",   1000, 5000, blk((1000, 5000)), False, "pgRNA_RT"),
        ("pgRNA_RT at exactly 3982", 0, 3982, blk((0, 3982)),   False, "pgRNA_RT"),
        ("spliced SP1",          1750, 3600, blk((1750, 2447), (3671, 3600 + 500)), False, "spliced"),
        ("precore TSS1730",      1730, 4400, blk((1730, 4400)), False, "precore"),
        ("precore TSS1815",      1815, 4500, blk((1815, 4500)), False, "precore"),
        ("pgRNA TSS1816",        1816, 4500, blk((1816, 4500)), False, "pgRNA"),
        ("pgRNA TSS1880",        1880, 4600, blk((1880, 4600)), False, "pgRNA"),
        ("precore too short",    1750, 3300, blk((1750, 3300)), False, "unclassified"),
        ("precore no junction",  1750, 3100, blk((1750, 3100)), False, "unclassified"),
        ("preS1 TSS2700",        2700, 4400, blk((2700, 4400)), False, "preS1"),
        ("preS1 TSS3099",        3099, 4800, blk((3099, 4800)), False, "preS1"),
        ("preS1 short -> uncl",  2700, 3500, blk((2700, 3500)), False, "unclassified"),
        ("preS2_S TSS3100",      3100, 4400, blk((3100, 4400)), False, "preS2_S"),
        ("preS2_S TSS150",        150, 1400, blk((150, 1400)),  False, "preS2_S"),
        ("X TSS1260 span300",    1260, 1560, blk((1260, 1560)), False, "X"),
        ("X TSS1450 span1000",   1450, 2450, blk((1450, 2450)), False, "X"),
        ("X span 299 -> uncl",   1260, 1559, blk((1260, 1559)), False, "unclassified"),
        ("X span 1001 -> uncl",  1260, 2261, blk((1260, 2261)), False, "unclassified"),
        ("small gap merged",     3100, 4400, blk((3100, 3700), (3800, 4400)), False, "preS2_S"),
    ]
    print(f"{VERSION}\n")
    fails = 0
    for label, s, e, b, rev, want in cases:
        got = classify(s, e, b, rev)["tx_class"]
        ok = got == want
        fails += (not ok)
        print(f"  [{'PASS' if ok else 'FAIL'}] {label:<28} -> {got:<14} (expected {want})")

    print("\n  splice junction annotation:")
    sp1 = classify(1750, 3600, blk((1750, 2447), (3671, 4100)), False)
    print(f"    SP1 case            splice_junction={sp1['splice_junction']!r} "
          f"max_intron={sp1['max_intron']}")
    fails += (sp1["splice_junction"] != "SP1")
    nc = classify(1750, 3600, blk((1750, 2100), (3400, 4100)), False)
    print(f"    non-canonical case  splice_junction={nc['splice_junction']!r} "
          f"max_intron={nc['max_intron']}")
    fails += (nc["splice_junction"] != "non_canonical")

    print("\n  poly-A status (single-copy 3' end):")
    for end2x, want in [(1936, "canonical"), (1856, "canonical"), (2016, "canonical"),
                        (1808, "alt_polya"), (1783, "alt_polya"), (1833, "alt_polya"),
                        (1700, "truncated"), (2500, "readthrough")]:
        got = polya_status_for(end2x, False)
        ok = got == want
        fails += (not ok)
        print(f"    [{'PASS' if ok else 'FAIL'}] 3' end {end2x:>5} -> {got:<12} (expected {want})")
    print(f"    [{'PASS' if polya_status_for(1936, True) == 'NA' else 'FAIL'}] antisense -> NA")

    print(f"\n  {'ALL PASSED' if fails == 0 else str(fails) + ' FAILURE(S)'}")
    raise SystemExit(1 if fails else 0)
