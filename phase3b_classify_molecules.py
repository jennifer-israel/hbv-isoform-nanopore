#!/usr/bin/env python3
"""
Phase 3b — transcript classification on DEDUPLICATED molecules, plus a per-class
duplication diagnostic.  (EXP26000892)

WHY THIS EXISTS
phase3_classify.py classifies hbv.umi.bam — every UMI-tagged READ, undeduplicated.
That was a deliberate choice in EXP26000559 (recovery reported in molecules,
composition in reads) and is sound under one assumption: that PCR duplication is
uniform across transcript classes.

This experiment likely violates that assumption. Several classes are gated on LENGTH
(X 300-1000 bp, preS2_S >=1200, precore/pgRNA >=2600), short amplicons duplicate
preferentially, and duplication here runs 92-99.7% with severe jackpotting at 0.1 ng.
Read-level composition gave X = 75% of bc04, which is what amplification bias alone
would produce.

WHAT IT DOES
  1. Regenerates the deduplicated BAM per library (umi_tools directional on
     well-formed 28-nt UMIs) and KEEPS it.
  2. Classifies those molecules with the shared hbv_transcript_classify module.
  3. Emits, per library:
        hbv_classified_molecules.tsv / .parquet   one row per molecule
  4. Emits three comparison tables under analysis/comparison/:
        phase3b_molecule_counts.tsv     class counts + proportions, molecule level
        phase3b_read_vs_molecule.tsv    read % vs molecule % side by side
        phase3b_class_duplication.tsv   reads per molecule WITHIN each class
                                        -- the diagnostic that says whether the
                                        read-level composition was trustworthy

READING phase3b_class_duplication.tsv
  If reads_per_molecule is similar across classes, duplication is uniform, the
  read-based composition is unbiased, and EXP26000559's approach carries over.
  If short-span classes (X) show markedly higher reads_per_molecule than long ones
  (precore/pgRNA/preS1), duplication is length-biased and only the molecule-level
  composition should be interpreted.

  conda activate hbv_lr
  python3 phase3b_classify_molecules.py                       # all libraries
  python3 phase3b_classify_molecules.py --sample-key <key>    # subset
  python3 phase3b_classify_molecules.py --keep-dedup-bam      # retain the BAMs
"""
import argparse, datetime, subprocess, sys, shutil
from pathlib import Path

import pysam
import pandas as pd

PROJECT_ROOT = Path("/data/EXP26000993")
ANALYSIS = PROJECT_ROOT / "analysis"
SAMPLES = PROJECT_ROOT / "config" / "samples.tsv"
SHARED = Path("/data/shared")
VERSION = "phase3b-v1"

sys.path.insert(0, str(SHARED))
try:
    from hbv_transcript_classify import classify, VERSION as CLASSIFY_VERSION
except ImportError:
    sys.exit(f"ERROR: hbv_transcript_classify.py not found in {SHARED}")

CLASS_ORDER = ["preS2_S", "preS1", "precore", "pgRNA", "pgRNA_RT", "X", "spliced",
               "unclassified", "antisense"]


def load_samples(keys=None):
    rows = []
    for line in SAMPLES.read_text().splitlines():
        if line.startswith("#") or line.startswith("barcode\t") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) >= 5 and (not keys or f[3] in keys):
            rows.append({"barcode": f[0], "lib_id": f[1], "sample_key": f[3],
                         "input_ng": f[4]})
    return rows


def make_dedup_bam(d: Path, threads=4):
    """hbv.umi.bam -> hbv.umi.dedup.bam (one representative read per molecule)."""
    src, umi28, ded = d / "hbv.umi.bam", d / "hbv.umi28.bam", d / "hbv.umi.dedup.bam"
    if not src.exists():
        return None
    if ded.exists():
        return ded
    subprocess.run(f"samtools view -b -e 'length([RX])==28' {src} "
                   f"| samtools sort -@ {threads} -o {umi28} -",
                   shell=True, check=True)
    pysam.index(str(umi28))
    r = subprocess.run(["umi_tools", "dedup", "-I", str(umi28), "-S", str(ded),
                        "--extract-umi-method=tag", "--umi-tag=RX",
                        "--method=directional", "-L", str(d / "hbv_umi_dedup.log")],
                       capture_output=True, text=True)
    if r.returncode != 0 or not ded.exists():
        sys.stderr.write(f"  ERROR umi_tools dedup failed in {d.name}:\n{r.stderr[:800]}\n")
        return None
    pysam.index(str(ded))
    umi28.unlink(missing_ok=True)
    Path(str(umi28) + ".bai").unlink(missing_ok=True)
    return ded


def classify_bam(bam: Path, sample_key, barcode):
    rows = []
    with pysam.AlignmentFile(bam) as af:
        for r in af.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            c = classify(r.reference_start, r.reference_end, r.get_blocks(), r.is_reverse)
            rows.append({
                "sample_key": sample_key, "barcode": barcode,
                "read_id": r.query_name,
                "umi": r.get_tag("RX") if r.has_tag("RX") else "",
                "strand": "-" if r.is_reverse else "+",
                "mapq": r.mapping_quality,
                "read_length": r.query_length or 0,
                "ref_start": r.reference_start, "ref_end": r.reference_end,
                "span_len": r.reference_end - r.reference_start,
                **c,
            })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-key", action="append")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--keep-dedup-bam", action="store_true")
    args = ap.parse_args()

    for t in ("samtools", "umi_tools"):
        if shutil.which(t) is None:
            sys.exit(f"ERROR: {t} not on PATH")

    today = datetime.date.today().isoformat()
    samples = load_samples(set(args.sample_key) if args.sample_key else None)
    comp = ANALYSIS / "comparison"; comp.mkdir(parents=True, exist_ok=True)

    print(f"=== Phase 3b: classification on deduplicated molecules ===")
    print(f"    classifier: {CLASSIFY_VERSION}\n")

    mol_frames, read_frames = [], []
    for s in samples:
        sk, bc = s["sample_key"], s["barcode"]
        d = ANALYSIS / "samples" / sk
        print(f"--- {sk} ({bc}, {s['input_ng']} ng)")

        ded = make_dedup_bam(d, args.threads)
        if ded is None:
            print("    skipped (no dedup BAM)"); continue

        dfm = classify_bam(ded, sk, bc)
        if dfm.empty:
            print("    0 molecules"); continue
        hdr = (f"# sample={sk}; pipeline_version={VERSION}; "
               f"classifier={CLASSIFY_VERSION}; level=MOLECULE; generated={today}")
        with open(d / "hbv_classified_molecules.tsv", "w") as fh:
            fh.write(hdr + "\n"); dfm.to_csv(fh, sep="\t", index=False)
        dfm.to_parquet(d / "hbv_classified_molecules.parquet", index=False)
        mol_frames.append(dfm)

        # read-level counts, for the comparison — reuse Phase 3 output if present
        rp = d / "hbv_classified.parquet"
        if rp.exists():
            read_frames.append(pd.read_parquet(rp).assign(sample_key=sk))

        vc = dfm.tx_class.value_counts()
        print(f"    {len(dfm):,} molecules — "
              + ", ".join(f"{k}={v}" for k, v in vc.items()))
        if not args.keep_dedup_bam:
            ded.unlink(missing_ok=True)
            Path(str(ded) + ".bai").unlink(missing_ok=True)

    if not mol_frames:
        sys.exit("no molecules classified")
    mols = pd.concat(mol_frames, ignore_index=True)

    # ---- molecule-level counts + proportions ----
    rows = []
    for sk, g in mols.groupby("sample_key"):
        tot = len(g)
        for cls in CLASS_ORDER:
            n = int((g.tx_class == cls).sum())
            rows.append({"sample_key": sk, "tx_class": cls, "n_molecules": n,
                         "total_molecules": tot,
                         "proportion": round(n / tot, 5) if tot else 0})
    mc = pd.DataFrame(rows)
    _write(mc, comp / "phase3b_molecule_counts.tsv", VERSION, today)

    print("\n=== molecule-level composition (% of classified molecules) ===")
    piv = mc.pivot_table(index="sample_key", columns="tx_class", values="proportion",
                         fill_value=0) * 100
    piv = piv[[c for c in CLASS_ORDER if c in piv.columns]]
    print(piv.round(2).to_string())

    # ---- read vs molecule, and per-class duplication ----
    if read_frames:
        reads = pd.concat(read_frames, ignore_index=True)
        cmp_rows, dup_rows = [], []
        for sk in mols.sample_key.unique():
            gm = mols[mols.sample_key == sk]
            gr = reads[reads.sample_key == sk]
            tm, tr = len(gm), len(gr)
            for cls in CLASS_ORDER:
                nm = int((gm.tx_class == cls).sum())
                nr = int((gr.tx_class == cls).sum())
                if nm == 0 and nr == 0:
                    continue
                cmp_rows.append({
                    "sample_key": sk, "tx_class": cls,
                    "reads": nr, "molecules": nm,
                    "pct_of_reads": round(100 * nr / tr, 3) if tr else 0,
                    "pct_of_molecules": round(100 * nm / tm, 3) if tm else 0,
                    "pct_shift": round((100 * nm / tm if tm else 0)
                                       - (100 * nr / tr if tr else 0), 3),
                })
                dup_rows.append({
                    "sample_key": sk, "tx_class": cls, "reads": nr, "molecules": nm,
                    "reads_per_molecule": round(nr / nm, 2) if nm else None,
                })
        cdf, ddf = pd.DataFrame(cmp_rows), pd.DataFrame(dup_rows)
        _write(cdf, comp / "phase3b_read_vs_molecule.tsv", VERSION, today)
        _write(ddf, comp / "phase3b_class_duplication.tsv", VERSION, today)

        print("\n=== reads per molecule WITHIN each class ===")
        print("    (uniform across classes -> read-level composition was unbiased;")
        print("     higher for short-span classes -> length-biased duplication)")
        dpiv = ddf.pivot_table(index="sample_key", columns="tx_class",
                               values="reads_per_molecule")
        dpiv = dpiv[[c for c in CLASS_ORDER if c in dpiv.columns]]
        print(dpiv.round(1).to_string())

        print("\n=== composition shift: molecule % minus read % ===")
        spiv = cdf.pivot_table(index="sample_key", columns="tx_class", values="pct_shift",
                               fill_value=0)
        spiv = spiv[[c for c in CLASS_ORDER if c in spiv.columns]]
        print(spiv.round(2).to_string())
    else:
        print("\n(no Phase 3 read-level parquet found; skipping the comparison tables)")

    print(f"\nWritten to {comp}/:")
    print("  phase3b_molecule_counts.tsv, phase3b_read_vs_molecule.tsv, "
          "phase3b_class_duplication.tsv")


def _write(df, path, ver, today):
    with open(path, "w") as fh:
        fh.write(f"# experiment=EXP26000892; pipeline_version={ver}; "
                 f"classifier={CLASSIFY_VERSION}; generated={today}\n")
        df.to_csv(fh, sep="\t", index=False)
    df.to_parquet(str(path).replace(".tsv", ".parquet"), index=False)


if __name__ == "__main__":
    main()
