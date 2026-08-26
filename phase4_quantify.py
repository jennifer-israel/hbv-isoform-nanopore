#!/usr/bin/env python3
"""
Phase 4 — baseline HBV abundance per condition, with confidence intervals.

Combines:
  - Phase 1 alignment summary  (composite primary-mapped denominator)
  - Phase 2 HBV-UMI summary    (HBV reads, unique molecules, duplication)
  - Phase 3 classification      (per-read transcript class)

Produces two tables (TSV + parquet) under analysis/comparison/:
  phase4_recovery.tsv          per condition: HBV reads, unique molecules,
                               HBV-per-million-mapped (raw + unique), duplication,
                               each with a CI (Poisson on counts, Wilson on rates)
  phase4_class_abundance.tsv   per condition × transcript class: read count,
                               proportion of classified HBV reads, Poisson CI on
                               the count, Wilson CI on the proportion

CIs: counts are single/double digits, so intervals are essential. Poisson exact
(chi-square) on absolute counts; Wilson on proportions/rates. Robust conclusions
survive the intervals (polyA ≫ NOpolyA); sparse per-class calls will not — by design.

Usage:  conda activate hbv_lr_analysis; python scripts/phase4_quantify.py
"""
import datetime
from pathlib import Path
import pandas as pd
from scipy.stats import chi2
from statsmodels.stats.proportion import proportion_confint

PROJECT_ROOT = Path("/data/EXP26000993")
ANALYSIS = PROJECT_ROOT / "analysis"
COMP = ANALYSIS / "comparison"
SAMPLES = PROJECT_ROOT / "config" / "samples.tsv"
VERSION = "phase4-v1"
CLASS_ORDER = ["preS2_S", "preS1", "precore", "pgRNA", "pgRNA_RT", "X", "spliced",
               "unclassified", "antisense"]


def poisson_ci(k, alpha=0.05):
    """Exact (chi-square) Poisson 95% CI for a count k."""
    lo = chi2.ppf(alpha / 2, 2 * k) / 2 if k > 0 else 0.0
    hi = chi2.ppf(1 - alpha / 2, 2 * (k + 1)) / 2
    return lo, hi


def wilson(k, n):
    if not n:
        return (float("nan"), float("nan"))
    return proportion_confint(k, n, alpha=0.05, method="wilson")


def load_samples():
    rows = []
    for line in SAMPLES.read_text().splitlines():
        if line.startswith("#") or line.startswith("barcode\t") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) >= 4:
            rows.append({"barcode": f[0], "lib_id": f[1], "sample_name": f[2],
                         "sample_key": f[3], "input_ng": f[4], "polya": f[5], "pcr": f[6]})
    return pd.DataFrame(rows)


def read_tsv(p):
    return pd.read_csv(p, sep="\t", comment="#") if p.exists() else pd.DataFrame()


def main():
    meta = load_samples()
    align = read_tsv(COMP / "phase1_align_summary.tsv")
    umi   = read_tsv(COMP / "phase2_hbv_umi_summary.tsv")
    today = datetime.date.today().isoformat()

    def composite_mapped(sk):
        """Primary-mapped from the per-sample flagstat (robust; summary may miss the probe lib)."""
        fs = ANALYSIS / "samples" / sk / "flagstat.txt"
        if fs.exists():
            for line in fs.read_text().splitlines():
                if "primary mapped" in line:
                    return int(line.split()[0])
        if len(align):
            row = align[align.sample_key == sk]
            if len(row):
                return int(row.primary_mapped.iloc[0])
        return 0

    # ---- per-condition recovery table ----
    rows = []
    for _, m in meta.iterrows():
        sk = m.sample_key
        u = umi[umi.sample_key == sk]
        hbv_reads = int(u.hbv_primary_reads.iloc[0]) if len(u) else 0
        uniq      = int(u.hbv_unique_molecules.iloc[0]) if len(u) else 0
        dup       = u.hbv_dup_rate.iloc[0] if len(u) else "NA"
        comp = composite_mapped(sk)
        rlo, rhi = poisson_ci(hbv_reads)
        ulo, uhi = poisson_ci(uniq)
        def perM(k):     # per million composite-mapped
            return round(k / comp * 1e6, 3) if comp else float("nan")
        rows.append({
            "sample_key": sk, "barcode": m.barcode, "input_ng": m.input_ng,
            "polya": m.polya, "pcr_cyc": m.pcr, "composite_mapped": comp,
            "hbv_reads": hbv_reads, "hbv_reads_ci": f"[{rlo:.1f}, {rhi:.1f}]",
            "hbv_unique": uniq,     "hbv_unique_ci": f"[{ulo:.1f}, {uhi:.1f}]",
            "dup_rate": dup,
            "hbv_reads_per_M": perM(hbv_reads), "hbv_unique_per_M": perM(uniq),
            "hbv_reads_per_M_ci": f"[{perM(rlo)}, {perM(rhi)}]",
            "hbv_unique_per_M_ci": f"[{perM(ulo)}, {perM(uhi)}]",
        })
    rec = pd.DataFrame(rows).sort_values("hbv_unique", ascending=False)
    _write(rec, COMP / "phase4_recovery.tsv", today)

    # ---- per-condition × class abundance table ----
    crows = []
    for _, m in meta.iterrows():
        sk = m.sample_key
        cp = ANALYSIS / "samples" / sk / "hbv_classified.parquet"
        if not cp.exists():
            continue
        df = pd.read_parquet(cp)
        total = len(df)
        vc = df.tx_class.value_counts().to_dict()
        for cls in CLASS_ORDER:
            n = int(vc.get(cls, 0))
            if n == 0 and cls not in vc:
                continue
            lo, hi = poisson_ci(n)
            plo, phi = wilson(n, total)
            crows.append({
                "sample_key": sk, "barcode": m.barcode, "polya": m.polya, "pcr_cyc": m.pcr,
                "tx_class": cls, "n_reads": n, "n_reads_ci": f"[{lo:.1f}, {hi:.1f}]",
                "total_hbv": total, "proportion": round(n / total, 4) if total else 0,
                "proportion_ci": f"[{plo:.3f}, {phi:.3f}]",
            })
    cls_df = pd.DataFrame(crows)
    _write(cls_df, COMP / "phase4_class_abundance.tsv", today)

    print("=== Per-condition recovery (sorted by unique HBV molecules) ===")
    print(rec[["sample_key", "polya", "pcr_cyc", "hbv_reads", "hbv_unique",
               "hbv_unique_ci", "hbv_unique_per_M", "dup_rate"]].to_string(index=False))
    print("\n=== Class abundance (conditions with HBV reads) ===")
    if len(cls_df):
        piv = cls_df.pivot_table(index="sample_key", columns="tx_class",
                                 values="n_reads", fill_value=0, aggfunc="sum")
        print(piv.to_string())
    print(f"\nWritten: phase4_recovery.tsv, phase4_class_abundance.tsv")


def _write(df, path, today):
    with open(path, "w") as fh:
        fh.write(f"# experiment=EXP26000892_cDNA003; pipeline_version={VERSION}; generated={today}\n")
        df.to_csv(fh, sep="\t", index=False)
    df.to_parquet(str(path).replace(".tsv", ".parquet"), index=False)


if __name__ == "__main__":
    main()
