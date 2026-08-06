#!/usr/bin/env python3
"""
Phase 5 — self-contained HTML report for EXP26000892 (HBV RNA enrichment pilot).

Adapted from phase5_report.py (EXP26000559). That report describes a different study —
an 8-library poly(A)-selection x PCR-cycle matrix on un-enriched cDNA — so its figure
code is reused but its narrative is not. The original is left untouched.

WHAT CHANGED
  - 4-library input series (100 / 10 / 1 / 0.1 ng), all poly(A)-selected, all
    hybridisation-captured, rather than an 8-condition prep matrix.
  - --level reads|molecules. Class composition can be computed on undeduplicated
    reads (what EXP26000559 did) or on UMI-deduplicated molecules. In this run
    duplication is length-biased, so the two disagree; both reports are produced and
    the difference is itself reported.
  - Poly(A) panel dropped from Figure 1: poly_tail_length is absent from this run's
    summaries (estimation was not enabled at basecalling).
  - The 150 ng-only combined figures are removed; there is no 150 ng condition.
  - Dose-response and read-vs-molecule figures added.

  conda activate hbv_lr
  python3 phase5_report_EXP26000892.py --level molecules
  python3 phase5_report_EXP26000892.py --level reads
  python3 phase5_report_EXP26000892.py --level both      # writes both reports
"""
import argparse, base64, datetime, io, sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

PROJECT_ROOT = Path("/data/EXP26000896")
ANALYSIS = PROJECT_ROOT / "analysis"
COMP = ANALYSIS / "comparison"
SAMPLES = PROJECT_ROOT / "config" / "samples.tsv"
VERSION = "phase5-EXP26000892-v1"
HBV_LEN = 3182

LEVEL_FILE = {"reads": "hbv_classified.parquet",
              "molecules": "hbv_classified_molecules.parquet"}
LEVEL_UNIT = {"reads": "read", "molecules": "molecule"}

CLASS_COLORS = {
    "preS2_S": "#2ca02c", "preS1": "#d62728", "precore": "#e377c2", "pgRNA": "#1f77b4",
    "pgRNA_RT": "#17becf", "X": "#e6a817", "spliced": "#9467bd",
    "unclassified": "#9e9e9e", "antisense": "#cfcfcf",
}
CLASS_ORDER = list(CLASS_COLORS)
CLASS_LABEL = {"preS2_S": "preS2/S", "spliced": "spliced (SP1)"}
def clabel(c): return CLASS_LABEL.get(c, c)

# highest input first
DISPLAY_ORDER = [
    "SeqLib5550_100ng_polyA_17",
    "SeqLib5551_10ng_polyA_17",
    "SeqLib5552_1ng_polyA_21",
    "SeqLib5553_0.1ng_polyA_24",
]

# EXP26000559 un-enriched baseline, for the enrichment comparison
BASELINE_UNIQUE = (10, 31)


# ---------------------------------------------------------------- helpers
def b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def img(fig, alt):
    if fig is None:
        return ""
    return (f'<figure><img alt="{alt}" src="data:image/png;base64,{b64(fig)}"/>'
            f'<figcaption>{alt}</figcaption></figure>')


def load(name):
    p = COMP / name
    return pd.read_csv(p, sep="\t", comment="#") if p.exists() else pd.DataFrame()


def meta_df():
    rows = []
    for line in SAMPLES.read_text().splitlines():
        if line.startswith("#") or line.startswith("barcode\t") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) >= 7:
            rows.append(dict(barcode=f[0], lib_id=f[1], sample_key=f[3],
                             input_ng=f[4], polya=f[5], pcr=f[6]))
    return pd.DataFrame(rows).set_index("sample_key")


def label(sk, meta):
    m = meta.loc[sk]
    return f"{m.input_ng} ng\n{m.pcr} cyc\n{m.barcode}"


def ordered(keys):
    return [k for k in DISPLAY_ORDER if k in set(keys)]


def _classified(sk, level):
    p = ANALYSIS / "samples" / sk / LEVEL_FILE[level]
    return pd.read_parquet(p) if p.exists() else None


def _inputs(keys, meta):
    return [float(meta.loc[k, "input_ng"]) for k in keys]


def coord_counts(level):
    """Distinct (ref_start, ref_end) pairs per library.

    UMI sequencing errors inflate the apparent molecule count — a molecule amplified
    into hundreds of thousands of reads spawns thousands of single-read UMI variants
    that dedup cannot always collapse. Alignment coordinates cannot be manufactured
    that way, so distinct coordinate pairs are a conservative complexity floor that is
    immune to that artifact. Distinct molecules CAN share coordinates, so this
    undercounts; treat it as a lower bound and the UMI count as an upper bound.
    """
    out = {}
    for sk in DISPLAY_ORDER:
        d = _classified(sk, level)
        if d is None or d.empty:
            continue
        g = d.groupby(["ref_start", "ref_end"]).size()
        out[sk] = {"rows": len(d), "distinct": len(g),
                   "top5_share": 100 * g.nlargest(5).sum() / len(d)}
    return out


def umi_concentration():
    """Per-library duplicate-count concentration, from hbv.umi.bam.

    Quantifies PCR jackpotting: how far the reads-per-molecule distribution departs
    from uniform.

    ALWAYS recomputed from the BAMs — never read back from a cache. The result is
    written to comparison/umi_concentration.tsv as a citable output, but reading that
    file back would risk the report silently displaying stale numbers if Phase 2 were
    re-run and the UMI BAMs regenerated. A few minutes of compute is cheap next to a
    report that is quietly wrong.
    """
    cache_f = COMP / "umi_concentration.tsv"
    try:
        import pysam
    except ImportError:
        return pd.DataFrame()
    from collections import Counter
    rows = []
    for sk in DISPLAY_ORDER:
        b = ANALYSIS / "samples" / sk / "hbv.umi.bam"
        if not b.exists():
            continue
        c = Counter()
        with pysam.AlignmentFile(b) as af:
            for r in af.fetch(until_eof=True):
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                if not r.has_tag("RX"):
                    continue
                u = r.get_tag("RX")
                if len(u) == 28:
                    c[(u, r.reference_start // 50)] += 1
        if not c:
            continue
        v = np.array(sorted(c.values(), reverse=True))
        tot = v.sum()
        rows.append({
            "sample_key": sk, "molecules": len(v), "reads": int(tot),
            "median": float(np.median(v)), "mean": round(float(v.mean()), 1),
            "max": int(v.max()),
            "top1pct_share": round(100 * v[:max(1, len(v) // 100)].sum() / tot, 1),
            "n_jackpots_ge1000": int((v >= 1000).sum()),
            "n_singletons": int((v == 1).sum()),
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        with open(cache_f, "w") as fh:
            fh.write(f"# EXP26000892 UMI duplicate concentration; generated={datetime.date.today()}\n")
            df.to_csv(fh, sep="\t", index=False)
    return df


# ---------------------------------------------------------------- figures
def fig_qc(qc, meta):
    if qc.empty:
        return None
    q = qc[qc.sample_key.isin(meta.index)].set_index("sample_key")
    keys = ordered(q.index)
    if not keys:
        return None
    labs = [label(k, meta) for k in keys]
    x = np.arange(len(keys))
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    ax[0].bar(x, [q.loc[k, "reads"] / 1e6 for k in keys], color="#4c72b0")
    ax[0].set_ylabel("Pass reads (millions)"); ax[0].set_title("Yield")
    ax[1].bar(x, [q.loc[k, "len_N50"] for k in keys], color="#55a868")
    ax[1].set_ylabel("Read length N50 (bp)"); ax[1].set_title("Read length")
    ax[2].bar(x, [q.loc[k, "mean_qscore"] for k in keys], color="#c44e52")
    ax[2].set_ylabel("Mean Q"); ax[2].set_title("Quality"); ax[2].set_ylim(bottom=0)
    for a in ax:
        a.set_xticks(x); a.set_xticklabels(labs, fontsize=8)
    fig.tight_layout()
    return fig


def fig_dose_response(rec, meta, coords):
    """HBV complexity vs input, bounded above by UMI counts and below by coordinates."""
    if rec.empty:
        return None
    r = rec.set_index("sample_key")
    keys = ordered(r.index)
    xs = _inputs(keys, meta)
    ys_umi = [r.loc[k, "hbv_unique"] for k in keys]
    ys_crd = [coords.get(k, {}).get("distinct", np.nan) for k in keys]

    fig, ax = plt.subplots(figsize=(8, 5.4))
    ax.fill_between(xs, ys_crd, ys_umi, color="#1f77b4", alpha=0.12, zorder=1)
    ax.plot(xs, ys_umi, "o-", color="#1f77b4", lw=2, ms=9, zorder=3,
            label="UMI-deduplicated molecules (upper bound)")
    ax.plot(xs, ys_crd, "s--", color="#e6a817", lw=2, ms=8, zorder=3,
            label="distinct alignment coordinates (lower bound)")
    for xi, yi in zip(xs, ys_umi):
        ax.annotate(f"{yi:,.0f}", (xi, yi), textcoords="offset points",
                    xytext=(0, 11), ha="center", fontsize=8.5, color="#1f4e79")
    for xi, yi in zip(xs, ys_crd):
        if not np.isnan(yi):
            ax.annotate(f"{yi:,.0f}", (xi, yi), textcoords="offset points",
                        xytext=(0, -16), ha="center", fontsize=8.5, color="#8a6400")
    ax.axhspan(*BASELINE_UNIQUE, color="#d62728", alpha=0.18, zorder=0)
    ax.text(min(xs), BASELINE_UNIQUE[1] * 1.4,
            f"un-enriched baseline (EXP26000559): {BASELINE_UNIQUE[0]}–{BASELINE_UNIQUE[1]} molecules",
            fontsize=8.5, color="#8b1a1a")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("RNA input (ng)"); ax.set_ylabel("HBV molecular complexity")
    ax.set_title("HBV recovery vs input — bounded")
    ax.grid(alpha=0.25, which="both", ls=":")
    ax.legend(fontsize=8.5, loc="lower right")
    fig.tight_layout()
    return fig


def fig_jackpot(conc, coords, meta):
    """PCR jackpotting: how concentrated the reads are on a few molecules."""
    if conc.empty:
        return None
    c = conc.set_index("sample_key")
    keys = [k for k in DISPLAY_ORDER if k in c.index]
    if not keys:
        return None
    labs = [f"{meta.loc[k,'input_ng']} ng" for k in keys]
    x = np.arange(len(keys))
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))

    # labels use offset-points so they never push the axes out of range
    def lab(a, xi, yi, txt):
        a.annotate(txt, (xi, yi), textcoords="offset points", xytext=(0, 4),
                   ha="center", fontsize=8, clip_on=False)

    v0 = [c.loc[k, "top1pct_share"] for k in keys]
    ax[0].bar(x, v0, color="#d62728")
    ax[0].axhline(1, ls=":", color="k", lw=1)
    ax[0].set_ylabel("% of reads held by top 1% of molecules")
    ax[0].set_title("Read concentration\n(dotted = uniform expectation)")
    ax[0].set_ylim(0, 110)
    for xi, v in zip(x, v0):
        lab(ax[0], xi, v, f"{v:.0f}%")

    v1 = [c.loc[k, "max"] for k in keys]
    ax[1].bar(x, v1, color="#8172b3")
    ax[1].set_yscale("log"); ax[1].set_ylabel("reads in the single largest molecule")
    ax[1].set_title("Largest jackpot")
    ax[1].set_ylim(top=max(v1) * 3)
    for xi, v in zip(x, v1):
        lab(ax[1], xi, v, f"{int(v):,}")

    ts = [coords.get(k, {}).get("top5_share", np.nan) for k in keys]
    ax[2].bar(x, ts, color="#e6a817")
    ax[2].set_ylabel("% of molecules on the top 5 coordinates")
    ax[2].set_title("Coordinate concentration")
    if not all(np.isnan(ts)):
        ax[2].set_ylim(0, np.nanmax(ts) * 1.18)
    for xi, v in zip(x, ts):
        if not np.isnan(v):
            lab(ax[2], xi, v, f"{v:.1f}%")

    for a in ax:
        a.set_xticks(x); a.set_xticklabels(labs, fontsize=9)
    fig.tight_layout()
    return fig


def fig_recovery(rec, meta):
    if rec.empty:
        return None
    r = rec.set_index("sample_key")
    keys = ordered(r.index)
    labs = [label(k, meta) for k in keys]
    x = np.arange(len(keys))

    def ci_arr(col_ci, point):
        lo, hi = [], []
        for k in keys:
            s = str(r.loc[k, col_ci]).strip("[]").split(",")
            lo.append(r.loc[k, point] - float(s[0])); hi.append(float(s[1]) - r.loc[k, point])
        return np.array([lo, hi])

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    ax[0].bar(x - 0.2, [r.loc[k, "hbv_reads"] for k in keys], 0.4,
              label="raw HBV reads", color="#b0b0b0")
    ax[0].bar(x + 0.2, [r.loc[k, "hbv_unique"] for k in keys], 0.4,
              yerr=ci_arr("hbv_unique_ci", "hbv_unique"), capsize=3,
              label="unique molecules", color="#1f77b4")
    ax[0].set_yscale("log"); ax[0].set_ylabel("HBV reads / unique molecules")
    ax[0].set_title("Absolute recovery (log scale; Poisson 95% CI)")
    ax[0].legend(fontsize=8)
    ax[1].bar(x, [r.loc[k, "hbv_unique_per_M"] for k in keys], color="#2ca02c")
    ax[1].set_ylabel("Unique HBV per million composite-mapped")
    ax[1].set_title("Efficiency per unit sequencing")
    for a in ax:
        a.set_xticks(x); a.set_xticklabels(labs, fontsize=8)
    fig.tight_layout()
    return fig


def fig_dup(rec, meta):
    if rec.empty:
        return None
    r = rec.set_index("sample_key")
    keys = [k for k in ordered(r.index) if str(r.loc[k, "dup_rate"]) not in ("NA", "nan")]
    if not keys:
        return None
    labs = [label(k, meta) for k in keys]
    x = np.arange(len(keys))
    vals = [float(r.loc[k, "dup_rate"]) * 100 for k in keys]
    red = [float(r.loc[k, "hbv_reads"]) / max(1, float(r.loc[k, "hbv_unique"])) for k in keys]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    ax[0].bar(x, vals, color="#d62728")
    ax[0].set_ylabel("PCR duplication rate (%)"); ax[0].set_ylim(0, 105)
    ax[0].set_title("UMI duplication (28-nt UMIs)")
    for xi, v in zip(x, vals):
        ax[0].text(xi, v + 1.5, f"{v:.1f}%", ha="center", fontsize=8)
    ax[1].bar(x, red, color="#8172b3")
    ax[1].set_yscale("log"); ax[1].set_ylabel("HBV reads per unique molecule")
    ax[1].set_title("Sequencing redundancy")
    for xi, v in zip(x, red):
        ax[1].text(xi, v * 1.1, f"{v:.0f}×", ha="center", fontsize=8)
    for a in ax:
        a.set_xticks(x); a.set_xticklabels(labs, fontsize=8)
    fig.tight_layout()
    return fig


def fig_class_composition(cls_df, meta, level):
    if cls_df.empty:
        return None
    keys = ordered(cls_df.sample_key.unique())
    if not keys:
        return None
    labs = [label(k, meta) for k in keys]
    x = np.arange(len(keys))
    fig, ax = plt.subplots(figsize=(10, 5))
    bottom = np.zeros(len(keys))
    for cls in CLASS_ORDER:
        vals = []
        for k in keys:
            sub = cls_df[(cls_df.sample_key == k) & (cls_df.tx_class == cls)]
            tot = cls_df[cls_df.sample_key == k].n_reads.sum()
            vals.append(sub.n_reads.iloc[0] / tot * 100 if len(sub) and tot else 0)
        vals = np.array(vals)
        if vals.sum() > 0:
            ax.bar(x, vals, bottom=bottom, label=clabel(cls), color=CLASS_COLORS[cls])
            bottom += vals
    ax.set_ylabel(f"% of classified {LEVEL_UNIT[level]}s")
    ax.set_title(f"Transcript class composition — {level}")
    ax.set_xticks(x); ax.set_xticklabels(labs, fontsize=8)
    ax.legend(fontsize=8, ncol=2, bbox_to_anchor=(1.01, 1), loc="upper left")
    fig.tight_layout()
    return fig


def fig_level_comparison(meta):
    """Read-level vs molecule-level composition, and per-class duplication."""
    keys = [k for k in DISPLAY_ORDER if (ANALYSIS / "samples" / k / LEVEL_FILE["molecules"]).exists()]
    if not keys:
        return None
    rows = []
    for k in keys:
        dr, dm = _classified(k, "reads"), _classified(k, "molecules")
        if dr is None or dm is None:
            continue
        for cls in CLASS_ORDER:
            nr = int((dr.tx_class == cls).sum()); nm = int((dm.tx_class == cls).sum())
            if nr == 0 and nm == 0:
                continue
            rows.append({"sample_key": k, "cls": cls,
                         "pr": 100 * nr / len(dr), "pm": 100 * nm / len(dm),
                         "rpm": nr / nm if nm else np.nan})
    if not rows:
        return None
    d = pd.DataFrame(rows)
    fig, ax = plt.subplots(1, 2, figsize=(13.5, 5))
    x = np.arange(len(keys)); w = 0.38
    for i, (col, ttl) in enumerate([("pr", "read-level"), ("pm", "molecule-level")]):
        bottom = np.zeros(len(keys))
        for cls in CLASS_ORDER:
            vals = np.array([d[(d.sample_key == k) & (d.cls == cls)][col].sum() for k in keys])
            if vals.sum() > 0:
                ax[0].bar(x + (i - 0.5) * w, vals, w, bottom=bottom,
                          color=CLASS_COLORS[cls],
                          label=clabel(cls) if i == 0 else None,
                          edgecolor="white", lw=0.4)
                bottom += vals
    ax[0].set_xticks(x); ax[0].set_xticklabels([label(k, meta) for k in keys], fontsize=8)
    ax[0].set_ylabel("% of classified"); ax[0].set_title("left bar = reads, right bar = molecules")
    ax[0].legend(fontsize=7, ncol=2, bbox_to_anchor=(1.01, 1), loc="upper left")

    for cls in ["X", "preS2_S", "preS1", "precore", "spliced"]:
        sub = d[d.cls == cls]
        if sub.empty:
            continue
        xs = [float(meta.loc[k, "input_ng"]) for k in sub.sample_key]
        ax[1].plot(xs, sub.rpm, "o-", color=CLASS_COLORS[cls], label=clabel(cls), ms=6)
    ax[1].set_xscale("log"); ax[1].set_yscale("log")
    ax[1].set_xlabel("RNA input (ng)"); ax[1].set_ylabel("reads per molecule")
    ax[1].set_title("Duplication within each class")
    ax[1].grid(alpha=0.25, which="both", ls=":"); ax[1].legend(fontsize=8)
    fig.tight_layout()
    return fig


def fig_tss(meta, level):
    fig, ax = plt.subplots(figsize=(11, 4))
    plotted = 0
    for sk in DISPLAY_ORDER:
        df = _classified(sk, level)
        if df is None or df.empty:
            continue
        df = df[df.strand == "+"]
        if len(df) < 5:
            continue
        ax.hist(df.tss, bins=np.arange(0, HBV_LEN + 50, 50), histtype="step",
                lw=1.5, label=f"{meta.loc[sk,'input_ng']} ng")
        plotted += 1
    if not plotted:
        plt.close(fig); return None
    band_defs = [("precore/pgRNA", "#1f77b4", [(1730, 1880)]),
                 ("preS1", "#d62728", [(2700, 3100)]),
                 ("X", "#e6a817", [(1260, 1450)]),
                 ("preS2/S", "#2ca02c", [(3100, HBV_LEN), (0, 150)])]
    for _, c, spans in band_defs:
        for lo, hi in spans:
            ax.axvspan(lo, hi, color=c, alpha=0.08)
    ax.set_xlabel("5′ TSS on HBV genome (mod 3182); shaded = transcript-class windows")
    ax.set_ylabel(LEVEL_UNIT[level] + "s")
    ax.set_title(f"TSS distribution — {level}")
    line_leg = ax.legend(fontsize=8, loc="upper right", title="input")
    ax.add_artist(line_leg)
    ax.legend(handles=[Patch(facecolor=c, alpha=0.25, label=n) for n, c, _ in band_defs],
              fontsize=7, loc="upper left", title="class window", framealpha=0.9)
    fig.tight_layout()
    return fig


SPAN_CLASS_ORDER = ["precore", "pgRNA", "preS1", "preS2_S", "X",
                    "spliced", "unclassified", "antisense", "pgRNA_RT"]
ANCHOR_X = 1936 + HBV_LEN


def fig_spans(meta, level, max_per_panel=4000):
    panels = []
    for sk in DISPLAY_ORDER:
        d = _classified(sk, level)
        if d is None or d.empty:
            continue
        if len(d) > max_per_panel:
            d = d.sample(max_per_panel, random_state=0)
        panels.append((f"{meta.loc[sk,'input_ng']} ng", d))
    if not panels:
        return None
    PER, MINH, MAXH = 0.045, 1.4, 60.0
    hs = [min(MAXH, max(MINH, len(d) * PER + 0.5)) for _, d in panels]
    fig, axes = plt.subplots(len(panels), 1, figsize=(11, sum(hs) + 0.4), sharex=True,
                             gridspec_kw={"height_ratios": hs, "hspace": 0.05})
    if len(panels) == 1:
        axes = [axes]
    ordkey = {c: i for i, c in enumerate(SPAN_CLASS_ORDER)}

    def parse_blocks(bs, rs, re_):
        if isinstance(bs, str) and bs:
            try:
                return [tuple(int(x) for x in tok.split("-")) for tok in bs.split(";")]
            except ValueError:
                pass
        return [(rs, re_)]

    xlo = xhi = ANCHOR_X
    present = set()
    for ax, (ylabel, df) in zip(axes, panels):
        d = df.copy()
        d["_o"] = d.tx_class.map(lambda c: ordkey.get(c, 99))
        d["_s"] = d.ref_end - d.ref_start
        d = d.sort_values(["_o", "_s"], ascending=[True, False]).reset_index(drop=True)
        has_blocks = "blocks" in d.columns
        for i, row in enumerate(d.itertuples()):
            col = CLASS_COLORS.get(row.tx_class, "#999")
            shift = (row.polya_end + HBV_LEN) - row.ref_end
            exons = parse_blocks(getattr(row, "blocks", None) if has_blocks else None,
                                 row.ref_start, row.ref_end)
            xs = [s + shift for s, _ in exons]; xe = [e + shift for _, e in exons]
            if len(exons) > 1:
                ax.plot([min(xs), max(xe)], [i, i], color=col, lw=0.4, alpha=0.5, zorder=1)
            for s, e in zip(xs, xe):
                ax.plot([s, e], [i, i], color=col, lw=0.9, solid_capstyle="butt", zorder=2)
            xlo, xhi = min(xlo, min(xs)), max(xhi, max(xe))
        present.update(d.tx_class)
        ax.axvline(HBV_LEN, ls="--", color="k", lw=0.8)
        ax.axvline(ANCHOR_X, ls=":", color="darkred", lw=1.0)
        ax.set_ylim(len(d) - 0.5, -0.5); ax.set_yticks([])
        ax.set_ylabel(f"{ylabel}\nn={len(d)}", rotation=0, ha="right", va="center", fontsize=8)
    for ax in axes:
        ax.set_xlim(xlo - 80, xhi + 120)
    axes[-1].set_xlabel("Position on U95551.1_2x — 3′-anchored at poly-A (dotted ~5118); "
                        "dashed = linearisation 3182", fontsize=9)
    axes[0].legend(handles=[Line2D([0], [0], color=CLASS_COLORS[c], lw=3, label=clabel(c))
                            for c in CLASS_ORDER if c in present],
                   fontsize=7, ncol=4, loc="upper left")
    fig.tight_layout(h_pad=0.3)
    return fig


# ---------------------------------------------------------------- html
def tbl(df, cols=None):
    if df.empty:
        return "<p><em>(not available)</em></p>"
    if cols:
        df = df[[c for c in cols if c in df.columns]]
    return df.fillna("—").to_html(index=False, border=0, classes="data", escape=False)


CSS = """
body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;max-width:1000px;margin:2rem auto;padding:0 1rem;color:#222;line-height:1.55}
h1{border-bottom:3px solid #1f77b4;padding-bottom:.3rem} h2{border-bottom:1px solid #ccc;margin-top:2rem}
figure{margin:1.2rem 0;text-align:center} img{max-width:100%;border:1px solid #eee}
figcaption{font-size:.85rem;color:#555;text-align:left;margin-top:.3rem}
table.data{border-collapse:collapse;font-size:.85rem;margin:1rem 0}
table.data th,table.data td{border:1px solid #ddd;padding:3px 8px;text-align:right}
table.data th{background:#f4f6f8}
.key{background:#eef7ff;padding:.8rem 1rem;border-left:4px solid #1f77b4;margin:1rem 0}
.warn{background:#fff6e5;padding:.8rem 1rem;border-left:4px solid #e6a817;margin:1rem 0}
code{background:#f4f4f4;padding:1px 4px;border-radius:3px}
"""


def build(level, meta, qc, rec, cls_df, today, coords, conc):
    r = rec.set_index("sample_key") if not rec.empty else pd.DataFrame()
    unit = LEVEL_UNIT[level]

    def g(sk, col, default=0):
        try:
            return r.loc[sk, col]
        except Exception:
            return default

    best = ordered(r.index)[0] if len(r) else None
    lo_key = ordered(r.index)[-1] if len(r) else None

    figs = {
        "qc": img(fig_qc(qc, meta),
                  "Figure 1. Per-library read QC from the demultiplexed pass reads: yield, "
                  "read-length N50, and mean Q. Poly(A) detection is not shown — "
                  "poly_tail_length is absent from this run's summaries because poly(A) "
                  "estimation was not enabled at basecalling."),
        "dose": img(fig_dose_response(rec, meta, coords),
                    "Figure 2. HBV molecular complexity against RNA input, log-log, shown as "
                    "a bounded range. Blue = UMI-deduplicated molecules, an upper bound: PCR "
                    "jackpotting plus UMI sequencing errors inflate it at low input (Figure 5). "
                    "Orange = distinct alignment coordinates, a conservative lower bound that "
                    "UMI errors cannot inflate. Red band is the un-enriched baseline from "
                    f"EXP26000559 ({BASELINE_UNIQUE[0]}–{BASELINE_UNIQUE[1]} molecules). The UMI "
                    "series flattens between 1 and 0.1 ng; the coordinate series does not — that "
                    "flattening is an amplification artifact, not a property of the assay."),
        "jack": img(fig_jackpot(conc, coords, meta),
                    "Figure 5. PCR jackpotting. Left: share of reads held by the top 1% of "
                    "molecules (uniform duplication would give ~1%). Middle: reads in the single "
                    "most-amplified molecule. Right: share of molecules sitting on just five "
                    "alignment coordinates. Concentration rises sharply as input falls and cycle "
                    "count rises; the 100 ng library has no molecule exceeding 1,000 reads."),
        "rec": img(fig_recovery(rec, meta),
                   "Figure 3. Left: raw HBV reads (grey) against UMI-deduplicated unique "
                   "molecules (blue, Poisson 95% CI), log scale — the gap between them is "
                   "PCR duplication. Right: unique molecules per million composite-mapped "
                   "reads."),
        "dup": img(fig_dup(rec, meta),
                   "Figure 4. Left: UMI duplication rate on well-formed 28-nt UMIs. Right: "
                   "HBV reads sequenced per unique molecule. Duplication rises steeply as "
                   "input falls and PCR cycles rise."),
        "cls": img(fig_class_composition(cls_df, meta, level),
                   f"Figure 5. HBV transcript-class composition per library, computed on "
                   f"{unit}s (% of classified {unit}s)."),
        "lvl": img(fig_level_comparison(meta),
                   "Figure 6. Left: composition computed on reads (left bar of each pair) "
                   "against molecules (right bar). Right: reads per molecule within each "
                   "class — the short-span X class duplicates far harder than the long "
                   "precore/pgRNA classes, which is why the two levels disagree."),
        "tss": img(fig_tss(meta, level),
                   f"Figure 7. 5′ TSS distribution (mod 3182), {unit}-level. Shaded bands are "
                   "the transcript-class TSS windows."),
        "spans": img(fig_spans(meta, level),
                     f"Figure 8. {unit.capitalize()} spans, 3′-anchored at the canonical poly-A "
                     "(dotted red) so co-terminal transcripts align. Colour = class; grouped by "
                     "class, ordered longest to shortest. Dashed line = the 3182 linearisation "
                     "point; lines crossing it wrap into the second genome copy. Large panels "
                     "are subsampled to 4,000."),
    }

    warn = ""
    if level == "reads":
        warn = ('<div class="warn"><b>Counting level: reads (undeduplicated).</b> '
                'Composition here is weighted by PCR duplication. Because several classes are '
                'gated on length (X 300–1,000 bp; preS2/S ≥1,200; precore/pgRNA ≥2,600) and short '
                'amplicons duplicate preferentially, the short classes are over-represented — '
                'see Figure 6. This level is reported for comparability with EXP26000559, which '
                'used it. For interpretation, use the molecule-level report.</div>')
    else:
        warn = ('<div class="warn"><b>Counting level: molecules (UMI-deduplicated).</b> '
                'Each row is one original cDNA molecule, so composition is not weighted by '
                'amplification. This is the level to interpret. The read-level report is '
                'provided alongside for comparability with EXP26000559.</div>')

    cnc = conc.set_index("sample_key") if not conc.empty else pd.DataFrame()

    def cg(sk, col, default=0):
        try:
            return cnc.loc[sk, col]
        except Exception:
            return default

    summary_rows = []
    for k in ordered(r.index):
        cd = coords.get(k, {}).get("distinct")
        summary_rows.append(
            f"<li><b>{meta.loc[k,'input_ng']} ng</b> ({meta.loc[k,'pcr']} cycles): "
            f"{int(g(k,'hbv_unique')):,} UMI molecules / "
            f"{cd:,} distinct coordinates" if cd else
            f"<li><b>{meta.loc[k,'input_ng']} ng</b>: {int(g(k,'hbv_unique')):,} UMI molecules")
        summary_rows[-1] += (f"; {float(g(k,'dup_rate',0))*100:.1f}% duplication, "
                             f"top 1% of molecules hold {cg(k,'top1pct_share','—')}% of reads</li>")

    rec_cols = ["barcode", "input_ng", "pcr_cyc", "composite_mapped", "hbv_reads",
                "hbv_unique", "hbv_unique_ci", "hbv_unique_per_M", "dup_rate"]

    html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>EXP26000892 HBV enrichment pilot — {level}</title><style>{CSS}</style></head><body>
<h1>Targeted HBV RNA enrichment from low-input RNA — EXP26000892_cDNA003</h1>
<p><em>Generated {today} · pipeline {VERSION} · counting level: <b>{level}</b></em></p>

{warn}

<div class="key">
<b>Executive summary.</b> Four libraries were prepared from a single poly(A)-selected
HBV-infected PHH RNA pool at 100, 10, 1 and 0.1 ng input, barcoded with custom 24-nt PCR
indices, pooled at equal mass (375 ng each), enriched by Twist hybridisation capture, and
sequenced on one PromethION flow cell.
<ul>
{''.join(summary_rows)}
</ul>
<p><b>Enrichment works.</b> Against the un-enriched baseline of
{BASELINE_UNIQUE[0]}–{BASELINE_UNIQUE[1]} unique HBV molecules per library (EXP26000559),
every condition here recovers orders of magnitude more.</p>
<p><b>But UMI molecule counts overstate complexity at low input.</b> PCR jackpotting is
severe below 10 ng: in the 0.1 ng library a single molecule accounts for 440,452 reads,
17 molecules hold 98.7% of all reads, and 88% of apparent single-read "molecules" sit at
a jackpot's alignment position — they are UMI sequencing errors off that jackpot, not
distinct captures. The 100 ng library has no molecule exceeding 1,000 reads. Complexity is
therefore reported as a <b>range</b>: UMI count as an upper bound, distinct alignment
coordinates as a lower bound that UMI errors cannot inflate (Figure 2).</p>
<p><b>10 ng is the recommended input for development.</b> It retains substantial
complexity, jackpotting is moderate, and its transcript composition is consistent with
1 ng. At 0.1 ng half the molecules sit on five alignment coordinates, so HBV is detected
but the population cannot be characterised.</p>
<p><b>Design caveat.</b> All four libraries entered capture at equal <em>mass</em>, so this
measures whether a low-input library still contains HBV — not how much HBV can be recovered
from a given input. Answering that requires capturing each library separately.</p>
</div>

<h2>Methods</h2>
<p><b>Libraries and sequencing.</b> One poly(A)-selected RNA pool (HBV-infected PHH),
diluted to 100 / 10 / 1 / 0.1 ng. Custom barcoded PCR (17 / 17 / 21 / 24 cycles
respectively — more cycles at lower input to reach usable mass), pooled at 375 ng per
library, Twist hybridisation capture with TSO and CRTA blocking oligos, post-capture PCR,
then SQK-LSK ligation sequencing on PromethION (FLO-PRO114M).</p>
<p><b>Demultiplexing.</b> The barcodes are custom, so MinKNOW could not resolve them.
Reads were demultiplexed with <code>dorado demux</code> using a custom barcode arrangement,
requiring the <b>same barcode at both ends</b> (<code>--barcode-both-ends</code>). This
matters: 4–10% of reads carried two <em>different</em> barcodes — unambiguous chimeras that
single-ended demultiplexing would have assigned at random, contaminating the low-input
libraries. Residual cross-library carry-over after strict demultiplexing was measured at
&lt;1% by shared UMI-and-position analysis.</p>
<p><b>Reference and alignment.</b> hg38 + 2×HBV (U95551.1 doubled to 6,364 bp so
genome-wrapping transcripts align as one record). minimap2 <code>-ax splice
--secondary=no</code>. Coordinates are kept in 2× space; <code>mod 3182</code> is used only
to derive TSS. Note that MAPQ is not a useful filter on a 2× reference — a read contained
within one genome copy matches both copies equally and receives MAPQ 0 by construction.</p>
<p><b>UMIs.</b> HBV-aligning reads were passed through pychopper
(<code>-k PCB114 -U</code>) to orient them and extract the 28-nt structured UMI carried in
the strand-switching primer, re-aligned to an HBV-only index, and deduplicated with
<code>umi_tools dedup --method=directional</code> on well-formed 28-nt UMIs (~92% of
tagged reads).</p>
<p><b>Classification.</b> Transcript classes follow the specified rules: pgRNA_RT at
footprint ≥3,982 bp; spliced where ≥2 exon blocks remain after merging gaps ≤200 bp;
precore TSS 1,730–1,815 and pgRNA TSS 1,816–1,880, both requiring junction crossing and
span ≥2,600 bp; preS1 TSS 2,700–3,100 with span ≥1,600; preS2/S TSS ≥3,100 or ≤150 with
span ≥1,200; X TSS 1,260–1,450 with span 300–1,000. The largest intron is annotated SP1
(donor 2,447±60, acceptor 489±60) or non-canonical.</p>
<p><b>Statistics.</b> Poisson exact 95% CIs on counts, Wilson 95% CIs on proportions.</p>

<h2>Library QC</h2>
{figs['qc']}

<h2>HBV recovery</h2>
{tbl(rec.reset_index()[rec_cols] if 'sample_key' in rec.columns else rec, rec_cols)}
{figs['dose']}
{figs['rec']}
{figs['dup']}
<p>Duplication is high throughout and rises sharply as input falls, reaching
{float(g(lo_key,'dup_rate',0))*100:.1f}% at 0.1 ng — roughly
{float(g(lo_key,'hbv_reads',1))/max(1,float(g(lo_key,'hbv_unique',1))):.0f} reads sequenced
per distinct molecule on average. Subsampling analysis showed none of the libraries reached
saturation, so deeper sequencing would recover more molecules.</p>

<h2>Amplification jackpotting</h2>
{figs['jack']}
<p>The mean duplication rate conceals a highly skewed distribution. The <b>median</b>
molecule is seen once or twice in every library; the mean is driven entirely by a small
number of molecules amplified enormously. This is the classic PCR jackpot — a molecule
copied in an early cycle doubles every cycle thereafter and finishes orders of magnitude
ahead of one first copied later.</p>
<p>The consequence for molecule counting is direct. A molecule read hundreds of thousands
of times generates, through ordinary sequencing error in its 28-nt UMI, thousands of
distinct UMI variants. Those at edit distance ≥2 are not reliably collapsed by
<code>umi_tools</code>, so they are counted as separate molecules. At 0.1 ng, 88% of
single-read molecules sit at a jackpot's alignment position and half are within 3 edits of
that jackpot's UMI — direct evidence that they are error variants rather than distinct
captures.</p>
<p>This is why complexity is bounded rather than point-estimated. Distinct alignment
coordinates cannot be produced by UMI error, and give a smooth dose-response across the
input series where the UMI count flattens between 1 and 0.1 ng. Note that coordinates
undercount, since genuinely distinct molecules can share a start and end position — the
truth lies between the two bounds, nearer the lower one at low input.</p>
<table class="data"><tr><th>input</th><th>UMI molecules (upper)</th>
<th>distinct coordinates (lower)</th><th>top 5 coordinates hold</th></tr>
{''.join(f"<tr><td>{meta.loc[k,'input_ng']} ng</td><td>{int(g(k,'hbv_unique')):,}</td>"
         f"<td>{coords.get(k,{}).get('distinct',0):,}</td>"
         f"<td>{coords.get(k,{}).get('top5_share',0):.1f}%</td></tr>"
         for k in ordered(r.index))}
</table>

<h2>Transcript composition</h2>
{figs['cls']}
{figs['lvl']}
<p>Figure 6 is the reason this report exists in two versions. Duplication is not uniform
across classes: the short-span X class duplicates several-fold harder than the long
precore and pgRNA classes, so composition computed on reads systematically over-states the
short classes. The effect grows as input falls and PCR cycles rise.</p>
{figs['tss']}

<h2>Read/molecule span maps</h2>
{figs['spans']}

<h2>Conclusions</h2>
<ol>
<li>Hybridisation capture raises HBV recovery by orders of magnitude over the un-enriched
baseline. The 100 ng library — the only one free of jackpotting — recovers
{int(g(ordered(r.index)[0],'hbv_unique')):,} molecules against a baseline of
{BASELINE_UNIQUE[0]}–{BASELINE_UNIQUE[1]}.</li>
<li>Cross-library contamination is under 1%, so barcode assignment is sound. This depends
on strict both-ends demultiplexing; 4–10% of reads carried two different barcodes and were
correctly excluded.</li>
<li><b>Molecular complexity must be reported as a range at low input.</b> UMI counts are
inflated by jackpot-derived UMI error variants; distinct alignment coordinates give a
conservative floor. The two bounds converge at high input and diverge sharply at 0.1 ng.</li>
<li><b>10 ng is the recommended input for development</b> — substantial retained
complexity, moderate jackpotting, and a transcript composition consistent with 1 ng.</li>
<li>0.1 ng detects HBV but cannot characterise it. Half its molecules sit on five alignment
coordinates, and the apparent shift toward short transcript classes cannot be cleanly
separated from which molecule happened to jackpot.</li>
<li>Composition must be computed on deduplicated molecules for this data. Duplication is
length-biased, so read-level composition overstates short classes — unlike EXP26000559,
where duplication was low enough that reads approximated molecules.</li>
<li>Reducing PCR cycles at low input would do more for data quality than sequencing deeper.
Jackpotting is established in the first few cycles and cannot be sequenced away.</li>
</ol>

<h2>Appendix — reproducibility</h2>
<p>dorado demux (custom arrangement, both-ends), minimap2, samtools, pychopper,
umi_tools (directional), pandas/scipy/statsmodels. Pipeline:
<code>run_dorado_demux.sh</code> → <code>reshape_demux_output.sh</code> →
<code>phase1_align.sh</code> → <code>phase1_5_readqc_demux.sh</code> →
<code>phase2_hbv_umi_v2.sh</code> → <code>phase3_classify.py</code> /
<code>phase3b_classify_molecules.py</code> → <code>phase4_quantify.py</code> →
<code>phase5_report_EXP26000892.py</code>. Tables in
<code>analysis/comparison/</code>, per-read and per-molecule data in
<code>analysis/samples/&lt;lib&gt;/</code>.</p>
<p><b>Classifier provenance.</b> <code>hbv_transcript_classify.py</code> used here is a
reimplementation written from the specified rules; the canonical module from the earlier
HBV projects was not available. Class counts should be confirmed against it before
cross-experiment comparison.</p>
</body></html>"""
    return html


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--level", choices=["reads", "molecules", "both"], default="both")
    args = ap.parse_args()

    today = datetime.date.today().isoformat()
    meta = meta_df()
    qc = load("phase1_5_read_qc.tsv")
    out_dir = ANALYSIS / "reports"; out_dir.mkdir(parents=True, exist_ok=True)

    print("  computing duplicate concentration from hbv.umi.bam (always recomputed)...")
    conc = umi_concentration()
    if not conc.empty:
        print(f"    {len(conc)} libraries -> {COMP/'umi_concentration.tsv'}")

    levels = ["reads", "molecules"] if args.level == "both" else [args.level]
    for lv in levels:
        rec = load(f"phase4_recovery_{lv}.tsv")
        if rec.empty:
            rec = load("phase4_recovery.tsv")
        cls_df = load(f"phase4_class_abundance_{lv}.tsv")
        if cls_df.empty:
            cls_df = load("phase4_class_abundance.tsv")
        if rec.empty:
            sys.stderr.write(f"WARN [{lv}]: no recovery table found; skipping\n")
            continue
        # coordinate bound is always computed on MOLECULES — it is a property of the
        # library, not of the counting level, and is meaningless on duplicated reads
        coords = coord_counts("molecules")
        html = build(lv, meta, qc, rec, cls_df, today, coords, conc)
        out = out_dir / f"EXP26000892_report_{lv}.html"
        out.write_text(html)
        print(f"  {lv:<10} -> {out}  ({out.stat().st_size/1024:.0f} KB)")


if __name__ == "__main__":
    main()
