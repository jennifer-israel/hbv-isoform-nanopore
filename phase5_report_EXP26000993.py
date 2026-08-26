#!/usr/bin/env python3
"""
Phase 5 — HTML report for EXP26000993 / cDNA005 (serum & plasma enrichment pilot).

Adapted from phase5_report_EXP26000892.py. Figure code is reused; the narrative is not —
that report describes an input titration on one RNA source, this one is a SAMPLE TYPE
comparison at fixed 1 ng input, plus a cross-run carryover control.

  conda activate hbv_lr
  python3 phase5_report_EXP26000993.py                 # both levels
  python3 phase5_report_EXP26000993.py --level molecules

WHAT CHANGED FROM THE EXP26000892 VERSION
  - 4 libraries, all 1 ng: human plasma pool, two Yecuris chimeric mice (g21-treated and
    PBS control), and SeqLib5552 carried over from cDNA003 as a cross-run control.
  - Dose-response replaced by a per-sample recovery figure — there is no input series.
  - Adds the short-fraction (phase2s) results: molecules recovered from reads the splice
    alignment discards, and their HBV depletion relative to the long fraction.
  - Narrative covers the two findings that actually drive decisions: pooling high with low
    viral load, and read length limiting transcript classification.
"""
import argparse, base64, datetime, io, os, sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", "/data/EXP26000993"))
ANALYSIS = PROJECT_ROOT / "analysis"
COMP = ANALYSIS / "comparison"
SAMPLES = PROJECT_ROOT / "config" / "samples.tsv"
VERSION = "phase5-EXP26000993-v1"
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

# control first, then treated, then the two low-load libraries
DISPLAY_ORDER = [
    "SeqLib5577_yecuris2_1ng_20",
    "SeqLib5576_yecuris1_1ng_20",
    "SeqLib5575_humanPlasma_1ng_20",
    "SeqLib5552_cDNA003_1ng_21",
]
NICE = {
    "SeqLib5577_yecuris2_1ng_20":    "yecuris 2\n(PBS control)",
    "SeqLib5576_yecuris1_1ng_20":    "yecuris 1\n(g21 treated)",
    "SeqLib5575_humanPlasma_1ng_20": "human plasma\npool",
    "SeqLib5552_cDNA003_1ng_21":     "cDNA003\ncarryover",
}
# HBV as % of mapped reads that SeqLib5552 achieved in its OWN experiment — the
# reference point for the probe-competition finding
CARRYOVER_PRIOR_PCT = 2.58
CARRYOVER_KEY = "SeqLib5552_cDNA003_1ng_21"

# Concatemer rate in cDNA003 (EXP26000896), measured retrospectively with the SAME
# filter_concatemers.py at the same settings, on HBV reads regenerated from
# aligned_sorted.bam. Validated by reproducing the cDNA005 phase-2 numbers to 0.03 pp.
# (sample_key, input, pct, total_hbv_reads, concatemers)
CDNA003_CONCAT = [
    ("SeqLib5550_100ng_polyA_17", "100 ng", 0.63, 1027904, 6444),
    ("SeqLib5551_10ng_polyA_17",  "10 ng",  0.91,  425032, 3875),
    ("SeqLib5552_1ng_polyA_21",   "1 ng",   0.97,  482567, 4669),
    ("SeqLib5553_0.1ng_polyA_24", "0.1 ng", 0.71, 1277869, 9114),
]
# the same physical library, measured in both runs
SHARED_LIB_CDNA003_PCT = 0.97
SHARED_LIB_CDNA005_PCT = 5.43


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
    if not SAMPLES.exists():
        sys.stderr.write(f"WARN: {SAMPLES} not found; labels will be sample keys\n")
        return pd.DataFrame()
    for line in SAMPLES.read_text().splitlines():
        if line.startswith("#") or line.startswith("barcode\t") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) >= 7:
            rows.append(dict(barcode=f[0], lib_id=f[1], sample_name=f[2], sample_key=f[3],
                             input_ng=f[4], polya=f[5], pcr=f[6]))
    return pd.DataFrame(rows).set_index("sample_key") if rows else pd.DataFrame()


def label(sk, meta):
    return NICE.get(sk, sk)


def ordered(keys):
    return [k for k in DISPLAY_ORDER if k in set(keys)]


def _classified(sk, level):
    p = ANALYSIS / "samples" / sk / LEVEL_FILE[level]
    return pd.read_parquet(p) if p.exists() else None


def coord_counts(level="molecules"):
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
    """Duplicate-count concentration per library, from hbv.umi.bam. Always recomputed —
    written out as a citable table but never read back, so it cannot go stale."""
    out_f = COMP / "umi_concentration.tsv"
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
        v = np.array(sorted(c.values(), reverse=True)); tot = v.sum()
        rows.append({"sample_key": sk, "molecules": len(v), "reads": int(tot),
                     "median": float(np.median(v)), "mean": round(float(v.mean()), 1),
                     "max": int(v.max()),
                     "top1pct_share": round(100 * v[:max(1, len(v)//100)].sum()/tot, 1),
                     "n_jackpots_ge1000": int((v >= 1000).sum())})
    df = pd.DataFrame(rows)
    if not df.empty:
        with open(out_f, "w") as fh:
            fh.write(f"# EXP26000993 UMI duplicate concentration; generated={datetime.date.today()}\n")
            df.to_csv(fh, sep="\t", index=False)
    return df


# ---------------------------------------------------------------- figures
def fig_qc(qc, meta):
    if qc.empty:
        return None
    q = qc.set_index("sample_key")
    keys = ordered(q.index)
    if not keys:
        return None
    labs = [label(k, meta) for k in keys]
    x = np.arange(len(keys))
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    ax[0].bar(x, [q.loc[k, "reads"]/1e6 for k in keys], color="#4c72b0")
    ax[0].set_ylabel("Pass reads (millions)"); ax[0].set_title("Yield")
    ax[1].bar(x, [q.loc[k, "len_N50"] for k in keys], color="#55a868")
    ax[1].set_ylabel("Read length N50 (bp)"); ax[1].set_title("Read length")
    ax[1].axhline(218, ls="--", color="k", lw=1)
    ax[1].annotate("218 bp = adapter construct,\nno insert", (0.02, 240),
                   fontsize=7.5, color="#444")
    ax[2].bar(x, [q.loc[k, "mean_qscore"] for k in keys], color="#c44e52")
    ax[2].set_ylabel("Mean Q"); ax[2].set_title("Quality"); ax[2].set_ylim(bottom=0)
    for a in ax:
        a.set_xticks(x); a.set_xticklabels(labs, fontsize=8)
    fig.tight_layout()
    return fig


def fig_recovery_by_sample(rec, short, meta):
    """Unique HBV molecules per library, long fraction plus short fraction."""
    if rec.empty:
        return None
    r = rec.set_index("sample_key")
    keys = ordered(r.index)
    labs = [label(k, meta) for k in keys]
    x = np.arange(len(keys))
    long_v = np.array([float(r.loc[k, "hbv_unique"]) for k in keys])
    s = short.set_index("sample_key") if not short.empty else pd.DataFrame()
    short_v = np.array([float(s.loc[k, "hbv_short_molecules"])
                        if (len(s) and k in s.index) else 0.0 for k in keys])

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    ax[0].bar(x, long_v, 0.6, label="long fraction (splice-aligned)", color="#1f77b4")
    ax[0].bar(x, short_v, 0.6, bottom=long_v, label="short fraction (sr-aligned)",
              color="#e6a817")
    ax[0].set_yscale("log"); ax[0].set_ylabel("Unique HBV molecules")
    ax[0].set_title("HBV recovery per library (log scale)")
    ax[0].legend(fontsize=8, loc="upper right")
    for xi, l, sv in zip(x, long_v, short_v):
        ax[0].annotate(f"{int(l+sv):,}", (xi, l+sv), textcoords="offset points",
                       xytext=(0, 5), ha="center", fontsize=8.5)

    pct = [100*float(r.loc[k, "hbv_reads"])/float(r.loc[k, "composite_mapped"])
           if float(r.loc[k, "composite_mapped"]) else 0 for k in keys]
    cols = ["#2ca02c" if p > 50 else "#b0b0b0" for p in pct]
    ax[1].bar(x, pct, color=cols)
    ax[1].set_ylabel("HBV as % of mapped reads")
    ax[1].set_title("On-target fraction (green = high viral load)")
    for xi, p in zip(x, pct):
        ax[1].annotate(f"{p:.2f}%", (xi, p), textcoords="offset points",
                       xytext=(0, 4), ha="center", fontsize=8.5)
    if CARRYOVER_KEY in keys:
        i = keys.index(CARRYOVER_KEY)
        ax[1].plot([i-0.35, i+0.35], [CARRYOVER_PRIOR_PCT]*2, ls="--", color="#8b1a1a", lw=1.6)
        ax[1].annotate(f"{CARRYOVER_PRIOR_PCT}% in its own run", (i, CARRYOVER_PRIOR_PCT),
                       textcoords="offset points", xytext=(0, 6), ha="center",
                       fontsize=8, color="#8b1a1a")
    for a in ax:
        a.set_xticks(x); a.set_xticklabels(labs, fontsize=8)
    fig.tight_layout()
    return fig


def fig_short_fraction(short, rec, meta):
    """How much of each library is short, and how HBV-depleted that fraction is."""
    if short.empty:
        return None
    s = short.set_index("sample_key")
    keys = [k for k in DISPLAY_ORDER if k in s.index]
    if not keys:
        return None
    labs = [label(k, meta) for k in keys]
    x = np.arange(len(keys))
    r = rec.set_index("sample_key") if not rec.empty else pd.DataFrame()

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    short_hbv_pct, long_hbv_pct = [], []
    for k in keys:
        py = float(s.loc[k, "pychopper_out"]) or 1
        short_hbv_pct.append(100*float(s.loc[k, "hbv_short_aligned"])/py)
        if len(r) and k in r.index and float(r.loc[k, "composite_mapped"]):
            long_hbv_pct.append(100*float(r.loc[k, "hbv_reads"])/float(r.loc[k, "composite_mapped"]))
        else:
            long_hbv_pct.append(0.0)
    w = 0.38
    ax[0].bar(x-w/2, long_hbv_pct, w, label="long fraction", color="#1f77b4")
    ax[0].bar(x+w/2, short_hbv_pct, w, label="short fraction", color="#e6a817")
    ax[0].set_yscale("log"); ax[0].set_ylabel("HBV as % of aligned reads")
    ax[0].set_title("HBV is strongly depleted in the short fraction")
    ax[0].legend(fontsize=8)

    dup_s = [float(s.loc[k, "short_dup_rate"]) if str(s.loc[k, "short_dup_rate"]) != "NA" else np.nan
             for k in keys]
    dup_l = [float(r.loc[k, "dup_rate"]) if (len(r) and k in r.index) else np.nan for k in keys]
    ax[1].bar(x-w/2, [d*100 for d in dup_l], w, label="long fraction", color="#1f77b4")
    ax[1].bar(x+w/2, [d*100 for d in dup_s], w, label="short fraction", color="#e6a817")
    ax[1].set_ylabel("duplication rate (%)"); ax[1].set_ylim(0, 105)
    ax[1].set_title("The short fraction was never capture-enriched,\nso it is far less duplicated")
    ax[1].legend(fontsize=8)
    for a in ax:
        a.set_xticks(x); a.set_xticklabels(labs, fontsize=8)
    fig.tight_layout()
    return fig


def fig_concatemers(p2, qc, meta):
    """Homologous concatemer rate per library, and its relationship to read length.

    Concatemers are two or more molecules from the SAME library ligated together. They
    carry a consistent barcode at both ends, so both-ends demultiplexing does not remove
    them — only an internal-adapter scan does. Rate is measured on HBV-aligned reads,
    before deduplication."""
    if p2.empty or "hbv_concatemer_pct" not in p2.columns:
        return None
    p = p2.set_index("sample_key")
    keys = [k for k in DISPLAY_ORDER if k in p.index]
    if not keys:
        return None
    labs = [label(k, meta) for k in keys]
    x = np.arange(len(keys))
    pct = [float(p.loc[k, "hbv_concatemer_pct"]) for k in keys]
    n = [int(p.loc[k, "hbv_concatemers_removed"]) for k in keys]

    fig, ax = plt.subplots(1, 2, figsize=(13.5, 4.6))
    ax[0].bar(x, pct, color="#8172b3")
    ax[0].set_ylabel("% of HBV reads that are concatemers")
    ax[0].set_title("cDNA005 — homologous concatemer rate\n(removed before deduplication)")
    ax[0].set_ylim(0, max(pct) * 1.25 if pct else 1)
    for xi, v, c in zip(x, pct, n):
        ax[0].annotate(f"{v:.2f}%\n({c:,})", (xi, v), textcoords="offset points",
                       xytext=(0, 4), ha="center", fontsize=8)
    ax[0].set_xticks(x); ax[0].set_xticklabels(labs, fontsize=8)

    # cross-experiment comparison
    c3_lab = [f"{inp}" for _, inp, _, _, _ in CDNA003_CONCAT]
    c3_pct = [p for _, _, p, _, _ in CDNA003_CONCAT]
    c5_lab = [labs[i].replace("\n", " ") for i in range(len(keys))]
    allv = c3_pct + pct
    x3 = np.arange(len(c3_pct))
    x5 = np.arange(len(pct)) + len(c3_pct) + 0.8
    ax[1].bar(x3, c3_pct, color="#4c72b0", label="cDNA003 (EXP26000896)")
    ax[1].bar(x5, pct, color="#8172b3", label="cDNA005 (EXP26000993)")
    for xi, v in zip(np.concatenate([x3, x5]), allv):
        ax[1].annotate(f"{v:.2f}", (xi, v), textcoords="offset points", xytext=(0, 3),
                       ha="center", fontsize=7.5)
    # mark the one library present in both runs
    try:
        i3 = [k for k, _, _, _, _ in CDNA003_CONCAT].index("SeqLib5552_1ng_polyA_21")
        i5 = keys.index(CARRYOVER_KEY)
        ax[1].annotate("", xy=(x3[i3], c3_pct[i3] + 0.9), xytext=(x5[i5], pct[i5] + 0.9),
                       arrowprops=dict(arrowstyle="<->", color="#8b1a1a", lw=1.4))
        ax[1].annotate("same physical library\n0.97% → 5.43%",
                       ((x3[i3] + x5[i5]) / 2, max(c3_pct[i3], pct[i5]) + 1.6),
                       ha="center", fontsize=8, color="#8b1a1a", fontweight="bold")
    except (ValueError, IndexError):
        pass
    ax[1].set_xticks(np.concatenate([x3, x5]))
    ax[1].set_xticklabels(c3_lab + c5_lab, fontsize=7.5, rotation=30, ha="right")
    ax[1].set_ylabel("% concatemers")
    ax[1].set_ylim(0, max(allv) * 1.3)
    ax[1].set_title("Same filter, same settings, both experiments")
    ax[1].legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    return fig


def fig_jackpot(conc, coords, meta):
    if conc.empty:
        return None
    c = conc.set_index("sample_key")
    keys = [k for k in DISPLAY_ORDER if k in c.index]
    if not keys:
        return None
    labs = [label(k, meta) for k in keys]
    x = np.arange(len(keys))
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))

    def lab(a, xi, yi, txt):
        a.annotate(txt, (xi, yi), textcoords="offset points", xytext=(0, 4),
                   ha="center", fontsize=8, clip_on=False)

    v0 = [c.loc[k, "top1pct_share"] for k in keys]
    ax[0].bar(x, v0, color="#d62728"); ax[0].axhline(1, ls=":", color="k", lw=1)
    ax[0].set_ylabel("% of reads held by top 1% of molecules")
    ax[0].set_title("Read concentration\n(dotted = uniform expectation)")
    ax[0].set_ylim(0, max(110, max(v0)*1.15))
    for xi, v in zip(x, v0): lab(ax[0], xi, v, f"{v:.0f}%")

    v1 = [c.loc[k, "mean"] for k in keys]
    ax[1].bar(x, v1, color="#8172b3")
    ax[1].set_ylabel("mean reads per molecule"); ax[1].set_title("Redundancy")
    ax[1].set_ylim(0, max(v1)*1.2)
    for xi, v in zip(x, v1): lab(ax[1], xi, v, f"{v:.1f}x")

    ts = [coords.get(k, {}).get("top5_share", np.nan) for k in keys]
    ax[2].bar(x, ts, color="#e6a817")
    ax[2].set_ylabel("% of molecules on the top 5 coordinates")
    ax[2].set_title("Coordinate concentration")
    if not all(np.isnan(ts)):
        ax[2].set_ylim(0, np.nanmax(ts)*1.2)
    for xi, v in zip(x, ts):
        if not np.isnan(v): lab(ax[2], xi, v, f"{v:.1f}%")
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
    fig, ax = plt.subplots(1, 2, figsize=(14, 5))
    for j, (logscale, ttl) in enumerate([(False, "all classes"),
                                         (True, "excluding unclassified (log scale)")]):
        bottom = np.zeros(len(keys))
        for cls in CLASS_ORDER:
            if logscale and cls == "unclassified":
                continue
            vals = []
            for k in keys:
                sub = cls_df[(cls_df.sample_key == k) & (cls_df.tx_class == cls)]
                tot = cls_df[cls_df.sample_key == k].n_reads.sum()
                vals.append(sub.n_reads.iloc[0]/tot*100 if len(sub) and tot else 0)
            vals = np.array(vals)
            if vals.sum() > 0:
                ax[j].bar(x, vals, bottom=bottom, label=clabel(cls) if j == 0 else None,
                          color=CLASS_COLORS[cls])
                bottom += vals
        ax[j].set_ylabel(f"% of classified {LEVEL_UNIT[level]}s")
        ax[j].set_title(ttl)
        ax[j].set_xticks(x); ax[j].set_xticklabels(labs, fontsize=8)
    ax[0].legend(fontsize=8, ncol=2, bbox_to_anchor=(1.02, 1), loc="upper left")
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
        ax.hist(df.tss, bins=np.arange(0, HBV_LEN+50, 50), histtype="step", lw=1.5,
                density=True, label=label(sk, meta).replace("\n", " "))
        plotted += 1
    if not plotted:
        plt.close(fig); return None
    bands = [("precore/pgRNA", "#1f77b4", [(1730, 1880)]), ("preS1", "#d62728", [(2700, 3100)]),
             ("X", "#e6a817", [(1260, 1450)]),
             ("preS2/S", "#2ca02c", [(3100, HBV_LEN), (0, 150)])]
    for _, c, spans in bands:
        for lo, hi in spans:
            ax.axvspan(lo, hi, color=c, alpha=0.08)
    ax.set_xlabel("5′ position on HBV genome (mod 3182); shaded = transcript-class TSS windows")
    ax.set_ylabel("density"); ax.set_title(f"5′ end distribution — {level}")
    leg = ax.legend(fontsize=8, loc="upper right", title="library"); ax.add_artist(leg)
    ax.legend(handles=[Patch(facecolor=c, alpha=0.25, label=n) for n, c, _ in bands],
              fontsize=7, loc="upper left", title="class window", framealpha=0.9)
    fig.tight_layout()
    return fig


SPAN_CLASS_ORDER = ["precore", "pgRNA", "preS1", "preS2_S", "X", "spliced",
                    "unclassified", "antisense", "pgRNA_RT"]
ANCHOR_X = 1936 + HBV_LEN


def fig_spans(meta, level, max_per_panel=3000):
    panels = []
    for sk in DISPLAY_ORDER:
        d = _classified(sk, level)
        if d is None or d.empty:
            continue
        if len(d) > max_per_panel:
            d = d.sample(max_per_panel, random_state=0)
        panels.append((label(sk, meta).replace("\n", " "), d))
    if not panels:
        return None
    PER, MINH, MAXH = 0.05, 1.4, 45.0
    hs = [min(MAXH, max(MINH, len(d)*PER + 0.5)) for _, d in panels]
    fig, axes = plt.subplots(len(panels), 1, figsize=(11, sum(hs)+0.4), sharex=True,
                             gridspec_kw={"height_ratios": hs, "hspace": 0.05})
    if len(panels) == 1:
        axes = [axes]
    ordk = {c: i for i, c in enumerate(SPAN_CLASS_ORDER)}

    def parse_blocks(bs, rs, re_):
        if isinstance(bs, str) and bs:
            try:
                return [tuple(int(x) for x in t.split("-")) for t in bs.split(";")]
            except ValueError:
                pass
        return [(rs, re_)]

    xlo = xhi = ANCHOR_X
    present = set()
    for ax, (ylabel, df) in zip(axes, panels):
        d = df.copy()
        d["_o"] = d.tx_class.map(lambda c: ordk.get(c, 99))
        d["_s"] = d.ref_end - d.ref_start
        d = d.sort_values(["_o", "_s"], ascending=[True, False]).reset_index(drop=True)
        has_b = "blocks" in d.columns
        for i, row in enumerate(d.itertuples()):
            col = CLASS_COLORS.get(row.tx_class, "#999")
            shift = (row.polya_end + HBV_LEN) - row.ref_end
            exons = parse_blocks(getattr(row, "blocks", None) if has_b else None,
                                 row.ref_start, row.ref_end)
            xs = [s+shift for s, _ in exons]; xe = [e+shift for _, e in exons]
            if len(exons) > 1:
                ax.plot([min(xs), max(xe)], [i, i], color=col, lw=0.4, alpha=0.5, zorder=1)
            for s, e in zip(xs, xe):
                ax.plot([s, e], [i, i], color=col, lw=0.9, solid_capstyle="butt", zorder=2)
            xlo, xhi = min(xlo, min(xs)), max(xhi, max(xe))
        present.update(d.tx_class)
        ax.axvline(HBV_LEN, ls="--", color="k", lw=0.8)
        ax.axvline(ANCHOR_X, ls=":", color="darkred", lw=1.0)
        ax.set_ylim(len(d)-0.5, -0.5); ax.set_yticks([])
        ax.set_ylabel(f"{ylabel}\nn={len(d)}", rotation=0, ha="right", va="center", fontsize=8)
    for ax in axes:
        ax.set_xlim(xlo-80, xhi+120)
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


def build(level, meta, qc, rec, cls_df, short, coords, conc, p2, today):
    r = rec.set_index("sample_key") if not rec.empty else pd.DataFrame()
    s = short.set_index("sample_key") if not short.empty else pd.DataFrame()
    unit = LEVEL_UNIT[level]

    def g(sk, col, default=0):
        try: return r.loc[sk, col]
        except Exception: return default

    def gs(sk, col, default=0):
        try: return s.loc[sk, col]
        except Exception: return default

    keys = ordered(r.index) if len(r) else DISPLAY_ORDER

    warn = ('<div class="warn"><b>Counting level: reads (undeduplicated).</b> '
            'Reported for comparability with earlier work. Use the molecule-level report '
            'for interpretation.</div>') if level == "reads" else \
           ('<div class="warn"><b>Counting level: molecules (UMI-deduplicated).</b> '
            'Each row is one original cDNA molecule. This is the level to interpret.</div>')

    rows_html = ""
    for k in keys:
        lf = int(g(k, "hbv_unique")); sf = int(gs(k, "hbv_short_molecules"))
        cm = float(g(k, "composite_mapped")) or 1
        rows_html += (f"<tr><td>{NICE.get(k,k).replace(chr(10),' ')}</td>"
                      f"<td>{int(g(k,'hbv_reads')):,}</td>"
                      f"<td>{100*float(g(k,'hbv_reads'))/cm:.2f}%</td>"
                      f"<td>{lf:,}</td><td>{sf:,}</td><td><b>{lf+sf:,}</b></td>"
                      f"<td>{float(g(k,'dup_rate',0))*100:.1f}%</td></tr>")

    p2i = p2.set_index("sample_key") if not p2.empty else pd.DataFrame()
    qci = qc.set_index("sample_key") if not qc.empty else pd.DataFrame()
    concat_rows = ""
    for sk, inp, p, tot, nc in CDNA003_CONCAT:
        star = ' <b>&#9733;</b>' if sk == "SeqLib5552_1ng_polyA_21" else ""
        concat_rows += (f'<tr style="background:#f7f9fc"><td>cDNA003</td>'
                        f'<td>{sk.split("_")[0]} ({inp}){star}</td><td>{tot:,}</td>'
                        f'<td>{nc:,}</td><td>{p:.2f}%</td>'
                        f'<td>no — filter post-dates that analysis</td></tr>')
    n_c5 = 0
    for k in keys:
        if len(p2i) and k in p2i.index and "hbv_concatemer_pct" in p2i.columns:
            n_c5 += 1
            star = ' <b>&#9733;</b>' if k == CARRYOVER_KEY else ""
            tot = p2i.loc[k, "hbv_reads_pre_filter"] if "hbv_reads_pre_filter" in p2i.columns else None
            tot_s = f"{int(tot):,}" if tot is not None and str(tot) != "nan" else "—"
            concat_rows += (f"<tr><td>cDNA005</td>"
                            f"<td>{NICE.get(k,k).replace(chr(10),' ')}{star}</td>"
                            f"<td>{tot_s}</td>"
                            f"<td>{int(p2i.loc[k,'hbv_concatemers_removed']):,}</td>"
                            f"<td>{float(p2i.loc[k,'hbv_concatemer_pct']):.2f}%</td>"
                            f"<td><b>yes</b></td></tr>")
    if not n_c5:
        concat_rows += ('<tr><td>cDNA005</td><td colspan="5"><em>phase 2 summary not found '
                        'in analysis/comparison/ — rerun to populate</em></td></tr>')
    concat_rows += ('<tr><td colspan="6" style="text-align:left;font-size:.8rem;'
                    'border:none;padding-top:6px">&#9733; the same physical library, '
                    'measured in both runs</td></tr>')

    figs = {
        "qc": img(fig_qc(qc, meta),
                  "Figure 1. Per-library read QC. The dashed line in the middle panel marks "
                  "218 bp — the length of the adapter, barcode, SSPII and CRTA construct with "
                  "no insert at all. Read lengths close to it mean very little biological "
                  "sequence per molecule."),
        "rec": img(fig_recovery_by_sample(rec, short, meta),
                   "Figure 2. Left: unique HBV molecules per library, long fraction (blue) plus "
                   "short fraction recovered by short-read realignment (orange), log scale. "
                   "Right: HBV as a percentage of mapped reads. The dashed red line on the "
                   f"carryover library marks the {CARRYOVER_PRIOR_PCT}% it achieved in its own "
                   "experiment — the same physical library, showing the effect of being pooled "
                   "with high-viral-load samples."),
        "short": img(fig_short_fraction(short, rec, meta),
                     "Figure 3. Left: HBV as a share of aligned reads in each fraction, log "
                     "scale. The short fraction is 40-600x depleted of HBV — fragments below "
                     "probe length are not captured. Right: duplication rate. The short "
                     "fraction was never enriched, so it was amplified far less."),
        "conc": img(fig_concatemers(p2, qc, meta),
                    "Figure 4. Homologous concatemer rate. Left: percentage of HBV-aligned "
                    "reads carrying an adapter or barcode motif more than 200 bp from both "
                    "ends, which is evidence the read is two or more molecules ligated "
                    "together. Because both fragments come from the same library they carry "
                    "the same barcode at both ends, so both-ends demultiplexing cannot "
                    "remove them. Right: the rate against read length N50."),
        "jack": img(fig_jackpot(conc, coords, meta),
                    "Figure 5. PCR jackpotting. Left: share of reads held by the top 1% of "
                    "molecules (uniform duplication gives ~1%). Middle: mean reads per "
                    "molecule. Right: share of molecules on the top five alignment "
                    "coordinates. Note the two low-load libraries are too shallowly sequenced "
                    "for these metrics to be meaningful."),
        "cls": img(fig_class_composition(cls_df, meta, level),
                   f"Figure 6. Transcript-class composition ({unit} level). Left: all classes — "
                   "unclassified dominates because reads are too short to satisfy the span "
                   "gates. Right: the same data excluding unclassified, log scale, so the "
                   "assigned classes are visible."),
        "tss": img(fig_tss(meta, level),
                   "Figure 7. 5′ end distribution, density-normalised so libraries of very "
                   "different size are comparable. Shaded bands are the transcript-class TSS "
                   "windows."),
        "spans": img(fig_spans(meta, level),
                     f"Figure 8. {unit.capitalize()} spans, 3′-anchored at the canonical poly-A. "
                     "Colour = class, grouped by class and ordered longest to shortest. Panels "
                     "are subsampled to 3,000."),
    }

    html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>EXP26000993 cDNA005 serum/plasma enrichment — {level}</title><style>{CSS}</style></head><body>
<h1>HBV enrichment from serum and plasma RNA — EXP26000993_cDNA005</h1>
<p><em>Generated {today} · pipeline {VERSION} · counting level: <b>{level}</b></em></p>

{warn}

<div class="key">
<b>Executive summary.</b> Four libraries at 1 ng input, pooled at equal mass (375 ng each)
and enriched together by Twist hybridisation capture: a human plasma pool from three HBV
patients, two Yecuris chimeric mice (g21-treated and PBS control), and SeqLib5552 carried
over from cDNA003 as a cross-run control. Unlike cDNA003 this is a <em>sample type</em>
comparison, not an input titration.
<ul>
<li><b>Where viral load was high, capture worked very well</b> — HBV is 86–95% of mapped
reads in the two Yecuris libraries.</li>
<li><b>High and low viral load samples should not be pooled together.</b> The Yecuris
libraries consumed the capture capacity. The carryover library — physically the same
material that ran in cDNA003 — gave 0.27% HBV here against {CARRYOVER_PRIOR_PCT}% there,
a ~10× suppression. Plasma reached 1.9%. HBV was detected in both low-load libraries, but
their yields understate what a balanced pool would give.</li>
<li><b>5.6× fewer HBV molecules in the treated animal</b>, consistent with the known viral
load difference. Both Yecuris libraries are ~90% HBV after capture, so the difference
appears as total reads recovered rather than as on-target percentage.</li>
<li><b>Read length limits transcript classification.</b> Median 250–350 bp in the
serum/plasma libraries against a 218 bp adapter construct. ~98% of HBV molecules fail the
span gates and come back unclassified — abundance is measurable, isoform usage is not.</li>
<li><b>Homologous concatemers range from 1.1% to 17.3% of HBV reads</b> and were removed
before deduplication, so nothing downstream includes them. cDNA003, measured retrospectively
with the same filter, sits at 0.63–0.97%. The same physical library reads 0.97% there and
5.43% here, so the rate is a property of the run rather than of the library.</li>
<li><b>Jackpotting is much improved where depth allows judging it</b> — 10–12× redundancy
in the Yecuris libraries against up to 387× in cDNA003. High viral load means high
complexity, so no single molecule dominates. The two low-load libraries are too shallow to
assess.</li>
</ul>
</div>

<h2>Methods</h2>
<p><b>Samples.</b> Human plasma pool (3 HBV patients), Yecuris chimeric mouse serum
(g21-treated and PBS control), and SeqLib5552 from cDNA003. All at 1 ng input. The three
new samples received <b>enzymatic poly(A) tailing</b> (PAP) because serum and plasma RNA is
fragmented and largely lacks native tails; the carryover library came from oligo-d(T)
poly(A)-selected intact cellular RNA. Those two routes are not equivalent and the
comparison between them is confounded on that axis.</p>
<p><b>Barcoding and demultiplexing.</b> The same four custom 24-nt PCR barcodes as cDNA003.
Demultiplexed with <code>dorado demux</code> using a custom arrangement requiring the same
barcode at both read ends.</p>
<p><b>Reference and alignment.</b> hg38 + 2×HBV (U95551.1 doubled to 6,364 bp so
genome-wrapping transcripts align contiguously). minimap2 <code>-ax splice</code>. Note the
Yecuris samples are chimeric <em>mice</em> and no mouse genome is present in the reference,
so a large share of their non-HBV reads cannot map.</p>
<p><b>Short fraction.</b> Reads that fail splice alignment — 56–91% of the serum/plasma
libraries — were realigned with <code>minimap2 -ax sr</code> against a <b>single-copy</b>
HBV reference (doubling gains nothing for a fragment that cannot wrap, and costs unique
mapping). Those molecules are reported separately and added to the totals.</p>
<p><b>UMIs and deduplication.</b> pychopper (<code>-k PCB114 -U</code>) to orient reads and
extract the 28-nt structured UMI, then <code>umi_tools dedup --method=directional</code> on
well-formed 28-nt UMIs. Homologous concatemers — fusions of two molecules from the same
library, which carry a consistent barcode at both ends and survive both-ends
demultiplexing — were filtered before deduplication (1.1% of HBV reads in the Yecuris
libraries, 5.4% in the carryover, 17.3% among the long reads of the plasma library).</p>
<p><b>Statistics.</b> Poisson 95% CIs on counts, Wilson 95% CIs on proportions.</p>

<h2>Library QC</h2>
{figs['qc']}

<h2>HBV recovery</h2>
<table class="data"><tr><th>library</th><th>HBV reads</th><th>% of mapped</th>
<th>molecules (long)</th><th>molecules (short)</th><th>total</th><th>duplication</th></tr>
{rows_html}</table>
{figs['rec']}
{figs['short']}
<p>The short fraction adds 2–8% to molecule counts and is 40–600× depleted of HBV relative
to the long fraction, which is what you would expect if fragments shorter than the ~120 nt
capture probes cannot form a stable enough duplex to be pulled down. Including it moves the
treated/control ratio from 5.46× to 5.63×, so excluding it was not biasing the comparison.</p>

<h2>Homologous concatemers</h2>

<h3>How they are identified</h3>
<p>A concatemer is two or more cDNA molecules ligated into a single read. When the two
fragments come from <em>different</em> libraries the read carries disagreeing barcodes at
its two ends and <code>dorado demux --barcode-both-ends</code> rejects it. When both come
from the <em>same</em> library — a <b>homologous</b> concatemer — the outer barcodes agree,
the read is assigned normally, and nothing downstream removes it.</p>
<p>They are detected by the one feature that cannot occur in an intact molecule: the library
construct places adapter, SSPII, CRTA and barcode sequences only at the termini, so any of
those motifs found <b>more than 200 bp from both ends</b> of a read marks a fusion junction.
Matching is by edit distance (adapters tolerate length/6, barcodes up to 8) because ONT
error rates make exact matching of a 27 bp motif fail more often than it succeeds. The scan
runs immediately before pychopper, which is the only place it can run — pychopper trims the
primers, so downstream of it the adapter sequence is gone whether or not the read was a
fusion.</p>
<p>Two consequences of leaving them in. A fused pair carries two UMIs, so deduplication sees
only whichever pychopper extracted and the second molecule is silently lost. And on a 2×
reference a fused pair of HBV molecules can align contiguously across the copy junction and
satisfy the ≥3,982 bp <code>pgRNA_RT</code> span gate — a PCR artifact scored as
tandem/concatemeric readthrough, i.e. as biology.</p>
<p>The reported rate is a <b>floor</b>. A fusion whose junction happens to fall within 200 bp
of either read end is invisible to the scan, as is one whose internal adapter has degraded
past the edit-distance threshold.</p>

<h3>Quantification and disposition</h3>
{figs['conc']}
<table class="data"><tr><th>experiment</th><th>library</th><th>HBV reads scanned</th>
<th>concatemers</th><th>% of HBV reads</th><th>excluded downstream?</th></tr>
{concat_rows}</table>
<p><b>In cDNA005 these reads were removed before deduplication and are excluded from every
downstream analysis in this report</b> — molecule counts, transcript classification, span
maps and the recovery table are all computed on the filtered set. <b>In cDNA003 they were
not.</b> The filter was written for this experiment; the cDNA003 figures above were measured
retrospectively, by regenerating HBV reads from its alignments and running the identical
filter at identical settings. Those 0.63–0.97% of reads were present in the cDNA003 analysis.
At under 1% the effect on its abundance estimates is negligible, but its
<code>pgRNA_RT</code> counts specifically should be treated as an upper bound, since that is
the class concatemers are most likely to be miscalled into.</p>
<p>The method was validated by round-tripping: regenerating cDNA005's HBV reads from its
alignments and re-running the filter reproduced the phase-2 numbers to within 0.03
percentage points (5.43 vs 5.4, 17.32 vs 17.3, 1.07 vs 1.1), confirming that the two
experiments are measured the same way and not merely with the same tool.</p>

<h3>What the comparison shows</h3>
<p>cDNA003 sits in a narrow band of 0.63–0.97% across a thousandfold input range, with no
trend in input amount. cDNA005 spans 1.07% to 17.32%. Read length does not explain the
difference — cDNA003's reads are considerably longer, and a longer read has <em>more</em>
room to be a detectable fusion, so if anything its rates are measured more sensitively.</p>
<p>The most informative single number is that <b>SeqLib5552 measures 0.97% in cDNA003 and
5.43% in cDNA005</b>. That is the same physical library scanned by the same code, so the
concatemer rate is not a fixed property of a library — it depends on the run it is measured
in. Two candidate explanations, not resolved here: the cDNA005 pool was dominated by the
high-load Yecuris libraries, so only ~3% as many HBV reads survived for SeqLib5552, and
concatemers are longer molecules that present more probe-binding surface and would be
preferentially retained when probes are limiting; alternatively the carryover material saw
additional handling between runs. Either way the practical reading is that concatemer rate
should be measured per run and not assumed to carry over.</p>

<h2>Amplification jackpotting</h2>
{figs['jack']}
<p>Redundancy of 10–12× in the Yecuris libraries is a different regime from cDNA003, where
the lowest-input library reached 387× and a single molecule held 440,452 reads. High viral
load produces high molecular complexity, so no molecule can dominate. The two low-load
libraries here sit at 2–5× redundancy, which is too shallow for these metrics to mean
anything — the carryover library shows this directly, since the same physical material
showed 90× redundancy and 82% of reads in its top 1% of molecules when sequenced deeply in
cDNA003.</p>

<h2>Transcript classification</h2>
{figs['cls']}
{figs['tss']}
<p>Roughly 98% of HBV molecules in the serum/plasma libraries are unclassified. Every
transcript class is gated on span — X requires 300–1,000 bp, preS2/S ≥1,200,
precore/pgRNA ≥2,600 — and a molecule with 30–130 bp of insert satisfies none of them. The
untreated Yecuris library still yields roughly 36,000 classified molecules in absolute
terms, which is ample for statistics, but they are the longest ~1.5% and therefore do not
represent the population. The full-length classes stay thin: 180 pgRNA and 4 precore.</p>
<p>The classifier used here is a reimplementation built from written specification, the
canonical module from cDNA001 not being available. Its logic was verified against a 30-case
test set, but class counts should be confirmed against the original before numbers are
compared across experiments.</p>

<h2>Read/molecule span maps</h2>
{figs['spans']}

<h2>Conclusions</h2>
<ol>
<li>Hybridisation capture is highly effective where viral load is high — 86–95% of mapped
reads on target in both Yecuris libraries.</li>
<li><b>Do not pool high and low viral load samples in one capture.</b> Probe capacity is
finite, and the carryover control quantifies the cost: 0.27% HBV against
{CARRYOVER_PRIOR_PCT}% for the same library in a balanced pool.</li>
<li>The treated animal yields 5.6× fewer HBV molecules than the PBS control, consistent
with the known viral load difference.</li>
<li>Read length, not depth or capture, is what prevents isoform-level analysis of serum and
plasma RNA. This is upstream of the library prep — it is a property of the input material.</li>
<li>The short fraction is real RNA but largely un-capturable and contributes only 2–8% of
molecules. It does not change any conclusion here, but it is a substantial share of the
sequencing spend.</li>
<li>Homologous concatemers are invisible to both-ends demultiplexing and require the
internal-adapter filter, now in the pipeline and applied here. They reach 17.3% of HBV reads
in the plasma library against 0.63–0.97% throughout cDNA003, and the same library measures
5.6× higher here than there — so the rate should be measured in every run rather than assumed
from a previous one.</li>
</ol>

<h2>Appendix — reproducibility</h2>
<p>dorado demux (custom arrangement, both-ends), minimap2, samtools, pychopper, umi_tools,
pandas/scipy/statsmodels. Pipeline: <code>run_dorado_demux.sh</code> →
<code>reshape_demux_output.sh</code> → <code>phase1_align.sh</code> →
<code>phase1_5_readqc_demux.sh</code> → <code>phase2_hbv_umi_v2.sh</code> →
<code>phase2s_short_fraction.sh</code> → <code>phase3_classify.py</code> /
<code>phase3b_classify_molecules.py</code> → <code>phase4_quantify.py</code> →
<code>phase5_report_EXP26000993.py</code>. Tables in <code>analysis/comparison/</code>.</p>
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
    short = load("phase2s_short_fraction_summary.tsv")
    p2 = pd.DataFrame()
    for cand in ("phase2_hbv_umi_summary.tsv", "phase2_summary.tsv",
                 "phase2_hbv_umi.tsv", "phase2_hbv_summary.tsv"):
        p2 = load(cand)
        if not p2.empty:
            print(f"  concatemer counts from {cand}")
            break
    if p2.empty:
        sys.stderr.write("WARN: no phase 2 summary found in analysis/comparison/ — "
                         "concatemer section will be empty. Check the filename:\n"
                         f"       ls {COMP}\n")
    out_dir = ANALYSIS / "reports"; out_dir.mkdir(parents=True, exist_ok=True)

    print("  computing duplicate concentration from hbv.umi.bam (always recomputed)...")
    conc = umi_concentration()
    if not conc.empty:
        print(f"    {len(conc)} libraries -> {COMP/'umi_concentration.tsv'}")
    coords = coord_counts("molecules")

    for lv in (["reads", "molecules"] if args.level == "both" else [args.level]):
        rec = load(f"phase4_recovery_{lv}.tsv")
        if rec.empty:
            rec = load("phase4_recovery.tsv")
        cls_df = load(f"phase4_class_abundance_{lv}.tsv")
        if cls_df.empty:
            cls_df = load("phase4_class_abundance.tsv")
        if rec.empty:
            sys.stderr.write(f"WARN [{lv}]: no recovery table; skipping\n"); continue
        html = build(lv, meta, qc, rec, cls_df, short, coords, conc, p2, today)
        out = out_dir / f"EXP26000993_report_{lv}.html"
        out.write_text(html)
        print(f"  {lv:<10} -> {out}  ({out.stat().st_size/1024:.0f} KB)")


if __name__ == "__main__":
    main()
