#!/usr/bin/env python3
"""
Phase 5 — self-contained HTML report for the cDNA HBV baseline (EXP26000559).

Manuscript style, modelled on the EXP26000465 report but adapted:
  - 8-condition library-prep matrix (input × polyA × PCR cycles), not treated/untreated
  - unique HBV molecules (UMI) alongside raw reads; PCR-duplication story
  - transcript classification on 2× coords (wrap kept visible), no m6A (cDNA)

Figures are matplotlib PNGs inlined as base64. Reads inputs from analysis/comparison/
and per-sample hbv_classified.parquet. Output: analysis/reports/EXP26000559_report.html

Usage:  conda activate hbv_lr_analysis; python scripts/phase5_report.py
"""
import base64, datetime, io
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

PROJECT_ROOT = Path("/home/ubuntu/matt_wolpert_claude_code_analysis/2026_07_02_EXP26000559_cDNA001")
ANALYSIS = PROJECT_ROOT / "analysis"
COMP = ANALYSIS / "comparison"
SAMPLES = PROJECT_ROOT / "config" / "samples.tsv"
OUT = ANALYSIS / "reports" / "EXP26000559_report.html"
VERSION = "phase5-v1"
HBV_LEN = 3182

CLASS_COLORS = {
    "preS2_S": "#2ca02c", "preS1": "#d62728", "precore": "#e377c2", "pgRNA": "#1f77b4",
    "pgRNA_RT": "#17becf", "X": "#e6a817", "spliced": "#9467bd",
    "unclassified": "#9e9e9e", "antisense": "#cfcfcf",
}
CLASS_ORDER = list(CLASS_COLORS)
# Display labels (internal key → figure label). Spliced reads are the canonical SP1 variant.
CLASS_LABEL = {"preS2_S": "preS2/S", "spliced": "spliced (SP1)"}
def clabel(c):
    return CLASS_LABEL.get(c, c)

# Display order: NOpolyA then polyA, 150ng then 10ng, 20cyc then 16cyc
DISPLAY_ORDER = [
    "SeqLib5543_150ng_NOpolyA_20", "SeqLib5540_150ng_NOpolyA_16",
    "SeqLib5544_10ng_NOpolyA_20", "SeqLib5541_10ng_NOpolyA_16",
    "SeqLib5542_150ng_polyA_20",  "SeqLib5539_150ng_PolyA_16",
    "NoName1_10ng_PolyA_16",      "NoName2_10ng_PolyA_20",
]


def b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def img(fig, alt):
    return f'<figure><img alt="{alt}" src="data:image/png;base64,{b64(fig)}"/><figcaption>{alt}</figcaption></figure>'


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
    df = pd.DataFrame(rows).set_index("sample_key")
    return df


def label(sk, meta):
    m = meta.loc[sk]
    pa = "polyA" if m.polya == "yes" else "NOpolyA"
    return f"{m.input_ng}ng {pa} {m.pcr}c\n{m.barcode}"


def ordered(keys):
    return [k for k in DISPLAY_ORDER if k in set(keys)]


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def fig_qc(qc, meta):
    q = qc[qc.sample_key.isin(meta.index)].set_index("sample_key")
    keys = ordered(q.index)
    labs = [label(k, meta) for k in keys]
    x = np.arange(len(keys))
    fig, ax = plt.subplots(2, 2, figsize=(12, 7))
    ax[0, 0].bar(x, [q.loc[k, "reads"] / 1e6 for k in keys], color="#4c72b0")
    ax[0, 0].set_ylabel("Pass reads (millions)"); ax[0, 0].set_yscale("log"); ax[0, 0].set_title("Yield")
    ax[0, 1].bar(x, [q.loc[k, "len_N50"] for k in keys], color="#55a868")
    ax[0, 1].set_ylabel("Read length N50 (bp)"); ax[0, 1].set_title("Read length")
    ax[1, 0].bar(x, [q.loc[k, "mean_qscore"] for k in keys], color="#c44e52")
    ax[1, 0].set_ylabel("Mean Q"); ax[1, 0].set_title("Quality")
    ax[1, 1].bar(x, [q.loc[k, "frac_with_polyA"] * 100 for k in keys], color="#8172b3")
    ax[1, 1].set_ylabel("% reads with poly(A) tail"); ax[1, 1].set_title("Poly(A) detection")
    for a in ax.flat:
        a.set_xticks(x); a.set_xticklabels(labs, rotation=45, ha="right", fontsize=7)
    fig.tight_layout()
    return fig


def fig_recovery(rec, meta):
    r = rec.set_index("sample_key")
    keys = ordered(r.index)
    labs = [label(k, meta) for k in keys]
    x = np.arange(len(keys))
    def ci_arr(col_ci, point):
        lo, hi = [], []
        for k in keys:
            s = r.loc[k, col_ci].strip("[]").split(",")
            lo.append(r.loc[k, point] - float(s[0])); hi.append(float(s[1]) - r.loc[k, point])
        return np.array([lo, hi])
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    uv = [r.loc[k, "hbv_unique"] for k in keys]
    rv = [r.loc[k, "hbv_reads"] for k in keys]
    ax[0].bar(x - 0.2, rv, 0.4, label="raw HBV reads", color="#b0b0b0")
    ax[0].bar(x + 0.2, uv, 0.4, yerr=ci_arr("hbv_unique_ci", "hbv_unique"),
              capsize=3, label="unique molecules", color="#1f77b4")
    ax[0].set_ylabel("HBV reads / unique molecules")
    ax[0].set_title("Absolute HBV recovery (Poisson 95% CI on unique)")
    ax[0].legend(fontsize=8)
    pm = [r.loc[k, "hbv_unique_per_M"] for k in keys]
    ax[1].bar(x, pm, color=["#2ca02c" if meta.loc[k, "polya"] == "yes" else "#b0b0b0" for k in keys])
    ax[1].set_ylabel("Unique HBV per million composite-mapped")
    ax[1].set_title("Efficiency (green = polyA-selected)")
    for a in ax:
        a.set_xticks(x); a.set_xticklabels(labs, rotation=45, ha="right", fontsize=7)
    fig.tight_layout()
    return fig


def fig_dup(rec, meta):
    r = rec.set_index("sample_key")
    keys = [k for k in ordered(r.index) if str(r.loc[k, "dup_rate"]) not in ("NA", "nan")]
    labs = [label(k, meta) for k in keys]
    x = np.arange(len(keys))
    vals = [float(r.loc[k, "dup_rate"]) * 100 for k in keys]
    col = ["#d62728" if meta.loc[k, "pcr"] == "20" else "#4c72b0" for k in keys]
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.bar(x, vals, color=col)
    ax.set_ylabel("PCR duplication rate (%)")
    ax.set_title("UMI duplication (red = 20 cycles, blue = 16 cycles)")
    ax.set_xticks(x); ax.set_xticklabels(labs, rotation=45, ha="right", fontsize=7)
    for xi, v in zip(x, vals):
        ax.text(xi, v + 1, f"{v:.0f}%", ha="center", fontsize=8)
    fig.tight_layout()
    return fig


def fig_class_composition(cls_df, meta):
    keys = ordered(cls_df.sample_key.unique())
    labs = [label(k, meta) for k in keys]
    x = np.arange(len(keys))
    fig, ax = plt.subplots(figsize=(11, 5))
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
    ax.set_ylabel("% of unique HBV molecules")
    ax.set_title("Transcript class composition per condition (UMI-deduplicated molecules)")
    ax.set_xticks(x); ax.set_xticklabels(labs, rotation=45, ha="right", fontsize=7)
    ax.legend(fontsize=8, ncol=2, bbox_to_anchor=(1.01, 1), loc="upper left")
    fig.tight_layout()
    return fig


def fig_class_composition_combined_150(cls_df, meta):
    """Fig 4b — 150 ng only, polyA + NOpolyA pooled within each PCR-cycle level."""
    groups = groups_150ng_bypcr(meta)
    labs = [g[0].replace("\n", " ") for g in groups]
    x = np.arange(len(groups))
    fig, ax = plt.subplots(figsize=(7, 5))
    bottom = np.zeros(len(groups))
    for cls in CLASS_ORDER:
        vals = []
        for _, sks in groups:
            sub = cls_df[(cls_df.sample_key.isin(sks)) & (cls_df.tx_class == cls)].n_reads.sum()
            tot = cls_df[cls_df.sample_key.isin(sks)].n_reads.sum()
            vals.append(sub / tot * 100 if tot else 0)
        vals = np.array(vals)
        if vals.sum() > 0:
            ax.bar(x, vals, bottom=bottom, label=clabel(cls), color=CLASS_COLORS[cls])
            bottom += vals
    ax.set_ylabel("% of unique HBV molecules")
    ax.set_title("150 ng transcript class composition — UMI-deduplicated (polyA + NOpolyA combined, per PCR level)")
    ax.set_xticks(x); ax.set_xticklabels(labs, fontsize=8)
    ax.legend(fontsize=8, ncol=2, bbox_to_anchor=(1.01, 1), loc="upper left")
    fig.tight_layout()
    return fig


# Top-to-bottom class grouping for the read-span map (user-specified).
SPAN_CLASS_ORDER = ["precore", "pgRNA", "preS1", "preS2_S", "X",
                    "spliced", "unclassified", "antisense", "pgRNA_RT"]
ANCHOR_X = 1936 + HBV_LEN          # canonical poly-A mapped into the copy-2 frame (~5118)


def _dedup_ids(sk):
    """UMI-dedup survivor read IDs (one per unique molecule) for a library, or None if the
    list is absent. Produced by phase2b_dedup_read_ids.sh; count == hbv_unique_molecules."""
    p = ANALYSIS / "samples" / sk / "hbv.dedup_read_ids.txt"
    if not p.exists():
        return None
    return {ln.strip() for ln in p.read_text().splitlines() if ln.strip()}


def _load_classified(sk):
    """Classified HBV reads for a library, restricted to UMI-deduplicated unique molecules
    (PCR repeats removed) so the transcript figures 4–6 match the unique-molecule counts in
    Fig 2. Falls back to all reads only if no dedup list exists."""
    p = ANALYSIS / "samples" / sk / "hbv_classified.parquet"
    if not p.exists():
        return None
    df = pd.read_parquet(p)
    ids = _dedup_ids(sk)
    if ids is not None:
        df = df[df.read_id.isin(ids)]
    return df


def dedup_class_abundance(meta):
    """Per-library transcript-class counts over UMI-deduplicated molecules — the dedup
    analogue of comparison/phase4_class_abundance.tsv (which counts all reads). Columns:
    sample_key, tx_class, n_reads (= unique molecules in that class)."""
    rows = []
    for sk in meta.index:
        df = _load_classified(sk)
        if df is None or df.empty:
            continue
        for cls, n in df.tx_class.value_counts().items():
            rows.append({"sample_key": sk, "tx_class": cls, "n_reads": int(n)})
    return pd.DataFrame(rows, columns=["sample_key", "tx_class", "n_reads"])


def _spans_from_panels(panels):
    """Render a multi-panel read-span figure. panels = [(ylabel_prefix, df), ...]. Reads are
    3′-anchored at the canonical poly-A (copy-2 frame ~5118), grouped by class and ordered
    longest→shortest within class. Panel height scales with read count (fixed per-read spacing
    so EVERY read is a separable line); tight y-limits mean the height is all reads, no margin."""
    panels = [(lbl, df) for lbl, df in panels if df is not None and not df.empty]
    if not panels:
        return None
    PER_READ, MIN_PANEL_H, MAX_PANEL_H = 0.06, 1.6, 200.0
    panel_h = [min(MAX_PANEL_H, max(MIN_PANEL_H, len(df) * PER_READ + 0.5)) for _, df in panels]
    fig, axes = plt.subplots(len(panels), 1, figsize=(11, sum(panel_h) + 0.3), sharex=True,
                             gridspec_kw={"height_ratios": panel_h, "hspace": 0.03})
    if len(panels) == 1:
        axes = [axes]
    ordkey = {c: i for i, c in enumerate(SPAN_CLASS_ORDER)}

    def parse_blocks(bs, ref_start, ref_end):
        if isinstance(bs, str) and bs:
            try:
                return [tuple(int(x) for x in tok.split("-")) for tok in bs.split(";")]
            except ValueError:
                pass
        return [(ref_start, ref_end)]

    xlo, xhi, present = ANCHOR_X, ANCHOR_X, set()
    for ax, (ylabel, df) in zip(axes, panels):
        d = df.copy()
        d["_ord"] = d.tx_class.map(lambda c: ordkey.get(c, 99))
        d["_span"] = d.ref_end - d.ref_start
        d = d.sort_values(["_ord", "_span"], ascending=[True, False]).reset_index(drop=True)
        has_blocks = "blocks" in d.columns
        for i, row in enumerate(d.itertuples()):
            col = CLASS_COLORS.get(row.tx_class, "#999")
            # 3′-anchor: shift so the poly-A end lands in the copy-2 frame (~5118).
            shift = (row.polya_end + HBV_LEN) - row.ref_end
            exons = parse_blocks(getattr(row, "blocks", None) if has_blocks else None,
                                 row.ref_start, row.ref_end)
            xs = [s + shift for s, _ in exons]
            xe = [e + shift for _, e in exons]
            if len(exons) > 1:                       # spliced: thin connector spans the intron(s)
                ax.plot([min(xs), max(xe)], [i, i], color=col, lw=0.5, alpha=0.55, zorder=1)
            for s, e in zip(xs, xe):                  # exon blocks (introns show as gaps)
                ax.plot([s, e], [i, i], color=col, lw=1.1, solid_capstyle="butt", zorder=2)
            xlo, xhi = min(xlo, min(xs)), max(xhi, max(xe))
        present.update(d.tx_class)
        ax.axvline(HBV_LEN, ls="--", color="k", lw=0.8)          # linearisation / copy boundary
        ax.axvline(ANCHOR_X, ls=":", color="darkred", lw=1.0)     # canonical poly-A anchor
        ax.set_ylim(len(d) - 0.5, -0.5)                           # tight (inverted): no blank margin above/below
        ax.set_yticks([])
        ax.set_ylabel(f"{ylabel}\nn={len(d)}", rotation=0, ha="right", va="center", fontsize=7)
    for ax in axes:
        ax.set_xlim(xlo - 80, xhi + 120)
    axes[-1].set_xlabel("Position on U95551.1_2x — 3′-anchored at poly-A (dotted ~5118); "
                        "dashed = linearisation 3182", fontsize=9)
    handles = [Line2D([0], [0], color=CLASS_COLORS[c], lw=3, label=clabel(c))
               for c in CLASS_ORDER if c in present]
    axes[0].legend(handles=handles, fontsize=7, ncol=min(4, len(handles)), loc="upper left")
    fig.tight_layout(h_pad=0.3)
    return fig


def fig_read_spans(meta):
    """Fig 6a — one panel per library (DISPLAY_ORDER)."""
    panels = [(f"{sk.split('_')[0]}\n{label(sk, meta).split(chr(10))[0]}", _load_classified(sk))
              for sk in DISPLAY_ORDER]
    return _spans_from_panels(panels)


def groups_150ng_bypcr(meta):
    """150 ng libraries only, combining polyA + NOpolyA within each PCR-cycle level.
    Returns [(group_label, [sample_key, ...]), ...] ordered by PCR cycles ascending."""
    m150 = meta[meta.input_ng == "150"]
    out = []
    for pcr in sorted(m150.pcr.unique(), key=lambda v: int(v)):
        sks = list(m150[m150.pcr == pcr].index)
        sks = sorted(sks, key=lambda k: meta.loc[k, "polya"] == "yes")   # NOpolyA before polyA
        out.append((f"150 ng {pcr}cyc\npolyA+NOpolyA", sks))
    return out


def fig_read_spans_combined_150(meta):
    """Fig 6b — 150 ng only, one panel per PCR-cycle level, pooling polyA + NOpolyA reads."""
    panels = []
    for glabel, sks in groups_150ng_bypcr(meta):
        dfs = [d for d in (_load_classified(sk) for sk in sks) if d is not None and not d.empty]
        panels.append((glabel.replace("\n", " "), pd.concat(dfs, ignore_index=True) if dfs else None))
    return _spans_from_panels(panels)


def fig_tss(meta):
    fig, ax = plt.subplots(figsize=(11, 4))
    keys = []
    for sk in DISPLAY_ORDER:
        df = _load_classified(sk)          # UMI-deduplicated unique molecules
        if df is None or df.empty:
            continue
        df = df[df.strand == "+"]
        if len(df) < 5:
            continue
        ax.hist(df.tss, bins=np.arange(0, HBV_LEN + 50, 50), histtype="step",
                lw=1.5, label=label(sk, meta).split("\n")[0])
        keys.append(sk)
    # shade class TSS windows (labels feed the band legend below)
    band_defs = [("precore/pgRNA", "#1f77b4", [(1730, 1880)]),
                 ("preS1", "#d62728", [(2700, 3100)]),
                 ("X", "#e6a817", [(1260, 1450)]),
                 ("preS2/S", "#2ca02c", [(3100, HBV_LEN), (0, 150)])]
    for _, c, spans in band_defs:
        for lo, hi in spans:
            ax.axvspan(lo, hi, color=c, alpha=0.08)
    ax.set_xlabel("5′ TSS on HBV genome (mod 3182); shaded = transcript-class TSS windows")
    ax.set_ylabel("unique molecules"); ax.set_title("TSS distribution (pychopper-oriented, UMI-deduplicated + molecules)")
    # two legends: per-library lines (upper right) + class-window colours (upper left)
    line_leg = ax.legend(fontsize=7, loc="upper right", title="library")
    ax.add_artist(line_leg)
    band_handles = [Patch(facecolor=c, alpha=0.25, label=name)
                    for name, c, _ in band_defs]
    ax.legend(handles=band_handles, fontsize=7, loc="upper left",
              title="class window", framealpha=0.9)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------

def tbl(df, cols=None):
    if cols:
        df = df[cols]
    return df.fillna("—").to_html(index=False, border=0, classes="data", escape=False)


def main():
    meta = meta_df()
    qc = load("phase1_5_read_qc.tsv")
    rec = load("phase4_recovery.tsv")
    # Transcript figures (4–6) are drawn on UMI-deduplicated unique molecules (PCR repeats
    # removed) so they match Fig 2's unique-molecule counts, not the raw-read class table.
    cls_df = dedup_class_abundance(meta)
    align = load("phase1_align_summary.tsv")
    today = datetime.date.today().isoformat()

    r = rec.set_index("sample_key")
    best_eff = rec.loc[rec.hbv_unique_per_M.idxmax()] if len(rec) else None
    best_abs = rec.loc[rec.hbv_unique.idxmax()] if len(rec) else None

    figs = {
        "qc": img(fig_qc(qc, meta), "Figure 1. Per-condition whole-library QC (all pass reads, not HBV-aligning only): yield (log-scale pass reads), read-length N50, mean Q, and % reads with a detected poly(A) tail."),
        "recovery": img(fig_recovery(rec, meta), "Figure 2. HBV recovery per condition. Left: raw HBV reads (grey) vs UMI-deduplicated unique molecules (blue, Poisson 95% CI). Right: unique HBV molecules per million composite-mapped reads (efficiency); green bars are poly(A)-selected."),
        "dup": img(fig_dup(rec, meta), "Figure 3. PCR-duplication rate from UMIs (well-formed 28-nt UMIs). Red = 20 PCR cycles, blue = 16. Higher input/cycle libraries duplicate more."),
        "cls": img(fig_class_composition(cls_df, meta), "Figure 4a. HBV transcript-class composition per condition, one bar per library (% of unique, UMI-deduplicated HBV molecules — PCR repeats removed). preS2/S dominates, as expected; the second-largest class is the canonical SP1 spliced variant (see text)."),
        "cls_combined": img(fig_class_composition_combined_150(cls_df, meta), "Figure 4b. Same composition (unique molecules) for the 150 ng libraries only, with polyA and NOpolyA pooled within each PCR-cycle level (16 cyc, 20 cyc) — the combined-fraction view."),
        "tss": img(fig_tss(meta), "Figure 5. 5′ TSS distribution over unique, UMI-deduplicated molecules (pychopper-oriented +, mod 3182). Shaded bands are the transcript-class TSS windows."),
    }
    span_caption = ("HBV read spans over unique, UMI-deduplicated molecules (PCR repeats removed), "
        "3′-anchored at the canonical poly-A (dotted red, copy-2 "
        "frame ~5118) so co-terminal transcripts line up. Each line is one unique molecule; colour = "
        "transcript class; reads are grouped by class (precore, pgRNA, preS1, preS2/S, X, spliced, "
        "unclassified) and ordered longest→shortest within class. Dashed line = 3182 linearisation "
        "(copy boundary): reads extending left of it wrap into the first genome copy; reads "
        "extending right of the poly-A anchor are 3′ readthrough.")
    span_figs = ""
    f = fig_read_spans(meta)
    if f is not None:
        span_figs = img(f, "Figure 6a. One panel per library. " + span_caption)
    span_figs_combined = ""
    fc = fig_read_spans_combined_150(meta)
    if fc is not None:
        span_figs_combined = img(fc, "Figure 6b. 150 ng libraries only, with polyA and NOpolyA "
                                 "reads pooled within each PCR-cycle level (16 cyc, 20 cyc). "
                                 + span_caption)

    css = """
    body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;max-width:1000px;margin:2rem auto;padding:0 1rem;color:#222;line-height:1.5}
    h1{border-bottom:3px solid #1f77b4;padding-bottom:.3rem} h2{border-bottom:1px solid #ccc;margin-top:2rem}
    figure{margin:1.2rem 0;text-align:center} img{max-width:100%;border:1px solid #eee}
    figcaption{font-size:.85rem;color:#555;text-align:left;margin-top:.3rem}
    table.data{border-collapse:collapse;font-size:.85rem;margin:1rem 0} table.data th,table.data td{border:1px solid #ddd;padding:3px 8px;text-align:right}
    table.data th{background:#f4f6f8} .key{background:#eef7ff;padding:.8rem 1rem;border-left:4px solid #1f77b4;margin:1rem 0}
    code{background:#f4f4f4;padding:1px 4px;border-radius:3px}
    """

    rec_cols = ["barcode", "input_ng", "polya", "pcr_cyc", "composite_mapped",
                "hbv_reads", "hbv_unique", "hbv_unique_ci", "hbv_unique_per_M", "dup_rate"]

    html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>EXP26000559 cDNA HBV Baseline</title><style>{css}</style></head><body>
<h1>HBV Detection Baseline with Standard cDNA (SQK-PCB114) — EXP26000559_cDNA001</h1>
<p><em>Generated {today} · pipeline {VERSION}</em></p>

<div class="key">
<b>Executive summary.</b> This run establishes the baseline HBV recovery of a conventional
Oxford Nanopore PCR-cDNA prep (no enrichment) across an 8-library matrix
(input × poly(A) selection × PCR cycles), on HBV-infected PHH total RNA (same source as
EXP26000465 untreated). Baseline HBV recovery is <b>low — on the order of 10–31 unique
molecules</b> even from the deepest libraries, establishing the "line in the sand" for a
future biotinylated-probe-capture workflow.
<ul>
<li><b>Poly(A) selection is the dominant factor</b>: best efficiency
{best_eff.hbv_unique_per_M:.1f} unique HBV/M ({best_eff.barcode}, poly(A)) vs 0.4–0.9/M for
matched non-selected libraries — a ~40–90× gap, far outside the confidence intervals.</li>
<li><b>16 cycles beats 20 for efficiency</b>: 20-cycle libraries reach 90% PCR duplication
(most reads are amplification, not new molecules); 16-cycle ~41%.</li>
<li><b>Best absolute recovery</b>: {best_abs.barcode} ({best_abs.hbv_unique} unique molecules,
CI {best_abs.hbv_unique_ci}) — but at 90% duplication, i.e. heavy over-sequencing.</li>
<li><b>Recommended standard prep</b>: 150 ng input, poly(A)-selected; 16 cycles is the more
molecule-efficient choice. 10 ng input failed (≤1 unique HBV molecule).</li>
<li><b>preS2/S dominates</b> the HBV transcriptome in every condition with data, consistent
with EXP26000465 and HBV biology.</li>
</ul>
</div>

<h2>Methods</h2>
<p><b>Samples & sequencing.</b> PromethION P2I, FLO-PRO114M, SQK-PCB114-24 (PCR-cDNA
barcoding, stranded, UMI in the strand-switching primer). 8 barcoded libraries, one biology
(HBV-infected PHH, untreated total RNA), varying input (150/10 ng) × poly(A) selection
(yes/no) × PCR cycles (16/20). On-instrument SUP basecalling (dorado 7.13.6, DNA e8.2 400K).</p>
<p><b>Reference & alignment (Phases 0–1).</b> hg38 + <b>2×HBV</b> (U95551.1 doubled to 6,364 bp so
genome-wrapping reads align as one continuous record). minimap2 2.30 <code>-ax splice
--secondary=no</code> (cDNA: no strand forcing; both orientations sequenced). Coordinates are
kept in the 2× space (0–6363) so the wrap stays visible; <code>mod 3182</code> is used only
to derive the genomic TSS for classification.</p>
<p><b>HBV reads & UMIs (Phases 1.5–2).</b> HBV-aligning primary reads (full sequence incl. soft-clipped
primer/UMI ends) were passed through pychopper (<code>-k PCB114 -U</code>) to orient and
extract the structured 28-nt UMI, re-aligned to the HBV-only 2× contig, and de-duplicated
with <code>umi_tools dedup --method=directional</code> on the well-formed 28-nt UMIs →
<b>unique HBV molecules</b>. HBV recovery is reported as raw reads, unique molecules, and
per-million-composite-mapped.</p>
<p><b>Classification (Phase 3).</b> On pychopper-oriented reads the 5′ TSS is reduced to
single-copy coordinates (<code>mod 3182</code>) and binned against fixed U95551.1 windows,
treated <em>circularly</em> because preS2/S starts just upstream of the linearization point and
wraps 0: precore/pgRNA TSS 1730–1880 (precore if ≤1815, upstream of the core promoter, else
pgRNA); preS1 2700–3100; preS2·S ≥3100 or ≤150; X 1260–1450. Each class also requires a
<em>plausible transcript length</em> (precore/pgRNA wrap the genome, span ≥2.6 kb + junction-
crossing; preS1 ≥1.6 kb; preS2/S ≥1.2 kb; X 0.3–1 kb) — a read whose 5′ lands in a window but is
far too short to be that transcript is a fragment and falls to <b>unclassified</b>. <b>pgRNA_RT</b>
is reserved for <b>tandem/concatemeric readthrough</b> (span &gt; 3,982 bp; a single terminally-
redundant pgRNA at ~3.3 kb is <em>not</em> readthrough). The major HBV mRNAs are unspliced;
reads carrying a <b>real intron</b> (aligned-block gap &gt; 200 bp — well above the 72–141 bp
minimap2 splice-mode artifacts, which are merged) are grouped as <b>spliced</b> regardless of
TSS and annotated by junction — essentially all match the canonical <b>SP1</b> junction
(donor ~2447 / acceptor ~489), the minor HBSP-encoding spliced variant, verified as genuine
(reads are ~1.2 kb shorter than pgRNA; the intron is absent from the molecule even without
splice-mode alignment). Remaining reads are <b>unclassified</b>.</p>
<p><b>Quantification (Phase 4).</b> Poisson 95% CIs on absolute counts, Wilson 95% CIs on proportions/
rates. Counts are single/double digits, so intervals are reported throughout; robust findings
(poly(A) effect) survive them, sparse per-class calls do not.</p>

<h2>Results</h2>
<h3>Library QC</h3>
{figs['qc']}
<p>Poly(A)-selected libraries yielded far less total data than non-selected at matched
input/cycles but had the longest, highest-quality reads (Q ~23, N50 ~2.3 kb).</p>

<h3>HBV recovery — the baseline</h3>
{tbl(rec.reset_index()[rec_cols] if 'sample_key' in rec.columns else rec[rec_cols])}
{figs['recovery']}
{figs['dup']}

<h3>Transcript classification</h3>
<p><em>Figures 4–6 are drawn on the <b>UMI-deduplicated unique molecules</b> (PCR repeats removed),
the same molecule set as the blue bars in Fig 2 — not the raw HBV reads. Per-molecule counts are
therefore small; class fractions are more meaningful than absolute heights.</em></p>
{figs['cls']}
{figs['cls_combined']}
{figs['tss']}
<p><b>pgRNA is full-length by definition; a spliced molecule is therefore not pgRNA.</b>
The major HBV mRNAs are unspliced, so pgRNA/precore/preS1/preS2·S/X count only full-length,
colinear transcripts (unspliced pgRNA here = 3 reads). The one spliced species is the
well-described minor <b>SP1</b> variant (donor ~2447 / acceptor ~489, encoding HBSP): 73 reads
carry a real ~1.2 kb intron, 72 at the exact canonical SP1 junction. The splice is genuine, not
a splice-mode artifact — the reads are ~1.2 kb shorter than pgRNA (the intron is absent from the
molecule even under non-splice alignment), and <b>all 72 have canonical GT…AG intron boundaries</b>,
the spliceosomal signature. Most SP1 reads have a pgRNA-region 5′ start (spliced pgRNA). That
full-length molecules are captured at all (the 30 precore reads are ~3.3 kb) shows the low
unspliced-pgRNA count is not a length dropout; the exact spliced:unspliced ratio, however, may be
modestly inflated by cDNA/PCR length bias and is not over-interpreted. Tightly-clustered
72–141 bp gaps in some preS2/S reads are minimap2 splice-mode artifacts, not introns.</p>
<h4>Read-span maps (2× reference — the wrap made visible)</h4>
{span_figs}
{span_figs_combined}

<h2>Conclusions</h2>
<p>Standard cDNA recovers HBV at a low absolute level (~10–31 unique molecules) even at depth.
Poly(A) selection is essential (~40–90× efficiency gain); 150 ng input is required (10 ng
failed); 16 PCR cycles is preferable to 20 (equivalent unique yield, far less duplication).
This is the baseline against which biotinylated-probe capture will be measured.</p>

<h2>Appendix — Reproducibility</h2>
<p>minimap2 2.30, samtools 1.23.1, pychopper 2.7.10, umi_tools 1.1.6, pandas/scipy/statsmodels.
Pipeline: <code>scripts/phase0_build_ref.sh</code> → <code>phase1_align.sh</code> →
<code>phase1_5_readqc.sh</code> / <code>phase1_5_umi.sh</code> →
<code>phase2_hbv_umi.sh</code> → <code>phase3_classify.py</code> →
<code>phase4_quantify.py</code> → <code>phase5_report.py</code>. Tables in
<code>analysis/comparison/</code>, per-read data in <code>analysis/samples/&lt;lib&gt;/</code>.</p>
</body></html>"""

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(html)
    print(f"Report written: {OUT} ({OUT.stat().st_size/1024:.0f} KB)")


if __name__ == "__main__":
    main()
