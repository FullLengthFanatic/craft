#!/usr/bin/env python
"""Bench 3 step 3: DE analysis + NMD-target enrichment + figure assembly.

Reads the 6 salmon quants under ``benchmarks/cache/bench3/quants/`` and the
CRAFT NMD label cache under ``benchmarks/cache/bench3/nmd_labels.tsv.gz``,
runs pydeseq2 (UPF1 KD vs scr control), tests whether CRAFT's NMD-sensitive
transcripts are enriched among UPF1-KD-upregulated transcripts, and emits a
1x3 plotly panel.

Run from the repo root:

    PYTHONPATH=benchmarks .venv/bin/python benchmarks/run_bench3_analysis.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import fisher_exact

REPO_ROOT = Path(__file__).resolve().parent.parent
CACHE = REPO_ROOT / "benchmarks/cache/bench3"
QUANT_DIR = CACHE / "quants"
NMD_TSV = CACHE / "nmd_labels.tsv.gz"
DE_OUT = CACHE / "de_results.tsv.gz"
ENRICH_TSV = REPO_ROOT / "benchmarks/figures/bench3_enrichment.tsv"
FIG_PNG = REPO_ROOT / "benchmarks/figures/bench3_enrichment_panel.png"
FIG_JSON = REPO_ROOT / "benchmarks/figures/bench3_enrichment_panel.json"

CONTROL = ("SRR4081222", "SRR4081223", "SRR4081224")
UPF1_KD = ("SRR4081225", "SRR4081226", "SRR4081227")
PADJ_THRESHOLD = 0.05
L2FC_UP_THRESHOLD = 1.0  # log2 fold-change cutoff for "upregulated"

# CRAFT NMD-sensitive call only makes sense for transcripts whose ORF outcome
# is propagated_intact or disrupted with a clean stop. Exclude STOP_AT_ALT_POLYA
# and STOP_NOT_OBSERVED rows; their nmd_status is forced to not_applicable.
KEEP_OUTCOMES = ("propagated_intact", "disrupted")


def progress(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def load_counts() -> pd.DataFrame:
    """Return samples-x-transcripts integer count matrix (NumReads rounded)."""
    frames: dict[str, pd.Series] = {}
    for srr in (*CONTROL, *UPF1_KD):
        quant = QUANT_DIR / srr / "quant.sf"
        sub = pd.read_csv(quant, sep="\t", usecols=["Name", "NumReads"])
        frames[srr] = sub.set_index("Name")["NumReads"]
    counts = pd.DataFrame(frames).fillna(0).round().astype("int64")
    # Salmon transcript IDs are full GENCODE names with version, e.g.
    # ENST00000456328.2|ENSG00000223972.5|...|; pydeseq2 just needs a unique
    # row label so the pipe-delimited string is fine.
    progress(f"[de] counts matrix: {counts.shape[0]} transcripts x {counts.shape[1]} samples")
    return counts


def stem_transcript_id(name: str) -> str:
    """Strip GENCODE's pipe-delimited extras, leaving just the ENST id."""
    return name.split("|", 1)[0]


def run_deseq(counts: pd.DataFrame) -> pd.DataFrame:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    metadata = pd.DataFrame(
        {"condition": ["control"] * len(CONTROL) + ["upf1_kd"] * len(UPF1_KD)},
        index=list(CONTROL) + list(UPF1_KD),
    )

    # pydeseq2 wants samples on rows.
    counts_T = counts.T
    keep_tx = (counts_T.sum(axis=0) >= 10) & ((counts_T > 0).sum(axis=0) >= 2)
    progress(f"[de] keeping {int(keep_tx.sum())} transcripts after low-count filter")
    counts_T = counts_T.loc[:, keep_tx]

    progress("[de] running DeseqDataSet.deseq2()...")
    t0 = time.time()
    inf = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(
        counts=counts_T,
        metadata=metadata,
        design="~condition",
        inference=inf,
        quiet=True,
    )
    dds.deseq2()
    progress(f"[de] deseq2 fit in {(time.time() - t0) / 60:.1f}m")

    ds = DeseqStats(
        dds,
        contrast=["condition", "upf1_kd", "control"],
        inference=inf,
        quiet=True,
    )
    ds.summary()
    res = ds.results_df.copy()
    res.index.name = "transcript_id_full"
    res.reset_index(inplace=True)
    res["transcript_id"] = res["transcript_id_full"].map(stem_transcript_id)
    progress(f"[de] significant @ padj<{PADJ_THRESHOLD}: {(res['padj'] < PADJ_THRESHOLD).sum()}")
    return res


def enrichment(merged: pd.DataFrame) -> dict:
    """Fisher's exact: NMD-sensitive vs upregulated in UPF1-KD."""
    eligible = merged[
        merged["orf_outcome"].isin(KEEP_OUTCOMES)
        & merged["nmd_status"].isin(["sensitive", "escaped"])
        & merged["log2FoldChange"].notna()
        & merged["padj"].notna()
    ].copy()
    eligible["is_sensitive"] = eligible["nmd_status"] == "sensitive"
    eligible["is_up"] = (
        (eligible["log2FoldChange"] >= L2FC_UP_THRESHOLD)
        & (eligible["padj"] < PADJ_THRESHOLD)
    )

    a = int(((eligible["is_sensitive"]) & (eligible["is_up"])).sum())
    b = int(((eligible["is_sensitive"]) & (~eligible["is_up"])).sum())
    c = int(((~eligible["is_sensitive"]) & (eligible["is_up"])).sum())
    d = int(((~eligible["is_sensitive"]) & (~eligible["is_up"])).sum())
    table = [[a, b], [c, d]]
    odds, pval = fisher_exact(table, alternative="greater")
    progress(
        f"[enrich] eligible tx = {len(eligible):,}; "
        f"sensitive+up={a}, sensitive+not={b}, other+up={c}, other+not={d}; "
        f"OR={odds:.3f}, one-sided p={pval:.3e}"
    )
    return {
        "table": table,
        "odds_ratio": float(odds),
        "p_value": float(pval),
        "n_sensitive": a + b,
        "n_not_sensitive": c + d,
        "n_up": a + c,
        "n_not_up": b + d,
        "n_total": int(len(eligible)),
        "l2fc_threshold": L2FC_UP_THRESHOLD,
        "padj_threshold": PADJ_THRESHOLD,
    }, eligible


def plot_panel(eligible: pd.DataFrame, result: dict) -> go.Figure:
    fig = make_subplots(
        rows=1,
        cols=3,
        subplot_titles=(
            f"Contingency (OR={result['odds_ratio']:.2f}, p={result['p_value']:.1e})",
            "Volcano (UPF1-KD vs control)",
            "log2FC CDF by CRAFT NMD label",
        ),
        column_widths=[0.30, 0.40, 0.30],
        horizontal_spacing=0.10,
    )

    a, b = result["table"][0]
    c, d = result["table"][1]
    z = [[a, b], [c, d]]
    fig.add_trace(
        go.Heatmap(
            z=z,
            x=["upregulated", "not upregulated"],
            y=["NMD-sensitive", "NMD-escaped"],
            text=[[f"{a:,}", f"{b:,}"], [f"{c:,}", f"{d:,}"]],
            texttemplate="%{text}",
            colorscale=[[0, "#eef2f6"], [1, "#5b7a9d"]],
            showscale=False,
            zmin=0,
        ),
        row=1,
        col=1,
    )

    for is_sens, color, name in (
        (False, "#cdd6df", "NMD-escaped"),
        (True, "#5b7a9d", "NMD-sensitive"),
    ):
        sub = eligible[eligible["is_sensitive"] == is_sens]
        fig.add_trace(
            go.Scattergl(
                x=sub["log2FoldChange"],
                y=-np.log10(sub["padj"].clip(lower=1e-300)),
                mode="markers",
                marker=dict(color=color, size=4, opacity=0.5),
                name=name,
                legendgroup=name,
            ),
            row=1,
            col=2,
        )
    fig.add_vline(x=L2FC_UP_THRESHOLD, line_dash="dot", line_color="#999", row=1, col=2)
    fig.add_hline(
        y=-np.log10(PADJ_THRESHOLD), line_dash="dot", line_color="#999", row=1, col=2
    )

    for is_sens, color, name in (
        (False, "#cdd6df", "NMD-escaped"),
        (True, "#5b7a9d", "NMD-sensitive"),
    ):
        sub = eligible[eligible["is_sensitive"] == is_sens]["log2FoldChange"].dropna()
        if sub.empty:
            continue
        x = np.sort(sub.values)
        y = np.linspace(0, 1, len(x), endpoint=False)
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                mode="lines",
                line=dict(color=color, width=2),
                name=name,
                showlegend=False,
                legendgroup=name,
            ),
            row=1,
            col=3,
        )
    fig.add_vline(x=0, line_dash="dot", line_color="#999", row=1, col=3)

    fig.update_xaxes(title_text="log2 FC (UPF1-KD / control)", row=1, col=2)
    fig.update_yaxes(title_text="-log10 padj", row=1, col=2)
    fig.update_xaxes(title_text="log2 FC", row=1, col=3, range=[-4, 6])
    fig.update_yaxes(title_text="cumulative fraction", row=1, col=3, range=[0, 1])
    fig.update_layout(
        height=420,
        width=1200,
        template="plotly_white",
        title_text=(
            f"Bench 3: CRAFT NMD-sensitive isoforms in UPF1-KD bulk RNA-seq "
            f"(GSE86148, HeLa, n=3v3; eligible tx={result['n_total']:,})"
        ),
        legend=dict(orientation="h", y=-0.18, x=0.5, xanchor="center"),
    )
    return fig


def main() -> int:
    counts = load_counts()
    if DE_OUT.exists():
        progress(f"[de] cached results at {DE_OUT}")
        de = pd.read_csv(DE_OUT, sep="\t")
    else:
        de = run_deseq(counts)
        DE_OUT.parent.mkdir(parents=True, exist_ok=True)
        de.to_csv(DE_OUT, sep="\t", index=False, compression="gzip")
        progress(f"[de] wrote {DE_OUT}")

    nmd = pd.read_csv(NMD_TSV, sep="\t")
    nmd["transcript_id_stem"] = nmd["transcript_id"].astype(str)
    de_join = de.copy()
    de_join["transcript_id_stem"] = de_join["transcript_id"]
    merged = de_join.merge(nmd, left_on="transcript_id_stem", right_on="transcript_id_stem", how="inner", suffixes=("_de", "_nmd"))
    progress(f"[merge] joined {len(merged):,} transcripts (DE x NMD universe)")

    result, eligible = enrichment(merged)

    ENRICH_TSV.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            ("contingency_a_sens_up", result["table"][0][0]),
            ("contingency_b_sens_notup", result["table"][0][1]),
            ("contingency_c_notsens_up", result["table"][1][0]),
            ("contingency_d_notsens_notup", result["table"][1][1]),
            ("odds_ratio", result["odds_ratio"]),
            ("p_value_one_sided", result["p_value"]),
            ("n_eligible", result["n_total"]),
            ("n_sensitive", result["n_sensitive"]),
            ("n_upregulated", result["n_up"]),
            ("l2fc_threshold", result["l2fc_threshold"]),
            ("padj_threshold", result["padj_threshold"]),
        ],
        columns=["metric", "value"],
    ).to_csv(ENRICH_TSV, sep="\t", index=False)
    progress(f"[enrich] wrote {ENRICH_TSV}")

    fig = plot_panel(eligible, result)
    fig.write_image(str(FIG_PNG), scale=2)
    FIG_JSON.write_text(fig.to_json())
    progress(f"[fig] PNG: {FIG_PNG} ({FIG_PNG.stat().st_size / 1024:.1f} KB)")
    progress(f"[fig] JSON: {FIG_JSON}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
