"""Plotly figure builders for the HTML report.

A single shared theme plus semantic colors: warm = NMD-sensitive / PTC /
intron-retained, green = intact / coding / escaped, grey = not-applicable /
novel / resolution-failed.
"""

from collections import Counter

import pandas as pd
import plotly.graph_objects as go

FONT = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif"

# Semantic palette.
_GREEN = "#2e8b57"
_RED = "#d1495b"
_ORANGE = "#e8743b"
_DARKRED = "#c0392b"
_PURPLE = "#7b6cd9"
_AMBER = "#d9b44a"
_GREY = "#9aa0a6"
_BLUE = "#3d7ea6"

COLOR = {
    # nmd_status
    "sensitive": _RED,
    "escaped": _GREEN,
    "not_applicable": _GREY,
    # resolved_orf_status
    "intact": _GREEN,
    "ptc_premature": _ORANGE,
    "ptc_intron_retained": _DARKRED,
    "cds_extension": _PURPLE,
    "no_stop_in_read": _AMBER,
    "resolution_failed": _GREY,
    # coding potential
    "coding": _GREEN,
    "noncoding": _GREY,
    # completeness
    "full_length": _GREEN,
    "alt_3prime_end": _BLUE,
    "truncated_5p": _AMBER,
    "truncated_3p": _AMBER,
    "truncated_both": _ORANGE,
    "internal_fragment": _ORANGE,
    "novel_no_match": _GREY,
}


def _theme(fig: go.Figure, title: str, height: int = 320) -> go.Figure:
    fig.update_layout(
        title=dict(text=title, font=dict(size=15, family=FONT, color="#222")),
        font=dict(family=FONT, size=12, color="#333"),
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=10, r=20, t=44, b=30),
        height=height,
        showlegend=False,
    )
    fig.update_xaxes(gridcolor="#eee", zeroline=False)
    fig.update_yaxes(gridcolor="#eee", zeroline=False)
    return fig


def _colors_for(categories: list[str]) -> list[str]:
    return [COLOR.get(c, _BLUE) for c in categories]


def category_bar(
    counts: dict[str, int] | Counter,
    title: str,
    order: list[str] | None = None,
) -> go.Figure:
    """Horizontal bar of category counts with count + percent labels."""
    if not counts:
        return _theme(go.Figure(), title, height=200)
    if order:
        items = [(k, counts.get(k, 0)) for k in order if counts.get(k, 0)]
    else:
        items = sorted(counts.items(), key=lambda kv: kv[1])
    cats = [str(k) for k, _ in items]
    vals = [int(v) for _, v in items]
    total = sum(vals) or 1
    labels = [f"{v:,} ({v / total * 100:.1f}%)" for v in vals]
    fig = go.Figure(
        go.Bar(
            x=vals,
            y=cats,
            orientation="h",
            marker_color=_colors_for(cats),
            text=labels,
            textposition="outside",
            cliponaxis=False,
        )
    )
    fig.update_xaxes(visible=False)
    return _theme(fig, title, height=max(200, 50 + 34 * len(cats)))


def funnel(stages: list[tuple[str, int]]) -> go.Figure:
    """Consequence cascade: total -> coding parent -> ORF resolved -> intact ..."""
    labels = [s for s, _ in stages]
    values = [v for _, v in stages]
    fig = go.Figure(
        go.Funnel(
            y=labels,
            x=values,
            textinfo="value+percent initial",
            marker=dict(color=[_BLUE, _BLUE, _GREEN, _GREEN, _RED][: len(stages)]),
        )
    )
    return _theme(fig, "Functional-consequence cascade", height=300)


def histogram(
    values: pd.Series,
    title: str,
    threshold: float | None = None,
    color: str = _BLUE,
    nbins: int = 40,
) -> go.Figure:
    """Histogram of a continuous quantity; optional vertical threshold line."""
    vals = pd.to_numeric(values, errors="coerce").dropna()
    fig = go.Figure(go.Histogram(x=vals, nbinsx=nbins, marker_color=color))
    if threshold is not None and len(vals):
        fig.add_vline(x=threshold, line_dash="dash", line_color="#444")
    return _theme(fig, title, height=300)
