"""Plotly figure builders for the HTML report."""

from collections import Counter

import pandas as pd
import plotly.graph_objects as go


def bar_chart(counts: dict[str, int] | Counter, title: str) -> go.Figure:
    """Horizontal bar chart of category counts, sorted descending."""
    if not counts:
        fig = go.Figure()
        fig.update_layout(title=title, height=240)
        return fig
    items = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)
    categories = [str(k) for k, _ in items]
    values = [int(v) for _, v in items]
    fig = go.Figure(data=[go.Bar(x=categories, y=values)])
    fig.update_layout(
        title=title,
        margin=dict(l=40, r=20, t=40, b=80),
        height=320,
        xaxis_tickangle=-30,
    )
    return fig


def isoform_track(per_isoform: pd.DataFrame, gene_id: str) -> go.Figure:
    """Per-gene isoform track view (exon + CDS structure). Deferred to v1.1."""
    raise NotImplementedError("Per-gene track view is not implemented in v1.")
