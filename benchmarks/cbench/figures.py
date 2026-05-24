"""Plotly figure helpers for Bench 1.

Colours and layout track scNoiseMeter / tecap / CRAFT's existing reports.
"""

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Two-tool palette (CRAFT propagation in the same muted slate blue used by the
# v1.2+ HTML report; orfipy in a warm complement).
COLOR_CRAFT = "#5b7a9d"
COLOR_ORFIPY = "#c97a5b"
COMPARATOR_COLORS = {"craft": COLOR_CRAFT, "orfipy": COLOR_ORFIPY}

ORIENTATIONS = ("5prime", "3prime", "both")


def recovery_panel(summary: pd.DataFrame) -> go.Figure:
    """2x2 panel: recovery rate by truncation rate, faceted by orientation, plus
    confidence calibration in the bottom-right if a column ``orf_confidence`` is
    present.

    Expects ``summary`` to have columns ``comparator``, ``rate``, ``orientation``,
    ``recovery_rate``, ``start_exact_rate``, ``mean_abs_length_error``. If a
    ``calibration`` DataFrame is added later, swap into the bottom-right panel.
    """
    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=(
            "Recovery rate by truncation %",
            "Start-codon exact-match rate",
            "Mean |ORF length error| (nt)",
            "Stop-codon exact-match rate",
        ),
        horizontal_spacing=0.12,
        vertical_spacing=0.18,
    )
    for comparator, color in COMPARATOR_COLORS.items():
        sub = summary[summary["comparator"] == comparator]
        if sub.empty:
            continue
        for orientation in ORIENTATIONS:
            cell = (
                sub[sub["orientation"] == orientation]
                .groupby("rate", as_index=False)
                .agg(
                    recovery_rate=("recovery_rate", "mean"),
                    start_exact_rate=("start_exact_rate", "mean"),
                    stop_exact_rate=("stop_exact_rate", "mean"),
                    mean_abs_length_error=("mean_abs_length_error", "mean"),
                )
                .sort_values("rate")
            )
            dash = {"5prime": "solid", "3prime": "dash", "both": "dot"}[orientation]
            fig.add_trace(
                go.Scatter(
                    x=cell["rate"] * 100,
                    y=cell["recovery_rate"],
                    mode="lines+markers",
                    line=dict(color=color, dash=dash),
                    name=f"{comparator} | {orientation}",
                    legendgroup=comparator,
                ),
                row=1,
                col=1,
            )
            fig.add_trace(
                go.Scatter(
                    x=cell["rate"] * 100,
                    y=cell["start_exact_rate"],
                    mode="lines+markers",
                    line=dict(color=color, dash=dash),
                    name=f"{comparator} | {orientation}",
                    showlegend=False,
                    legendgroup=comparator,
                ),
                row=1,
                col=2,
            )
            fig.add_trace(
                go.Scatter(
                    x=cell["rate"] * 100,
                    y=cell["mean_abs_length_error"],
                    mode="lines+markers",
                    line=dict(color=color, dash=dash),
                    name=f"{comparator} | {orientation}",
                    showlegend=False,
                    legendgroup=comparator,
                ),
                row=2,
                col=1,
            )
            fig.add_trace(
                go.Scatter(
                    x=cell["rate"] * 100,
                    y=cell["stop_exact_rate"],
                    mode="lines+markers",
                    line=dict(color=color, dash=dash),
                    name=f"{comparator} | {orientation}",
                    showlegend=False,
                    legendgroup=comparator,
                ),
                row=2,
                col=2,
            )

    fig.update_xaxes(title_text="Truncation rate (%)")
    fig.update_yaxes(title_text="Rate", row=1, col=1, range=[0, 1])
    fig.update_yaxes(title_text="Rate", row=1, col=2, range=[0, 1])
    fig.update_yaxes(title_text="nt", row=2, col=1)
    fig.update_yaxes(title_text="Rate", row=2, col=2, range=[0, 1])
    fig.update_layout(
        height=640,
        width=900,
        template="plotly_white",
        title_text="Bench 1: simulated truncation - CRAFT propagation vs orfipy",
    )
    return fig


def calibration_bar(calibration: pd.DataFrame) -> go.Figure:
    """Standalone bar chart for the CRAFT confidence calibration table.

    Expects columns ``orf_confidence`` and one or more rate columns (e.g.,
    ``recovery_rate``, ``start_exact_rate``).
    """
    order = ["high", "medium", "low", "none", "missing"]
    calibration = calibration.copy()
    calibration["__order"] = calibration["orf_confidence"].apply(
        lambda v: order.index(v) if v in order else len(order)
    )
    calibration = calibration.sort_values("__order").drop(columns="__order")
    fig = go.Figure()
    for metric, marker in (
        ("recovery_rate", COLOR_CRAFT),
        ("start_exact_rate", "#7c97b3"),
        ("stop_exact_rate", "#9eb3c8"),
    ):
        if metric not in calibration.columns:
            continue
        fig.add_trace(
            go.Bar(
                x=calibration["orf_confidence"],
                y=calibration[metric],
                name=metric.replace("_", " "),
                marker_color=marker,
            )
        )
    fig.update_layout(
        barmode="group",
        height=380,
        width=620,
        template="plotly_white",
        title_text="CRAFT confidence calibration on truncated transcripts",
        xaxis_title="CRAFT orf_confidence",
        yaxis_title="Rate",
        yaxis_range=[0, 1],
    )
    return fig
