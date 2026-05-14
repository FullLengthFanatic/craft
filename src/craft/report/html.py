"""Top-level HTML report assembly.

Generates a self-contained interactive HTML summarising the per-isoform
annotations produced by :func:`craft.pipeline.run_annotate`. The report has
three sections: summary cards, distribution bar charts (plotly), and a
per-isoform table.
"""

from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from craft import __version__
from craft.report.plots import bar_chart

_SUMMARY_FIELDS = (
    ("completeness", "Completeness"),
    ("orf_outcome", "ORF outcome"),
    ("nmd_status", "NMD status"),
    ("orf_confidence", "ORF confidence"),
)

_HIDDEN_TABLE_COLUMNS = ("propagated_cds_intervals", "denovo_cds_intervals")

_DEFAULT_TABLE_ROW_CAP = 1000

_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>CRAFT annotation report</title>
<style>
  body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    max-width: 1400px;
    margin: 2em auto;
    padding: 0 1em;
    color: #222;
  }}
  h1, h2 {{ color: #111; }}
  h1 {{ margin-bottom: 0.25em; }}
  .meta {{ color: #666; font-size: 0.9em; margin-top: 0; }}
  .summary-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
    gap: 1em;
    margin: 1em 0 2em 0;
  }}
  .summary-card {{ background: #f7f7f7; padding: 0.8em 1em; border-radius: 4px; }}
  .summary-card h3 {{ margin: 0 0 0.4em 0; font-size: 0.95em; }}
  .summary-list {{ margin: 0; padding-left: 1.2em; font-size: 0.9em; }}
  .figures-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(420px, 1fr));
    gap: 1em;
    margin: 1em 0 2em 0;
  }}
  .table-container {{
    max-height: 600px;
    overflow: auto;
    border: 1px solid #ddd;
    border-radius: 4px;
    margin-top: 1em;
  }}
  table.craft-table {{ border-collapse: collapse; font-size: 0.83em; width: 100%; }}
  table.craft-table th, table.craft-table td {{
    padding: 4px 8px;
    border-bottom: 1px solid #eee;
    text-align: left;
    white-space: nowrap;
  }}
  table.craft-table th {{ background: #f0f0f0; position: sticky; top: 0; z-index: 1; }}
  .footer {{ color: #999; font-size: 0.8em; margin-top: 3em; }}
</style>
</head>
<body>
<h1>CRAFT annotation report</h1>
<p class="meta">Generated {timestamp} UTC | craft v{version} | {n_isoforms} isoforms</p>

<h2>Summary</h2>
<div class="summary-grid">
{summary_cards}
</div>

<h2>Distributions</h2>
<div class="figures-grid">
{figure_blocks}
</div>

<h2>Per-isoform annotations</h2>
{table_note}
<div class="table-container">
{table_html}
</div>

<p class="footer">CRAFT (Coding Region Annotation From Templates).</p>
</body>
</html>
"""


def _value_counts(df: pd.DataFrame, col: str) -> dict[str, int]:
    if col not in df.columns or df.empty:
        return {}
    series = df[col].dropna()
    return {str(k): int(v) for k, v in series.value_counts().items()}


def _summary_card(title: str, counts: dict[str, int], total: int) -> str:
    if not counts:
        return f'<div class="summary-card"><h3>{title}</h3><p>No data.</p></div>'
    items = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)
    lines = []
    for label, count in items:
        pct = (count / total * 100.0) if total else 0.0
        lines.append(f"<li>{label}: {count} ({pct:.1f}%)</li>")
    return (
        f'<div class="summary-card"><h3>{title}</h3>'
        f'<ul class="summary-list">{"".join(lines)}</ul></div>'
    )


def _figure_block(fig, *, include_js: bool) -> str:
    return fig.to_html(
        include_plotlyjs="inline" if include_js else False,
        full_html=False,
    )


def _table_html(df: pd.DataFrame, row_cap: int) -> tuple[str, int]:
    drop = [c for c in _HIDDEN_TABLE_COLUMNS if c in df.columns]
    display = df.drop(columns=drop)
    truncated_n = max(0, len(display) - row_cap)
    if truncated_n:
        display = display.head(row_cap)
    return (
        display.to_html(
            index=False,
            classes="craft-table",
            border=0,
            na_rep="",
        ),
        truncated_n,
    )


def render(
    per_isoform: pd.DataFrame,
    output: Path,
    row_cap: int = _DEFAULT_TABLE_ROW_CAP,
) -> None:
    """Render the interactive HTML report.

    Args:
        per_isoform: Per-isoform annotation table (output of
            :func:`craft.pipeline.run_annotate`).
        output: Destination ``.html`` path. Parent directories are created.
        row_cap: Maximum rows shown in the per-isoform table. Excess rows are
            truncated with a note; the full data is in the companion TSV/JSON.
    """
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)

    n_isoforms = len(per_isoform)
    summary = {col: _value_counts(per_isoform, col) for col, _ in _SUMMARY_FIELDS}

    summary_cards = "\n".join(
        _summary_card(title, summary[col], n_isoforms) for col, title in _SUMMARY_FIELDS
    )

    figs = [(title, bar_chart(summary[col], title)) for col, title in _SUMMARY_FIELDS]
    figure_blocks = "\n".join(
        f"<div>{_figure_block(fig, include_js=(i == 0))}</div>"
        for i, (_, fig) in enumerate(figs)
    )

    table_html, truncated = _table_html(per_isoform, row_cap)
    table_note = (
        f"<p class='meta'>Showing first {row_cap} of {n_isoforms} isoforms. "
        f"See the TSV/JSON for the full table.</p>"
        if truncated
        else ""
    )

    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M")
    html = _TEMPLATE.format(
        timestamp=timestamp,
        version=__version__,
        n_isoforms=n_isoforms,
        summary_cards=summary_cards,
        figure_blocks=figure_blocks,
        table_note=table_note,
        table_html=table_html,
    )
    output.write_text(html)
