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

_NOTABLE_TOP_N = 10

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
  table.craft-table {{
    border-collapse: collapse;
    font-size: 0.85em;
    width: auto;
    margin: 0.5em 0 1.5em 0;
  }}
  table.craft-table th, table.craft-table td {{
    padding: 4px 10px;
    border-bottom: 1px solid #eee;
    text-align: left;
    white-space: nowrap;
  }}
  table.craft-table th {{ background: #f0f0f0; }}
  table.craft-table-compact td {{ font-variant-numeric: tabular-nums; }}
  h3 .count {{ color: #888; font-weight: normal; font-size: 0.85em; }}
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

<h2>Notable findings</h2>
{notable_findings}

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


def _small_table_html(df: pd.DataFrame, columns: list[str]) -> str:
    cols_present = [c for c in columns if c in df.columns]
    return df[cols_present].to_html(
        index=False,
        classes="craft-table craft-table-compact",
        border=0,
        na_rep="-",
    )


def _top_nmd_sensitive(df: pd.DataFrame) -> tuple[str, int]:
    if "nmd_status" not in df.columns or "nmd_confidence" not in df.columns:
        return "", 0
    sub = df[(df["nmd_status"] == "sensitive") & (df["nmd_confidence"] == "high")]
    n_total = len(sub)
    if n_total == 0:
        return "", 0
    sort_col = "orf_confidence_score" if "orf_confidence_score" in sub.columns else None
    top = sub.nlargest(_NOTABLE_TOP_N, sort_col) if sort_col else sub.head(_NOTABLE_TOP_N)
    cols = [
        "transcript_id",
        "parent_gene_name",
        "parent_tx_id",
        "completeness",
        "nmd_rule",
        "stop_to_last_junction_nt",
        "orf_confidence_score",
    ]
    return _small_table_html(top, cols), n_total


def _top_disrupted(df: pd.DataFrame) -> tuple[str, int]:
    if "orf_outcome" not in df.columns or "orf_confidence" not in df.columns:
        return "", 0
    sub = df[(df["orf_outcome"] == "disrupted") & (df["orf_confidence"] == "high")].copy()
    n_total = len(sub)
    if n_total == 0:
        return "", 0
    if {"parent_cds_bp", "propagated_cds_bp"}.issubset(sub.columns):
        sub["cds_bp_lost"] = sub["parent_cds_bp"] - sub["propagated_cds_bp"]
        top = sub.nlargest(_NOTABLE_TOP_N, "cds_bp_lost")
    else:
        top = sub.head(_NOTABLE_TOP_N)
    cols = [
        "transcript_id",
        "parent_gene_name",
        "parent_tx_id",
        "completeness",
        "propagated_cds_bp",
        "parent_cds_bp",
        "cds_bp_lost",
        "pfam_lost",
    ]
    return _small_table_html(top, cols), n_total


def _top_isoform_diversity(df: pd.DataFrame) -> tuple[str, int]:
    """Distinct functional variants per parent gene.

    Counts unique ``(parent_tx_id, orf_outcome)`` pairs per
    ``parent_gene_id``, restricted to isoforms with high or medium ORF
    confidence. This collapses the PacBio-collapse over-fragmentation noise
    (50 PB.X.Y entries that all map to the same parent transcript with the
    same outcome collapse to a single functional variant).
    """
    required = {"parent_gene_id", "parent_tx_id", "orf_outcome", "orf_confidence"}
    if not required.issubset(df.columns):
        return "", 0
    sub = df[df["orf_confidence"].isin(["high", "medium"])].copy()
    sub["parent_gene_id"] = sub["parent_gene_id"].astype(str)
    sub = sub[sub["parent_gene_id"] != ""]
    if sub.empty:
        return "", 0

    distinct_pairs = sub[
        ["parent_gene_id", "parent_tx_id", "orf_outcome"]
    ].drop_duplicates()
    diversity = (
        distinct_pairs.groupby("parent_gene_id")
        .size()
        .rename("n_functional_variants")
        .reset_index()
    )

    if "parent_gene_name" in sub.columns:
        names = (
            sub.drop_duplicates("parent_gene_id")
            .set_index("parent_gene_id")["parent_gene_name"]
            .to_dict()
        )
        diversity["parent_gene_name"] = diversity["parent_gene_id"].map(names).fillna("")

    diversity = diversity.sort_values("n_functional_variants", ascending=False).head(
        _NOTABLE_TOP_N
    )
    cols = ["parent_gene_name", "parent_gene_id", "n_functional_variants"]
    total = (
        int(diversity["n_functional_variants"].sum()) if not diversity.empty else 0
    )
    return _small_table_html(diversity, cols), total


def _notable_findings_html(df: pd.DataFrame) -> str:
    sections: list[str] = []

    nmd_html, nmd_total = _top_nmd_sensitive(df)
    if nmd_html:
        sections.append(
            f'<h3>Top NMD-sensitive isoforms <span class="count">({nmd_total} '
            f'total high-confidence)</span></h3>{nmd_html}'
        )

    disrupted_html, disrupted_total = _top_disrupted(df)
    if disrupted_html:
        sections.append(
            f'<h3>Top ORF-disrupted isoforms <span class="count">({disrupted_total} '
            f'total high-confidence)</span></h3>{disrupted_html}'
        )

    diversity_html, variant_total = _top_isoform_diversity(df)
    if diversity_html:
        sections.append(
            f'<h3>Genes with most functional isoform diversity '
            f'<span class="count">(distinct (parent_tx, orf_outcome) pairs '
            f'among high/medium-confidence isoforms; '
            f'{variant_total} variants in top {_NOTABLE_TOP_N})</span></h3>'
            f'{diversity_html}'
        )

    if not sections:
        return (
            '<p class="meta">No notable findings to surface. See '
            '<code>per_isoform.tsv</code> for the full per-isoform table.</p>'
        )
    return "\n".join(sections)


def render(per_isoform: pd.DataFrame, output: Path) -> None:
    """Render the interactive HTML report.

    Args:
        per_isoform: Per-isoform annotation table (output of
            :func:`craft.pipeline.run_annotate`).
        output: Destination ``.html`` path. Parent directories are created.
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

    notable_findings = _notable_findings_html(per_isoform)

    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M")
    html = _TEMPLATE.format(
        timestamp=timestamp,
        version=__version__,
        n_isoforms=n_isoforms,
        summary_cards=summary_cards,
        figure_blocks=figure_blocks,
        notable_findings=notable_findings,
    )
    output.write_text(html)
