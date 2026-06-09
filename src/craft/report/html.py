"""Top-level HTML report assembly.

A self-contained interactive report: a KPI strip, a functional-consequence
cascade, semantic distribution bars (completeness, resolved ORF status, NMD),
coding-potential and UTR-delta histograms, and notable-findings tables.
"""

from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from craft import __version__
from craft.report.plots import category_bar, funnel, histogram

_NOTABLE_TOP_N = 10

_COMPLETENESS_ORDER = [
    "full_length",
    "alt_3prime_end",
    "truncated_5p",
    "truncated_3p",
    "truncated_both",
    "internal_fragment",
    "novel_no_match",
]
_RESOLVED_ORDER = [
    "intact",
    "ptc_premature",
    "ptc_intron_retained",
    "cds_extension",
    "no_stop_in_read",
    "resolution_failed",
]
_NMD_ORDER = ["sensitive", "escaped", "not_applicable"]
_RESOLVED_WITH_STOP = {"intact", "ptc_premature", "ptc_intron_retained", "cds_extension"}

_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>CRAFT annotation report</title>
<style>
  body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    max-width: 1320px; margin: 2em auto; padding: 0 1.2em; color: #222;
  }}
  h1 {{ margin-bottom: 0.1em; font-weight: 650; }}
  h2 {{ color: #111; margin-top: 1.6em; border-bottom: 1px solid #eee; padding-bottom: 0.2em; }}
  .meta {{ color: #888; font-size: 0.85em; margin-top: 0; }}
  .kpis {{ display: flex; flex-wrap: wrap; gap: 0.8em; margin: 1.2em 0; }}
  .kpi {{
    flex: 1 1 150px; background: #f7f8fa; border: 1px solid #eef0f3;
    border-radius: 8px; padding: 0.8em 1em;
  }}
  .kpi .val {{ font-size: 1.7em; font-weight: 650; color: #1a1a1a; }}
  .kpi .lab {{ font-size: 0.82em; color: #666; margin-top: 0.2em; }}
  .grid {{
    display: grid; grid-template-columns: repeat(auto-fit, minmax(440px, 1fr)); gap: 0.6em;
  }}
  table.craft-table {{ border-collapse: collapse; font-size: 0.85em; margin: 0.4em 0 1.4em 0; }}
  table.craft-table th, table.craft-table td {{
    padding: 4px 10px; border-bottom: 1px solid #eee; text-align: left; white-space: nowrap;
  }}
  table.craft-table th {{ background: #f0f2f5; }}
  h3 .count {{ color: #999; font-weight: normal; font-size: 0.85em; }}
  .footer {{ color: #aaa; font-size: 0.8em; margin-top: 3em; }}
</style>
</head>
<body>
<h1>CRAFT annotation report</h1>
<p class="meta">Generated {timestamp} UTC &middot; craft v{version} &middot; {n_isoforms:,} iso</p>

<div class="kpis">{kpis}</div>

<h2>Functional-consequence cascade</h2>
<div class="grid">{cascade}</div>

<h2>Distributions</h2>
<div class="grid">{distributions}</div>

<h2>Notable findings</h2>
{notable}

<p class="footer">CRAFT (Coding Region Annotation From Templates). See docs/features.md.</p>
</body>
</html>
"""


def _pct(n: int, total: int) -> str:
    return f"{n / total * 100:.0f}%" if total else "0%"


def _kpis(df: pd.DataFrame) -> str:
    n = len(df)
    tiles: list[tuple[str, str]] = [(f"{n:,}", "isoforms")]
    if "nmd_status" in df.columns:
        tiles.append((_pct(int((df["nmd_status"] == "sensitive").sum()), n), "NMD-sensitive"))
    if "ptc_introduced" in df.columns:
        tiles.append((_pct(int(df["ptc_introduced"].sum()), n), "premature stop"))
    if "intron_retained_in_cds" in df.columns:
        tiles.append((_pct(int(df["intron_retained_in_cds"].sum()), n), "retained CDS intron"))
    if "completeness" in df.columns:
        alt = int((df["completeness"] == "alt_3prime_end").sum())
        tiles.append((_pct(alt, n), "alt 3' end (APA)"))
    cp = df.get("coding_potential_label")
    if cp is not None and (cp == "coding").any():
        tiles.append((_pct(int((cp == "coding").sum()), n), "coding (potential)"))
    return "".join(
        f'<div class="kpi"><div class="val">{v}</div><div class="lab">{lab}</div></div>'
        for v, lab in tiles
    )


def _figure_block(fig, *, include_js: bool) -> str:
    inner = fig.to_html(
        include_plotlyjs="inline" if include_js else False, full_html=False
    )
    return f"<div>{inner}</div>"


def _value_counts(df: pd.DataFrame, col: str) -> dict[str, int]:
    if col not in df.columns or df.empty:
        return {}
    return {str(k): int(v) for k, v in df[col].dropna().value_counts().items()}


def _cascade_fig(df: pd.DataFrame):
    n = len(df)
    coding_parent = int((~df["orf_outcome"].isin(["no_parent", "no_parent_cds"])).sum())
    resolved = int(df["resolved_orf_status"].isin(_RESOLVED_WITH_STOP).sum())
    intact = int((df["resolved_orf_status"] == "intact").sum())
    return funnel(
        [
            ("all isoforms", n),
            ("coding parent", coding_parent),
            ("ORF reconstructed", resolved),
            ("ORF intact", intact),
        ]
    )


def _small_table_html(df: pd.DataFrame, columns: list[str]) -> str:
    cols_present = [c for c in columns if c in df.columns]
    return df[cols_present].to_html(
        index=False, classes="craft-table", border=0, na_rep="-"
    )


def _top_nmd_sensitive(df: pd.DataFrame) -> tuple[str, int]:
    if "nmd_status" not in df.columns or "nmd_confidence" not in df.columns:
        return "", 0
    sub = df[(df["nmd_status"] == "sensitive") & (df["nmd_confidence"] == "high")]
    if sub.empty:
        return "", 0
    sort_col = "orf_confidence_score" if "orf_confidence_score" in sub.columns else None
    top = sub.nlargest(_NOTABLE_TOP_N, sort_col) if sort_col else sub.head(_NOTABLE_TOP_N)
    cols = [
        "transcript_id", "parent_gene_name", "completeness", "resolved_orf_status",
        "nmd_rule", "nmd_basis", "stop_to_last_junction_nt",
    ]
    return _small_table_html(top, cols), len(sub)


def _top_intron_retained(df: pd.DataFrame) -> tuple[str, int]:
    if "intron_retained_in_cds" not in df.columns:
        return "", 0
    sub = df[df["intron_retained_in_cds"]].copy()
    if sub.empty:
        return "", 0
    top = sub.head(_NOTABLE_TOP_N)
    cols = [
        "transcript_id", "parent_gene_name", "completeness", "resolved_orf_status",
        "nmd_status", "resolved_aa_length",
    ]
    return _small_table_html(top, cols), len(sub)


def _top_isoform_diversity(df: pd.DataFrame) -> tuple[str, int]:
    required = {"parent_gene_id", "parent_tx_id", "resolved_orf_status", "orf_confidence"}
    if not required.issubset(df.columns):
        return "", 0
    sub = df[df["orf_confidence"].isin(["high", "medium"])].copy()
    sub["parent_gene_id"] = sub["parent_gene_id"].astype(str)
    sub = sub[sub["parent_gene_id"] != ""]
    if sub.empty:
        return "", 0
    distinct = sub[["parent_gene_id", "parent_tx_id", "resolved_orf_status"]].drop_duplicates()
    diversity = (
        distinct.groupby("parent_gene_id").size().rename("n_functional_variants").reset_index()
    )
    if "parent_gene_name" in sub.columns:
        names = sub.drop_duplicates("parent_gene_id").set_index("parent_gene_id")[
            "parent_gene_name"
        ].to_dict()
        diversity["parent_gene_name"] = diversity["parent_gene_id"].map(names).fillna("")
    diversity = diversity.sort_values("n_functional_variants", ascending=False).head(_NOTABLE_TOP_N)
    total = int(diversity["n_functional_variants"].sum()) if not diversity.empty else 0
    cols = ["parent_gene_name", "parent_gene_id", "n_functional_variants"]
    return _small_table_html(diversity, cols), total


def _notable_findings_html(df: pd.DataFrame) -> str:
    sections: list[str] = []
    nmd_html, nmd_total = _top_nmd_sensitive(df)
    if nmd_html:
        sections.append(
            f'<h3>Top NMD-sensitive isoforms <span class="count">({nmd_total} '
            f'high-confidence)</span></h3>{nmd_html}'
        )
    ir_html, ir_total = _top_intron_retained(df)
    if ir_html:
        sections.append(
            f'<h3>Intron-retained-in-CDS isoforms <span class="count">({ir_total} '
            f'total)</span></h3>{ir_html}'
        )
    div_html, div_total = _top_isoform_diversity(df)
    if div_html:
        sections.append(
            f'<h3>Genes with most functional isoform diversity '
            f'<span class="count">(distinct parent_tx x resolved status; '
            f'{div_total} variants in top {_NOTABLE_TOP_N})</span></h3>{div_html}'
        )
    if not sections:
        return '<p class="meta">No notable findings to surface. See per_isoform.tsv.</p>'
    return "\n".join(sections)


def render(per_isoform: pd.DataFrame, output: Path) -> None:
    """Render the interactive HTML report."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    n_isoforms = len(per_isoform)

    figs = []  # (fig, needs_js)
    if not per_isoform.empty and "resolved_orf_status" in per_isoform.columns:
        figs.append(_cascade_fig(per_isoform))
    cascade_block = (
        _figure_block(figs[0], include_js=True) if figs else "<p class='meta'>No data.</p>"
    )

    dist_figs = [
        category_bar(
            _value_counts(per_isoform, "completeness"), "Completeness", _COMPLETENESS_ORDER
        ),
        category_bar(
            _value_counts(per_isoform, "resolved_orf_status"),
            "Resolved ORF status",
            _RESOLVED_ORDER,
        ),
        category_bar(_value_counts(per_isoform, "nmd_status"), "NMD status", _NMD_ORDER),
    ]
    cp = per_isoform.get("coding_potential_score")
    if cp is not None and pd.to_numeric(cp, errors="coerce").notna().any():
        dist_figs.append(histogram(cp, "Coding-potential score", threshold=0.5))
    u3 = per_isoform.get("utr3_length_delta_nt")
    if u3 is not None and pd.to_numeric(u3, errors="coerce").notna().any():
        dist_figs.append(histogram(u3, "3'UTR length delta vs parent (nt)", threshold=0))

    distributions = "\n".join(
        _figure_block(f, include_js=(figs == [] and i == 0)) for i, f in enumerate(dist_figs)
    )

    html = _TEMPLATE.format(
        timestamp=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M"),
        version=__version__,
        n_isoforms=n_isoforms,
        kpis=_kpis(per_isoform) if n_isoforms else "",
        cascade=cascade_block,
        distributions=distributions,
        notable=_notable_findings_html(per_isoform),
    )
    output.write_text(html)
