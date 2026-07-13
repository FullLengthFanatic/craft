"""Tests for craft.report.html and craft.report.plots."""

from pathlib import Path

import pandas as pd
import plotly.graph_objects as go

from craft.report.html import render
from craft.report.plots import category_bar, funnel, histogram


def _row(tx, completeness, outcome, resolved, nmd, **kw):
    base = {
        "transcript_id": tx,
        "completeness": completeness,
        "parent_tx_id": kw.get("parent_tx_id", "t_ref"),
        "parent_gene_id": kw.get("parent_gene_id", ""),
        "parent_gene_name": kw.get("parent_gene_name", ""),
        "orf_outcome": outcome,
        "resolved_orf_status": resolved,
        "orf_confidence": kw.get("orf_confidence", "high"),
        "orf_confidence_score": kw.get("orf_confidence_score", 1.0),
        "nmd_status": nmd,
        "nmd_rule": kw.get("nmd_rule", ""),
        "nmd_confidence": kw.get("nmd_confidence", "high"),
        "nmd_basis": kw.get("nmd_basis", "resolved"),
        "stop_to_last_junction_nt": kw.get("stop_to_last_junction_nt", 0),
        "ptc_introduced": kw.get("ptc_introduced", False),
        "intron_retained_in_cds": kw.get("intron_retained_in_cds", False),
        "resolved_aa_length": kw.get("resolved_aa_length", 100),
        "coding_potential_score": kw.get("coding_potential_score", 0.9),
        "coding_potential_label": kw.get("coding_potential_label", "coding"),
        "utr3_length_delta_nt": kw.get("utr3_length_delta_nt", 0),
    }
    return base


def _example_df() -> pd.DataFrame:
    return pd.DataFrame(
        [
            _row("t1", "full_length", "propagated_intact", "intact", "escaped"),
            _row(
                "t2", "truncated_5p", "disrupted", "ptc_premature", "sensitive",
                nmd_rule="ptc_50nt_rule", ptc_introduced=True, stop_to_last_junction_nt=120,
            ),
            _row(
                "t3", "full_length", "disrupted", "ptc_intron_retained", "sensitive",
                intron_retained_in_cds=True, ptc_introduced=True, coding_potential_label="coding",
            ),
            _row(
                "t_novel", "novel_no_match", "no_parent", "resolution_failed", "not_applicable",
                parent_tx_id="", orf_confidence="low", nmd_confidence="none", nmd_basis="none",
                coding_potential_label="noncoding", coding_potential_score=0.1,
                utr3_length_delta_nt=None,
            ),
        ]
    )


def test_category_bar_returns_figure_with_semantic_color() -> None:
    fig = category_bar({"sensitive": 10, "escaped": 5}, "NMD")
    assert isinstance(fig, go.Figure)
    assert fig.layout.title.text == "NMD"
    # 'sensitive' should be the warm/red semantic color, not a flat default.
    assert "#d1495b" in list(fig.data[0].marker.color)


def test_category_bar_handles_empty() -> None:
    assert isinstance(category_bar({}, "Empty"), go.Figure)


def test_funnel_and_histogram_return_figures() -> None:
    assert isinstance(funnel([("all", 100), ("coding parent", 60)]), go.Figure)
    assert isinstance(histogram(pd.Series([0.1, 0.9, 0.5]), "scores", threshold=0.5), go.Figure)


def test_render_creates_file_and_sections(tmp_path: Path) -> None:
    output = tmp_path / "report.html"
    render(_example_df(), output)
    assert output.exists() and output.stat().st_size > 0
    html = output.read_text()
    assert "CRAFT annotation report" in html
    assert "Functional-consequence cascade" in html
    assert "Distributions" in html
    assert "Notable findings" in html
    assert "Plotly" in html  # inlined JS, self-contained


def test_render_kpis_and_distribution_categories(tmp_path: Path) -> None:
    output = tmp_path / "report.html"
    render(_example_df(), output)
    html = output.read_text()
    assert "predicted NMD candidate" in html  # KPI label
    assert "full_length" in html
    assert "resolution_failed" in html


def test_render_notable_nmd_and_intron_retained(tmp_path: Path) -> None:
    output = tmp_path / "report.html"
    render(_example_df(), output)
    html = output.read_text()
    assert "Top predicted NMD candidates" in html
    assert "Intron-retained-in-CDS isoforms" in html
    assert "t3" in html  # the IR isoform


def test_render_gene_diversity_when_gene_id_present(tmp_path: Path) -> None:
    df = _example_df()
    df["parent_gene_id"] = ["g1", "g1", "g1", ""]
    df["parent_gene_name"] = ["GENE1", "GENE1", "GENE1", ""]
    output = tmp_path / "report.html"
    render(df, output)
    html = output.read_text()
    assert "Genes with most functional isoform diversity" in html
    assert "GENE1" in html


def test_render_celltype_as_nmd_panel(tmp_path: Path) -> None:
    celltype = pd.DataFrame(
        {
            "cell_group": ["neuron", "glia"],
            "transcript_id": ["t2", "t3"],
            "parent_gene_name": ["GENE2", "GENE3"],
            "nmd_rule": ["ptc_50nt_rule", "long_exon"],
            "recurrence_score": [0.98, 0.96],
            "molecules_in_group": [120.0, 80.0],
            "frac_of_group": [0.4, 0.3],
        }
    )
    output = tmp_path / "report.html"
    render(_example_df(), output, celltype_as_nmd=celltype)
    html = output.read_text()
    assert "Cell-type AS-NMD map" in html
    assert "neuron" in html and "GENE2" in html


def test_render_celltype_panel_absent_when_no_data(tmp_path: Path) -> None:
    output = tmp_path / "report.html"
    render(_example_df(), output)  # no celltype arg
    assert "Cell-type AS-NMD map" not in output.read_text()


def test_render_nested_dir_and_fallback(tmp_path: Path) -> None:
    out = tmp_path / "a" / "b" / "report.html"
    render(
        pd.DataFrame([_row("t1", "full_length", "propagated_intact", "intact", "escaped")]),
        out,
    )
    assert out.exists()
    html = out.read_text()
    # single intact/escaped isoform, no gene_id -> no notable sections.
    assert "No notable findings" in html
