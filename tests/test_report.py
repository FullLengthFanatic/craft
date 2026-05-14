"""Tests for craft.report.html and craft.report.plots."""

from pathlib import Path

import pandas as pd
import plotly.graph_objects as go

from craft.report.html import render
from craft.report.plots import bar_chart


def _example_df() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "transcript_id": "t1",
                "completeness": "full_length",
                "parent_tx_id": "t_ref",
                "orf_outcome": "propagated_intact",
                "orf_confidence": "high",
                "orf_confidence_score": 1.0,
                "nmd_status": "escaped",
                "nmd_rule": "stop_in_last_exon",
                "stop_to_last_junction_nt": 0,
                "last_exon_length_nt": 100,
                "nmd_confidence": "high",
                "iso_utr3_length_nt": 50,
                "parent_utr3_length_nt": 50,
                "utr3_length_delta_nt": 0,
                "utr3_length_delta_pct": 0.0,
                "polya_signal_motif": "",
                "polya_signal_distance_nt": None,
                "denovo_orf_found": False,
                "denovo_cds_bp": 0,
                "denovo_orf_aa_length": 0,
                "denovo_start_codon": "",
                "denovo_stop_codon": "",
                "denovo_cds_intervals": [],
                "propagated_cds_intervals": [("chr1", 150, 200, "+")],
                "shared_junctions": 2,
                "parent_overlap_bp": 300,
                "propagated_cds_bp": 200,
                "parent_cds_bp": 200,
                "start_codon_covered": True,
                "stop_codon_covered": True,
            },
            {
                "transcript_id": "t2",
                "completeness": "truncated_5p",
                "parent_tx_id": "t_ref",
                "orf_outcome": "propagated_intact",
                "orf_confidence": "high",
                "orf_confidence_score": 0.9,
                "nmd_status": "sensitive",
                "nmd_rule": "ptc_50nt_rule",
                "stop_to_last_junction_nt": 120,
                "last_exon_length_nt": 100,
                "nmd_confidence": "high",
                "iso_utr3_length_nt": 50,
                "parent_utr3_length_nt": 50,
                "utr3_length_delta_nt": 0,
                "utr3_length_delta_pct": 0.0,
                "polya_signal_motif": "AATAAA",
                "polya_signal_distance_nt": 20,
                "denovo_orf_found": False,
                "denovo_cds_bp": 0,
                "denovo_orf_aa_length": 0,
                "denovo_start_codon": "",
                "denovo_stop_codon": "",
                "denovo_cds_intervals": [],
                "propagated_cds_intervals": [("chr1", 200, 350, "+")],
                "shared_junctions": 1,
                "parent_overlap_bp": 200,
                "propagated_cds_bp": 150,
                "parent_cds_bp": 200,
                "start_codon_covered": True,
                "stop_codon_covered": True,
            },
            {
                "transcript_id": "t_novel",
                "completeness": "novel_no_match",
                "parent_tx_id": "",
                "orf_outcome": "no_parent",
                "orf_confidence": "low",
                "orf_confidence_score": 0.25,
                "nmd_status": "not_applicable",
                "nmd_rule": "",
                "stop_to_last_junction_nt": None,
                "last_exon_length_nt": None,
                "nmd_confidence": "none",
                "iso_utr3_length_nt": None,
                "parent_utr3_length_nt": None,
                "utr3_length_delta_nt": None,
                "utr3_length_delta_pct": None,
                "polya_signal_motif": "",
                "polya_signal_distance_nt": None,
                "denovo_orf_found": True,
                "denovo_cds_bp": 180,
                "denovo_orf_aa_length": 60,
                "denovo_start_codon": "ATG",
                "denovo_stop_codon": "TAA",
                "denovo_cds_intervals": [("chr2", 0, 180, "+")],
                "propagated_cds_intervals": [],
                "shared_junctions": 0,
                "parent_overlap_bp": 0,
                "propagated_cds_bp": 0,
                "parent_cds_bp": 0,
                "start_codon_covered": False,
                "stop_codon_covered": False,
            },
        ]
    )


def test_bar_chart_returns_figure() -> None:
    fig = bar_chart({"full_length": 10, "truncated_5p": 5}, "Test")
    assert isinstance(fig, go.Figure)
    assert fig.layout.title.text == "Test"


def test_bar_chart_handles_empty_counts() -> None:
    fig = bar_chart({}, "Empty")
    assert isinstance(fig, go.Figure)
    assert fig.layout.title.text == "Empty"


def test_bar_chart_sorts_descending() -> None:
    fig = bar_chart({"a": 5, "b": 20, "c": 10}, "Sort")
    xs = list(fig.data[0].x)
    ys = list(fig.data[0].y)
    assert xs == ["b", "c", "a"]
    assert ys == [20, 10, 5]


def test_render_creates_html_file(tmp_path: Path) -> None:
    output = tmp_path / "report.html"
    render(_example_df(), output)
    assert output.exists()
    assert output.stat().st_size > 0


def test_render_contains_expected_sections(tmp_path: Path) -> None:
    output = tmp_path / "report.html"
    render(_example_df(), output)
    html = output.read_text()
    assert "CRAFT annotation report" in html
    assert "Summary" in html
    assert "Distributions" in html
    assert "Per-isoform annotations" in html


def test_render_includes_transcript_ids_in_table(tmp_path: Path) -> None:
    output = tmp_path / "report.html"
    render(_example_df(), output)
    html = output.read_text()
    assert "t1" in html
    assert "t2" in html
    assert "t_novel" in html


def test_render_summary_counts_appear_in_html(tmp_path: Path) -> None:
    output = tmp_path / "report.html"
    render(_example_df(), output)
    html = output.read_text()
    # Two isoforms are full_length / truncated_5p; one is novel_no_match.
    assert "full_length" in html
    assert "truncated_5p" in html
    assert "novel_no_match" in html


def test_render_drops_list_columns_from_table(tmp_path: Path) -> None:
    output = tmp_path / "report.html"
    render(_example_df(), output)
    html = output.read_text()
    # The list columns are dropped from the rendered table, but the column
    # headers should not appear at all.
    assert "propagated_cds_intervals" not in html
    assert "denovo_cds_intervals" not in html


def test_render_creates_nested_output_dir(tmp_path: Path) -> None:
    output = tmp_path / "a" / "b" / "report.html"
    render(_example_df(), output)
    assert output.exists()


def test_render_includes_plotly_js(tmp_path: Path) -> None:
    """Plotly JS should be inlined in the first figure block (self-contained)."""
    output = tmp_path / "report.html"
    render(_example_df(), output)
    html = output.read_text()
    # Plotly's inlined JS contains the global Plotly namespace symbol.
    assert "Plotly" in html


def test_render_row_cap_truncates_large_tables(tmp_path: Path) -> None:
    # Make a 20-row df, cap at 5.
    base = _example_df().iloc[0:1]
    big = pd.concat([base.assign(transcript_id=f"t{i}") for i in range(20)], ignore_index=True)
    output = tmp_path / "report.html"
    render(big, output, row_cap=5)
    html = output.read_text()
    assert "Showing first 5 of 20" in html
