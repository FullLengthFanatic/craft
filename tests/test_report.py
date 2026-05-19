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
    assert "Notable findings" in html


def test_render_summary_counts_appear_in_html(tmp_path: Path) -> None:
    output = tmp_path / "report.html"
    render(_example_df(), output)
    html = output.read_text()
    assert "full_length" in html
    assert "truncated_5p" in html
    assert "novel_no_match" in html


def test_render_creates_nested_output_dir(tmp_path: Path) -> None:
    output = tmp_path / "a" / "b" / "report.html"
    render(_example_df(), output)
    assert output.exists()


def test_render_includes_plotly_js(tmp_path: Path) -> None:
    """Plotly JS should be inlined in the first figure block (self-contained)."""
    output = tmp_path / "report.html"
    render(_example_df(), output)
    html = output.read_text()
    assert "Plotly" in html


def test_render_includes_top_nmd_sensitive_section(tmp_path: Path) -> None:
    """t2 is sensitive + high confidence in the example fixture."""
    output = tmp_path / "report.html"
    render(_example_df(), output)
    html = output.read_text()
    assert "Top NMD-sensitive isoforms" in html
    assert "t2" in html  # should appear in the NMD table


def test_render_skips_empty_sections(tmp_path: Path) -> None:
    """No 'disrupted' rows in the fixture; that section should be absent."""
    output = tmp_path / "report.html"
    render(_example_df(), output)
    html = output.read_text()
    assert "Top ORF-disrupted isoforms" not in html
    # No gene_id in the fixture; diversity table should also be absent.
    assert "Genes with most isoform diversity" not in html


def test_render_includes_gene_diversity_section_when_gene_id_present(
    tmp_path: Path,
) -> None:
    df = _example_df()
    df["parent_gene_id"] = ["g1", "g1", ""]
    df["parent_gene_name"] = ["GENE1", "GENE1", ""]
    output = tmp_path / "report.html"
    render(df, output)
    html = output.read_text()
    assert "Genes with most functional isoform diversity" in html
    assert "GENE1" in html


def test_render_gene_diversity_collapses_duplicate_functional_variants(
    tmp_path: Path,
) -> None:
    """Two isoforms with the same (parent_tx_id, orf_outcome) count as one
    functional variant; differing on either dimension counts as two."""
    base = {
        "completeness": "full_length",
        "nmd_status": "escaped",
        "nmd_confidence": "high",
        "orf_confidence": "high",
        "orf_confidence_score": 1.0,
    }
    df = pd.DataFrame(
        [
            # gALPHA: two isos with identical (parent_tx_id, orf_outcome) -> 1 variant.
            {**base, "transcript_id": "t1a", "parent_gene_id": "gALPHA",
             "parent_gene_name": "GENEALPHA", "parent_tx_id": "t_ref_a",
             "orf_outcome": "propagated_intact"},
            {**base, "transcript_id": "t1b", "parent_gene_id": "gALPHA",
             "parent_gene_name": "GENEALPHA", "parent_tx_id": "t_ref_a",
             "orf_outcome": "propagated_intact"},
            # gBETA: two isos with different parent_tx_ids -> 2 variants.
            {**base, "transcript_id": "t2a", "parent_gene_id": "gBETA",
             "parent_gene_name": "GENEBETA", "parent_tx_id": "t_ref_b",
             "orf_outcome": "propagated_intact"},
            {**base, "transcript_id": "t2b", "parent_gene_id": "gBETA",
             "parent_gene_name": "GENEBETA", "parent_tx_id": "t_ref_c",
             "orf_outcome": "propagated_intact"},
        ]
    )
    output = tmp_path / "report.html"
    render(df, output)
    html = output.read_text()
    # Search inside the diversity table's HTML (after its heading) to avoid
    # matches in plotly's inlined JS, which is dumped earlier in the file.
    section_start = html.find("Genes with most functional isoform diversity")
    assert section_start != -1
    section = html[section_start:]
    alpha_pos = section.find("GENEALPHA")
    beta_pos = section.find("GENEBETA")
    assert alpha_pos != -1 and beta_pos != -1
    assert beta_pos < alpha_pos, "Gene with more functional variants should rank first"


def test_render_gene_diversity_filters_out_low_confidence(tmp_path: Path) -> None:
    """LOW-confidence isoforms should NOT count toward gene diversity."""
    df = pd.DataFrame(
        [
            {
                "transcript_id": "t1", "parent_gene_id": "g1",
                "parent_gene_name": "ONLYLOW", "parent_tx_id": "t_ref",
                "orf_outcome": "no_parent_cds", "orf_confidence": "low",
                "completeness": "full_length", "nmd_status": "not_applicable",
                "nmd_confidence": "none",
            },
        ]
    )
    output = tmp_path / "report.html"
    render(df, output)
    html = output.read_text()
    # No high/medium confidence isoforms -> gene diversity section absent.
    assert "Genes with most functional isoform diversity" not in html
    assert "ONLYLOW" not in html


def test_render_shows_fallback_when_no_findings(tmp_path: Path) -> None:
    """No NMD-sensitive, no disrupted, no gene_id -> fallback message."""
    df = pd.DataFrame(
        [
            {
                "transcript_id": "t1",
                "completeness": "full_length",
                "orf_outcome": "propagated_intact",
                "orf_confidence": "high",
                "nmd_status": "escaped",
                "nmd_confidence": "high",
            }
        ]
    )
    output = tmp_path / "report.html"
    render(df, output)
    html = output.read_text()
    assert "No notable findings" in html
    assert "per_isoform.tsv" in html


def test_bar_chart_uses_muted_color() -> None:
    """Bar charts should not use plotly's default bright blue."""
    fig = bar_chart({"a": 5, "b": 10}, "Test")
    marker_color = fig.data[0].marker.color
    assert marker_color == "#5b7a9d"
