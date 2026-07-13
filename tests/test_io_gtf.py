"""Tests for craft.io.gtf."""

import gzip
from pathlib import Path

import pyranges as pr
import pytest

from craft.io.gtf import load_isoforms, load_reference


@pytest.fixture
def synthetic_gtf(tmp_path: Path) -> Path:
    """Small GTF with two transcripts, one + strand on chr1 (3 exons + 3 CDS),
    one - strand on chr2 (2 exons + 1 CDS), plus assorted non-exon/CDS rows
    (gene, transcript, start_codon) to verify filtering."""
    rows = [
        'chr1\tCRAFT\tgene\t101\t600\t.\t+\t.\tgene_id "g1";',
        'chr1\tCRAFT\ttranscript\t101\t600\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
        'chr1\tCRAFT\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
        'chr1\tCRAFT\texon\t301\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
        'chr1\tCRAFT\texon\t501\t600\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
        'chr1\tCRAFT\tCDS\t151\t200\t.\t+\t0\tgene_id "g1"; transcript_id "t1";',
        'chr1\tCRAFT\tCDS\t301\t400\t.\t+\t0\tgene_id "g1"; transcript_id "t1";',
        'chr1\tCRAFT\tCDS\t501\t550\t.\t+\t0\tgene_id "g1"; transcript_id "t1";',
        'chr1\tCRAFT\tstart_codon\t151\t153\t.\t+\t0\tgene_id "g1"; transcript_id "t1";',
        'chr2\tCRAFT\tgene\t101\t400\t.\t-\t.\tgene_id "g2";',
        'chr2\tCRAFT\ttranscript\t101\t400\t.\t-\t.\tgene_id "g2"; transcript_id "t2";',
        'chr2\tCRAFT\texon\t101\t200\t.\t-\t.\tgene_id "g2"; transcript_id "t2";',
        'chr2\tCRAFT\texon\t301\t400\t.\t-\t.\tgene_id "g2"; transcript_id "t2";',
        'chr2\tCRAFT\tCDS\t301\t380\t.\t-\t0\tgene_id "g2"; transcript_id "t2";',
    ]
    path = tmp_path / "test.gtf"
    path.write_text("\n".join(rows) + "\n")
    return path


def test_load_isoforms_returns_pyranges(synthetic_gtf: Path) -> None:
    iso = load_isoforms(synthetic_gtf)
    assert isinstance(iso, pr.PyRanges)


def test_load_isoforms_keeps_only_exon_rows(synthetic_gtf: Path) -> None:
    iso = load_isoforms(synthetic_gtf)
    # 3 exons on t1 + 2 exons on t2 = 5
    assert len(iso) == 5


def test_load_isoforms_coordinates_are_zero_based_half_open(synthetic_gtf: Path) -> None:
    iso = load_isoforms(synthetic_gtf)
    df = iso.df.sort_values(["transcript_id", "Start"]).reset_index(drop=True)
    # First exon of t1: GTF 101-200 (1-based inclusive) -> PyRanges 100-200
    first = df.iloc[0]
    assert first["transcript_id"] == "t1"
    assert int(first["Start"]) == 100
    assert int(first["End"]) == 200


def test_load_isoforms_preserves_strand(synthetic_gtf: Path) -> None:
    iso = load_isoforms(synthetic_gtf)
    df = iso.df
    assert iso.stranded
    assert (df.loc[df["transcript_id"] == "t1", "Strand"] == "+").all()
    assert (df.loc[df["transcript_id"] == "t2", "Strand"] == "-").all()


def test_load_isoforms_has_minimum_columns(synthetic_gtf: Path) -> None:
    iso = load_isoforms(synthetic_gtf)
    required = {"Chromosome", "Start", "End", "Strand", "transcript_id"}
    assert required.issubset(set(iso.columns))


def test_load_reference_keeps_exon_and_cds_rows(synthetic_gtf: Path) -> None:
    ref = load_reference(synthetic_gtf)
    df = ref.df
    features = set(df["Feature"].unique())
    assert features == {"exon", "CDS", "start_codon"}
    # 5 exons + 4 CDS + explicit start codon = 10 rows
    assert len(ref) == 10


def test_load_reference_feature_column_preserved(synthetic_gtf: Path) -> None:
    ref = load_reference(synthetic_gtf)
    assert "Feature" in ref.columns


def test_load_reference_coordinates_match_gtf_conversion(synthetic_gtf: Path) -> None:
    ref = load_reference(synthetic_gtf)
    df = ref.df
    # CDS for t1 first segment: GTF 151-200 -> PyRanges 150-200
    cds = df[(df["transcript_id"] == "t1") & (df["Feature"] == "CDS")].sort_values("Start")
    first = cds.iloc[0]
    assert int(first["Start"]) == 150
    assert int(first["End"]) == 200


def test_loaded_isoforms_run_through_completeness_pipeline(synthetic_gtf: Path) -> None:
    """Integration: loaded PyRanges work with the existing core modules."""
    from craft.core.completeness import Completeness, classify

    iso = load_isoforms(synthetic_gtf)
    ref_full = load_reference(synthetic_gtf)
    # Build exon-only reference PyRanges with str-typed Strand
    ref_exons_df = (
        ref_full.df[ref_full.df["Feature"] == "exon"].drop(columns="Feature").copy()
    )
    ref_exons_df["Strand"] = ref_exons_df["Strand"].astype(str)
    ref_exons = pr.PyRanges(ref_exons_df)

    result = classify(iso, ref_exons)
    df = result.df
    # Iso == ref for both t1 and t2; both classified as FULL_LENGTH.
    t1_class = df.loc[df["transcript_id"] == "t1", "completeness"].iloc[0]
    t2_class = df.loc[df["transcript_id"] == "t2", "completeness"].iloc[0]
    assert Completeness(t1_class) == Completeness.FULL_LENGTH
    assert Completeness(t2_class) == Completeness.FULL_LENGTH


def test_load_isoforms_handles_gzipped_gtf(tmp_path: Path) -> None:
    gtf_text = 'chr1\tCRAFT\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    path = tmp_path / "test.gtf.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(gtf_text)
    iso = load_isoforms(path)
    assert len(iso) == 1
    df = iso.df
    assert df.iloc[0]["transcript_id"] == "t1"
    assert int(df.iloc[0]["Start"]) == 100
    assert int(df.iloc[0]["End"]) == 200
