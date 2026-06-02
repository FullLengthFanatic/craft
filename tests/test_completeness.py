"""Tests for craft.core.completeness."""

import pandas as pd
import pyranges as pr

from craft.core.completeness import Completeness, classify


def _exons(records: list[tuple]) -> pr.PyRanges:
    cols = ["Chromosome", "Start", "End", "Strand", "transcript_id"]
    df = pd.DataFrame(records, columns=cols)
    return pr.PyRanges(df)


def _class_of(result: pr.PyRanges, tx_id: str) -> Completeness:
    rows = result.df[result.df["transcript_id"] == tx_id]
    cats = rows["completeness"].unique()
    assert len(cats) == 1, f"Inconsistent classification for {tx_id}: {cats}"
    return Completeness(cats[0])


def _parent_of(result: pr.PyRanges, tx_id: str) -> str:
    rows = result.df[result.df["transcript_id"] == tx_id]
    parents = rows["parent_tx_id"].unique()
    assert len(parents) == 1
    return parents[0]


def test_identical_isoform_classified_as_full_length_plus_strand() -> None:
    ref = _exons(
        [
            ("chr1", 100, 200, "+", "t_ref"),
            ("chr1", 300, 400, "+", "t_ref"),
            ("chr1", 500, 600, "+", "t_ref"),
        ]
    )
    iso = _exons(
        [
            ("chr1", 100, 200, "+", "t1"),
            ("chr1", 300, 400, "+", "t1"),
            ("chr1", 500, 600, "+", "t1"),
        ]
    )
    result = classify(iso, ref)
    assert _class_of(result, "t1") == Completeness.FULL_LENGTH
    assert _parent_of(result, "t1") == "t_ref"


def test_identical_isoform_classified_as_full_length_minus_strand() -> None:
    ref = _exons(
        [
            ("chr1", 100, 200, "-", "t_ref"),
            ("chr1", 300, 400, "-", "t_ref"),
            ("chr1", 500, 600, "-", "t_ref"),
        ]
    )
    iso = _exons(
        [
            ("chr1", 100, 200, "-", "t1"),
            ("chr1", 300, 400, "-", "t1"),
            ("chr1", 500, 600, "-", "t1"),
        ]
    )
    result = classify(iso, ref)
    assert _class_of(result, "t1") == Completeness.FULL_LENGTH


def test_5prime_truncated_plus_strand() -> None:
    ref = _exons(
        [
            ("chr1", 100, 200, "+", "t_ref"),
            ("chr1", 300, 400, "+", "t_ref"),
            ("chr1", 500, 600, "+", "t_ref"),
        ]
    )
    iso = _exons(
        [
            ("chr1", 350, 400, "+", "t1"),
            ("chr1", 500, 600, "+", "t1"),
        ]
    )
    result = classify(iso, ref)
    assert _class_of(result, "t1") == Completeness.TRUNCATED_5P


def test_3prime_truncated_plus_strand() -> None:
    ref = _exons(
        [
            ("chr1", 100, 200, "+", "t_ref"),
            ("chr1", 300, 400, "+", "t_ref"),
            ("chr1", 500, 600, "+", "t_ref"),
        ]
    )
    iso = _exons(
        [
            ("chr1", 100, 200, "+", "t1"),
            ("chr1", 300, 400, "+", "t1"),
        ]
    )
    result = classify(iso, ref)
    assert _class_of(result, "t1") == Completeness.TRUNCATED_3P


def test_5prime_truncated_minus_strand() -> None:
    ref = _exons(
        [
            ("chr1", 100, 200, "-", "t_ref"),
            ("chr1", 300, 400, "-", "t_ref"),
            ("chr1", 500, 600, "-", "t_ref"),
        ]
    )
    iso = _exons(
        [
            ("chr1", 100, 200, "-", "t1"),
            ("chr1", 300, 400, "-", "t1"),
        ]
    )
    result = classify(iso, ref)
    assert _class_of(result, "t1") == Completeness.TRUNCATED_5P


def test_3prime_truncated_minus_strand() -> None:
    ref = _exons(
        [
            ("chr1", 100, 200, "-", "t_ref"),
            ("chr1", 300, 400, "-", "t_ref"),
            ("chr1", 500, 600, "-", "t_ref"),
        ]
    )
    iso = _exons(
        [
            ("chr1", 300, 400, "-", "t1"),
            ("chr1", 500, 600, "-", "t1"),
        ]
    )
    result = classify(iso, ref)
    assert _class_of(result, "t1") == Completeness.TRUNCATED_3P


def test_truncated_both_close_to_parent_boundaries() -> None:
    ref = _exons(
        [
            ("chr1", 100, 200, "+", "t_ref"),
            ("chr1", 300, 400, "+", "t_ref"),
            ("chr1", 500, 600, "+", "t_ref"),
            ("chr1", 700, 800, "+", "t_ref"),
            ("chr1", 900, 1000, "+", "t_ref"),
        ]
    )
    iso = _exons(
        [
            ("chr1", 180, 200, "+", "t1"),
            ("chr1", 300, 400, "+", "t1"),
            ("chr1", 500, 600, "+", "t1"),
            ("chr1", 700, 800, "+", "t1"),
            ("chr1", 900, 920, "+", "t1"),
        ]
    )
    result = classify(iso, ref)
    assert _class_of(result, "t1") == Completeness.TRUNCATED_BOTH


def test_internal_fragment_well_inside_parent() -> None:
    ref = _exons(
        [
            ("chr1", 100, 200, "+", "t_ref"),
            ("chr1", 300, 400, "+", "t_ref"),
            ("chr1", 500, 600, "+", "t_ref"),
            ("chr1", 700, 800, "+", "t_ref"),
            ("chr1", 900, 1000, "+", "t_ref"),
        ]
    )
    iso = _exons(
        [
            ("chr1", 500, 600, "+", "t1"),
            ("chr1", 700, 800, "+", "t1"),
        ]
    )
    result = classify(iso, ref)
    assert _class_of(result, "t1") == Completeness.INTERNAL_FRAGMENT


def test_no_overlap_classified_as_novel_no_match() -> None:
    ref = _exons(
        [
            ("chr1", 100, 200, "+", "t_ref"),
            ("chr1", 300, 400, "+", "t_ref"),
        ]
    )
    iso = _exons(
        [
            ("chr2", 100, 200, "+", "t_novel"),
            ("chr2", 300, 400, "+", "t_novel"),
        ]
    )
    result = classify(iso, ref)
    assert _class_of(result, "t_novel") == Completeness.NOVEL_NO_MATCH
    assert _parent_of(result, "t_novel") == ""


def test_picks_parent_with_most_shared_junctions() -> None:
    ref = _exons(
        [
            ("chr1", 100, 200, "+", "ref_a"),
            ("chr1", 300, 400, "+", "ref_a"),
            ("chr1", 500, 600, "+", "ref_a"),
            ("chr1", 100, 200, "+", "ref_b"),
            ("chr1", 300, 350, "+", "ref_b"),
            ("chr1", 500, 600, "+", "ref_b"),
        ]
    )
    iso = _exons(
        [
            ("chr1", 100, 200, "+", "t1"),
            ("chr1", 300, 400, "+", "t1"),
            ("chr1", 500, 600, "+", "t1"),
        ]
    )
    result = classify(iso, ref)
    assert _parent_of(result, "t1") == "ref_a"


def test_emits_per_isoform_metadata_columns() -> None:
    ref = _exons(
        [
            ("chr1", 100, 200, "+", "t_ref"),
            ("chr1", 300, 400, "+", "t_ref"),
            ("chr1", 500, 600, "+", "t_ref"),
        ]
    )
    iso = _exons(
        [
            ("chr1", 100, 200, "+", "t1"),
            ("chr1", 300, 400, "+", "t1"),
            ("chr1", 500, 600, "+", "t1"),
        ]
    )
    result = classify(iso, ref).df
    for col in ("completeness", "parent_tx_id", "shared_junctions", "parent_overlap_bp"):
        assert col in result.columns, f"Missing column: {col}"
    row = result[result["transcript_id"] == "t1"].iloc[0]
    assert int(row["shared_junctions"]) == 2
    assert int(row["parent_overlap_bp"]) > 0


def test_coding_aware_parent_tiebreak_prefers_cds_bearing() -> None:
    # Two reference transcripts with identical structure tie on shared junctions
    # and exon overlap; only t_cds carries a CDS.
    ref = _exons(
        [
            ("chr1", 100, 200, "+", "t_cds"),
            ("chr1", 300, 400, "+", "t_cds"),
            ("chr1", 100, 200, "+", "t_noncds"),
            ("chr1", 300, 400, "+", "t_noncds"),
        ]
    )
    iso = _exons([("chr1", 100, 200, "+", "t1"), ("chr1", 300, 400, "+", "t1")])
    result = classify(iso, ref, cds_tx_ids={"t_cds"}, prefer_coding_parent=True)
    assert _parent_of(result, "t1") == "t_cds"
    row = result.df[result.df["transcript_id"] == "t1"].iloc[0]
    assert bool(row["has_cds_bearing_parent"])


def test_has_cds_bearing_parent_false_without_cds_set() -> None:
    ref = _exons([("chr1", 100, 200, "+", "t_ref"), ("chr1", 300, 400, "+", "t_ref")])
    iso = _exons([("chr1", 100, 200, "+", "t1"), ("chr1", 300, 400, "+", "t1")])
    result = classify(iso, ref)
    row = result.df[result.df["transcript_id"] == "t1"].iloc[0]
    assert not bool(row["has_cds_bearing_parent"])
