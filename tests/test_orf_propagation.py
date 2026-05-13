"""Tests for craft.core.orf.propagation."""

import pandas as pd
import pyranges as pr

from craft.core.completeness import classify
from craft.core.orf.propagation import ORFOutcome, propagate


def _exons(records: list[tuple]) -> pr.PyRanges:
    cols = ["Chromosome", "Start", "End", "Strand", "transcript_id"]
    df = pd.DataFrame(records, columns=cols)
    return pr.PyRanges(df)


def _reference(exon_records: list[tuple], cds_records: list[tuple]) -> pr.PyRanges:
    cols = ["Chromosome", "Start", "End", "Strand", "transcript_id", "Feature"]
    rows = [(*r, "exon") for r in exon_records] + [(*r, "CDS") for r in cds_records]
    df = pd.DataFrame(rows, columns=cols)
    return pr.PyRanges(df)


def _classified(iso: pr.PyRanges, ref: pr.PyRanges) -> pr.PyRanges:
    ref_exon_df = ref.df[ref.df["Feature"] == "exon"].drop(columns="Feature").copy()
    ref_exon_df["Strand"] = ref_exon_df["Strand"].astype(str)
    return classify(iso, pr.PyRanges(ref_exon_df))


def _outcome(result: pd.DataFrame, tx: str) -> ORFOutcome:
    row = result[result["transcript_id"] == tx].iloc[0]
    return ORFOutcome(row["orf_outcome"])


def test_identical_iso_and_parent_propagated_intact_plus_strand() -> None:
    ref = _reference(
        exon_records=[
            ("chr1", 100, 200, "+", "t_ref"),
            ("chr1", 300, 400, "+", "t_ref"),
            ("chr1", 500, 600, "+", "t_ref"),
        ],
        cds_records=[
            ("chr1", 150, 200, "+", "t_ref"),
            ("chr1", 300, 400, "+", "t_ref"),
            ("chr1", 500, 550, "+", "t_ref"),
        ],
    )
    iso = _exons(
        [
            ("chr1", 100, 200, "+", "t1"),
            ("chr1", 300, 400, "+", "t1"),
            ("chr1", 500, 600, "+", "t1"),
        ]
    )
    result = propagate(_classified(iso, ref), ref)
    row = result[result["transcript_id"] == "t1"].iloc[0]
    assert ORFOutcome(row["orf_outcome"]) == ORFOutcome.PROPAGATED_INTACT
    assert row["propagated_cds_bp"] == 200
    assert row["parent_cds_bp"] == 200
    assert row["start_codon_covered"]
    assert row["stop_codon_covered"]
    assert row["propagated_cds_intervals"] == [
        ("chr1", 150, 200, "+"),
        ("chr1", 300, 400, "+"),
        ("chr1", 500, 550, "+"),
    ]


def test_identical_iso_and_parent_propagated_intact_minus_strand() -> None:
    ref = _reference(
        exon_records=[
            ("chr1", 100, 200, "-", "t_ref"),
            ("chr1", 300, 400, "-", "t_ref"),
            ("chr1", 500, 600, "-", "t_ref"),
        ],
        cds_records=[
            ("chr1", 150, 200, "-", "t_ref"),
            ("chr1", 300, 400, "-", "t_ref"),
            ("chr1", 500, 550, "-", "t_ref"),
        ],
    )
    iso = _exons(
        [
            ("chr1", 100, 200, "-", "t1"),
            ("chr1", 300, 400, "-", "t1"),
            ("chr1", 500, 600, "-", "t1"),
        ]
    )
    result = propagate(_classified(iso, ref), ref)
    assert _outcome(result, "t1") == ORFOutcome.PROPAGATED_INTACT


def test_5prime_truncation_with_start_codon_preserved_is_intact() -> None:
    ref = _reference(
        exon_records=[
            ("chr1", 100, 200, "+", "t_ref"),
            ("chr1", 300, 400, "+", "t_ref"),
            ("chr1", 500, 600, "+", "t_ref"),
        ],
        cds_records=[
            ("chr1", 150, 200, "+", "t_ref"),
            ("chr1", 300, 400, "+", "t_ref"),
            ("chr1", 500, 550, "+", "t_ref"),
        ],
    )
    iso = _exons(
        [
            ("chr1", 140, 200, "+", "t1"),
            ("chr1", 300, 400, "+", "t1"),
            ("chr1", 500, 600, "+", "t1"),
        ]
    )
    result = propagate(_classified(iso, ref), ref)
    assert _outcome(result, "t1") == ORFOutcome.PROPAGATED_INTACT


def test_5prime_truncation_cuts_off_start_codon_yields_start_lost() -> None:
    ref = _reference(
        exon_records=[
            ("chr1", 100, 200, "+", "t_ref"),
            ("chr1", 300, 400, "+", "t_ref"),
            ("chr1", 500, 600, "+", "t_ref"),
        ],
        cds_records=[
            ("chr1", 150, 200, "+", "t_ref"),
            ("chr1", 300, 400, "+", "t_ref"),
            ("chr1", 500, 550, "+", "t_ref"),
        ],
    )
    iso = _exons(
        [
            ("chr1", 170, 200, "+", "t1"),
            ("chr1", 300, 400, "+", "t1"),
            ("chr1", 500, 600, "+", "t1"),
        ]
    )
    result = propagate(_classified(iso, ref), ref)
    assert _outcome(result, "t1") == ORFOutcome.START_LOST


def test_3prime_truncation_past_stop_codon_yields_stop_not_observed() -> None:
    ref = _reference(
        exon_records=[
            ("chr1", 100, 200, "+", "t_ref"),
            ("chr1", 300, 400, "+", "t_ref"),
            ("chr1", 500, 600, "+", "t_ref"),
        ],
        cds_records=[
            ("chr1", 150, 200, "+", "t_ref"),
            ("chr1", 300, 400, "+", "t_ref"),
            ("chr1", 500, 550, "+", "t_ref"),
        ],
    )
    iso = _exons(
        [
            ("chr1", 100, 200, "+", "t1"),
            ("chr1", 300, 400, "+", "t1"),
            ("chr1", 500, 530, "+", "t1"),
        ]
    )
    result = propagate(_classified(iso, ref), ref)
    assert _outcome(result, "t1") == ORFOutcome.STOP_NOT_OBSERVED


def test_exon_skip_in_cds_region_yields_disrupted() -> None:
    ref = _reference(
        exon_records=[
            ("chr1", 100, 200, "+", "t_ref"),
            ("chr1", 300, 400, "+", "t_ref"),
            ("chr1", 500, 600, "+", "t_ref"),
            ("chr1", 700, 800, "+", "t_ref"),
        ],
        cds_records=[
            ("chr1", 150, 200, "+", "t_ref"),
            ("chr1", 300, 400, "+", "t_ref"),
            ("chr1", 500, 600, "+", "t_ref"),
            ("chr1", 700, 750, "+", "t_ref"),
        ],
    )
    iso = _exons(
        [
            ("chr1", 100, 200, "+", "t1"),
            ("chr1", 500, 600, "+", "t1"),
            ("chr1", 700, 800, "+", "t1"),
        ]
    )
    result = propagate(_classified(iso, ref), ref)
    assert _outcome(result, "t1") == ORFOutcome.DISRUPTED
    row = result[result["transcript_id"] == "t1"].iloc[0]
    assert row["propagated_cds_bp"] < row["parent_cds_bp"]


def test_novel_no_match_yields_no_parent_outcome() -> None:
    ref = _reference(
        exon_records=[("chr1", 100, 200, "+", "t_ref")],
        cds_records=[("chr1", 150, 200, "+", "t_ref")],
    )
    iso = _exons([("chr2", 100, 200, "+", "t_novel")])
    result = propagate(_classified(iso, ref), ref)
    assert _outcome(result, "t_novel") == ORFOutcome.NO_PARENT
    row = result[result["transcript_id"] == "t_novel"].iloc[0]
    assert row["propagated_cds_bp"] == 0
    assert row["propagated_cds_intervals"] == []


def test_parent_without_cds_records_yields_no_parent_cds() -> None:
    ref = _reference(
        exon_records=[
            ("chr1", 100, 200, "+", "t_ref"),
            ("chr1", 300, 400, "+", "t_ref"),
        ],
        cds_records=[],
    )
    iso = _exons(
        [
            ("chr1", 100, 200, "+", "t1"),
            ("chr1", 300, 400, "+", "t1"),
        ]
    )
    result = propagate(_classified(iso, ref), ref)
    assert _outcome(result, "t1") == ORFOutcome.NO_PARENT_CDS


def test_minus_strand_loss_of_5prime_exon_cuts_off_start_codon() -> None:
    # On - strand, 5' is at high coord. Start codon is at parent_cds["End"].max() - 1.
    ref = _reference(
        exon_records=[
            ("chr1", 100, 200, "-", "t_ref"),
            ("chr1", 300, 400, "-", "t_ref"),
            ("chr1", 500, 600, "-", "t_ref"),
        ],
        cds_records=[
            ("chr1", 150, 200, "-", "t_ref"),
            ("chr1", 300, 400, "-", "t_ref"),
            ("chr1", 500, 550, "-", "t_ref"),  # start codon area at 549 on - strand
        ],
    )
    iso = _exons(
        [
            ("chr1", 100, 200, "-", "t1"),
            ("chr1", 300, 400, "-", "t1"),
        ]
    )
    result = propagate(_classified(iso, ref), ref)
    assert _outcome(result, "t1") == ORFOutcome.START_LOST
