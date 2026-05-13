"""Tests for craft.core.nmd."""

import pandas as pd
import pyranges as pr

from craft.core.nmd import NMDStatus, predict
from craft.core.orf.confidence import ORFConfidence
from craft.core.orf.propagation import ORFOutcome


def _iso_pr(records: list[tuple]) -> pr.PyRanges:
    cols = ["Chromosome", "Start", "End", "Strand", "transcript_id"]
    df = pd.DataFrame(records, columns=cols)
    return pr.PyRanges(df)


def _prop_row(
    tx_id: str,
    outcome: ORFOutcome,
    intervals: list[tuple],
    cds_bp: int = 0,
    parent_bp: int = 0,
    start_covered: bool = True,
    stop_covered: bool = True,
    parent_tx_id: str = "t_ref",
) -> dict:
    return {
        "transcript_id": tx_id,
        "parent_tx_id": parent_tx_id,
        "orf_outcome": outcome.value,
        "propagated_cds_bp": cds_bp,
        "parent_cds_bp": parent_bp,
        "start_codon_covered": start_covered,
        "stop_codon_covered": stop_covered,
        "propagated_cds_intervals": intervals,
    }


def _row(result: pd.DataFrame, tx: str) -> pd.Series:
    return result[result["transcript_id"] == tx].iloc[0]


def test_stop_in_last_exon_plus_strand_escapes() -> None:
    iso = _iso_pr(
        [
            ("chr1", 100, 200, "+", "t1"),
            ("chr1", 300, 400, "+", "t1"),
        ]
    )
    prop = pd.DataFrame(
        [
            _prop_row(
                "t1",
                ORFOutcome.PROPAGATED_INTACT,
                intervals=[("chr1", 150, 200, "+"), ("chr1", 300, 380, "+")],
                cds_bp=130,
                parent_bp=130,
            )
        ]
    )
    row = _row(predict(iso, prop), "t1")
    assert row["nmd_status"] == NMDStatus.ESCAPED.value
    assert row["nmd_rule"] == "stop_in_last_exon"


def test_stop_far_upstream_of_last_junction_is_sensitive() -> None:
    iso = _iso_pr(
        [
            ("chr1", 100, 300, "+", "t1"),
            ("chr1", 400, 500, "+", "t1"),
            ("chr1", 600, 700, "+", "t1"),
        ]
    )
    # CDS bp 150 (avoids start_proximal). Stop at 249 (in first exon).
    # distance to last junction = (300-249) + (500-400) = 51 + 100 = 151
    prop = pd.DataFrame(
        [
            _prop_row(
                "t1",
                ORFOutcome.PROPAGATED_INTACT,
                intervals=[("chr1", 100, 250, "+")],
                cds_bp=150,
                parent_bp=150,
            )
        ]
    )
    row = _row(predict(iso, prop), "t1")
    assert row["nmd_status"] == NMDStatus.SENSITIVE.value
    assert row["nmd_rule"] == "ptc_50nt_rule"
    assert int(row["stop_to_last_junction_nt"]) == 151


def test_stop_within_50nt_of_last_junction_escapes() -> None:
    iso = _iso_pr(
        [
            ("chr1", 100, 300, "+", "t1"),
            ("chr1", 400, 500, "+", "t1"),
            ("chr1", 600, 700, "+", "t1"),
        ]
    )
    # Stop at 495 (in middle exon). distance = 500-495 = 5 <= 50 -> escape.
    prop = pd.DataFrame(
        [
            _prop_row(
                "t1",
                ORFOutcome.PROPAGATED_INTACT,
                intervals=[("chr1", 200, 300, "+"), ("chr1", 400, 496, "+")],
                cds_bp=196,
                parent_bp=196,
            )
        ]
    )
    row = _row(predict(iso, prop), "t1")
    assert row["nmd_status"] == NMDStatus.ESCAPED.value
    assert row["nmd_rule"] == "within_50nt_of_last_junction"


def test_single_exon_iso_escapes() -> None:
    iso = _iso_pr([("chr1", 100, 500, "+", "t1")])
    prop = pd.DataFrame(
        [
            _prop_row(
                "t1",
                ORFOutcome.PROPAGATED_INTACT,
                intervals=[("chr1", 100, 400, "+")],
                cds_bp=300,
                parent_bp=300,
            )
        ]
    )
    row = _row(predict(iso, prop), "t1")
    assert row["nmd_status"] == NMDStatus.ESCAPED.value
    assert row["nmd_rule"] == "stop_in_last_exon"


def test_start_proximal_short_orf_escapes() -> None:
    iso = _iso_pr(
        [
            ("chr1", 100, 200, "+", "t1"),
            ("chr1", 300, 400, "+", "t1"),
            ("chr1", 500, 600, "+", "t1"),
            ("chr1", 700, 800, "+", "t1"),
        ]
    )
    # CDS bp 100 < 150. Stop at 349 (middle exon idx 1). distance > 50.
    prop = pd.DataFrame(
        [
            _prop_row(
                "t1",
                ORFOutcome.PROPAGATED_INTACT,
                intervals=[("chr1", 150, 200, "+"), ("chr1", 300, 350, "+")],
                cds_bp=100,
                parent_bp=100,
            )
        ]
    )
    row = _row(predict(iso, prop), "t1")
    assert row["nmd_status"] == NMDStatus.ESCAPED.value
    assert row["nmd_rule"] == "start_proximal"


def test_long_last_exon_escapes() -> None:
    iso = _iso_pr(
        [
            ("chr1", 100, 300, "+", "t1"),
            ("chr1", 400, 500, "+", "t1"),
            ("chr1", 600, 1200, "+", "t1"),  # last exon 600 bp
        ]
    )
    # Stop at 449 (middle exon). distance = 500-449 = 51 > 50. CDS bp 250 >= 150.
    # Last exon length 600 > 400 -> long_last_exon escape.
    prop = pd.DataFrame(
        [
            _prop_row(
                "t1",
                ORFOutcome.PROPAGATED_INTACT,
                intervals=[("chr1", 100, 300, "+"), ("chr1", 400, 450, "+")],
                cds_bp=250,
                parent_bp=250,
            )
        ]
    )
    row = _row(predict(iso, prop), "t1")
    assert row["nmd_status"] == NMDStatus.ESCAPED.value
    assert row["nmd_rule"] == "long_last_exon"


def test_minus_strand_stop_in_last_exon_escapes() -> None:
    iso = _iso_pr(
        [
            ("chr1", 100, 200, "-", "t1"),
            ("chr1", 300, 400, "-", "t1"),
            ("chr1", 500, 600, "-", "t1"),
        ]
    )
    # - strand: last exon (transcript order) = (100, 200). Stop = min(Start) = 150 -> in last exon.
    prop = pd.DataFrame(
        [
            _prop_row(
                "t1",
                ORFOutcome.PROPAGATED_INTACT,
                intervals=[
                    ("chr1", 150, 200, "-"),
                    ("chr1", 300, 400, "-"),
                    ("chr1", 500, 550, "-"),
                ],
                cds_bp=200,
                parent_bp=200,
            )
        ]
    )
    row = _row(predict(iso, prop), "t1")
    assert row["nmd_status"] == NMDStatus.ESCAPED.value
    assert row["nmd_rule"] == "stop_in_last_exon"


def test_minus_strand_stop_far_upstream_is_sensitive() -> None:
    iso = _iso_pr(
        [
            ("chr1", 100, 200, "-", "t1"),
            ("chr1", 300, 400, "-", "t1"),
            ("chr1", 500, 700, "-", "t1"),
            ("chr1", 800, 900, "-", "t1"),
        ]
    )
    # transcript order: (800,900) -> (500,700) -> (300,400) -> (100,200) (last)
    # stop = min(Start of CDS) = 600 in (500, 700) idx 2
    # distance = (600-500) + length((300,400)) = 100 + 100 = 200
    prop = pd.DataFrame(
        [
            _prop_row(
                "t1",
                ORFOutcome.PROPAGATED_INTACT,
                intervals=[
                    ("chr1", 600, 700, "-"),
                    ("chr1", 800, 850, "-"),
                ],
                cds_bp=150,
                parent_bp=150,
            )
        ]
    )
    row = _row(predict(iso, prop), "t1")
    assert row["nmd_status"] == NMDStatus.SENSITIVE.value
    assert row["nmd_rule"] == "ptc_50nt_rule"
    assert int(row["stop_to_last_junction_nt"]) == 200


def test_stop_not_observed_is_not_applicable() -> None:
    iso = _iso_pr(
        [
            ("chr1", 100, 200, "+", "t1"),
            ("chr1", 300, 400, "+", "t1"),
        ]
    )
    prop = pd.DataFrame(
        [
            _prop_row(
                "t1",
                ORFOutcome.STOP_NOT_OBSERVED,
                intervals=[("chr1", 150, 200, "+")],
                cds_bp=50,
                parent_bp=200,
                stop_covered=False,
            )
        ]
    )
    row = _row(predict(iso, prop), "t1")
    assert row["nmd_status"] == NMDStatus.NOT_APPLICABLE.value
    assert row["nmd_confidence"] == ORFConfidence.NONE.value


def test_start_lost_is_not_applicable() -> None:
    iso = _iso_pr([("chr1", 100, 200, "+", "t1")])
    prop = pd.DataFrame(
        [
            _prop_row(
                "t1",
                ORFOutcome.START_LOST,
                intervals=[],
                cds_bp=0,
                parent_bp=0,
                start_covered=False,
            )
        ]
    )
    row = _row(predict(iso, prop), "t1")
    assert row["nmd_status"] == NMDStatus.NOT_APPLICABLE.value


def test_no_parent_is_not_applicable() -> None:
    iso = _iso_pr([("chr1", 100, 200, "+", "t1")])
    prop = pd.DataFrame(
        [
            _prop_row(
                "t1",
                ORFOutcome.NO_PARENT,
                intervals=[],
                cds_bp=0,
                parent_bp=0,
                start_covered=False,
                stop_covered=False,
                parent_tx_id="",
            )
        ]
    )
    row = _row(predict(iso, prop), "t1")
    assert row["nmd_status"] == NMDStatus.NOT_APPLICABLE.value


def test_disrupted_outcome_yields_medium_confidence() -> None:
    iso = _iso_pr(
        [
            ("chr1", 100, 200, "+", "t1"),
            ("chr1", 300, 400, "+", "t1"),
        ]
    )
    prop = pd.DataFrame(
        [
            _prop_row(
                "t1",
                ORFOutcome.DISRUPTED,
                intervals=[("chr1", 150, 200, "+"), ("chr1", 300, 380, "+")],
                cds_bp=130,
                parent_bp=200,
            )
        ]
    )
    row = _row(predict(iso, prop), "t1")
    assert row["nmd_status"] == NMDStatus.ESCAPED.value
    assert row["nmd_confidence"] == ORFConfidence.MEDIUM.value


def test_propagated_intact_yields_high_confidence() -> None:
    iso = _iso_pr(
        [
            ("chr1", 100, 200, "+", "t1"),
            ("chr1", 300, 400, "+", "t1"),
        ]
    )
    prop = pd.DataFrame(
        [
            _prop_row(
                "t1",
                ORFOutcome.PROPAGATED_INTACT,
                intervals=[("chr1", 150, 200, "+"), ("chr1", 300, 380, "+")],
                cds_bp=130,
                parent_bp=130,
            )
        ]
    )
    row = _row(predict(iso, prop), "t1")
    assert row["nmd_confidence"] == ORFConfidence.HIGH.value


def test_multiple_isoforms_classified_independently() -> None:
    iso = _iso_pr(
        [
            ("chr1", 100, 200, "+", "t1"),
            ("chr1", 300, 400, "+", "t1"),
            ("chr1", 100, 300, "+", "t2"),
            ("chr1", 400, 500, "+", "t2"),
            ("chr1", 600, 700, "+", "t2"),
        ]
    )
    prop = pd.DataFrame(
        [
            _prop_row(
                "t1",
                ORFOutcome.PROPAGATED_INTACT,
                intervals=[("chr1", 150, 200, "+"), ("chr1", 300, 380, "+")],
                cds_bp=130,
                parent_bp=130,
            ),
            _prop_row(
                "t2",
                ORFOutcome.PROPAGATED_INTACT,
                intervals=[("chr1", 100, 250, "+")],
                cds_bp=150,
                parent_bp=150,
            ),
        ]
    )
    result = predict(iso, prop)
    assert _row(result, "t1")["nmd_status"] == NMDStatus.ESCAPED.value
    assert _row(result, "t2")["nmd_status"] == NMDStatus.SENSITIVE.value
