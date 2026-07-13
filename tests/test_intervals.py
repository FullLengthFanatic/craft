"""Tests for craft.core.intervals."""

import pandas as pd
import pyranges as pr

from craft.core.intervals import (
    genomic_position_at_transcript_coordinate,
    splice_junctions,
    transcript_coordinate,
)


def _exons(records: list[tuple]) -> pr.PyRanges:
    cols = ["Chromosome", "Start", "End", "Strand", "transcript_id"]
    df = pd.DataFrame(records, columns=cols)
    return pr.PyRanges(df)


def test_three_exons_yield_two_junctions() -> None:
    exons = _exons(
        [
            ("chr1", 100, 150, "+", "t1"),
            ("chr1", 200, 280, "+", "t1"),
            ("chr1", 350, 420, "+", "t1"),
        ]
    )
    j = splice_junctions(exons).df.sort_values("Start").reset_index(drop=True)
    assert len(j) == 2
    assert list(j["Start"]) == [150, 280]
    assert list(j["End"]) == [200, 350]
    assert (j["transcript_id"] == "t1").all()


def test_single_exon_transcript_yields_no_junctions() -> None:
    exons = _exons([("chr1", 100, 200, "+", "t1")])
    assert len(splice_junctions(exons)) == 0


def test_empty_input_yields_empty_output() -> None:
    empty = _exons([])
    assert len(splice_junctions(empty)) == 0


def test_multiple_transcripts_grouped_independently() -> None:
    exons = _exons(
        [
            ("chr1", 100, 150, "+", "t1"),
            ("chr1", 200, 280, "+", "t1"),
            ("chr2", 500, 580, "-", "t2"),
            ("chr2", 700, 800, "-", "t2"),
            ("chr2", 900, 1000, "-", "t2"),
        ]
    )
    j = splice_junctions(exons).df
    t1 = j[j["transcript_id"] == "t1"].sort_values("Start").reset_index(drop=True)
    t2 = j[j["transcript_id"] == "t2"].sort_values("Start").reset_index(drop=True)
    assert len(t1) == 1
    assert (t1["Start"].iloc[0], t1["End"].iloc[0]) == (150, 200)
    assert len(t2) == 2
    assert list(t2["Start"]) == [580, 800]
    assert list(t2["End"]) == [700, 900]
    assert list(t2["Strand"]) == ["-", "-"]


def test_unsorted_input_still_produces_correct_junctions() -> None:
    exons = _exons(
        [
            ("chr1", 350, 420, "+", "t1"),
            ("chr1", 100, 150, "+", "t1"),
            ("chr1", 200, 280, "+", "t1"),
        ]
    )
    j = splice_junctions(exons).df.sort_values("Start").reset_index(drop=True)
    assert list(j["Start"]) == [150, 280]
    assert list(j["End"]) == [200, 350]


def test_junction_index_in_genomic_order() -> None:
    exons = _exons(
        [
            ("chr1", 100, 150, "+", "t1"),
            ("chr1", 200, 280, "+", "t1"),
            ("chr1", 350, 420, "+", "t1"),
        ]
    )
    j = splice_junctions(exons).df.sort_values("Start").reset_index(drop=True)
    assert list(j["junction_index"]) == [0, 1]


def test_strand_preserved_on_minus() -> None:
    exons = _exons(
        [
            ("chr1", 100, 150, "-", "t1"),
            ("chr1", 200, 280, "-", "t1"),
        ]
    )
    j = splice_junctions(exons).df
    assert j["Strand"].iloc[0] == "-"


def test_transcript_coordinate_round_trip_across_junctions() -> None:
    exons = pd.DataFrame({"Start": [100, 300], "End": [105, 310]})
    cases = (
        ("+", [(0, 100), (4, 104), (5, 300), (14, 309)]),
        ("-", [(0, 309), (9, 300), (10, 104), (14, 100)]),
    )
    for strand, positions in cases:
        for expected, genomic in positions:
            assert transcript_coordinate(exons, genomic, strand) == expected
            assert genomic_position_at_transcript_coordinate(exons, expected, strand) == genomic
