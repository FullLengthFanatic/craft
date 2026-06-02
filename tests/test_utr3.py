"""Tests for craft.core.utr3."""

from pathlib import Path

import pandas as pd
import pyranges as pr
import pysam
import pytest

from craft.core.orf.propagation import ORFOutcome
from craft.core.utr3 import annotate, polya_near_3prime_end, polya_signal


def _iso_pr(records: list[tuple]) -> pr.PyRanges:
    cols = ["Chromosome", "Start", "End", "Strand", "transcript_id"]
    df = pd.DataFrame(records, columns=cols)
    return pr.PyRanges(df)


def _reference(exon_records: list[tuple], cds_records: list[tuple]) -> pr.PyRanges:
    cols = ["Chromosome", "Start", "End", "Strand", "transcript_id", "Feature"]
    rows = [(*r, "exon") for r in exon_records] + [(*r, "CDS") for r in cds_records]
    df = pd.DataFrame(rows, columns=cols)
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


@pytest.fixture
def synthetic_genome(tmp_path: Path) -> Path:
    # chr1: 60 bp; positions 0-29 all N; positions 30-59 = UTR for + strand test
    #   "GCGC" (4) + "AATAAA" (6) + "GCGCGCGCGCGCGCGCGCGC" (20) = 30
    seq1 = "N" * 30 + "GCGC" + "AATAAA" + "GCGCGCGCGCGCGCGCGCGC"
    # chr2: 60 bp; positions 0-29 forward-strand sequence, positions 30-59 all N
    #   forward 0-29 = "GCGCGCGCGCGCGCGCGCGCTTTATTGCGC"
    #   revcomp of that = "GCGCAATAAAGCGCGCGCGCGCGCGCGCGC"
    #   So on - strand iso, the UTR (genomic 0-29) reverse-complements to a sequence
    #   containing AATAAA in transcript order.
    seq2 = "GCGCGCGCGCGCGCGCGCGCTTTATTGCGC" + "N" * 30
    fasta_path = tmp_path / "genome.fa"
    fasta_path.write_text(f">chr1\n{seq1}\n>chr2\n{seq2}\n")
    pysam.faidx(str(fasta_path))
    return fasta_path


def test_polya_signal_finds_canonical_aataaa() -> None:
    sig = polya_signal("GCGCAATAAAGC")
    assert sig["motif"] == "AATAAA"
    assert sig["distance_from_3p_end"] == 2


def test_polya_signal_finds_variant_when_canonical_absent() -> None:
    sig = polya_signal("GCGCATTAAAGC")
    assert sig["motif"] == "ATTAAA"
    assert sig["distance_from_3p_end"] == 2


def test_polya_signal_prefers_canonical_when_both_present() -> None:
    sig = polya_signal("ATTAAAGCGCAATAAAGCGC")
    assert sig["motif"] == "AATAAA"


def test_polya_signal_picks_rightmost_occurrence_for_same_motif() -> None:
    sig = polya_signal("AATAAAGCGCAATAAAGC")
    assert sig["motif"] == "AATAAA"
    assert sig["distance_from_3p_end"] == 2


def test_polya_signal_returns_empty_when_no_motif() -> None:
    sig = polya_signal("GCGCGCGCGCGCGC")
    assert sig["motif"] == ""
    assert sig["distance_from_3p_end"] == -1


def test_polya_signal_is_case_insensitive() -> None:
    sig = polya_signal("gcgcaataaagc")
    assert sig["motif"] == "AATAAA"


def test_identical_iso_and_parent_yield_zero_utr_delta() -> None:
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
    iso = _iso_pr(
        [
            ("chr1", 100, 200, "+", "t1"),
            ("chr1", 300, 400, "+", "t1"),
            ("chr1", 500, 600, "+", "t1"),
        ]
    )
    prop = pd.DataFrame(
        [
            _prop_row(
                "t1",
                ORFOutcome.PROPAGATED_INTACT,
                intervals=[
                    ("chr1", 150, 200, "+"),
                    ("chr1", 300, 400, "+"),
                    ("chr1", 500, 550, "+"),
                ],
                cds_bp=200,
                parent_bp=200,
            )
        ]
    )
    result = annotate(iso, prop, ref)
    row = result.iloc[0]
    assert row["transcript_id"] == "t1"
    assert int(row["iso_utr3_length_nt"]) == 50
    assert int(row["parent_utr3_length_nt"]) == 50
    assert int(row["utr3_length_delta_nt"]) == 0


def test_iso_3p_truncated_has_shorter_utr_than_parent() -> None:
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
    iso = _iso_pr(
        [
            ("chr1", 100, 200, "+", "t1"),
            ("chr1", 300, 400, "+", "t1"),
            ("chr1", 500, 580, "+", "t1"),
        ]
    )
    prop = pd.DataFrame(
        [
            _prop_row(
                "t1",
                ORFOutcome.PROPAGATED_INTACT,
                intervals=[
                    ("chr1", 150, 200, "+"),
                    ("chr1", 300, 400, "+"),
                    ("chr1", 500, 550, "+"),
                ],
                cds_bp=200,
                parent_bp=200,
            )
        ]
    )
    result = annotate(iso, prop, ref)
    row = result.iloc[0]
    assert int(row["iso_utr3_length_nt"]) == 30
    assert int(row["parent_utr3_length_nt"]) == 50
    assert int(row["utr3_length_delta_nt"]) == -20


def test_minus_strand_utr_length() -> None:
    ref = _reference(
        exon_records=[
            ("chr1", 50, 200, "-", "t_ref"),
            ("chr1", 300, 400, "-", "t_ref"),
            ("chr1", 500, 600, "-", "t_ref"),
        ],
        cds_records=[
            ("chr1", 100, 200, "-", "t_ref"),
            ("chr1", 300, 400, "-", "t_ref"),
            ("chr1", 500, 550, "-", "t_ref"),
        ],
    )
    iso = _iso_pr(
        [
            ("chr1", 50, 200, "-", "t1"),
            ("chr1", 300, 400, "-", "t1"),
            ("chr1", 500, 600, "-", "t1"),
        ]
    )
    prop = pd.DataFrame(
        [
            _prop_row(
                "t1",
                ORFOutcome.PROPAGATED_INTACT,
                intervals=[
                    ("chr1", 100, 200, "-"),
                    ("chr1", 300, 400, "-"),
                    ("chr1", 500, 550, "-"),
                ],
                cds_bp=200,
                parent_bp=200,
            )
        ]
    )
    result = annotate(iso, prop, ref)
    row = result.iloc[0]
    # - strand stop = min(Start of CDS) = 100. UTR = exonic positions < 100.
    # (50, 200) contributes: min(200, 100) - 50 = 50. Others: 0.
    assert int(row["iso_utr3_length_nt"]) == 50
    assert int(row["parent_utr3_length_nt"]) == 50


def test_no_parent_yields_none_metrics() -> None:
    ref = _reference([], [])
    iso = _iso_pr([("chr2", 100, 200, "+", "t_novel")])
    prop = pd.DataFrame(
        [
            _prop_row(
                "t_novel",
                ORFOutcome.NO_PARENT,
                intervals=[],
                parent_tx_id="",
                start_covered=False,
                stop_covered=False,
            )
        ]
    )
    result = annotate(iso, prop, ref)
    row = result.iloc[0]
    assert row["transcript_id"] == "t_novel"
    assert pd.isna(row["iso_utr3_length_nt"])
    assert pd.isna(row["utr3_length_delta_nt"])
    assert row["polya_signal_motif"] == ""


def test_stop_not_observed_yields_none_metrics() -> None:
    ref = _reference(
        exon_records=[
            ("chr1", 100, 200, "+", "t_ref"),
            ("chr1", 300, 400, "+", "t_ref"),
        ],
        cds_records=[
            ("chr1", 150, 200, "+", "t_ref"),
            ("chr1", 300, 400, "+", "t_ref"),
        ],
    )
    iso = _iso_pr(
        [
            ("chr1", 100, 200, "+", "t1"),
            ("chr1", 300, 350, "+", "t1"),
        ]
    )
    prop = pd.DataFrame(
        [
            _prop_row(
                "t1",
                ORFOutcome.STOP_NOT_OBSERVED,
                intervals=[("chr1", 150, 200, "+"), ("chr1", 300, 350, "+")],
                cds_bp=100,
                parent_bp=150,
                stop_covered=False,
            )
        ]
    )
    result = annotate(iso, prop, ref)
    row = result.iloc[0]
    assert pd.isna(row["iso_utr3_length_nt"])


def test_polya_signal_extracted_from_plus_strand_utr(synthetic_genome: Path) -> None:
    ref = _reference(
        exon_records=[("chr1", 0, 60, "+", "t_ref")],
        cds_records=[("chr1", 0, 30, "+", "t_ref")],
    )
    iso = _iso_pr([("chr1", 0, 60, "+", "t1")])
    prop = pd.DataFrame(
        [
            _prop_row(
                "t1",
                ORFOutcome.PROPAGATED_INTACT,
                intervals=[("chr1", 0, 30, "+")],
                cds_bp=30,
                parent_bp=30,
            )
        ]
    )
    result = annotate(iso, prop, ref, genome_fasta=synthetic_genome)
    row = result.iloc[0]
    assert int(row["iso_utr3_length_nt"]) == 30
    assert row["polya_signal_motif"] == "AATAAA"
    # UTR = "GCGCAATAAAGCGCGCGCGCGCGCGCGCGC" (length 30). AATAAA at idx 4, ends at 10.
    # Distance from 3' end = 30 - 10 = 20.
    assert int(row["polya_signal_distance_nt"]) == 20


def test_polya_signal_extracted_from_minus_strand_utr(synthetic_genome: Path) -> None:
    ref = _reference(
        exon_records=[("chr2", 0, 60, "-", "t_ref")],
        cds_records=[("chr2", 30, 60, "-", "t_ref")],
    )
    iso = _iso_pr([("chr2", 0, 60, "-", "t1")])
    prop = pd.DataFrame(
        [
            _prop_row(
                "t1",
                ORFOutcome.PROPAGATED_INTACT,
                intervals=[("chr2", 30, 60, "-")],
                cds_bp=30,
                parent_bp=30,
            )
        ]
    )
    result = annotate(iso, prop, ref, genome_fasta=synthetic_genome)
    row = result.iloc[0]
    assert int(row["iso_utr3_length_nt"]) == 30
    assert row["polya_signal_motif"] == "AATAAA"


def test_polya_near_3prime_end_plus_strand_with_signal(
    synthetic_genome: Path,
) -> None:
    # chr1 positions 30-59 contain "GCGCAATAAAGCGCGCGCGCGCGCGCGCGC".
    # + strand iso (0, 60) has 3' end at 60; last 50 bp window is positions 10-60,
    # which contains AATAAA.
    exons = pd.DataFrame(
        [{"Chromosome": "chr1", "Start": 0, "End": 60, "Strand": "+", "transcript_id": "t1"}]
    )
    result = polya_near_3prime_end(exons, "+", synthetic_genome)
    assert result["found"] is True
    assert result["motif"] == "AATAAA"


def test_polya_near_3prime_end_minus_strand_with_signal(
    synthetic_genome: Path,
) -> None:
    # chr2 positions 0-29 forward = "GCGCGCGCGCGCGCGCGCGCTTTATTGCGC",
    # revcomp contains AATAAA. - strand iso (0, 60) sees this in transcript orientation.
    exons = pd.DataFrame(
        [{"Chromosome": "chr2", "Start": 0, "End": 60, "Strand": "-", "transcript_id": "t1"}]
    )
    result = polya_near_3prime_end(exons, "-", synthetic_genome)
    assert result["found"] is True
    assert result["motif"] == "AATAAA"


def test_polya_near_3prime_end_no_signal(synthetic_genome: Path) -> None:
    # chr1 positions 0-29 are all N (no motif). - strand iso ending there has no signal.
    exons = pd.DataFrame(
        [{"Chromosome": "chr1", "Start": 0, "End": 30, "Strand": "-", "transcript_id": "t1"}]
    )
    result = polya_near_3prime_end(exons, "-", synthetic_genome)
    assert result["found"] is False
    assert result["motif"] == ""


def test_polya_near_3prime_end_accepts_path_argument(synthetic_genome: Path) -> None:
    """Passing a Path works (function opens FASTA internally and closes it)."""
    exons = pd.DataFrame(
        [{"Chromosome": "chr1", "Start": 0, "End": 60, "Strand": "+", "transcript_id": "t1"}]
    )
    # Pass the Path object directly, not an open pysam.FastaFile handle.
    result = polya_near_3prime_end(exons, "+", synthetic_genome)
    assert result["found"] is True


def test_polya_signal_motif_empty_when_utr_lacks_signal(synthetic_genome: Path) -> None:
    # On chr1, position 0-29 is all Ns. Use - strand iso with CDS at 30-60.
    ref = _reference(
        exon_records=[("chr1", 0, 60, "-", "t_ref")],
        cds_records=[("chr1", 30, 60, "-", "t_ref")],
    )
    iso = _iso_pr([("chr1", 0, 60, "-", "t1")])
    prop = pd.DataFrame(
        [
            _prop_row(
                "t1",
                ORFOutcome.PROPAGATED_INTACT,
                intervals=[("chr1", 30, 60, "-")],
                cds_bp=30,
                parent_bp=30,
            )
        ]
    )
    result = annotate(iso, prop, ref, genome_fasta=synthetic_genome)
    row = result.iloc[0]
    assert row["polya_signal_motif"] == ""


# ---- resolved 3'UTR + 5'UTR (annotate_resolved) ------------------------------

from craft.core.utr3 import _utr5_length, annotate_resolved  # noqa: E402


def _classified(iso_records: list[tuple], parent_tx: str, tx_id: str = "t1") -> pr.PyRanges:
    cols = ["Chromosome", "Start", "End", "Strand", "transcript_id", "parent_tx_id"]
    rows = [(c, s, e, st, tx_id, parent_tx) for (c, s, e, st) in iso_records]
    return pr.PyRanges(pd.DataFrame(rows, columns=cols))


def _res_df(tx_id, status, intervals, stop=True):
    return pd.DataFrame(
        [
            {
                "transcript_id": tx_id,
                "resolved_orf_status": status,
                "resolved_stop_pos": None,
                "resolved_cds_bp": 100,
                "resolved_aa_length": 33,
                "resolved_cds_intervals": intervals,
                "ptc_introduced": False,
                "intron_retained_in_cds": False,
                "frame_consistent": True,
                "stop_in_transcript": stop,
                "uorf_count": 0,
                "uorf_triggers_nmd": False,
            }
        ]
    )


def test_utr5_length_plus_and_minus_strand() -> None:
    exons = pd.DataFrame(
        {"Start": [100, 300], "End": [200, 400]},
    )
    # + strand: bases upstream of start 150 -> exon1 [100,150) = 50.
    assert _utr5_length(exons, 150, "+") == 50
    # - strand: 5' is high coord; start 350 -> bases strictly above it (351..399) = 49.
    assert _utr5_length(exons, 350, "-") == 49


def test_annotate_resolved_full_length_has_zero_deltas() -> None:
    iso = [("chr1", 100, 200, "+"), ("chr1", 300, 400, "+")]
    classified = _classified(iso, "t_ref")
    ref = _reference(
        exon_records=[("chr1", 100, 200, "+", "t_ref"), ("chr1", 300, 400, "+", "t_ref")],
        cds_records=[("chr1", 150, 200, "+", "t_ref"), ("chr1", 300, 350, "+", "t_ref")],
    )
    res = _res_df("t1", "intact", [("chr1", 150, 200, "+"), ("chr1", 300, 350, "+")])
    out = annotate_resolved(classified, res, ref)
    row = out[out["transcript_id"] == "t1"].iloc[0]
    assert row["iso_utr5_length_nt"] == 50
    assert row["parent_utr5_length_nt"] == 50
    assert row["utr5_length_delta_nt"] == 0
    assert row["iso_utr3_length_resolved_nt"] == 50  # exon2 [350,400)
    assert not row["long_utr3_triggers_nmd"]


def test_annotate_resolved_long_utr3_triggers_flag() -> None:
    iso = [("chr1", 100, 200, "+"), ("chr1", 300, 2000, "+")]
    classified = _classified(iso, "t_ref")
    ref = _reference(
        exon_records=[("chr1", 100, 200, "+", "t_ref"), ("chr1", 300, 2000, "+", "t_ref")],
        cds_records=[("chr1", 150, 200, "+", "t_ref"), ("chr1", 300, 350, "+", "t_ref")],
    )
    res = _res_df("t1", "intact", [("chr1", 150, 200, "+"), ("chr1", 300, 350, "+")])
    out = annotate_resolved(classified, res, ref, long_utr3_nt=1000)
    row = out[out["transcript_id"] == "t1"].iloc[0]
    assert row["iso_utr3_length_resolved_nt"] == 1650
    assert row["long_utr3_triggers_nmd"]
