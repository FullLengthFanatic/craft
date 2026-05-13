"""Tests for craft.core.orf.denovo."""

from pathlib import Path

import pandas as pd
import pyranges as pr
import pysam
import pytest

from craft.core.orf.denovo import predict

_RC_TABLE = str.maketrans("ACGTN", "TGCAN")


def _rc(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]


def _iso_pr(records: list[tuple]) -> pr.PyRanges:
    cols = ["Chromosome", "Start", "End", "Strand", "transcript_id"]
    df = pd.DataFrame(records, columns=cols)
    return pr.PyRanges(df)


def _row(result: pd.DataFrame, tx: str) -> pd.Series:
    return result[result["transcript_id"] == tx].iloc[0]


@pytest.fixture
def synthetic_genome(tmp_path: Path) -> Path:
    """FASTA with two contigs designed for + and - strand ORF tests, plus a
    multi-exon contig for genomic-mapping checks."""

    # chr1 (200 bp): ATG at 0, then 59 GCC codons, TAA at 180, then filler.
    # An iso single-exon (0, 200) on + strand gets ORF at transcript [0, 180).
    fwd_chr1 = "ATG" + ("GCC" * 59) + "TAA" + ("G" * 17)
    assert len(fwd_chr1) == 200

    # chr2 (200 bp): same transcript content as chr1, but the iso is on - strand,
    # so the FORWARD genomic is the reverse-complement of that transcript.
    fwd_chr2 = _rc(fwd_chr1)
    assert len(fwd_chr2) == 200

    # chr3 (600 bp): three "exon" slots at [100,200), [300,400), [500,600).
    # Transcript = concatenation of those slots in genomic order = 300 bp.
    # In the transcript: 50 G's, ATG, 59 GCC codons, TAA, 67 G's.
    tx_chr3 = ("G" * 50) + "ATG" + ("GCC" * 59) + "TAA" + ("G" * 67)
    assert len(tx_chr3) == 300
    fwd_chr3 = (
        ("N" * 100) + tx_chr3[0:100]
        + ("N" * 100) + tx_chr3[100:200]
        + ("N" * 100) + tx_chr3[200:300]
    )
    assert len(fwd_chr3) == 600

    # chr4 (200 bp): sequence with no ATG/stop pattern -- pure filler.
    fwd_chr4 = "GCG" * 66 + "GC"
    assert len(fwd_chr4) == 200

    # chr5 (300 bp): two ORFs of different lengths; predict should pick the longest.
    # Transcript: ATG + 9 codons (27 bp) + TAA + filler + ATG + 59 codons + TAA + filler.
    short_orf = "ATG" + ("GCC" * 9) + "TAA"  # len 3 + 27 + 3 = 33
    long_orf = "ATG" + ("GCC" * 59) + "TAA"  # len 3 + 177 + 3 = 183
    fwd_chr5 = short_orf + ("G" * 24) + long_orf + ("G" * 60)
    assert len(fwd_chr5) == 300

    fasta_path = tmp_path / "genome.fa"
    fasta_path.write_text(
        f">chr1\n{fwd_chr1}\n"
        f">chr2\n{fwd_chr2}\n"
        f">chr3\n{fwd_chr3}\n"
        f">chr4\n{fwd_chr4}\n"
        f">chr5\n{fwd_chr5}\n"
    )
    pysam.faidx(str(fasta_path))
    return fasta_path


def test_single_exon_plus_strand_finds_orf(synthetic_genome: Path) -> None:
    iso = _iso_pr([("chr1", 0, 200, "+", "t1")])
    result = predict(iso, synthetic_genome)
    row = _row(result, "t1")
    assert bool(row["denovo_orf_found"]) is True
    assert int(row["denovo_cds_bp"]) == 180
    assert int(row["denovo_orf_aa_length"]) == 60
    assert row["denovo_start_codon"] == "ATG"
    assert row["denovo_stop_codon"] == "TAA"
    assert row["denovo_cds_intervals"] == [("chr1", 0, 180, "+")]


def test_single_exon_minus_strand_finds_orf(synthetic_genome: Path) -> None:
    iso = _iso_pr([("chr2", 0, 200, "-", "t1")])
    result = predict(iso, synthetic_genome)
    row = _row(result, "t1")
    assert bool(row["denovo_orf_found"]) is True
    assert int(row["denovo_cds_bp"]) == 180
    # transcript [0, 180) maps to genomic [200-180, 200-0) = [20, 200) on chr2
    assert row["denovo_cds_intervals"] == [("chr2", 20, 200, "-")]


def test_multi_exon_plus_strand_orf_spans_junctions(synthetic_genome: Path) -> None:
    iso = _iso_pr(
        [
            ("chr3", 100, 200, "+", "t1"),
            ("chr3", 300, 400, "+", "t1"),
            ("chr3", 500, 600, "+", "t1"),
        ]
    )
    result = predict(iso, synthetic_genome)
    row = _row(result, "t1")
    assert bool(row["denovo_orf_found"]) is True
    # transcript ORF at [50, 230). Mapping to genomic:
    # exon1 overlap [50, 100) -> [100+50, 100+100) = [150, 200)
    # exon2 overlap [100, 200) -> [300+0, 300+100) = [300, 400)
    # exon3 overlap [200, 230) -> [500+0, 500+30) = [500, 530)
    assert row["denovo_cds_intervals"] == [
        ("chr3", 150, 200, "+"),
        ("chr3", 300, 400, "+"),
        ("chr3", 500, 530, "+"),
    ]
    assert int(row["denovo_cds_bp"]) == 180
    assert int(row["denovo_orf_aa_length"]) == 60


def test_no_orf_returns_empty(synthetic_genome: Path) -> None:
    iso = _iso_pr([("chr4", 0, 200, "+", "t1")])
    result = predict(iso, synthetic_genome)
    row = _row(result, "t1")
    assert bool(row["denovo_orf_found"]) is False
    assert int(row["denovo_cds_bp"]) == 0
    assert row["denovo_cds_intervals"] == []
    assert row["denovo_start_codon"] == ""


def test_below_min_length_returns_empty(synthetic_genome: Path) -> None:
    # chr1 has an ORF of 60 aa. With min_orf_aa=100 it should be filtered out.
    iso = _iso_pr([("chr1", 0, 200, "+", "t1")])
    result = predict(iso, synthetic_genome, min_orf_aa=100)
    row = _row(result, "t1")
    assert bool(row["denovo_orf_found"]) is False


def test_picks_longest_orf_when_multiple_present(synthetic_genome: Path) -> None:
    # chr5 has a 10 aa ORF followed by a 60 aa ORF (default min 50 keeps only the long one).
    iso = _iso_pr([("chr5", 0, 300, "+", "t1")])
    result = predict(iso, synthetic_genome)
    row = _row(result, "t1")
    assert int(row["denovo_orf_aa_length"]) == 60


def test_picks_longest_orf_even_when_short_orf_passes_min(synthetic_genome: Path) -> None:
    # Lower min to 5 aa so both ORFs qualify; ensure the 60 aa one still wins.
    iso = _iso_pr([("chr5", 0, 300, "+", "t1")])
    result = predict(iso, synthetic_genome, min_orf_aa=5)
    row = _row(result, "t1")
    assert int(row["denovo_orf_aa_length"]) == 60


def test_empty_input_returns_empty_frame(synthetic_genome: Path) -> None:
    iso = _iso_pr([])
    result = predict(iso, synthetic_genome)
    assert result.empty
    assert list(result.columns) == [
        "transcript_id",
        "denovo_orf_found",
        "denovo_cds_bp",
        "denovo_cds_intervals",
        "denovo_orf_aa_length",
        "denovo_start_codon",
        "denovo_stop_codon",
    ]


def test_multiple_isoforms_handled_independently(synthetic_genome: Path) -> None:
    iso = _iso_pr(
        [
            ("chr1", 0, 200, "+", "t_with_orf"),
            ("chr4", 0, 200, "+", "t_no_orf"),
        ]
    )
    result = predict(iso, synthetic_genome)
    assert bool(_row(result, "t_with_orf")["denovo_orf_found"]) is True
    assert bool(_row(result, "t_no_orf")["denovo_orf_found"]) is False
