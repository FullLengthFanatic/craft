"""Tests for craft.core.pfam."""

from pathlib import Path

import pandas as pd
import pyhmmer
import pyranges as pr
import pysam
import pytest

from craft.core.pfam import scan, translate

# 26 aa protein and its in-frame DNA encoding (one canonical codon per aa).
# M  V  K  L  L  A  A  A  K  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W
_PROTEIN = "MVKLLAAAKDEFGHIKLMNPQRSTVW"
_DNA = (
    "ATG" "GTT" "AAA" "CTT" "CTT" "GCT" "GCT" "GCT" "AAA" "GAT"
    "GAA" "TTT" "GGT" "CAT" "ATT" "AAA" "CTT" "ATG" "AAT" "CCT"
    "CAA" "CGT" "TCT" "ACT" "GTT" "TGG"
)
assert len(_DNA) == 78


def _rc(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _write_test_hmm(path: Path, name: str = "TEST_DOMAIN", sequence: str = _PROTEIN) -> None:
    """Build a tiny single-sequence HMM and write it to ``path``."""
    alphabet = pyhmmer.easel.Alphabet.amino()
    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    seq = pyhmmer.easel.TextSequence(
        name=name.encode(), sequence=sequence
    ).digitize(alphabet)
    hmm, _, _ = builder.build(seq, background)
    # pyhmmer writes a 'COM   [1]' line with no command; HMMer's parser
    # rejects it. Setting command_line gives the line a body.
    hmm.command_line = "hmmbuild test"
    with open(path, "wb") as fh:
        hmm.write(fh)


@pytest.fixture
def test_hmm(tmp_path: Path) -> Path:
    p = tmp_path / "test.hmm"
    _write_test_hmm(p)
    return p


@pytest.fixture
def genome_with_orf(tmp_path: Path) -> Path:
    """FASTA with a 78-bp CDS at chr1:100-178 (matches _PROTEIN) and the
    reverse-complement of that at chr2:100-178 (for minus-strand tests)."""
    # chr1 (+ strand): N(100) + DNA(78) + N(22) = 200
    chr1 = ("N" * 100) + _DNA + ("N" * 22)
    assert len(chr1) == 200
    # chr2 (- strand): the forward genomic sequence is _rc(_DNA) so that
    # reading on - strand from 100-178 gives the original DNA -> protein.
    chr2 = ("N" * 100) + _rc(_DNA) + ("N" * 22)
    assert len(chr2) == 200
    fasta = tmp_path / "genome.fa"
    fasta.write_text(f">chr1\n{chr1}\n>chr2\n{chr2}\n")
    pysam.faidx(str(fasta))
    return fasta


def _reference(rows: list[tuple]) -> pr.PyRanges:
    cols = ["Chromosome", "Start", "End", "Strand", "transcript_id", "Feature"]
    df = pd.DataFrame(rows, columns=cols)
    return pr.PyRanges(df)


def _per_isoform(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows)


def test_translate_basic() -> None:
    assert translate(_DNA) == _PROTEIN


def test_translate_stops_at_stop_codon() -> None:
    # ATG GCT TAA -> M A stop
    assert translate("ATGGCTTAA") == "MA"


def test_translate_unknown_codon_becomes_x() -> None:
    # ATG NGN ATG -> M X M
    assert translate("ATGNGNATG") == "MXM"


def test_translate_empty() -> None:
    assert translate("") == ""


def test_scan_empty_per_isoform_returns_empty(
    test_hmm: Path, genome_with_orf: Path
) -> None:
    ref = _reference(
        [("chr1", 100, 178, "+", "t_ref", "CDS")]
    )
    result = scan(
        pd.DataFrame(
            columns=[
                "transcript_id",
                "parent_tx_id",
                "propagated_cds_intervals",
                "denovo_cds_intervals",
            ]
        ),
        ref,
        test_hmm,
        genome_with_orf,
    )
    assert result.empty
    assert list(result.columns) == [
        "transcript_id",
        "iso_pfam_domains",
        "parent_pfam_domains",
        "pfam_preserved",
        "pfam_lost",
        "pfam_gained",
    ]


def test_scan_full_match_yields_preserved_domain(
    test_hmm: Path, genome_with_orf: Path
) -> None:
    ref = _reference(
        [
            ("chr1", 100, 200, "+", "t_ref", "exon"),
            ("chr1", 100, 178, "+", "t_ref", "CDS"),
        ]
    )
    df = _per_isoform(
        [
            {
                "transcript_id": "t_intact",
                "parent_tx_id": "t_ref",
                "propagated_cds_intervals": [("chr1", 100, 178, "+")],
                "denovo_cds_intervals": [],
            }
        ]
    )
    result = scan(df, ref, test_hmm, genome_with_orf)
    row = result.iloc[0]
    assert row["iso_pfam_domains"] == ["TEST_DOMAIN"]
    assert row["parent_pfam_domains"] == ["TEST_DOMAIN"]
    assert row["pfam_preserved"] == ["TEST_DOMAIN"]
    assert row["pfam_lost"] == []
    assert row["pfam_gained"] == []


def test_scan_iso_truncation_yields_lost_domain(
    test_hmm: Path, genome_with_orf: Path
) -> None:
    ref = _reference(
        [
            ("chr1", 100, 200, "+", "t_ref", "exon"),
            ("chr1", 100, 178, "+", "t_ref", "CDS"),
        ]
    )
    # iso CDS is only 30 bp = 10 aa -> too short for HMM hit.
    df = _per_isoform(
        [
            {
                "transcript_id": "t_truncated",
                "parent_tx_id": "t_ref",
                "propagated_cds_intervals": [("chr1", 100, 130, "+")],
                "denovo_cds_intervals": [],
            }
        ]
    )
    result = scan(df, ref, test_hmm, genome_with_orf)
    row = result.iloc[0]
    assert row["iso_pfam_domains"] == []
    assert row["parent_pfam_domains"] == ["TEST_DOMAIN"]
    assert row["pfam_lost"] == ["TEST_DOMAIN"]
    assert row["pfam_preserved"] == []
    assert row["pfam_gained"] == []


def test_scan_no_parent_with_denovo_yields_gained_domain(
    test_hmm: Path, genome_with_orf: Path
) -> None:
    ref = _reference([])
    df = _per_isoform(
        [
            {
                "transcript_id": "t_novel",
                "parent_tx_id": "",
                "propagated_cds_intervals": [],
                "denovo_cds_intervals": [("chr1", 100, 178, "+")],
            }
        ]
    )
    result = scan(df, ref, test_hmm, genome_with_orf)
    row = result.iloc[0]
    assert row["iso_pfam_domains"] == ["TEST_DOMAIN"]
    assert row["parent_pfam_domains"] == []
    assert row["pfam_gained"] == ["TEST_DOMAIN"]
    assert row["pfam_preserved"] == []
    assert row["pfam_lost"] == []


def test_scan_minus_strand_translates_via_reverse_complement(
    test_hmm: Path, genome_with_orf: Path
) -> None:
    ref = _reference(
        [
            ("chr2", 100, 200, "-", "t_ref_minus", "exon"),
            ("chr2", 100, 178, "-", "t_ref_minus", "CDS"),
        ]
    )
    df = _per_isoform(
        [
            {
                "transcript_id": "t_intact_minus",
                "parent_tx_id": "t_ref_minus",
                "propagated_cds_intervals": [("chr2", 100, 178, "-")],
                "denovo_cds_intervals": [],
            }
        ]
    )
    result = scan(df, ref, test_hmm, genome_with_orf)
    row = result.iloc[0]
    assert row["iso_pfam_domains"] == ["TEST_DOMAIN"]
    assert row["parent_pfam_domains"] == ["TEST_DOMAIN"]
    assert row["pfam_preserved"] == ["TEST_DOMAIN"]


def test_scan_caches_repeated_proteins(
    test_hmm: Path, genome_with_orf: Path
) -> None:
    """Same protein on multiple isoforms should still classify consistently."""
    ref = _reference(
        [
            ("chr1", 100, 200, "+", "t_ref", "exon"),
            ("chr1", 100, 178, "+", "t_ref", "CDS"),
        ]
    )
    df = _per_isoform(
        [
            {
                "transcript_id": "t1",
                "parent_tx_id": "t_ref",
                "propagated_cds_intervals": [("chr1", 100, 178, "+")],
                "denovo_cds_intervals": [],
            },
            {
                "transcript_id": "t2",
                "parent_tx_id": "t_ref",
                "propagated_cds_intervals": [("chr1", 100, 178, "+")],
                "denovo_cds_intervals": [],
            },
        ]
    )
    result = scan(df, ref, test_hmm, genome_with_orf)
    assert len(result) == 2
    for _, row in result.iterrows():
        assert row["pfam_preserved"] == ["TEST_DOMAIN"]
