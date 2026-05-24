"""Tests for cbench.truncation."""

import pytest
from cbench.truncation import truncate_exons, truncate_sequence


def _total_length(exons: list[tuple[int, int]]) -> int:
    return sum(e - s for s, e in exons)


def test_truncate_sequence_no_op_when_rate_zero() -> None:
    assert truncate_sequence("ACGTACGT", 0.0, "5prime") == "ACGTACGT"


def test_truncate_sequence_3prime() -> None:
    assert truncate_sequence("ACGTACGT", 0.25, "3prime") == "ACGTAC"


def test_truncate_sequence_5prime() -> None:
    assert truncate_sequence("ACGTACGT", 0.25, "5prime") == "GTACGT"


def test_truncate_sequence_both_splits_loss_evenly() -> None:
    # 8 bp * 0.5 = 4 bp total removed; 2 from each end.
    assert truncate_sequence("ACGTACGT", 0.5, "both") == "GTAC"


def test_truncate_sequence_both_odd_loss_takes_extra_from_3prime() -> None:
    # 9 bp * 0.5 = 4.5 -> int = 4. half=2 from 5', remainder 2 from 3'.
    assert truncate_sequence("ACGTACGTA", 0.5, "both") == "GTACG"


def test_truncate_sequence_full_trim_returns_empty() -> None:
    assert truncate_sequence("ACGTACGT", 1.0, "3prime") == ""


def test_truncate_exons_invalid_orientation_raises() -> None:
    with pytest.raises(ValueError, match="orientation"):
        truncate_exons([(0, 100)], "+", 0.25, "middle")


def test_truncate_exons_invalid_strand_raises() -> None:
    with pytest.raises(ValueError, match="strand"):
        truncate_exons([(0, 100)], ".", 0.25, "5prime")


def test_truncate_exons_no_op_when_rate_zero() -> None:
    exons = [(100, 200), (300, 400)]
    assert truncate_exons(exons, "+", 0.0, "5prime") == exons


def test_truncate_exons_full_trim_returns_empty() -> None:
    assert truncate_exons([(100, 200), (300, 400)], "+", 1.0, "3prime") == []


def test_truncate_exons_single_exon_3prime_plus_strand() -> None:
    # 100 bp exon, trim 25% from 3' (genomic right on + strand).
    result = truncate_exons([(1000, 1100)], "+", 0.25, "3prime")
    assert result == [(1000, 1075)]
    assert _total_length(result) == 75


def test_truncate_exons_single_exon_5prime_plus_strand() -> None:
    # 100 bp exon, trim 25% from 5' (genomic left on + strand).
    result = truncate_exons([(1000, 1100)], "+", 0.25, "5prime")
    assert result == [(1025, 1100)]


def test_truncate_exons_single_exon_3prime_minus_strand_trims_left_in_genomic() -> None:
    # On - strand, the 3' transcript end is the genomic-leftmost base.
    # Trimming 25% from 3' should remove 25 bp from the left genomic edge.
    result = truncate_exons([(1000, 1100)], "-", 0.25, "3prime")
    assert result == [(1025, 1100)]


def test_truncate_exons_single_exon_5prime_minus_strand_trims_right_in_genomic() -> None:
    result = truncate_exons([(1000, 1100)], "-", 0.25, "5prime")
    assert result == [(1000, 1075)]


def test_truncate_exons_multi_exon_3prime_removes_last_exon_when_deep_enough() -> None:
    # Three 100 bp exons = 300 bp total. Trim 40% (120 bp) from 3':
    # last exon (100 bp) fully removed; trim 20 bp from end of middle exon.
    exons = [(1000, 1100), (2000, 2100), (3000, 3100)]
    result = truncate_exons(exons, "+", 0.40, "3prime")
    assert result == [(1000, 1100), (2000, 2080)]
    assert _total_length(result) == 180


def test_truncate_exons_multi_exon_5prime_removes_first_exon_when_deep_enough() -> None:
    exons = [(1000, 1100), (2000, 2100), (3000, 3100)]
    result = truncate_exons(exons, "+", 0.40, "5prime")
    assert result == [(2020, 2100), (3000, 3100)]
    assert _total_length(result) == 180


def test_truncate_exons_both_orientation_total_loss_matches_rate() -> None:
    exons = [(1000, 1100), (2000, 2100), (3000, 3100)]  # 300 bp total
    result = truncate_exons(exons, "+", 0.50, "both")
    # 50% of 300 = 150 bp trimmed total. half = 75 from 5', 75 from 3'.
    # 5' end of exon1 moves +75; 3' end of exon3 moves -75.
    assert result == [(1075, 1100), (2000, 2100), (3000, 3025)]
    assert _total_length(result) == 150


def test_truncate_exons_minus_strand_multi_exon_3prime() -> None:
    # On - strand, transcript order is reversed in genomic coords.
    # The 3' end of the transcript is the genomic-leftmost base of the leftmost exon.
    exons = [(1000, 1100), (2000, 2100), (3000, 3100)]
    result = truncate_exons(exons, "-", 0.40, "3prime")
    # Trim 120 bp from transcript 3' end = remove leftmost 120 bp in genomic coords.
    # That removes exon1 (100 bp) entirely and 20 bp from the left of exon2.
    assert result == [(2020, 2100), (3000, 3100)]
    assert _total_length(result) == 180


def test_truncate_exons_preserves_genomic_sort_order() -> None:
    exons = [(1000, 1100), (2000, 2100), (3000, 3100)]
    for strand in ("+", "-"):
        for orientation in ("5prime", "3prime", "both"):
            result = truncate_exons(exons, strand, 0.30, orientation)
            sorted_starts = sorted(s for s, _ in result)
            assert [s for s, _ in result] == sorted_starts, (
                f"strand={strand} orientation={orientation}: not genome-sorted"
            )


def test_truncate_exons_total_length_matches_truncate_sequence() -> None:
    """The same (rate, orientation) on exons and on a transcript sequence must
    leave the same number of bases."""
    exons = [(1000, 1050), (2000, 2070), (3000, 3080)]  # 50 + 70 + 80 = 200 bp
    seq = "A" * 200
    for orientation in ("5prime", "3prime", "both"):
        for rate in (0.05, 0.10, 0.25, 0.50):
            tx_seq = truncate_sequence(seq, rate, orientation)
            for strand in ("+", "-"):
                truncated_exons = truncate_exons(exons, strand, rate, orientation)
                assert _total_length(truncated_exons) == len(tx_seq), (
                    f"strand={strand} rate={rate} orient={orientation}: "
                    f"exon length {_total_length(truncated_exons)} != "
                    f"seq length {len(tx_seq)}"
                )
