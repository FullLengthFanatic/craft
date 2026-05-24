"""Tests for cbench.metrics."""


from cbench.comparators import ORFCall
from cbench.metrics import score_one, summarize_cells, to_dataframe


def _truth(start: int = 30, end: int = 300) -> ORFCall:
    return ORFCall(
        transcript_id="t1",
        found=True,
        tx_start=start,
        tx_end=end,
        start_codon="ATG",
        stop_codon="TAA",
    )


def _pred(found: bool, start: int = 0, end: int = 0) -> ORFCall:
    return ORFCall(
        transcript_id="t1",
        found=found,
        tx_start=start,
        tx_end=end,
        start_codon="ATG" if found else "",
        stop_codon="TAA" if found else "",
    )


def test_score_one_no_call_marks_full_length_loss() -> None:
    row = score_one(
        _truth(30, 300), _pred(False), "orfipy", 0.25, "3prime", 0, "t1"
    )
    assert row.found is False
    assert row.start_exact is False
    assert row.stop_exact is False
    assert row.length_error_nt == -(300 - 30)


def test_score_one_exact_match() -> None:
    row = score_one(
        _truth(30, 300), _pred(True, 30, 300), "craft", 0.10, "5prime", 0, "t1"
    )
    assert row.found is True
    assert row.start_exact is True
    assert row.stop_exact is True
    assert row.start_in_frame is True
    assert row.stop_in_frame is True
    assert row.length_error_nt == 0


def test_score_one_off_by_one_breaks_frame() -> None:
    row = score_one(
        _truth(30, 300), _pred(True, 31, 300), "craft", 0.10, "5prime", 0, "t1"
    )
    assert row.start_exact is False
    assert row.start_in_frame is False
    assert row.stop_exact is True


def test_score_one_off_by_three_stays_in_frame() -> None:
    row = score_one(
        _truth(30, 300), _pred(True, 33, 300), "craft", 0.10, "5prime", 0, "t1"
    )
    assert row.start_exact is False
    assert row.start_in_frame is True


def test_summarize_cells_recovery_and_length_error() -> None:
    rows = [
        score_one(_truth(30, 300), _pred(True, 30, 300), "craft", 0.1, "5prime", 0, "t1"),
        score_one(_truth(30, 300), _pred(False), "craft", 0.1, "5prime", 0, "t2"),
        score_one(_truth(30, 300), _pred(True, 30, 297), "craft", 0.1, "5prime", 0, "t3"),
    ]
    df = to_dataframe(rows)
    summary = summarize_cells(df)
    assert len(summary) == 1
    cell = summary.iloc[0]
    assert cell["n"] == 3
    assert cell["recovery_rate"] == pytest_approx(2 / 3)
    # |length_error|: 0 + 270 + 3 = 273; mean = 91.0
    assert cell["mean_abs_length_error"] == pytest_approx(273 / 3)


def test_to_dataframe_empty_returns_typed_empty() -> None:
    df = to_dataframe([])
    assert df.empty
    for col in ("transcript_id", "comparator", "rate", "orientation"):
        assert col in df.columns


def pytest_approx(v):
    import pytest

    return pytest.approx(v, rel=1e-6)
