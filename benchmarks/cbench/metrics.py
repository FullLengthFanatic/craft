"""Scoring functions for Bench 1.

Compares ``ORFCall`` predictions against GENCODE-truth ORFs and aggregates per
benchmark cell (rate, orientation, comparator). Outputs are pandas DataFrames
ready to feed into the figure helpers.
"""

from dataclasses import dataclass

import pandas as pd

from cbench.comparators import ORFCall


@dataclass(frozen=True)
class ORFScoreRow:
    transcript_id: str
    comparator: str
    rate: float
    orientation: str
    seed: int
    found: bool
    start_exact: bool
    start_in_frame: bool
    stop_exact: bool
    stop_in_frame: bool
    length_error_nt: int  # predicted_length - truth_length
    truth_length_nt: int


_NO_SCORE = ORFScoreRow("", "", 0.0, "", 0, False, False, False, False, False, 0, 0)


def score_one(
    truth: ORFCall,
    pred: ORFCall,
    comparator: str,
    rate: float,
    orientation: str,
    seed: int,
    transcript_id: str,
) -> ORFScoreRow:
    """Score a single prediction against the truth ORF for one transcript.

    ``truth.tx_start`` and ``pred.tx_start`` are in the GENCODE transcript
    coordinate system; the truncation simulator preserves the same frame by
    slicing the transcript, so ``pred`` from CRAFT or orfipy is also expressed
    against that same coordinate system after we map back from the truncated iso.
    """
    truth_length = truth.tx_end - truth.tx_start
    if not pred.found:
        return ORFScoreRow(
            transcript_id=transcript_id,
            comparator=comparator,
            rate=rate,
            orientation=orientation,
            seed=seed,
            found=False,
            start_exact=False,
            start_in_frame=False,
            stop_exact=False,
            stop_in_frame=False,
            length_error_nt=-truth_length,
            truth_length_nt=truth_length,
        )
    start_exact = pred.tx_start == truth.tx_start
    stop_exact = pred.tx_end == truth.tx_end
    start_in_frame = (pred.tx_start - truth.tx_start) % 3 == 0
    stop_in_frame = (pred.tx_end - truth.tx_end) % 3 == 0
    return ORFScoreRow(
        transcript_id=transcript_id,
        comparator=comparator,
        rate=rate,
        orientation=orientation,
        seed=seed,
        found=True,
        start_exact=start_exact,
        start_in_frame=start_in_frame,
        stop_exact=stop_exact,
        stop_in_frame=stop_in_frame,
        length_error_nt=(pred.tx_end - pred.tx_start) - truth_length,
        truth_length_nt=truth_length,
    )


def to_dataframe(rows: list[ORFScoreRow]) -> pd.DataFrame:
    """Pack a list of ``ORFScoreRow`` into a pandas DataFrame."""
    if not rows:
        return pd.DataFrame(columns=list(_NO_SCORE.__dataclass_fields__.keys()))
    return pd.DataFrame([r.__dict__ for r in rows])


def summarize_cells(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate per-row scores into (comparator, rate, orientation, seed) cells.

    Returns a long DataFrame with columns: ``comparator``, ``rate``, ``orientation``,
    ``seed``, ``n``, ``recovery_rate``, ``start_exact_rate``, ``stop_exact_rate``,
    ``start_in_frame_rate``, ``stop_in_frame_rate``, ``mean_abs_length_error``.
    """
    if df.empty:
        return pd.DataFrame()
    grouped = df.groupby(["comparator", "rate", "orientation", "seed"], dropna=False)
    out = grouped.agg(
        n=("found", "size"),
        recovery_rate=("found", "mean"),
        start_exact_rate=("start_exact", "mean"),
        start_in_frame_rate=("start_in_frame", "mean"),
        stop_exact_rate=("stop_exact", "mean"),
        stop_in_frame_rate=("stop_in_frame", "mean"),
        mean_abs_length_error=(
            "length_error_nt",
            lambda s: float(s.abs().mean()),
        ),
    )
    return out.reset_index()


def confidence_calibration(
    score_df: pd.DataFrame,
    confidence_df: pd.DataFrame,
    on: str = "transcript_id",
) -> pd.DataFrame:
    """Join CRAFT confidence labels onto the score table and bin by HIGH/MED/LOW.

    Args:
        score_df: per-row scores from ``to_dataframe`` (only rows with comparator
            == "craft" make sense here).
        confidence_df: DataFrame with at least columns ``transcript_id`` and
            ``orf_confidence`` (from CRAFT's per_isoform.tsv).
        on: join key.

    Returns:
        DataFrame with columns ``orf_confidence``, ``n``, ``recovery_rate``,
        ``start_exact_rate``, ``stop_exact_rate``.
    """
    craft_rows = score_df[score_df["comparator"] == "craft"].copy()
    if craft_rows.empty:
        return pd.DataFrame()
    joined = craft_rows.merge(
        confidence_df[[on, "orf_confidence"]], on=on, how="left"
    )
    joined["orf_confidence"] = joined["orf_confidence"].fillna("missing")
    grouped = joined.groupby("orf_confidence")
    return grouped.agg(
        n=("found", "size"),
        recovery_rate=("found", "mean"),
        start_exact_rate=("start_exact", "mean"),
        stop_exact_rate=("stop_exact", "mean"),
    ).reset_index()
