"""Tests for craft.core.orf.confidence."""

import pytest

from craft.core.completeness import Completeness
from craft.core.orf.confidence import ORFConfidence, score
from craft.core.orf.propagation import ORFOutcome


def test_full_length_intact_is_high_with_max_score() -> None:
    category, value = score(Completeness.FULL_LENGTH, ORFOutcome.PROPAGATED_INTACT)
    assert category == ORFConfidence.HIGH
    assert value == pytest.approx(1.0)


def test_truncated_5p_intact_is_high() -> None:
    category, value = score(Completeness.TRUNCATED_5P, ORFOutcome.PROPAGATED_INTACT)
    assert category == ORFConfidence.HIGH
    assert 0.85 <= value < 1.0


def test_truncated_3p_intact_is_high() -> None:
    category, value = score(Completeness.TRUNCATED_3P, ORFOutcome.PROPAGATED_INTACT)
    assert category == ORFConfidence.HIGH


def test_truncated_both_intact_is_medium() -> None:
    category, _ = score(Completeness.TRUNCATED_BOTH, ORFOutcome.PROPAGATED_INTACT)
    assert category == ORFConfidence.MEDIUM


def test_internal_fragment_intact_is_medium() -> None:
    category, _ = score(Completeness.INTERNAL_FRAGMENT, ORFOutcome.PROPAGATED_INTACT)
    assert category == ORFConfidence.MEDIUM


def test_full_length_stop_not_observed_is_medium() -> None:
    category, _ = score(Completeness.FULL_LENGTH, ORFOutcome.STOP_NOT_OBSERVED)
    assert category == ORFConfidence.MEDIUM


def test_truncated_3p_stop_not_observed_is_low() -> None:
    category, _ = score(Completeness.TRUNCATED_3P, ORFOutcome.STOP_NOT_OBSERVED)
    assert category == ORFConfidence.LOW


def test_full_length_disrupted_is_low() -> None:
    category, _ = score(Completeness.FULL_LENGTH, ORFOutcome.DISRUPTED)
    assert category == ORFConfidence.LOW


def test_full_length_start_lost_is_low() -> None:
    category, value = score(Completeness.FULL_LENGTH, ORFOutcome.START_LOST)
    assert category == ORFConfidence.LOW
    assert value < 0.5


def test_no_parent_is_none_with_zero_score() -> None:
    category, value = score(Completeness.NOVEL_NO_MATCH, ORFOutcome.NO_PARENT)
    assert category == ORFConfidence.NONE
    assert value == 0.0


def test_no_parent_cds_is_none_with_zero_score() -> None:
    category, value = score(Completeness.FULL_LENGTH, ORFOutcome.NO_PARENT_CDS)
    assert category == ORFConfidence.NONE
    assert value == 0.0


def test_accepts_string_inputs() -> None:
    category, value = score("full_length", "propagated_intact")
    assert category == ORFConfidence.HIGH
    assert value == pytest.approx(1.0)


def test_score_is_in_unit_interval_for_all_combinations() -> None:
    for comp in Completeness:
        for outcome in ORFOutcome:
            _, value = score(comp, outcome)
            assert 0.0 <= value <= 1.0, f"Score out of range for {comp}/{outcome}: {value}"


def test_high_threshold_strictly_above_medium_threshold() -> None:
    # Sanity: PROPAGATED_INTACT scores are always >= STOP_NOT_OBSERVED scores
    # for the same completeness level.
    for comp in (
        Completeness.FULL_LENGTH,
        Completeness.TRUNCATED_5P,
        Completeness.TRUNCATED_3P,
        Completeness.TRUNCATED_BOTH,
        Completeness.INTERNAL_FRAGMENT,
    ):
        _, intact_value = score(comp, ORFOutcome.PROPAGATED_INTACT)
        _, stop_lost_value = score(comp, ORFOutcome.STOP_NOT_OBSERVED)
        _, disrupted_value = score(comp, ORFOutcome.DISRUPTED)
        _, start_lost_value = score(comp, ORFOutcome.START_LOST)
        assert intact_value >= stop_lost_value
        assert intact_value >= disrupted_value
        assert intact_value >= start_lost_value


def test_full_length_scores_at_least_as_high_as_truncated_for_same_outcome() -> None:
    for outcome in (
        ORFOutcome.PROPAGATED_INTACT,
        ORFOutcome.STOP_NOT_OBSERVED,
        ORFOutcome.DISRUPTED,
        ORFOutcome.START_LOST,
    ):
        _, fl_value = score(Completeness.FULL_LENGTH, outcome)
        for comp in (
            Completeness.TRUNCATED_5P,
            Completeness.TRUNCATED_3P,
            Completeness.TRUNCATED_BOTH,
            Completeness.INTERNAL_FRAGMENT,
        ):
            _, partial_value = score(comp, outcome)
            assert fl_value >= partial_value, (
                f"FL/{outcome} score ({fl_value}) "
                f"< {comp}/{outcome} score ({partial_value})"
            )
