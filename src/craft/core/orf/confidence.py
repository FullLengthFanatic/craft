"""ORF confidence scoring as a function of completeness and propagation outcome."""

from enum import Enum

from craft.core.completeness import Completeness
from craft.core.orf.propagation import ORFOutcome

HIGH_THRESHOLD = 0.85
MEDIUM_THRESHOLD = 0.5

_BASE_BY_OUTCOME: dict[ORFOutcome, float] = {
    ORFOutcome.PROPAGATED_INTACT: 1.0,
    ORFOutcome.STOP_NOT_OBSERVED: 0.55,
    ORFOutcome.DISRUPTED: 0.45,
    ORFOutcome.START_LOST: 0.2,
}

_COMPLETENESS_FACTOR: dict[Completeness, float] = {
    Completeness.FULL_LENGTH: 1.0,
    Completeness.TRUNCATED_5P: 0.9,
    Completeness.TRUNCATED_3P: 0.9,
    Completeness.TRUNCATED_BOTH: 0.65,
    Completeness.INTERNAL_FRAGMENT: 0.5,
    Completeness.NOVEL_NO_MATCH: 0.0,
}


class ORFConfidence(str, Enum):
    """Per-isoform ORF call confidence."""

    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"
    NONE = "none"


def score(
    completeness: Completeness | str,
    orf_outcome: ORFOutcome | str,
) -> tuple[ORFConfidence, float]:
    """Compute a categorical and numeric (0-1) ORF confidence.

    The numeric score combines the propagation outcome (PROPAGATED_INTACT is the
    strongest signal; START_LOST the weakest among attempted calls) with a
    completeness factor that penalises truncation and internal-fragment cases.
    Outcomes ``NO_PARENT`` and ``NO_PARENT_CDS`` short-circuit to
    ``ORFConfidence.NONE`` with score 0.0 regardless of completeness.

    Args:
        completeness: A :class:`Completeness` enum value or its string equivalent.
        orf_outcome: An :class:`ORFOutcome` enum value or its string equivalent.

    Returns:
        Tuple of (confidence category, numeric score in [0.0, 1.0]).
    """
    comp = Completeness(completeness) if isinstance(completeness, str) else completeness
    outcome = ORFOutcome(orf_outcome) if isinstance(orf_outcome, str) else orf_outcome

    if outcome in (ORFOutcome.NO_PARENT, ORFOutcome.NO_PARENT_CDS):
        return ORFConfidence.NONE, 0.0
    if comp == Completeness.NOVEL_NO_MATCH:
        return ORFConfidence.NONE, 0.0

    value = _BASE_BY_OUTCOME[outcome] * _COMPLETENESS_FACTOR[comp]

    if value >= HIGH_THRESHOLD:
        category = ORFConfidence.HIGH
    elif value >= MEDIUM_THRESHOLD:
        category = ORFConfidence.MEDIUM
    else:
        category = ORFConfidence.LOW
    return category, value
