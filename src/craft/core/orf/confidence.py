"""ORF confidence scoring as a function of completeness and propagation success."""

from enum import Enum


class ORFConfidence(str, Enum):
    """Per-isoform ORF call confidence."""

    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


def score(
    completeness: str,
    propagation_success: bool,
    stop_codon_observed: bool,
) -> tuple[ORFConfidence, float]:
    """Compute a categorical and numeric (0-1) ORF confidence."""
    raise NotImplementedError
