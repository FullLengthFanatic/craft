"""3' UTR feature consequences: length delta and poly(A) signal motif scan."""

from pathlib import Path

import pyranges as pr

POLYA_SIGNALS: tuple[str, ...] = (
    "AATAAA",
    "ATTAAA",
    "AGTAAA",
    "TATAAA",
    "CATAAA",
    "GATAAA",
    "AATATA",
    "AATACA",
    "AATAGA",
    "AAAAAG",
    "ACTAAA",
)


def annotate(
    isoforms: pr.PyRanges,
    reference: pr.PyRanges,
    genome_fasta: Path,
) -> pr.PyRanges:
    """Emit UTR length delta vs parent and the strongest poly(A) signal per isoform."""
    raise NotImplementedError
