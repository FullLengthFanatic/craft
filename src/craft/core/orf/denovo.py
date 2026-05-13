"""De novo ORF prediction for genuinely novel isoforms (no usable parent)."""

from pathlib import Path

import pyranges as pr


def predict(isoforms: pr.PyRanges, genome_fasta: Path) -> pr.PyRanges:
    """Predict ORFs de novo for isoforms with no usable reference parent."""
    raise NotImplementedError
