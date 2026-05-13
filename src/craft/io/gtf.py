"""GTF parsing for isoform and reference annotations."""

from pathlib import Path

import pyranges as pr


def load_isoforms(path: Path) -> pr.PyRanges:
    """Load an isoform GTF into a PyRanges with exon and transcript records."""
    raise NotImplementedError


def load_reference(path: Path) -> pr.PyRanges:
    """Load a reference annotation GTF (GENCODE / Ensembl) with CDS records."""
    raise NotImplementedError
