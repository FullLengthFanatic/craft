"""BAM access wrappers for optional read-end position checks."""

from pathlib import Path

import pysam


def open_alignments(path: Path) -> pysam.AlignmentFile:
    """Open a coordinate-sorted BAM with index."""
    raise NotImplementedError
