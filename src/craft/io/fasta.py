"""FASTA access wrappers around pysam."""

from pathlib import Path

import pysam


def open_genome(path: Path) -> pysam.FastaFile:
    """Open an indexed genome FASTA; build the .fai index if missing."""
    raise NotImplementedError
