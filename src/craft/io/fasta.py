"""FASTA access wrappers around pysam."""

from pathlib import Path

import pysam


def open_genome(path: Path) -> pysam.FastaFile:
    """Open an indexed genome FASTA, building the ``.fai`` index if missing.

    Args:
        path: Path to a FASTA file. Compressed files (``.fa.gz``) require a
            BGZF-compressed FASTA with a ``.gzi`` index alongside the ``.fai``;
            pysam will raise on that case.

    Returns:
        Open :class:`pysam.FastaFile` handle. The caller is responsible for
        closing it (``handle.close()`` or use it as a context manager).
    """
    fai = Path(f"{path}.fai")
    if not fai.exists():
        pysam.faidx(str(path))
    return pysam.FastaFile(str(path))
