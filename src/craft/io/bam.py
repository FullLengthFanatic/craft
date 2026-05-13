"""BAM access wrappers for optional read-end position checks."""

from pathlib import Path

import pysam


def open_alignments(path: Path) -> pysam.AlignmentFile:
    """Open a coordinate-sorted BAM, building the ``.bai`` index if missing.

    The BAM must already be coordinate-sorted; ``pysam.index`` will raise if it
    is not.

    Args:
        path: Path to a coordinate-sorted BAM.

    Returns:
        Open :class:`pysam.AlignmentFile` handle (read mode). The caller is
        responsible for closing it.
    """
    bai_alongside = Path(f"{path}.bai")
    bai_renamed = path.with_suffix(".bai")
    if not bai_alongside.exists() and not bai_renamed.exists():
        pysam.index(str(path))
    return pysam.AlignmentFile(str(path), "rb")
