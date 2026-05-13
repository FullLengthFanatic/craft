"""PyRanges helpers for splice-junction extraction and interval ops."""

import pyranges as pr


def splice_junctions(exons: pr.PyRanges) -> pr.PyRanges:
    """Compute splice junctions from a set of exon intervals grouped by transcript."""
    raise NotImplementedError
