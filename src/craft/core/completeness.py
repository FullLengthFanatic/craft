"""Completeness classification: full-length vs truncated isoforms vs reference."""

from enum import Enum

import pyranges as pr


class Completeness(str, Enum):
    """Structural completeness of a novel isoform relative to a reference parent."""

    FULL_LENGTH = "full_length"
    TRUNCATED_5P = "truncated_5p"
    TRUNCATED_3P = "truncated_3p"
    TRUNCATED_BOTH = "truncated_both"
    INTERNAL_FRAGMENT = "internal_fragment"
    NOVEL_NO_MATCH = "novel_no_match"


def classify(isoforms: pr.PyRanges, reference: pr.PyRanges) -> pr.PyRanges:
    """Classify each isoform's completeness vs its best-matching reference transcript."""
    raise NotImplementedError
