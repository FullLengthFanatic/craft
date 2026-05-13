"""Nonsense-mediated decay (NMD) susceptibility prediction.

Rule set: 50nt PTC rule, long-last-exon escape (>400 nt), start-proximal escape
(<150 nt from start codon), and explicit handling of intron-retention stop codons.
"""

from enum import Enum

import pyranges as pr


class NMDStatus(str, Enum):
    """NMD susceptibility class."""

    SENSITIVE = "sensitive"
    ESCAPED = "escaped"
    NOT_APPLICABLE = "not_applicable"


def predict(orfs: pr.PyRanges) -> pr.PyRanges:
    """Apply NMD rules with confidence flags that drop on 3' truncation."""
    raise NotImplementedError
