"""Pfam domain disruption via local pyhmmer scanning.

Hits are cached by protein-sequence SHA256 so repeated proteins across cells or samples
scan only once.
"""

from pathlib import Path

import pyranges as pr


def scan(orfs: pr.PyRanges, pfam_hmm: Path, genome_fasta: Path) -> pr.PyRanges:
    """Translate ORFs and scan against Pfam-A.hmm; emit per-domain preservation status."""
    raise NotImplementedError
