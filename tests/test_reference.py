"""Tests for reference quality and curated-parent priority metadata."""

import pandas as pd
import pyranges as pr

from craft.core.reference import transcript_metadata


def test_mane_is_prioritized_and_incomplete_cds_is_flagged() -> None:
    columns = [
        "Chromosome", "Start", "End", "Strand", "transcript_id", "Feature",
        "gene_id", "Frame", "tag",
    ]
    rows = [
        ("chr1", 0, 100, "+", "mane", "exon", "g1", ".", "MANE_Select,basic"),
        ("chr1", 10, 91, "+", "mane", "CDS", "g1", "0", "MANE_Select,basic"),
        ("chr1", 10, 13, "+", "mane", "start_codon", "g1", "0", "MANE_Select,basic"),
        ("chr1", 91, 94, "+", "mane", "stop_codon", "g1", "0", "MANE_Select,basic"),
        ("chr1", 0, 100, "+", "partial", "exon", "g1", ".", "cds_start_NF"),
        ("chr1", 20, 90, "+", "partial", "CDS", "g1", "0", "cds_start_NF"),
    ]
    meta = transcript_metadata(pr.PyRanges(pd.DataFrame(rows, columns=columns))).set_index(
        "transcript_id"
    )
    assert meta.loc["mane", "reference_priority"] > meta.loc["partial", "reference_priority"]
    assert bool(meta.loc["mane", "reference_cds_complete"])
    assert not bool(meta.loc["partial", "reference_cds_complete"])


def test_cds_phase_continuity_is_checked_in_transcript_order() -> None:
    columns = [
        "Chromosome", "Start", "End", "Strand", "transcript_id", "Feature",
        "gene_id", "Frame",
    ]
    rows = [
        ("chr1", 0, 5, "+", "valid", "CDS", "g1", "0"),
        ("chr1", 10, 16, "+", "valid", "CDS", "g1", "1"),
        ("chr1", 0, 5, "+", "invalid", "CDS", "g1", "0"),
        ("chr1", 10, 16, "+", "invalid", "CDS", "g1", "0"),
    ]
    meta = transcript_metadata(pr.PyRanges(pd.DataFrame(rows, columns=columns))).set_index(
        "transcript_id"
    )
    assert bool(meta.loc["valid", "reference_cds_phase_valid"])
    assert not bool(meta.loc["invalid", "reference_cds_phase_valid"])
