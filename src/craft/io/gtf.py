"""GTF parsing for isoform and reference annotations.

Wraps :func:`pyranges.read_gtf` and normalises the output to the minimal column
set the downstream core modules expect (Chromosome, Start, End, Strand,
transcript_id, plus Feature for the reference loader). Coordinates are returned
as PyRanges 0-based half-open, matching how :func:`pyranges.read_gtf` already
converts from 1-based-inclusive GTF.

``Strand`` is normalised to plain string dtype to avoid the PyRanges 0.1.4
quirk where a categorical column with unused levels marks the PyRanges as
unstranded after a downstream `.df` round-trip.
"""

from pathlib import Path

import pandas as pd
import pyranges as pr


def _normalise(df, columns: list[str]) -> pr.PyRanges:
    df = df[columns].copy()
    df["Strand"] = df["Strand"].astype(str)
    if "Feature" in df.columns:
        df["Feature"] = df["Feature"].astype(str)
    return pr.PyRanges(df)


def _present(df: pd.DataFrame, required: list[str], optional: list[str]) -> list[str]:
    """Return required columns plus optional annotation metadata that are present."""
    return [*required, *(column for column in optional if column in df.columns)]


_ISOFORM_OPTIONAL = ["gene_id", "gene_name"]

# Keep the reference fields needed to judge whether a CDS is complete and to
# prioritize curated parents.  PyRanges/read_gtf exposes GTF column 8 as Frame.
_REFERENCE_OPTIONAL = [
    "gene_name",
    "Frame",
    "transcript_type",
    "transcript_biotype",
    "gene_type",
    "gene_biotype",
    "transcript_support_level",
    "tag",
    "ccdsid",
    "havana_transcript",
]


def load_isoforms(path: Path) -> pr.PyRanges:
    """Load isoform exon records from a GTF.

    Filters to ``Feature == "exon"`` rows and returns a PyRanges with columns
    ``Chromosome``, ``Start``, ``End``, ``Strand``, ``transcript_id``.

    Args:
        path: GTF file path (plain or gzipped; pyranges handles both).

    Returns:
        PyRanges of exon intervals.
    """
    gr = pr.read_gtf(str(path))
    df = gr.df
    df = df[df["Feature"] == "exon"]
    required = ["Chromosome", "Start", "End", "Strand", "transcript_id"]
    return _normalise(df, _present(df, required, _ISOFORM_OPTIONAL))


def load_reference(path: Path) -> pr.PyRanges:
    """Load reference annotation records from a GTF.

    Keeps exon, CDS, explicit start-codon and explicit stop-codon rows.  It also
    preserves CDS phase and transcript-quality fields when present; downstream
    code must not silently treat an incomplete reference CDS as complete.
    columns ``Chromosome``, ``Start``, ``End``, ``Strand``, ``transcript_id``,
    ``Feature``, ``gene_id``, and (if the GTF has it) ``gene_name``.

    Args:
        path: GTF file path (plain or gzipped).

    Returns:
        PyRanges of exon and CDS intervals.
    """
    gr = pr.read_gtf(str(path))
    df = gr.df
    df = df[df["Feature"].isin(["exon", "CDS", "start_codon", "stop_codon"])]
    required = [
        "Chromosome", "Start", "End", "Strand", "transcript_id", "Feature", "gene_id"
    ]
    columns = _present(df, required, _REFERENCE_OPTIONAL)
    return _normalise(df, columns)
