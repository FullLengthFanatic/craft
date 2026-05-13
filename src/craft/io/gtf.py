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

import pyranges as pr


def _normalise(df, columns: list[str]) -> pr.PyRanges:
    df = df[columns].copy()
    df["Strand"] = df["Strand"].astype(str)
    if "Feature" in df.columns:
        df["Feature"] = df["Feature"].astype(str)
    return pr.PyRanges(df)


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
    return _normalise(df, ["Chromosome", "Start", "End", "Strand", "transcript_id"])


def load_reference(path: Path) -> pr.PyRanges:
    """Load reference annotation records from a GTF.

    Filters to ``Feature in {"exon", "CDS"}`` rows and returns a PyRanges with
    columns ``Chromosome``, ``Start``, ``End``, ``Strand``, ``transcript_id``,
    ``Feature``.

    Args:
        path: GTF file path (plain or gzipped).

    Returns:
        PyRanges of exon and CDS intervals.
    """
    gr = pr.read_gtf(str(path))
    df = gr.df
    df = df[df["Feature"].isin(["exon", "CDS"])]
    return _normalise(
        df, ["Chromosome", "Start", "End", "Strand", "transcript_id", "Feature"]
    )
