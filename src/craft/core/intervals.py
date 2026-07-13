"""PyRanges helpers for splice-junction extraction and interval ops."""

import pandas as pd
import pyranges as pr


def transcript_coordinate(exons: pd.DataFrame, genomic_pos: int, strand: str) -> int | None:
    """Map a genomic base to a 0-based coordinate in spliced transcript order."""
    if strand not in {"+", "-"}:
        raise ValueError(f"Unsupported strand: {strand!r}")
    ordered = exons.sort_values("Start", ascending=strand == "+", kind="stable")
    offset = 0
    for exon in ordered.itertuples(index=False):
        start = int(exon.Start)
        end = int(exon.End)
        if start <= genomic_pos < end:
            within = genomic_pos - start if strand == "+" else end - 1 - genomic_pos
            return offset + within
        offset += end - start
    return None


def genomic_position_at_transcript_coordinate(
    exons: pd.DataFrame, coordinate: int, strand: str
) -> int | None:
    """Map one 0-based spliced-transcript coordinate back to the genome."""
    if coordinate < 0:
        return None
    if strand not in {"+", "-"}:
        raise ValueError(f"Unsupported strand: {strand!r}")
    ordered = exons.sort_values("Start", ascending=strand == "+", kind="stable")
    offset = 0
    for exon in ordered.itertuples(index=False):
        start = int(exon.Start)
        end = int(exon.End)
        length = end - start
        if coordinate < offset + length:
            within = coordinate - offset
            return start + within if strand == "+" else end - 1 - within
        offset += length
    return None


def spliced_length(exons: pd.DataFrame) -> int:
    """Total exonic length of one transcript."""
    return int((exons["End"] - exons["Start"]).sum())


def splice_junctions(exons: pr.PyRanges) -> pr.PyRanges:
    """Compute splice junctions from exon intervals grouped by transcript.

    For each transcript, junctions are the introns between consecutive exons in
    genomic order. Single-exon transcripts contribute no junctions. Output coordinates
    are half-open intron coordinates: ``Start`` is the upstream exon's ``End``,
    ``End`` is the downstream exon's ``Start``.

    Args:
        exons: PyRanges of exons with ``transcript_id`` and ``Strand`` columns.

    Returns:
        PyRanges with columns ``Chromosome``, ``Start``, ``End``, ``Strand``,
        ``transcript_id``, ``junction_index`` (0-based, in genomic order within the
        transcript). Empty PyRanges if no junctions are present.
    """
    if len(exons) == 0:
        return pr.PyRanges()

    df = exons.df.sort_values(["transcript_id", "Start"], kind="stable").reset_index(drop=True)
    df["next_start"] = df.groupby("transcript_id")["Start"].shift(-1)
    df["junction_index"] = df.groupby("transcript_id").cumcount()

    junctions = df.dropna(subset=["next_start"]).copy()
    junctions = junctions[junctions["next_start"] > junctions["End"]]

    if junctions.empty:
        return pr.PyRanges()

    out = pd.DataFrame(
        {
            "Chromosome": junctions["Chromosome"].to_numpy(),
            "Start": junctions["End"].astype("int64").to_numpy(),
            "End": junctions["next_start"].astype("int64").to_numpy(),
            "Strand": junctions["Strand"].to_numpy(),
            "transcript_id": junctions["transcript_id"].to_numpy(),
            "junction_index": junctions["junction_index"].astype("int64").to_numpy(),
        }
    )
    return pr.PyRanges(out)
