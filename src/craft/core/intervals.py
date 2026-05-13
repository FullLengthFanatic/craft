"""PyRanges helpers for splice-junction extraction and interval ops."""

import pandas as pd
import pyranges as pr


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
