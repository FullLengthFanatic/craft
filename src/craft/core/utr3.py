"""3' UTR feature consequences: length delta vs parent and poly(A) signal motif scan.

Internal priming detection stays in `tecap`; this module only reports structural and
sequence features of the 3' UTR.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pyranges as pr
import pysam

from craft.core.orf.propagation import ORFOutcome

POLYA_SIGNALS: tuple[str, ...] = (
    "AATAAA",
    "ATTAAA",
    "AGTAAA",
    "TATAAA",
    "CATAAA",
    "GATAAA",
    "AATATA",
    "AATACA",
    "AATAGA",
    "AAAAAG",
    "ACTAAA",
)

_RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _reverse_complement(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]


def polya_signal(sequence: str) -> dict[str, int | str]:
    """Find the strongest poly(A) signal motif in a 3' UTR sequence.

    Motifs are searched in priority order (canonical AATAAA first, then known
    variants). For each motif, the rightmost (most 3'-proximal) occurrence is
    chosen. The first motif with any occurrence wins, so priority over distance.

    Args:
        sequence: 3' UTR sequence in transcript orientation (5' to 3').

    Returns:
        Dict with keys ``motif`` (empty string if none) and
        ``distance_from_3p_end`` (-1 if none; otherwise nt from motif end to
        sequence end).
    """
    upper = sequence.upper()
    for motif in POLYA_SIGNALS:
        idx = upper.rfind(motif)
        if idx >= 0:
            distance = len(upper) - (idx + len(motif))
            return {"motif": motif, "distance_from_3p_end": distance}
    return {"motif": "", "distance_from_3p_end": -1}


def _utr3_length(exons: pd.DataFrame, stop_pos: int, strand: str) -> int:
    """Total mRNA bp downstream of the stop codon (in transcript order)."""
    starts = exons["Start"].to_numpy()
    ends = exons["End"].to_numpy()
    if strand == "+":
        contrib = np.maximum(0, ends - np.maximum(starts, stop_pos + 1))
    elif strand == "-":
        contrib = np.maximum(0, np.minimum(ends, stop_pos) - starts)
    else:
        raise ValueError(f"Unsupported strand: {strand!r}")
    return int(contrib.sum())


def _parent_stop_pos(cds_df: pd.DataFrame, strand: str) -> int:
    if strand == "+":
        return int(cds_df["End"].max()) - 1
    if strand == "-":
        return int(cds_df["Start"].min())
    raise ValueError(f"Unsupported strand: {strand!r}")


def _iso_stop_pos(intervals: list[tuple], strand: str) -> int:
    if strand == "+":
        return max(end - 1 for _, _, end, _ in intervals)
    if strand == "-":
        return min(start for _, start, _, _ in intervals)
    raise ValueError(f"Unsupported strand: {strand!r}")


def _extract_utr3_sequence(
    exons: pd.DataFrame,
    stop_pos: int,
    strand: str,
    genome: pysam.FastaFile,
) -> str:
    """Extract the 3' UTR sequence in transcript orientation (5' to 3')."""
    chrom = str(exons["Chromosome"].iloc[0])
    sorted_exons = exons.sort_values("Start").reset_index(drop=True)
    parts: list[str] = []
    for _, ex in sorted_exons.iterrows():
        ex_start = int(ex["Start"])
        ex_end = int(ex["End"])
        if strand == "+":
            utr_start = max(ex_start, stop_pos + 1)
            utr_end = ex_end
        else:
            utr_start = ex_start
            utr_end = min(ex_end, stop_pos)
        if utr_end > utr_start:
            parts.append(genome.fetch(chrom, utr_start, utr_end))
    sequence = "".join(parts).upper()
    if strand == "-":
        sequence = _reverse_complement(sequence)
    return sequence


def annotate(
    classified: pr.PyRanges,
    propagated: pd.DataFrame,
    reference: pr.PyRanges,
    genome_fasta: Path | None = None,
) -> pd.DataFrame:
    """3' UTR feature consequences per isoform.

    Computes the isoform's 3' UTR length, its parent's 3' UTR length (when a
    parent transcript and parent CDS are available), the absolute and percent
    deltas, and, if ``genome_fasta`` is provided, the strongest poly(A) signal
    motif in the isoform's 3' UTR plus its distance from the 3' end.

    Args:
        classified: PyRanges of isoform exons with ``transcript_id``, ``Strand``.
        propagated: DataFrame returned by :func:`craft.core.orf.propagation.propagate`.
        reference: Reference PyRanges with a ``Feature`` column (exon / CDS rows).
        genome_fasta: Optional path to an indexed genome FASTA. If omitted, the
            poly(A) scan is skipped.

    Returns:
        DataFrame with one row per isoform and columns: ``transcript_id``,
        ``iso_utr3_length_nt``, ``parent_utr3_length_nt``, ``utr3_length_delta_nt``,
        ``utr3_length_delta_pct``, ``polya_signal_motif``, ``polya_signal_distance_nt``.
    """
    cols = [
        "transcript_id",
        "iso_utr3_length_nt",
        "parent_utr3_length_nt",
        "utr3_length_delta_nt",
        "utr3_length_delta_pct",
        "polya_signal_motif",
        "polya_signal_distance_nt",
    ]
    if propagated.empty or len(classified) == 0:
        return pd.DataFrame(columns=cols)

    iso_df = classified.df
    iso_strand = iso_df.groupby("transcript_id")["Strand"].first().to_dict()
    iso_exons_by_tx = {tx: g for tx, g in iso_df.groupby("transcript_id", sort=False)}

    ref_df = reference.df
    if "Feature" in ref_df.columns and "transcript_id" in ref_df.columns:
        parent_exons_all = ref_df[ref_df["Feature"] == "exon"]
        parent_cds_all = ref_df[ref_df["Feature"] == "CDS"]
        parent_exons_by_tx = {
            tx: g for tx, g in parent_exons_all.groupby("transcript_id", sort=False)
        }
        parent_cds_by_tx = {
            tx: g for tx, g in parent_cds_all.groupby("transcript_id", sort=False)
        }
    else:
        parent_exons_by_tx = {}
        parent_cds_by_tx = {}

    genome = pysam.FastaFile(str(genome_fasta)) if genome_fasta is not None else None

    try:
        rows: list[dict] = []
        for _, prop_row in propagated.iterrows():
            tx_id = prop_row["transcript_id"]
            outcome = ORFOutcome(prop_row["orf_outcome"])
            intervals = prop_row["propagated_cds_intervals"]
            stop_covered = bool(prop_row["stop_codon_covered"])
            parent_tx = prop_row["parent_tx_id"]

            applicable = (
                outcome in (ORFOutcome.PROPAGATED_INTACT, ORFOutcome.DISRUPTED)
                and stop_covered
                and intervals
            )
            if not applicable:
                rows.append(
                    {
                        "transcript_id": tx_id,
                        "iso_utr3_length_nt": None,
                        "parent_utr3_length_nt": None,
                        "utr3_length_delta_nt": None,
                        "utr3_length_delta_pct": None,
                        "polya_signal_motif": "",
                        "polya_signal_distance_nt": None,
                    }
                )
                continue

            strand = str(iso_strand[tx_id])
            iso_exons = iso_exons_by_tx[tx_id]
            iso_stop = _iso_stop_pos(intervals, strand)
            iso_utr_len = _utr3_length(iso_exons, iso_stop, strand)

            parent_utr_len: int | None = None
            if (
                parent_tx
                and parent_tx in parent_exons_by_tx
                and parent_tx in parent_cds_by_tx
            ):
                parent_exons = parent_exons_by_tx[parent_tx]
                parent_cds = parent_cds_by_tx[parent_tx]
                parent_stop = _parent_stop_pos(parent_cds, strand)
                parent_utr_len = _utr3_length(parent_exons, parent_stop, strand)

            delta: int | None = None
            delta_pct: float | None = None
            if parent_utr_len is not None:
                delta = iso_utr_len - parent_utr_len
                if parent_utr_len > 0:
                    delta_pct = delta / parent_utr_len * 100.0

            polya_motif = ""
            polya_dist: int | None = None
            if genome is not None and iso_utr_len > 0:
                seq = _extract_utr3_sequence(iso_exons, iso_stop, strand, genome)
                sig = polya_signal(seq)
                polya_motif = str(sig["motif"])
                if polya_motif:
                    polya_dist = int(sig["distance_from_3p_end"])

            rows.append(
                {
                    "transcript_id": tx_id,
                    "iso_utr3_length_nt": iso_utr_len,
                    "parent_utr3_length_nt": parent_utr_len,
                    "utr3_length_delta_nt": delta,
                    "utr3_length_delta_pct": delta_pct,
                    "polya_signal_motif": polya_motif,
                    "polya_signal_distance_nt": polya_dist,
                }
            )
        return pd.DataFrame(rows, columns=cols)
    finally:
        if genome is not None:
            genome.close()
