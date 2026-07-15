"""3' and 5' UTR feature consequences: length deltas vs parent + poly(A) motif scan.

UTR lengths are measured from the **resolved** ORF (the real in-frame start and
stop), so the columns are single, canonical names (no geometric/resolved split).
Internal priming detection stays in `tecap`; this module only reports structural
and sequence features of the UTRs.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pyranges as pr
import pysam

from craft.core.intervals import (
    genomic_position_at_transcript_coordinate,
    spliced_length,
    transcript_coordinate,
)

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

# Resolved-ORF statuses (from craft.core.orf.resolve) that carry a real stop.
_RESOLVED_WITH_STOP = frozenset(
    {"intact", "ptc_premature", "ptc_intron_retained", "cds_extension", "start_rescued"}
)

LONG_UTR3_NT = 1000


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


def _utr3_length(exons: pd.DataFrame, stop_codon_pos: int, strand: str) -> int:
    """Total mRNA bp after the complete three-base stop codon."""
    stop_tx = transcript_coordinate(exons, stop_codon_pos, strand)
    if stop_tx is None:
        return 0
    return max(spliced_length(exons) - (stop_tx + 3), 0)


def _parent_stop_codon_pos(parent_df: pd.DataFrame, strand: str) -> int | None:
    stop = parent_df[parent_df["Feature"] == "stop_codon"]
    if not stop.empty:
        if strand == "+":
            return int(stop["Start"].min())
        return int(stop["End"].max()) - 1
    cds_df = parent_df[parent_df["Feature"] == "CDS"]
    exons = parent_df[parent_df["Feature"] == "exon"]
    last_sense = (
        int(cds_df["End"].max()) - 1 if strand == "+" else int(cds_df["Start"].min())
    )
    last_sense_tx = transcript_coordinate(exons, last_sense, strand)
    # None when the parent CDS ends at the transcript 3' end and no stop_codon is
    # annotated (e.g. GENCODE cds_end_NF): the stop base falls off the exon chain.
    # Return None so the caller skips the parent-relative 3'UTR fields for this
    # isoform rather than crashing the whole run.
    return (
        genomic_position_at_transcript_coordinate(exons, last_sense_tx + 1, strand)
        if last_sense_tx is not None else None
    )


def _start_pos(parent_df: pd.DataFrame, strand: str) -> int:
    """Genomic 0-based position of the start codon's first base."""
    starts = parent_df[parent_df["Feature"] == "start_codon"]
    if not starts.empty:
        if strand == "+":
            return int(starts["Start"].min())
        return int(starts["End"].max()) - 1
    cds_df = parent_df[parent_df["Feature"] == "CDS"]
    if strand == "+":
        return int(cds_df["Start"].min())
    if strand == "-":
        return int(cds_df["End"].max()) - 1
    raise ValueError(f"Unsupported strand: {strand!r}")


def _utr5_length(exons: pd.DataFrame, start_pos: int, strand: str) -> int:
    """Total mRNA bp upstream of the start codon (in transcript order)."""
    starts = exons["Start"].to_numpy()
    ends = exons["End"].to_numpy()
    if strand == "+":
        contrib = np.maximum(0, np.minimum(ends, start_pos) - starts)
    elif strand == "-":
        contrib = np.maximum(0, ends - np.maximum(starts, start_pos + 1))
    else:
        raise ValueError(f"Unsupported strand: {strand!r}")
    return int(contrib.sum())


def _stop_codon_from_cds_intervals(
    intervals: list[tuple], exons: pd.DataFrame, strand: str
) -> int:
    """Infer stop-codon first base from CDS intervals that exclude the stop."""
    last_sense = (
        max(int(end) for _, _, end, _ in intervals) - 1
        if strand == "+"
        else min(int(start) for _, start, _, _ in intervals)
    )
    last_sense_tx = transcript_coordinate(exons, last_sense, strand)
    stop_pos = (
        genomic_position_at_transcript_coordinate(exons, last_sense_tx + 1, strand)
        if last_sense_tx is not None else None
    )
    if stop_pos is None:
        raise ValueError("Could not place the inferred stop codon on the exon chain")
    return stop_pos


def polya_near_3prime_end(
    exons: pd.DataFrame,
    strand: str,
    genome: pysam.FastaFile | Path,
    window: int = 50,
) -> dict:
    """Scan the last ``window`` bp of the isoform (transcript orientation) for a poly(A) signal.

    For oligo-dT primed long-read data, the iso's 3' end *is* the polyadenylation
    site (the polyA tail was the priming substrate). A canonical poly(A) signal
    motif sitting within ~10-30 nt upstream of the cleavage site is strong
    biological evidence that the iso's 3' end reflects alternative
    polyadenylation rather than technical truncation.

    Args:
        exons: PyRanges-style DataFrame of the isoform's exons (single
            transcript_id).
        strand: "+" or "-".
        genome: either a path to an indexed FASTA, or an already-open
            :class:`pysam.FastaFile` (preferred when calling per-isoform in a
            loop to avoid repeated open/close).
        window: nt window upstream of the iso's 3' end to scan. Default 50.

    Returns:
        Dict with keys ``motif`` (empty string if no signal), ``distance_from_3p_end``
        (-1 if no signal), ``found`` (bool).
    """
    empty = {"motif": "", "distance_from_3p_end": -1, "found": False}
    if len(exons) == 0:
        return empty

    if isinstance(genome, str | Path):
        with pysam.FastaFile(str(genome)) as fh:
            return polya_near_3prime_end(exons, strand, fh, window)

    chrom = str(exons["Chromosome"].iloc[0])
    sorted_exons = exons.sort_values("Start").reset_index(drop=True)
    if strand == "+":
        last_exon = sorted_exons.iloc[-1]
        ex_start = int(last_exon["Start"])
        ex_end = int(last_exon["End"])
        fetch_start = max(ex_start, ex_end - window)
        fetch_end = ex_end
    elif strand == "-":
        last_exon = sorted_exons.iloc[0]
        ex_start = int(last_exon["Start"])
        ex_end = int(last_exon["End"])
        fetch_start = ex_start
        fetch_end = min(ex_end, ex_start + window)
    else:
        raise ValueError(f"Unsupported strand: {strand!r}")

    if fetch_end <= fetch_start:
        return empty

    seq = genome.fetch(chrom, fetch_start, fetch_end).upper()
    if strand == "-":
        seq = _reverse_complement(seq)
    sig = polya_signal(seq)
    return {**sig, "found": bool(sig["motif"])}


def _extract_utr3_sequence(
    exons: pd.DataFrame,
    stop_codon_pos: int,
    strand: str,
    genome: pysam.FastaFile,
) -> str:
    """Extract the 3' UTR sequence in transcript orientation (5' to 3')."""
    chrom = str(exons["Chromosome"].iloc[0])
    sorted_exons = exons.sort_values("Start", kind="stable")
    sequence = "".join(
        genome.fetch(chrom, int(exon.Start), int(exon.End))
        for exon in sorted_exons.itertuples(index=False)
    ).upper()
    if strand == "-":
        sequence = _reverse_complement(sequence)
    stop_tx = transcript_coordinate(exons, stop_codon_pos, strand)
    return sequence[stop_tx + 3 :] if stop_tx is not None else ""


_COLUMNS = [
    "transcript_id",
    "iso_utr3_length_nt",
    "parent_utr3_length_nt",
    "utr3_length_delta_nt",
    "utr3_length_delta_pct",
    "iso_utr5_length_nt",
    "parent_utr5_length_nt",
    "utr5_length_delta_nt",
    "utr5_length_delta_pct",
    "long_utr3_triggers_nmd",
    "polya_signal_motif",
    "polya_signal_distance_nt",
]


def annotate(
    classified: pr.PyRanges,
    resolved: pd.DataFrame,
    reference: pr.PyRanges,
    genome_fasta: Path | None = None,
    long_utr3_nt: int = LONG_UTR3_NT,
) -> pd.DataFrame:
    """UTR consequences per isoform, measured from the resolved ORF.

    3'UTR length/delta are measured from the resolved (real in-frame) stop; 5'UTR
    length/delta from the start codon. When ``genome_fasta`` is given, the iso's
    3'UTR is scanned for the strongest canonical poly(A) signal.

    Args:
        classified: PyRanges of isoform exons (``transcript_id``, ``Strand``,
            ``parent_tx_id``).
        resolved: DataFrame from :func:`craft.core.orf.resolve.resolve`.
        reference: Reference PyRanges with a ``Feature`` column (exon / CDS rows).
        genome_fasta: Optional indexed genome FASTA; if omitted the poly(A) scan
            is skipped.
        long_utr3_nt: 3'UTR length above which ``long_utr3_triggers_nmd`` is set.

    Returns:
        DataFrame with ``transcript_id`` plus ``iso_utr3_length_nt``,
        ``parent_utr3_length_nt``, ``utr3_length_delta_nt``, ``utr3_length_delta_pct``,
        ``iso_utr5_length_nt``, ``parent_utr5_length_nt``, ``utr5_length_delta_nt``,
        ``utr5_length_delta_pct``, ``long_utr3_triggers_nmd``, ``polya_signal_motif``,
        ``polya_signal_distance_nt``.
    """
    if len(classified) == 0:
        return pd.DataFrame(columns=_COLUMNS)

    iso_df = classified.df
    iso_strand = iso_df.groupby("transcript_id")["Strand"].first().to_dict()
    iso_parent = iso_df.groupby("transcript_id")["parent_tx_id"].first().to_dict()
    iso_exons_by_tx = {tx: g for tx, g in iso_df.groupby("transcript_id", sort=False)}
    resolved_by_tx: dict = {}
    if resolved is not None and not resolved.empty:
        resolved_by_tx = {r["transcript_id"]: r for _, r in resolved.iterrows()}

    ref_df = reference.df
    parent_exons_by_tx = {
        tx: g for tx, g in ref_df[ref_df["Feature"] == "exon"].groupby("transcript_id", sort=False)
    }
    parent_records_by_tx = {tx: g for tx, g in ref_df.groupby("transcript_id", sort=False)}

    genome = pysam.FastaFile(str(genome_fasta)) if genome_fasta is not None else None
    try:
        rows: list[dict] = []
        for tx_id, iso_exons in iso_exons_by_tx.items():
            res_row = resolved_by_tx.get(tx_id)
            strand = str(iso_strand.get(tx_id, "+"))
            parent_tx = iso_parent.get(tx_id, "")
            parent_records = parent_records_by_tx.get(parent_tx)
            parent_exons = parent_exons_by_tx.get(parent_tx)
            status = (
                str(res_row["resolved_orf_status"]) if res_row is not None else "resolution_failed"
            )

            row = dict.fromkeys(_COLUMNS)
            row["transcript_id"] = tx_id
            row["long_utr3_triggers_nmd"] = False
            row["polya_signal_motif"] = ""

            has_stop = (
                res_row is not None
                and status in _RESOLVED_WITH_STOP
                and bool(res_row["stop_in_transcript"])
                and res_row["resolved_cds_intervals"]
            )
            if has_stop and iso_exons is not None:
                iso_stop = (
                    int(res_row["resolved_stop_codon_pos"])
                    if pd.notna(res_row.get("resolved_stop_codon_pos"))
                    else _stop_codon_from_cds_intervals(
                        res_row["resolved_cds_intervals"], iso_exons, strand
                    )
                )
                iso_utr3 = _utr3_length(iso_exons, iso_stop, strand)
                row["iso_utr3_length_nt"] = iso_utr3
                row["long_utr3_triggers_nmd"] = iso_utr3 > long_utr3_nt
                if parent_records is not None and parent_exons is not None:
                    parent_stop = _parent_stop_codon_pos(parent_records, strand)
                    if parent_stop is not None:
                        parent_utr3 = _utr3_length(parent_exons, parent_stop, strand)
                        row["parent_utr3_length_nt"] = parent_utr3
                        row["utr3_length_delta_nt"] = iso_utr3 - parent_utr3
                        if parent_utr3 > 0:
                            row["utr3_length_delta_pct"] = (
                                (iso_utr3 - parent_utr3) / parent_utr3 * 100.0
                            )
                if genome is not None and iso_utr3 > 0:
                    seq = _extract_utr3_sequence(iso_exons, iso_stop, strand, genome)
                    sig = polya_signal(seq)
                    motif = str(sig["motif"])
                    row["polya_signal_motif"] = motif
                    if motif:
                        row["polya_signal_distance_nt"] = int(sig["distance_from_3p_end"])

            # 5'UTR: only when the start codon is observed (resolution succeeded).
            if (
                parent_records is not None
                and iso_exons is not None
                and status not in {"resolution_failed", "left_censored"}
            ):
                start_pos = (
                    int(res_row["resolved_start_pos"])
                    if res_row is not None and pd.notna(res_row.get("resolved_start_pos"))
                    else _start_pos(parent_records, strand)
                )
                iso_utr5 = _utr5_length(iso_exons, start_pos, strand)
                row["iso_utr5_length_nt"] = iso_utr5
                if parent_exons is not None:
                    parent_utr5 = _utr5_length(parent_exons, start_pos, strand)
                    row["parent_utr5_length_nt"] = parent_utr5
                    row["utr5_length_delta_nt"] = iso_utr5 - parent_utr5
                    if parent_utr5 > 0:
                        row["utr5_length_delta_pct"] = (
                            (iso_utr5 - parent_utr5) / parent_utr5 * 100.0
                        )
            rows.append(row)
        return pd.DataFrame(rows, columns=_COLUMNS)
    finally:
        if genome is not None:
            genome.close()
