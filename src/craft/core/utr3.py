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


def _start_pos(cds_df: pd.DataFrame, strand: str) -> int:
    """Genomic 0-based position of the start codon's first base."""
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


def _iso_stop_pos(intervals: list[tuple], strand: str) -> int:
    if strand == "+":
        return max(end - 1 for _, _, end, _ in intervals)
    if strand == "-":
        return min(start for _, start, _, _ in intervals)
    raise ValueError(f"Unsupported strand: {strand!r}")


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
    parent_cds_by_tx = {
        tx: g for tx, g in ref_df[ref_df["Feature"] == "CDS"].groupby("transcript_id", sort=False)
    }

    genome = pysam.FastaFile(str(genome_fasta)) if genome_fasta is not None else None
    try:
        rows: list[dict] = []
        for tx_id, iso_exons in iso_exons_by_tx.items():
            res_row = resolved_by_tx.get(tx_id)
            strand = str(iso_strand.get(tx_id, "+"))
            parent_tx = iso_parent.get(tx_id, "")
            parent_cds = parent_cds_by_tx.get(parent_tx)
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
                iso_stop = _iso_stop_pos(res_row["resolved_cds_intervals"], strand)
                iso_utr3 = _utr3_length(iso_exons, iso_stop, strand)
                row["iso_utr3_length_nt"] = iso_utr3
                row["long_utr3_triggers_nmd"] = iso_utr3 > long_utr3_nt
                if parent_cds is not None and parent_exons is not None:
                    parent_stop = _parent_stop_pos(parent_cds, strand)
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
            if parent_cds is not None and iso_exons is not None and status != "resolution_failed":
                start_pos = _start_pos(parent_cds, strand)
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
