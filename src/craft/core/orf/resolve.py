"""Sequence-level ORF resolution.

Where :mod:`craft.core.orf.propagation` projects the parent CDS coordinates onto
the isoform purely geometrically, this module reconstructs the *real* ORF by
translating the isoform's own spliced sequence. For every isoform whose parent
start codon is observed, it walks the isoform transcript in 3-nt codons from the
projected start to the first in-frame stop. This naturally captures the
consequences the geometric path cannot see:

* a frameshift from an alternative splice site, which moves the stop;
* an exon skip that introduces a premature termination codon (PTC);
* an intron retained inside the CDS, which the parent-CDS ∩ iso-exon
  intersection silently drops (and which almost always carries a PTC).

The geometric columns are left untouched; everything here lands in new
``resolved_*`` columns so the v1.4 outputs stay reproducible.
"""

from enum import Enum
from pathlib import Path

import pandas as pd
import pyranges as pr
import pysam

from craft.core.orf.denovo import _transcript_sequence, _transcript_to_genomic_intervals
from craft.core.orf.propagation import (
    ORFOutcome,
    _start_codon_pos,
    _stop_codon_pos,
)

_STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})

# Default 50nt PTC window, shared with NMD. uORF NMD reuses the same threshold.
UORF_PTC_THRESHOLD_NT = 50

# Outcomes that carry an observed parent start codon and are therefore resolvable.
_RESOLVABLE_OUTCOMES = frozenset(
    {
        ORFOutcome.PROPAGATED_INTACT.value,
        ORFOutcome.DISRUPTED.value,
        ORFOutcome.STOP_NOT_OBSERVED.value,
        ORFOutcome.STOP_AT_ALT_POLYA.value,
    }
)


class ResolvedORFStatus(str, Enum):
    """Sequence-resolved ORF outcome."""

    INTACT = "intact"
    PTC_PREMATURE = "ptc_premature"
    PTC_INTRON_RETAINED = "ptc_intron_retained"
    CDS_EXTENSION = "cds_extension"
    NO_STOP_IN_READ = "no_stop_in_read"
    RESOLUTION_FAILED = "resolution_failed"


_COLUMNS = [
    "transcript_id",
    "resolved_orf_status",
    "resolved_stop_pos",
    "resolved_cds_bp",
    "resolved_aa_length",
    "resolved_cds_intervals",
    "ptc_introduced",
    "intron_retained_in_cds",
    "frame_consistent",
    "stop_in_transcript",
    "uorf_count",
    "uorf_triggers_nmd",
]


def _genomic_pos_to_tx_coord(pos: int, exons: pd.DataFrame, strand: str) -> int | None:
    """Transcript coordinate (0-based, 5'→3') of a genomic position, or None.

    For ``+`` strand the transcript runs in genomic order; for ``-`` strand it
    runs in reverse, so within an exon the 5'-most (smallest transcript) base is
    the highest genomic coordinate.
    """
    sorted_exons = exons.sort_values("Start").reset_index(drop=True)
    if strand == "+":
        walk = list(sorted_exons.itertuples(index=False))
    elif strand == "-":
        walk = list(sorted_exons.itertuples(index=False))[::-1]
    else:
        raise ValueError(f"Unsupported strand: {strand!r}")

    cum = 0
    for ex in walk:
        ex_start = int(ex.Start)
        ex_end = int(ex.End)
        length = ex_end - ex_start
        if ex_start <= pos < ex_end:
            if strand == "+":
                return cum + (pos - ex_start)
            return cum + (ex_end - 1 - pos)
        cum += length
    return None


def _walk_to_stop(seq: str, start_tx: int) -> tuple[int, bool]:
    """Walk codons from ``start_tx``; return (tx coord of stop-codon first base, found).

    If no in-frame stop appears before the transcript end, returns
    ``(len(seq), False)``.
    """
    n = len(seq)
    i = start_tx
    while i + 3 <= n:
        if seq[i : i + 3] in _STOP_CODONS:
            return i, True
        i += 3
    return n, False


def _introns_of(exons: pd.DataFrame) -> list[tuple[int, int]]:
    """Genomic introns (gaps between consecutive exons), sorted by start."""
    s = exons.sort_values("Start")
    starts = s["Start"].to_numpy()
    ends = s["End"].to_numpy()
    return [
        (int(ends[i]), int(starts[i + 1]))
        for i in range(len(starts) - 1)
        if starts[i + 1] > ends[i]
    ]


def _intron_retained_in_cds(
    iso_exons: pd.DataFrame, parent_introns: list[tuple[int, int]], cds_lo: int, cds_hi: int
) -> bool:
    """True if a parent intron inside the CDS span is engulfed by one iso exon.

    This is the precise intron-retention test: a CDS-region intron that the
    parent splices out but the isoform carries as continuous exonic sequence.
    Using engulfment (rather than a junction-count difference) avoids mislabeling
    an exon skip, which also lowers the junction count but is not retention.
    """
    if not parent_introns:
        return False
    iso = iso_exons[["Start", "End"]].to_numpy()
    for js, je in parent_introns:
        if js < cds_lo or je > cds_hi:
            continue
        for es, ee in iso:
            if int(es) <= js and je <= int(ee):
                return True
    return False


def _stop_pos_from_intervals(intervals: list[tuple], strand: str) -> int:
    """Genomic position of the last coding base of the resolved CDS."""
    if strand == "+":
        return max(end - 1 for _, _, end, _ in intervals)
    return min(start for _, start, _, _ in intervals)


def _scan_uorfs(
    seq: str, start_tx: int, exon_lengths: list[int], ptc_threshold: int
) -> tuple[int, bool]:
    """Count upstream ORFs in the 5'UTR and flag uORF-triggered NMD.

    A uORF is an ATG in ``seq[:start_tx]`` (the 5'UTR) with an in-frame stop codon
    before the main start. The NMD flag is the standard heuristic: a uORF stop
    sitting more than ``ptc_threshold`` mRNA-nt upstream of the transcript's last
    exon-exon junction is recognised as premature. Both are advisory.
    """
    utr5 = seq[:start_tx]
    stops: list[int] = []
    for p in range(0, len(utr5) - 2):
        if utr5[p : p + 3] != "ATG":
            continue
        j = p
        while j + 3 <= start_tx:
            if seq[j : j + 3] in _STOP_CODONS:
                stops.append(j)
                break
            j += 3
    if not stops:
        return 0, False

    last_junction_tx = sum(exon_lengths[:-1]) if len(exon_lengths) > 1 else None
    triggers = (
        last_junction_tx is not None
        and (last_junction_tx - min(stops)) > ptc_threshold
    )
    return len(stops), bool(triggers)


def resolve(
    classified: pr.PyRanges,
    propagated: pd.DataFrame,
    reference: pr.PyRanges,
    genome_fasta,
    ptc_threshold_nt: int = UORF_PTC_THRESHOLD_NT,
) -> pd.DataFrame:
    """Reconstruct the true ORF for each resolvable isoform.

    Args:
        classified: PyRanges of isoform exons (``transcript_id``, ``Strand``,
            ``Chromosome``).
        propagated: DataFrame from :func:`craft.core.orf.propagation.propagate`.
        reference: Reference PyRanges with a ``Feature`` column (exon / CDS rows).
        genome_fasta: Path to an indexed genome FASTA, or an open
            :class:`pysam.FastaFile`.

    Returns:
        DataFrame with one row per isoform and the ``resolved_*`` columns plus
        ``ptc_introduced``, ``intron_retained_in_cds``, ``frame_consistent``,
        ``stop_in_transcript``.
    """
    if propagated.empty or len(classified) == 0:
        return pd.DataFrame(columns=_COLUMNS)

    iso_df = classified.df
    iso_strand = iso_df.groupby("transcript_id")["Strand"].first().to_dict()
    iso_exons_by_tx = {tx: g for tx, g in iso_df.groupby("transcript_id", sort=False)}

    ref_df = reference.df
    parent_cds_by_tx = {
        tx: g for tx, g in ref_df[ref_df["Feature"] == "CDS"].groupby("transcript_id", sort=False)
    }
    parent_exons_by_tx = {
        tx: g for tx, g in ref_df[ref_df["Feature"] == "exon"].groupby("transcript_id", sort=False)
    }
    parent_introns_by_tx = {tx: _introns_of(g) for tx, g in parent_exons_by_tx.items()}

    own_genome = isinstance(genome_fasta, str | Path)
    genome = pysam.FastaFile(str(genome_fasta)) if own_genome else genome_fasta
    try:
        rows: list[dict] = []
        for _, prop_row in propagated.iterrows():
            tx_id = prop_row["transcript_id"]
            outcome = str(prop_row["orf_outcome"])
            parent_tx = prop_row["parent_tx_id"]

            if outcome not in _RESOLVABLE_OUTCOMES or parent_tx not in parent_cds_by_tx:
                rows.append(_failed_row(tx_id))
                continue

            strand = str(iso_strand[tx_id])
            iso_exons = iso_exons_by_tx[tx_id]
            parent_cds = parent_cds_by_tx[parent_tx]

            start_pos = _start_codon_pos(parent_cds, strand)
            start_tx = _genomic_pos_to_tx_coord(start_pos, iso_exons, strand)
            if start_tx is None:
                rows.append(_failed_row(tx_id))
                continue

            seq = _transcript_sequence(iso_exons, strand, genome)
            stop_tx, found = _walk_to_stop(seq, start_tx)

            sorted_exons = iso_exons.sort_values("Start")
            exon_lengths = (sorted_exons["End"] - sorted_exons["Start"]).astype(int).tolist()
            if strand == "-":
                exon_lengths = exon_lengths[::-1]
            uorf_count, uorf_triggers = _scan_uorfs(
                seq, start_tx, exon_lengths, ptc_threshold_nt
            )

            chrom = str(iso_exons["Chromosome"].iloc[0])
            intervals = _transcript_to_genomic_intervals(
                start_tx, stop_tx, iso_exons, strand, chrom
            )
            cds_bp = stop_tx - start_tx
            resolved_stop = (
                _stop_pos_from_intervals(intervals, strand) if intervals else None
            )

            # Where does the parent's stop fall in the iso's transcript?
            parent_stop_pos = _stop_codon_pos(parent_cds, strand)
            parent_stop_tx = _genomic_pos_to_tx_coord(parent_stop_pos, iso_exons, strand)

            # Intron retention anywhere in the parent CDS span (independent of
            # where the resolved stop lands, so an early PTC inside the retained
            # intron is still labelled as retention).
            cds_lo = int(parent_cds["Start"].min())
            cds_hi = int(parent_cds["End"].max())
            parent_introns = parent_introns_by_tx.get(parent_tx, [])
            ir = _intron_retained_in_cds(iso_exons, parent_introns, cds_lo, cds_hi)

            status, ptc, frame_ok = _classify(found, stop_tx, parent_stop_tx, ir)

            rows.append(
                {
                    "transcript_id": tx_id,
                    "resolved_orf_status": status.value,
                    "resolved_stop_pos": resolved_stop,
                    "resolved_cds_bp": cds_bp if found else 0,
                    "resolved_aa_length": cds_bp // 3 if found else 0,
                    "resolved_cds_intervals": intervals if found else [],
                    "ptc_introduced": ptc,
                    "intron_retained_in_cds": ir,
                    "frame_consistent": frame_ok,
                    "stop_in_transcript": found,
                    "uorf_count": uorf_count,
                    "uorf_triggers_nmd": uorf_triggers,
                }
            )
        return pd.DataFrame(rows, columns=_COLUMNS)
    finally:
        if own_genome:
            genome.close()


def _classify(
    found: bool,
    stop_tx: int,
    parent_stop_tx: int | None,
    ir: bool,
) -> tuple[ResolvedORFStatus, bool, bool]:
    """Map the walk result to (status, ptc_introduced, frame_consistent)."""
    if not found:
        return ResolvedORFStatus.NO_STOP_IN_READ, False, False

    # Transcript coord of the resolved last coding base.
    resolved_last_coding_tx = stop_tx - 1

    if parent_stop_tx is not None and resolved_last_coding_tx == parent_stop_tx:
        return ResolvedORFStatus.INTACT, False, True

    if parent_stop_tx is not None and resolved_last_coding_tx > parent_stop_tx:
        # Read past the parent stop: stop-loss / 3' extension.
        return ResolvedORFStatus.CDS_EXTENSION, False, False

    # Either the parent stop is not in the iso, or the resolved stop is upstream:
    # a premature stop relative to the parent ORF.
    if ir:
        return ResolvedORFStatus.PTC_INTRON_RETAINED, True, False
    return ResolvedORFStatus.PTC_PREMATURE, True, False


def _failed_row(tx_id: str) -> dict:
    return {
        "transcript_id": tx_id,
        "resolved_orf_status": ResolvedORFStatus.RESOLUTION_FAILED.value,
        "resolved_stop_pos": None,
        "resolved_cds_bp": 0,
        "resolved_aa_length": 0,
        "resolved_cds_intervals": [],
        "ptc_introduced": False,
        "intron_retained_in_cds": False,
        "frame_consistent": False,
        "stop_in_transcript": False,
        "uorf_count": 0,
        "uorf_triggers_nmd": False,
    }
