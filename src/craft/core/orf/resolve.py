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

Implementation note: this runs on *every* isoform, so the per-isoform work is
done on precomputed numpy arrays rather than per-row pandas operations. The
coordinate helpers mirror the (DataFrame-based) ones in
:mod:`craft.core.orf.denovo`; they are duplicated here in array form purely for
speed on full-genome inputs.
"""

from enum import Enum
from pathlib import Path

import numpy as np
import pandas as pd
import pyranges as pr
import pysam

from craft.core.orf.propagation import ORFOutcome

_STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})
_RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")

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
    START_RESCUED = "start_rescued"
    LEFT_CENSORED = "left_censored"
    RIGHT_CENSORED = "right_censored"
    RESOLUTION_FAILED = "resolution_failed"


_COLUMNS = [
    "transcript_id",
    "resolved_orf_status",
    "resolved_stop_pos",
    "resolved_start_pos",
    "resolved_stop_codon_pos",
    "resolved_cds_bp",
    "resolved_aa_length",
    "resolved_cds_intervals",
    "ptc_introduced",
    "intron_retained_in_cds",
    "frame_consistent",
    "stop_in_transcript",
    "uorf_count",
    "uorf_triggers_nmd",
    "orf_start_observed",
    "orf_stop_observed",
    "orf_censoring",
    "partial_cds_bp",
    "partial_cds_intervals",
    "alternative_start_inferred",
]


def _reverse_complement(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]


def _genomic_pos_to_tx_coord(
    pos: int, starts: np.ndarray, ends: np.ndarray, strand: str
) -> int | None:
    """Transcript coordinate (0-based, 5'→3') of a genomic position, or None.

    ``starts``/``ends`` are the isoform's exons sorted ascending by start. For
    ``+`` strand the transcript runs in genomic order; for ``-`` strand it runs
    in reverse, so within an exon the 5'-most base is the highest coordinate.
    """
    n = len(starts)
    cum = 0
    if strand == "+":
        rng = range(n)
    elif strand == "-":
        rng = range(n - 1, -1, -1)
    else:
        raise ValueError(f"Unsupported strand: {strand!r}")
    for i in rng:
        s = int(starts[i])
        e = int(ends[i])
        if s <= pos < e:
            return cum + (pos - s) if strand == "+" else cum + (e - 1 - pos)
        cum += e - s
    return None


def _spliced_sequence(
    chrom: str, starts: np.ndarray, ends: np.ndarray, strand: str, genome: pysam.FastaFile
) -> str:
    """Spliced transcript sequence (5'→3'); exons given sorted ascending by start."""
    parts = [genome.fetch(chrom, int(s), int(e)) for s, e in zip(starts, ends, strict=True)]
    seq = "".join(parts).upper()
    return _reverse_complement(seq) if strand == "-" else seq


def _tx_to_genomic_intervals(
    tx_start: int, tx_end: int, chrom: str, starts: np.ndarray, ends: np.ndarray, strand: str
) -> list[tuple[str, int, int, str]]:
    """Map [tx_start, tx_end) in transcript coords back to genomic intervals."""
    rng = range(len(starts)) if strand == "+" else range(len(starts) - 1, -1, -1)
    out: list[tuple[str, int, int, str]] = []
    cum = 0
    for i in rng:
        s = int(starts[i])
        e = int(ends[i])
        ex_tx_end = cum + (e - s)
        ov_s = max(cum, tx_start)
        ov_e = min(ex_tx_end, tx_end)
        if ov_e > ov_s:
            if strand == "+":
                out.append((chrom, s + (ov_s - cum), s + (ov_e - cum), strand))
            else:
                out.append((chrom, e - (ov_e - cum), e - (ov_s - cum), strand))
        cum = ex_tx_end
    out.sort(key=lambda x: x[1])
    return out


def _tx_coord_to_genomic_pos(
    coord: int, starts: np.ndarray, ends: np.ndarray, strand: str
) -> int | None:
    """Map one transcript coordinate to its genomic base position."""
    if coord < 0:
        return None
    rng = range(len(starts)) if strand == "+" else range(len(starts) - 1, -1, -1)
    cum = 0
    for i in rng:
        start = int(starts[i])
        end = int(ends[i])
        length = end - start
        if coord < cum + length:
            offset = coord - cum
            return start + offset if strand == "+" else end - 1 - offset
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


def _introns_of(starts: np.ndarray, ends: np.ndarray) -> list[tuple[int, int]]:
    """Genomic introns (gaps between consecutive exons), exons sorted ascending."""
    return [
        (int(ends[i]), int(starts[i + 1]))
        for i in range(len(starts) - 1)
        if starts[i + 1] > ends[i]
    ]


def _intron_retained_in_cds(
    starts: np.ndarray,
    ends: np.ndarray,
    parent_introns: list[tuple[int, int]],
    cds_lo: int,
    cds_hi: int,
) -> bool:
    """True if a parent intron inside the CDS span is engulfed by one iso exon.

    This is the precise intron-retention test: a CDS-region intron that the
    parent splices out but the isoform carries as continuous exonic sequence.
    Using engulfment (rather than a junction-count difference) avoids mislabeling
    an exon skip, which also lowers the junction count but is not retention.
    """
    if not parent_introns:
        return False
    for js, je in parent_introns:
        if js < cds_lo or je > cds_hi:
            continue
        for s, e in zip(starts, ends, strict=True):
            if int(s) <= js and je <= int(e):
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


def _parent_phase(g: int, cds_intervals: list[tuple[int, int]], strand: str) -> int | None:
    """Reading-frame phase (0/1/2) of genomic position ``g`` in the parent CDS.

    The phase is the count of parent CDS bases 5' of ``g`` in transcript order,
    modulo 3. Returns None if ``g`` is not inside the parent CDS.
    """
    ivs = sorted((int(s), int(e)) for s, e in cds_intervals)
    cum = 0
    if strand == "+":
        for s, e in ivs:
            if g < s:
                return None
            if g < e:
                return (cum + (g - s)) % 3
            cum += e - s
        return None
    for s, e in reversed(ivs):
        if g >= e:
            return None
        if g >= s:
            return (cum + (e - 1 - g)) % 3
        cum += e - s
    return None


def _first_inframe_atg(seq: str, frame_offset: int) -> int | None:
    """Transcript coord of the first ATG on the ``frame_offset`` codon frame, or None."""
    i = frame_offset % 3
    n = len(seq)
    while i + 3 <= n:
        if seq[i : i + 3] == "ATG":
            return i
        i += 3
    return None


def _rescue_start_lost(
    tx_id: str,
    chrom: str,
    strand: str,
    starts: np.ndarray,
    ends: np.ndarray,
    seq: str,
    prop_intervals: list[tuple] | None,
    parent_cds_intervals: list[tuple[int, int]],
    parent_introns: list[tuple[int, int]],
    cds_span: tuple[int, int],
    ptc_threshold_nt: int,
) -> dict | None:
    """Frame-aware ORF rescue for a 5'-truncated isoform whose start codon is gone.

    Keeps the parent CDS reading frame (anchored on the 5'-most parent-CDS base
    still present in the isoform) and takes the first in-frame ATG from the
    isoform's 5' end, translating to the first in-frame stop. This extends the
    reference-trust to the ``start_lost`` case instead of surrendering it to a
    blind de-novo longest-ORF search. Returns a resolved row dict, or None when no
    in-frame ATG + stop exists (the caller then falls back to de novo).
    """
    if not prop_intervals:
        return None
    if strand == "+":
        anchor_g = min(int(s) for _, s, _, _ in prop_intervals)
    else:
        anchor_g = max(int(e) for _, _, e, _ in prop_intervals) - 1
    phase = _parent_phase(anchor_g, parent_cds_intervals, strand)
    anchor_t = _genomic_pos_to_tx_coord(anchor_g, starts, ends, strand)
    if phase is None or anchor_t is None:
        return None

    frame_offset = (anchor_t - phase) % 3
    rescued_start = _first_inframe_atg(seq, frame_offset)
    if rescued_start is None:
        return None
    stop_tx, found = _walk_to_stop(seq, rescued_start)
    if not found:
        return None

    intervals = _tx_to_genomic_intervals(rescued_start, stop_tx, chrom, starts, ends, strand)
    cds_bp = stop_tx - rescued_start
    resolved_stop = _stop_pos_from_intervals(intervals, strand) if intervals else None
    cds_lo, cds_hi = cds_span
    ir = _intron_retained_in_cds(starts, ends, parent_introns, cds_lo, cds_hi)
    lengths = (ends - starts).astype(int)
    exon_lengths = lengths.tolist() if strand == "+" else lengths[::-1].tolist()
    uorf_count, uorf_triggers = _scan_uorfs(seq, rescued_start, exon_lengths, ptc_threshold_nt)
    return {
        "transcript_id": tx_id,
        "resolved_orf_status": ResolvedORFStatus.START_RESCUED.value,
        "resolved_stop_pos": resolved_stop,
        "resolved_start_pos": _tx_coord_to_genomic_pos(rescued_start, starts, ends, strand),
        "resolved_stop_codon_pos": _tx_coord_to_genomic_pos(stop_tx, starts, ends, strand),
        "resolved_cds_bp": cds_bp,
        "resolved_aa_length": cds_bp // 3,
        "resolved_cds_intervals": intervals,
        "ptc_introduced": False,
        "intron_retained_in_cds": ir,
        "frame_consistent": False,
        "stop_in_transcript": True,
        "uorf_count": uorf_count,
        "uorf_triggers_nmd": uorf_triggers,
        "orf_start_observed": False,
        "orf_stop_observed": True,
        "orf_censoring": "left",
        "partial_cds_bp": cds_bp,
        "partial_cds_intervals": intervals,
        "alternative_start_inferred": True,
    }


def _left_censored_row(
    tx_id: str,
    chrom: str,
    strand: str,
    starts: np.ndarray,
    ends: np.ndarray,
    prop_intervals: list[tuple] | None,
) -> dict:
    """Represent an observed CDS fragment without inventing a translation start."""
    intervals = list(prop_intervals or [])
    bp = sum(int(end) - int(start) for _, start, end, _ in intervals)
    row = _failed_row(tx_id)
    row.update(
        {
            "resolved_orf_status": ResolvedORFStatus.LEFT_CENSORED.value,
            "resolved_start_pos": None,
            "resolved_stop_codon_pos": None,
            "resolved_cds_intervals": [],
            "orf_start_observed": False,
            "orf_stop_observed": False,
            "orf_censoring": "left",
            "partial_cds_bp": bp,
            "partial_cds_intervals": intervals,
        }
    )
    return row


def resolve(
    classified: pr.PyRanges,
    propagated: pd.DataFrame,
    reference: pr.PyRanges,
    genome_fasta,
    ptc_threshold_nt: int = UORF_PTC_THRESHOLD_NT,
    allow_start_rescue: bool = False,
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

    # Per-isoform exon arrays, sorted ascending by start, built once.
    iso_df = classified.df.sort_values(["transcript_id", "Start"], kind="stable")
    iso_by_tx: dict[str, tuple[str, str, np.ndarray, np.ndarray]] = {}
    for tx, g in iso_df.groupby("transcript_id", sort=False):
        iso_by_tx[tx] = (
            str(g["Chromosome"].iat[0]),
            str(g["Strand"].iat[0]),
            g["Start"].to_numpy(),
            g["End"].to_numpy(),
        )

    # Parent start/stop positions, CDS span, and introns, all precomputed once.
    ref_df = reference.df
    parent_start: dict[str, int] = {}
    parent_stop: dict[str, int] = {}
    parent_cds_span: dict[str, tuple[int, int]] = {}
    parent_cds_intervals: dict[str, list[tuple[int, int]]] = {}
    for tx, g in ref_df[ref_df["Feature"] == "CDS"].groupby("transcript_id", sort=False):
        strand = str(g["Strand"].iat[0])
        smin = int(g["Start"].min())
        emax = int(g["End"].max())
        parent_cds_span[tx] = (smin, emax)
        parent_cds_intervals[tx] = [
            (int(s), int(e)) for s, e in zip(g["Start"], g["End"], strict=True)
        ]
        if strand == "+":
            parent_start[tx], parent_stop[tx] = smin, emax - 1
        else:
            parent_start[tx], parent_stop[tx] = emax - 1, smin
    parent_introns: dict[str, list[tuple[int, int]]] = {}
    for tx, g in ref_df[ref_df["Feature"] == "exon"].groupby("transcript_id", sort=False):
        s = g.sort_values("Start")
        parent_introns[tx] = _introns_of(s["Start"].to_numpy(), s["End"].to_numpy())

    tx_ids = propagated["transcript_id"].to_numpy()
    outcomes = propagated["orf_outcome"].astype(str).to_numpy()
    parents = propagated["parent_tx_id"].astype(str).to_numpy()
    prop_iv = dict(
        zip(propagated["transcript_id"], propagated["propagated_cds_intervals"], strict=True)
    )

    own_genome = isinstance(genome_fasta, str | Path)
    genome = pysam.FastaFile(str(genome_fasta)) if own_genome else genome_fasta
    try:
        rows: list[dict] = []
        for tx_id, outcome, parent_tx in zip(tx_ids, outcomes, parents, strict=True):
            if parent_tx not in parent_start:
                rows.append(_failed_row(tx_id))
                continue

            chrom, strand, starts, ends = iso_by_tx[tx_id]
            seq = _spliced_sequence(chrom, starts, ends, strand, genome)

            if outcome == ORFOutcome.START_LOST.value:
                rescued = None
                if allow_start_rescue:
                    rescued = _rescue_start_lost(
                        tx_id, chrom, strand, starts, ends, seq,
                        prop_iv.get(tx_id), parent_cds_intervals.get(parent_tx, []),
                        parent_introns.get(parent_tx, []), parent_cds_span[parent_tx],
                        ptc_threshold_nt,
                    )
                rows.append(
                    rescued if rescued is not None else _left_censored_row(
                        tx_id, chrom, strand, starts, ends, prop_iv.get(tx_id)
                    )
                )
                continue

            if outcome not in _RESOLVABLE_OUTCOMES:
                rows.append(_failed_row(tx_id))
                continue
            start_tx = _genomic_pos_to_tx_coord(parent_start[parent_tx], starts, ends, strand)
            if start_tx is None:
                rows.append(_failed_row(tx_id))
                continue

            stop_tx, found = _walk_to_stop(seq, start_tx)

            lengths = (ends - starts).astype(int)
            exon_lengths = lengths.tolist() if strand == "+" else lengths[::-1].tolist()
            uorf_count, uorf_triggers = _scan_uorfs(
                seq, start_tx, exon_lengths, ptc_threshold_nt
            )

            intervals = _tx_to_genomic_intervals(start_tx, stop_tx, chrom, starts, ends, strand)
            cds_bp = stop_tx - start_tx
            resolved_stop = (
                _stop_pos_from_intervals(intervals, strand) if intervals else None
            )

            parent_stop_tx = _genomic_pos_to_tx_coord(
                parent_stop[parent_tx], starts, ends, strand
            )
            cds_lo, cds_hi = parent_cds_span[parent_tx]
            ir = _intron_retained_in_cds(
                starts, ends, parent_introns.get(parent_tx, []), cds_lo, cds_hi
            )

            status, ptc, frame_ok = _classify(found, stop_tx, parent_stop_tx, ir)
            if not found:
                status = ResolvedORFStatus.RIGHT_CENSORED
            partial_intervals = (
                _tx_to_genomic_intervals(start_tx, len(seq), chrom, starts, ends, strand)
                if not found else intervals
            )
            partial_bp = len(seq) - start_tx if not found else cds_bp

            rows.append(
                {
                    "transcript_id": tx_id,
                    "resolved_orf_status": status.value,
                    "resolved_stop_pos": resolved_stop,
                    "resolved_start_pos": _tx_coord_to_genomic_pos(start_tx, starts, ends, strand),
                    "resolved_stop_codon_pos": (
                        _tx_coord_to_genomic_pos(stop_tx, starts, ends, strand) if found else None
                    ),
                    "resolved_cds_bp": cds_bp if found else 0,
                    "resolved_aa_length": cds_bp // 3 if found else 0,
                    "resolved_cds_intervals": intervals if found else [],
                    "ptc_introduced": ptc,
                    "intron_retained_in_cds": ir,
                    "frame_consistent": frame_ok,
                    "stop_in_transcript": found,
                    "uorf_count": uorf_count,
                    "uorf_triggers_nmd": uorf_triggers,
                    "orf_start_observed": True,
                    "orf_stop_observed": found,
                    "orf_censoring": "none" if found else "right",
                    "partial_cds_bp": partial_bp,
                    "partial_cds_intervals": partial_intervals,
                    "alternative_start_inferred": False,
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
        "resolved_start_pos": None,
        "resolved_stop_codon_pos": None,
        "resolved_cds_bp": 0,
        "resolved_aa_length": 0,
        "resolved_cds_intervals": [],
        "ptc_introduced": False,
        "intron_retained_in_cds": False,
        "frame_consistent": False,
        "stop_in_transcript": False,
        "uorf_count": 0,
        "uorf_triggers_nmd": False,
        "orf_start_observed": False,
        "orf_stop_observed": False,
        "orf_censoring": "unknown",
        "partial_cds_bp": 0,
        "partial_cds_intervals": [],
        "alternative_start_inferred": False,
    }
