"""Nonsense-mediated decay (NMD) susceptibility prediction.

NMD is predicted from the **sequence-resolved** ORF stop (the real in-frame stop
found by translating the isoform's own spliced CDS). For orphan isoforms that
have no reference-anchored ORF, the call falls back to the de-novo ORF's stop.
There is a single NMD answer per isoform; ``nmd_basis`` records which ORF it came
from (``resolved`` / ``denovo`` / ``none``).

Rule cascade (each is a sufficient condition for escape, evaluated in order):
1. Stop codon in the transcript's last exon (or single-exon transcript).
2. Stop codon within ``PTC_THRESHOLD_NT`` (default 50) mRNA-nt of the last junction.
3. Start-proximal: CDS shorter than ``START_PROXIMAL_NT`` (default 150) bp.
4. Last exon longer than ``LONG_LAST_EXON_NT`` (default 400) bp.
Otherwise the transcript is NMD-sensitive (50nt rule violated).

Confidence is ``high`` for a resolved intact ORF, ``medium`` for a resolved but
altered ORF (PTC / intron-retained / extension), ``low`` for a de-novo call (the
stop is from a predicted ORF, not a reference), and ``none`` when not applicable.
"""

from enum import Enum

import pandas as pd
import pyranges as pr

from craft.core.orf.confidence import ORFConfidence

PTC_THRESHOLD_NT = 50
START_PROXIMAL_NT = 150
LONG_LAST_EXON_NT = 400

# Resolved-ORF statuses (from craft.core.orf.resolve) that carry a real stop.
_RESOLVED_WITH_STOP = frozenset(
    {"intact", "ptc_premature", "ptc_intron_retained", "cds_extension"}
)

_COLUMNS = [
    "transcript_id",
    "nmd_status",
    "nmd_rule",
    "nmd_confidence",
    "nmd_basis",
    "stop_to_last_junction_nt",
    "last_exon_length_nt",
]


class NMDStatus(str, Enum):
    """NMD susceptibility class."""

    SENSITIVE = "sensitive"
    ESCAPED = "escaped"
    NOT_APPLICABLE = "not_applicable"


def _stop_codon_genome(intervals: list[tuple], strand: str) -> int:
    """Genomic 0-based position of the last CDS base in transcript order."""
    if strand == "+":
        return max(end - 1 for _, _, end, _ in intervals)
    if strand == "-":
        return min(start for _, start, _, _ in intervals)
    raise ValueError(f"Unsupported strand: {strand!r}")


def _distance_stop_to_last_junction(
    stop_pos: int, exons: pd.DataFrame, strand: str
) -> tuple[int, bool]:
    """Distance (mRNA bp) from the stop codon to the last exon-exon junction.

    Returns ``(distance, is_in_last_exon)``. For single-exon transcripts and for
    stops in the transcript's last exon, returns ``(0, True)``.
    """
    sorted_exons = exons.sort_values("Start").reset_index(drop=True)
    n_exons = len(sorted_exons)
    if n_exons <= 1:
        return 0, True

    mask = (sorted_exons["Start"] <= stop_pos) & (stop_pos < sorted_exons["End"])
    if not mask.any():
        return -1, False
    stop_exon_idx = int(sorted_exons.index[mask][0])

    if strand == "+":
        last_exon_idx = n_exons - 1
    elif strand == "-":
        last_exon_idx = 0
    else:
        raise ValueError(f"Unsupported strand: {strand!r}")

    if stop_exon_idx == last_exon_idx:
        return 0, True

    stop_exon = sorted_exons.iloc[stop_exon_idx]
    if strand == "+":
        distance = int(stop_exon["End"]) - stop_pos
        intermediate = range(stop_exon_idx + 1, last_exon_idx)
    else:
        distance = stop_pos - int(stop_exon["Start"])
        intermediate = range(1, stop_exon_idx)
    for i in intermediate:
        ex = sorted_exons.iloc[i]
        distance += int(ex["End"]) - int(ex["Start"])
    return distance, False


def _last_exon_length(exons: pd.DataFrame, strand: str) -> int:
    sorted_exons = exons.sort_values("Start").reset_index(drop=True)
    if strand == "+":
        last = sorted_exons.iloc[-1]
    elif strand == "-":
        last = sorted_exons.iloc[0]
    else:
        raise ValueError(f"Unsupported strand: {strand!r}")
    return int(last["End"]) - int(last["Start"])


def _cascade(
    distance: int,
    in_last: bool,
    cds_bp: int,
    last_exon_len: int,
    ptc_threshold_nt: int,
    start_proximal_nt: int,
    long_last_exon_nt: int,
) -> tuple[NMDStatus, str]:
    if in_last:
        return NMDStatus.ESCAPED, "stop_in_last_exon"
    if distance <= ptc_threshold_nt:
        return NMDStatus.ESCAPED, "within_50nt_of_last_junction"
    if cds_bp < start_proximal_nt:
        return NMDStatus.ESCAPED, "start_proximal"
    if last_exon_len > long_last_exon_nt:
        return NMDStatus.ESCAPED, "long_last_exon"
    return NMDStatus.SENSITIVE, "ptc_50nt_rule"


def _not_applicable(tx_id: str) -> dict:
    return {
        "transcript_id": tx_id,
        "nmd_status": NMDStatus.NOT_APPLICABLE.value,
        "nmd_rule": "",
        "nmd_confidence": ORFConfidence.NONE.value,
        "nmd_basis": "none",
        "stop_to_last_junction_nt": None,
        "last_exon_length_nt": None,
    }


def predict(
    classified: pr.PyRanges,
    resolved: pd.DataFrame,
    denovo: pd.DataFrame | None = None,
    ptc_threshold_nt: int = PTC_THRESHOLD_NT,
    start_proximal_nt: int = START_PROXIMAL_NT,
    long_last_exon_nt: int = LONG_LAST_EXON_NT,
) -> pd.DataFrame:
    """Predict per-isoform NMD from the resolved ORF (de-novo fallback for orphans).

    Args:
        classified: PyRanges of isoform exons (``transcript_id``, ``Strand``).
        resolved: DataFrame from :func:`craft.core.orf.resolve.resolve`.
        denovo: DataFrame from :func:`craft.core.orf.denovo.predict` (optional;
            supplies the orphan fallback stop).
        ptc_threshold_nt: 50nt PTC rule threshold (mRNA distance to last junction).
        start_proximal_nt: start-proximal escape window.
        long_last_exon_nt: long-last-exon escape threshold.

    Returns:
        DataFrame with ``transcript_id``, ``nmd_status``, ``nmd_rule``,
        ``nmd_confidence``, ``nmd_basis``, ``stop_to_last_junction_nt``,
        ``last_exon_length_nt``.
    """
    if len(classified) == 0:
        return pd.DataFrame(columns=_COLUMNS)

    iso_df = classified.df
    iso_strand = iso_df.groupby("transcript_id")["Strand"].first().to_dict()
    iso_exons_by_tx = {tx: g for tx, g in iso_df.groupby("transcript_id", sort=False)}

    resolved_by_tx: dict = {}
    if resolved is not None and not resolved.empty:
        resolved_by_tx = {r["transcript_id"]: r for _, r in resolved.iterrows()}
    denovo_by_tx: dict = {}
    if denovo is not None and not denovo.empty:
        denovo_by_tx = {r["transcript_id"]: r for _, r in denovo.iterrows()}

    rows: list[dict] = []
    for tx_id in iso_exons_by_tx:
        intervals = None
        cds_bp = 0
        basis = "none"
        confidence = ORFConfidence.NONE

        res = resolved_by_tx.get(tx_id)
        if (
            res is not None
            and str(res["resolved_orf_status"]) in _RESOLVED_WITH_STOP
            and bool(res["stop_in_transcript"])
            and res["resolved_cds_intervals"]
        ):
            intervals = res["resolved_cds_intervals"]
            cds_bp = int(res["resolved_cds_bp"])
            basis = "resolved"
            confidence = (
                ORFConfidence.HIGH
                if str(res["resolved_orf_status"]) == "intact"
                else ORFConfidence.MEDIUM
            )
        else:
            dn = denovo_by_tx.get(tx_id)
            if (
                dn is not None
                and bool(dn["denovo_orf_found"])
                and isinstance(dn["denovo_cds_intervals"], list)
                and dn["denovo_cds_intervals"]
            ):
                intervals = dn["denovo_cds_intervals"]
                cds_bp = int(dn["denovo_cds_bp"])
                basis = "denovo"
                confidence = ORFConfidence.LOW

        if intervals is None:
            rows.append(_not_applicable(tx_id))
            continue

        strand = str(iso_strand[tx_id])
        iso_exons = iso_exons_by_tx[tx_id]
        stop_pos = _stop_codon_genome(intervals, strand)
        distance, in_last = _distance_stop_to_last_junction(stop_pos, iso_exons, strand)
        last_exon_len = _last_exon_length(iso_exons, strand)
        status, rule = _cascade(
            distance, in_last, cds_bp, last_exon_len,
            ptc_threshold_nt, start_proximal_nt, long_last_exon_nt,
        )
        rows.append(
            {
                "transcript_id": tx_id,
                "nmd_status": status.value,
                "nmd_rule": rule,
                "nmd_confidence": confidence.value,
                "nmd_basis": basis,
                "stop_to_last_junction_nt": int(distance),
                "last_exon_length_nt": int(last_exon_len),
            }
        )
    return pd.DataFrame(rows, columns=_COLUMNS)
