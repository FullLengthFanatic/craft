"""Nonsense-mediated decay (NMD) susceptibility prediction.

Rule cascade (each is a sufficient condition for escape, evaluated in order):
1. Stop codon in the transcript's last exon (or single-exon transcript).
2. Stop codon within `PTC_THRESHOLD_NT` (default 50) mRNA-nt of the last
   exon-exon junction.
3. Start-proximal: CDS shorter than `START_PROXIMAL_NT` (default 150) bp,
   permitting re-initiation past the PTC.
4. Last exon longer than `LONG_LAST_EXON_NT` (default 400) bp.
Otherwise the transcript is NMD-sensitive (50nt rule violated).

Confidence drops to MEDIUM when the propagation outcome is DISRUPTED and to
NONE for outcomes where NMD cannot be evaluated (no parent, no parent CDS,
start lost, stop not observed).
"""

from enum import Enum

import pandas as pd
import pyranges as pr

from craft.core.orf.confidence import ORFConfidence
from craft.core.orf.propagation import ORFOutcome

PTC_THRESHOLD_NT = 50
START_PROXIMAL_NT = 150
LONG_LAST_EXON_NT = 400


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


def predict(
    classified: pr.PyRanges,
    propagated: pd.DataFrame,
    ptc_threshold_nt: int = PTC_THRESHOLD_NT,
    start_proximal_nt: int = START_PROXIMAL_NT,
    long_last_exon_nt: int = LONG_LAST_EXON_NT,
) -> pd.DataFrame:
    """Predict per-isoform NMD susceptibility from propagated ORFs.

    Args:
        classified: PyRanges of isoform exons (with ``transcript_id`` and ``Strand``).
        propagated: DataFrame returned by :func:`craft.core.orf.propagation.propagate`.
        ptc_threshold_nt: 50nt PTC rule threshold (mRNA distance to last junction).
        start_proximal_nt: start-proximal escape window (CDS length below this escapes).
        long_last_exon_nt: long-last-exon escape threshold.

    Returns:
        DataFrame with one row per isoform and columns: ``transcript_id``,
        ``nmd_status``, ``nmd_rule``, ``stop_to_last_junction_nt``,
        ``last_exon_length_nt``, ``nmd_confidence``.
    """
    cols = [
        "transcript_id",
        "nmd_status",
        "nmd_rule",
        "stop_to_last_junction_nt",
        "last_exon_length_nt",
        "nmd_confidence",
    ]
    if propagated.empty or len(classified) == 0:
        return pd.DataFrame(columns=cols)

    iso_df = classified.df
    iso_strand = iso_df.groupby("transcript_id")["Strand"].first().to_dict()
    iso_exons_by_tx = {tx: g for tx, g in iso_df.groupby("transcript_id", sort=False)}

    rows: list[dict] = []
    for _, prop_row in propagated.iterrows():
        tx_id = prop_row["transcript_id"]
        outcome = ORFOutcome(prop_row["orf_outcome"])
        intervals = prop_row["propagated_cds_intervals"]
        stop_covered = bool(prop_row["stop_codon_covered"])

        applicable = (
            outcome in (ORFOutcome.PROPAGATED_INTACT, ORFOutcome.DISRUPTED)
            and stop_covered
            and intervals
        )
        if not applicable:
            rows.append(
                {
                    "transcript_id": tx_id,
                    "nmd_status": NMDStatus.NOT_APPLICABLE.value,
                    "nmd_rule": "",
                    "stop_to_last_junction_nt": None,
                    "last_exon_length_nt": None,
                    "nmd_confidence": ORFConfidence.NONE.value,
                }
            )
            continue

        strand = str(iso_strand[tx_id])
        iso_exons = iso_exons_by_tx[tx_id]
        stop_pos = _stop_codon_genome(intervals, strand)
        distance, in_last = _distance_stop_to_last_junction(stop_pos, iso_exons, strand)
        last_exon_len = _last_exon_length(iso_exons, strand)
        cds_bp = int(prop_row["propagated_cds_bp"])

        if in_last:
            status, rule = NMDStatus.ESCAPED, "stop_in_last_exon"
        elif distance <= ptc_threshold_nt:
            status, rule = NMDStatus.ESCAPED, "within_50nt_of_last_junction"
        elif cds_bp < start_proximal_nt:
            status, rule = NMDStatus.ESCAPED, "start_proximal"
        elif last_exon_len > long_last_exon_nt:
            status, rule = NMDStatus.ESCAPED, "long_last_exon"
        else:
            status, rule = NMDStatus.SENSITIVE, "ptc_50nt_rule"

        confidence = (
            ORFConfidence.HIGH
            if outcome == ORFOutcome.PROPAGATED_INTACT
            else ORFConfidence.MEDIUM
        )
        rows.append(
            {
                "transcript_id": tx_id,
                "nmd_status": status.value,
                "nmd_rule": rule,
                "stop_to_last_junction_nt": int(distance),
                "last_exon_length_nt": int(last_exon_len),
                "nmd_confidence": confidence.value,
            }
        )
    return pd.DataFrame(rows, columns=cols)


# Resolved-ORF statuses (from craft.core.orf.resolve) for which NMD can be evaluated:
# a real in-frame stop was found in the isoform's own sequence.
_RESOLVED_NMD_APPLICABLE = frozenset(
    {"intact", "ptc_premature", "ptc_intron_retained", "cds_extension"}
)


def predict_resolved(
    classified: pr.PyRanges,
    resolved: pd.DataFrame,
    ptc_threshold_nt: int = PTC_THRESHOLD_NT,
    start_proximal_nt: int = START_PROXIMAL_NT,
    long_last_exon_nt: int = LONG_LAST_EXON_NT,
) -> pd.DataFrame:
    """Predict NMD from the sequence-resolved stop instead of the geometric one.

    Mirrors :func:`predict` but reads the true stop position from
    :func:`craft.core.orf.resolve.resolve` output. Emits the parallel
    ``nmd_status_resolved`` / ``nmd_rule_resolved`` / ``nmd_confidence_resolved``
    columns and leaves the geometric columns untouched.

    Args:
        classified: PyRanges of isoform exons (``transcript_id``, ``Strand``).
        resolved: DataFrame from :func:`craft.core.orf.resolve.resolve`.
        ptc_threshold_nt: 50nt PTC rule threshold (mRNA distance to last junction).
        start_proximal_nt: start-proximal escape window.
        long_last_exon_nt: long-last-exon escape threshold.

    Returns:
        DataFrame with ``transcript_id``, ``nmd_status_resolved``,
        ``nmd_rule_resolved``, ``nmd_confidence_resolved``.
    """
    cols = ["transcript_id", "nmd_status_resolved", "nmd_rule_resolved", "nmd_confidence_resolved"]
    if resolved.empty or len(classified) == 0:
        return pd.DataFrame(columns=cols)

    iso_df = classified.df
    iso_strand = iso_df.groupby("transcript_id")["Strand"].first().to_dict()
    iso_exons_by_tx = {tx: g for tx, g in iso_df.groupby("transcript_id", sort=False)}

    rows: list[dict] = []
    for _, res_row in resolved.iterrows():
        tx_id = res_row["transcript_id"]
        status_label = str(res_row["resolved_orf_status"])
        intervals = res_row["resolved_cds_intervals"]
        applicable = (
            status_label in _RESOLVED_NMD_APPLICABLE
            and bool(res_row["stop_in_transcript"])
            and intervals
        )
        if not applicable:
            rows.append(
                {
                    "transcript_id": tx_id,
                    "nmd_status_resolved": NMDStatus.NOT_APPLICABLE.value,
                    "nmd_rule_resolved": "",
                    "nmd_confidence_resolved": ORFConfidence.NONE.value,
                }
            )
            continue

        strand = str(iso_strand[tx_id])
        iso_exons = iso_exons_by_tx[tx_id]
        stop_pos = _stop_codon_genome(intervals, strand)
        distance, in_last = _distance_stop_to_last_junction(stop_pos, iso_exons, strand)
        last_exon_len = _last_exon_length(iso_exons, strand)
        cds_bp = int(res_row["resolved_cds_bp"])

        if in_last:
            status, rule = NMDStatus.ESCAPED, "stop_in_last_exon"
        elif distance <= ptc_threshold_nt:
            status, rule = NMDStatus.ESCAPED, "within_50nt_of_last_junction"
        elif cds_bp < start_proximal_nt:
            status, rule = NMDStatus.ESCAPED, "start_proximal"
        elif last_exon_len > long_last_exon_nt:
            status, rule = NMDStatus.ESCAPED, "long_last_exon"
        else:
            status, rule = NMDStatus.SENSITIVE, "ptc_50nt_rule"

        confidence = ORFConfidence.HIGH if status_label == "intact" else ORFConfidence.MEDIUM
        rows.append(
            {
                "transcript_id": tx_id,
                "nmd_status_resolved": status.value,
                "nmd_rule_resolved": rule,
                "nmd_confidence_resolved": confidence.value,
            }
        )
    return pd.DataFrame(rows, columns=cols)


def predict_denovo(
    classified: pr.PyRanges,
    denovo: pd.DataFrame,
    ptc_threshold_nt: int = PTC_THRESHOLD_NT,
    start_proximal_nt: int = START_PROXIMAL_NT,
    long_last_exon_nt: int = LONG_LAST_EXON_NT,
) -> pd.DataFrame:
    """Predict NMD from a de novo ORF, for isoforms with no usable parent.

    Orphan isoforms (``no_parent`` / ``no_parent_cds`` / ``start_lost``) carry no
    reference-anchored stop, so :func:`predict` and :func:`predict_resolved` leave
    them ``not_applicable``. This applies the same escape-rule cascade to the de
    novo ORF's stop instead, so a novel isoform with a predicted ORF still gets an
    NMD call. Confidence is always ``low``: the stop comes from a predicted ORF,
    not a curated reference.

    Args:
        classified: PyRanges of isoform exons (``transcript_id``, ``Strand``).
        denovo: DataFrame from :func:`craft.core.orf.denovo.predict` (needs
            ``denovo_orf_found``, ``denovo_cds_intervals``, ``denovo_cds_bp``).
        ptc_threshold_nt: 50nt PTC rule threshold (mRNA distance to last junction).
        start_proximal_nt: start-proximal escape window.
        long_last_exon_nt: long-last-exon escape threshold.

    Returns:
        DataFrame with ``transcript_id``, ``nmd_status_denovo``,
        ``nmd_rule_denovo``, ``nmd_confidence_denovo``.
    """
    cols = ["transcript_id", "nmd_status_denovo", "nmd_rule_denovo", "nmd_confidence_denovo"]
    if denovo is None or denovo.empty or len(classified) == 0:
        return pd.DataFrame(columns=cols)

    iso_df = classified.df
    iso_strand = iso_df.groupby("transcript_id")["Strand"].first().to_dict()
    iso_exons_by_tx = {tx: g for tx, g in iso_df.groupby("transcript_id", sort=False)}

    rows: list[dict] = []
    for _, dn_row in denovo.iterrows():
        tx_id = dn_row["transcript_id"]
        intervals = dn_row["denovo_cds_intervals"]
        found = bool(dn_row["denovo_orf_found"]) if pd.notna(dn_row["denovo_orf_found"]) else False
        applicable = (
            found and isinstance(intervals, list) and len(intervals) > 0
            and tx_id in iso_exons_by_tx
        )
        if not applicable:
            rows.append(
                {
                    "transcript_id": tx_id,
                    "nmd_status_denovo": NMDStatus.NOT_APPLICABLE.value,
                    "nmd_rule_denovo": "",
                    "nmd_confidence_denovo": ORFConfidence.NONE.value,
                }
            )
            continue

        strand = str(iso_strand[tx_id])
        iso_exons = iso_exons_by_tx[tx_id]
        stop_pos = _stop_codon_genome(intervals, strand)
        distance, in_last = _distance_stop_to_last_junction(stop_pos, iso_exons, strand)
        last_exon_len = _last_exon_length(iso_exons, strand)
        cds_bp = int(dn_row["denovo_cds_bp"])

        if in_last:
            status, rule = NMDStatus.ESCAPED, "stop_in_last_exon"
        elif distance <= ptc_threshold_nt:
            status, rule = NMDStatus.ESCAPED, "within_50nt_of_last_junction"
        elif cds_bp < start_proximal_nt:
            status, rule = NMDStatus.ESCAPED, "start_proximal"
        elif last_exon_len > long_last_exon_nt:
            status, rule = NMDStatus.ESCAPED, "long_last_exon"
        else:
            status, rule = NMDStatus.SENSITIVE, "ptc_50nt_rule"

        rows.append(
            {
                "transcript_id": tx_id,
                "nmd_status_denovo": status.value,
                "nmd_rule_denovo": rule,
                "nmd_confidence_denovo": ORFConfidence.LOW.value,
            }
        )
    return pd.DataFrame(rows, columns=cols)
