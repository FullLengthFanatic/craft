"""Rule-based RNA-surveillance susceptibility annotation.

NMD is predicted from the **sequence-resolved** ORF stop (the real in-frame stop
found by translating the isoform's own spliced CDS). For orphan isoforms that
have no reference-anchored ORF, the call falls back to the de-novo ORF's stop.
There is a single structural NMD hypothesis per isoform; ``nmd_basis`` records which ORF it came
from (``resolved`` / ``denovo`` / ``none``).

Rule cascade (each is a sufficient condition for escape, evaluated in order):
1. Stop codon in the transcript's last exon (or single-exon transcript).
2. Stop codon within ``PTC_THRESHOLD_NT`` (default 50) mRNA-nt of the last junction.
3. Start-proximal: CDS shorter than ``START_PROXIMAL_NT`` (default 150) bp.
4. The exon containing the PTC is longer than ``LONG_LAST_EXON_NT`` (default 400) bp
   (the long-exon rule; Lindeboom et al. 2016, on the PTC-bearing exon).
Otherwise the transcript is structurally compatible with NMD (50nt rule violated).

Confidence is ``high`` for a resolved intact ORF, ``medium`` for a resolved but
altered ORF (PTC / intron-retained / extension), ``low`` for a de-novo call (the
stop is from a predicted ORF, not a reference), and ``none`` when not applicable.
"""

from enum import Enum

import pandas as pd
import pyranges as pr

from craft.core.intervals import (
    genomic_position_at_transcript_coordinate,
    spliced_length,
    transcript_coordinate,
)
from craft.core.orf.confidence import ORFConfidence

PTC_THRESHOLD_NT = 50
START_PROXIMAL_NT = 150
LONG_LAST_EXON_NT = 400

# Resolved-ORF statuses (from craft.core.orf.resolve) that carry a real stop.
_RESOLVED_WITH_STOP = frozenset(
    {"intact", "ptc_premature", "ptc_intron_retained", "cds_extension", "start_rescued"}
)

_COLUMNS = [
    "transcript_id",
    "nmd_status",
    "nmd_rule",
    "nmd_confidence",
    "nmd_basis",
    "stop_to_last_junction_nt",
    "last_exon_length_nt",
    "ptc_exon_length_nt",
    "nmd_susceptibility",
    "nmd_rule_score",
    "nmd_evidence_tier",
    "surveillance_status",
    "surveillance_mechanism",
    "nonstop_decay_candidate",
    "surveillance_limitations",
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
    stop_codon_pos: int, exons: pd.DataFrame, strand: str
) -> tuple[int, bool]:
    """Distance (mRNA bp) from the stop codon to the last exon-exon junction.

    Returns ``(distance, is_in_last_exon)``. For single-exon transcripts and for
    stops in the transcript's last exon, returns ``(0, True)``.
    """
    ordered = exons.sort_values("Start", ascending=strand == "+", kind="stable")
    n_exons = len(ordered)
    if n_exons <= 1:
        return 0, True
    stop_tx = transcript_coordinate(exons, stop_codon_pos, strand)
    if stop_tx is None:
        return -1, False
    last_exon_start_tx = spliced_length(ordered.iloc[:-1])
    if stop_tx >= last_exon_start_tx:
        return 0, True
    # Transcript coordinates handle stop codons that cross an exon junction.
    return last_exon_start_tx - (stop_tx + 3), False


def _last_exon_length(exons: pd.DataFrame, strand: str) -> int:
    sorted_exons = exons.sort_values("Start").reset_index(drop=True)
    if strand == "+":
        last = sorted_exons.iloc[-1]
    elif strand == "-":
        last = sorted_exons.iloc[0]
    else:
        raise ValueError(f"Unsupported strand: {strand!r}")
    return int(last["End"]) - int(last["Start"])


def _ptc_exon_length(stop_pos: int, exons: pd.DataFrame) -> int:
    """Length (bp) of the exon that contains the stop codon; 0 if none does.

    The long-exon NMD-escape rule applies to the exon carrying the PTC (a long
    exon deposits fewer EJCs per unit length), not the transcript's terminal exon.
    """
    mask = (exons["Start"] <= stop_pos) & (stop_pos < exons["End"])
    if not mask.any():
        return 0
    ex = exons[mask].iloc[0]
    return int(ex["End"]) - int(ex["Start"])


def _cascade(
    distance: int,
    in_last: bool,
    cds_bp: int,
    ptc_exon_len: int,
    ptc_threshold_nt: int,
    start_proximal_nt: int,
    long_last_exon_nt: int,
) -> tuple[NMDStatus, str, float]:
    if in_last:
        return NMDStatus.ESCAPED, "stop_in_last_exon", 0.05
    if distance <= ptc_threshold_nt:
        return NMDStatus.ESCAPED, "within_50nt_of_last_junction", 0.10
    if cds_bp < start_proximal_nt:
        return NMDStatus.ESCAPED, "start_proximal", 0.20
    if ptc_exon_len > long_last_exon_nt:
        return NMDStatus.ESCAPED, "long_exon", 0.25
    return NMDStatus.SENSITIVE, "ptc_50nt_rule", 0.80


def _not_applicable(tx_id: str) -> dict:
    return {
        "transcript_id": tx_id,
        "nmd_status": NMDStatus.NOT_APPLICABLE.value,
        "nmd_rule": "",
        "nmd_confidence": ORFConfidence.NONE.value,
        "nmd_basis": "none",
        "stop_to_last_junction_nt": None,
        "last_exon_length_nt": None,
        "ptc_exon_length_nt": None,
        "nmd_susceptibility": "indeterminate",
        "nmd_rule_score": None,
        "nmd_evidence_tier": "none",
        "surveillance_status": "indeterminate",
        "surveillance_mechanism": "none",
        "nonstop_decay_candidate": False,
        "surveillance_limitations": "no_resolved_termination_event",
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
    per_tx_meta = iso_df.groupby("transcript_id").first()

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
        evidence_tier = "none"
        limitation = ""

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
            status_str = str(res["resolved_orf_status"])
            if status_str == "intact":
                confidence = ORFConfidence.HIGH
            elif status_str == "start_rescued":
                # The start is inferred (frame-anchored), so cap confidence low.
                confidence = ORFConfidence.LOW
            else:
                confidence = ORFConfidence.MEDIUM
            meta = per_tx_meta.loc[tx_id]
            parent_ambiguous = bool(meta.get("parent_ambiguous", False))
            reference_complete = bool(meta.get("reference_cds_complete", True))
            reference_phase_valid = bool(meta.get("reference_cds_phase_valid", True))
            reference_issues: list[str] = []
            if parent_ambiguous:
                reference_issues.append("ambiguous_reference_parent")
            if not reference_complete:
                reference_issues.append("incomplete_reference_cds")
            if not reference_phase_valid:
                reference_issues.append("inconsistent_reference_cds_phase")
            if reference_issues:
                evidence_tier = "moderate" if status_str == "intact" else "limited"
                limitation = ";".join(reference_issues)
            else:
                evidence_tier = "strong" if status_str == "intact" else "moderate"
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
                evidence_tier = "limited"
                limitation = "de_novo_orf"

        if intervals is None:
            row = _not_applicable(tx_id)
            res = resolved_by_tx.get(tx_id)
            resolved_status = str(res["resolved_orf_status"]) if res is not None else ""
            completeness = str(per_tx_meta.loc[tx_id].get("completeness", ""))
            if resolved_status == "right_censored" and completeness == "alt_3prime_end":
                row.update(
                    {
                        "surveillance_status": "candidate",
                        "surveillance_mechanism": "nonstop_decay",
                        "nonstop_decay_candidate": True,
                        "surveillance_limitations": (
                            "polyadenylated_transcript_without_observed_stop"
                        ),
                    }
                )
            elif resolved_status in {"left_censored", "right_censored"}:
                row["surveillance_limitations"] = f"{resolved_status}_orf"
            rows.append(row)
            continue

        strand = str(iso_strand[tx_id])
        iso_exons = iso_exons_by_tx[tx_id]
        res = resolved_by_tx.get(tx_id)
        if res is not None and pd.notna(res.get("resolved_stop_codon_pos")):
            stop_pos = int(res["resolved_stop_codon_pos"])
        else:
            # De-novo intervals omit the stop codon; advance one spliced base
            # from the last sense base, including across an exon junction.
            last_sense = _stop_codon_genome(intervals, strand)
            last_sense_tx = transcript_coordinate(iso_exons, last_sense, strand)
            stop_pos = (
                genomic_position_at_transcript_coordinate(
                    iso_exons, last_sense_tx + 1, strand
                )
                if last_sense_tx is not None else None
            )
            if stop_pos is None:
                rows.append(_not_applicable(tx_id))
                continue
        distance, in_last = _distance_stop_to_last_junction(stop_pos, iso_exons, strand)
        last_exon_len = _last_exon_length(iso_exons, strand)
        ptc_exon_len = _ptc_exon_length(stop_pos, iso_exons)
        status, rule, rule_score = _cascade(
            distance, in_last, cds_bp, ptc_exon_len,
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
                "ptc_exon_length_nt": int(ptc_exon_len),
                "nmd_susceptibility": (
                    "likely_sensitive" if status == NMDStatus.SENSITIVE else "likely_escape"
                ),
                "nmd_rule_score": rule_score,
                "nmd_evidence_tier": evidence_tier,
                "surveillance_status": "candidate" if status == NMDStatus.SENSITIVE else "none",
                "surveillance_mechanism": "nmd" if status == NMDStatus.SENSITIVE else "none",
                "nonstop_decay_candidate": False,
                "surveillance_limitations": limitation,
            }
        )
    return pd.DataFrame(rows, columns=_COLUMNS)
