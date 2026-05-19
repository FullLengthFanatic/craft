"""Reference-isoform ORF propagation.

Project the parent transcript's CDS coordinates onto the novel isoform where structure
is preserved. The v1 algorithm intersects the parent's CDS intervals with the
isoform's exons and classifies the outcome by start- and stop-codon coverage and
total preserved CDS length.
"""

from enum import Enum

import pandas as pd
import pyranges as pr


class ORFOutcome(str, Enum):
    """Per-isoform ORF propagation outcome."""

    PROPAGATED_INTACT = "propagated_intact"
    DISRUPTED = "disrupted"
    START_LOST = "start_lost"
    STOP_NOT_OBSERVED = "stop_not_observed"
    STOP_AT_ALT_POLYA = "stop_at_alt_polya"
    NO_PARENT = "no_parent"
    NO_PARENT_CDS = "no_parent_cds"


def _start_codon_pos(parent_cds: pd.DataFrame, strand: str) -> int:
    """Genomic 0-based position of the start codon's first base."""
    if strand == "+":
        return int(parent_cds["Start"].min())
    if strand == "-":
        return int(parent_cds["End"].max()) - 1
    raise ValueError(f"Unsupported strand: {strand!r}")


def _stop_codon_pos(parent_cds: pd.DataFrame, strand: str) -> int:
    """Genomic 0-based position of the last CDS base (the stop codon area)."""
    if strand == "+":
        return int(parent_cds["End"].max()) - 1
    if strand == "-":
        return int(parent_cds["Start"].min())
    raise ValueError(f"Unsupported strand: {strand!r}")


def _position_in_exons(pos: int, exons: pd.DataFrame) -> bool:
    return bool(((exons["Start"] <= pos) & (pos < exons["End"])).any())


def _intersect_intervals(
    a_exons: pd.DataFrame,
    b_cds: pd.DataFrame,
    chrom: str,
    strand: str,
) -> list[tuple[str, int, int, str]]:
    """Pairwise interval intersection on same chromosome and strand."""
    if a_exons.empty or b_cds.empty:
        return []
    out: list[tuple[str, int, int, str]] = []
    a_arr = a_exons[["Start", "End"]].to_numpy()
    b_arr = b_cds[["Start", "End"]].to_numpy()
    for a_start, a_end in a_arr:
        for b_start, b_end in b_arr:
            start = max(int(a_start), int(b_start))
            end = min(int(a_end), int(b_end))
            if end > start:
                out.append((chrom, start, end, strand))
    out.sort(key=lambda x: x[1])
    return out


def _propagate_one(
    iso_exons: pd.DataFrame,
    parent_cds: pd.DataFrame,
    strand: str,
) -> tuple[ORFOutcome, list[tuple[str, int, int, str]], int, bool, bool]:
    chrom = str(iso_exons["Chromosome"].iloc[0])
    start_pos = _start_codon_pos(parent_cds, strand)
    stop_pos = _stop_codon_pos(parent_cds, strand)
    start_covered = _position_in_exons(start_pos, iso_exons)
    stop_covered = _position_in_exons(stop_pos, iso_exons)
    propagated = _intersect_intervals(iso_exons, parent_cds, chrom, strand)
    propagated_bp = sum(e - s for _, s, e, _ in propagated)
    parent_bp = int((parent_cds["End"] - parent_cds["Start"]).sum())

    if not start_covered:
        outcome = ORFOutcome.START_LOST
    elif not stop_covered:
        outcome = ORFOutcome.STOP_NOT_OBSERVED
    elif propagated_bp == parent_bp:
        outcome = ORFOutcome.PROPAGATED_INTACT
    else:
        outcome = ORFOutcome.DISRUPTED
    return outcome, propagated, propagated_bp, start_covered, stop_covered


def propagate(
    classified: pr.PyRanges,
    reference: pr.PyRanges,
) -> pd.DataFrame:
    """Propagate parent CDS coordinates onto novel isoforms.

    For each isoform with an identified parent reference transcript (``parent_tx_id``
    column on ``classified``), intersect the parent's CDS intervals with the
    isoform's exons. The per-isoform outcome is determined by whether the parent's
    start and stop codon genomic positions are observed in the isoform and by how
    much of the parent's CDS is preserved.

    Args:
        classified: Output of :func:`craft.core.completeness.classify` (isoform
            exons with ``transcript_id`` and ``parent_tx_id`` columns).
        reference: Reference annotation PyRanges with a ``Feature`` column
            distinguishing ``exon`` from ``CDS`` rows.

    Returns:
        DataFrame with one row per isoform transcript and columns:
        ``transcript_id``, ``parent_tx_id``, ``orf_outcome``, ``propagated_cds_bp``,
        ``parent_cds_bp``, ``start_codon_covered``, ``stop_codon_covered``,
        ``propagated_cds_intervals`` (list of ``(Chromosome, Start, End, Strand)``).
    """
    columns = [
        "transcript_id",
        "parent_tx_id",
        "orf_outcome",
        "propagated_cds_bp",
        "parent_cds_bp",
        "start_codon_covered",
        "stop_codon_covered",
        "propagated_cds_intervals",
    ]
    if len(classified) == 0:
        return pd.DataFrame(columns=columns)

    iso_df = classified.df
    ref_df = reference.df
    parent_cds_all = ref_df[ref_df["Feature"] == "CDS"]

    iso_strand = iso_df.groupby("transcript_id")["Strand"].first().to_dict()
    iso_parent = iso_df.groupby("transcript_id")["parent_tx_id"].first().to_dict()
    iso_exons_by_tx = {tx: g for tx, g in iso_df.groupby("transcript_id", sort=False)}
    parent_cds_by_tx = {
        tx: g for tx, g in parent_cds_all.groupby("transcript_id", sort=False)
    }

    rows: list[dict] = []
    for tx_id, iso_exons in iso_exons_by_tx.items():
        parent_tx = iso_parent.get(tx_id, "")
        strand = iso_strand.get(tx_id, "+")

        if not parent_tx:
            rows.append(
                {
                    "transcript_id": tx_id,
                    "parent_tx_id": "",
                    "orf_outcome": ORFOutcome.NO_PARENT.value,
                    "propagated_cds_bp": 0,
                    "parent_cds_bp": 0,
                    "start_codon_covered": False,
                    "stop_codon_covered": False,
                    "propagated_cds_intervals": [],
                }
            )
            continue

        if parent_tx not in parent_cds_by_tx:
            rows.append(
                {
                    "transcript_id": tx_id,
                    "parent_tx_id": parent_tx,
                    "orf_outcome": ORFOutcome.NO_PARENT_CDS.value,
                    "propagated_cds_bp": 0,
                    "parent_cds_bp": 0,
                    "start_codon_covered": False,
                    "stop_codon_covered": False,
                    "propagated_cds_intervals": [],
                }
            )
            continue

        parent_cds = parent_cds_by_tx[parent_tx]
        outcome, intervals, prop_bp, start_cov, stop_cov = _propagate_one(
            iso_exons, parent_cds, strand
        )
        parent_bp = int((parent_cds["End"] - parent_cds["Start"]).sum())
        rows.append(
            {
                "transcript_id": tx_id,
                "parent_tx_id": parent_tx,
                "orf_outcome": outcome.value,
                "propagated_cds_bp": prop_bp,
                "parent_cds_bp": parent_bp,
                "start_codon_covered": start_cov,
                "stop_codon_covered": stop_cov,
                "propagated_cds_intervals": intervals,
            }
        )

    return pd.DataFrame(rows, columns=columns)
