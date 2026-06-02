"""Completeness classification: full-length vs truncated isoforms vs reference."""

from enum import Enum

import pandas as pd
import pyranges as pr

from craft.core.intervals import splice_junctions


class Completeness(str, Enum):
    """Structural completeness of a novel isoform relative to a reference parent."""

    FULL_LENGTH = "full_length"
    TRUNCATED_5P = "truncated_5p"
    TRUNCATED_3P = "truncated_3p"
    TRUNCATED_BOTH = "truncated_both"
    INTERNAL_FRAGMENT = "internal_fragment"
    NOVEL_NO_MATCH = "novel_no_match"
    ALT_3PRIME_END = "alt_3prime_end"


def _transcript_spans(exons: pr.PyRanges) -> pd.DataFrame:
    """Return per-transcript span (one row per transcript_id)."""
    df = exons.df
    return df.groupby("transcript_id", as_index=False).agg(
        Chromosome=("Chromosome", "first"),
        Strand=("Strand", "first"),
        Start=("Start", "min"),
        End=("End", "max"),
    )


def _shared_junction_counts(iso_jx: pr.PyRanges, ref_jx: pr.PyRanges) -> pd.DataFrame:
    """Per (iso_tx, ref_tx) pair, count of exactly matching splice junctions."""
    empty = pd.DataFrame(columns=["iso_tx", "ref_tx", "shared_jx"])
    if len(iso_jx) == 0 or len(ref_jx) == 0:
        return empty
    iso_df = iso_jx.df.rename(columns={"transcript_id": "iso_tx"})
    ref_df = ref_jx.df.rename(columns={"transcript_id": "ref_tx"})
    merged = iso_df.merge(
        ref_df,
        on=["Chromosome", "Start", "End", "Strand"],
        how="inner",
        suffixes=("", "_ref"),
    )
    if merged.empty:
        return empty
    counts = merged.groupby(["iso_tx", "ref_tx"], as_index=False).size()
    return counts.rename(columns={"size": "shared_jx"})[["iso_tx", "ref_tx", "shared_jx"]]


def _exon_overlap_bp(iso_exons: pr.PyRanges, ref_exons: pr.PyRanges) -> pd.DataFrame:
    """Per (iso_tx, ref_tx) pair, total bp of stranded exon overlap."""
    empty = pd.DataFrame(columns=["iso_tx", "ref_tx", "overlap_bp"])
    if len(iso_exons) == 0 or len(ref_exons) == 0:
        return empty
    joined = iso_exons.join(ref_exons, strandedness="same")
    if len(joined) == 0:
        return empty
    df = joined.df
    df["overlap_bp"] = (
        df[["End", "End_b"]].min(axis=1) - df[["Start", "Start_b"]].max(axis=1)
    ).clip(lower=0)
    total = df.groupby(["transcript_id", "transcript_id_b"], as_index=False)["overlap_bp"].sum()
    total = total[total["overlap_bp"] > 0].rename(
        columns={"transcript_id": "iso_tx", "transcript_id_b": "ref_tx"}
    )
    return total[["iso_tx", "ref_tx", "overlap_bp"]]


def _select_parent(
    jx_counts: pd.DataFrame,
    overlap_bp: pd.DataFrame,
    cds_tx_ids: set[str] | None = None,
    prefer_coding_parent: bool = False,
) -> pd.DataFrame:
    """Per iso_tx, pick the ref_tx with max shared junctions, tiebreak by overlap_bp.

    With ``prefer_coding_parent=True``, a CDS-bearing reference transcript is
    preferred over a non-coding one as the lowest-priority tiebreaker (only when
    shared junctions and exon overlap are equal). The default is off, so the
    parent selection is byte-identical to v1.4.
    """
    if jx_counts.empty and overlap_bp.empty:
        return pd.DataFrame(columns=["iso_tx", "ref_tx", "shared_jx", "overlap_bp"])
    scored = jx_counts.merge(overlap_bp, on=["iso_tx", "ref_tx"], how="outer")
    scored["shared_jx"] = scored["shared_jx"].fillna(0).astype("int64")
    scored["overlap_bp"] = scored["overlap_bp"].fillna(0).astype("int64")
    sort_cols = ["iso_tx", "shared_jx", "overlap_bp"]
    ascending = [True, False, False]
    if prefer_coding_parent and cds_tx_ids:
        scored["ref_has_cds"] = scored["ref_tx"].isin(cds_tx_ids).astype("int64")
        sort_cols.append("ref_has_cds")
        ascending.append(False)
    scored = scored.sort_values(sort_cols, ascending=ascending, kind="stable")
    return scored.groupby("iso_tx", as_index=False).first()


def _classify_one(
    iso_start: int,
    iso_end: int,
    parent_start: int,
    parent_end: int,
    strand: str,
    tolerance: int,
) -> Completeness:
    """Classify a single isoform/parent pair by end-position comparison."""
    if strand == "+":
        five_prime_complete = iso_start <= parent_start + tolerance
        three_prime_complete = iso_end >= parent_end - tolerance
    elif strand == "-":
        five_prime_complete = iso_end >= parent_end - tolerance
        three_prime_complete = iso_start <= parent_start + tolerance
    else:
        raise ValueError(f"Unsupported strand: {strand!r}")

    if five_prime_complete and three_prime_complete:
        return Completeness.FULL_LENGTH
    if three_prime_complete:
        return Completeness.TRUNCATED_5P
    if five_prime_complete:
        return Completeness.TRUNCATED_3P

    interior_pad = max(tolerance * 2, 100)
    if iso_start > parent_start + interior_pad and iso_end < parent_end - interior_pad:
        return Completeness.INTERNAL_FRAGMENT
    return Completeness.TRUNCATED_BOTH


def classify(
    isoforms: pr.PyRanges,
    reference: pr.PyRanges,
    tolerance: int = 50,
    cds_tx_ids: set[str] | None = None,
    prefer_coding_parent: bool = False,
) -> pr.PyRanges:
    """Classify each isoform's completeness vs its best-matching reference transcript.

    Best match is the reference transcript with the most exactly shared splice
    junctions; ties are broken by total stranded exon-overlap bp. An isoform with
    no shared junctions and no overlapping exons is classified ``NOVEL_NO_MATCH``.

    Args:
        isoforms: PyRanges of isoform exons with ``transcript_id`` and ``Strand``.
        reference: PyRanges of reference transcript exons (same columns).
        tolerance: Allowed slack (bp) on each end before considering it truncated.

    Returns:
        Input isoform PyRanges with added per-row columns ``completeness``,
        ``parent_tx_id``, ``shared_junctions``, ``parent_overlap_bp``.
    """
    if len(isoforms) == 0:
        return pr.PyRanges()

    iso_spans = _transcript_spans(isoforms).set_index("transcript_id")
    ref_spans = _transcript_spans(reference).set_index("transcript_id")
    parents = _select_parent(
        _shared_junction_counts(splice_junctions(isoforms), splice_junctions(reference)),
        _exon_overlap_bp(isoforms, reference),
        cds_tx_ids=cds_tx_ids,
        prefer_coding_parent=prefer_coding_parent,
    )

    parent_lookup = parents.set_index("iso_tx") if not parents.empty else parents

    rows: list[dict] = []
    for iso_tx, span in iso_spans.iterrows():
        if iso_tx not in parent_lookup.index:
            rows.append(
                {
                    "transcript_id": iso_tx,
                    "completeness": Completeness.NOVEL_NO_MATCH.value,
                    "parent_tx_id": "",
                    "shared_junctions": 0,
                    "parent_overlap_bp": 0,
                    "has_cds_bearing_parent": False,
                }
            )
            continue
        match = parent_lookup.loc[iso_tx]
        ref_tx = str(match["ref_tx"])
        ref_span = ref_spans.loc[ref_tx]
        category = _classify_one(
            iso_start=int(span["Start"]),
            iso_end=int(span["End"]),
            parent_start=int(ref_span["Start"]),
            parent_end=int(ref_span["End"]),
            strand=str(span["Strand"]),
            tolerance=tolerance,
        )
        rows.append(
            {
                "transcript_id": iso_tx,
                "completeness": category.value,
                "parent_tx_id": ref_tx,
                "shared_junctions": int(match["shared_jx"]),
                "parent_overlap_bp": int(match["overlap_bp"]),
                "has_cds_bearing_parent": bool(cds_tx_ids and ref_tx in cds_tx_ids),
            }
        )

    per_tx = pd.DataFrame(rows)
    annotated = isoforms.df.merge(per_tx, on="transcript_id", how="left")
    return pr.PyRanges(annotated)
