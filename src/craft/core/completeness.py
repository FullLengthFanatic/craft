"""Completeness classification and ambiguity-aware reference-parent selection."""

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
    iso_junction_totals: dict[str, int] | None = None,
    ref_junction_totals: dict[str, int] | None = None,
    iso_lengths: dict[str, int] | None = None,
    ref_lengths: dict[str, int] | None = None,
    iso_gene_ids: dict[str, str] | None = None,
    ref_gene_ids: dict[str, str] | None = None,
    reference_priority: dict[str, float] | None = None,
    parent_hints: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Rank candidate parents and retain ambiguity rather than hiding ties.

    Candidates are restricted to the input gene when both annotations carry a
    gene id.  An upstream associated-transcript hint is strongest, followed by
    exact/contained intron chains, junction F1, normalized exon overlap and
    curated-reference priority.  The returned margin and ambiguity flag make it
    explicit when several parents remain equally plausible.
    """
    if jx_counts.empty and overlap_bp.empty:
        return pd.DataFrame(columns=["iso_tx", "ref_tx", "shared_jx", "overlap_bp"])
    scored = jx_counts.merge(overlap_bp, on=["iso_tx", "ref_tx"], how="outer")
    scored["shared_jx"] = (
        pd.to_numeric(scored["shared_jx"], errors="coerce").fillna(0).astype("int64")
    )
    scored["overlap_bp"] = (
        pd.to_numeric(scored["overlap_bp"], errors="coerce").fillna(0).astype("int64")
    )
    iso_junction_totals = iso_junction_totals or {}
    ref_junction_totals = ref_junction_totals or {}
    iso_lengths = iso_lengths or {}
    ref_lengths = ref_lengths or {}
    iso_gene_ids = iso_gene_ids or {}
    ref_gene_ids = ref_gene_ids or {}
    reference_priority = reference_priority or {}
    parent_hints = parent_hints or {}

    scored["iso_junctions"] = scored["iso_tx"].map(iso_junction_totals).fillna(0).astype(int)
    scored["ref_junctions"] = scored["ref_tx"].map(ref_junction_totals).fillna(0).astype(int)
    scored["junction_precision"] = scored["shared_jx"] / scored["iso_junctions"].replace(0, 1)
    scored["junction_recall"] = scored["shared_jx"] / scored["ref_junctions"].replace(0, 1)
    denom = scored["junction_precision"] + scored["junction_recall"]
    scored["junction_f1"] = (
        2 * scored["junction_precision"] * scored["junction_recall"] / denom.replace(0, 1)
    )
    scored["exact_intron_chain"] = (
        (scored["iso_junctions"] > 0)
        & (scored["iso_junctions"] == scored["ref_junctions"])
        & (scored["shared_jx"] == scored["iso_junctions"])
    )
    scored["iso_chain_contained"] = (
        (scored["iso_junctions"] > 0) & (scored["shared_jx"] == scored["iso_junctions"])
    )
    min_len = pd.concat(
        [scored["iso_tx"].map(iso_lengths), scored["ref_tx"].map(ref_lengths)], axis=1
    ).min(axis=1).replace(0, 1)
    scored["overlap_fraction"] = (scored["overlap_bp"] / min_len).clip(upper=1.0)
    scored["ref_has_cds"] = scored["ref_tx"].isin(cds_tx_ids or set()).astype(int)
    scored["reference_priority"] = scored["ref_tx"].map(reference_priority).fillna(0.0)
    scored["hint_match"] = scored.apply(
        lambda row: parent_hints.get(str(row["iso_tx"]), "") == str(row["ref_tx"]), axis=1
    )

    iso_gene = scored["iso_tx"].map(iso_gene_ids).fillna("").astype(str)
    ref_gene = scored["ref_tx"].map(ref_gene_ids).fillna("").astype(str)
    scored["gene_match"] = (iso_gene.ne("") & ref_gene.ne("") & iso_gene.eq(ref_gene))
    scored["gene_conflict"] = (iso_gene.ne("") & ref_gene.ne("") & iso_gene.ne(ref_gene))
    # If a gene-consistent candidate exists, a conflicting overlap cannot win.
    has_gene_match = scored.groupby("iso_tx")["gene_match"].transform("any")
    scored = scored[~(has_gene_match & scored["gene_conflict"])].copy()

    scored["parent_match_score"] = (
        5.0 * scored["hint_match"].astype(float)
        + 3.0 * scored["gene_match"].astype(float)
        + 2.0 * scored["exact_intron_chain"].astype(float)
        + 1.5 * scored["iso_chain_contained"].astype(float)
        + scored["junction_f1"]
        + scored["overlap_fraction"]
        + scored["reference_priority"].clip(-100, 100) / 1000.0
        + (0.01 * scored["ref_has_cds"] if prefer_coding_parent else 0.0)
    )
    scored = scored.sort_values(
        [
            "iso_tx", "parent_match_score", "shared_jx", "overlap_bp",
            "reference_priority", "ref_tx",
        ],
        ascending=[True, False, False, False, False, True],
        kind="stable",
    )
    scored["parent_candidate_rank"] = scored.groupby("iso_tx").cumcount() + 1
    scored["parent_candidate_count"] = scored.groupby("iso_tx")["ref_tx"].transform("size")
    scored["second_parent_score"] = scored.groupby("iso_tx")["parent_match_score"].transform(
        lambda values: values.iloc[1] if len(values) > 1 else float("nan")
    )
    scored["parent_match_margin"] = (
        scored["parent_match_score"] - scored["second_parent_score"]
    )
    scored["parent_ambiguous"] = (
        scored["parent_candidate_count"].gt(1)
        & scored["parent_match_margin"].fillna(float("inf")).lt(0.05)
    )
    scored["parent_selection_reason"] = scored.apply(_selection_reason, axis=1)
    return scored[scored["parent_candidate_rank"] == 1].reset_index(drop=True)


def _selection_reason(row: pd.Series) -> str:
    reasons: list[str] = []
    if bool(row.get("hint_match", False)):
        reasons.append("upstream_hint")
    if bool(row.get("gene_match", False)):
        reasons.append("gene_match")
    if bool(row.get("exact_intron_chain", False)):
        reasons.append("exact_intron_chain")
    elif bool(row.get("iso_chain_contained", False)):
        reasons.append("contained_intron_chain")
    elif float(row.get("shared_jx", 0)) > 0:
        reasons.append("partial_junction_match")
    else:
        reasons.append("exon_overlap_only")
    if float(row.get("reference_priority", 0)) > 0:
        reasons.append("curated_reference")
    return ";".join(reasons)


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
    reference_metadata: pd.DataFrame | None = None,
    parent_hints: dict[str, str] | None = None,
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
    iso_jx = splice_junctions(isoforms)
    ref_jx = splice_junctions(reference)
    iso_df = isoforms.df
    ref_df = reference.df
    iso_gene_ids = (
        iso_df.groupby("transcript_id")["gene_id"].first().fillna("").astype(str).to_dict()
        if "gene_id" in iso_df.columns else {}
    )
    ref_gene_ids = (
        ref_df.groupby("transcript_id")["gene_id"].first().fillna("").astype(str).to_dict()
        if "gene_id" in ref_df.columns else {}
    )
    priority = {}
    meta_lookup: dict[str, dict] = {}
    if reference_metadata is not None and not reference_metadata.empty:
        priority = reference_metadata.set_index("transcript_id")["reference_priority"].to_dict()
        meta_lookup = reference_metadata.set_index("transcript_id").to_dict(orient="index")
    parents = _select_parent(
        _shared_junction_counts(iso_jx, ref_jx),
        _exon_overlap_bp(isoforms, reference),
        cds_tx_ids=cds_tx_ids,
        prefer_coding_parent=prefer_coding_parent,
        iso_junction_totals=(
            iso_jx.df.groupby("transcript_id").size().to_dict() if len(iso_jx) else {}
        ),
        ref_junction_totals=(
            ref_jx.df.groupby("transcript_id").size().to_dict() if len(ref_jx) else {}
        ),
        iso_lengths=(iso_df.assign(_length=iso_df["End"] - iso_df["Start"])
                     .groupby("transcript_id")["_length"].sum().to_dict()),
        ref_lengths=(ref_df.assign(_length=ref_df["End"] - ref_df["Start"])
                     .groupby("transcript_id")["_length"].sum().to_dict()),
        iso_gene_ids=iso_gene_ids,
        ref_gene_ids=ref_gene_ids,
        reference_priority=priority,
        parent_hints=parent_hints,
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
                    "parent_candidate_count": 0,
                    "parent_ambiguous": False,
                    "parent_match_score": 0.0,
                    "parent_match_margin": None,
                    "parent_selection_reason": "no_overlap",
                    "junction_precision": 0.0,
                    "junction_recall": 0.0,
                    "junction_f1": 0.0,
                    "exact_intron_chain": False,
                    "iso_chain_contained": False,
                    "parent_reference_priority": 0.0,
                    "reference_cds_complete": False,
                    "reference_has_explicit_start": False,
                    "reference_has_explicit_stop": False,
                    "reference_cds_phase_valid": False,
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
                "parent_candidate_count": int(match["parent_candidate_count"]),
                "parent_ambiguous": bool(match["parent_ambiguous"]),
                "parent_match_score": float(match["parent_match_score"]),
                "parent_match_margin": (
                    float(match["parent_match_margin"])
                    if pd.notna(match["parent_match_margin"]) else None
                ),
                "parent_selection_reason": str(match["parent_selection_reason"]),
                "junction_precision": float(match["junction_precision"]),
                "junction_recall": float(match["junction_recall"]),
                "junction_f1": float(match["junction_f1"]),
                "exact_intron_chain": bool(match["exact_intron_chain"]),
                "iso_chain_contained": bool(match["iso_chain_contained"]),
                "parent_reference_priority": float(match["reference_priority"]),
                "reference_cds_complete": bool(
                    meta_lookup.get(ref_tx, {}).get("reference_cds_complete", True)
                ),
                "reference_has_explicit_start": bool(
                    meta_lookup.get(ref_tx, {}).get("reference_has_explicit_start", False)
                ),
                "reference_has_explicit_stop": bool(
                    meta_lookup.get(ref_tx, {}).get("reference_has_explicit_stop", False)
                ),
                "reference_cds_phase_valid": bool(
                    meta_lookup.get(ref_tx, {}).get("reference_cds_phase_valid", True)
                ),
            }
        )

    per_tx = pd.DataFrame(rows)
    annotated = isoforms.df.merge(per_tx, on="transcript_id", how="left")
    return pr.PyRanges(annotated)
