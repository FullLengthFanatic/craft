"""End-to-end CRAFT annotation pipeline.

Orchestrates the existing core modules into a single ``run_annotate`` entry
point and emits per-isoform TSV and JSON outputs to a directory.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyranges as pr
import pysam

from craft.core.coding_potential import score_isoforms as coding_potential_score
from craft.core.completeness import Completeness, classify
from craft.core.nmd import predict as nmd_predict
from craft.core.nmd import predict_denovo as nmd_predict_denovo
from craft.core.nmd import predict_resolved as nmd_predict_resolved
from craft.core.orf.confidence import ORFConfidence, score
from craft.core.orf.denovo import predict as denovo_predict
from craft.core.orf.propagation import ORFOutcome, propagate
from craft.core.orf.resolve import resolve as resolve_orf
from craft.core.pfam import scan as pfam_scan
from craft.core.polya_atlas import build_atlas_index, load_atlas, match_iso_end
from craft.core.utr3 import annotate as utr3_annotate
from craft.core.utr3 import annotate_resolved as utr3_annotate_resolved
from craft.core.utr3 import polya_near_3prime_end
from craft.export.anndata import to_anndata, write_h5ad
from craft.export.celltype import aggregate_consequences
from craft.io.counts import load_counts
from craft.io.gtf import load_isoforms, load_reference
from craft.report.html import render as render_report

_LIST_COLUMNS = (
    "propagated_cds_intervals",
    "denovo_cds_intervals",
    "resolved_cds_intervals",
    "iso_pfam_domains",
    "parent_pfam_domains",
    "pfam_preserved",
    "pfam_lost",
    "pfam_gained",
)

_PFAM_COLUMNS = (
    "iso_pfam_domains",
    "parent_pfam_domains",
    "pfam_preserved",
    "pfam_lost",
    "pfam_gained",
)

_OUTPUT_COLUMNS = [
    "transcript_id",
    "completeness",
    "parent_tx_id",
    "parent_gene_id",
    "parent_gene_name",
    "shared_junctions",
    "parent_overlap_bp",
    "orf_outcome",
    "propagated_cds_bp",
    "parent_cds_bp",
    "start_codon_covered",
    "stop_codon_covered",
    "propagated_cds_intervals",
    "denovo_orf_found",
    "denovo_cds_bp",
    "denovo_orf_aa_length",
    "denovo_start_codon",
    "denovo_stop_codon",
    "denovo_cds_intervals",
    "orf_confidence",
    "orf_confidence_score",
    "nmd_status",
    "nmd_rule",
    "stop_to_last_junction_nt",
    "last_exon_length_nt",
    "nmd_confidence",
    "iso_utr3_length_nt",
    "parent_utr3_length_nt",
    "utr3_length_delta_nt",
    "utr3_length_delta_pct",
    "polya_signal_motif",
    "polya_signal_distance_nt",
    "polya_evidence_source",
    "polya_db_site_id",
    "iso_pfam_domains",
    "parent_pfam_domains",
    "pfam_preserved",
    "pfam_lost",
    "pfam_gained",
    # v1.5 additive columns (existing columns above are unchanged).
    "has_cds_bearing_parent",
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
    "nmd_status_resolved",
    "nmd_rule_resolved",
    "nmd_confidence_resolved",
    "iso_utr3_length_resolved_nt",
    "utr3_length_delta_resolved_nt",
    "utr3_length_delta_pct_resolved",
    "long_utr3_triggers_nmd",
    "iso_utr5_length_nt",
    "parent_utr5_length_nt",
    "utr5_length_delta_nt",
    "utr5_length_delta_pct",
    "nmd_status_denovo",
    "nmd_rule_denovo",
    "nmd_confidence_denovo",
    "coding_potential_score",
    "coding_potential_label",
    "coding_potential_orf_source",
]

_RESOLVE_COLUMNS = (
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
)

_DENOVO_TRIGGER_OUTCOMES = frozenset(
    {
        ORFOutcome.NO_PARENT.value,
        ORFOutcome.NO_PARENT_CDS.value,
        ORFOutcome.START_LOST.value,
    }
)


def _filter_isoforms_by_genome_contigs(
    isoforms: pr.PyRanges, genome_path: Path
) -> pr.PyRanges:
    """Drop isoforms whose chromosome is absent from the genome FASTA.

    PacBio collapse outputs frequently reference random/alt contigs
    (``chr1_KI270706v1_random``, ``chrUn_*``, etc.) that the GRCh38
    primary_assembly FASTA omits. Those isoforms have no parent annotation
    in standard GENCODE anyway (NOVEL_NO_MATCH downstream) and any FASTA
    fetch on them raises ``KeyError`` mid-run. Drop them up front with a
    one-line warning to stderr.
    """
    with pysam.FastaFile(str(genome_path)) as genome:
        available = set(genome.references)
    iso_df = isoforms.df
    iso_contigs = set(iso_df["Chromosome"].astype(str).unique())
    missing = iso_contigs - available
    if not missing:
        return isoforms
    mask = ~iso_df["Chromosome"].astype(str).isin(missing)
    keep = iso_df[mask].copy()
    dropped_tx = (
        iso_df.loc[~mask, "transcript_id"].nunique()
        if "transcript_id" in iso_df.columns
        else 0
    )
    sample = sorted(missing)[:5]
    suffix = ", ..." if len(missing) > 5 else ""
    print(
        f"[craft] Skipping {dropped_tx} isoforms on {len(missing)} contigs "
        f"not in the genome FASTA: {sample}{suffix}",
        file=sys.stderr,
    )
    if keep.empty:
        return pr.PyRanges()
    keep["Strand"] = keep["Strand"].astype(str)
    return pr.PyRanges(keep)


def _exon_only_reference(reference: pr.PyRanges) -> pr.PyRanges:
    df = reference.df
    exons_df = df[df["Feature"] == "exon"].drop(columns="Feature").copy()
    exons_df["Strand"] = exons_df["Strand"].astype(str)
    return pr.PyRanges(exons_df)


def _orphan_isoforms(isoforms: pr.PyRanges, tx_ids: list[str]) -> pr.PyRanges:
    df = isoforms.df
    orphans = df[df["transcript_id"].isin(tx_ids)].copy()
    orphans["Strand"] = orphans["Strand"].astype(str)
    return pr.PyRanges(orphans)


def _confidence_for(
    completeness_value: str,
    orf_outcome: str,
    denovo_found: bool,
    high_threshold: float,
    medium_threshold: float,
) -> tuple[str, float]:
    category, score_value = score(
        completeness_value, orf_outcome, high_threshold, medium_threshold
    )
    if category == ORFConfidence.NONE and denovo_found:
        return ORFConfidence.LOW.value, 0.25
    return category.value, score_value


def _iso_3prime_pos(exons: pd.DataFrame, strand: str) -> tuple[str, int]:
    """Genomic position of the iso's 3' end (transcript orientation)."""
    sorted_exons = exons.sort_values("Start").reset_index(drop=True)
    chrom = str(sorted_exons.iloc[0]["Chromosome"])
    if strand == "+":
        pos = int(sorted_exons.iloc[-1]["End"]) - 1
    elif strand == "-":
        pos = int(sorted_exons.iloc[0]["Start"])
    else:
        raise ValueError(f"Unsupported strand: {strand!r}")
    return chrom, pos


def _compute_polya_evidence(
    isoforms: pr.PyRanges,
    genome_path: Path,
    polya_atlas_path: Path | None = None,
) -> dict[str, dict]:
    """For every iso transcript_id, return poly(A) evidence:
    {found: bool, source: "polya_db"|"canonical_motif"|"none", pas_id: str}.

    Atlas matching runs first when a BED path is provided; misses (and the
    no-atlas case) fall back to the canonical motif scan in the iso's last 50 nt.
    """
    atlas_index = (
        build_atlas_index(load_atlas(polya_atlas_path))
        if polya_atlas_path is not None
        else {}
    )
    iso_df = isoforms.df
    iso_strand = iso_df.groupby("transcript_id")["Strand"].first().to_dict()
    iso_exons_by_tx = {tx: g for tx, g in iso_df.groupby("transcript_id", sort=False)}
    evidence: dict[str, dict] = {}
    with pysam.FastaFile(str(genome_path)) as genome:
        for tx_id, exons in iso_exons_by_tx.items():
            strand = str(iso_strand[tx_id])

            if atlas_index:
                chrom, iso_3p = _iso_3prime_pos(exons, strand)
                hit = match_iso_end(iso_3p, chrom, strand, atlas_index)
                if hit["matched"]:
                    evidence[str(tx_id)] = {
                        "found": True,
                        "source": "polya_db",
                        "pas_id": str(hit["pas_id"]),
                    }
                    continue

            sig = polya_near_3prime_end(exons, strand, genome)
            if sig["found"]:
                evidence[str(tx_id)] = {
                    "found": True,
                    "source": "canonical_motif",
                    "pas_id": "",
                }
            else:
                evidence[str(tx_id)] = {
                    "found": False,
                    "source": "none",
                    "pas_id": "",
                }
    return evidence


def _reclassify_with_polya(
    classified: pr.PyRanges,
    propagated: pd.DataFrame,
    polya_evidence: dict[str, dict],
) -> tuple[pr.PyRanges, pd.DataFrame]:
    """Use poly(A) evidence (atlas hit OR canonical motif) to relabel
    TRUNCATED_3P -> ALT_3PRIME_END and STOP_NOT_OBSERVED -> STOP_AT_ALT_POLYA."""

    def _found(tx: str) -> bool:
        return bool(polya_evidence.get(tx, {}).get("found", False))

    cls_df = classified.df.copy()

    def _new_completeness(row: pd.Series) -> str:
        if row["completeness"] == Completeness.TRUNCATED_3P.value and _found(
            str(row["transcript_id"])
        ):
            return Completeness.ALT_3PRIME_END.value
        return str(row["completeness"])

    cls_df["completeness"] = cls_df.apply(_new_completeness, axis=1)
    classified = pr.PyRanges(cls_df)

    propagated = propagated.copy()

    def _new_outcome(row: pd.Series) -> str:
        if row["orf_outcome"] == ORFOutcome.STOP_NOT_OBSERVED.value and _found(
            str(row["transcript_id"])
        ):
            return ORFOutcome.STOP_AT_ALT_POLYA.value
        return str(row["orf_outcome"])

    propagated["orf_outcome"] = propagated.apply(_new_outcome, axis=1)
    return classified, propagated


def _polya_evidence_columns(polya_evidence: dict[str, dict]) -> pd.DataFrame:
    """Build a DataFrame of (transcript_id, polya_evidence_source, polya_db_site_id)
    rows ready to merge into the per-isoform output."""
    if not polya_evidence:
        return pd.DataFrame(
            columns=["transcript_id", "polya_evidence_source", "polya_db_site_id"]
        )
    return pd.DataFrame(
        [
            {
                "transcript_id": tx,
                "polya_evidence_source": info.get("source", "none"),
                "polya_db_site_id": info.get("pas_id", ""),
            }
            for tx, info in polya_evidence.items()
        ]
    )


def _write_coding_model(model: dict, path: Path) -> None:
    """Write a small, readable summary of the coding-potential model."""
    summary = {
        "features": ["hexamer_llr", "log10_orf_len", "orf_coverage", "fickett"],
        "weights": [float(w) for w in model["weights"][:-1]],
        "intercept": float(model["weights"][-1]),
        "feature_mean": [float(v) for v in model["mean"]],
        "feature_std": [float(v) for v in model["std"]],
        "threshold": float(model["threshold"]),
        "n_train_coding": int(model["n_coding"]),
        "n_train_noncoding": int(model["n_noncoding"]),
        "heldout_auc": None if np.isnan(model["heldout_auc"]) else float(model["heldout_auc"]),
    }
    with open(path, "w") as fh:
        json.dump(summary, fh, indent=2)


def _fill_resolved_defaults(df: pd.DataFrame) -> None:
    """Fill defaults for the resolved columns in place (covers the no-genome path)."""
    str_defaults = {
        "resolved_orf_status": "resolution_failed",
        "nmd_status_resolved": "not_applicable",
        "nmd_rule_resolved": "",
        "nmd_confidence_resolved": "none",
        "nmd_status_denovo": "not_applicable",
        "nmd_rule_denovo": "",
        "nmd_confidence_denovo": "none",
    }
    for col, default in str_defaults.items():
        if col in df.columns:
            df[col] = df[col].fillna(default)
    for col in (
        "ptc_introduced",
        "intron_retained_in_cds",
        "frame_consistent",
        "stop_in_transcript",
        "uorf_triggers_nmd",
        "long_utr3_triggers_nmd",
    ):
        if col in df.columns:
            df[col] = df[col].fillna(False).astype(bool)
    for col in ("resolved_cds_bp", "resolved_aa_length", "uorf_count"):
        if col in df.columns:
            df[col] = df[col].fillna(0).astype("int64")
    if "resolved_cds_intervals" in df.columns:
        df["resolved_cds_intervals"] = df["resolved_cds_intervals"].apply(
            lambda v: v if isinstance(v, list) else []
        )


def _empty_denovo() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "transcript_id",
            "denovo_orf_found",
            "denovo_cds_bp",
            "denovo_cds_intervals",
            "denovo_orf_aa_length",
            "denovo_start_codon",
            "denovo_stop_codon",
        ]
    )


def run_annotate(
    isoforms_path: Path,
    reference_path: Path,
    output_dir: Path,
    genome_path: Path | None = None,
    counts_path: Path | None = None,
    pfam_hmm_path: Path | None = None,
    polya_atlas_path: Path | None = None,
    tolerance: int = 50,
    ptc_threshold_nt: int = 50,
    start_proximal_nt: int = 150,
    long_last_exon_nt: int = 400,
    min_orf_aa: int = 50,
    orf_high_confidence: float = 0.85,
    orf_medium_confidence: float = 0.5,
    long_utr3_nt: int = 1000,
    prefer_coding_parent: bool = False,
    coding_potential: bool = True,
    group_by: str | None = None,
) -> pd.DataFrame:
    """Run the full CRAFT annotation pipeline.

    Args:
        isoforms_path: GTF of isoform exons from FLAIR / IsoQuant / Bambu / etc.
        reference_path: GTF of reference annotation with both exon and CDS records.
        output_dir: Directory for ``per_isoform.tsv`` and ``per_isoform.json``.
        genome_path: Optional indexed genome FASTA. Required for de novo ORF
            prediction on orphans and for poly(A) signal scanning.
        counts_path: Optional per-cell count matrix (h5ad or 10x MTX dir).
            Currently parsed elsewhere; not used in this pipeline yet.

    Returns:
        Per-isoform DataFrame with the columns written to ``per_isoform.tsv``.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    isoforms = load_isoforms(isoforms_path)
    if genome_path is not None and len(isoforms) > 0:
        isoforms = _filter_isoforms_by_genome_contigs(isoforms, genome_path)
    reference = load_reference(reference_path)

    ref_df = reference.df
    cds_tx_ids = set(
        ref_df.loc[ref_df["Feature"] == "CDS", "transcript_id"].astype(str).unique()
    )
    classified = classify(
        isoforms,
        _exon_only_reference(reference),
        tolerance=tolerance,
        cds_tx_ids=cds_tx_ids,
        prefer_coding_parent=prefer_coding_parent,
    )
    propagated = propagate(classified, reference)

    polya_evidence: dict[str, dict] = {}
    if genome_path is not None:
        polya_evidence = _compute_polya_evidence(
            isoforms, genome_path, polya_atlas_path=polya_atlas_path
        )
        classified, propagated = _reclassify_with_polya(
            classified, propagated, polya_evidence
        )

    orphan_ids = propagated.loc[
        propagated["orf_outcome"].isin(_DENOVO_TRIGGER_OUTCOMES),
        "transcript_id",
    ].tolist()
    if genome_path is not None and orphan_ids:
        denovo_df = denovo_predict(
            _orphan_isoforms(isoforms, orphan_ids), genome_path, min_orf_aa=min_orf_aa
        )
    else:
        denovo_df = _empty_denovo()

    nmd_df = nmd_predict(
        classified,
        propagated,
        ptc_threshold_nt=ptc_threshold_nt,
        start_proximal_nt=start_proximal_nt,
        long_last_exon_nt=long_last_exon_nt,
    )
    utr3_df = utr3_annotate(classified, propagated, reference, genome_fasta=genome_path)

    # Sequence-level ORF resolution + its resolved consumers (additive columns).
    if genome_path is not None:
        resolve_df = resolve_orf(
            classified, propagated, reference, genome_path, ptc_threshold_nt=ptc_threshold_nt
        )
        nmd_resolved_df = nmd_predict_resolved(
            classified,
            resolve_df,
            ptc_threshold_nt=ptc_threshold_nt,
            start_proximal_nt=start_proximal_nt,
            long_last_exon_nt=long_last_exon_nt,
        )
        utr3_resolved_df = utr3_annotate_resolved(
            classified, resolve_df, reference, long_utr3_nt=long_utr3_nt
        )
    else:
        resolve_df = pd.DataFrame(columns=["transcript_id", *_RESOLVE_COLUMNS])
        nmd_resolved_df = pd.DataFrame(
            columns=[
                "transcript_id",
                "nmd_status_resolved",
                "nmd_rule_resolved",
                "nmd_confidence_resolved",
            ]
        )
        utr3_resolved_df = pd.DataFrame(
            columns=[
                "transcript_id",
                "iso_utr3_length_resolved_nt",
                "utr3_length_delta_resolved_nt",
                "utr3_length_delta_pct_resolved",
                "long_utr3_triggers_nmd",
                "iso_utr5_length_nt",
                "parent_utr5_length_nt",
                "utr5_length_delta_nt",
                "utr5_length_delta_pct",
            ]
        )

    # NMD on the de novo ORF, for orphan isoforms with no reference-anchored stop.
    nmd_denovo_df = nmd_predict_denovo(
        classified,
        denovo_df,
        ptc_threshold_nt=ptc_threshold_nt,
        start_proximal_nt=start_proximal_nt,
        long_last_exon_nt=long_last_exon_nt,
    )

    per_tx = classified.df.groupby("transcript_id").first().reset_index()[
        ["transcript_id", "completeness"]
    ]

    merged = propagated.merge(per_tx, on="transcript_id", how="left")
    if not denovo_df.empty:
        merged = merged.merge(denovo_df, on="transcript_id", how="left")
    else:
        n = len(merged)
        merged["denovo_orf_found"] = [False] * n
        merged["denovo_cds_bp"] = [0] * n
        merged["denovo_orf_aa_length"] = [0] * n
        merged["denovo_cds_intervals"] = [[] for _ in range(n)]
        merged["denovo_start_codon"] = [""] * n
        merged["denovo_stop_codon"] = [""] * n
    merged = merged.merge(nmd_df, on="transcript_id", how="left")
    merged = merged.merge(utr3_df, on="transcript_id", how="left")
    merged = merged.merge(resolve_df, on="transcript_id", how="left")
    merged = merged.merge(nmd_resolved_df, on="transcript_id", how="left")
    merged = merged.merge(utr3_resolved_df, on="transcript_id", how="left")
    merged = merged.merge(nmd_denovo_df, on="transcript_id", how="left")
    _fill_resolved_defaults(merged)

    completeness_default = Completeness.NOVEL_NO_MATCH.value
    merged["completeness"] = merged["completeness"].fillna(completeness_default)
    merged["denovo_orf_found"] = merged["denovo_orf_found"].map(
        lambda v: bool(v) if pd.notna(v) else False
    )

    confidence = merged.apply(
        lambda r: _confidence_for(
            r["completeness"],
            r["orf_outcome"],
            bool(r["denovo_orf_found"]),
            orf_high_confidence,
            orf_medium_confidence,
        ),
        axis=1,
        result_type="expand",
    )
    confidence.columns = ["orf_confidence", "orf_confidence_score"]
    merged = pd.concat([merged, confidence], axis=1)

    # Add classify-derived columns from per-exon rows
    classify_meta = classified.df.groupby("transcript_id").first().reset_index()[
        ["transcript_id", "shared_junctions", "parent_overlap_bp", "has_cds_bearing_parent"]
    ]
    merged = merged.merge(classify_meta, on="transcript_id", how="left")
    merged["has_cds_bearing_parent"] = merged["has_cds_bearing_parent"].fillna(False).astype(bool)

    # Look up parent gene from reference for the report's notable-findings panel.
    ref_df = reference.df
    gene_id_map: dict[str, str] = {}
    gene_name_map: dict[str, str] = {}
    if "gene_id" in ref_df.columns:
        gene_id_map = (
            ref_df[["transcript_id", "gene_id"]]
            .dropna()
            .drop_duplicates()
            .set_index("transcript_id")["gene_id"]
            .to_dict()
        )
    if "gene_name" in ref_df.columns:
        gene_name_map = (
            ref_df[["transcript_id", "gene_name"]]
            .dropna()
            .drop_duplicates()
            .set_index("transcript_id")["gene_name"]
            .to_dict()
        )
    merged["parent_gene_id"] = merged["parent_tx_id"].map(gene_id_map).fillna("")
    merged["parent_gene_name"] = merged["parent_tx_id"].map(gene_name_map).fillna("")

    evidence_df = _polya_evidence_columns(polya_evidence)
    if not evidence_df.empty:
        merged = merged.merge(evidence_df, on="transcript_id", how="left")
    if "polya_evidence_source" not in merged.columns:
        merged["polya_evidence_source"] = "none"
    if "polya_db_site_id" not in merged.columns:
        merged["polya_db_site_id"] = ""
    merged["polya_evidence_source"] = merged["polya_evidence_source"].fillna("none")
    merged["polya_db_site_id"] = merged["polya_db_site_id"].fillna("")

    if pfam_hmm_path is not None and genome_path is not None:
        pfam_df = pfam_scan(merged, reference, pfam_hmm_path, genome_path)
        merged = merged.merge(pfam_df, on="transcript_id", how="left")
    for col in _PFAM_COLUMNS:
        if col not in merged.columns:
            merged[col] = [[] for _ in range(len(merged))]
        else:
            merged[col] = merged[col].apply(lambda v: v if isinstance(v, list) else [])

    if coding_potential and genome_path is not None:
        cp_df, cp_model = coding_potential_score(merged, classified, reference, genome_path)
        merged = merged.merge(cp_df, on="transcript_id", how="left")
        if cp_model is not None:
            _write_coding_model(cp_model, output_dir / "coding_potential_model.json")
        else:
            print(
                "[craft] Reference has no non-coding transcripts; coding-potential "
                "score skipped (columns left empty).",
                file=sys.stderr,
            )
    if "coding_potential_score" not in merged.columns:
        merged["coding_potential_score"] = float("nan")
        merged["coding_potential_label"] = ""
        merged["coding_potential_orf_source"] = "none"

    merged = merged.reindex(columns=_OUTPUT_COLUMNS)

    counts_adata = load_counts(counts_path) if counts_path is not None else None
    adata = to_anndata(merged, counts=counts_adata)

    if group_by is not None and counts_adata is not None:
        try:
            aggregate_consequences(
                adata, merged, group_by, output_dir / "per_celltype_consequence.tsv"
            )
        except ValueError as exc:
            print(f"[craft] Skipping per-cell-type aggregation: {exc}", file=sys.stderr)
    elif group_by is not None:
        print(
            "[craft] --group-by given but no --counts; skipping per-cell-type aggregation.",
            file=sys.stderr,
        )

    write_h5ad(adata, output_dir / "annotated.h5ad")

    _write_outputs(merged, output_dir)
    return merged


def _write_outputs(merged: pd.DataFrame, output_dir: Path) -> None:
    tsv_df = merged.copy()
    for col in _LIST_COLUMNS:
        if col in tsv_df.columns:
            tsv_df[col] = tsv_df[col].apply(
                lambda v: json.dumps(v) if isinstance(v, list) else ("[]" if v is None else str(v))
            )
    tsv_df.to_csv(output_dir / "per_isoform.tsv", sep="\t", index=False)

    records = merged.to_dict(orient="records")
    with open(output_dir / "per_isoform.json", "w") as fh:
        json.dump(records, fh, default=str, indent=2)

    render_report(merged, output_dir / "report.html")
