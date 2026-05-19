"""End-to-end CRAFT annotation pipeline.

Orchestrates the existing core modules into a single ``run_annotate`` entry
point and emits per-isoform TSV and JSON outputs to a directory.
"""

import json
from pathlib import Path

import pandas as pd
import pyranges as pr
import pysam

from craft.core.completeness import Completeness, classify
from craft.core.nmd import predict as nmd_predict
from craft.core.orf.confidence import ORFConfidence, score
from craft.core.orf.denovo import predict as denovo_predict
from craft.core.orf.propagation import ORFOutcome, propagate
from craft.core.pfam import scan as pfam_scan
from craft.core.utr3 import annotate as utr3_annotate
from craft.core.utr3 import polya_near_3prime_end
from craft.export.anndata import to_anndata, write_h5ad
from craft.io.counts import load_counts
from craft.io.gtf import load_isoforms, load_reference
from craft.report.html import render as render_report

_LIST_COLUMNS = (
    "propagated_cds_intervals",
    "denovo_cds_intervals",
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
    "iso_pfam_domains",
    "parent_pfam_domains",
    "pfam_preserved",
    "pfam_lost",
    "pfam_gained",
]

_DENOVO_TRIGGER_OUTCOMES = frozenset(
    {
        ORFOutcome.NO_PARENT.value,
        ORFOutcome.NO_PARENT_CDS.value,
        ORFOutcome.START_LOST.value,
    }
)


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
) -> tuple[str, float]:
    category, score_value = score(completeness_value, orf_outcome)
    if category == ORFConfidence.NONE and denovo_found:
        return ORFConfidence.LOW.value, 0.25
    return category.value, score_value


def _compute_polya_evidence(
    isoforms: pr.PyRanges, genome_path: Path
) -> dict[str, bool]:
    """For every iso transcript_id, True iff a canonical poly(A) signal sits
    in the last 50 nt of its 3' end (transcript orientation)."""
    iso_df = isoforms.df
    iso_strand = iso_df.groupby("transcript_id")["Strand"].first().to_dict()
    iso_exons_by_tx = {tx: g for tx, g in iso_df.groupby("transcript_id", sort=False)}
    evidence: dict[str, bool] = {}
    with pysam.FastaFile(str(genome_path)) as genome:
        for tx_id, exons in iso_exons_by_tx.items():
            sig = polya_near_3prime_end(exons, str(iso_strand[tx_id]), genome)
            evidence[str(tx_id)] = bool(sig["found"])
    return evidence


def _reclassify_with_polya(
    classified: pr.PyRanges,
    propagated: pd.DataFrame,
    polya_evidence: dict[str, bool],
) -> tuple[pr.PyRanges, pd.DataFrame]:
    """Use poly(A) signal evidence to relabel APA isoforms previously called
    truncations: TRUNCATED_3P -> ALT_3PRIME_END, STOP_NOT_OBSERVED -> STOP_AT_ALT_POLYA."""

    cls_df = classified.df.copy()

    def _new_completeness(row: pd.Series) -> str:
        if (
            row["completeness"] == Completeness.TRUNCATED_3P.value
            and polya_evidence.get(str(row["transcript_id"]), False)
        ):
            return Completeness.ALT_3PRIME_END.value
        return str(row["completeness"])

    cls_df["completeness"] = cls_df.apply(_new_completeness, axis=1)
    classified = pr.PyRanges(cls_df)

    propagated = propagated.copy()

    def _new_outcome(row: pd.Series) -> str:
        if (
            row["orf_outcome"] == ORFOutcome.STOP_NOT_OBSERVED.value
            and polya_evidence.get(str(row["transcript_id"]), False)
        ):
            return ORFOutcome.STOP_AT_ALT_POLYA.value
        return str(row["orf_outcome"])

    propagated["orf_outcome"] = propagated.apply(_new_outcome, axis=1)
    return classified, propagated


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
    reference = load_reference(reference_path)

    classified = classify(isoforms, _exon_only_reference(reference))
    propagated = propagate(classified, reference)

    if genome_path is not None:
        polya_evidence = _compute_polya_evidence(isoforms, genome_path)
        classified, propagated = _reclassify_with_polya(
            classified, propagated, polya_evidence
        )

    orphan_ids = propagated.loc[
        propagated["orf_outcome"].isin(_DENOVO_TRIGGER_OUTCOMES),
        "transcript_id",
    ].tolist()
    if genome_path is not None and orphan_ids:
        denovo_df = denovo_predict(_orphan_isoforms(isoforms, orphan_ids), genome_path)
    else:
        denovo_df = _empty_denovo()

    nmd_df = nmd_predict(classified, propagated)
    utr3_df = utr3_annotate(classified, propagated, reference, genome_fasta=genome_path)

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

    completeness_default = Completeness.NOVEL_NO_MATCH.value
    merged["completeness"] = merged["completeness"].fillna(completeness_default)
    merged["denovo_orf_found"] = merged["denovo_orf_found"].map(
        lambda v: bool(v) if pd.notna(v) else False
    )

    confidence = merged.apply(
        lambda r: _confidence_for(
            r["completeness"], r["orf_outcome"], bool(r["denovo_orf_found"])
        ),
        axis=1,
        result_type="expand",
    )
    confidence.columns = ["orf_confidence", "orf_confidence_score"]
    merged = pd.concat([merged, confidence], axis=1)

    # Add classify-derived columns from per-exon rows
    classify_meta = classified.df.groupby("transcript_id").first().reset_index()[
        ["transcript_id", "shared_junctions", "parent_overlap_bp"]
    ]
    merged = merged.merge(classify_meta, on="transcript_id", how="left")

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

    if pfam_hmm_path is not None and genome_path is not None:
        pfam_df = pfam_scan(merged, reference, pfam_hmm_path, genome_path)
        merged = merged.merge(pfam_df, on="transcript_id", how="left")
    for col in _PFAM_COLUMNS:
        if col not in merged.columns:
            merged[col] = [[] for _ in range(len(merged))]
        else:
            merged[col] = merged[col].apply(lambda v: v if isinstance(v, list) else [])

    merged = merged.reindex(columns=_OUTPUT_COLUMNS)

    counts_adata = load_counts(counts_path) if counts_path is not None else None
    adata = to_anndata(merged, counts=counts_adata)
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
