"""End-to-end CRAFT annotation pipeline.

Orchestrates the existing core modules into a single ``run_annotate`` entry
point and emits per-isoform TSV and JSON outputs to a directory.
"""

import json
from pathlib import Path

import pandas as pd
import pyranges as pr

from craft.core.completeness import Completeness, classify
from craft.core.nmd import predict as nmd_predict
from craft.core.orf.confidence import ORFConfidence, score
from craft.core.orf.denovo import predict as denovo_predict
from craft.core.orf.propagation import ORFOutcome, propagate
from craft.core.utr3 import annotate as utr3_annotate
from craft.export.anndata import to_anndata, write_h5ad
from craft.io.counts import load_counts
from craft.io.gtf import load_isoforms, load_reference
from craft.report.html import render as render_report

_LIST_COLUMNS = ("propagated_cds_intervals", "denovo_cds_intervals")

_OUTPUT_COLUMNS = [
    "transcript_id",
    "completeness",
    "parent_tx_id",
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
        for col in _empty_denovo().columns:
            if col == "transcript_id":
                continue
            merged[col] = (
                False
                if col == "denovo_orf_found"
                else 0 if col in {"denovo_cds_bp", "denovo_orf_aa_length"}
                else [] if col == "denovo_cds_intervals"
                else ""
            )
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
