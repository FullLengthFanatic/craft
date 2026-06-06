"""Click-based command-line interface for CRAFT."""

from pathlib import Path

import click

from craft.pipeline import run_annotate


@click.group()
@click.version_option(package_name="craft")
def cli() -> None:
    """CRAFT: long-read isoform functional-consequence annotator."""


@cli.command()
@click.option(
    "--isoforms",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Isoform GTF from FLAIR / IsoQuant / Bambu / FLAMES / SQANTI3.",
)
@click.option(
    "--reference",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Reference annotation GTF with CDS records (GENCODE or Ensembl).",
)
@click.option(
    "--genome",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Genome FASTA (indexed).",
)
@click.option(
    "--counts",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="Optional per-cell isoform count matrix (h5ad / MTX) for sc workflows.",
)
@click.option(
    "--pfam-hmm",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional Pfam-A.hmm file. When provided, isoforms are scanned for "
    "Pfam domains and compared to the parent transcript's domain set.",
)
@click.option(
    "--polya-atlas",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional BED file of polyadenylation sites (e.g. PolyASite v3.0 or "
    "PolyA_DB v4). When provided, atlas matches drive ALT_3PRIME_END / "
    "STOP_AT_ALT_POLYA reclassification; the canonical poly(A) motif scan stays "
    "as the fallback. See docs/user_guide.md for the expected BED format.",
)
@click.option(
    "--output-dir",
    type=click.Path(file_okay=False, path_type=Path),
    required=True,
    help="Output directory (created if missing).",
)
@click.option(
    "--group-by",
    type=str,
    default=None,
    help="Column in the counts obs to aggregate functional consequences by "
    "(e.g. cell_type, leiden). Writes per_celltype_consequence.tsv. Requires --counts.",
)
@click.option(
    "--tolerance",
    type=int,
    default=50,
    show_default=True,
    help="Slack (bp) on each transcript end before calling it truncated.",
)
@click.option(
    "--ptc-threshold-nt",
    type=int,
    default=50,
    show_default=True,
    help="PTC rule: stop more than this many mRNA-nt upstream of the last "
    "exon-exon junction is NMD-sensitive.",
)
@click.option(
    "--start-proximal-nt",
    type=int,
    default=150,
    show_default=True,
    help="CDS shorter than this (bp) escapes NMD via re-initiation.",
)
@click.option(
    "--long-last-exon-nt",
    type=int,
    default=400,
    show_default=True,
    help="Last exon longer than this (bp) escapes NMD.",
)
@click.option(
    "--min-orf-aa",
    type=int,
    default=50,
    show_default=True,
    help="Minimum de novo ORF length in amino acids.",
)
@click.option(
    "--orf-high-confidence",
    type=float,
    default=0.85,
    show_default=True,
    help="ORF confidence score at or above this is HIGH.",
)
@click.option(
    "--orf-medium-confidence",
    type=float,
    default=0.5,
    show_default=True,
    help="ORF confidence score at or above this (and below high) is MEDIUM.",
)
@click.option(
    "--long-utr3-nt",
    type=int,
    default=1000,
    show_default=True,
    help="Resolved 3'UTR longer than this sets long_utr3_triggers_nmd.",
)
@click.option(
    "--prefer-coding-parent",
    is_flag=True,
    default=False,
    help="On ties (equal shared junctions and exon overlap), prefer a CDS-bearing "
    "reference parent. Off by default to keep parent selection reproducible.",
)
@click.option(
    "--coding-potential/--no-coding-potential",
    default=True,
    show_default=True,
    help="Score each isoform's ORF for coding potential, using a hexamer + ORF "
    "model self-calibrated to the reference. Skipped if the reference has no "
    "non-coding transcripts.",
)
@click.option(
    "--classification",
    "classification_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional SQANTI3/pigeon (or any) classification TSV, keyed by isoform "
    "id. Selected columns are joined onto the per-isoform output by transcript_id "
    "(e.g. to carry structural_category for novel-boundary x consequence analysis).",
)
@click.option(
    "--classification-columns",
    default="structural_category",
    show_default=True,
    help="Comma-separated column names to carry from --classification.",
)
def annotate(
    isoforms: Path,
    reference: Path,
    genome: Path,
    counts: Path | None,
    pfam_hmm: Path | None,
    polya_atlas: Path | None,
    output_dir: Path,
    group_by: str | None,
    tolerance: int,
    ptc_threshold_nt: int,
    start_proximal_nt: int,
    long_last_exon_nt: int,
    min_orf_aa: int,
    orf_high_confidence: float,
    orf_medium_confidence: float,
    long_utr3_nt: int,
    prefer_coding_parent: bool,
    coding_potential: bool,
    classification_path: Path | None,
    classification_columns: str,
) -> None:
    """Annotate isoforms with functional consequences (ORF, NMD, Pfam, 3' UTR)."""
    result = run_annotate(
        isoforms_path=isoforms,
        reference_path=reference,
        output_dir=output_dir,
        genome_path=genome,
        counts_path=counts,
        pfam_hmm_path=pfam_hmm,
        polya_atlas_path=polya_atlas,
        tolerance=tolerance,
        ptc_threshold_nt=ptc_threshold_nt,
        start_proximal_nt=start_proximal_nt,
        long_last_exon_nt=long_last_exon_nt,
        min_orf_aa=min_orf_aa,
        orf_high_confidence=orf_high_confidence,
        orf_medium_confidence=orf_medium_confidence,
        long_utr3_nt=long_utr3_nt,
        prefer_coding_parent=prefer_coding_parent,
        coding_potential=coding_potential,
        classification_path=classification_path,
        classification_columns=classification_columns,
        group_by=group_by,
    )
    click.echo(f"Annotated {len(result)} isoforms -> {output_dir}/")
