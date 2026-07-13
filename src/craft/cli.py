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
    "--cell-whitelist",
    "cell_whitelist_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional text file of called-cell barcodes (one per line). With --counts, "
    "total_count / n_cells_detected are computed over these cells only; otherwise "
    "every barcode in the matrix is used (includes ambient droplets).",
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
    "exon-exon junction is predicted NMD-susceptible.",
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
    help="Exon containing the PTC longer than this (bp) escapes NMD (long-exon rule).",
)
@click.option(
    "--min-orf-aa",
    type=int,
    default=50,
    show_default=True,
    help="Minimum de novo ORF length in amino acids.",
)
@click.option(
    "--infer-alternative-start",
    is_flag=True,
    default=False,
    help="For a 5'-censored CDS, optionally infer the first downstream in-frame ATG. "
    "Off by default because an internal ATG is not evidence of the biological start.",
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
    help="Add a small explicit preference for a CDS-bearing reference parent in "
    "otherwise close candidate rankings. Off by default.",
)
@click.option(
    "--coding-potential/--no-coding-potential",
    default=True,
    show_default=True,
    help="Score each isoform's ORF for coding potential, using a hexamer + ORF "
    "classifier trained on the reference. Skipped if the reference has no "
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
@click.option(
    "--evidence-table",
    "evidence_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional per-isoform molecule/read evidence TSV from scNoiseMeter, tecap, "
    "or a caller. Adds an explicit uncalibrated evidence score and tier; never filters rows.",
)
@click.option(
    "--orf-comparator-gtf",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional independent CDS-annotated GTF (for example ORFanage output). "
    "Reports start/stop/CDS agreement without allowing either caller to define truth.",
)
@click.option(
    "--recurrence-null",
    type=click.Choice(["none", "occupancy", "betabinom"]),
    default="none",
    show_default=True,
    help="Compare per-cell recurrence with an exploratory null (needs --counts): "
    "'occupancy' uses a depth-aware Poissonized occupancy null; 'betabinom' fits "
    "an empirical beta-binomial "
    "(stratified by structural_category when --classification is given). Emits "
    "recurrence_pvalue / recurrence_score. These are exploratory dispersion "
    "statistics, not probabilities that an isoform is real. 'none' leaves them empty.",
)
def annotate(
    isoforms: Path,
    reference: Path,
    genome: Path,
    counts: Path | None,
    cell_whitelist_path: Path | None,
    pfam_hmm: Path | None,
    polya_atlas: Path | None,
    output_dir: Path,
    group_by: str | None,
    tolerance: int,
    ptc_threshold_nt: int,
    start_proximal_nt: int,
    long_last_exon_nt: int,
    min_orf_aa: int,
    infer_alternative_start: bool,
    orf_high_confidence: float,
    orf_medium_confidence: float,
    long_utr3_nt: int,
    prefer_coding_parent: bool,
    coding_potential: bool,
    classification_path: Path | None,
    classification_columns: str,
    evidence_path: Path | None,
    orf_comparator_gtf: Path | None,
    recurrence_null: str,
) -> None:
    """Annotate isoforms with structure evidence, ORFs, surveillance, and UTRs."""
    result = run_annotate(
        isoforms_path=isoforms,
        reference_path=reference,
        output_dir=output_dir,
        genome_path=genome,
        counts_path=counts,
        cell_whitelist_path=cell_whitelist_path,
        pfam_hmm_path=pfam_hmm,
        polya_atlas_path=polya_atlas,
        tolerance=tolerance,
        ptc_threshold_nt=ptc_threshold_nt,
        start_proximal_nt=start_proximal_nt,
        long_last_exon_nt=long_last_exon_nt,
        min_orf_aa=min_orf_aa,
        infer_alternative_start=infer_alternative_start,
        orf_high_confidence=orf_high_confidence,
        orf_medium_confidence=orf_medium_confidence,
        long_utr3_nt=long_utr3_nt,
        prefer_coding_parent=prefer_coding_parent,
        coding_potential=coding_potential,
        classification_path=classification_path,
        classification_columns=classification_columns,
        evidence_path=evidence_path,
        orf_comparator_gtf=orf_comparator_gtf,
        group_by=group_by,
        recurrence_null=recurrence_null,
    )
    click.echo(f"Annotated {len(result)} isoforms -> {output_dir}/")
