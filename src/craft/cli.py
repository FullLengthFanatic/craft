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
    "--output-dir",
    type=click.Path(file_okay=False, path_type=Path),
    required=True,
    help="Output directory (created if missing).",
)
def annotate(
    isoforms: Path,
    reference: Path,
    genome: Path,
    counts: Path | None,
    pfam_hmm: Path | None,
    output_dir: Path,
) -> None:
    """Annotate isoforms with functional consequences (ORF, NMD, Pfam, 3' UTR)."""
    result = run_annotate(
        isoforms_path=isoforms,
        reference_path=reference,
        output_dir=output_dir,
        genome_path=genome,
        counts_path=counts,
        pfam_hmm_path=pfam_hmm,
    )
    click.echo(f"Annotated {len(result)} isoforms -> {output_dir}/")
