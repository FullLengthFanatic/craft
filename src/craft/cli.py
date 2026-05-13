"""Click-based command-line interface for CRAFT."""

from pathlib import Path

import click


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
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Optional per-cell isoform count matrix (h5ad / MTX) for sc workflows.",
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
    output_dir: Path,
) -> None:
    """Annotate isoforms with functional consequences (ORF, NMD, Pfam, 3' UTR)."""
    raise NotImplementedError("annotate command not yet implemented")
