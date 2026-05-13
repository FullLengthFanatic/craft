"""Per-cell count matrix loaders.

Two input shapes are supported:

* A single ``.h5ad`` file: read directly via :func:`anndata.read_h5ad`.
* A directory in 10x-style MTX layout: ``matrix.mtx[.gz]`` plus
  ``barcodes.tsv[.gz]`` plus ``features.tsv[.gz]`` or ``genes.tsv[.gz]``.
  The MTX is read with :func:`anndata.io.read_mtx`, transposed to put cells in
  ``obs`` and features in ``var``, and labelled with the barcode and feature
  IDs (only the first whitespace-separated column of the TSV is used, so
  10x's ``id, name, type`` files work as-is).
"""

import gzip
from pathlib import Path

import anndata as ad
import anndata.io
import pandas as pd


def _read_first_column(path: Path) -> list[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return [line.rstrip("\n").split("\t")[0] for line in fh if line.strip()]


def _find_first(directory: Path, candidates: list[str]) -> Path:
    for name in candidates:
        p = directory / name
        if p.exists():
            return p
    raise FileNotFoundError(
        f"Missing required file in {directory}; expected one of: {candidates}"
    )


def _load_10x_mtx(directory: Path) -> ad.AnnData:
    matrix_path = _find_first(directory, ["matrix.mtx.gz", "matrix.mtx"])
    barcodes_path = _find_first(directory, ["barcodes.tsv.gz", "barcodes.tsv"])
    features_path = _find_first(
        directory,
        ["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"],
    )

    barcodes = _read_first_column(barcodes_path)
    features = _read_first_column(features_path)

    raw = anndata.io.read_mtx(matrix_path)
    if raw.shape == (len(features), len(barcodes)):
        adata = raw.T
    elif raw.shape == (len(barcodes), len(features)):
        adata = raw
    else:
        raise ValueError(
            f"MTX shape {raw.shape} does not match barcodes ({len(barcodes)}) "
            f"x features ({len(features)})"
        )

    adata.obs_names = pd.Index(barcodes)
    adata.var_names = pd.Index(features)
    return adata


def load_counts(path: Path) -> ad.AnnData:
    """Load per-cell isoform counts from h5ad or 10x-style MTX.

    Args:
        path: Either a ``.h5ad`` file or a directory containing 10x-style MTX
            files (matrix.mtx[.gz], barcodes.tsv[.gz], and features.tsv[.gz]
            or genes.tsv[.gz]).

    Returns:
        AnnData with cells in ``obs`` and isoforms in ``var``.
    """
    path = Path(path)
    if path.is_file() and path.suffix == ".h5ad":
        return ad.read_h5ad(path)
    if path.is_dir():
        return _load_10x_mtx(path)
    raise ValueError(
        f"Unsupported counts path: {path}. "
        "Expected a .h5ad file or a directory with 10x-style MTX files."
    )
