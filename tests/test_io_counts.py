"""Tests for craft.io.counts."""

import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from craft.io.counts import load_counts

# 3 features x 2 cells, 4 non-zero entries:
# features (rows): iso1, iso2, iso3
# cells   (cols): cell1, cell2
# cell1: iso1=5, iso2=3
# cell2: iso1=10, iso3=7
_MTX_TEXT = (
    "%%MatrixMarket matrix coordinate integer general\n"
    "%\n"
    "3 2 4\n"
    "1 1 5\n"
    "2 1 3\n"
    "1 2 10\n"
    "3 2 7\n"
)


def _write_mtx_dir(
    directory: Path,
    *,
    features_name: str = "features.tsv",
    matrix_name: str = "matrix.mtx",
    barcodes_name: str = "barcodes.tsv",
    gzip_files: bool = False,
) -> None:
    def _put(name: str, content: str) -> None:
        path = directory / name
        if name.endswith(".gz") or gzip_files:
            out = path if name.endswith(".gz") else path.with_suffix(path.suffix + ".gz")
            with gzip.open(out, "wt") as fh:
                fh.write(content)
        else:
            path.write_text(content)

    _put(matrix_name, _MTX_TEXT)
    _put(barcodes_name, "cell1\ncell2\n")
    _put(features_name, "iso1\niso2\niso3\n")


def test_load_h5ad(tmp_path: Path) -> None:
    adata = ad.AnnData(
        X=np.array([[1, 2], [3, 4]], dtype=np.float32),
        obs=pd.DataFrame(index=["cell1", "cell2"]),
        var=pd.DataFrame(index=["iso1", "iso2"]),
    )
    h5 = tmp_path / "test.h5ad"
    adata.write_h5ad(h5)
    loaded = load_counts(h5)
    assert loaded.shape == (2, 2)
    assert list(loaded.obs_names) == ["cell1", "cell2"]
    assert list(loaded.var_names) == ["iso1", "iso2"]


def test_load_10x_mtx_directory_plain(tmp_path: Path) -> None:
    _write_mtx_dir(tmp_path)
    loaded = load_counts(tmp_path)
    # After transpose: 2 cells x 3 features
    assert loaded.shape == (2, 3)
    assert list(loaded.obs_names) == ["cell1", "cell2"]
    assert list(loaded.var_names) == ["iso1", "iso2", "iso3"]
    # Spot-check a value: cell1 / iso1 = 5
    cell1_iso1 = loaded[loaded.obs_names == "cell1", loaded.var_names == "iso1"].X
    if hasattr(cell1_iso1, "toarray"):
        cell1_iso1 = cell1_iso1.toarray()
    assert int(cell1_iso1.flatten()[0]) == 5


def test_load_10x_mtx_directory_gzipped(tmp_path: Path) -> None:
    _write_mtx_dir(
        tmp_path,
        matrix_name="matrix.mtx.gz",
        barcodes_name="barcodes.tsv.gz",
        features_name="features.tsv.gz",
    )
    loaded = load_counts(tmp_path)
    assert loaded.shape == (2, 3)
    assert list(loaded.obs_names) == ["cell1", "cell2"]
    assert list(loaded.var_names) == ["iso1", "iso2", "iso3"]


def test_load_10x_mtx_directory_with_genes_tsv(tmp_path: Path) -> None:
    """Older 10x convention uses genes.tsv instead of features.tsv."""
    _write_mtx_dir(tmp_path, features_name="genes.tsv")
    loaded = load_counts(tmp_path)
    assert loaded.shape == (2, 3)
    assert list(loaded.var_names) == ["iso1", "iso2", "iso3"]


def test_load_10x_mtx_features_with_extra_columns(tmp_path: Path) -> None:
    """features.tsv may have id\tname\ttype columns; only id is needed."""
    (tmp_path / "matrix.mtx").write_text(_MTX_TEXT)
    (tmp_path / "barcodes.tsv").write_text("cell1\ncell2\n")
    (tmp_path / "features.tsv").write_text(
        "iso1\tname1\tIsoform\niso2\tname2\tIsoform\niso3\tname3\tIsoform\n"
    )
    loaded = load_counts(tmp_path)
    assert list(loaded.var_names) == ["iso1", "iso2", "iso3"]


def test_load_directory_missing_matrix_raises(tmp_path: Path) -> None:
    (tmp_path / "barcodes.tsv").write_text("cell1\n")
    (tmp_path / "features.tsv").write_text("iso1\n")
    with pytest.raises(FileNotFoundError):
        load_counts(tmp_path)


def test_load_directory_missing_features_raises(tmp_path: Path) -> None:
    (tmp_path / "matrix.mtx").write_text(_MTX_TEXT)
    (tmp_path / "barcodes.tsv").write_text("cell1\ncell2\n")
    with pytest.raises(FileNotFoundError):
        load_counts(tmp_path)


def test_load_path_neither_h5ad_nor_directory_raises(tmp_path: Path) -> None:
    bogus = tmp_path / "something.txt"
    bogus.write_text("nope")
    with pytest.raises(ValueError):
        load_counts(bogus)
