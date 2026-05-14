"""Tests for craft.export.anndata."""

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from craft.export.anndata import to_anndata, write_h5ad


def _per_isoform(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows)


def test_to_anndata_without_counts_has_empty_obs() -> None:
    df = _per_isoform(
        [
            {"transcript_id": "t1", "orf_outcome": "propagated_intact"},
            {"transcript_id": "t2", "orf_outcome": "no_parent"},
        ]
    )
    adata = to_anndata(df)
    assert adata.n_obs == 0
    assert adata.n_vars == 2
    assert list(adata.var_names) == ["t1", "t2"]
    assert "orf_outcome" in adata.var.columns


def test_to_anndata_with_matching_counts_preserves_X() -> None:
    df = _per_isoform([{"transcript_id": "t1", "orf_outcome": "x"}])
    counts = ad.AnnData(
        X=np.array([[5.0], [3.0]], dtype=np.float32),
        obs=pd.DataFrame(index=["cell1", "cell2"]),
        var=pd.DataFrame(index=["t1"]),
    )
    adata = to_anndata(df, counts=counts)
    assert adata.shape == (2, 1)
    assert list(adata.obs_names) == ["cell1", "cell2"]
    assert adata.X.tolist() == [[5.0], [3.0]]


def test_to_anndata_unmeasured_isoforms_get_zero_columns() -> None:
    df = _per_isoform(
        [
            {"transcript_id": "t1", "orf_outcome": "intact"},
            {"transcript_id": "t2", "orf_outcome": "novel"},  # not in counts
        ]
    )
    counts = ad.AnnData(
        X=np.array([[5.0], [3.0]], dtype=np.float32),
        obs=pd.DataFrame(index=["cell1", "cell2"]),
        var=pd.DataFrame(index=["t1"]),
    )
    adata = to_anndata(df, counts=counts)
    assert adata.shape == (2, 2)
    assert list(adata.var_names) == ["t1", "t2"]
    np.testing.assert_array_equal(np.asarray(adata.X[:, 1]).flatten(), [0.0, 0.0])
    np.testing.assert_array_equal(np.asarray(adata.X[:, 0]).flatten(), [5.0, 3.0])


def test_to_anndata_supports_sparse_counts() -> None:
    df = _per_isoform(
        [
            {"transcript_id": "t1"},
            {"transcript_id": "t2"},
        ]
    )
    counts = ad.AnnData(
        X=sp.csr_matrix(np.array([[5.0, 0.0], [0.0, 3.0]], dtype=np.float32)),
        obs=pd.DataFrame(index=["cell1", "cell2"]),
        var=pd.DataFrame(index=["t1", "t2"]),
    )
    adata = to_anndata(df, counts=counts)
    assert sp.issparse(adata.X)
    assert adata.shape == (2, 2)
    np.testing.assert_array_equal(adata.X.toarray(), [[5.0, 0.0], [0.0, 3.0]])


def test_to_anndata_serializes_list_columns_to_json() -> None:
    df = _per_isoform(
        [
            {
                "transcript_id": "t1",
                "propagated_cds_intervals": [("chr1", 100, 200, "+")],
                "denovo_cds_intervals": [],
            },
        ]
    )
    adata = to_anndata(df)
    raw = adata.var.iloc[0]["propagated_cds_intervals"]
    assert isinstance(raw, str)
    assert json.loads(raw) == [["chr1", 100, 200, "+"]]
    assert json.loads(adata.var.iloc[0]["denovo_cds_intervals"]) == []


def test_write_h5ad_roundtrip(tmp_path: Path) -> None:
    df = _per_isoform([{"transcript_id": "t1", "orf_outcome": "intact"}])
    counts = ad.AnnData(
        X=np.array([[5.0]], dtype=np.float32),
        obs=pd.DataFrame(index=["cell1"]),
        var=pd.DataFrame(index=["t1"]),
    )
    adata = to_anndata(df, counts=counts)
    out = tmp_path / "deep" / "nested" / "test.h5ad"
    write_h5ad(adata, out)
    assert out.exists()

    loaded = ad.read_h5ad(out)
    assert loaded.shape == adata.shape
    assert list(loaded.var_names) == ["t1"]
    assert loaded.var["orf_outcome"].iloc[0] == "intact"
    assert loaded.X.tolist() == [[5.0]]


def test_to_anndata_empty_per_isoform_yields_zero_var() -> None:
    df = pd.DataFrame(columns=["transcript_id", "orf_outcome"])
    adata = to_anndata(df)
    assert adata.n_obs == 0
    assert adata.n_vars == 0
