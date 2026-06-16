"""Unit tests for craft.core.recurrence."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from craft.core.recurrence import (
    compute_recurrence,
    load_cell_whitelist,
    within_gene_fraction,
)


def _adata(matrix, isoforms, cells):
    a = ad.AnnData(X=matrix)
    a.obs_names = pd.Index(cells)
    a.var_names = pd.Index(isoforms)
    return a


# cells x isoforms: A in all 3 cells (sum 6), B in 1 cell (sum 5), C in none.
_X = np.array([[2, 0, 0], [1, 5, 0], [3, 0, 0]], dtype=float)


def test_compute_recurrence_dense():
    df = compute_recurrence(_adata(_X, ["A", "B", "C"], ["c1", "c2", "c3"])).set_index(
        "transcript_id"
    )
    assert df.loc["A", "total_count"] == 6
    assert df.loc["A", "n_cells_detected"] == 3
    assert df.loc["B", "total_count"] == 5
    assert df.loc["B", "n_cells_detected"] == 1
    assert df.loc["C", "total_count"] == 0
    assert df.loc["C", "n_cells_detected"] == 0
    assert df["total_count"].dtype == "int64"


def test_compute_recurrence_sparse_matches_dense():
    iso, cells = ["A", "B", "C"], ["c1", "c2", "c3"]
    dense = compute_recurrence(_adata(_X, iso, cells))
    sparse = compute_recurrence(_adata(sp.csr_matrix(_X), iso, cells))
    pd.testing.assert_frame_equal(dense, sparse)


def test_compute_recurrence_whitelist_subsets_cells():
    df = compute_recurrence(
        _adata(_X, ["A", "B", "C"], ["c1", "c2", "c3"]), cell_whitelist=["c1", "c3"]
    ).set_index("transcript_id")
    assert df.loc["A", "total_count"] == 5  # 2 + 3, c2 excluded
    assert df.loc["A", "n_cells_detected"] == 2
    assert df.loc["B", "total_count"] == 0  # B only lived in c2
    assert df.loc["B", "n_cells_detected"] == 0


def test_compute_recurrence_empty_whitelist_falls_back_to_all(capsys):
    df = compute_recurrence(
        _adata(_X[:2, :2], ["A", "B"], ["c1", "c2"]), cell_whitelist=["not-a-barcode"]
    ).set_index("transcript_id")
    assert df.loc["A", "total_count"] == 3  # all cells used
    assert "matched 0" in capsys.readouterr().err


def test_within_gene_fraction_basic():
    total = pd.Series([10.0, 30.0, 5.0])
    gene = pd.Series(["g1", "g1", "g2"])
    frac = within_gene_fraction(total, gene)
    assert frac[0] == pytest.approx(0.25)
    assert frac[1] == pytest.approx(0.75)
    assert frac[2] == pytest.approx(1.0)


def test_within_gene_fraction_orphan_and_nan():
    total = pd.Series([10.0, np.nan, 5.0])
    gene = pd.Series(["", "g2", "g2"])
    frac = within_gene_fraction(total, gene)
    assert np.isnan(frac[0])  # orphan: empty gene id
    assert np.isnan(frac[1])  # unmeasured isoform
    assert frac[2] == pytest.approx(1.0)  # only measured isoform of g2


def test_load_cell_whitelist(tmp_path):
    p = tmp_path / "wl.txt"
    p.write_text("cell_barcode\nAAA-BBB\nCCC-DDD\n\nEEE-FFF\textra\n")
    assert load_cell_whitelist(p) == ["AAA-BBB", "CCC-DDD", "EEE-FFF"]
