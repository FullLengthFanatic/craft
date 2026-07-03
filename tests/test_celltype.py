"""Tests for craft.export.celltype."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from craft.export.celltype import aggregate_consequences, celltype_as_nmd


def _per_isoform() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "transcript_id": ["iso1", "iso2", "iso3"],
            "completeness": ["full_length", "truncated_5p", "alt_3prime_end"],
            "nmd_status": ["sensitive", "escaped", "not_applicable"],
            "ptc_introduced": [True, False, False],
            "intron_retained_in_cds": [False, False, False],
            "pfam_lost": [[], ["PF00001"], []],
        }
    )


def _adata(x: np.ndarray | sp.spmatrix) -> ad.AnnData:
    obs = pd.DataFrame({"cell_type": ["A", "A", "B", "B"]}, index=["c1", "c2", "c3", "c4"])
    var = pd.DataFrame(index=["iso1", "iso2", "iso3"])
    return ad.AnnData(X=x, obs=obs, var=var)


def _grp(result: pd.DataFrame, label: str) -> pd.Series:
    return result[result["cell_group"] == label].iloc[0]


_COUNTS = np.array(
    [
        [10.0, 0.0, 5.0],
        [10.0, 0.0, 5.0],
        [0.0, 8.0, 2.0],
        [0.0, 8.0, 2.0],
    ]
)


def test_happy_path_molecule_weighted_fractions(tmp_path) -> None:
    adata = _adata(_COUNTS)
    out = aggregate_consequences(adata, _per_isoform(), "cell_type", tmp_path / "ct.tsv")
    a = _grp(out, "A")
    assert a["n_cells"] == 2
    assert a["total_molecules"] == 30
    assert a["n_isoforms"] == 2
    assert a["frac_nmd_sensitive"] == pytest.approx(20 / 30)
    assert a["frac_ptc_introduced"] == pytest.approx(20 / 30)
    assert a["frac_alt_3prime_end"] == pytest.approx(10 / 30)
    assert a["frac_truncated_5p"] == pytest.approx(0.0)
    b = _grp(out, "B")
    assert b["frac_truncated_5p"] == pytest.approx(16 / 20)
    assert b["frac_domain_lost"] == pytest.approx(16 / 20)
    assert b["frac_nmd_sensitive"] == pytest.approx(0.0)
    assert (tmp_path / "ct.tsv").exists()
    assert adata.uns["celltype_consequences"]["group_by"] == "cell_type"


def test_sparse_matches_dense() -> None:
    dense = aggregate_consequences(_adata(_COUNTS), _per_isoform(), "cell_type")
    sparse = aggregate_consequences(_adata(sp.csr_matrix(_COUNTS)), _per_isoform(), "cell_type")
    pd.testing.assert_frame_equal(
        dense.set_index("cell_group").sort_index(),
        sparse.set_index("cell_group").sort_index(),
    )


def test_missing_group_column_raises() -> None:
    with pytest.raises(ValueError, match="not in counts obs"):
        aggregate_consequences(_adata(_COUNTS), _per_isoform(), "leiden")


def test_no_isoform_overlap_raises() -> None:
    per = _per_isoform()
    per["transcript_id"] = ["x1", "x2", "x3"]
    with pytest.raises(ValueError, match="No overlap"):
        aggregate_consequences(_adata(_COUNTS), per, "cell_type")


def test_zero_count_group_yields_nan_fractions() -> None:
    counts = _COUNTS.copy()
    counts[2:, :] = 0.0  # group B has no molecules
    out = aggregate_consequences(_adata(counts), _per_isoform(), "cell_type")
    b = _grp(out, "B")
    assert b["total_molecules"] == 0
    assert b["n_cells"] == 2
    assert np.isnan(b["frac_nmd_sensitive"])


def test_isoform_in_counts_not_in_per_isoform_counts_in_total() -> None:
    # iso3 dropped from per_isoform: still counted in totals, never in numerators.
    per = _per_isoform().iloc[:2].copy()
    out = aggregate_consequences(_adata(_COUNTS), per, "cell_type")
    a = _grp(out, "A")
    assert a["total_molecules"] == 30  # iso3 molecules still in denominator
    assert a["frac_alt_3prime_end"] == pytest.approx(0.0)  # iso3 not annotated


def test_as_nmd_lists_recurrent_sensitive_isoforms_by_group(tmp_path) -> None:
    # iso1 is the only NMD-sensitive isoform; it lives in group A (20 mol) not B.
    out = celltype_as_nmd(_adata(_COUNTS), _per_isoform(), "cell_type", tmp_path / "as_nmd.tsv")
    assert set(out["cell_group"]) == {"A"}
    row = out.iloc[0]
    assert row["transcript_id"] == "iso1"
    assert row["molecules_in_group"] == 20
    assert row["frac_of_group"] == pytest.approx(20 / 30)
    assert (tmp_path / "as_nmd.tsv").exists()


def test_as_nmd_recurrence_score_gates_membership() -> None:
    per = _per_isoform()
    per["recurrence_score"] = [0.99, 0.99, 0.99]  # iso1 sensitive + recurrent -> kept
    kept = celltype_as_nmd(_adata(_COUNTS), per, "cell_type", recurrence_score_min=0.95)
    assert "iso1" in set(kept["transcript_id"])

    per["recurrence_score"] = [0.5, 0.99, 0.99]  # iso1 now below threshold -> dropped
    dropped = celltype_as_nmd(_adata(_COUNTS), per, "cell_type", recurrence_score_min=0.95)
    assert dropped.empty


def test_as_nmd_missing_group_raises() -> None:
    with pytest.raises(ValueError, match="not in counts obs"):
        celltype_as_nmd(_adata(_COUNTS), _per_isoform(), "leiden")


def test_as_nmd_without_nmd_status_is_empty() -> None:
    per = _per_isoform().drop(columns=["nmd_status"])
    out = celltype_as_nmd(_adata(_COUNTS), per, "cell_type")
    assert out.empty
    assert list(out.columns)[0] == "cell_group"
