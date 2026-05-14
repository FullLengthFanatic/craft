"""AnnData and MuData writers.

Builds an AnnData where isoforms live in ``var`` (indexed by ``transcript_id``)
and all per-isoform functional annotations from ``run_annotate`` are var columns.
When per-cell counts are supplied, they are reindexed to the per-isoform table
so unmeasured isoforms get zero-filled count columns rather than being dropped.
"""

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def _serialize_value(v: object) -> str:
    if isinstance(v, list | tuple):
        return json.dumps(v)
    if pd.isna(v):
        return "[]"
    return str(v)


def _serialize_list_columns(df: pd.DataFrame) -> pd.DataFrame:
    """JSON-encode columns containing list/tuple values (h5ad-incompatible)."""
    df = df.copy()
    for col in df.columns:
        is_listlike = df[col].apply(lambda v: isinstance(v, list | tuple))
        if is_listlike.any():
            df[col] = df[col].apply(_serialize_value)
    return df


def _normalize_var_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    """Coerce object-dtype columns to a single homogeneous dtype.

    h5ad cannot serialise object arrays with mixed numeric/None values. For each
    remaining object column we attempt numeric conversion of the non-null
    values; if it succeeds for all of them, the whole column becomes float
    (with NaN replacing None). Otherwise we cast to string with empty-string
    standing in for null.
    """
    df = df.copy()
    for col in df.columns:
        if df[col].dtype != object:
            continue
        non_null = df[col].dropna()
        if non_null.empty:
            df[col] = df[col].fillna("").astype(str)
            continue
        numeric = pd.to_numeric(non_null, errors="coerce")
        if numeric.notna().all():
            df[col] = pd.to_numeric(df[col], errors="coerce")
        else:
            df[col] = df[col].fillna("").astype(str)
    return df


def _reindex_counts(
    counts: ad.AnnData, target: pd.Index
) -> np.ndarray | sp.spmatrix:
    """Build (n_obs, len(target)) X with counts' values for shared isoforms, zeros elsewhere."""
    src_to_idx = {tx: i for i, tx in enumerate(counts.var_names)}
    n_obs = counts.n_obs
    n_target = len(target)

    if sp.issparse(counts.X):
        out = sp.lil_matrix((n_obs, n_target), dtype=counts.X.dtype)
        for j, tx in enumerate(target):
            src = src_to_idx.get(tx)
            if src is not None:
                out[:, j] = counts.X[:, src]
        return out.tocsr()

    dtype = counts.X.dtype if counts.X is not None else np.float32
    out = np.zeros((n_obs, n_target), dtype=dtype)
    for j, tx in enumerate(target):
        src = src_to_idx.get(tx)
        if src is not None:
            col = counts.X[:, src]
            out[:, j] = np.asarray(col).flatten()
    return out


def to_anndata(
    per_isoform: pd.DataFrame,
    counts: ad.AnnData | None = None,
) -> ad.AnnData:
    """Build an AnnData of isoform annotations, optionally with per-cell counts.

    Args:
        per_isoform: Per-isoform annotation table from
            :func:`craft.pipeline.run_annotate`.
        counts: Optional per-cell counts AnnData with cells in obs and
            isoforms in var. Reindexed to the per-isoform table: unmeasured
            isoforms become zero-filled count columns rather than being dropped.

    Returns:
        AnnData with:
        - var: per_isoform, indexed by transcript_id (list-valued columns
          like ``propagated_cds_intervals`` are JSON-encoded so h5ad write works).
        - obs: counts.obs when counts provided; empty (0 rows) otherwise.
        - X: reindexed counts when provided; shape (0, n_var) otherwise.
    """
    var = per_isoform.set_index("transcript_id").copy()
    var = _serialize_list_columns(var)
    var = _normalize_var_dtypes(var)

    if counts is None:
        x = np.zeros((0, len(var)), dtype=np.float32)
        return ad.AnnData(X=x, var=var)

    x = _reindex_counts(counts, var.index)
    return ad.AnnData(X=x, obs=counts.obs.copy(), var=var)


def write_h5ad(adata: ad.AnnData, output: Path) -> None:
    """Serialize the AnnData to disk (creates parent dirs if needed)."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output)
