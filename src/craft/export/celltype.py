"""Per-cell-type functional-consequence aggregation.

Given per-cell isoform counts (an AnnData with cells in ``obs`` and isoforms in
``var``) and CRAFT's per-isoform consequence calls, this collapses the calls to
expression-weighted fractions per cell group. For a group ``g`` and consequence
class ``c`` the fraction is

    sum of molecules in g over isoforms with class c
    -------------------------------------------------
    sum of all molecules in g

so a highly expressed isoform contributes proportionally to its read support.

This is the layer the single-cell long-read callers (FLAMES, scisorseqr,
Isosceles) do not provide: they stop at abundance; this turns abundance plus
CRAFT's annotations into "which cell types express more NMD-targeted /
truncated / domain-disrupted isoforms".
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

# Consequence class -> (required column, predicate over that column's Series).
# Order here is the output column order.
_CONSEQUENCES: dict[str, tuple[str, callable]] = {
    "nmd_sensitive": ("nmd_status", lambda s: s == "sensitive"),
    "ptc_introduced": ("ptc_introduced", lambda s: s.fillna(False).astype(bool)),
    "intron_retained_in_cds": ("intron_retained_in_cds", lambda s: s.fillna(False).astype(bool)),
    "truncated_5p": ("completeness", lambda s: s == "truncated_5p"),
    "truncated_3p": ("completeness", lambda s: s == "truncated_3p"),
    "truncated_both": ("completeness", lambda s: s == "truncated_both"),
    "internal_fragment": ("completeness", lambda s: s == "internal_fragment"),
    "alt_3prime_end": ("completeness", lambda s: s == "alt_3prime_end"),
    "domain_lost": ("pfam_lost", lambda s: s.apply(lambda v: isinstance(v, list) and len(v) > 0)),
}


def _build_consequence_flags(per_isoform: pd.DataFrame) -> pd.DataFrame:
    """Boolean DataFrame (transcript_id index, one column per consequence class)."""
    flags = pd.DataFrame(index=per_isoform["transcript_id"].astype(str))
    for name, (col, predicate) in _CONSEQUENCES.items():
        if col in per_isoform.columns:
            flags[name] = predicate(per_isoform[col]).fillna(False).to_numpy(dtype=bool)
        else:
            flags[name] = False
    return flags


def aggregate_consequences(
    adata: ad.AnnData,
    per_isoform: pd.DataFrame,
    group_by: str,
    output_tsv: Path | None = None,
) -> pd.DataFrame:
    """Expression-weighted per-group consequence fractions.

    Args:
        adata: AnnData with per-cell counts (cells in ``obs``, isoforms in
            ``var`` indexed by ``transcript_id``) and ``X`` the count matrix.
        per_isoform: CRAFT per-isoform table (must include ``transcript_id`` and
            the consequence columns).
        group_by: A column in ``adata.obs`` to group cells by.
        output_tsv: Optional path; the result is written there as TSV.

    Returns:
        DataFrame with one row per group and columns ``cell_group``, ``n_cells``,
        ``total_molecules``, ``n_isoforms`` and ``frac_<consequence>``.

    Raises:
        ValueError: if ``group_by`` is not a column of ``adata.obs`` or there is
            no isoform overlap between ``adata.var`` and ``per_isoform``.
    """
    if group_by not in adata.obs.columns:
        raise ValueError(
            f"--group-by column {group_by!r} not in counts obs; "
            f"available: {list(adata.obs.columns)}"
        )

    flags = _build_consequence_flags(per_isoform)
    var_names = adata.var_names.astype(str)
    shared = var_names.intersection(flags.index)
    if len(shared) == 0:
        raise ValueError(
            "No overlap between counts var_names and per-isoform transcript_id; "
            "cannot aggregate."
        )

    # Align flags to the count matrix column order; isoforms absent from the
    # per-isoform table contribute no consequence but still count in the total.
    flags = flags.reindex(var_names, fill_value=False)
    flag_arrays = {name: flags[name].to_numpy(dtype=bool) for name in _CONSEQUENCES}

    labels = adata.obs[group_by]
    rows: list[dict] = []
    for label in pd.unique(labels.dropna()):
        mask = (labels == label).to_numpy()
        xg = adata.X[mask, :]
        col_sums = np.asarray(xg.sum(axis=0)).ravel() if sp.issparse(xg) else xg.sum(axis=0)
        col_sums = np.asarray(col_sums).ravel()
        total = float(col_sums.sum())
        row = {
            "cell_group": str(label),
            "n_cells": int(mask.sum()),
            "total_molecules": total,
            "n_isoforms": int((col_sums > 0).sum()),
        }
        for name, arr in flag_arrays.items():
            row[f"frac_{name}"] = float(col_sums[arr].sum()) / total if total > 0 else float("nan")
        rows.append(row)

    result = pd.DataFrame(rows)
    if output_tsv is not None:
        result.to_csv(Path(output_tsv), sep="\t", index=False, float_format="%.4f")
    adata.uns["celltype_consequences"] = {"group_by": group_by, "aggregates": result}
    return result
