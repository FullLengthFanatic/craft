"""AnnData and MuData writers."""

from pathlib import Path

import anndata as ad
import pandas as pd


def to_anndata(per_isoform: pd.DataFrame, counts: ad.AnnData | None = None) -> ad.AnnData:
    """Build an AnnData with isoforms as `var` and functional annotations as `var` columns."""
    raise NotImplementedError


def write_h5ad(adata: ad.AnnData, output: Path) -> None:
    """Serialize the AnnData to disk."""
    raise NotImplementedError
