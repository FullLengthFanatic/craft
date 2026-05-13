"""Per-cell count matrix loaders (h5ad and MTX)."""

from pathlib import Path

import anndata as ad


def load_counts(path: Path) -> ad.AnnData:
    """Load per-cell isoform counts from h5ad or MTX into an AnnData."""
    raise NotImplementedError
