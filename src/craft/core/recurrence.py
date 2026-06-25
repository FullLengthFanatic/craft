"""Per-isoform recurrence from a per-cell count matrix.

``n_cells_detected`` and ``total_count`` are depth-stable filtering signals: an
isoform observed in many independent cells is supported regardless of per-cell
sequencing depth, unlike a raw read count, which scales with how deeply each cell was sequenced.
``isoform_fraction_within_gene`` puts an isoform's abundance on a per-gene scale,
which cancels the depth term as well.

Restrict to called cells with an optional barcode whitelist. The count matrices
that ship with pigeon / isoseq carry every observed barcode, most of which are
ambient droplets, so counting cells over the raw matrix inflates recurrence.
"""

import sys
from collections.abc import Iterable
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

_HEADER_TOKENS = {"barcode", "barcodes", "cell", "cell_barcode", "cellbarcode", "cb"}


def load_cell_whitelist(path: Path) -> list[str]:
    """Read a called-cell barcode list (one barcode per line, first column).

    A leading header line (e.g. ``barcode`` / ``cell_barcode``) is skipped. Blank
    lines are ignored. Only the first tab/whitespace-separated token is taken, so a
    plain list or a wider TSV both work.
    """
    barcodes: list[str] = []
    with open(path) as fh:
        for i, line in enumerate(fh):
            tok = line.strip().split()[0] if line.strip() else ""
            if not tok:
                continue
            if i == 0 and tok.lower() in _HEADER_TOKENS:
                continue
            barcodes.append(tok)
    return barcodes


def compute_recurrence(
    counts: ad.AnnData, cell_whitelist: Iterable[str] | None = None
) -> pd.DataFrame:
    """Per-isoform molecule total and cell recurrence over (optionally whitelisted) cells.

    Args:
        counts: AnnData with cells in ``obs`` and isoforms in ``var`` (``var_names``
            are transcript ids matching the per-isoform table).
        cell_whitelist: Optional called-cell barcodes. ``obs`` is subset to the
            intersection; if nothing matches, all cells are used and a note is logged.

    Returns:
        DataFrame with ``transcript_id``, ``n_cells_detected``, ``total_count``.
    """
    if cell_whitelist is not None:
        wl = set(cell_whitelist)
        mask = np.fromiter((b in wl for b in counts.obs_names), bool, counts.n_obs)
        n = int(mask.sum())
        if n == 0:
            print(
                f"[craft] cell whitelist matched 0 of {counts.n_obs} matrix barcodes; "
                "recurrence computed over all cells.",
                file=sys.stderr,
            )
        else:
            counts = counts[mask]
            print(
                f"[craft] recurrence over {n} whitelisted cells "
                f"(of {len(mask)} barcodes in the count matrix).",
                file=sys.stderr,
            )

    x = counts.X
    if sp.issparse(x):
        x = x.tocsc()
        total = np.asarray(x.sum(axis=0)).ravel()
        n_cells = np.asarray((x > 0).sum(axis=0)).ravel()
    else:
        x = np.asarray(x)
        total = x.sum(axis=0).ravel()
        n_cells = (x > 0).sum(axis=0).ravel()

    return pd.DataFrame(
        {
            "transcript_id": np.asarray(counts.var_names),
            "n_cells_detected": np.rint(n_cells).astype("int64"),
            "total_count": np.rint(total).astype("int64"),
        }
    )


def within_gene_fraction(
    total_count: pd.Series, parent_gene_id: pd.Series
) -> pd.Series:
    """Each isoform's share of its parent gene's molecules.

    Isoforms with no parent gene (``""`` / NaN) or whose gene has zero total return
    NaN, as do isoforms with no measured count (NaN ``total_count``).
    """
    gene = parent_gene_id.where(parent_gene_id.astype(str) != "", other=np.nan)
    gene_total = total_count.groupby(gene).transform("sum")
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = total_count / gene_total
    return frac.where(gene.notna() & (gene_total > 0), np.nan)
