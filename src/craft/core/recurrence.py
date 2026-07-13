"""Per-isoform recurrence from a per-cell count matrix.

``n_cells_detected`` and ``total_count`` are descriptive support features.  They
are not truth labels and are not depth invariant: both depend on capture,
sequencing saturation, cell number and cell-type abundance.

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
from scipy.stats import betabinom, binom, norm

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


def _whitelist_mask(
    counts: ad.AnnData, cell_whitelist: Iterable[str] | None
) -> np.ndarray | None:
    """Boolean cell mask for the whitelist, or None if absent or matching nothing."""
    if cell_whitelist is None:
        return None
    wl = set(cell_whitelist)
    mask = np.fromiter((b in wl for b in counts.obs_names), bool, counts.n_obs)
    return mask if mask.any() else None


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
        total_barcodes = counts.n_obs
        mask = _whitelist_mask(counts, cell_whitelist)
        if mask is None:
            print(
                f"[craft] cell whitelist matched 0 of {total_barcodes} matrix barcodes; "
                "recurrence computed over all cells.",
                file=sys.stderr,
            )
        else:
            counts = counts[mask]
            print(
                f"[craft] recurrence over {int(mask.sum())} whitelisted cells "
                f"(of {total_barcodes} barcodes in the count matrix).",
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
            "n_cells_total": int(counts.n_obs),
            "detection_fraction": n_cells / max(int(counts.n_obs), 1),
            "molecules_per_detected_cell": np.divide(
                total,
                n_cells,
                out=np.full(total.shape, np.nan, dtype=float),
                where=n_cells > 0,
            ),
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


_CONFIDENCE_COLUMNS = ["transcript_id", "recurrence_pvalue", "recurrence_score"]


def _occupancy_pvalue(counts: ad.AnnData) -> pd.Series:
    """Upper-tail p-value that an isoform occupies >= its observed cell count.

    Null: independent Poisson molecule counts have mean ``T * p_c`` per cell,
    where ``p_c`` is the cell's fraction of total library depth and ``T`` is the
    observed isoform total. Occupancy indicators are then independent with
    probability ``1 - exp(-T * p_c)``. The occupied-cell count is Poisson-binomial;
    its upper tail uses a normal approximation with continuity correction.

    This is a coherent Poissonized approximation to fixed-total multinomial
    allocation. Treat the result as an exploratory dispersion statistic: it does
    not model expression heterogeneity and is not an isoform-validity probability.
    Returns a Series indexed by ``var_names`` (NaN for isoforms with no molecules).
    """
    x = counts.X
    if sp.issparse(x):
        x = x.tocsc()
        depth = np.asarray(x.sum(axis=1)).ravel()
        total = np.asarray(x.sum(axis=0)).ravel()
        occ = np.asarray((x > 0).sum(axis=0)).ravel()
    else:
        x = np.asarray(x)
        depth = x.sum(axis=1).ravel()
        total = x.sum(axis=0).ravel()
        occ = (x > 0).sum(axis=0).ravel()

    var_names = np.asarray(counts.var_names)
    n_cells = int(counts.n_obs)
    depth_total = float(depth.sum())
    if depth_total <= 0 or n_cells == 0:
        return pd.Series(np.nan, index=var_names)

    p = depth / depth_total
    totals = np.rint(total).astype(np.int64)

    moments: dict[int, tuple[float, float]] = {}
    for molecules in np.unique(totals):
        occupied_probability = 1.0 - np.exp(-int(molecules) * p)
        moments[int(molecules)] = (
            float(occupied_probability.sum()),
            float((occupied_probability * (1.0 - occupied_probability)).sum()),
        )

    mean = np.array([moments[int(t)][0] for t in totals], dtype=float)
    variance = np.array([moments[int(t)][1] for t in totals], dtype=float)
    variance = np.clip(variance, 1e-12, None)
    z = (occ - 0.5 - mean) / np.sqrt(variance)
    pval = norm.sf(z)
    pval = np.where(totals <= 0, np.nan, pval)
    return pd.Series(pval, index=var_names)


def _betabinom_pvalue(
    k: np.ndarray, n_cells: int, classes: np.ndarray | None = None
) -> np.ndarray:
    """Upper-tail p-value under a beta-binomial fitted (by moments) to the data.

    The observed (cells-detected, n-cells) counts define an empirical detection
    distribution; a beta-binomial fitted to it gives ``P(K >= k)`` per isoform.
    When ``classes`` is given, a separate fit is made per class. A degenerate class
    (no variance) falls back to a binomial with the class mean rate.
    """
    k = np.asarray(k, dtype=float)
    pval = np.full(k.shape, np.nan)
    if n_cells <= 0:
        return pval
    if classes is None:
        groups = [np.ones(k.shape, dtype=bool)]
    else:
        cl = pd.Series(classes).astype(str).to_numpy()
        groups = [cl == g for g in pd.unique(cl)]

    for g in groups:
        valid = g & ~np.isnan(k)
        if valid.sum() < 2:
            continue
        kk = k[valid]
        f = kk / n_cells
        m = float(np.mean(f))
        v = float(np.var(f))
        if v <= 0 or m <= 0 or m >= 1:
            rate = min(max(m, 1e-9), 1 - 1e-9)
            pval[valid] = binom.sf(kk - 1, n_cells, rate)
            continue
        common = m * (1 - m) / v - 1
        a = max(m * common, 1e-6)
        b = max((1 - m) * common, 1e-6)
        pval[valid] = betabinom.sf(kk - 1, n_cells, a, b)
    return pval


def recurrence_confidence(
    counts: ad.AnnData,
    recurrence: pd.DataFrame,
    method: str = "none",
    cell_whitelist: Iterable[str] | None = None,
    classes: pd.Series | None = None,
) -> pd.DataFrame:
    """Exploratory recurrence: is detection broader than a chosen null predicts?

    Returns ``transcript_id``, ``recurrence_pvalue`` (upper-tail probability of the
    observed cell count under the null) and ``recurrence_score`` (``1 - pvalue``;
    higher = broader detection relative to that null). Neither statistic is a
    probability that an isoform is real. Both are NaN when ``method == "none"``
    (the default, which leaves the prior output unchanged) and for isoforms with no
    molecules.

    Args:
        counts: the per-cell count matrix (whitelist applied here if given, matching
            :func:`compute_recurrence`).
        recurrence: output of :func:`compute_recurrence` (defines isoform order).
        method: ``occupancy`` (depth-aware Poissonized occupancy null), ``betabinom``
            (empirical beta-binomial on the cells-detected counts), or ``none``.
        classes: optional per-isoform structural class, aligned to ``recurrence``
            rows, used to stratify the ``betabinom`` fit.
    """
    tx = recurrence["transcript_id"].to_numpy()
    if method == "none":
        return pd.DataFrame(
            {"transcript_id": tx, "recurrence_pvalue": np.nan, "recurrence_score": np.nan}
        )[_CONFIDENCE_COLUMNS]

    mask = _whitelist_mask(counts, cell_whitelist)
    if mask is not None:
        counts = counts[mask]

    if method == "occupancy":
        by_var = _occupancy_pvalue(counts)
        pval = by_var.reindex(pd.Index(tx)).to_numpy()
    elif method == "betabinom":
        k = recurrence["n_cells_detected"].to_numpy(dtype=float)
        cls = classes.to_numpy() if classes is not None else None
        pval = _betabinom_pvalue(k, int(counts.n_obs), cls)
    else:
        raise ValueError(f"Unknown recurrence null: {method!r}")

    return pd.DataFrame(
        {"transcript_id": tx, "recurrence_pvalue": pval, "recurrence_score": 1.0 - pval}
    )[_CONFIDENCE_COLUMNS]
