"""User-supplied polyA atlas loader and matching.

Accepts a BED file of curated polyadenylation sites (PAS), indexes it by
chromosome+strand, and exposes a per-isoform 3'-end match function. This is the
direct-evidence alternative to the v1.1 canonical poly(A) signal motif scan in
:mod:`craft.core.utr3`.

Recommended sources (user-supplied; CRAFT does not auto-download):

- PolyASite v3.0 (https://polyasite.unibas.ch/) -- multi-species, scRNA-seq-derived,
  ships at three stringency levels.
- PolyA_DB v4 (https://exon.apps.wistar.org/PolyA_DB/v4/) -- human + mouse, from
  3' end and long-read sequencing.

Expected file format: BED 6-column at minimum. Lines starting with ``#``,
``track``, or ``browser`` are skipped. Extra columns are tolerated and ignored.
``.bed`` and ``.bed.gz`` both work.
"""

import gzip
from pathlib import Path

import pandas as pd
import pyranges as pr


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def load_atlas(bed_path: Path) -> pr.PyRanges:
    """Read a BED file of polyA sites into a stranded PyRanges.

    Args:
        bed_path: path to a BED file (plain or ``.bed.gz``). Must have at least
            6 tab-separated columns: chrom, chromStart, chromEnd, name, score,
            strand. Extra columns are ignored. Header lines starting with
            ``#``, ``track``, or ``browser`` are skipped.

    Returns:
        PyRanges with columns ``Chromosome``, ``Start``, ``End``, ``Name``,
        ``Score``, ``Strand``. ``Strand`` is plain ``str`` (not categorical),
        which keeps the PyRanges flagged stranded after downstream slicing.
    """
    rows: list[tuple[str, int, int, str, float, str]] = []
    with _open_text(bed_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3] if parts[3] else ""
            try:
                score = float(parts[4])
            except (ValueError, IndexError):
                score = 0.0
            strand = parts[5]
            if strand not in ("+", "-"):
                continue
            rows.append((chrom, start, end, name, score, strand))

    if not rows:
        return pr.PyRanges(
            pd.DataFrame(
                columns=["Chromosome", "Start", "End", "Name", "Score", "Strand"]
            )
        )

    df = pd.DataFrame(
        rows, columns=["Chromosome", "Start", "End", "Name", "Score", "Strand"]
    )
    df["Strand"] = df["Strand"].astype(str)
    return pr.PyRanges(df)


def match_iso_end(
    iso_3prime_pos: int,
    chrom: str,
    strand: str,
    atlas: pr.PyRanges,
    tolerance: int = 24,
) -> dict[str, int | str | bool]:
    """Look up an iso's 3' end against the atlas.

    Args:
        iso_3prime_pos: genomic 0-based position of the iso's 3' end. For ``+``
            strand: the last exon's End - 1 (last exonic base). For ``-`` strand:
            the leftmost exon's Start (first exonic base in genome coordinates).
        chrom: chromosome name. Must match atlas chromosome naming.
        strand: ``+`` or ``-``.
        atlas: PyRanges from :func:`load_atlas`.
        tolerance: nt window around the iso 3' end. The PAS midpoint
            ``(Start + End) // 2`` must fall within ``[iso_3prime_pos - tol,
            iso_3prime_pos + tol]`` for a hit.

    Returns:
        Dict with keys ``matched`` (bool), ``pas_id`` (str, atlas Name column;
        empty when no match), ``distance_nt`` (int, signed nt from iso 3' end to
        matched PAS midpoint; -1 when no match).
    """
    empty: dict[str, int | str | bool] = {
        "matched": False,
        "pas_id": "",
        "distance_nt": -1,
    }
    if len(atlas) == 0:
        return empty

    df = atlas.df
    candidates = df[(df["Chromosome"] == chrom) & (df["Strand"] == strand)]
    if candidates.empty:
        return empty

    midpoints = ((candidates["Start"] + candidates["End"]) // 2).to_numpy()
    distances = midpoints - iso_3prime_pos
    abs_distances = abs(distances)
    if abs_distances.size == 0:
        return empty

    best_idx = int(abs_distances.argmin())
    best_dist = int(abs_distances[best_idx])
    if best_dist > tolerance:
        return empty

    matched_row = candidates.iloc[best_idx]
    return {
        "matched": True,
        "pas_id": str(matched_row["Name"]) if pd.notna(matched_row["Name"]) else "",
        "distance_nt": int(distances[best_idx]),
    }
