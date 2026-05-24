"""In-silico truncation of GENCODE transcripts for Bench 1.

Truncate a fraction of the transcript length from the 5' end, 3' end, or both.
"Transcript orientation" matters: for the ``-`` strand the genomic-rightmost
base is the transcript's 5' end, so the same ``rate`` + ``orientation`` arguments
trim opposite sides in genomic coordinates depending on strand.

The exon truncator and the sequence truncator are kept in lock-step: the same
(rate, orientation) call on a transcript's exons and on its transcript-orientation
sequence produces a coherent pair (the resulting exon spans cover exactly the
returned sequence).
"""

from collections.abc import Sequence

Exon = tuple[int, int]  # half-open [start, end) genomic coordinates


def truncate_exons(
    exons: Sequence[Exon],
    strand: str,
    rate: float,
    orientation: str,
) -> list[Exon]:
    """Trim ``rate`` of the transcript length from the requested end(s).

    Args:
        exons: list of ``(start, end)`` half-open genomic intervals, sorted by
            ``start``. Order is genomic, not transcript.
        strand: ``"+"`` or ``"-"``.
        rate: fraction of transcript length to remove. ``0`` returns input
            unchanged; ``1`` returns ``[]``.
        orientation: ``"5prime"``, ``"3prime"``, or ``"both"``. For ``both``,
            half of ``rate`` is removed from each end (the 5' side gets the
            integer half when ``rate * total_length`` is odd).

    Returns:
        New exon intervals in genomic order. First and/or last exon may be
        shortened; exons fully within the trimmed portion are removed. Returns
        ``[]`` if the trim consumes the whole transcript.
    """
    if rate <= 0:
        return [tuple(e) for e in exons]
    total_length = sum(e - s for s, e in exons)
    if rate >= 1.0 or total_length == 0:
        return []
    if strand not in ("+", "-"):
        raise ValueError(f"strand must be + or -, got {strand!r}")

    trim_bp = int(total_length * rate)
    if orientation == "5prime":
        keep_start_tx, keep_end_tx = trim_bp, total_length
    elif orientation == "3prime":
        keep_start_tx, keep_end_tx = 0, total_length - trim_bp
    elif orientation == "both":
        half = trim_bp // 2
        keep_start_tx = half
        keep_end_tx = total_length - (trim_bp - half)
    else:
        raise ValueError(
            f"orientation must be 5prime/3prime/both, got {orientation!r}"
        )

    if keep_end_tx <= keep_start_tx:
        return []

    return _slice_exons_by_transcript_range(
        exons, strand, keep_start_tx, keep_end_tx
    )


def truncate_sequence(seq: str, rate: float, orientation: str) -> str:
    """Trim ``rate`` of the sequence from the 5' end (start), 3' end (end), or both.

    ``seq`` is expected in transcript orientation (already reverse-complemented
    for ``-`` strand transcripts).
    """
    if rate <= 0:
        return seq
    n = len(seq)
    if rate >= 1.0 or n == 0:
        return ""

    trim = int(n * rate)
    if orientation == "5prime":
        return seq[trim:]
    if orientation == "3prime":
        return seq[: n - trim]
    if orientation == "both":
        half = trim // 2
        return seq[half : n - (trim - half)]
    raise ValueError(
        f"orientation must be 5prime/3prime/both, got {orientation!r}"
    )


def _slice_exons_by_transcript_range(
    exons: Sequence[Exon], strand: str, tx_start: int, tx_end: int
) -> list[Exon]:
    """Return exons covering ``[tx_start, tx_end)`` in transcript coordinates."""
    ordered = list(exons) if strand == "+" else list(reversed(exons))
    out_tx_order: list[Exon] = []
    cursor = 0
    for g_start, g_end in ordered:
        exon_len = g_end - g_start
        exon_tx_start = cursor
        exon_tx_end = cursor + exon_len
        overlap_start = max(exon_tx_start, tx_start)
        overlap_end = min(exon_tx_end, tx_end)
        if overlap_start < overlap_end:
            offset_left = overlap_start - exon_tx_start
            offset_right = overlap_end - exon_tx_start
            if strand == "+":
                new_g_start = g_start + offset_left
                new_g_end = g_start + offset_right
            else:
                new_g_end = g_end - offset_left
                new_g_start = g_end - offset_right
            out_tx_order.append((new_g_start, new_g_end))
        cursor = exon_tx_end
        if cursor >= tx_end:
            break

    return out_tx_order if strand == "+" else list(reversed(out_tx_order))
