"""ORF prediction baselines for Bench 1.

CRAFT's propagation lives in the main ``craft`` package; this module wraps the
de novo comparator (``orfipy``) plus the GENCODE-truth ORF accessor.
"""

import json
import shutil
import subprocess
import tempfile
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from cbench.data import Transcript


@dataclass(frozen=True)
class ORFCall:
    """A single ORF prediction with transcript coordinates."""

    transcript_id: str
    found: bool
    tx_start: int = -1  # 0-based, inclusive
    tx_end: int = -1  # exclusive
    start_codon: str = ""
    stop_codon: str = ""

    @property
    def length(self) -> int:
        return max(0, self.tx_end - self.tx_start)


def truth_orf(tx: Transcript) -> ORFCall:
    """The GENCODE-truth ORF for ``tx`` (in transcript coordinates)."""
    seq = tx.transcript_seq
    return ORFCall(
        transcript_id=tx.transcript_id,
        found=True,
        tx_start=tx.cds_tx_start,
        tx_end=tx.cds_tx_end,
        start_codon=seq[tx.cds_tx_start : tx.cds_tx_start + 3],
        stop_codon=seq[max(0, tx.cds_tx_end - 3) : tx.cds_tx_end],
    )


def orfipy_predict(
    seqs: dict[str, str],
    min_orf_nt: int = 75,
    starts: tuple[str, ...] = ("ATG",),
    stops: tuple[str, ...] = ("TAA", "TAG", "TGA"),
    extra_args: tuple[str, ...] = (),
    orfipy_bin: str | None = None,
) -> dict[str, ORFCall]:
    """Run orfipy on a dict of ``{transcript_id: transcript_sequence}``.

    Returns the longest ATG-initiated ORF on the forward strand per sequence.
    Sequences with no qualifying ORF get an ``ORFCall(found=False)``.

    Args:
        seqs: transcript-orientation sequences keyed by transcript id.
        min_orf_nt: minimum ORF length to report (orfipy ``--min``).
        starts: allowed start codons. Defaults to ``("ATG",)``.
        stops: allowed stop codons.
        extra_args: extra flags passed straight to ``orfipy`` (e.g. ``("--partial-3",)``).
    """
    if not seqs:
        return {}

    with tempfile.TemporaryDirectory(prefix="cbench_orfipy_") as td:
        td_path = Path(td)
        fasta = td_path / "in.fa"
        with fasta.open("w") as fh:
            for tx_id, seq in seqs.items():
                fh.write(f">{tx_id}\n{seq}\n")
        out_dir = td_path / "out"
        out_dir.mkdir()
        binpath = orfipy_bin or shutil.which("orfipy") or _orfipy_in_venv()
        if binpath is None:
            raise FileNotFoundError(
                "orfipy not found on PATH; pass orfipy_bin or activate the venv"
            )
        cmd = [
            binpath,
            str(fasta),
            "--min",
            str(min_orf_nt),
            "--start",
            ",".join(starts),
            "--stop",
            ",".join(stops),
            "--strand",
            "f",
            "--longest",
            "--bed",
            "orf.bed",
            "--outdir",
            str(out_dir),
            *extra_args,
        ]
        subprocess.run(cmd, check=True, capture_output=True)
        bed = out_dir / "orf.bed"
        return _parse_orfipy_bed(bed, seqs)


def _orfipy_in_venv() -> str | None:
    """Best-effort: locate ``orfipy`` next to the running Python interpreter."""
    import sys

    candidate = Path(sys.executable).parent / "orfipy"
    return str(candidate) if candidate.exists() else None


def _parse_orfipy_bed(
    bed_path: Path, seqs: dict[str, str]
) -> dict[str, ORFCall]:
    out: dict[str, ORFCall] = {tx_id: ORFCall(tx_id, False) for tx_id in seqs}
    if not bed_path.exists() or bed_path.stat().st_size == 0:
        return out
    df = pd.read_csv(
        bed_path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "info", "score", "strand"],
        dtype={"chrom": "string"},
    )
    df = df[df["strand"] == "+"]
    df["length"] = df["end"] - df["start"]
    df = df.sort_values("length", ascending=False).drop_duplicates(
        "chrom", keep="first"
    )
    for row in df.itertuples(index=False):
        seq = seqs.get(row.chrom, "")
        start = int(row.start)
        end = int(row.end)
        out[row.chrom] = ORFCall(
            transcript_id=row.chrom,
            found=True,
            tx_start=start,
            tx_end=end,
            start_codon=seq[start : start + 3] if seq else "",
            stop_codon=seq[max(0, end - 3) : end] if seq else "",
        )
    return out


def craft_propagated_calls(
    per_isoform_tsv: Path,
    transcripts: Iterable[Transcript],
) -> dict[str, ORFCall]:
    """Read a CRAFT ``per_isoform.tsv`` and convert each row to an ``ORFCall``.

    CRAFT's propagated CDS coordinates live in ``propagated_cds_intervals`` as a
    JSON-like list of genomic intervals. We map them back to transcript
    coordinates using the *truncated* iso's exon set (passed in ``transcripts``).
    Rows whose ``orf_outcome`` is one of the no-call categories (``no_parent``,
    ``no_parent_cds``, ``start_lost``) become ``found=False``.
    """
    from cbench.data import _genomic_to_tx_coord  # local to avoid circular hints

    tx_by_id = {t.transcript_id: t for t in transcripts}
    df = pd.read_csv(per_isoform_tsv, sep="\t", dtype={"transcript_id": "string"})
    no_call_outcomes = {"no_parent", "no_parent_cds", "start_lost"}
    out: dict[str, ORFCall] = {}
    for row in df.itertuples(index=False):
        tx_id = str(row.transcript_id)
        outcome = str(getattr(row, "orf_outcome", ""))
        intervals_raw = getattr(row, "propagated_cds_intervals", "")
        intervals = _parse_intervals(intervals_raw)
        if outcome in no_call_outcomes or not intervals:
            out[tx_id] = ORFCall(tx_id, False)
            continue
        tx = tx_by_id.get(tx_id)
        if tx is None or not tx.exons:
            out[tx_id] = ORFCall(tx_id, False)
            continue
        try:
            tx_positions = []
            for g_start, g_end in intervals:
                for pos in (g_start, g_end - 1):
                    tx_positions.append(_genomic_to_tx_coord(pos, tx.exons, tx.strand))
        except ValueError:
            out[tx_id] = ORFCall(tx_id, False)
            continue
        if not tx_positions:
            out[tx_id] = ORFCall(tx_id, False)
            continue
        tx_lo = min(tx_positions)
        tx_hi = max(tx_positions) + 1
        seq = tx.transcript_seq
        out[tx_id] = ORFCall(
            transcript_id=tx_id,
            found=True,
            tx_start=tx_lo,
            tx_end=tx_hi,
            start_codon=seq[tx_lo : tx_lo + 3] if 0 <= tx_lo < len(seq) else "",
            stop_codon=seq[max(0, tx_hi - 3) : tx_hi] if tx_hi <= len(seq) else "",
        )
    return out


def _parse_intervals(raw) -> list[tuple[int, int]]:
    """Parse CRAFT's ``propagated_cds_intervals`` field.

    Format is a JSON-encoded list of ``[chrom, start, end, strand]`` quads, one
    per CDS exon segment. Empty CDS is represented as ``"[]"`` (or NaN when
    CRAFT didn't propagate anything).
    """
    if raw is None or pd.isna(raw):
        return []
    s = str(raw).strip()
    if not s or s == "[]":
        return []
    try:
        parsed = json.loads(s)
    except json.JSONDecodeError:
        return []
    if not isinstance(parsed, list):
        return []
    out: list[tuple[int, int]] = []
    for item in parsed:
        if not isinstance(item, (list, tuple)) or len(item) < 3:
            continue
        try:
            out.append((int(item[1]), int(item[2])))
        except (TypeError, ValueError):
            continue
    return out
