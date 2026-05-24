"""Orchestration for Bench 1: simulated truncation vs CRAFT propagation + orfipy.

A single cell of the benchmark is one ``(rate, orientation, seed)`` combination.
``run_cell`` does the per-cell work end-to-end:

1. Optionally subsample ``n_per_cell`` transcripts deterministically (by ``seed``).
2. Truncate each transcript's exons and transcript-orientation sequence.
3. Write a minimal iso GTF (truncated exons) and reference GTF (full GENCODE
   transcript with CDS + start/stop codons) to a temp dir.
4. Invoke ``craft annotate`` on the pair.
5. Invoke ``orfipy`` on the truncated transcript sequences.
6. Map both predictions to the truncated iso's transcript coordinates and
   score against the GENCODE-truth ORF (also re-expressed in those coordinates).

The ground-truth ORF is only considered ``intact`` when truncation didn't clip
either the start codon or the stop codon. Scoring is restricted to those rows;
transcripts whose ORF was partially clipped are reported as a separate count so
sampling bias is visible.
"""

import random
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

from cbench.comparators import (
    ORFCall,
    craft_propagated_calls,
    orfipy_predict,
)
from cbench.data import Transcript
from cbench.metrics import ORFScoreRow, score_one
from cbench.truncation import truncate_exons, truncate_sequence


@dataclass
class CellResult:
    rate: float
    orientation: str
    seed: int
    rows: list[ORFScoreRow] = field(default_factory=list)
    n_input: int = 0
    n_intact_truth: int = 0
    n_skipped_no_exons: int = 0


def _truth_in_truncated_coords(
    tx: Transcript, rate: float, orientation: str
) -> tuple[int, int] | None:
    n = tx.transcript_length
    trim = int(n * rate)
    if orientation == "5prime":
        offset_5p = trim
    elif orientation == "3prime":
        offset_5p = 0
    elif orientation == "both":
        offset_5p = trim // 2
    else:
        raise ValueError(orientation)
    truncated_length = n - trim
    new_start = tx.cds_tx_start - offset_5p
    new_end = tx.cds_tx_end - offset_5p
    if new_start < 0 or new_end > truncated_length or new_end <= new_start:
        return None
    return new_start, new_end


def _iso_id(tx: Transcript, rate: float, orientation: str, seed: int) -> str:
    return f"{tx.transcript_id}__r{int(rate * 100):02d}_{orientation}_s{seed}"


def _write_iso_gtf(
    path: Path,
    isos: list[tuple[str, Transcript, list[tuple[int, int]]]],
) -> None:
    """Write one transcript + exon set per iso, no CDS records (CRAFT propagates them)."""
    with path.open("w") as fh:
        for iso_id, tx, exons in isos:
            if not exons:
                continue
            attrs = f'transcript_id "{iso_id}"; gene_id "{tx.gene_id}";'
            start = min(s for s, _ in exons)
            end = max(e for _, e in exons)
            fh.write(
                f"{tx.chrom}\tcbench\ttranscript\t{start + 1}\t{end}\t.\t"
                f"{tx.strand}\t.\t{attrs}\n"
            )
            for e_start, e_end in sorted(exons):
                fh.write(
                    f"{tx.chrom}\tcbench\texon\t{e_start + 1}\t{e_end}\t.\t"
                    f"{tx.strand}\t.\t{attrs}\n"
                )


def _write_reference_gtf(path: Path, transcripts: list[Transcript]) -> None:
    """Write full GENCODE records (transcript + exons + CDS + start + stop) per tx."""
    with path.open("w") as fh:
        for tx in transcripts:
            attrs = (
                f'gene_id "{tx.gene_id}"; transcript_id "{tx.transcript_id}";'
            )
            t_start = min(s for s, _ in tx.exons)
            t_end = max(e for _, e in tx.exons)
            fh.write(
                f"{tx.chrom}\tcbench_ref\ttranscript\t{t_start + 1}\t{t_end}\t.\t"
                f"{tx.strand}\t.\t{attrs}\n"
            )
            for e_start, e_end in tx.exons:
                fh.write(
                    f"{tx.chrom}\tcbench_ref\texon\t{e_start + 1}\t{e_end}\t.\t"
                    f"{tx.strand}\t.\t{attrs}\n"
                )
            for c_start, c_end in tx.cds_genomic_intervals:
                fh.write(
                    f"{tx.chrom}\tcbench_ref\tCDS\t{c_start + 1}\t{c_end}\t.\t"
                    f"{tx.strand}\t0\t{attrs}\n"
                )
            s_start, s_end = tx.start_codon_genomic
            fh.write(
                f"{tx.chrom}\tcbench_ref\tstart_codon\t{s_start + 1}\t{s_end}\t.\t"
                f"{tx.strand}\t0\t{attrs}\n"
            )
            p_start, p_end = tx.stop_codon_genomic
            fh.write(
                f"{tx.chrom}\tcbench_ref\tstop_codon\t{p_start + 1}\t{p_end}\t.\t"
                f"{tx.strand}\t0\t{attrs}\n"
            )


def run_cell(
    transcripts: list[Transcript],
    rate: float,
    orientation: str,
    seed: int,
    genome_path: Path,
    workdir: Path,
    n_per_cell: int | None = None,
    craft_bin: str = "craft",
) -> CellResult:
    """Execute one cell of Bench 1 end-to-end and return scored rows.

    Args:
        transcripts: pool of GENCODE-truth transcripts to sample from.
        rate: truncation rate (0-1).
        orientation: ``"5prime"``, ``"3prime"``, or ``"both"``.
        seed: random seed; controls the per-cell subsample.
        genome_path: indexed genome FASTA.
        workdir: temp dir for GTFs and CRAFT output (created if missing).
        n_per_cell: subsample size per cell. ``None`` uses all transcripts.
        craft_bin: path to the ``craft`` CLI (default assumes it's on PATH).
    """
    workdir.mkdir(parents=True, exist_ok=True)
    pool = transcripts
    if n_per_cell is not None and len(pool) > n_per_cell:
        rng = random.Random(seed)
        pool = rng.sample(pool, n_per_cell)

    truncated_isos: list[tuple[str, Transcript, list[tuple[int, int]]]] = []
    iso_seqs: dict[str, str] = {}
    intact_truth: dict[str, tuple[int, int]] = {}
    truncated_tx_records: list[Transcript] = []
    n_skipped_no_exons = 0

    for tx in pool:
        truncated_exons = truncate_exons(tx.exons, tx.strand, rate, orientation)
        if not truncated_exons:
            n_skipped_no_exons += 1
            continue
        truth = _truth_in_truncated_coords(tx, rate, orientation)
        if truth is None:
            continue
        iso_id = _iso_id(tx, rate, orientation, seed)
        truncated_isos.append((iso_id, tx, truncated_exons))
        iso_seqs[iso_id] = truncate_sequence(tx.transcript_seq, rate, orientation)
        intact_truth[iso_id] = truth
        truncated_tx_records.append(
            Transcript(
                transcript_id=iso_id,
                gene_id=tx.gene_id,
                chrom=tx.chrom,
                strand=tx.strand,
                exons=truncated_exons,
                cds_genomic_intervals=[],
                start_codon_genomic=tx.start_codon_genomic,
                stop_codon_genomic=tx.stop_codon_genomic,
                transcript_seq=iso_seqs[iso_id],
                cds_tx_start=truth[0],
                cds_tx_end=truth[1],
                extras=dict(tx.extras),
            )
        )

    result = CellResult(
        rate=rate,
        orientation=orientation,
        seed=seed,
        n_input=len(pool),
        n_intact_truth=len(intact_truth),
        n_skipped_no_exons=n_skipped_no_exons,
    )
    if not intact_truth:
        return result

    iso_gtf = workdir / "iso.gtf"
    ref_gtf = workdir / "ref.gtf"
    craft_out = workdir / "craft_out"
    _write_iso_gtf(iso_gtf, truncated_isos)
    _write_reference_gtf(ref_gtf, pool)

    subprocess.run(
        [
            craft_bin,
            "annotate",
            "--isoforms",
            str(iso_gtf),
            "--reference",
            str(ref_gtf),
            "--genome",
            str(genome_path),
            "--output-dir",
            str(craft_out),
        ],
        check=True,
        capture_output=True,
    )

    craft_calls = craft_propagated_calls(
        craft_out / "per_isoform.tsv", truncated_tx_records
    )
    orfipy_calls = orfipy_predict(iso_seqs)

    rows: list[ORFScoreRow] = []
    for iso_id, (truth_start, truth_end) in intact_truth.items():
        seq = iso_seqs[iso_id]
        truth_call = ORFCall(
            transcript_id=iso_id,
            found=True,
            tx_start=truth_start,
            tx_end=truth_end,
            start_codon=seq[truth_start : truth_start + 3],
            stop_codon=seq[max(0, truth_end - 3) : truth_end],
        )
        for comparator, calls in (("craft", craft_calls), ("orfipy", orfipy_calls)):
            pred = calls.get(iso_id, ORFCall(iso_id, False))
            rows.append(
                score_one(
                    truth_call,
                    pred,
                    comparator,
                    rate,
                    orientation,
                    seed,
                    iso_id,
                )
            )
    result.rows = rows
    return result


def progress(*args, **kwargs) -> None:
    print(*args, file=sys.stderr, flush=True, **kwargs)
