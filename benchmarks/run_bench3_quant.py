#!/usr/bin/env python
"""Bench 3 step 2: salmon transcript-level quant of GSE86148 UPF1-KD vs control.

Samples (GSE86148, SRP083135, HeLa, siRNA, Lykke-Andersen 2017):

    Control (scrambled):  SRR4081222, SRR4081223, SRR4081224
    UPF1 knockdown:       SRR4081225, SRR4081226, SRR4081227

Pipeline per sample:
    prefetch -> fasterq-dump -> salmon quant -> delete fastq

Peak disk: ~7 GB (one sample's FASTQ at a time). Salmon index ~3 GB sticks
around for reuse. Total wall: ~1.5 h on this VM.

Run from the repo root:

    PYTHONPATH=benchmarks .venv/bin/python benchmarks/run_bench3_quant.py
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SALMON_BIN = str(REPO_ROOT / ".tools/salmon-linux-x86_64/bin/salmon")
CACHE = REPO_ROOT / "benchmarks/cache/bench3"
TX_FASTA = CACHE / "gencode.v45.transcripts.fa.gz"
INDEX_DIR = CACHE / "salmon_index_v45"
QUANT_DIR = CACHE / "quants"
SRA_DIR = CACHE / "sra"
FASTQ_DIR = CACHE / "fastq"

TX_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/"
    "gencode.v45.transcripts.fa.gz"
)

SAMPLES = {
    "control": ["SRR4081222", "SRR4081223", "SRR4081224"],
    "upf1_kd": ["SRR4081225", "SRR4081226", "SRR4081227"],
}


def progress(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def fetch_transcriptome() -> None:
    if TX_FASTA.exists():
        progress(f"[quant] transcriptome cached: {TX_FASTA}")
        return
    progress(f"[quant] downloading {TX_URL}")
    TX_FASTA.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with urllib.request.urlopen(TX_URL) as resp, TX_FASTA.open("wb") as out:
        shutil.copyfileobj(resp, out, length=1 << 20)
    progress(f"[quant] saved {TX_FASTA.stat().st_size / 1e6:.1f} MB in {time.time() - t0:.1f}s")


def build_index() -> None:
    sentinel = INDEX_DIR / "info.json"
    if sentinel.exists():
        progress(f"[quant] salmon index cached: {INDEX_DIR}")
        return
    progress(f"[quant] building salmon index -> {INDEX_DIR}")
    INDEX_DIR.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    subprocess.run(
        [
            SALMON_BIN,
            "index",
            "-t",
            str(TX_FASTA),
            "-i",
            str(INDEX_DIR),
            "-k",
            "31",
            "--threads",
            "8",
            "--gencode",
        ],
        check=True,
    )
    progress(f"[quant] index built in {(time.time() - t0) / 60:.1f}m")


def fastqs_for(srr: str) -> tuple[Path, Path] | tuple[Path]:
    pe1 = FASTQ_DIR / f"{srr}_1.fastq"
    pe2 = FASTQ_DIR / f"{srr}_2.fastq"
    se = FASTQ_DIR / f"{srr}.fastq"
    if pe1.exists() and pe2.exists():
        return pe1, pe2
    if se.exists():
        return (se,)
    return ()


def fetch_fastq(srr: str) -> None:
    """Download SRA + dump FASTQ. fasterq-dump auto-splits PE into _1/_2."""
    existing = fastqs_for(srr)
    if existing:
        progress(f"[quant] {srr}: fastq cached")
        return
    SRA_DIR.mkdir(parents=True, exist_ok=True)
    FASTQ_DIR.mkdir(parents=True, exist_ok=True)
    progress(f"[quant] {srr}: prefetch")
    t0 = time.time()
    subprocess.run(
        ["prefetch", "--output-directory", str(SRA_DIR), "--max-size", "30g", srr],
        check=True,
    )
    progress(f"[quant] {srr}: fasterq-dump")
    subprocess.run(
        [
            "fasterq-dump",
            "--split-files",
            "--threads",
            "8",
            "--outdir",
            str(FASTQ_DIR),
            str(SRA_DIR / srr / f"{srr}.sra"),
        ],
        check=True,
    )
    shutil.rmtree(SRA_DIR / srr, ignore_errors=True)
    progress(f"[quant] {srr}: fastq ready in {(time.time() - t0) / 60:.1f}m")


def quant_one(srr: str) -> None:
    quant_out = QUANT_DIR / srr
    if (quant_out / "quant.sf").exists():
        progress(f"[quant] {srr}: quant cached")
        return
    fastqs = fastqs_for(srr)
    if not fastqs:
        raise RuntimeError(f"no fastq found for {srr}")
    QUANT_DIR.mkdir(parents=True, exist_ok=True)
    progress(f"[quant] {srr}: salmon quant")
    t0 = time.time()
    cmd = [
        SALMON_BIN,
        "quant",
        "-i",
        str(INDEX_DIR),
        "-l",
        "A",
        "-o",
        str(quant_out),
        "--threads",
        "8",
        "--validateMappings",
        "--gcBias",
    ]
    if len(fastqs) == 2:
        cmd += ["-1", str(fastqs[0]), "-2", str(fastqs[1])]
    else:
        cmd += ["-r", str(fastqs[0])]
    subprocess.run(cmd, check=True, capture_output=True)
    for fq in fastqs:
        fq.unlink(missing_ok=True)
    progress(f"[quant] {srr}: quant done in {(time.time() - t0) / 60:.1f}m")


def main() -> int:
    fetch_transcriptome()
    build_index()
    flat = [(srr, cond) for cond, runs in SAMPLES.items() for srr in runs]
    total = len(flat)
    t0 = time.time()
    for idx, (srr, cond) in enumerate(flat, start=1):
        progress(f"[quant] sample {idx}/{total}: {srr} ({cond})")
        fetch_fastq(srr)
        quant_one(srr)
        elapsed = (time.time() - t0) / 60
        progress(f"[quant] sample {idx}/{total} done; cumulative {elapsed:.1f}m")
    progress(f"[quant] all samples done in {(time.time() - t0) / 60:.1f}m")
    progress(f"[quant] quants under {QUANT_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
