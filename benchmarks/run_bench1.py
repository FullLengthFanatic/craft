#!/usr/bin/env python
"""Run the full Bench 1 grid (4 rates x 3 orientations x 3 seeds = 36 cells).

Saves per-cell scored DataFrames to ``benchmarks/cache/bench1_scores/`` so the
script can be resumed without re-running completed cells. The GENCODE pool is
cached as a pickle on first run.

Run from the repo root:

    PYTHONPATH=benchmarks .venv/bin/python benchmarks/run_bench1.py
"""

from __future__ import annotations

import pickle
import shutil
import sys
import time
from pathlib import Path

from cbench.bench1 import progress, run_cell
from cbench.data import load_protein_coding_transcripts
from cbench.metrics import to_dataframe

REPO_ROOT = Path(__file__).resolve().parent.parent
GTF_PATH = Path("/data/scNoiseMeter/gencode.v45.annotation.gtf.gz")
GENOME_PATH = Path("/data/scNoiseMeter/GRCh38.primary_assembly.genome.fa.gz")
CRAFT_BIN = str(REPO_ROOT / ".venv/bin/craft")
POOL_CACHE = REPO_ROOT / "benchmarks/cache/gencode_v45_pool.pkl"
SCORES_DIR = REPO_ROOT / "benchmarks/cache/bench1_scores"
WORKDIR = Path("/tmp/bench1_workdir")

RATES = (0.05, 0.10, 0.25, 0.50)
ORIENTATIONS = ("5prime", "3prime", "both")
SEEDS = (0, 1, 2)
N_PER_CELL = 3000


def load_or_cache_pool():
    if POOL_CACHE.exists():
        progress(f"[bench1] loading cached pool: {POOL_CACHE}")
        with POOL_CACHE.open("rb") as fh:
            return pickle.load(fh)
    progress("[bench1] parsing GENCODE v45 protein-coding transcripts (slow)...")
    t0 = time.time()
    pool = load_protein_coding_transcripts(
        GTF_PATH,
        GENOME_PATH,
        min_cds_bp=150,
        max_n=None,
        seed=0,
    )
    progress(f"[bench1] parsed {len(pool)} transcripts in {time.time() - t0:.1f}s")
    POOL_CACHE.parent.mkdir(parents=True, exist_ok=True)
    with POOL_CACHE.open("wb") as fh:
        pickle.dump(pool, fh)
    progress(f"[bench1] cached pool to {POOL_CACHE}")
    return pool


def main() -> int:
    pool = load_or_cache_pool()
    SCORES_DIR.mkdir(parents=True, exist_ok=True)

    cells = [
        (rate, orientation, seed)
        for rate in RATES
        for orientation in ORIENTATIONS
        for seed in SEEDS
    ]
    n_cells = len(cells)
    t0 = time.time()
    for idx, (rate, orientation, seed) in enumerate(cells, start=1):
        tag = f"r{int(rate * 100):02d}_{orientation}_s{seed}"
        out_path = SCORES_DIR / f"{tag}.tsv"
        if out_path.exists():
            progress(f"[bench1] {idx}/{n_cells} {tag}: SKIP (cached)")
            continue
        cell_start = time.time()
        shutil.rmtree(WORKDIR, ignore_errors=True)
        result = run_cell(
            transcripts=pool,
            rate=rate,
            orientation=orientation,
            seed=seed,
            genome_path=GENOME_PATH,
            workdir=WORKDIR,
            n_per_cell=N_PER_CELL,
            craft_bin=CRAFT_BIN,
        )
        df = to_dataframe(result.rows)
        df.to_csv(out_path, sep="\t", index=False)
        cell_elapsed = time.time() - cell_start
        total_elapsed = time.time() - t0
        progress(
            f"[bench1] {idx}/{n_cells} {tag}: "
            f"n_input={result.n_input} n_intact={result.n_intact_truth} "
            f"rows={len(result.rows)} cell={cell_elapsed:.1f}s "
            f"total={total_elapsed / 60:.1f}m"
        )

    progress(f"[bench1] DONE in {(time.time() - t0) / 60:.1f}m; scores under {SCORES_DIR}")
    shutil.rmtree(WORKDIR, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
