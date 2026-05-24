#!/usr/bin/env python
"""Bench 2: real-data ORF concordance on bcM0003 sc PacBio Iso-Seq.

For each iso in CRAFT's bcM0003 output that has a GENCODE parent with a
complete CDS, project the parent's start/stop codon genomic positions onto
the iso's exon structure (= the "bulk truth" ORF in iso transcript coords).
Run orfipy on the iso's transcript-orientation sequence to get a de novo
ORF call. Compare CRAFT's propagated ORF (from per_isoform.tsv) and
orfipy's de novo call against the truth, restricted to isos whose truth
ORF is fully intact in the truncated iso.

Stratify the resulting concordance rate by CRAFT's completeness category
so the figure shows where propagation buys you something over de novo on
real-world long-read truncation patterns.

Run from the repo root:

    PYTHONPATH=benchmarks .venv/bin/python benchmarks/run_bench2.py
"""

from __future__ import annotations

import gzip
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pysam
from plotly.subplots import make_subplots

from cbench.comparators import ORFCall, orfipy_predict
from cbench.data import _genomic_to_tx_coord, _rc
from cbench.metrics import ORFScoreRow, score_one, to_dataframe

REPO_ROOT = Path(__file__).resolve().parent.parent
ISO_GFF = REPO_ROOT / "test_dataset/bcM0003_bcM0003.collapsed.sorted.gff"
PER_ISO_TSV = (
    REPO_ROOT
    / "test_dataset/bcM0003_full_genome/craft_out_atlas_filtered/per_isoform.tsv"
)
GENOME = Path("/data/scNoiseMeter/GRCh38.primary_assembly.genome.fa.gz")
POOL_PKL = REPO_ROOT / "benchmarks/cache/gencode_v45_pool.pkl"

CACHE = REPO_ROOT / "benchmarks/cache/bench2"
EXONS_PKL = CACHE / "bcm0003_iso_exons.pkl"
SCORES_TSV = CACHE / "bench2_scores.tsv.gz"
FIG_PNG = REPO_ROOT / "benchmarks/figures/bench2_concordance_panel.png"
FIG_JSON = REPO_ROOT / "benchmarks/figures/bench2_concordance_panel.json"


def progress(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def parse_iso_exons() -> dict[str, dict]:
    """Stream the bcM0003 GFF, yield exon lists keyed by transcript_id.

    Each value is {chrom, strand, exons: [(start, end), ...]} with 0-based
    half-open coordinates, sorted by start.
    """
    if EXONS_PKL.exists():
        progress(f"[exons] cached: {EXONS_PKL}")
        with EXONS_PKL.open("rb") as fh:
            return pickle.load(fh)
    progress(f"[exons] parsing {ISO_GFF}")
    t0 = time.time()
    iso: dict[str, dict] = {}
    with ISO_GFF.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            attrs = cols[8]
            tx_id = ""
            for chunk in attrs.split(";"):
                chunk = chunk.strip()
                if chunk.startswith("transcript_id"):
                    tx_id = chunk.split('"', 2)[1] if '"' in chunk else chunk.split(" ", 1)[1]
                    break
            if not tx_id:
                continue
            chrom = cols[0]
            strand = cols[6]
            start = int(cols[3]) - 1
            end = int(cols[4])
            rec = iso.get(tx_id)
            if rec is None:
                iso[tx_id] = {"chrom": chrom, "strand": strand, "exons": [(start, end)]}
            else:
                rec["exons"].append((start, end))
    for rec in iso.values():
        rec["exons"].sort()
    progress(f"[exons] {len(iso):,} isoforms parsed in {time.time() - t0:.1f}s")
    CACHE.mkdir(parents=True, exist_ok=True)
    with EXONS_PKL.open("wb") as fh:
        pickle.dump(iso, fh)
    return iso


def truth_in_iso_coords(
    parent_start_codon: tuple[int, int],
    parent_stop_codon: tuple[int, int],
    iso_exons: list[tuple[int, int]],
    iso_strand: str,
) -> tuple[int, int] | None:
    """Project parent's start_codon + stop_codon genomic intervals onto the iso.

    Returns ``(cds_tx_start, cds_tx_end)`` in iso transcript coords (CDS
    excludes the stop codon to match GENCODE), or ``None`` if either codon is
    missing from the iso's exons.

    Parent codons are 3 bp half-open ``(g_start, g_end)`` intervals in
    genomic coords. On ``+`` strand the A of ATG sits at ``g_start``; on
    ``-`` strand it sits at ``g_end - 1``.
    """

    def _genomic_inside(pos: int) -> bool:
        return any(s <= pos < e for s, e in iso_exons)

    if iso_strand == "+":
        start_pos = parent_start_codon[0]
        stop_pos = parent_stop_codon[0]
    elif iso_strand == "-":
        start_pos = parent_start_codon[1] - 1
        stop_pos = parent_stop_codon[1] - 1
    else:
        return None
    if not _genomic_inside(start_pos) or not _genomic_inside(stop_pos):
        return None
    try:
        cds_tx_start = _genomic_to_tx_coord(start_pos, iso_exons, iso_strand)
        cds_tx_end = _genomic_to_tx_coord(stop_pos, iso_exons, iso_strand)
    except ValueError:
        return None
    if cds_tx_end <= cds_tx_start:
        return None
    return cds_tx_start, cds_tx_end


def iso_transcript_sequence(
    chrom: str, strand: str, exons: list[tuple[int, int]], genome: pysam.FastaFile
) -> str:
    if chrom not in genome.references:
        return ""
    pieces = [genome.fetch(chrom, s, e).upper() for s, e in exons]
    seq = "".join(pieces)
    return _rc(seq) if strand == "-" else seq


def parse_craft_intervals(raw) -> list[tuple[int, int]]:
    if raw is None or pd.isna(raw):
        return []
    s = str(raw).strip()
    if not s or s == "[]":
        return []
    try:
        parsed = json.loads(s)
    except json.JSONDecodeError:
        return []
    out: list[tuple[int, int]] = []
    for item in parsed:
        if isinstance(item, (list, tuple)) and len(item) >= 3:
            try:
                out.append((int(item[1]), int(item[2])))
            except (TypeError, ValueError):
                continue
    return out


def craft_orf_call(
    iso_id: str,
    craft_row: pd.Series,
    iso_exons: list[tuple[int, int]],
    iso_strand: str,
    iso_seq: str,
) -> ORFCall:
    outcome = str(craft_row.get("orf_outcome", ""))
    if outcome in ("no_parent", "no_parent_cds", "start_lost"):
        return ORFCall(iso_id, False)
    intervals = parse_craft_intervals(craft_row.get("propagated_cds_intervals", ""))
    if not intervals:
        return ORFCall(iso_id, False)
    try:
        tx_positions: list[int] = []
        for g_start, g_end in intervals:
            for pos in (g_start, g_end - 1):
                tx_positions.append(_genomic_to_tx_coord(pos, iso_exons, iso_strand))
    except ValueError:
        return ORFCall(iso_id, False)
    if not tx_positions:
        return ORFCall(iso_id, False)
    tx_lo, tx_hi = min(tx_positions), max(tx_positions) + 1
    return ORFCall(
        transcript_id=iso_id,
        found=True,
        tx_start=tx_lo,
        tx_end=tx_hi,
        start_codon=iso_seq[tx_lo : tx_lo + 3] if 0 <= tx_lo < len(iso_seq) else "",
        stop_codon=iso_seq[max(0, tx_hi - 3) : tx_hi] if tx_hi <= len(iso_seq) else "",
    )


def main() -> int:
    progress("[bench2] loading GENCODE pool")
    with POOL_PKL.open("rb") as fh:
        pool = pickle.load(fh)
    pool_by_id = {t.transcript_id: t for t in pool}
    progress(f"[bench2] {len(pool_by_id):,} GENCODE transcripts in pool")

    iso_exons_by_id = parse_iso_exons()

    progress(f"[bench2] reading {PER_ISO_TSV}")
    craft = pd.read_csv(
        PER_ISO_TSV,
        sep="\t",
        usecols=[
            "transcript_id",
            "completeness",
            "parent_tx_id",
            "orf_outcome",
            "orf_confidence",
            "propagated_cds_intervals",
        ],
        dtype={"transcript_id": "string", "parent_tx_id": "string"},
    ).set_index("transcript_id")
    progress(f"[bench2] {len(craft):,} CRAFT rows")

    progress("[bench2] building truth + extracting sequences (eligible isos)")
    genome = pysam.FastaFile(str(GENOME))
    t0 = time.time()
    eligible: list[dict] = []
    seqs: dict[str, str] = {}
    for iso_id, craft_row in craft.iterrows():
        parent_id = craft_row["parent_tx_id"]
        if pd.isna(parent_id):
            continue
        parent_id = str(parent_id)
        if not parent_id:
            continue
        parent = pool_by_id.get(parent_id)
        if parent is None:
            continue
        iso_rec = iso_exons_by_id.get(iso_id)
        if iso_rec is None:
            continue
        iso_strand = iso_rec["strand"]
        if iso_strand != parent.strand or iso_rec["chrom"] != parent.chrom:
            continue
        truth = truth_in_iso_coords(
            parent.start_codon_genomic,
            parent.stop_codon_genomic,
            iso_rec["exons"],
            iso_strand,
        )
        if truth is None:
            continue
        seq = iso_transcript_sequence(iso_rec["chrom"], iso_strand, iso_rec["exons"], genome)
        if not seq or len(seq) < 60:
            continue
        eligible.append(
            {
                "iso_id": iso_id,
                "iso_strand": iso_strand,
                "iso_exons": iso_rec["exons"],
                "iso_seq": seq,
                "truth_start": truth[0],
                "truth_end": truth[1],
                "completeness": str(craft_row["completeness"]),
                "orf_confidence": str(craft_row["orf_confidence"]),
                "craft_row": craft_row,
            }
        )
        seqs[iso_id] = seq
    genome.close()
    progress(
        f"[bench2] {len(eligible):,} eligible isos (intact GENCODE truth) "
        f"prepared in {(time.time() - t0) / 60:.1f}m"
    )

    progress("[bench2] running orfipy on eligible iso sequences (batched)")
    t0 = time.time()
    orfipy_calls = orfipy_predict(seqs, min_orf_nt=75)
    progress(f"[bench2] orfipy done in {(time.time() - t0) / 60:.1f}m")

    progress("[bench2] scoring")
    rows: list[ORFScoreRow] = []
    for rec in eligible:
        truth_call = ORFCall(
            transcript_id=rec["iso_id"],
            found=True,
            tx_start=rec["truth_start"],
            tx_end=rec["truth_end"],
            start_codon=rec["iso_seq"][rec["truth_start"] : rec["truth_start"] + 3],
            stop_codon=rec["iso_seq"][max(0, rec["truth_end"] - 3) : rec["truth_end"]],
        )
        craft_call = craft_orf_call(
            rec["iso_id"], rec["craft_row"], rec["iso_exons"], rec["iso_strand"], rec["iso_seq"]
        )
        orfipy_call = orfipy_calls.get(rec["iso_id"], ORFCall(rec["iso_id"], False))
        for comparator, pred in (("craft", craft_call), ("orfipy", orfipy_call)):
            row = score_one(
                truth_call,
                pred,
                comparator,
                rate=0.0,
                orientation="",
                seed=0,
                transcript_id=rec["iso_id"],
            )
            rows.append(row)

    df = to_dataframe(rows)
    completeness_lookup = {rec["iso_id"]: rec["completeness"] for rec in eligible}
    confidence_lookup = {rec["iso_id"]: rec["orf_confidence"] for rec in eligible}
    df["completeness"] = df["transcript_id"].map(completeness_lookup)
    df["orf_confidence"] = df["transcript_id"].map(confidence_lookup)
    CACHE.mkdir(parents=True, exist_ok=True)
    df.to_csv(SCORES_TSV, sep="\t", index=False, compression="gzip")
    progress(f"[bench2] wrote {SCORES_TSV} ({len(df):,} rows)")

    summary = (
        df.groupby(["completeness", "comparator"], dropna=False)
        .agg(
            n=("found", "size"),
            recovery_rate=("found", "mean"),
            start_exact_rate=("start_exact", "mean"),
            stop_exact_rate=("stop_exact", "mean"),
            mean_abs_length_error=("length_error_nt", lambda s: float(s.abs().mean())),
        )
        .reset_index()
    )
    pivot_start = summary.pivot(
        index="completeness", columns="comparator", values="start_exact_rate"
    )
    pivot_n = summary.pivot(index="completeness", columns="comparator", values="n")
    progress("[bench2] per-completeness start-exact rate:")
    progress(str(pivot_start.round(3)))
    progress("[bench2] per-completeness n:")
    progress(str(pivot_n.astype(int)))

    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=(
            "Start-codon exact match by CRAFT completeness",
            "CRAFT start-exact rate by orf_confidence",
        ),
        column_widths=[0.62, 0.38],
        horizontal_spacing=0.16,
    )
    ordered_cats = [
        c
        for c in (
            "full_length",
            "alt_3prime_end",
            "truncated_5p",
            "truncated_3p",
            "truncated_both",
            "internal_fragment",
        )
        if c in pivot_start.index
    ]
    fig.add_trace(
        go.Bar(
            x=ordered_cats,
            y=[pivot_start.loc[c, "craft"] for c in ordered_cats],
            text=[f"n={int(pivot_n.loc[c, 'craft'])}" for c in ordered_cats],
            textposition="outside",
            marker_color="#5b7a9d",
            name="CRAFT propagation",
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Bar(
            x=ordered_cats,
            y=[pivot_start.loc[c, "orfipy"] for c in ordered_cats],
            marker_color="#c97a5b",
            name="orfipy de novo",
        ),
        row=1,
        col=1,
    )
    fig.update_xaxes(title_text="CRAFT completeness", tickangle=-25, row=1, col=1)
    fig.update_yaxes(title_text="Start-codon exact-match rate", range=[0, 1.05], row=1, col=1)

    craft_only = df[df["comparator"] == "craft"]
    by_conf = (
        craft_only.groupby("orf_confidence")
        .agg(start_exact_rate=("start_exact", "mean"), n=("found", "size"))
        .reset_index()
    )
    conf_order = [c for c in ("high", "medium", "low", "none") if c in by_conf["orf_confidence"].values]
    by_conf = by_conf.set_index("orf_confidence").loc[conf_order].reset_index()
    fig.add_trace(
        go.Bar(
            x=by_conf["orf_confidence"],
            y=by_conf["start_exact_rate"],
            text=[f"n={int(n)}" for n in by_conf["n"]],
            textposition="outside",
            marker_color="#5b7a9d",
            showlegend=False,
        ),
        row=1,
        col=2,
    )
    fig.update_xaxes(title_text="CRAFT orf_confidence", row=1, col=2)
    fig.update_yaxes(title_text="Start-codon exact-match rate", range=[0, 1.05], row=1, col=2)

    fig.update_layout(
        height=480,
        width=1200,
        barmode="group",
        template="plotly_white",
        title_text=(
            f"Bench 2: real-data ORF concordance on bcM0003 sc PacBio Iso-Seq "
            f"(eligible isos: {len(eligible):,})"
        ),
        legend=dict(orientation="h", y=-0.22, x=0.5, xanchor="center"),
    )
    fig.write_image(str(FIG_PNG), scale=2)
    FIG_JSON.write_text(fig.to_json())
    progress(f"[bench2] PNG: {FIG_PNG} ({FIG_PNG.stat().st_size / 1024:.1f} KB)")
    progress(f"[bench2] JSON: {FIG_JSON}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
