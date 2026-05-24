#!/usr/bin/env python
"""Bench 3 step 1: build the CRAFT NMD label universe.

Filter GENCODE v45 to transcripts with type ``protein_coding`` or
``nonsense_mediated_decay`` AND a complete CDS (both start and stop codons
annotated). Write an exon-only iso GTF; run ``craft annotate`` against the
full GENCODE GTF as reference; cache per-transcript NMD labels for the
downstream UPF1-KD enrichment analysis.

Run from the repo root:

    PYTHONPATH=benchmarks .venv/bin/python benchmarks/run_bench3_universe.py
"""

from __future__ import annotations

import gzip
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
GTF_PATH = Path("/data/scNoiseMeter/gencode.v45.annotation.gtf.gz")
GENOME_PATH = Path("/data/scNoiseMeter/GRCh38.primary_assembly.genome.fa.gz")
CRAFT_BIN = str(REPO_ROOT / ".venv/bin/craft")

CACHE_DIR = REPO_ROOT / "benchmarks/cache/bench3"
ISO_GTF = CACHE_DIR / "gencode_v45_nmd_universe.iso.gtf"
CRAFT_OUT = CACHE_DIR / "craft_out_universe"
LABEL_TSV = CACHE_DIR / "nmd_labels.tsv.gz"

KEEP_TYPES = {"protein_coding", "nonsense_mediated_decay"}


def progress(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def parse_attrs(attr_str: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for chunk in attr_str.strip().rstrip(";").split(";"):
        chunk = chunk.strip()
        if not chunk or " " not in chunk:
            continue
        k, v = chunk.split(" ", 1)
        out[k] = v.strip().strip('"')
    return out


def build_iso_gtf() -> dict[str, int]:
    """Stream GENCODE -> filtered exon-only iso GTF.

    First pass: collect transcript_ids whose transcript record has
    ``transcript_type`` in ``KEEP_TYPES`` AND that have both ``start_codon`` and
    ``stop_codon`` features. Second pass: emit exon records for those tx.
    """
    keep: set[str] = set()
    has_start: set[str] = set()
    has_stop: set[str] = set()
    typed_in_scope: set[str] = set()

    t0 = time.time()
    with gzip.open(GTF_PATH, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            feature = cols[2]
            if feature not in ("transcript", "start_codon", "stop_codon"):
                continue
            attrs = parse_attrs(cols[8])
            tx_id = attrs.get("transcript_id")
            if tx_id is None:
                continue
            if feature == "transcript":
                if attrs.get("transcript_type") in KEEP_TYPES:
                    typed_in_scope.add(tx_id)
            elif feature == "start_codon":
                has_start.add(tx_id)
            elif feature == "stop_codon":
                has_stop.add(tx_id)
    keep = typed_in_scope & has_start & has_stop
    progress(
        f"[universe] candidates in scope: {len(typed_in_scope)}, "
        f"with start_codon: {len(has_start)}, with stop_codon: {len(has_stop)}, "
        f"kept (intersection): {len(keep)}; first pass {time.time() - t0:.1f}s"
    )

    counts = {"transcripts_in": 0, "exons_out": 0}
    ISO_GTF.parent.mkdir(parents=True, exist_ok=True)
    t1 = time.time()
    with gzip.open(GTF_PATH, "rt") as fh, ISO_GTF.open("w") as out:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            attrs = parse_attrs(cols[8])
            tx_id = attrs.get("transcript_id")
            if tx_id is None or tx_id not in keep:
                continue
            # Minimal exon record: keep transcript_id + gene_id for CRAFT loader.
            slim_attrs = f'transcript_id "{tx_id}";'
            gene_id = attrs.get("gene_id")
            if gene_id:
                slim_attrs += f' gene_id "{gene_id}";'
            out.write(
                f"{cols[0]}\tcbench\texon\t{cols[3]}\t{cols[4]}\t.\t"
                f"{cols[6]}\t.\t{slim_attrs}\n"
            )
            counts["exons_out"] += 1
    counts["transcripts_in"] = len(keep)
    progress(
        f"[universe] wrote {counts['exons_out']} exon records "
        f"({counts['transcripts_in']} transcripts) to {ISO_GTF}; "
        f"second pass {time.time() - t1:.1f}s"
    )
    return counts


def run_craft() -> None:
    progress(f"[universe] launching craft annotate -> {CRAFT_OUT}")
    t0 = time.time()
    subprocess.run(
        [
            CRAFT_BIN,
            "annotate",
            "--isoforms",
            str(ISO_GTF),
            "--reference",
            str(GTF_PATH),
            "--genome",
            str(GENOME_PATH),
            "--output-dir",
            str(CRAFT_OUT),
        ],
        check=True,
    )
    progress(f"[universe] craft annotate done in {(time.time() - t0) / 60:.1f}m")


def extract_labels() -> None:
    import pandas as pd

    tsv = CRAFT_OUT / "per_isoform.tsv"
    progress(f"[universe] parsing {tsv}")
    df = pd.read_csv(
        tsv,
        sep="\t",
        usecols=[
            "transcript_id",
            "parent_tx_id",
            "parent_gene_id",
            "parent_gene_name",
            "completeness",
            "orf_outcome",
            "orf_confidence",
            "nmd_status",
            "nmd_rule",
            "nmd_confidence",
            "stop_to_last_junction_nt",
            "last_exon_length_nt",
        ],
        dtype={"transcript_id": "string"},
    )
    progress(
        f"[universe] {len(df)} rows; nmd_status: "
        f"{df.nmd_status.value_counts().to_dict()}"
    )
    df.to_csv(LABEL_TSV, sep="\t", index=False, compression="gzip")
    progress(f"[universe] wrote {LABEL_TSV}")


def main() -> int:
    if not LABEL_TSV.exists():
        if not ISO_GTF.exists():
            build_iso_gtf()
        if not (CRAFT_OUT / "per_isoform.tsv").exists():
            run_craft()
        extract_labels()
    else:
        progress(f"[universe] cache exists at {LABEL_TSV}; nothing to do")
    return 0


if __name__ == "__main__":
    sys.exit(main())
