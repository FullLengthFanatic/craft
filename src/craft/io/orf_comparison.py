"""Comparison with an independent GTF-native ORF caller such as ORFanage."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pyranges as pr

COLUMNS = [
    "transcript_id",
    "comparator_orf_present",
    "comparator_start_pos",
    "comparator_stop_codon_pos",
    "comparator_cds_bp",
    "comparator_start_agrees",
    "comparator_stop_agrees",
    "comparator_cds_bp_delta",
]


def _positions(group: pd.DataFrame) -> tuple[int, int, int]:
    strand = str(group["Strand"].iloc[0])
    cds = group[group["Feature"] == "CDS"]
    starts = group[group["Feature"] == "start_codon"]
    stops = group[group["Feature"] == "stop_codon"]
    if strand == "+":
        start = int(starts["Start"].min()) if not starts.empty else int(cds["Start"].min())
        stop = int(stops["Start"].min()) if not stops.empty else int(cds["End"].max())
    else:
        start = int(starts["End"].max()) - 1 if not starts.empty else int(cds["End"].max()) - 1
        stop = int(stops["End"].max()) - 1 if not stops.empty else int(cds["Start"].min()) - 1
    bp = int((cds["End"] - cds["Start"]).sum())
    return start, stop, bp


def compare_orf_gtf(per_isoform: pd.DataFrame, path: Path) -> pd.DataFrame:
    """Compare CRAFT calls with CDS/start/stop features in an independent GTF."""
    df = pr.read_gtf(str(path)).df
    if "transcript_id" not in df.columns:
        raise ValueError("ORF comparator GTF needs transcript_id attributes")
    df = df[df["Feature"].isin(["CDS", "start_codon", "stop_codon"])]
    calls: dict[str, tuple[int, int, int]] = {}
    for tx, group in df.groupby("transcript_id", sort=False):
        if (group["Feature"] == "CDS").any():
            calls[str(tx)] = _positions(group)

    rows: list[dict] = []
    for _, record in per_isoform.iterrows():
        tx = str(record["transcript_id"])
        if tx not in calls:
            rows.append({"transcript_id": tx, "comparator_orf_present": False})
            continue
        start, stop, bp = calls[tx]
        craft_start = record.get("resolved_start_pos")
        craft_stop = record.get("resolved_stop_codon_pos")
        craft_bp = record.get("resolved_cds_bp")
        rows.append(
            {
                "transcript_id": tx,
                "comparator_orf_present": True,
                "comparator_start_pos": start,
                "comparator_stop_codon_pos": stop,
                "comparator_cds_bp": bp,
                "comparator_start_agrees": (
                    bool(pd.notna(craft_start) and int(craft_start) == start)
                ),
                "comparator_stop_agrees": bool(pd.notna(craft_stop) and int(craft_stop) == stop),
                "comparator_cds_bp_delta": (
                    int(craft_bp) - bp if pd.notna(craft_bp) else None
                ),
            }
        )
    return pd.DataFrame(rows).reindex(columns=COLUMNS)

