"""Reference-transcript quality and parent-priority metadata.

Reference annotations contain incomplete CDS models and many biologically
equivalent transcript choices.  CRAFT keeps those facts explicit instead of
assuming that the first/longest overlapping CDS is a complete ground truth.
"""

from __future__ import annotations

import re

import pandas as pd
import pyranges as pr

_INCOMPLETE_TAGS = ("cds_start_NF", "cds_end_NF", "mRNA_start_NF", "mRNA_end_NF")


def _joined_values(group: pd.DataFrame, column: str) -> str:
    if column not in group.columns:
        return ""
    return ";".join(sorted({str(v) for v in group[column].dropna() if str(v) not in {"", "."}}))


def _has_token(text: str, token: str) -> bool:
    return bool(re.search(rf"(?:^|[;,\s]){re.escape(token)}(?:$|[;,\s])", text))


def _phase_consistent(cds: pd.DataFrame, strand: str) -> bool:
    """Check GTF CDS phase values and their continuity in transcript order."""
    if cds.empty or "Frame" not in cds.columns:
        return True
    ordered = cds.sort_values("Start", ascending=strand == "+", kind="stable")
    phases = pd.to_numeric(ordered["Frame"].replace(".", pd.NA), errors="coerce")
    if phases.isna().any() or not phases.isin([0, 1, 2]).all():
        return False
    phase_values = phases.astype(int).tolist()
    lengths = (ordered["End"] - ordered["Start"]).astype(int).tolist()
    initial_phase = phase_values[0]
    cumulative = lengths[0] - initial_phase
    for phase, length in zip(phase_values[1:], lengths[1:], strict=True):
        if phase != (3 - (cumulative % 3)) % 3:
            return False
        cumulative += length
    return True


def transcript_metadata(reference: pr.PyRanges) -> pd.DataFrame:
    """One row per reference transcript with CDS completeness and curation priority."""
    df = reference.df
    columns = [
        "transcript_id",
        "reference_gene_id",
        "reference_has_cds",
        "reference_cds_complete",
        "reference_has_explicit_start",
        "reference_has_explicit_stop",
        "reference_cds_phase_valid",
        "reference_priority",
        "reference_priority_reason",
        "reference_transcript_type",
        "reference_tags",
    ]
    if df.empty:
        return pd.DataFrame(columns=columns)

    rows: list[dict] = []
    for tx, group in df.groupby("transcript_id", sort=False):
        features = set(group["Feature"].astype(str))
        cds = group[group["Feature"] == "CDS"]
        tags = _joined_values(group, "tag")
        tx_type = _joined_values(group, "transcript_type") or _joined_values(
            group, "transcript_biotype"
        )
        incomplete = any(token in tags for token in _INCOMPLETE_TAGS)
        has_start = "start_codon" in features
        has_stop = "stop_codon" in features

        phase_valid = _phase_consistent(cds, str(group["Strand"].iloc[0]))

        priority = 0
        reasons: list[str] = []
        if _has_token(tags, "MANE_Select"):
            priority += 100
            reasons.append("MANE_Select")
        if _has_token(tags, "MANE_Plus_Clinical"):
            priority += 90
            reasons.append("MANE_Plus_Clinical")
        appris = next(
            (t for t in re.split(r"[;,\s]+", tags) if t.startswith("appris_principal")),
            "",
        )
        if appris:
            priority += 70
            reasons.append(appris)
        if _has_token(tags, "CCDS") or ("ccdsid" in group and group["ccdsid"].notna().any()):
            priority += 40
            reasons.append("CCDS")
        if _has_token(tags, "basic"):
            priority += 10
            reasons.append("GENCODE_basic")
        if incomplete:
            priority -= 100
            reasons.append("incomplete_CDS")
        if not phase_valid:
            priority -= 50
            reasons.append("phase_inconsistent")
        if not cds.empty and not has_start:
            priority -= 20
            reasons.append("no_explicit_start")
        if not cds.empty and not has_stop:
            priority -= 20
            reasons.append("no_explicit_stop")

        gene_id = ""
        if "gene_id" in group.columns and group["gene_id"].notna().any():
            gene_id = str(group["gene_id"].dropna().iloc[0])
        rows.append(
            {
                "transcript_id": str(tx),
                "reference_gene_id": gene_id,
                "reference_has_cds": not cds.empty,
                "reference_cds_complete": bool(
                    not cds.empty
                    and not incomplete
                    and has_start
                    and has_stop
                    and phase_valid
                ),
                "reference_has_explicit_start": has_start,
                "reference_has_explicit_stop": has_stop,
                "reference_cds_phase_valid": phase_valid,
                "reference_priority": priority,
                "reference_priority_reason": ";".join(reasons) or "unranked",
                "reference_transcript_type": tx_type,
                "reference_tags": tags,
            }
        )
    return pd.DataFrame(rows, columns=columns)
