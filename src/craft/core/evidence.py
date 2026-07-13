"""Transparent molecule/read evidence aggregation for isoform plausibility.

The score here is deliberately called an evidence score, not a probability.
It combines only measured features and reports how many were available.  A
future dataset-specific calibrator can replace this model without changing the
input/output contract.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

POSITIVE_WEIGHTS = {
    "unique_molecule_fraction": 2.0,
    "canonical_junction_fraction": 1.5,
    "short_read_junction_support_fraction": 1.0,
    "full_length_fraction": 1.0,
    "five_prime_adapter_fraction": 0.75,
    "polya_tail_fraction": 1.0,
    "mapq_fraction": 0.5,
    "replicate_support_fraction": 1.0,
}

NEGATIVE_WEIGHTS = {
    "ambiguous_molecule_fraction": 1.25,
    "internal_priming_fraction": 2.0,
    "template_switch_fraction": 1.5,
    "chimera_fraction": 2.0,
}

OUTPUT_COLUMNS = [
    "transcript_id",
    "isoform_evidence_score",
    "isoform_evidence_tier",
    "evidence_feature_count",
    "evidence_model",
    "evidence_warnings",
    *POSITIVE_WEIGHTS,
    *NEGATIVE_WEIGHTS,
]


def _read_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python")
    id_col = next(
        (c for c in ("transcript_id", "isoform", "pbid", "id") if c in df.columns),
        None,
    )
    if id_col is None:
        raise ValueError("Evidence table needs transcript_id/isoform/pbid/id")
    return df.rename(columns={id_col: "transcript_id"}).drop_duplicates("transcript_id")


def _numeric(df: pd.DataFrame, names: tuple[str, ...]) -> pd.Series:
    column = next((name for name in names if name in df.columns), None)
    if column is None:
        return pd.Series(np.nan, index=df.index, dtype=float)
    return pd.to_numeric(df[column], errors="coerce")


def _fraction(series: pd.Series) -> pd.Series:
    # Accept either [0, 1] fractions or [0, 100] percentages.
    converted = series.where(series.le(1.0), series / 100.0)
    return converted.clip(0.0, 1.0)


def load_evidence(path: Path) -> pd.DataFrame:
    """Load scNoiseMeter/tecap/caller-derived per-isoform evidence columns."""
    raw = _read_table(path)
    out = pd.DataFrame({"transcript_id": raw["transcript_id"].astype(str)})
    total = _numeric(raw, ("total_molecules", "total_count", "molecules")).replace(0, np.nan)
    unique = _numeric(raw, ("unique_molecules", "unique_support", "unique_reads"))
    ambiguous = _numeric(raw, ("ambiguous_molecules", "ambiguous_support", "ambiguous_reads"))
    out["unique_molecule_fraction"] = _fraction(
        _numeric(raw, ("unique_molecule_fraction", "unique_fraction")).fillna(unique / total)
    )
    out["ambiguous_molecule_fraction"] = _fraction(
        _numeric(raw, ("ambiguous_molecule_fraction", "ambiguous_fraction")).fillna(
            ambiguous / total
        )
    )
    aliases = {
        "canonical_junction_fraction": ("canonical_junction_fraction", "all_canonical_fraction"),
        "short_read_junction_support_fraction": (
            "short_read_junction_support_fraction", "short_read_support_fraction"
        ),
        "full_length_fraction": ("full_length_fraction", "fl_fraction"),
        "five_prime_adapter_fraction": (
            "five_prime_adapter_fraction", "tso_fraction", "five_prime_complete_fraction"
        ),
        "polya_tail_fraction": ("polya_tail_fraction", "polya_fraction"),
        "internal_priming_fraction": (
            "internal_priming_fraction", "intrapriming_fraction", "oligodt_mispriming_fraction"
        ),
        "template_switch_fraction": (
            "template_switch_fraction", "rts_fraction", "strand_invasion_fraction"
        ),
        "chimera_fraction": ("chimera_fraction", "chimeric_fraction"),
    }
    for target, names in aliases.items():
        out[target] = _fraction(_numeric(raw, names))
    mapq_fraction = _fraction(_numeric(raw, ("mapq_fraction",)))
    mapq = _numeric(raw, ("median_mapq", "mean_mapq", "mapq"))
    out["mapq_fraction"] = mapq_fraction.fillna((mapq / 60.0).clip(0.0, 1.0))
    replicate_fraction = _fraction(_numeric(raw, ("replicate_support_fraction",)))
    n_samples = _numeric(raw, ("n_samples", "n_replicates", "replicate_count"))
    out["replicate_support_fraction"] = replicate_fraction.fillna(
        (n_samples / 2.0).clip(0.0, 1.0)
    )
    return score_evidence(out)


def score_evidence(features: pd.DataFrame) -> pd.DataFrame:
    """Combine available evidence into a transparent, uncalibrated 0-1 score."""
    rows: list[dict] = []
    for _, record in features.iterrows():
        contributions: list[tuple[float, float]] = []
        warnings: list[str] = []
        row = {column: record.get(column, np.nan) for column in OUTPUT_COLUMNS}
        row["transcript_id"] = str(record["transcript_id"])
        for feature, weight in POSITIVE_WEIGHTS.items():
            value = record.get(feature, np.nan)
            if pd.notna(value):
                contributions.append((weight, float(np.clip(value, 0, 1))))
        for feature, weight in NEGATIVE_WEIGHTS.items():
            value = record.get(feature, np.nan)
            if pd.notna(value):
                contributions.append((weight, 1.0 - float(np.clip(value, 0, 1))))
        score = (
            sum(weight * value for weight, value in contributions)
            / sum(weight for weight, _ in contributions)
            if contributions else float("nan")
        )
        n_features = len(contributions)
        if n_features < 3:
            tier = "insufficient_evidence"
            warnings.append("fewer_than_three_evidence_features")
        elif (
            float(record.get("internal_priming_fraction", 0) or 0) >= 0.5
            or float(record.get("chimera_fraction", 0) or 0) >= 0.25
            or float(record.get("template_switch_fraction", 0) or 0) >= 0.4
        ):
            tier = "artifact_likely"
            warnings.append("strong_artifact_signal")
        elif score >= 0.8:
            tier = "strong"
        elif score >= 0.65:
            tier = "supported"
        elif score >= 0.45:
            tier = "limited"
        else:
            tier = "artifact_likely"
        row.update(
            {
                "isoform_evidence_score": score,
                "isoform_evidence_tier": tier,
                "evidence_feature_count": n_features,
                "evidence_model": "transparent_uncalibrated_v1",
                "evidence_warnings": json.dumps(warnings),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
