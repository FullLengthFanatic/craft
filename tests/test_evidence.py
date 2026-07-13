"""Tests for transparent molecule/read evidence aggregation."""

from pathlib import Path

import pandas as pd

from craft.core.evidence import load_evidence, score_evidence


def test_strong_and_artifact_tiers() -> None:
    features = pd.DataFrame(
        [
            {
                "transcript_id": "good",
                "unique_molecule_fraction": 0.95,
                "canonical_junction_fraction": 1.0,
                "polya_tail_fraction": 0.9,
                "internal_priming_fraction": 0.01,
                "chimera_fraction": 0.0,
            },
            {
                "transcript_id": "artifact",
                "unique_molecule_fraction": 0.2,
                "canonical_junction_fraction": 0.5,
                "polya_tail_fraction": 0.1,
                "internal_priming_fraction": 0.8,
                "chimera_fraction": 0.4,
            },
        ]
    )
    out = score_evidence(features).set_index("transcript_id")
    assert out.loc["good", "isoform_evidence_tier"] == "strong"
    assert out.loc["artifact", "isoform_evidence_tier"] == "artifact_likely"
    assert out.loc["good", "isoform_evidence_score"] > out.loc["artifact", "isoform_evidence_score"]
    assert "uncalibrated" in out.loc["good", "evidence_model"]


def test_load_evidence_accepts_percentages_and_derives_unique_fraction(tmp_path: Path) -> None:
    path = tmp_path / "evidence.tsv"
    path.write_text(
        "transcript_id\ttotal_molecules\tunique_molecules\tcanonical_junction_fraction"
        "\tinternal_priming_fraction\n"
        "tx1\t10\t8\t100\t5\n"
    )
    row = load_evidence(path).iloc[0]
    assert row["unique_molecule_fraction"] == 0.8
    assert row["canonical_junction_fraction"] == 1.0
    assert row["internal_priming_fraction"] == 0.05


def test_sparse_feature_set_is_not_overinterpreted() -> None:
    out = score_evidence(
        pd.DataFrame([{"transcript_id": "tx", "polya_tail_fraction": 1.0}])
    ).iloc[0]
    assert out["isoform_evidence_tier"] == "insufficient_evidence"


def test_canonical_mapq_and_replicate_fractions_are_accepted(tmp_path: Path) -> None:
    path = tmp_path / "evidence.tsv"
    path.write_text(
        "transcript_id\tmapq_fraction\treplicate_support_fraction\tpolya_tail_fraction\n"
        "tx1\t0.75\t50\t1\n"
    )
    row = load_evidence(path).iloc[0]
    assert row["mapq_fraction"] == 0.75
    assert row["replicate_support_fraction"] == 0.5
