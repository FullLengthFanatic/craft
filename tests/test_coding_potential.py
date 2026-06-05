"""Tests for craft.core.coding_potential."""

from pathlib import Path

import numpy as np
import pandas as pd
import pyranges as pr
import pysam

from craft.core.coding_potential import (
    _auc,
    _hexamer_score,
    _log_ratio_table,
    build_model,
    score_isoforms,
)

# Two distinct hexamer regimes: coding transcripts are GCT-rich (Ala repeats),
# non-coding are AT-rich. The model should separate them.
_CODING_UNIT = "GCTGCAGCCGCG"   # all Ala codons, varied 3rd base
_NONCODING_UNIT = "ATATATATTAAT"


def _coding_seq(n: int) -> str:
    return (_CODING_UNIT * n)[: n * 3]


def _noncoding_seq(n: int) -> str:
    return (_NONCODING_UNIT * n)[: n * 3]


def _write_genome(tmp_path: Path, seqs: dict[str, str]) -> Path:
    fa = tmp_path / "g.fa"
    fa.write_text("".join(f">{c}\n{s}\n" for c, s in seqs.items()))
    pysam.faidx(str(fa))
    return fa


def _reference(coding: dict[str, tuple], noncoding: dict[str, tuple]) -> pr.PyRanges:
    """coding/noncoding: {tx_id: (chrom, start, end, strand)} single-exon transcripts."""
    cols = ["Chromosome", "Start", "End", "Strand", "transcript_id", "Feature"]
    rows = []
    for tx, (c, s, e, st) in coding.items():
        rows.append((c, s, e, st, tx, "exon"))
        rows.append((c, s, e, st, tx, "CDS"))
    for tx, (c, s, e, st) in noncoding.items():
        rows.append((c, s, e, st, tx, "exon"))
    return pr.PyRanges(pd.DataFrame(rows, columns=cols))


def test_auc_perfect_and_random() -> None:
    assert _auc(np.array([0.9, 0.8, 0.2, 0.1]), np.array([1, 1, 0, 0])) == 1.0
    assert _auc(np.array([0.1, 0.9, 0.2, 0.8]), np.array([1, 0, 1, 0])) == 0.0
    assert _auc(np.array([0.5, 0.5]), np.array([1, 0])) == 0.5


def test_hexamer_score_separates_regimes() -> None:
    coding = _log_ratio_table(
        _hex_counts(_coding_seq(40)), _hex_counts(_noncoding_seq(40))
    )
    s_cod = _hexamer_score(_coding_seq(20), coding)
    s_non = _hexamer_score(_noncoding_seq(20), coding)
    assert s_cod > s_non


def _hex_counts(seq: str) -> np.ndarray:
    from craft.core.coding_potential import _hexamer_counts

    return _hexamer_counts(seq)


def test_build_model_returns_none_without_noncoding(tmp_path: Path) -> None:
    chrom = _coding_seq(100)
    genome = _write_genome(tmp_path, {"chr1": chrom})
    ref = _reference(coding={"c1": ("chr1", 0, 120, "+")}, noncoding={})
    with pysam.FastaFile(str(genome)) as g:
        assert build_model(ref, g) is None


def test_score_isoforms_coding_beats_noncoding(tmp_path: Path) -> None:
    # Genome: chr1 coding-like (300 nt), chr2 noncoding-like (300 nt).
    genome = _write_genome(
        tmp_path, {"chr1": _coding_seq(100), "chr2": _noncoding_seq(100)}
    )
    # Reference: several coding transcripts on chr1, several non-coding on chr2.
    coding = {f"c{i}": ("chr1", 0, 300, "+") for i in range(6)}
    noncoding = {f"n{i}": ("chr2", 0, 300, "+") for i in range(6)}
    ref = _reference(coding, noncoding)

    classified = pr.PyRanges(
        pd.DataFrame(
            [
                ("chr1", 0, 300, "+", "iso_cod"),
                ("chr2", 0, 300, "+", "iso_non"),
            ],
            columns=["Chromosome", "Start", "End", "Strand", "transcript_id"],
        )
    )
    per_iso = pd.DataFrame(
        {
            "transcript_id": ["iso_cod", "iso_non"],
            "resolved_cds_intervals": [[("chr1", 0, 300, "+")], []],
            "propagated_cds_intervals": [[], []],
            "denovo_cds_intervals": [[], [("chr2", 0, 300, "+")]],
        }
    )
    scored, model = score_isoforms(per_iso, classified, ref, genome)
    assert model is not None
    cod = scored.set_index("transcript_id").loc["iso_cod"]
    non = scored.set_index("transcript_id").loc["iso_non"]
    assert cod["coding_potential_score"] > non["coding_potential_score"]
    assert cod["coding_potential_orf_source"] == "resolved"
    assert non["coding_potential_orf_source"] == "denovo"


def test_score_isoforms_no_orf_is_noncoding(tmp_path: Path) -> None:
    genome = _write_genome(
        tmp_path, {"chr1": _coding_seq(100), "chr2": _noncoding_seq(100)}
    )
    ref = _reference(
        {f"c{i}": ("chr1", 0, 300, "+") for i in range(4)},
        {f"n{i}": ("chr2", 0, 300, "+") for i in range(4)},
    )
    classified = pr.PyRanges(
        pd.DataFrame(
            [("chr1", 0, 300, "+", "iso_x")],
            columns=["Chromosome", "Start", "End", "Strand", "transcript_id"],
        )
    )
    per_iso = pd.DataFrame(
        {
            "transcript_id": ["iso_x"],
            "resolved_cds_intervals": [[]],
            "propagated_cds_intervals": [[]],
            "denovo_cds_intervals": [[]],
        }
    )
    scored, _ = score_isoforms(per_iso, classified, ref, genome)
    row = scored.iloc[0]
    assert row["coding_potential_orf_source"] == "none"
    assert row["coding_potential_label"] == "noncoding"
    assert np.isnan(row["coding_potential_score"])
