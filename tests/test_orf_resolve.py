"""Tests for craft.core.orf.resolve (sequence-level ORF resolution)."""

from pathlib import Path

import pandas as pd
import pyranges as pr
import pysam

from craft.core.completeness import classify
from craft.core.orf.denovo import (
    _transcript_to_genomic_intervals,
)
from craft.core.orf.propagation import ORFOutcome, propagate
from craft.core.orf.resolve import ResolvedORFStatus, resolve

_RC = str.maketrans("ACGTN", "TGCAN")


def _rc(seq: str) -> str:
    return seq.translate(_RC)[::-1]


def _layout(exons: list[tuple[int, int]], strand: str, tx: str) -> list[tuple[int, str]]:
    """Genome fills so that the spliced transcript over ``exons`` equals ``tx``.

    For ``-`` strand the forward-genome sequence is the reverse complement of the
    transcript, split across exons in genomic order.
    """
    fwd = tx if strand == "+" else _rc(tx)
    assert sum(e - s for s, e in exons) == len(fwd), (sum(e - s for s, e in exons), len(fwd))
    fills, pos = [], 0
    for s, e in sorted(exons):
        length = e - s
        fills.append((s, fwd[pos : pos + length]))
        pos += length
    return fills


def _write_genome(tmp_path: Path, length: int, fills: list[tuple[int, str]]) -> Path:
    ba = bytearray(b"A" * length)
    for start, seq in fills:
        ba[start : start + len(seq)] = seq.encode()
    fa = tmp_path / "genome.fa"
    fa.write_text(f">chr1\n{ba.decode()}\n")
    pysam.faidx(str(fa))
    return fa


def _exons(records: list[tuple], strand: str, tx: str) -> pr.PyRanges:
    cols = ["Chromosome", "Start", "End", "Strand", "transcript_id"]
    rows = [("chr1", s, e, strand, tx) for s, e in records]
    return pr.PyRanges(pd.DataFrame(rows, columns=cols))


def _reference(
    exons: list[tuple[int, int]], cds_intervals: list[tuple], strand: str, tx: str
) -> pr.PyRanges:
    cols = ["Chromosome", "Start", "End", "Strand", "transcript_id", "Feature"]
    rows = [("chr1", s, e, strand, tx, "exon") for s, e in exons]
    rows += [(c, s, e, st, tx, "CDS") for c, s, e, st in cds_intervals]
    return pr.PyRanges(pd.DataFrame(rows, columns=cols))


def _ref_exons_only(ref: pr.PyRanges) -> pr.PyRanges:
    df = ref.df[ref.df["Feature"] == "exon"].drop(columns="Feature").copy()
    df["Strand"] = df["Strand"].astype(str)
    return pr.PyRanges(df)


def _cds_for(exons: list[tuple[int, int]], strand: str, cstart: int, cend: int) -> list[tuple]:
    ref_df = pd.DataFrame(
        [("chr1", s, e, strand, "ref") for s, e in exons],
        columns=["Chromosome", "Start", "End", "Strand", "transcript_id"],
    )
    return _transcript_to_genomic_intervals(cstart, cend, ref_df, strand, "chr1")


def _run(iso: pr.PyRanges, ref: pr.PyRanges, genome: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    classified = classify(iso, _ref_exons_only(ref))
    propagated = propagate(classified, ref)
    resolved = resolve(classified, propagated, ref, genome)
    return propagated, resolved


def _row(df: pd.DataFrame, tx: str) -> pd.Series:
    return df[df["transcript_id"] == tx].iloc[0]


# transcript: 12nt 5'UTR, 48nt coding (16 Ala), stop, 3'UTR. Coding tx is [12, 60).
_TX_INTACT = "C" * 12 + "GCT" * 16 + "TAA" + "A" * 7  # len 70


def test_intact_plus_strand(tmp_path: Path) -> None:
    exons = [(100, 130), (200, 240)]
    cds = _cds_for(exons, "+", 12, 60)
    ref = _reference(exons, cds, "+", "ref")
    iso = _exons(exons, "+", "t1")
    genome = _write_genome(tmp_path, 240, _layout(exons, "+", _TX_INTACT))
    _, resolved = _run(iso, ref, genome)
    r = _row(resolved, "t1")
    assert r["resolved_orf_status"] == ResolvedORFStatus.INTACT.value
    assert r["stop_in_transcript"]
    assert not r["ptc_introduced"]
    assert r["frame_consistent"]
    assert not r["intron_retained_in_cds"]
    assert r["resolved_aa_length"] == 16
    assert r["resolved_stop_pos"] == 229


def test_intact_minus_strand(tmp_path: Path) -> None:
    exons = [(100, 130), (200, 240)]
    cds = _cds_for(exons, "-", 12, 60)
    ref = _reference(exons, cds, "-", "ref")
    iso = _exons(exons, "-", "t1")
    genome = _write_genome(tmp_path, 240, _layout(exons, "-", _TX_INTACT))
    _, resolved = _run(iso, ref, genome)
    r = _row(resolved, "t1")
    assert r["resolved_orf_status"] == ResolvedORFStatus.INTACT.value
    assert r["resolved_aa_length"] == 16
    assert r["frame_consistent"]


def test_in_frame_ptc_caught_when_geometric_says_intact(tmp_path: Path) -> None:
    # Parent CDS is annotated to [12,60) but the sequence has a stop at tx30.
    exons = [(100, 130), (200, 240)]
    cds = _cds_for(exons, "+", 12, 60)
    ref = _reference(exons, cds, "+", "ref")
    iso = _exons(exons, "+", "t1")
    tx = "C" * 12 + "GCT" * 6 + "TAA" + "A" * 37  # stop at tx30; len 70
    genome = _write_genome(tmp_path, 240, _layout(exons, "+", tx))
    propagated, resolved = _run(iso, ref, genome)
    assert _row(propagated, "t1")["orf_outcome"] == ORFOutcome.PROPAGATED_INTACT.value
    r = _row(resolved, "t1")
    assert r["resolved_orf_status"] == ResolvedORFStatus.PTC_PREMATURE.value
    assert r["ptc_introduced"]
    assert not r["intron_retained_in_cds"]
    assert r["resolved_aa_length"] == 6


def test_intron_retention_in_cds_introduces_ptc(tmp_path: Path) -> None:
    parent_exons = [(100, 130), (200, 240)]
    cds = _cds_for(parent_exons, "+", 12, 60)
    ref = _reference(parent_exons, cds, "+", "ref")
    iso = _exons([(100, 240)], "+", "t1")  # single exon: retains intron [130,200)
    # parent exonic content + a stop in the retained intron, in frame from the start.
    fills = _layout(parent_exons, "+", _TX_INTACT) + [(130, "TAA")]
    genome = _write_genome(tmp_path, 240, fills)
    propagated, resolved = _run(iso, ref, genome)
    # geometric propagation cannot see the retained intron:
    assert _row(propagated, "t1")["orf_outcome"] == ORFOutcome.PROPAGATED_INTACT.value
    r = _row(resolved, "t1")
    assert r["intron_retained_in_cds"]
    assert r["ptc_introduced"]
    assert r["resolved_orf_status"] == ResolvedORFStatus.PTC_INTRON_RETAINED.value


def test_no_stop_in_read_when_3prime_truncated(tmp_path: Path) -> None:
    parent_exons = [(100, 130), (200, 240)]
    cds = _cds_for(parent_exons, "+", 12, 60)
    ref = _reference(parent_exons, cds, "+", "ref")
    iso = _exons([(100, 130), (200, 225)], "+", "t1")  # drops the stop region
    genome = _write_genome(tmp_path, 240, _layout(parent_exons, "+", _TX_INTACT))
    propagated, resolved = _run(iso, ref, genome)
    assert _row(propagated, "t1")["orf_outcome"] == ORFOutcome.STOP_NOT_OBSERVED.value
    r = _row(resolved, "t1")
    assert not r["stop_in_transcript"]
    assert r["resolved_orf_status"] == ResolvedORFStatus.NO_STOP_IN_READ.value
    assert r["resolved_aa_length"] == 0


def test_no_parent_yields_resolution_failed(tmp_path: Path) -> None:
    exons = [(100, 130), (200, 240)]
    cds = _cds_for(exons, "+", 12, 60)
    ref = _reference(exons, cds, "+", "ref")
    iso = pr.PyRanges(
        pd.DataFrame(
            [("chr2", 100, 130, "+", "t_novel")],
            columns=["Chromosome", "Start", "End", "Strand", "transcript_id"],
        )
    )
    genome = _write_genome(tmp_path, 240, _layout(exons, "+", _TX_INTACT))
    propagated, resolved = _run(iso, ref, genome)
    assert _row(propagated, "t_novel")["orf_outcome"] == ORFOutcome.NO_PARENT.value
    status = _row(resolved, "t_novel")["resolved_orf_status"]
    assert status == ResolvedORFStatus.RESOLUTION_FAILED.value


def test_uorf_detected_and_flags_nmd(tmp_path: Path) -> None:
    exons = [(100, 250), (400, 440)]  # long first exon -> distant last junction
    cds = _cds_for(exons, "+", 12, 60)
    ref = _reference(exons, cds, "+", "ref")
    iso = _exons(exons, "+", "t1")
    # 5'UTR carries one uORF (ATG AAA TGA), then the main CDS at tx12.
    tx = "ATGAAATGA" + "CCC" + "GCT" * 16 + "TAA" + "A" * (190 - 9 - 3 - 48 - 3)
    assert len(tx) == 190
    genome = _write_genome(tmp_path, 440, _layout(exons, "+", tx))
    _, resolved = _run(iso, ref, genome)
    r = _row(resolved, "t1")
    assert r["uorf_count"] == 1
    assert r["uorf_triggers_nmd"]
    assert r["resolved_orf_status"] == ResolvedORFStatus.INTACT.value


def test_no_uorf_when_5utr_clean(tmp_path: Path) -> None:
    exons = [(100, 130), (200, 240)]
    cds = _cds_for(exons, "+", 12, 60)
    ref = _reference(exons, cds, "+", "ref")
    iso = _exons(exons, "+", "t1")
    genome = _write_genome(tmp_path, 240, _layout(exons, "+", _TX_INTACT))
    _, resolved = _run(iso, ref, genome)
    r = _row(resolved, "t1")
    assert r["uorf_count"] == 0
    assert not r["uorf_triggers_nmd"]


def test_empty_inputs_return_empty_frame(tmp_path: Path) -> None:
    genome = _write_genome(tmp_path, 240, [])
    out = resolve(pr.PyRanges(), pd.DataFrame(), pr.PyRanges(), genome)
    assert out.empty
    assert list(out.columns)[0] == "transcript_id"
