"""Tests for craft.core.utr3 (resolved UTR metrics + poly(A) motif scan)."""

from pathlib import Path

import pandas as pd
import pyranges as pr
import pysam

from craft.core.utr3 import _utr5_length, annotate, polya_near_3prime_end, polya_signal


def _write_genome(tmp_path: Path, length: int, fills: list[tuple[int, str]]) -> Path:
    ba = bytearray(b"C" * length)
    for start, seq in fills:
        ba[start : start + len(seq)] = seq.encode()
    fa = tmp_path / "g.fa"
    fa.write_text(f">chr1\n{ba.decode()}\n")
    pysam.faidx(str(fa))
    return fa


def _classified(iso_records: list[tuple], parent_tx: str, tx_id: str = "t1") -> pr.PyRanges:
    cols = ["Chromosome", "Start", "End", "Strand", "transcript_id", "parent_tx_id"]
    rows = [(c, s, e, st, tx_id, parent_tx) for (c, s, e, st) in iso_records]
    return pr.PyRanges(pd.DataFrame(rows, columns=cols))


def _reference(exons: list[tuple], cds: list[tuple]) -> pr.PyRanges:
    cols = ["Chromosome", "Start", "End", "Strand", "transcript_id", "Feature"]
    rows = [(*r, "exon") for r in exons] + [(*r, "CDS") for r in cds]
    return pr.PyRanges(pd.DataFrame(rows, columns=cols))


def _res(tx_id, status, intervals, stop=True):
    return pd.DataFrame(
        [
            {
                "transcript_id": tx_id,
                "resolved_orf_status": status,
                "resolved_stop_pos": None,
                "resolved_cds_bp": 100,
                "resolved_aa_length": 33,
                "resolved_cds_intervals": intervals,
                "ptc_introduced": False,
                "intron_retained_in_cds": False,
                "frame_consistent": True,
                "stop_in_transcript": stop,
                "uorf_count": 0,
                "uorf_triggers_nmd": False,
            }
        ]
    )


# ---- pure functions ----------------------------------------------------------

def test_polya_signal_finds_canonical_motif() -> None:
    sig = polya_signal("GGGGAATAAACCC")
    assert sig["motif"] == "AATAAA"
    assert sig["distance_from_3p_end"] == 3


def test_polya_signal_none() -> None:
    assert polya_signal("GGGGCCCCGGGG")["motif"] == ""


def test_utr5_length_plus_and_minus() -> None:
    exons = pd.DataFrame({"Start": [100, 300], "End": [200, 400]})
    assert _utr5_length(exons, 150, "+") == 50
    assert _utr5_length(exons, 350, "-") == 49


# ---- annotate ----------------------------------------------------------------

def test_annotate_full_length_zero_deltas_and_motif(tmp_path: Path) -> None:
    iso = [("chr1", 100, 200, "+"), ("chr1", 300, 400, "+")]
    classified = _classified(iso, "t_ref")
    ref = _reference(
        exons=[("chr1", 100, 200, "+", "t_ref"), ("chr1", 300, 400, "+", "t_ref")],
        cds=[("chr1", 150, 200, "+", "t_ref"), ("chr1", 300, 350, "+", "t_ref")],
    )
    res = _res("t1", "intact", [("chr1", 150, 200, "+"), ("chr1", 300, 350, "+")])
    genome = _write_genome(tmp_path, 400, [(360, "AATAAA")])  # motif in the 3'UTR
    out = annotate(classified, res, ref, genome_fasta=genome)
    r = out[out["transcript_id"] == "t1"].iloc[0]
    assert r["iso_utr3_length_nt"] == 47  # stop [350,353), UTR [353,400)
    assert r["parent_utr3_length_nt"] == 47
    assert r["utr3_length_delta_nt"] == 0
    assert r["iso_utr5_length_nt"] == 50  # exon1 [100,150)
    assert r["utr5_length_delta_nt"] == 0
    assert r["polya_signal_motif"] == "AATAAA"
    assert not r["long_utr3_triggers_nmd"]


def test_annotate_long_utr3_flag(tmp_path: Path) -> None:
    iso = [("chr1", 100, 200, "+"), ("chr1", 300, 2000, "+")]
    classified = _classified(iso, "t_ref")
    ref = _reference(
        exons=[("chr1", 100, 200, "+", "t_ref"), ("chr1", 300, 2000, "+", "t_ref")],
        cds=[("chr1", 150, 200, "+", "t_ref"), ("chr1", 300, 350, "+", "t_ref")],
    )
    res = _res("t1", "intact", [("chr1", 150, 200, "+"), ("chr1", 300, 350, "+")])
    out = annotate(classified, res, ref, genome_fasta=None, long_utr3_nt=1000)
    r = out[out["transcript_id"] == "t1"].iloc[0]
    assert r["iso_utr3_length_nt"] == 1647
    assert r["long_utr3_triggers_nmd"]
    assert r["polya_signal_motif"] == ""  # no genome -> no scan


def test_utr3_excludes_split_stop_codon_across_junction(tmp_path: Path) -> None:
    iso = [("chr1", 100, 105, "+"), ("chr1", 300, 310, "+")]
    classified = _classified(iso, "t_ref")
    ref = _reference(
        exons=[("chr1", 100, 105, "+", "t_ref"), ("chr1", 300, 310, "+", "t_ref")],
        cds=[("chr1", 100, 104, "+", "t_ref")],
    )
    res = _res("t1", "intact", [("chr1", 100, 104, "+")])
    res["resolved_stop_codon_pos"] = 103
    out = annotate(classified, res, ref, genome_fasta=None)
    row = out.iloc[0]
    # Stop bases are genomic 103, 104 and 300; UTR is exon-2 positions 301..309.
    assert row["iso_utr3_length_nt"] == 9


def test_annotate_resolution_failed_yields_nulls(tmp_path: Path) -> None:
    iso = [("chr1", 100, 200, "+"), ("chr1", 300, 400, "+")]
    classified = _classified(iso, "t_ref")
    ref = _reference(
        exons=[("chr1", 100, 200, "+", "t_ref"), ("chr1", 300, 400, "+", "t_ref")],
        cds=[("chr1", 150, 200, "+", "t_ref"), ("chr1", 300, 350, "+", "t_ref")],
    )
    res = _res("t1", "resolution_failed", [], stop=False)
    out = annotate(classified, res, ref, genome_fasta=None)
    r = out[out["transcript_id"] == "t1"].iloc[0]
    assert pd.isna(r["iso_utr3_length_nt"])
    assert pd.isna(r["iso_utr5_length_nt"])
    assert r["polya_signal_motif"] == ""


def test_annotate_parent_stop_unplaceable_skips_parent_delta(tmp_path: Path) -> None:
    # Parent CDS ends exactly at the transcript 3' end and has no stop_codon record
    # (e.g. GENCODE cds_end_NF): the stop codon cannot be placed on its exon chain.
    # This must not crash the run; the parent-relative 3'UTR fields are left empty
    # while the isoform's own 3'UTR is still reported.
    iso = [("chr1", 100, 200, "+"), ("chr1", 300, 400, "+")]
    classified = _classified(iso, "t_ref")
    ref = _reference(
        exons=[("chr1", 100, 200, "+", "t_ref")],
        cds=[("chr1", 150, 200, "+", "t_ref")],
    )
    res = _res("t1", "intact", [("chr1", 150, 200, "+"), ("chr1", 300, 350, "+")])
    res["resolved_stop_codon_pos"] = 350
    out = annotate(classified, res, ref, genome_fasta=None)
    r = out[out["transcript_id"] == "t1"].iloc[0]
    assert r["iso_utr3_length_nt"] == 47
    assert pd.isna(r["parent_utr3_length_nt"])
    assert pd.isna(r["utr3_length_delta_nt"])


def test_polya_near_3prime_end_plus_strand(tmp_path: Path) -> None:
    genome = _write_genome(tmp_path, 200, [(140, "AATAAA")])
    exons = pd.DataFrame(
        [{"Chromosome": "chr1", "Start": 100, "End": 160, "Strand": "+", "transcript_id": "t1"}]
    )
    res = polya_near_3prime_end(exons, "+", genome)
    assert res["found"] is True
    assert res["motif"] == "AATAAA"
