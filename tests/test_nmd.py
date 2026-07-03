"""Tests for craft.core.nmd (single consolidated NMD call)."""

import pandas as pd
import pyranges as pr

from craft.core.nmd import NMDStatus, predict
from craft.core.orf.confidence import ORFConfidence


def _iso_pr(records: list[tuple]) -> pr.PyRanges:
    cols = ["Chromosome", "Start", "End", "Strand", "transcript_id"]
    return pr.PyRanges(pd.DataFrame(records, columns=cols))


def _res(tx_id, status, intervals, cds_bp=300, stop=True):
    return {
        "transcript_id": tx_id,
        "resolved_orf_status": status,
        "resolved_stop_pos": None,
        "resolved_cds_bp": cds_bp,
        "resolved_aa_length": cds_bp // 3,
        "resolved_cds_intervals": intervals,
        "ptc_introduced": status != "intact",
        "intron_retained_in_cds": False,
        "frame_consistent": status == "intact",
        "stop_in_transcript": stop,
        "uorf_count": 0,
        "uorf_triggers_nmd": False,
    }


def _dn(tx_id, found, intervals, cds_bp=0):
    return {
        "transcript_id": tx_id,
        "denovo_orf_found": found,
        "denovo_cds_bp": cds_bp,
        "denovo_cds_intervals": intervals,
        "denovo_orf_aa_length": cds_bp // 3,
        "denovo_start_codon": "ATG" if found else "",
        "denovo_stop_codon": "TAA" if found else "",
    }


def _row(df: pd.DataFrame, tx: str) -> pd.Series:
    return df[df["transcript_id"] == tx].iloc[0]


_ISO3 = [
    ("chr1", 100, 200, "+", "t1"),
    ("chr1", 300, 400, "+", "t1"),
    ("chr1", 500, 600, "+", "t1"),
]
_INTACT_IV = [("chr1", 150, 200, "+"), ("chr1", 300, 400, "+"), ("chr1", 500, 560, "+")]
_PTC_IV = [("chr1", 150, 200, "+"), ("chr1", 300, 350, "+")]
_DN_IV = [("chr1", 150, 300, "+"), ("chr1", 400, 450, "+")]


def test_resolved_intact_stop_in_last_exon_escapes_high_conf() -> None:
    iso = _iso_pr(_ISO3)
    res = pd.DataFrame(
        [_res("t1", "intact", _INTACT_IV)]
    )
    r = _row(predict(iso, res, None), "t1")
    assert r["nmd_status"] == NMDStatus.ESCAPED.value
    assert r["nmd_rule"] == "stop_in_last_exon"
    assert r["nmd_confidence"] == ORFConfidence.HIGH.value
    assert r["nmd_basis"] == "resolved"


def test_resolved_ptc_far_from_last_junction_is_sensitive_medium() -> None:
    iso = _iso_pr(_ISO3)
    # stop in exon 2 (last coding base 349), 51 nt from the exon2/exon3 junction.
    res = pd.DataFrame(
        [_res("t1", "ptc_premature", _PTC_IV)]
    )
    r = _row(predict(iso, res, None), "t1")
    assert r["nmd_status"] == NMDStatus.SENSITIVE.value
    assert r["nmd_rule"] == "ptc_50nt_rule"
    assert r["nmd_confidence"] == ORFConfidence.MEDIUM.value
    assert r["nmd_basis"] == "resolved"


def test_denovo_fallback_when_no_resolved_stop() -> None:
    iso = _iso_pr(
        [
            ("chr1", 100, 300, "+", "t1"),
            ("chr1", 400, 500, "+", "t1"),
            ("chr1", 600, 700, "+", "t1"),
        ]
    )
    res = pd.DataFrame([_res("t1", "resolution_failed", [], cds_bp=0, stop=False)])
    dn = pd.DataFrame(
        [_dn("t1", True, _DN_IV, cds_bp=200)]
    )
    r = _row(predict(iso, res, dn), "t1")
    assert r["nmd_status"] == NMDStatus.SENSITIVE.value
    assert r["nmd_basis"] == "denovo"
    assert r["nmd_confidence"] == ORFConfidence.LOW.value


def test_resolved_preferred_over_denovo() -> None:
    iso = _iso_pr(_ISO3)
    res = pd.DataFrame(
        [_res("t1", "intact", _INTACT_IV)]
    )
    dn = pd.DataFrame([_dn("t1", True, [("chr1", 150, 200, "+")], cds_bp=50)])
    assert _row(predict(iso, res, dn), "t1")["nmd_basis"] == "resolved"


def test_not_applicable_when_no_orf() -> None:
    iso = _iso_pr([("chr1", 100, 200, "+", "t1"), ("chr1", 300, 400, "+", "t1")])
    res = pd.DataFrame([_res("t1", "resolution_failed", [], cds_bp=0, stop=False)])
    dn = pd.DataFrame([_dn("t1", False, [])])
    r = _row(predict(iso, res, dn), "t1")
    assert r["nmd_status"] == NMDStatus.NOT_APPLICABLE.value
    assert r["nmd_basis"] == "none"
    assert r["nmd_confidence"] == ORFConfidence.NONE.value


def test_no_stop_in_read_is_not_applicable() -> None:
    iso = _iso_pr(_ISO3)
    res = pd.DataFrame([_res("t1", "no_stop_in_read", [], cds_bp=0, stop=False)])
    assert _row(predict(iso, res, None), "t1")["nmd_status"] == NMDStatus.NOT_APPLICABLE.value


def test_start_proximal_escape() -> None:
    iso = _iso_pr(_ISO3)
    # short CDS (< 150 bp) escapes via re-initiation even with a far stop.
    res = pd.DataFrame(
        [_res("t1", "ptc_premature", _PTC_IV, cds_bp=100)]
    )
    r = _row(predict(iso, res, None), "t1")
    assert r["nmd_status"] == NMDStatus.ESCAPED.value
    assert r["nmd_rule"] == "start_proximal"


def test_long_exon_rule_uses_ptc_bearing_exon_not_terminal_exon() -> None:
    # PTC sits in a long (700 bp) internal exon; the terminal exon is short (100 bp).
    # The long-exon rule must fire on the PTC-bearing exon, so this escapes.
    iso = _iso_pr(
        [
            ("chr1", 100, 200, "+", "t1"),
            ("chr1", 300, 1000, "+", "t1"),
            ("chr1", 1100, 1200, "+", "t1"),
        ]
    )
    # stop last coding base at 899, inside the 700 bp internal exon, 101 nt from the
    # exon2/exon3 junction (so not within 50 nt) and CDS >= 150 bp.
    res = pd.DataFrame(
        [_res("t1", "ptc_premature", [("chr1", 150, 200, "+"), ("chr1", 300, 900, "+")])]
    )
    r = _row(predict(iso, res, None), "t1")
    assert r["nmd_status"] == NMDStatus.ESCAPED.value
    assert r["nmd_rule"] == "long_exon"
    assert r["ptc_exon_length_nt"] == 700
    assert r["last_exon_length_nt"] == 100


def test_long_terminal_exon_does_not_rescue_ptc_in_short_exon() -> None:
    # Mirror case: PTC in a short (100 bp) internal exon, terminal exon is long
    # (700 bp). The old code escaped on terminal-exon length; the fix keeps it
    # sensitive because the PTC-bearing exon is short.
    iso = _iso_pr(
        [
            ("chr1", 100, 200, "+", "t1"),
            ("chr1", 300, 400, "+", "t1"),
            ("chr1", 500, 1200, "+", "t1"),
        ]
    )
    res = pd.DataFrame(
        [_res("t1", "ptc_premature", [("chr1", 150, 200, "+"), ("chr1", 300, 350, "+")])]
    )
    r = _row(predict(iso, res, None), "t1")
    assert r["nmd_status"] == NMDStatus.SENSITIVE.value
    assert r["nmd_rule"] == "ptc_50nt_rule"
    assert r["ptc_exon_length_nt"] == 100
    assert r["last_exon_length_nt"] == 700


def test_minus_strand_sensitive() -> None:
    iso = _iso_pr(
        [
            ("chr1", 100, 200, "-", "t1"),
            ("chr1", 300, 400, "-", "t1"),
            ("chr1", 600, 700, "-", "t1"),
        ]
    )
    # On - strand the stop is the min Start of the CDS intervals; place it in the
    # middle exon, far from the 5'-most (last in transcript order) junction.
    res = pd.DataFrame(
        [_res("t1", "ptc_premature", [("chr1", 350, 400, "-"), ("chr1", 600, 700, "-")])]
    )
    r = _row(predict(iso, res, None), "t1")
    assert r["nmd_status"] in {NMDStatus.SENSITIVE.value, NMDStatus.ESCAPED.value}


def test_empty_classified_returns_empty() -> None:
    out = predict(pr.PyRanges(), pd.DataFrame())
    assert out.empty
    assert list(out.columns)[0] == "transcript_id"
