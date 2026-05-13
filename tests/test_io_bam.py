"""Tests for craft.io.bam."""

from pathlib import Path

import pysam
import pytest

from craft.io.bam import open_alignments


def _write_sorted_bam(path: Path) -> None:
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"LN": 1000, "SN": "chr1"}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for i, start in enumerate([100, 200, 300]):
            a = pysam.AlignedSegment(bam.header)
            a.query_name = f"read{i}"
            a.query_sequence = "ACGT"
            a.flag = 0
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigartuples = [(0, 4)]
            a.query_qualities = pysam.qualitystring_to_array("IIII")
            bam.write(a)


def test_open_alignments_returns_pysam_handle(tmp_path: Path) -> None:
    bam = tmp_path / "test.bam"
    _write_sorted_bam(bam)
    pysam.index(str(bam))
    handle = open_alignments(bam)
    try:
        assert isinstance(handle, pysam.AlignmentFile)
    finally:
        handle.close()


def test_open_alignments_builds_index_when_missing(tmp_path: Path) -> None:
    bam = tmp_path / "test.bam"
    _write_sorted_bam(bam)
    bai = Path(str(bam) + ".bai")
    assert not bai.exists()
    handle = open_alignments(bam)
    try:
        assert bai.exists(), ".bai should be created if missing"
    finally:
        handle.close()


def test_open_alignments_can_fetch_reads(tmp_path: Path) -> None:
    bam = tmp_path / "test.bam"
    _write_sorted_bam(bam)
    handle = open_alignments(bam)
    try:
        names = sorted(r.query_name for r in handle.fetch("chr1"))
        assert names == ["read0", "read1", "read2"]
    finally:
        handle.close()


def test_open_alignments_reuses_existing_index(tmp_path: Path) -> None:
    bam = tmp_path / "test.bam"
    _write_sorted_bam(bam)
    pysam.index(str(bam))
    bai = Path(str(bam) + ".bai")
    mtime_before = bai.stat().st_mtime
    handle = open_alignments(bam)
    try:
        assert bai.stat().st_mtime == mtime_before, "Should not rebuild existing .bai"
    finally:
        handle.close()


def test_open_alignments_raises_for_missing_file(tmp_path: Path) -> None:
    missing = tmp_path / "nope.bam"
    with pytest.raises((FileNotFoundError, OSError, ValueError)):
        open_alignments(missing)
