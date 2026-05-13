"""Tests for craft.io.fasta."""

from pathlib import Path

import pysam

from craft.io.fasta import open_genome


def _write_fasta(path: Path, contigs: dict[str, str]) -> None:
    text = "".join(f">{name}\n{seq}\n" for name, seq in contigs.items())
    path.write_text(text)


def test_open_genome_returns_pysam_handle(tmp_path: Path) -> None:
    fasta = tmp_path / "g.fa"
    _write_fasta(fasta, {"chr1": "ACGTACGTACGT"})
    pysam.faidx(str(fasta))
    handle = open_genome(fasta)
    try:
        assert isinstance(handle, pysam.FastaFile)
    finally:
        handle.close()


def test_open_genome_builds_index_when_missing(tmp_path: Path) -> None:
    fasta = tmp_path / "g.fa"
    _write_fasta(fasta, {"chr1": "ACGTACGTACGT"})
    fai = Path(str(fasta) + ".fai")
    assert not fai.exists()
    handle = open_genome(fasta)
    try:
        assert fai.exists(), ".fai should be created if missing"
    finally:
        handle.close()


def test_open_genome_can_fetch_sequence(tmp_path: Path) -> None:
    fasta = tmp_path / "g.fa"
    _write_fasta(fasta, {"chr1": "ACGTACGTACGT"})
    handle = open_genome(fasta)
    try:
        assert handle.fetch("chr1", 0, 4) == "ACGT"
        assert handle.fetch("chr1", 4, 8) == "ACGT"
    finally:
        handle.close()


def test_open_genome_reuses_existing_index(tmp_path: Path) -> None:
    fasta = tmp_path / "g.fa"
    _write_fasta(fasta, {"chr1": "ACGTACGTACGT"})
    pysam.faidx(str(fasta))
    fai = Path(str(fasta) + ".fai")
    mtime_before = fai.stat().st_mtime
    handle = open_genome(fasta)
    try:
        assert fai.stat().st_mtime == mtime_before, "Should not rebuild existing .fai"
    finally:
        handle.close()
