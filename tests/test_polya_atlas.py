"""Tests for craft.core.polya_atlas."""

import gzip
from pathlib import Path

import pyranges as pr

from craft.core.polya_atlas import load_atlas, match_iso_end

_BED = "\n".join(
    [
        "# header comment",
        "track name=PAS",
        "browser hide all",
        "chr1\t100\t110\tpas_a\t100\t+",
        "chr1\t500\t510\tpas_b\t80\t+",
        "chr1\t800\t810\tpas_c\t50\t-",
        "chr2\t200\t210\tpas_d\t90\t+",
        "",
        "chr_unstranded\t300\t310\tpas_skip\t10\t.",
    ]
)


def test_load_atlas_plain(tmp_path: Path) -> None:
    p = tmp_path / "atlas.bed"
    p.write_text(_BED)
    atlas = load_atlas(p)
    assert isinstance(atlas, pr.PyRanges)
    df = atlas.df
    # 5 valid rows; 4 with +/-, 1 with "." -> skipped
    assert len(df) == 4
    assert set(df["Name"]) == {"pas_a", "pas_b", "pas_c", "pas_d"}
    assert set(df["Strand"].unique()) <= {"+", "-"}


def test_load_atlas_gzipped(tmp_path: Path) -> None:
    p = tmp_path / "atlas.bed.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(_BED)
    atlas = load_atlas(p)
    assert len(atlas) == 4


def test_load_atlas_handles_empty_or_invalid(tmp_path: Path) -> None:
    p = tmp_path / "empty.bed"
    p.write_text("# only header\ntrack thing\n")
    atlas = load_atlas(p)
    assert len(atlas) == 0
    # match_iso_end on an empty atlas should produce a clean no-match.
    assert match_iso_end(100, "chr1", "+", atlas)["matched"] is False


def test_match_iso_end_plus_strand_hit(tmp_path: Path) -> None:
    p = tmp_path / "atlas.bed"
    p.write_text(_BED)
    atlas = load_atlas(p)
    # PAS at chr1:100-110, midpoint=105. iso 3' end at 100 -> distance 5, within tol.
    result = match_iso_end(100, "chr1", "+", atlas, tolerance=24)
    assert result["matched"] is True
    assert result["pas_id"] == "pas_a"
    assert result["distance_nt"] == 5  # midpoint - iso = 105 - 100


def test_match_iso_end_minus_strand_hit(tmp_path: Path) -> None:
    p = tmp_path / "atlas.bed"
    p.write_text(_BED)
    atlas = load_atlas(p)
    # PAS at chr1:800-810 on -, midpoint=805. iso 3' end at 810.
    result = match_iso_end(810, "chr1", "-", atlas, tolerance=24)
    assert result["matched"] is True
    assert result["pas_id"] == "pas_c"


def test_match_iso_end_outside_tolerance(tmp_path: Path) -> None:
    p = tmp_path / "atlas.bed"
    p.write_text(_BED)
    atlas = load_atlas(p)
    # PAS midpoints at 105, 505, 805, 205. iso at 200 on + closest is 205 (dist 5)
    # but on chr2 not chr1.
    result = match_iso_end(200, "chr1", "+", atlas, tolerance=24)
    assert result["matched"] is False
    assert result["pas_id"] == ""
    assert result["distance_nt"] == -1


def test_match_iso_end_wrong_chromosome(tmp_path: Path) -> None:
    p = tmp_path / "atlas.bed"
    p.write_text(_BED)
    atlas = load_atlas(p)
    result = match_iso_end(105, "chr_other", "+", atlas, tolerance=24)
    assert result["matched"] is False


def test_match_iso_end_wrong_strand(tmp_path: Path) -> None:
    p = tmp_path / "atlas.bed"
    p.write_text(_BED)
    atlas = load_atlas(p)
    # pas_a is on +; querying - at the same position should miss.
    result = match_iso_end(105, "chr1", "-", atlas, tolerance=24)
    assert result["matched"] is False


def test_match_iso_end_empty_atlas() -> None:
    empty_atlas = pr.PyRanges()
    result = match_iso_end(100, "chr1", "+", empty_atlas)
    assert result["matched"] is False
