"""Integration tests for craft.pipeline."""

import json
from pathlib import Path

import pandas as pd
import pysam
import pytest

from craft.pipeline import run_annotate

_RC_TABLE = str.maketrans("ACGTN", "TGCAN")


def _rc(s: str) -> str:
    return s.translate(_RC_TABLE)[::-1]


@pytest.fixture
def pipeline_inputs(tmp_path: Path) -> dict[str, Path]:
    # Reference GTF: one parent transcript on chr1 + strand, 3 exons + 3 CDS.
    ref_rows = [
        'chr1\tCRAFT\tgene\t101\t600\t.\t+\t.\tgene_id "g1";',
        'chr1\tCRAFT\ttranscript\t101\t600\t.\t+\t.\tgene_id "g1"; transcript_id "t_parent";',
        'chr1\tCRAFT\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t_parent";',
        'chr1\tCRAFT\texon\t301\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t_parent";',
        'chr1\tCRAFT\texon\t501\t600\t.\t+\t.\tgene_id "g1"; transcript_id "t_parent";',
        'chr1\tCRAFT\tCDS\t151\t200\t.\t+\t0\tgene_id "g1"; transcript_id "t_parent";',
        'chr1\tCRAFT\tCDS\t301\t400\t.\t+\t0\tgene_id "g1"; transcript_id "t_parent";',
        'chr1\tCRAFT\tCDS\t501\t550\t.\t+\t0\tgene_id "g1"; transcript_id "t_parent";',
    ]
    ref_path = tmp_path / "reference.gtf"
    ref_path.write_text("\n".join(ref_rows) + "\n")

    # Isoforms GTF:
    # - t_intact: identical exons to t_parent -> PROPAGATED_INTACT, FULL_LENGTH
    # - t_novel: on chr2, no parent in reference -> NOVEL_NO_MATCH (de novo path)
    iso_rows = [
        'chr1\tFLAIR\texon\t101\t200\t.\t+\t.\ttranscript_id "t_intact";',
        'chr1\tFLAIR\texon\t301\t400\t.\t+\t.\ttranscript_id "t_intact";',
        'chr1\tFLAIR\texon\t501\t600\t.\t+\t.\ttranscript_id "t_intact";',
        'chr2\tFLAIR\texon\t1\t200\t.\t+\t.\ttranscript_id "t_novel";',
    ]
    iso_path = tmp_path / "isoforms.gtf"
    iso_path.write_text("\n".join(iso_rows) + "\n")

    # FASTA:
    # chr1: 600 bp. exon1 region (100-199) needs AATAAA-free filler in the UTR portion.
    #   The iso's 3' UTR for t_intact = positions 550-599. Make it "GCGCGCGC..." (50 chars).
    # chr2: 200 bp with a clear ORF for the de novo path: ATG + 59 GCC + TAA + 17 G.
    chr1_seq = (
        ("N" * 100)
        + ("G" * 100)  # exon1
        + ("N" * 100)  # intron 1
        + ("G" * 100)  # exon2
        + ("N" * 100)  # intron 2
        + ("G" * 50)   # CDS portion of exon3 (positions 500-549)
        + ("GCGC" * 12 + "GC")  # UTR portion of exon3 (positions 550-599); no AATAAA
    )
    assert len(chr1_seq) == 600
    chr2_seq = "ATG" + ("GCC" * 59) + "TAA" + ("G" * 17)
    assert len(chr2_seq) == 200

    fasta_path = tmp_path / "genome.fa"
    fasta_path.write_text(f">chr1\n{chr1_seq}\n>chr2\n{chr2_seq}\n")
    pysam.faidx(str(fasta_path))

    return {
        "isoforms": iso_path,
        "reference": ref_path,
        "genome": fasta_path,
        "output_dir": tmp_path / "out",
    }


def test_pipeline_runs_end_to_end(pipeline_inputs: dict[str, Path]) -> None:
    result = run_annotate(
        isoforms_path=pipeline_inputs["isoforms"],
        reference_path=pipeline_inputs["reference"],
        output_dir=pipeline_inputs["output_dir"],
        genome_path=pipeline_inputs["genome"],
    )
    assert len(result) == 2
    assert {"t_intact", "t_novel"} == set(result["transcript_id"])


def test_pipeline_writes_tsv_and_json(pipeline_inputs: dict[str, Path]) -> None:
    run_annotate(
        isoforms_path=pipeline_inputs["isoforms"],
        reference_path=pipeline_inputs["reference"],
        output_dir=pipeline_inputs["output_dir"],
        genome_path=pipeline_inputs["genome"],
    )
    tsv = pipeline_inputs["output_dir"] / "per_isoform.tsv"
    js = pipeline_inputs["output_dir"] / "per_isoform.json"
    html = pipeline_inputs["output_dir"] / "report.html"
    h5ad = pipeline_inputs["output_dir"] / "annotated.h5ad"
    assert tsv.exists()
    assert js.exists()
    assert html.exists()
    assert h5ad.exists()

    loaded_tsv = pd.read_csv(tsv, sep="\t")
    assert len(loaded_tsv) == 2
    for col in (
        "transcript_id",
        "completeness",
        "orf_outcome",
        "orf_confidence",
        "nmd_status",
        "iso_utr3_length_nt",
        "polya_signal_motif",
        "denovo_orf_found",
    ):
        assert col in loaded_tsv.columns

    with open(js) as fh:
        records = json.load(fh)
    assert len(records) == 2
    assert {r["transcript_id"] for r in records} == {"t_intact", "t_novel"}


def test_pipeline_intact_isoform_is_propagated_with_high_confidence(
    pipeline_inputs: dict[str, Path],
) -> None:
    result = run_annotate(
        isoforms_path=pipeline_inputs["isoforms"],
        reference_path=pipeline_inputs["reference"],
        output_dir=pipeline_inputs["output_dir"],
        genome_path=pipeline_inputs["genome"],
    )
    intact = result[result["transcript_id"] == "t_intact"].iloc[0]
    assert intact["completeness"] == "full_length"
    assert intact["orf_outcome"] == "propagated_intact"
    assert intact["orf_confidence"] == "high"
    assert float(intact["orf_confidence_score"]) == pytest.approx(1.0)
    assert intact["parent_tx_id"] == "t_parent"


def test_pipeline_novel_isoform_uses_denovo_path(
    pipeline_inputs: dict[str, Path],
) -> None:
    result = run_annotate(
        isoforms_path=pipeline_inputs["isoforms"],
        reference_path=pipeline_inputs["reference"],
        output_dir=pipeline_inputs["output_dir"],
        genome_path=pipeline_inputs["genome"],
    )
    novel = result[result["transcript_id"] == "t_novel"].iloc[0]
    assert novel["completeness"] == "novel_no_match"
    assert novel["orf_outcome"] == "no_parent"
    assert bool(novel["denovo_orf_found"]) is True
    assert int(novel["denovo_orf_aa_length"]) == 60
    assert novel["denovo_start_codon"] == "ATG"
    assert novel["denovo_stop_codon"] == "TAA"
    # Confidence is downgraded to LOW because we have a denovo ORF but no parent.
    assert novel["orf_confidence"] == "low"


def test_pipeline_creates_output_directory(pipeline_inputs: dict[str, Path]) -> None:
    nested = pipeline_inputs["output_dir"] / "deeply" / "nested"
    assert not nested.exists()
    run_annotate(
        isoforms_path=pipeline_inputs["isoforms"],
        reference_path=pipeline_inputs["reference"],
        output_dir=nested,
        genome_path=pipeline_inputs["genome"],
    )
    assert nested.exists()
    assert (nested / "per_isoform.tsv").exists()


def test_pipeline_polya_atlas_takes_precedence_over_motif(tmp_path: Path) -> None:
    """When --polya-atlas hits, polya_evidence_source should be 'polya_db' and
    polya_db_site_id should carry the BED's name column."""
    import pysam as _pysam

    # chr1: 600 bp, no AATAAA motif anywhere (so any 'found' must come from the atlas).
    chr1 = ("N" * 100) + ("G" * 100) + ("N" * 100) + ("G" * 100) + ("N" * 100) + ("G" * 100)
    assert len(chr1) == 600
    fasta_path = tmp_path / "genome.fa"
    fasta_path.write_text(f">chr1\n{chr1}\n")
    _pysam.faidx(str(fasta_path))

    # Reference: parent with 3 exons, last exon to 600.
    ref_rows = [
        'chr1\tCRAFT\tgene\t101\t600\t.\t+\t.\tgene_id "g1";',
        'chr1\tCRAFT\ttranscript\t101\t600\t.\t+\t.\tgene_id "g1"; transcript_id "t_parent";',
        'chr1\tCRAFT\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t_parent";',
        'chr1\tCRAFT\texon\t301\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t_parent";',
        'chr1\tCRAFT\texon\t501\t600\t.\t+\t.\tgene_id "g1"; transcript_id "t_parent";',
        'chr1\tCRAFT\tCDS\t151\t200\t.\t+\t0\tgene_id "g1"; transcript_id "t_parent";',
        'chr1\tCRAFT\tCDS\t301\t400\t.\t+\t0\tgene_id "g1"; transcript_id "t_parent";',
        'chr1\tCRAFT\tCDS\t501\t530\t.\t+\t0\tgene_id "g1"; transcript_id "t_parent";',
    ]
    ref_path = tmp_path / "reference.gtf"
    ref_path.write_text("\n".join(ref_rows) + "\n")

    # Iso: 3'-truncated at GTF 549 (PyRanges End=549, iso 3' end=548). Parent ends
    # at GTF 600 (PyRanges 600); 3'-complete threshold is iso.End >= 550, so 549
    # is TRUNCATED_3P. Iso still covers the parent's stop codon (PyRanges 529).
    iso_rows = [
        'chr1\tCRAFT\texon\t101\t200\t.\t+\t.\ttranscript_id "t_apa";',
        'chr1\tCRAFT\texon\t301\t400\t.\t+\t.\ttranscript_id "t_apa";',
        'chr1\tCRAFT\texon\t501\t549\t.\t+\t.\ttranscript_id "t_apa";',
    ]
    iso_path = tmp_path / "isoforms.gtf"
    iso_path.write_text("\n".join(iso_rows) + "\n")

    # PAS at chr1:540-555 on +; midpoint 547, iso 3' end 548, distance 1 -> match.
    bed_path = tmp_path / "polya.bed"
    bed_path.write_text("chr1\t540\t555\tPAS_TEST_ID\t100\t+\n")

    result = run_annotate(
        isoforms_path=iso_path,
        reference_path=ref_path,
        output_dir=tmp_path / "out",
        genome_path=fasta_path,
        polya_atlas_path=bed_path,
    )
    row = result[result["transcript_id"] == "t_apa"].iloc[0]
    assert row["completeness"] == "alt_3prime_end"
    assert row["orf_outcome"] == "propagated_intact"
    assert row["polya_evidence_source"] == "polya_db"
    assert row["polya_db_site_id"] == "PAS_TEST_ID"


def test_pipeline_polya_atlas_miss_falls_back_to_motif(tmp_path: Path) -> None:
    """When the atlas is provided but the iso 3' end has no PAS hit, the motif
    fallback runs; polya_evidence_source should be either 'canonical_motif' or
    'none' depending on whether a motif is found."""
    import pysam as _pysam

    # chr1 layout matches test_pipeline_reclassifies_alt_polya_isoforms (with AATAAA
    # at positions 565-570 inside iso's last 50 nt window).
    exon3 = "C" * 10 + "AATAAA" + "C" * 49 + "AATAAA" + "C" * 129
    chr1 = ("N" * 100) + ("A" * 100) + ("N" * 100) + ("A" * 100) + ("N" * 100) + exon3
    assert len(chr1) == 700
    fasta_path = tmp_path / "genome.fa"
    fasta_path.write_text(f">chr1\n{chr1}\n")
    _pysam.faidx(str(fasta_path))

    ref_rows = [
        'chr1\tCRAFT\tgene\t101\t700\t.\t+\t.\tgene_id "g1";',
        'chr1\tCRAFT\ttranscript\t101\t700\t.\t+\t.\tgene_id "g1"; transcript_id "t_parent";',
        'chr1\tCRAFT\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t_parent";',
        'chr1\tCRAFT\texon\t301\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t_parent";',
        'chr1\tCRAFT\texon\t501\t700\t.\t+\t.\tgene_id "g1"; transcript_id "t_parent";',
        'chr1\tCRAFT\tCDS\t151\t200\t.\t+\t0\tgene_id "g1"; transcript_id "t_parent";',
        'chr1\tCRAFT\tCDS\t301\t400\t.\t+\t0\tgene_id "g1"; transcript_id "t_parent";',
        'chr1\tCRAFT\tCDS\t501\t530\t.\t+\t0\tgene_id "g1"; transcript_id "t_parent";',
    ]
    ref_path = tmp_path / "reference.gtf"
    ref_path.write_text("\n".join(ref_rows) + "\n")

    iso_rows = [
        'chr1\tCRAFT\texon\t101\t200\t.\t+\t.\ttranscript_id "t_apa";',
        'chr1\tCRAFT\texon\t301\t400\t.\t+\t.\ttranscript_id "t_apa";',
        'chr1\tCRAFT\texon\t501\t580\t.\t+\t.\ttranscript_id "t_apa";',
    ]
    iso_path = tmp_path / "isoforms.gtf"
    iso_path.write_text("\n".join(iso_rows) + "\n")

    # Atlas BED has a PAS far from the iso's 3' end (position 200, not 580).
    bed_path = tmp_path / "polya.bed"
    bed_path.write_text("chr1\t190\t210\tFAR_PAS\t100\t+\n")

    result = run_annotate(
        isoforms_path=iso_path,
        reference_path=ref_path,
        output_dir=tmp_path / "out",
        genome_path=fasta_path,
        polya_atlas_path=bed_path,
    )
    row = result[result["transcript_id"] == "t_apa"].iloc[0]
    # Atlas misses, motif (AATAAA at 565-570) hits -> source should be canonical_motif.
    assert row["polya_evidence_source"] == "canonical_motif"
    assert row["polya_db_site_id"] == ""
    assert row["completeness"] == "alt_3prime_end"


def test_pipeline_reclassifies_alt_polya_isoforms(tmp_path: Path) -> None:
    """When the iso has a canonical poly(A) signal in its last 50 bp, the
    pipeline should reclassify TRUNCATED_3P -> ALT_3PRIME_END and
    STOP_NOT_OBSERVED -> STOP_AT_ALT_POLYA."""
    import pysam as _pysam

    # chr1 is 700 bp. Exon3 region (500-699) carries AATAAA at positions
    # 510-516 (for the premature-stop iso) and 565-571 (for the post-stop iso).
    exon3 = (
        "C" * 10
        + "AATAAA"
        + "C" * 49
        + "AATAAA"
        + "C" * 129
    )
    assert len(exon3) == 200
    chr1 = (
        ("N" * 100)
        + ("A" * 100)
        + ("N" * 100)
        + ("A" * 100)
        + ("N" * 100)
        + exon3
    )
    assert len(chr1) == 700

    fasta_path = tmp_path / "genome.fa"
    fasta_path.write_text(f">chr1\n{chr1}\n")
    _pysam.faidx(str(fasta_path))

    ref_rows = [
        'chr1\tCRAFT\tgene\t101\t700\t.\t+\t.\tgene_id "g1";',
        (
            'chr1\tCRAFT\ttranscript\t101\t700\t.\t+\t.\tgene_id "g1"; '
            'transcript_id "t_parent"; gene_name "GENE1";'
        ),
        (
            'chr1\tCRAFT\texon\t101\t200\t.\t+\t.\tgene_id "g1"; '
            'transcript_id "t_parent"; gene_name "GENE1";'
        ),
        (
            'chr1\tCRAFT\texon\t301\t400\t.\t+\t.\tgene_id "g1"; '
            'transcript_id "t_parent"; gene_name "GENE1";'
        ),
        (
            'chr1\tCRAFT\texon\t501\t700\t.\t+\t.\tgene_id "g1"; '
            'transcript_id "t_parent"; gene_name "GENE1";'
        ),
        (
            'chr1\tCRAFT\tCDS\t151\t200\t.\t+\t0\tgene_id "g1"; '
            'transcript_id "t_parent"; gene_name "GENE1";'
        ),
        (
            'chr1\tCRAFT\tCDS\t301\t400\t.\t+\t0\tgene_id "g1"; '
            'transcript_id "t_parent"; gene_name "GENE1";'
        ),
        (
            'chr1\tCRAFT\tCDS\t501\t530\t.\t+\t0\tgene_id "g1"; '
            'transcript_id "t_parent"; gene_name "GENE1";'
        ),
    ]
    ref_path = tmp_path / "reference.gtf"
    ref_path.write_text("\n".join(ref_rows) + "\n")

    iso_rows = [
        'chr1\tCRAFT\texon\t101\t200\t.\t+\t.\ttranscript_id "t_apa";',
        'chr1\tCRAFT\texon\t301\t400\t.\t+\t.\ttranscript_id "t_apa";',
        'chr1\tCRAFT\texon\t501\t580\t.\t+\t.\ttranscript_id "t_apa";',
        'chr1\tCRAFT\texon\t101\t200\t.\t+\t.\ttranscript_id "t_premature";',
        'chr1\tCRAFT\texon\t301\t400\t.\t+\t.\ttranscript_id "t_premature";',
        'chr1\tCRAFT\texon\t501\t525\t.\t+\t.\ttranscript_id "t_premature";',
    ]
    iso_path = tmp_path / "isoforms.gtf"
    iso_path.write_text("\n".join(iso_rows) + "\n")

    result = run_annotate(
        isoforms_path=iso_path,
        reference_path=ref_path,
        output_dir=tmp_path / "out",
        genome_path=fasta_path,
    )

    t_apa = result[result["transcript_id"] == "t_apa"].iloc[0]
    assert t_apa["completeness"] == "alt_3prime_end"
    assert t_apa["orf_outcome"] == "propagated_intact"
    assert t_apa["orf_confidence"] == "high"
    assert t_apa["parent_gene_name"] == "GENE1"

    t_prem = result[result["transcript_id"] == "t_premature"].iloc[0]
    assert t_prem["completeness"] == "alt_3prime_end"
    assert t_prem["orf_outcome"] == "stop_at_alt_polya"
    assert t_prem["orf_confidence"] == "high"


def test_pipeline_pfam_columns_populate_when_hmm_provided(
    pipeline_inputs: dict[str, Path], tmp_path: Path
) -> None:
    """When --pfam-hmm is supplied, the per-isoform table includes domain columns."""
    from tests.test_pfam import _DNA, _rc, _write_test_hmm

    # Replace chr1 sequence so positions 100-178 encode the test protein, and
    # set t_parent's CDS to that range.
    hmm_path = tmp_path / "test.hmm"
    _write_test_hmm(hmm_path)

    # Rebuild the FASTA with a chr1 that actually encodes the test protein in
    # the 100-178 window so propagation -> protein matches the HMM.
    new_chr1 = ("N" * 100) + _DNA + ("N" * 22)
    chr2_seq = "ATG" + ("GCC" * 59) + "TAA" + ("G" * 17)
    fasta_path = pipeline_inputs["genome"]
    fasta_path.write_text(f">chr1\n{new_chr1}\n>chr2\n{chr2_seq}\n")
    fai = Path(str(fasta_path) + ".fai")
    if fai.exists():
        fai.unlink()
    import pysam
    pysam.faidx(str(fasta_path))

    # Adjust the reference + iso GTFs so chr1 t_parent encodes the matching
    # protein at positions 101..178 (1-based GTF).
    ref_rows = [
        'chr1\tCRAFT\tgene\t101\t200\t.\t+\t.\tgene_id "g1";',
        'chr1\tCRAFT\ttranscript\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t_parent";',
        'chr1\tCRAFT\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t_parent";',
        'chr1\tCRAFT\tCDS\t101\t178\t.\t+\t0\tgene_id "g1"; transcript_id "t_parent";',
    ]
    pipeline_inputs["reference"].write_text("\n".join(ref_rows) + "\n")

    iso_rows = [
        'chr1\tFLAIR\texon\t101\t200\t.\t+\t.\ttranscript_id "t_intact";',
        'chr2\tFLAIR\texon\t1\t200\t.\t+\t.\ttranscript_id "t_novel";',
    ]
    pipeline_inputs["isoforms"].write_text("\n".join(iso_rows) + "\n")

    result = run_annotate(
        isoforms_path=pipeline_inputs["isoforms"],
        reference_path=pipeline_inputs["reference"],
        output_dir=pipeline_inputs["output_dir"],
        genome_path=pipeline_inputs["genome"],
        pfam_hmm_path=hmm_path,
    )
    intact = result[result["transcript_id"] == "t_intact"].iloc[0]
    assert intact["pfam_preserved"] == ["TEST_DOMAIN"]
    assert intact["pfam_lost"] == []

    # Unused import marker.
    _ = _rc


def test_pipeline_skips_isoforms_on_contigs_missing_from_fasta(
    pipeline_inputs: dict[str, Path], capsys: pytest.CaptureFixture[str]
) -> None:
    """Isoforms on contigs absent from the genome FASTA must be dropped up front.

    PacBio collapse outputs reference random/alt contigs that the
    primary_assembly FASTA omits; without this filter any FASTA fetch on those
    isoforms aborts the whole run.
    """
    rows = pipeline_inputs["isoforms"].read_text().splitlines()
    rows.append(
        'chr_missing\tFLAIR\texon\t1\t200\t.\t+\t.\ttranscript_id "t_orphan";'
    )
    pipeline_inputs["isoforms"].write_text("\n".join(rows) + "\n")

    result = run_annotate(
        isoforms_path=pipeline_inputs["isoforms"],
        reference_path=pipeline_inputs["reference"],
        output_dir=pipeline_inputs["output_dir"],
        genome_path=pipeline_inputs["genome"],
    )

    assert set(result["transcript_id"]) == {"t_intact", "t_novel"}
    stderr = capsys.readouterr().err
    assert "Skipping 1 isoforms" in stderr
    assert "chr_missing" in stderr


def test_pipeline_cli_smoke(pipeline_inputs: dict[str, Path]) -> None:
    """Smoke-test the CLI entry point with click's testing harness."""
    from click.testing import CliRunner

    from craft.cli import cli

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "annotate",
            "--isoforms", str(pipeline_inputs["isoforms"]),
            "--reference", str(pipeline_inputs["reference"]),
            "--genome", str(pipeline_inputs["genome"]),
            "--output-dir", str(pipeline_inputs["output_dir"]),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "Annotated 2 isoforms" in result.output
    assert (pipeline_inputs["output_dir"] / "per_isoform.tsv").exists()
