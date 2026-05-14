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
    assert tsv.exists()
    assert js.exists()
    assert html.exists()

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
