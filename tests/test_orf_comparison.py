"""Tests for independent ORF-GTF comparison."""

from pathlib import Path

import pandas as pd

from craft.io.orf_comparison import compare_orf_gtf


def test_compare_orf_gtf_reports_exact_agreement(tmp_path: Path) -> None:
    path = tmp_path / "orfs.gtf"
    path.write_text(
        'chr1\tORFanage\tCDS\t11\t30\t.\t+\t0\tgene_id "g"; transcript_id "tx";\n'
        'chr1\tORFanage\tstart_codon\t11\t13\t.\t+\t0\tgene_id "g"; transcript_id "tx";\n'
        'chr1\tORFanage\tstop_codon\t31\t33\t.\t+\t0\tgene_id "g"; transcript_id "tx";\n'
    )
    craft = pd.DataFrame(
        [{
            "transcript_id": "tx",
            "resolved_start_pos": 10,
            "resolved_stop_codon_pos": 30,
            "resolved_cds_bp": 20,
        }]
    )
    row = compare_orf_gtf(craft, path).iloc[0]
    assert bool(row["comparator_start_agrees"])
    assert bool(row["comparator_stop_agrees"])
    assert row["comparator_cds_bp_delta"] == 0

