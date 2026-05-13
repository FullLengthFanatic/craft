"""De novo ORF prediction for genuinely novel isoforms (no usable parent).

Uses orfipy with ATG-only start codons. For each isoform, builds the transcript
sequence in 5'-to-3' orientation (reverse-complementing minus-strand exon
sequence), scans for ORFs, picks the longest one above the minimum-length
threshold, and maps the ORF's transcript coordinates back to genomic intervals.
"""

from pathlib import Path

import orfipy_core
import pandas as pd
import pyranges as pr
import pysam

_RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_COLUMNS = (
    "transcript_id",
    "denovo_orf_found",
    "denovo_cds_bp",
    "denovo_cds_intervals",
    "denovo_orf_aa_length",
    "denovo_start_codon",
    "denovo_stop_codon",
)


def _reverse_complement(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]


def _transcript_sequence(
    exons: pd.DataFrame,
    strand: str,
    genome: pysam.FastaFile,
) -> str:
    chrom = str(exons["Chromosome"].iloc[0])
    sorted_exons = exons.sort_values("Start").reset_index(drop=True)
    parts: list[str] = []
    for _, ex in sorted_exons.iterrows():
        parts.append(genome.fetch(chrom, int(ex["Start"]), int(ex["End"])))
    sequence = "".join(parts).upper()
    if strand == "-":
        sequence = _reverse_complement(sequence)
    return sequence


def _transcript_to_genomic_intervals(
    orf_start: int,
    orf_end: int,
    exons: pd.DataFrame,
    strand: str,
    chrom: str,
) -> list[tuple[str, int, int, str]]:
    """Map [orf_start, orf_end) in transcript coords back to genomic intervals."""
    sorted_exons = exons.sort_values("Start").reset_index(drop=True)
    if strand == "+":
        walk = list(sorted_exons.iterrows())
    elif strand == "-":
        walk = list(sorted_exons.iterrows())[::-1]
    else:
        raise ValueError(f"Unsupported strand: {strand!r}")

    intervals: list[tuple[str, int, int, str]] = []
    cum_tx = 0
    for _, ex in walk:
        ex_start = int(ex["Start"])
        ex_end = int(ex["End"])
        ex_tx_end = cum_tx + (ex_end - ex_start)
        overlap_start = max(cum_tx, orf_start)
        overlap_end = min(ex_tx_end, orf_end)
        if overlap_end > overlap_start:
            if strand == "+":
                gen_start = ex_start + (overlap_start - cum_tx)
                gen_end = ex_start + (overlap_end - cum_tx)
            else:
                rel_start = overlap_start - cum_tx
                rel_end = overlap_end - cum_tx
                gen_start = ex_end - rel_end
                gen_end = ex_end - rel_start
            intervals.append((chrom, gen_start, gen_end, strand))
        cum_tx = ex_tx_end

    intervals.sort(key=lambda x: x[1])
    return intervals


def _parse_codon(description: str, field: str) -> str:
    marker = f"{field}:"
    if marker not in description:
        return ""
    rest = description.split(marker, 1)[1]
    return rest.split(";", 1)[0].strip()


def _empty_row(tx_id: str) -> dict:
    return {
        "transcript_id": tx_id,
        "denovo_orf_found": False,
        "denovo_cds_bp": 0,
        "denovo_cds_intervals": [],
        "denovo_orf_aa_length": 0,
        "denovo_start_codon": "",
        "denovo_stop_codon": "",
    }


def predict(
    isoforms: pr.PyRanges,
    genome_fasta: Path,
    min_orf_aa: int = 50,
) -> pd.DataFrame:
    """Predict ORFs de novo for isoforms (typically those with no usable parent).

    Args:
        isoforms: PyRanges of isoform exons (with ``transcript_id`` and ``Strand``).
        genome_fasta: Path to indexed genome FASTA (``.fai`` expected alongside).
        min_orf_aa: Minimum ORF length in amino acids. Default 50.

    Returns:
        DataFrame with one row per isoform and columns ``transcript_id``,
        ``denovo_orf_found``, ``denovo_cds_bp``, ``denovo_cds_intervals``,
        ``denovo_orf_aa_length``, ``denovo_start_codon``, ``denovo_stop_codon``.
    """
    if len(isoforms) == 0:
        return pd.DataFrame(columns=list(_COLUMNS))

    iso_df = isoforms.df
    iso_strand = iso_df.groupby("transcript_id")["Strand"].first().to_dict()
    iso_chrom = iso_df.groupby("transcript_id")["Chromosome"].first().to_dict()
    iso_exons_by_tx = {tx: g for tx, g in iso_df.groupby("transcript_id", sort=False)}
    min_bp = min_orf_aa * 3

    rows: list[dict] = []
    with pysam.FastaFile(str(genome_fasta)) as genome:
        for tx_id, iso_exons in iso_exons_by_tx.items():
            strand = str(iso_strand[tx_id])
            chrom = str(iso_chrom[tx_id])
            seq = _transcript_sequence(iso_exons, strand, genome)

            candidates = list(
                orfipy_core.orfs(
                    seq,
                    minlen=min_bp,
                    starts=["ATG"],
                    strand="f",
                )
            )
            if not candidates:
                rows.append(_empty_row(tx_id))
                continue

            candidates.sort(key=lambda o: o[1] - o[0], reverse=True)
            orf_start, orf_stop, _, description = candidates[0]
            cds_bp = orf_stop - orf_start
            aa_length = cds_bp // 3
            start_codon = _parse_codon(description, "Start")
            stop_codon = _parse_codon(description, "Stop")
            intervals = _transcript_to_genomic_intervals(
                orf_start, orf_stop, iso_exons, strand, chrom
            )
            rows.append(
                {
                    "transcript_id": tx_id,
                    "denovo_orf_found": True,
                    "denovo_cds_bp": cds_bp,
                    "denovo_cds_intervals": intervals,
                    "denovo_orf_aa_length": aa_length,
                    "denovo_start_codon": start_codon,
                    "denovo_stop_codon": stop_codon,
                }
            )

    return pd.DataFrame(rows, columns=list(_COLUMNS))
