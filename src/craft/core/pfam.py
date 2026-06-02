"""Pfam domain disruption via local pyhmmer scanning.

For each isoform, translates the propagated (or de novo, when no propagation
exists) CDS to protein, scans it against a Pfam-format HMM database, and
compares the resulting domain set against the parent reference transcript's
domain set. Hits are cached by SHA256 of the protein sequence so repeated
proteins (common across cells in single-cell data) scan only once.
"""

import hashlib
from pathlib import Path

import pandas as pd
import pyhmmer
import pyranges as pr
import pysam

_RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")

_CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

_COLUMNS = [
    "transcript_id",
    "iso_pfam_domains",
    "parent_pfam_domains",
    "pfam_preserved",
    "pfam_lost",
    "pfam_gained",
]


def _reverse_complement(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]


def translate(nt_seq: str) -> str:
    """Translate a nucleotide sequence to single-letter amino-acid codes.

    Reads codons in frame 0; stops at the first stop codon. Unknown codons
    (containing N or any non-standard base) become ``X``.
    """
    aa: list[str] = []
    seq = nt_seq.upper()
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        if len(codon) < 3:
            break
        c = _CODON_TABLE.get(codon, "X")
        if c == "*":
            break
        aa.append(c)
    return "".join(aa)


def _extract_protein(
    intervals: list[tuple[str, int, int, str]],
    genome: pysam.FastaFile,
) -> str:
    if not intervals:
        return ""
    chrom = intervals[0][0]
    strand = intervals[0][3]
    sorted_intervals = sorted(intervals, key=lambda x: x[1])
    parts = [genome.fetch(chrom, int(s), int(e)) for _, s, e, _ in sorted_intervals]
    nt = "".join(parts).upper()
    if strand == "-":
        nt = _reverse_complement(nt)
    return translate(nt)


class _PfamScanner:
    """Cached HMM scanner. Caches per-protein hit sets by SHA256."""

    def __init__(self, hmms: list[pyhmmer.plan7.HMM]):
        self._alphabet = pyhmmer.easel.Alphabet.amino()
        self._hmms = hmms
        self._cache: dict[str, frozenset[str]] = {}

    def scan(self, protein: str) -> frozenset[str]:
        if not protein:
            return frozenset()
        key = hashlib.sha256(protein.encode()).hexdigest()
        if key in self._cache:
            return self._cache[key]
        seq = pyhmmer.easel.TextSequence(
            name=b"query", sequence=protein
        ).digitize(self._alphabet)
        domain_names: set[str] = set()
        for top in pyhmmer.hmmsearch(self._hmms, [seq]):
            for hit in top:
                if hit.included:
                    domain_names.add(str(top.query.name))
                    break
        result = frozenset(domain_names)
        self._cache[key] = result
        return result


def _parent_cds_lookup(reference: pr.PyRanges) -> dict[str, list[tuple]]:
    df = reference.df
    if "Feature" not in df.columns or "transcript_id" not in df.columns:
        return {}
    cds = df[df["Feature"] == "CDS"]
    out: dict[str, list[tuple]] = {}
    for tx, group in cds.groupby("transcript_id", sort=False):
        out[str(tx)] = [
            (str(r["Chromosome"]), int(r["Start"]), int(r["End"]), str(r["Strand"]))
            for _, r in group.iterrows()
        ]
    return out


def scan(
    per_isoform: pd.DataFrame,
    reference: pr.PyRanges,
    pfam_hmm: Path,
    genome_fasta: Path,
) -> pd.DataFrame:
    """Scan iso and parent CDS for Pfam domains and emit per-isoform comparison.

    Args:
        per_isoform: DataFrame with at least ``transcript_id``, ``parent_tx_id``,
            ``propagated_cds_intervals``, ``denovo_cds_intervals``. Typically
            the merged output of the propagation + denovo steps from
            :func:`craft.pipeline.run_annotate`.
        reference: Reference annotation PyRanges with a ``Feature`` column
            (CDS rows are used).
        pfam_hmm: Path to a ``.hmm`` file (e.g., Pfam-A.hmm).
        genome_fasta: Path to indexed genome FASTA.

    Returns:
        Per-isoform DataFrame with ``transcript_id``, ``iso_pfam_domains``,
        ``parent_pfam_domains``, ``pfam_preserved``, ``pfam_lost``,
        ``pfam_gained`` (each domain column is a sorted list of HMM names).
    """
    if per_isoform.empty:
        return pd.DataFrame(columns=_COLUMNS)

    with pyhmmer.plan7.HMMFile(str(pfam_hmm)) as f:
        hmms = list(f)
    scanner = _PfamScanner(hmms)

    parent_cds_by_tx = _parent_cds_lookup(reference)

    rows: list[dict] = []
    with pysam.FastaFile(str(genome_fasta)) as genome:
        for _, prop_row in per_isoform.iterrows():
            tx_id = prop_row["transcript_id"]
            parent_tx = prop_row.get("parent_tx_id", "") or ""

            # Prefer the sequence-resolved CDS (frameshift- and intron-retention-
            # aware) when available; fall back to the geometric propagated CDS,
            # then the de novo ORF.
            iso_intervals = prop_row.get("resolved_cds_intervals") or []
            if not iso_intervals:
                iso_intervals = prop_row.get("propagated_cds_intervals") or []
            if not iso_intervals:
                iso_intervals = prop_row.get("denovo_cds_intervals") or []

            iso_protein = _extract_protein(iso_intervals, genome) if iso_intervals else ""
            iso_domains = scanner.scan(iso_protein)

            parent_domains: frozenset[str] = frozenset()
            if parent_tx and parent_tx in parent_cds_by_tx:
                parent_protein = _extract_protein(parent_cds_by_tx[parent_tx], genome)
                parent_domains = scanner.scan(parent_protein)

            rows.append(
                {
                    "transcript_id": tx_id,
                    "iso_pfam_domains": sorted(iso_domains),
                    "parent_pfam_domains": sorted(parent_domains),
                    "pfam_preserved": sorted(iso_domains & parent_domains),
                    "pfam_lost": sorted(parent_domains - iso_domains),
                    "pfam_gained": sorted(iso_domains - parent_domains),
                }
            )

    return pd.DataFrame(rows, columns=_COLUMNS)
