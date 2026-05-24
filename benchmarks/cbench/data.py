"""GENCODE transcript loader for Bench 1.

Stream-parses a GENCODE GTF + genome FASTA, yields ``Transcript`` records for
protein-coding transcripts whose CDS is complete (both start and stop codons
present, no ``cds_start_NF`` / ``cds_end_NF`` tag). Each record carries enough
state to (a) write the original transcript out as a GTF for CRAFT, (b) extract
the transcript-orientation sequence for orfipy, and (c) score predictions
against the GENCODE-truth ORF.
"""

import gzip
import random
from collections import defaultdict
from collections.abc import Callable, Iterator
from dataclasses import dataclass, field
from pathlib import Path

import pysam

_RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _rc(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]


@dataclass
class Transcript:
    """A GENCODE-truth protein-coding transcript with extracted sequence + ORF coords."""

    transcript_id: str
    gene_id: str
    chrom: str
    strand: str
    exons: list[tuple[int, int]]  # half-open, genomic, sorted by start
    cds_genomic_intervals: list[tuple[int, int]]
    start_codon_genomic: tuple[int, int]  # 3 bp interval
    stop_codon_genomic: tuple[int, int]  # 3 bp interval
    transcript_seq: str  # transcript orientation, uppercase
    cds_tx_start: int  # transcript coord (0-based) of the A in ATG
    cds_tx_end: int  # transcript coord (exclusive) past the last base of stop codon
    extras: dict[str, str] = field(default_factory=dict)

    @property
    def transcript_length(self) -> int:
        return len(self.transcript_seq)

    @property
    def cds_length(self) -> int:
        return self.cds_tx_end - self.cds_tx_start


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _parse_attrs(attr_str: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for chunk in attr_str.strip().rstrip(";").split(";"):
        chunk = chunk.strip()
        if not chunk or " " not in chunk:
            continue
        key, val = chunk.split(" ", 1)
        out[key] = val.strip().strip('"')
    return out


def _genomic_to_tx_coord(
    pos: int, exons: list[tuple[int, int]], strand: str
) -> int:
    """Map a 0-based genomic position to the 0-based transcript coordinate.

    ``pos`` must fall within one of the exons. Raises ``ValueError`` otherwise.
    """
    ordered = exons if strand == "+" else list(reversed(exons))
    cursor = 0
    for start, end in ordered:
        if start <= pos < end:
            return cursor + (pos - start if strand == "+" else end - 1 - pos)
        cursor += end - start
    raise ValueError(f"position {pos} not in any exon")


def _protein_coding_with_complete_cds(buf: dict) -> bool:
    tx = buf.get("transcript")
    if tx is None:
        return False
    attrs = tx["attrs"]
    if attrs.get("transcript_type") != "protein_coding":
        return False
    if not buf["exon"] or not buf["CDS"]:
        return False
    if not buf["start_codon"] or not buf["stop_codon"]:
        return False
    tag_str = attrs.get("tag", "")
    tags = {t.strip() for t in tag_str.split(",")}
    if "cds_start_NF" in tags or "cds_end_NF" in tags:
        return False
    return True


def _iter_transcript_buffers(
    gtf_path: Path,
    predicate: Callable[[dict], bool] | None = None,
) -> Iterator[dict]:
    """Stream a GTF, yield per-transcript dicts of feature lists.

    GENCODE GTFs are grouped by transcript; we lean on that ordering instead of
    a full PyRanges load. Each yielded dict has keys ``transcript`` (single
    record) plus ``exon``, ``CDS``, ``start_codon``, ``stop_codon`` (lists).
    """
    current_tx: str | None = None
    buffer: dict = _empty_buffer()
    with _open_text(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            feature = cols[2]
            if feature not in ("exon", "CDS", "start_codon", "stop_codon", "transcript"):
                continue
            attrs = _parse_attrs(cols[8])
            tx_id = attrs.get("transcript_id")
            if tx_id is None:
                continue
            if tx_id != current_tx:
                if current_tx is not None and (predicate is None or predicate(buffer)):
                    yield buffer
                current_tx = tx_id
                buffer = _empty_buffer()
            record = {
                "chrom": cols[0],
                "start": int(cols[3]) - 1,
                "end": int(cols[4]),
                "strand": cols[6],
                "attrs": attrs,
            }
            if feature == "transcript":
                buffer["transcript"] = record
            else:
                buffer[feature].append(record)
        if current_tx is not None and (predicate is None or predicate(buffer)):
            yield buffer


def _empty_buffer() -> dict:
    return {
        "transcript": None,
        "exon": [],
        "CDS": [],
        "start_codon": [],
        "stop_codon": [],
    }


def _buffer_to_transcript(buf: dict, genome: pysam.FastaFile) -> Transcript | None:
    tx = buf["transcript"]
    attrs = tx["attrs"]
    chrom = tx["chrom"]
    strand = tx["strand"]

    if chrom not in genome.references:
        return None

    exons = sorted((e["start"], e["end"]) for e in buf["exon"])
    cds = sorted((c["start"], c["end"]) for c in buf["CDS"])

    pieces = [genome.fetch(chrom, s, e).upper() for s, e in exons]
    full = "".join(pieces)
    if strand == "-":
        full = _rc(full)
    if "N" in full[:10] or len(full) < 60:
        return None

    starts = sorted((c["start"], c["end"]) for c in buf["start_codon"])
    stops = sorted((c["start"], c["end"]) for c in buf["stop_codon"])
    start_codon = (starts[0][0], starts[-1][1])
    stop_codon = (stops[0][0], stops[-1][1])

    # GENCODE convention: CDS records exclude the stop codon (it's a separate
    # ``stop_codon`` feature). ``cds_tx_end`` is the first base of the stop
    # codon in transcript orientation -> CDS spans [cds_tx_start, cds_tx_end).
    # Both CRAFT (which inherits GENCODE's intervals) and orfipy (default flag
    # set) report the same convention, so truth must match.
    if strand == "+":
        cds_tx_start = _genomic_to_tx_coord(start_codon[0], exons, strand)
        cds_tx_end = _genomic_to_tx_coord(stop_codon[0], exons, strand)
    else:
        cds_tx_start = _genomic_to_tx_coord(start_codon[1] - 1, exons, strand)
        cds_tx_end = _genomic_to_tx_coord(stop_codon[1] - 1, exons, strand)

    if cds_tx_end <= cds_tx_start or cds_tx_end > len(full):
        return None

    return Transcript(
        transcript_id=attrs["transcript_id"],
        gene_id=attrs.get("gene_id", ""),
        chrom=chrom,
        strand=strand,
        exons=exons,
        cds_genomic_intervals=cds,
        start_codon_genomic=start_codon,
        stop_codon_genomic=stop_codon,
        transcript_seq=full,
        cds_tx_start=cds_tx_start,
        cds_tx_end=cds_tx_end,
        extras={
            "gene_name": attrs.get("gene_name", ""),
            "transcript_name": attrs.get("transcript_name", ""),
        },
    )


def load_protein_coding_transcripts(
    gtf_path: Path,
    fasta_path: Path,
    min_cds_bp: int = 150,
    max_n: int | None = None,
    seed: int = 0,
) -> list[Transcript]:
    """Return up to ``max_n`` protein-coding transcripts with complete CDS.

    Args:
        gtf_path: GENCODE GTF (plain or ``.gz``).
        fasta_path: indexed genome FASTA (must have ``.fai`` sibling).
        min_cds_bp: minimum CDS length to keep (default 150 bp = 50 aa).
        max_n: if given, return a random sample of this size (seeded). If
            ``None``, return all matches.
        seed: random seed for sampling when ``max_n`` is given.

    Returns:
        List of ``Transcript`` records.
    """
    genome = pysam.FastaFile(str(fasta_path))
    try:
        kept: list[Transcript] = []
        for buf in _iter_transcript_buffers(gtf_path, _protein_coding_with_complete_cds):
            tx = _buffer_to_transcript(buf, genome)
            if tx is None or tx.cds_length < min_cds_bp:
                continue
            kept.append(tx)
    finally:
        genome.close()
    if max_n is not None and len(kept) > max_n:
        rng = random.Random(seed)
        kept = rng.sample(kept, max_n)
    return kept


def transcripts_by_gene(
    transcripts: list[Transcript],
) -> dict[str, list[Transcript]]:
    out: dict[str, list[Transcript]] = defaultdict(list)
    for tx in transcripts:
        out[tx.gene_id].append(tx)
    return dict(out)
