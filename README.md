# CRAFT

**Coding Region Annotation From Templates**

![status](https://img.shields.io/badge/status-pre--alpha-orange)
![license](https://img.shields.io/badge/license-MIT-blue)
![python](https://img.shields.io/badge/python-3.10%2B-blue)

A Python toolkit for long-read isoform functional-consequence annotation. Takes
the output of any long-read isoform caller (FLAIR, IsoQuant, Bambu, FLAMES,
SQANTI3, isoseq+pigeon) and emits per-isoform structural completeness, ORF
status, NMD susceptibility, 3' UTR features, and (optionally) Pfam domain
preservation.

The methods novelty is **reference-isoform ORF propagation with truncation-aware
confidence**. When an iso is truncated but still covers the parent transcript's
CDS region, CRAFT projects the parent's ORF coordinates onto the iso and flags
exactly where the call becomes uncertain (start codon outside the read, stop
codon outside the read, structural divergence in CDS). De-novo ORF prediction is
used only as a fallback for genuinely novel isoforms.

For the full method, every category definition, and the rationale for every
design choice, read [`docs/methods.md`](docs/methods.md).

## Install

```bash
pip install -e ".[dev]"
```

Requires Python ≥ 3.10. Runtime dependencies: pysam, pyranges, pandas, plotly,
pyhmmer, orfipy, anndata, click. No R, no Java, no scanpy.

## Quick start

```bash
craft annotate \
    --isoforms  /path/to/iso.gtf \
    --reference /path/to/gencode.v45.annotation.gtf \
    --genome    /path/to/GRCh38.fa \
    --output-dir out/
```

Optional flags:

- `--counts h5ad_file_or_10x_mtx_dir` per-cell counts; populates `annotated.h5ad`.
- `--pfam-hmm Pfam-A.hmm` enables Pfam domain preservation analysis (slow with
  full Pfam; v1.5 will switch to `hmmscan` against a pressed database).
- `--polya-atlas sites.bed` provides curated polyA sites (PolyASite v3.0,
  PolyA_DB v4, or any user-supplied BED 6+). When supplied, atlas hits drive
  the ALT_3PRIME_END / STOP_AT_ALT_POLYA reclassification with the canonical
  poly(A) motif scan as fallback. **Pre-filter the atlas by usage score**
  (`awk '$5 >= 0.01'` for PolyASite v3.0) before passing — the unfiltered
  atlas is too dense (~one PAS every 200 bp) and produces uninformative
  98% match rates. See [`docs/user_guide.md`](docs/user_guide.md) for the
  BED format spec, recommended sources, and the full stringency story.

Runtime on chr22 of a real PacBio Iso-Seq sample (~13k isoforms): ~1 minute.
Full-genome scale (~600k iso rows) is roughly 10-15 minutes without `--pfam-hmm`.

## Inputs

| Input            | Required | Format                                           |
| ---------------- | -------- | ------------------------------------------------ |
| Isoform GTF      | yes      | GTF with `exon` rows and a `transcript_id` attr  |
| Reference GTF    | yes      | GTF with `exon` AND `CDS` rows (GENCODE/Ensembl) |
| Genome FASTA     | yes      | indexed (`.fai` built on-the-fly if missing)     |
| Per-cell counts  | no       | `.h5ad` or 10x-style MTX directory               |
| Pfam HMM         | no       | pyhmmer-compatible `.hmm` (e.g. `Pfam-A.hmm`)    |

All three required inputs must use the same chromosome naming (`chr1` vs `1`).

## Outputs

| File                          | Description                                                                 |
| ----------------------------- | --------------------------------------------------------------------------- |
| `per_isoform.tsv`             | per-iso annotation table, 62 columns, list columns JSON-encoded             |
| `per_isoform.json`            | same content as records; list columns stay as lists                         |
| `report.html`                 | self-contained interactive report (summary cards + plotly + table)          |
| `annotated.h5ad`              | AnnData with iso annotations in `var`, per-cell counts in `X` (if given)    |
| `per_celltype_consequence.tsv`| with `--counts --group-by`: consequence fractions per cell group            |

Every column is documented in [`docs/features.md`](docs/features.md).

CRAFT reports each isoform's ORF two ways: the original **geometric** projection
(`orf_outcome`, `nmd_status`, ...) and a v1.5 **sequence-resolved** view
(`resolved_orf_status`, `nmd_status_resolved`, `intron_retained_in_cds`, ...)
that reads the spliced sequence to find the real stop. Prefer the resolved
columns for functional-consequence calls.

## Filter recipes

```python
import pandas as pd
df = pd.read_csv("out/per_isoform.tsv", sep="\t")

# trustworthy ORF calls
df[df["orf_confidence"].isin(["high", "medium"])]

# biological NMD substrates (not truncation artefacts)
df[(df["nmd_status"] == "sensitive") & (df["nmd_confidence"] == "high")]

# alternative 3' UTR isoforms with a canonical poly(A) signal
df[(df["utr3_length_delta_nt"] != 0) & (df["polya_signal_motif"] != "")]

# novel coding isoforms supported by a de-novo ORF
df[(df["orf_outcome"] == "no_parent") & (df["denovo_orf_found"])]

# APA isoforms (alternative 3' end backed by a poly(A) signal)
df[df["completeness"] == "alt_3prime_end"]
```

## How it works (one paragraph)

For each iso: pick the best-matching reference parent by maximal splice-junction
sharing. Classify the iso's structural completeness from its end positions vs
the parent's. Propagate the parent's CDS coordinates onto the iso, flagging
cases where the start or stop codon falls outside the read. Apply NMD rules to
the resulting stop position (50nt PTC rule + start-proximal, long-last-exon,
and last-exon escapes). Compute 3' UTR length delta and scan for canonical
poly(A) signals. Optionally translate the propagated CDS and scan against a
Pfam HMM database. Emit TSV, JSON, HTML report, and AnnData.

For the full algorithm, threshold defaults, and design rationale, see
[`docs/methods.md`](docs/methods.md).

## What CRAFT does *not* do

- It does not do structural QC. The iso GTF is assumed to be post-QC (e.g.,
  pigeon, SQANTI3-curated). CRAFT describes what's there, it doesn't filter.
- It does not call cell types. With `--group-by` it summarises functional
  consequences per existing cell grouping, but clustering/cell-typing is upstream.
- It does not harmonise chromosome naming. All three inputs must agree.

As of v1.5 it does detect intron retention inside the CDS and the resulting
premature stops (`intron_retained_in_cds`, `resolved_orf_status`).

## Status

Pre-alpha (v0.1.0). The pipeline composes end-to-end on real long-read isoform
GTFs. Smoke-tested on a PacBio Iso-Seq sample (chr22 subset, ~13k isoforms).
Methods paper in preparation: benchmarking reference-isoform ORF propagation vs
de-novo prediction on simulated truncated reads.

## Citation

See [`CITATION.cff`](CITATION.cff).

## License

MIT. See [`LICENSE`](LICENSE).
