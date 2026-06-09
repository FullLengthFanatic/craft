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
- `--group-by obs_column` with `--counts`, aggregates functional consequences per
  cell group into `per_celltype_consequence.tsv`.
- `--pfam-hmm Pfam-A.hmm` enables Pfam domain preservation analysis (slow with
  full Pfam; a future release will switch to `hmmscan` against a pressed database).
- `--polya-atlas sites.bed` provides curated polyA sites (PolyASite v3.0,
  PolyA_DB v4, or any user-supplied BED 6+). When supplied, atlas hits drive
  the ALT_3PRIME_END / STOP_AT_ALT_POLYA reclassification with the canonical
  poly(A) motif scan as fallback. **Pre-filter the atlas by usage score**
  (`awk '$5 >= 0.01'` for PolyASite v3.0) before passing — the unfiltered
  atlas is too dense (~one PAS every 200 bp) and produces uninformative
  98% match rates. See [`docs/user_guide.md`](docs/user_guide.md) for the
  BED format spec, recommended sources, and the full stringency story.
- `--no-coding-potential` turns off the reference-calibrated coding-potential
  score (on by default). Threshold flags (`--ptc-threshold-nt`, `--min-orf-aa`,
  ...) are listed in [`docs/features.md`](docs/features.md#command-line-options).
- `--classification sqanti_classification.txt` joins SQANTI3/pigeon columns
  (default `structural_category`) onto the output by transcript_id, so you can
  cut novel-boundary classes against CRAFT's consequence calls.

Runtime on chr22 of a real PacBio Iso-Seq sample (~13k isoforms): ~2 minutes.
Full-genome scale runs end-to-end: the bcM0003 PacBio Iso-Seq sample (698,049
isoforms, with `--polya-atlas` and coding potential on, no `--pfam-hmm`) takes
~1h45m wall and ~19 GB RAM on a 32-core VM.

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
| `per_isoform.tsv`             | per-iso annotation table, 60 columns, list columns JSON-encoded             |
| `per_isoform.json`            | same content as records; list columns stay as lists                         |
| `report.html`                 | self-contained interactive report (summary cards + plotly + table)          |
| `annotated.h5ad`              | AnnData with iso annotations in `var`, per-cell counts in `X` (if given)    |
| `per_celltype_consequence.tsv`| with `--counts --group-by`: consequence fractions per cell group            |

Every column is documented in [`docs/features.md`](docs/features.md).

CRAFT classifies each ORF by geometric propagation (`orf_outcome`,
`propagated_cds_*`) and then reconstructs the spliced CDS to find the real stop
(`resolved_orf_status`, `intron_retained_in_cds`, ...). NMD and UTR consequences
are computed once from the resolved ORF (single `nmd_status`, with a de-novo
fallback for orphans recorded in `nmd_basis`).

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
cases where the start or stop codon falls outside the read. Then reconstruct the
iso's own spliced CDS and walk it to the real in-frame stop, which catches
frameshifts, exon-skip premature stops, and introns retained in the CDS
(`resolved_*` columns). Apply NMD rules once to the resolved stop
(50nt PTC rule + start-proximal, long-last-exon, and last-exon escapes), falling
back to the de-novo ORF for orphan isoforms (`nmd_basis`). Compute 3'/5' UTR length deltas and scan
for poly(A) signals. Score each ORF for coding potential against a model
self-calibrated to the reference. Optionally scan the translated CDS against a
Pfam HMM database. Emit TSV, JSON, HTML report, and AnnData.

For the full algorithm, threshold defaults, and design rationale, see
[`docs/methods.md`](docs/methods.md).

## What CRAFT does *not* do

- It does not do structural QC. The iso GTF is assumed to be post-QC (e.g.,
  pigeon, SQANTI3-curated). CRAFT describes what's there, it doesn't filter.
- It does not call cell types. With `--group-by` it summarises functional
  consequences per existing cell grouping, but clustering/cell-typing is upstream.
- It does not harmonise chromosome naming. All three inputs must agree.

Since v1.5 it detects intron retention inside the CDS and the resulting
premature stops (`intron_retained_in_cds`, `resolved_orf_status`).

## Status

v1.7.0. The pipeline runs end-to-end on real long-read isoform GTFs at
full-genome scale. Validated on a PacBio Iso-Seq sample (chr22 subset and the
full bcM0003 sample, ~698k isoforms). Methods paper in preparation: benchmarking
reference-isoform ORF propagation vs de-novo prediction on simulated truncated
reads.

## Citation

See [`CITATION.cff`](CITATION.cff).

## License

MIT. See [`LICENSE`](LICENSE).
