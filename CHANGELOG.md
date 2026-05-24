# Changelog

All notable changes to CRAFT will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- `core/polya_atlas.match_iso_end` now takes a pre-built index from `build_atlas_index` and uses `numpy.searchsorted` for O(log n) lookup per isoform instead of the previous linear pandas filter. Chr22 runtime with a filtered PolyASite v3.0 atlas drops from ~3 min to <30 s; the perf fix is what makes the full-genome run practical (commit `ce02461`).

### Fixed
- `pipeline.run_annotate` now filters isoforms whose chromosome is absent from the genome FASTA before any per-isoform processing. PacBio collapse outputs reference random/alt contigs (`chr*_KI270*_random`, `chrUn_*`) that the GRCh38 primary_assembly FASTA omits; any downstream `pysam.fetch` on those isoforms raised `KeyError` and aborted the run mid-pipeline. The filter logs the count + first 5 contig names to stderr and proceeds (commit `3a25d9c`).

### Verified on full bcM0003 sample (2026-05-22)

Two full-genome runs on PacBio Iso-Seq sample bcM0003 (698,049 isoforms after skipping 673 isos on 22 contigs missing from the FASTA). Both runs took ~71 min wall, ~18 GB peak RSS.

| Completeness | Baseline (no atlas) | Atlas-filtered (PolyASite v3.0, score ≥ 0.01) |
|---|---|---|
| full_length | 25.1% | 25.1% |
| truncated_5p | 12.1% | 12.1% |
| truncated_3p | **20.5%** | **3.4%** |
| truncated_both | 4.4% | 4.4% |
| internal_fragment | 12.7% | 12.7% |
| alt_3prime_end | **13.9%** | **31.0%** |
| novel_no_match | 11.3% | 11.3% |

The atlas reclassifies 119,255 isoforms from `truncated_3p` to `alt_3prime_end`; structural categories are unchanged.

ORF confidence: HIGH rises from 26.4% to 35.4%, LOW drops from 63.4% to 53.5%. `polya_evidence_source` is 85.1% `polya_db` / 4.1% `canonical_motif` / 10.8% `none` with the atlas, vs 0% / 52.9% / 47.1% without. NMD distribution is unchanged (4.1% sensitive, 32.6% escaped, 63.4% not_applicable) because NMD depends on whether the propagated stop is observed in the read, not on PAS evidence.

### Benchmarks

- New dev-only `benchmarks/cbench/` package: shared truncation simulator, GENCODE loader, orfipy wrapper, scoring, and plotly figure helpers for the three methods-paper benchmarks.
- 26 cbench tests (truncation simulator + metrics), main suite still at 164 passing, ruff clean.
- `benchmarks/run_bench1.py` runs the full 4 × 3 × 3 grid of (rate, orientation, seed) cells over GENCODE v45 protein-coding transcripts.

### Bench 1: simulated truncation vs orfipy (2026-05-24)

Pool: 63,332 GENCODE v45 protein-coding transcripts with complete CDS. 36 cells = {5, 10, 25, 50}% truncation × {5', 3', both} orientation × 3 seeds; 3,000 transcripts sampled per cell. Total wall time 8.6 min on the dev VM, 92,718 scored rows. Per-cell scores cached under `benchmarks/cache/bench1_scores/` (regenerable from the runner; gitignored).

Scoring is restricted to transcripts whose ground-truth ORF survives the truncation intact (both start and stop codons retained in the truncated sequence); the per-cell `n_intact` counts are reported.

| rate | orientation | mean n | CRAFT start_exact | orfipy start_exact | CRAFT \|len_err\| | orfipy \|len_err\| |
|------|-------------|--------|---|---|---|---|
| 5%   | 3'          | 2,701  | **1.000** | 0.947 |  0.0 nt | 11.0 nt |
| 5%   | 5'          | 1,796  | **0.999** | 0.948 |  1.2 nt |  8.1 nt |
| 5%   | both        | 2,263  | **0.999** | 0.948 |  0.6 nt | 11.6 nt |
| 10%  | 3'          | 2,445  | **1.000** | 0.944 |  0.0 nt | 10.9 nt |
| 10%  | 5'          | 1,047  | **0.996** | 0.945 |  3.4 nt |  8.1 nt |
| 10%  | both        | 1,635  | **0.999** | 0.945 |  1.1 nt | 11.8 nt |
| 25%  | 3'          | 1,748  | **1.000** | 0.945 |  0.0 nt | 11.5 nt |
| 25%  | 5'          |   253  | **0.990** | 0.936 |  5.8 nt |  7.9 nt |
| 25%  | both        |   604  | **0.994** | 0.938 |  3.5 nt | 11.5 nt |
| 50%  | 3'          |   842  | **1.000** | 0.945 |  0.0 nt | 12.2 nt |
| 50%  | 5'          |    35  | 0.983     | 0.962 |  4.6 nt |  3.4 nt |
| 50%  | both        |    84  | **0.985** | 0.937 |  5.7 nt |  9.4 nt |

CRAFT propagation hits the start codon **0.98-1.00** across the grid; orfipy bottoms out at **0.94-0.95** consistently. CRAFT's mean ORF length error is **0.0 nt for every 3'-truncated cell** and 1-6 nt elsewhere; orfipy sits at 8-12 nt across all conditions. The single cell where orfipy edges out CRAFT (50% / 5', n=35) is too small to read into - one CRAFT failure flips the rate by 2.8 ppt.

Figure committed at `benchmarks/figures/bench1_recovery_panel.{png,json}`.

### Bench 3: NMD-target enrichment in UPF1-KD bulk RNA-seq (2026-05-24)

Three-step pipeline:

1. **NMD universe.** `benchmarks/run_bench3_universe.py` filters GENCODE v45 to 80,441 transcripts (`protein_coding` + `nonsense_mediated_decay`, complete CDS) and runs CRAFT on the filtered iso GTF against the full GENCODE reference. CRAFT labels 7,282 transcripts (9.1%) as NMD-sensitive, 72,906 as escaped, 253 as not_applicable. ~16 min wall.
2. **Salmon quant.** `benchmarks/run_bench3_quant.py` pulls 6 GSE86148 (HeLa, Lykke-Andersen lab, SRP083135) samples - 3 scr controls (SRR4081222-224) + 3 UPF1 KDs (SRR4081225-227) - and runs salmon transcript-level quant against the GENCODE v45 transcriptome. Mapping rates 83.9-85.9% across all 6 samples, no batch outliers. ~70 min wall.
3. **DE + enrichment.** `benchmarks/run_bench3_analysis.py` runs pydeseq2 (Wald test, `~condition`, UPF1-KD vs control), joins per-transcript log2FC + padj against the CRAFT NMD universe, and tests whether NMD-sensitive transcripts are enriched among UPF1-KD-upregulated (log2FC >= 1, padj < 0.05). 2x2 contingency on 47,378 eligible transcripts (propagated_intact or disrupted, NMD label sensitive or escaped, complete DE results):

|              | upregulated | not upregulated |
|---|---|---|
| **NMD-sensitive** |   620 |  3,228 |
| **NMD-escaped**   | 5,070 | 38,460 |

Fisher's exact (one-sided, alternative=greater): **odds ratio 1.457, p = 2.4e-15**.

NMD-sensitive transcripts are upregulated under UPF1 KD at 16.1% vs 11.6% for NMD-escaped transcripts. The effect size is modest (the rule cascade is structural, not biochemical, and 48h siRNA KD has secondary effects beyond NMD), but the directionality and significance match the methods-paper claim that CRAFT's NMD labels track real NMD biology.

Figure committed at `benchmarks/figures/bench3_enrichment_panel.{png,json}`; raw contingency + odds ratio in `benchmarks/figures/bench3_enrichment.tsv`. Quant summary in `benchmarks/figures/bench3_quant_summary.tsv`.

## [0.1.0] - 2026-05-13

### Added
- Initial repository scaffold: hatchling-based `pyproject.toml`, MIT license, CITATION.cff, .zenodo.json, CHANGELOG.
- Package skeleton under `src/craft/` with typed-signature `NotImplementedError` stubs for every planned module: `completeness`, `orf` (propagation + denovo + confidence), `nmd`, `pfam`, `utr3`, `report` (html + plots), `export` (anndata).
- Click-based CLI with `craft annotate` subcommand stub.
- pytest smoke test asserting package import and `__version__`.
- ruff configuration (line length 100; rule set E, F, W, I, UP, B).
- GitHub Actions test workflow with Python matrix 3.10 / 3.11 / 3.12 and ruff lint.
