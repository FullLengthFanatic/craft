# Changelog

All notable changes to CRAFT will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.6.0] - 2026-06-05

Two additive functional-annotation features for novel/orphan isoforms, plus
their documentation. The existing columns are unchanged; 68 columns total.

### Added
- Coding-potential scoring (`core/coding_potential.py`), self-calibrated to the
  supplied reference. CRAFT trains a hexamer coding/non-coding log-likelihood
  table plus a logistic regression on four features (hexamer LLR, log10 ORF
  length, ORF coverage, Fickett TESTCODE) using the reference's own CDS-bearing
  transcripts as the coding set and CDS-less transcripts as the non-coding set,
  then scores every isoform's best ORF (resolved -> propagated -> de novo). No
  model is shipped and no external tool is required. New columns
  `coding_potential_score`, `coding_potential_label`, `coding_potential_orf_source`;
  the fitted model and a held-out AUC are written to `coding_potential_model.json`.
  On chr22 the held-out AUC is 0.88; de-novo orphan ORFs are 13% coding (the gate
  for trusting their NMD calls vs flagging lncRNA). Default on;
  `--no-coding-potential` disables it, and it auto-skips when the reference has no
  non-coding transcripts.
- De novo NMD: `nmd.predict_denovo` applies the escape-rule cascade to the de
  novo ORF, so orphan isoforms (`no_parent` / `no_parent_cds` / `start_lost`)
  with a predicted ORF now get an NMD call instead of `not_applicable`. New
  columns `nmd_status_denovo`, `nmd_rule_denovo`, `nmd_confidence_denovo`
  (confidence always `low`: the stop is from a predicted ORF, not a reference).
  This is the call to use for the orphan tail of novel-not-in-catalog (NNC)
  isoforms. On chr22, 592 orphans become NMD-sensitive that previously had no NMD
  call; orphans with no de novo ORF stay `not_applicable` (lncRNA candidates).
  `docs/features.md` documents the columns and a recipe to coalesce resolved +
  de novo into a single `nmd_call`.

### Packaging
- Version bumped `1.5.1` -> `1.6.0` in `src/craft/__init__.py` and `CITATION.cff`.

## [1.5.1] - 2026-06-03

### Performance
- `core/orf/resolve.py`: rewrote the per-isoform hot path to operate on
  precomputed numpy exon arrays instead of per-row pandas operations
  (`iterrows`/`itertuples`, repeated `sort_values`, `.iloc`). On chr22 the
  resolve step drops from ~103 s to ~6 s (17x) and the end-to-end run from
  ~2m53s to ~1m50s. Output is byte-identical across all 62 columns; this is a
  pure speedup with no behaviour change, which matters for full-genome inputs.

### Packaging
- Version bumped `1.5.0` -> `1.5.1` in `src/craft/__init__.py` and `CITATION.cff`.

## [1.5.0] - 2026-06-02

Sequence-aware functional consequences. All additions are additive: the original
39 `per_isoform` columns keep identical values (verified byte/float-identical
against the chr22 `craft_out_v13_atlas` baseline), so prior benchmark numbers
stay reproducible. 23 new columns and one new output file.

### Added
- `core/orf/resolve.py`: sequence-level ORF resolution. For each isoform with an
  observed parent start, CRAFT reconstructs the isoform's own spliced CDS and
  walks it in 3-nt codons to the first in-frame stop, so it finds the real stop
  through frameshifts, exon-skip premature stops, and retained CDS introns.
  New columns: `resolved_orf_status`, `resolved_stop_pos`, `resolved_cds_bp`,
  `resolved_aa_length`, `resolved_cds_intervals`, `ptc_introduced`,
  `intron_retained_in_cds`, `frame_consistent`, `stop_in_transcript`.
- `nmd.predict_resolved`: NMD on the resolved stop -> `nmd_status_resolved`,
  `nmd_rule_resolved`, `nmd_confidence_resolved`. On chr22, 739 isoforms flip to
  NMD-sensitive vs the geometric call (499 from `escaped`, 240 from
  `not_applicable`), the frameshift/PTC cases the geometric method could not see.
- uORF detection and a long-3'UTR flag as advisory NMD branches: `uorf_count`,
  `uorf_triggers_nmd`, `long_utr3_triggers_nmd`.
- 5'UTR length deltas: `iso_utr5_length_nt`, `parent_utr5_length_nt`,
  `utr5_length_delta_nt`, `utr5_length_delta_pct`; resolved 3'UTR deltas:
  `iso_utr3_length_resolved_nt`, `utr3_length_delta_resolved_nt`,
  `utr3_length_delta_pct_resolved`.
- `export/celltype.py` + `--group-by`: expression-weighted per-cell-group
  consequence fractions written to `per_celltype_consequence.tsv` and stored in
  `annotated.h5ad` `uns['celltype_consequences']`.
- `has_cds_bearing_parent` column and an opt-in `--prefer-coding-parent` flag
  that breaks parent-selection ties toward CDS-bearing transcripts (off by
  default to preserve reproducibility).
- CLI options exposing previously hardcoded thresholds: `--tolerance`,
  `--ptc-threshold-nt`, `--start-proximal-nt`, `--long-last-exon-nt`,
  `--min-orf-aa`, `--orf-high-confidence`, `--orf-medium-confidence`,
  `--long-utr3-nt`. Defaults equal the prior constants.
- `docs/features.md`: complete column dictionary covering every output column.

### Changed
- `pfam.scan` now takes the isoform protein from the resolved CDS when available
  (frameshift- and intron-retention-aware), falling back to the propagated then
  de novo CDS. Domain comparison columns are unchanged in name.
- `docs/methods.md`: the v1 "Known limitations" items now implemented (intron
  retention, frame tracking, 5'UTR, cell-type aggregation) moved into a new
  sequence-resolution section.

### Tests
- New: `tests/test_orf_resolve.py` (9), `tests/test_celltype.py` (6), plus
  resolved-NMD, resolved/5'UTR, and parent-tiebreak cases. Suite: 187 passing,
  ruff clean.

## [1.4.1] - 2026-05-25

Docs-only patch release. Addresses review findings on the methods-paper draft committed in v1.4.

### Changed
- `docs/methods_paper.md`:
  - Verified GSE86148 dataset citation. Replaced wrong attribution (Boehm 2021 / Lykke-Andersen lab) with Colombo M, Karousis ED, Bourquin J, Bruggmann R, Mühlemann O (2017) *RNA* 23:189-201; Mühlemann lab at Bern. PMC5238794.
  - Merged duplicate refs 5 and 6 (both de la Fuente 2020 *Genome Biology*) into a single ref covering tappAS and the IsoAnnotLite companion tool. Renumbered in-text refs 6-24 down by 1; total references now 23.
  - Added a Bench 3 positive-control comparison using GENCODE's curated `nonsense_mediated_decay` annotation: GENCODE NMD OR=1.415 (p=1.3e-24) vs CRAFT NMD-sensitive OR=1.457 on a 4× smaller positive set, with 3,724-transcript overlap (96.8% of CRAFT-sensitive transcripts are also GENCODE-NMD). New Table 3b.
  - Added drop-out accounting paragraph: 80,441 NMD universe → 47,378 eligible (28,261 transcripts dropped by salmon low-count filter, 4,802 dropped because NMD rules don't apply).
  - Annotated Table 1 caption as "8 of 12 cells" with pointer to the full per-cell cache.
  - Added markdown image embeds for all three figures (`![Figure N](../benchmarks/figures/...)`).

### Packaging
- Package version bumped from `1.4` to `1.4.1` in `src/craft/__init__.py` and `CITATION.cff`. No source code changes; tests + ruff unchanged from v1.4.

## [1.4] - 2026-05-25

### Packaging
- Package version bumped from `1.3` to `1.4` in `src/craft/__init__.py` and `CITATION.cff`. `craft --version` now reports `1.4`.

### Added
- `benchmarks/cbench/` dev-only package (truncation simulator, GENCODE loader, orfipy wrapper, scoring, plotly figure helpers). 26 tests.
- Three runnable benchmark scripts under `benchmarks/`: `run_bench1.py` (simulated truncation, full grid), `run_bench2.py` (real-data ORF concordance on bcM0003), `run_bench3_universe.py` + `run_bench3_quant.py` + `run_bench3_analysis.py` (NMD-target enrichment in GSE86148 UPF1-KD bulk RNA-seq).
- All benchmark figures committed under `benchmarks/figures/` as PNG + plotly JSON.
- `bench` optional-dependency group (`jupyter`, `ipykernel`, `nbformat`).

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

### Bench 2: real-data ORF concordance on bcM0003 (2026-05-24)

Real PacBio sc Iso-Seq stand-in for the originally-planned MAS-Seq paired bulk + sc design (the MAS-Seq paper turned out to use tumor T cells, not WTC11 - the plan was based on a wrong premise). Bench 2 instead asks: "on the user's own bcM0003 sc Iso-Seq (698,049 isoforms), when each iso retains a recoverable ORF from its GENCODE parent, does CRAFT's propagation outperform orfipy de novo?"

`benchmarks/run_bench2.py` reuses the cached bcM0003 CRAFT output (`craft_out_atlas_filtered/per_isoform.tsv`) and the GENCODE v45 pool pickle from Bench 1. For each iso it projects the parent's `start_codon` + `stop_codon` genomic positions onto the iso's exon structure; isos where either codon falls outside the iso's exons are dropped (truth not recoverable). 223,976 eligible isos remain. Runs orfipy on iso transcript-orientation sequences (batched, 18 s), parses CRAFT's `propagated_cds_intervals` and maps both calls into iso transcript coordinates, scores against truth, stratifies by `completeness`.

Start-codon exact-match rate by completeness category:

| completeness | n | CRAFT | orfipy | gap |
|---|---|---|---|---|
| full_length        | 95,395 | 1.000 | 0.727 | +27.3 ppt |
| alt_3prime_end     | 79,511 | 1.000 | 0.765 | +23.5 ppt |
| truncated_5p       | 23,576 | 1.000 | 0.740 | +26.0 ppt |
| truncated_3p       |  4,410 | 1.000 | 0.851 | +14.9 ppt |
| truncated_both     |  9,506 | 1.000 | 0.790 | +21.0 ppt |
| internal_fragment  | 11,578 | 1.000 | 0.812 | +18.8 ppt |

Caveat. CRAFT's near-perfect score is partially tautological: truth is defined as "parent's GENCODE CDS projected onto the iso", and CRAFT's propagated coordinates come from the same parent. The honest claim from this bench is the orfipy bar: **on real long-read sc isoforms, de novo prediction misses the start codon 15-28% of the time across every completeness category, even when the ORF is fully recoverable**. Propagation closes that gap. The right panel of the figure (CRAFT accuracy by `orf_confidence`) is uninformative for the same reason and is kept only as a sanity check.

Figure committed at `benchmarks/figures/bench2_concordance_panel.{png,json}`; per-row scored DataFrame cached at `benchmarks/cache/bench2/bench2_scores.tsv.gz` (gitignored). Total wall: 3 min (parse 22 s, truth projection 2.4 min, orfipy 18 s, score + figure < 30 s).

## [0.1.0] - 2026-05-13

### Added
- Initial repository scaffold: hatchling-based `pyproject.toml`, MIT license, CITATION.cff, .zenodo.json, CHANGELOG.
- Package skeleton under `src/craft/` with typed-signature `NotImplementedError` stubs for every planned module: `completeness`, `orf` (propagation + denovo + confidence), `nmd`, `pfam`, `utr3`, `report` (html + plots), `export` (anndata).
- Click-based CLI with `craft annotate` subcommand stub.
- pytest smoke test asserting package import and `__version__`.
- ruff configuration (line length 100; rule set E, F, W, I, UP, B).
- GitHub Actions test workflow with Python matrix 3.10 / 3.11 / 3.12 and ruff lint.
