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

## [0.1.0] - 2026-05-13

### Added
- Initial repository scaffold: hatchling-based `pyproject.toml`, MIT license, CITATION.cff, .zenodo.json, CHANGELOG.
- Package skeleton under `src/craft/` with typed-signature `NotImplementedError` stubs for every planned module: `completeness`, `orf` (propagation + denovo + confidence), `nmd`, `pfam`, `utr3`, `report` (html + plots), `export` (anndata).
- Click-based CLI with `craft annotate` subcommand stub.
- pytest smoke test asserting package import and `__version__`.
- ruff configuration (line length 100; rule set E, F, W, I, UP, B).
- GitHub Actions test workflow with Python matrix 3.10 / 3.11 / 3.12 and ruff lint.
