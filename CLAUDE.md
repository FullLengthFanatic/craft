# CLAUDE.md

Context for working on CRAFT that isn't obvious from the code. Keep this short; the
docs below carry the detail.

## What CRAFT is

Long-read isoform functional-consequence annotator. Sits downstream of an isoform
caller (isoseq+pigeon, FLAIR, IsoQuant, Bambu, FLAMES, SQANTI3). Per isoform it emits
structural completeness, ORF status (propagated + sequence-resolved + de-novo
fallback), NMD, UTR features, coding potential, optional Pfam, and per-cell
recurrence. Output is a 63-column table + HTML report + AnnData. Current version
v1.8.0. MIT. Public repo: https://github.com/FullLengthFanatic/craft.

## Two design commitments (don't violate these)

- **Describe, don't filter.** CRAFT never deletes isoform rows. The input GTF is
  assumed structurally QC'd; CRAFT describes what's there and hands the user
  confidence + recurrence columns to filter on. A dropped row can't be recovered
  downstream; a column can always be filtered.
- **Trust the reference over de-novo.** Propagate a curated parent CDS onto truncated
  novel isoforms; de-novo ORF prediction is only the orphan fallback. The benchmarks
  exist to justify this (de-novo misses the start codon on truncated reads).

## Project decisions / state

- **GitHub-only, no paper** (decided at v1.7.0, see CHANGELOG). Don't try to revive
  the manuscript.
- `docs/craft_explained.md` is the technical reference (how/why per stage, design
  rationale, interpretation FAQ, code anchors). The old `methods.md` and the
  academic-manuscript draft `methods_paper.md` were removed: `methods.md`'s content
  folded into `craft_explained.md`, and CRAFT is GitHub-only (no paper planned).
- Development is complete and validated (bench1/2/3 + the BD70 single-cell run).

## Doc map

Four docs, one job each:
- `docs/whitepaper.md` - plain-language primer for newcomers (the "quick view").
- `docs/user_guide.md` - how to run CRAFT end to end (operational reference).
- `docs/craft_explained.md` - detailed, code-anchored design guide: how/why per
  stage + snippets + file:line + permalinks pinned to the v1.8.0 commit, the 18 CLI
  knobs vs source-only params, strengths/limitations vs tools, interpretation FAQ.
- `docs/features.md` - per-column dictionary (authoritative for the 63 columns).

## Hard rule on numbers

Any number that appears in docs must trace to a committed script/output:
`benchmarks/figures/` (bench1/2/3) or `test_dataset/input_files/analysis/`
(BD70 recurrence/recovery: `recurrence_v18.py`, `recurrence_pigeon_min20.py`,
`recovery_gencode.py`). Verify before writing; never invent or eyeball a figure.

## Code layout

`src/craft/`: `cli.py` (Click entry, all flags) -> `pipeline.py::run_annotate`
(orchestration, `_OUTPUT_COLUMNS` = the canonical 63-column order). Stages under
`core/`: `intervals`, `completeness`, `orf/{propagation,denovo,resolve,confidence}`,
`nmd`, `utr3`, `coding_potential`, `pfam`, `polya_atlas`, `recurrence`. IO in `io/`,
writers in `export/` and `report/`.

## Parameters

- CLI flags (the user-tunable surface) are all in `cli.py`. Numeric defaults:
  `--tolerance` 50, `--ptc-threshold-nt` 50, `--start-proximal-nt` 150,
  `--long-last-exon-nt` 400, `--min-orf-aa` 50, `--orf-high-confidence` 0.85,
  `--orf-medium-confidence` 0.5, `--long-utr3-nt` 1000.
- Source-only (no flag): confidence base/factor tables (`orf/confidence.py`),
  poly(A) atlas window 24nt (`polya_atlas.py`), motif window 50nt + the 11
  `POLYA_SIGNALS` (`utr3.py`), ATG-only de-novo starts (`orf/denovo.py`),
  coding-potential cutoff 0.5 and train cap 4000 (`coding_potential.py`),
  Pfam evalue<=0.01 (`pfam.py`).
- Full knob inventory with line numbers: `docs/craft_explained.md` section 4.

## Commands

```bash
pip install -e ".[dev]"     # install (needs Python >= 3.10)
pytest                       # main test suite
ruff check .                 # lint (config in pyproject.toml; line-length 100)
craft annotate --isoforms iso.gtf --reference gencode.gtf --genome genome.fa --output-dir out/
```

Inputs must share chromosome naming (`chr1` vs `1`); CRAFT does not harmonise it.
Pre-filter a poly(A) atlas by usage score (`awk '$5 >= 0.01'` for PolyASite v3.0)
before passing `--polya-atlas`.

## Open questions (where work would help)

Calibrated recurrence threshold (vs the fixed `n_cells_detected >= 3`); unique vs
ambiguous read support (to collapse the recovery floor/ceiling range); validating
intron-retention NMD calls against IR-PSI tables; NMD beyond the 4 structural rules.
