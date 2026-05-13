# Changelog

All notable changes to CRAFT will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2026-05-13

### Added
- Initial repository scaffold: hatchling-based `pyproject.toml`, MIT license, CITATION.cff, .zenodo.json, CHANGELOG.
- Package skeleton under `src/craft/` with typed-signature `NotImplementedError` stubs for every planned module: `completeness`, `orf` (propagation + denovo + confidence), `nmd`, `pfam`, `utr3`, `report` (html + plots), `export` (anndata).
- Click-based CLI with `craft annotate` subcommand stub.
- pytest smoke test asserting package import and `__version__`.
- ruff configuration (line length 100; rule set E, F, W, I, UP, B).
- GitHub Actions test workflow with Python matrix 3.10 / 3.11 / 3.12 and ruff lint.
