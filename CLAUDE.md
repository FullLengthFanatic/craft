# CRAFT contributor context

CRAFT is an evidence-aware long-read isoform functional annotator. Current
version: `2.0.0`.

## Scientific contracts

- Describe and annotate; never silently delete input isoforms.
- Keep structure evidence, abundance/recurrence, ORF inference, and RNA-surveillance
  prediction separate.
- Treat transcript ends and CDS boundaries as observed or censored. Do not infer a
  biological start from an internal ATG unless the user explicitly opts in.
- Reference projection is preferred when the parent assignment and CDS are adequate,
  but reference annotation is evidence rather than truth.
- Report close parent ties, incomplete CDS metadata, and phase problems.
- Never describe recurrence, coding-potential output, the evidence score, or the NMD
  rule score as a calibrated probability unless a calibration study has established it.
- Any empirical number in the docs must trace to a committed analysis and output.

## Code map

- `src/craft/cli.py`: public options.
- `src/craft/pipeline.py`: orchestration and canonical column order.
- `src/craft/core/reference.py`: reference CDS/tag quality metadata.
- `src/craft/core/completeness.py`: parent ranking and transcript-end classification.
- `src/craft/core/orf/`: projection, sequence resolution, de novo fallback, confidence.
- `src/craft/core/evidence.py`: transparent external molecule/read evidence aggregation.
- `src/craft/core/nmd.py`: rule-based surveillance susceptibility.
- `src/craft/io/orf_comparison.py`: independent ORF-caller comparison.
- `docs/v2_design.md`: authoritative scientific semantics.
- `docs/features.md`: field dictionary.
- `docs/migration_v2.md`: v1-to-v2 behavior changes.

## Commands

```bash
pip install -e ".[dev]"
pytest
ruff check .
craft annotate --help
```

Python is 3.10+. Keep dependencies bounded where upstream format/API changes can
break runtime compatibility. Inputs must use the same genome assembly and contig naming.

## Development priorities

The largest remaining scientific task is prospective calibration on FLIGHT-seq:
positive controls, negative/artifact controls, held-out genes/samples, technical
replicates, and orthogonal junction/end/translation evidence. Until that exists,
preserve transparent component features and avoid learned probability language.
