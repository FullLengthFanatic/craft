# CRAFT

**Coding Region Annotation From Templates**

A Python toolkit for long-read isoform functional-consequence annotation with truncation-aware ORF propagation.

## Status

Pre-alpha (v0.1.0). Repository scaffold only; modules are typed stubs awaiting implementation. See `/home/simone.picelli/.claude/plans/i-want-you-to-foamy-micali.md` for the full v1 plan.

## Overview

CRAFT takes long-read isoform calls (FLAIR, IsoQuant, Bambu, FLAMES, SQANTI3 output) and emits per-isoform functional consequences: ORF prediction (with reference-isoform propagation under truncation), NMD susceptibility, Pfam domain disruption, and 3' UTR feature changes. Outputs include a per-isoform TSV, an AnnData/MuData file for downstream single-cell analysis, and an interactive HTML report.

The methods novelty is **reference-isoform ORF propagation with truncation-aware confidence**: novel single-cell long-read isoforms commonly arrive truncated, and CRAFT propagates ORF coordinates from the parent annotated transcript where structure is preserved, with explicit confidence flags that drop when key structural regions (start codon, stop codon, last junction) are not observed in the read.

## Installation

```bash
pip install -e ".[dev]"
```

(Conda and Docker installs to come.)

## Quick start

```bash
craft annotate \
    --isoforms isoforms.gtf \
    --reference gencode.v44.annotation.gtf \
    --genome hg38.fa \
    --output-dir out/
```

## Outputs

- `per_isoform.tsv` — every functional annotation with confidence flags
- `per_isoform.json` — same data, programmatic consumers
- `annotated.h5ad` — AnnData with isoforms as `var`, functional features as `var` columns
- `report.html` — self-contained interactive plotly report

## Citation

See `CITATION.cff`.

## License

MIT. See `LICENSE`.
