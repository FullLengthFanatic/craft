# CRAFT

**Coding Region Annotation From Templates**

![status](https://img.shields.io/badge/status-v2_development-orange)
![license](https://img.shields.io/badge/license-MIT-blue)
![python](https://img.shields.io/badge/python-3.10%2B-blue)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20844836.svg)](https://doi.org/10.5281/zenodo.20844836)

CRAFT is an evidence-aware functional annotation toolkit for long-read
single-cell transcriptomes. It accepts isoforms from Pigeon/SQANTI3, IsoQuant,
Bambu, FLAIR, FLAMES, or another caller and keeps four questions separate:

1. How well is the transcript structure supported by molecules and orthogonal evidence?
2. Which transcript ends—and therefore which CDS boundaries—were actually observed?
3. What coding sequence is supported by the reference and by the observed sequence?
4. Which RNA-surveillance rules could apply if that isoform exists in the cell?

That separation is the central v2 design decision. Cell recurrence is descriptive,
not a probability that an isoform is real. An internal ATG is not silently promoted
to the biological start. An NMD rule match is reported as predicted susceptibility,
not observed degradation.

CRAFT never deletes input rows. It exposes evidence, ambiguity, censoring, and
limitations so that filtering can be calibrated on positive and negative controls.

## Install

```bash
pip install -e ".[dev]"
```

Python 3.10 or newer is required.

## Quick start

```bash
craft annotate \
  --isoforms isoforms.gtf \
  --reference gencode.annotation.gtf \
  --genome genome.fa \
  --classification pigeon_classification.txt \
  --evidence-table molecule_evidence.tsv \
  --output-dir craft_out
```

Useful optional inputs:

- `--counts counts.h5ad` adds per-cell abundance and detection summaries.
- `--classification FILE` imports SQANTI3/Pigeon structural class and parent/gene hints.
- `--evidence-table FILE` imports molecule-level and artifact evidence and emits a
  transparent, explicitly uncalibrated evidence score and tier.
- `--orf-comparator-gtf FILE` compares CRAFT with an independent CDS caller such as
  ORFanage without treating either answer as truth.
- `--polya-atlas FILE` adds curated cleavage/polyadenylation-site evidence.
- `--pfam-hmm Pfam-A.hmm` compares predicted protein-domain content.
- `--infer-alternative-start` opts into downstream in-frame ATG inference for
  5'-censored transcripts. It is off by default.

Run `craft annotate --help` for the complete option list.

## Inputs

| Input | Required | Purpose |
| --- | --- | --- |
| Isoform GTF | yes | Exon structures keyed by `transcript_id` |
| Reference GTF | yes | Exons, CDS, and preferably explicit start/stop codons and tags |
| Genome FASTA | yes | Sequence reconstruction; indexed automatically when possible |
| Classification table | no | Structural categories plus parent/gene hints |
| Evidence table | no | Unique/ambiguous support, junction, end, artifact, and replicate evidence |
| Per-cell counts | no | Abundance, detection fraction, and group summaries |
| Independent ORF GTF | no | Caller agreement/disagreement analysis |

All genomic inputs must use the same assembly and chromosome naming convention.

## Outputs

| File | Description |
| --- | --- |
| `per_isoform.tsv` | Versioned per-isoform annotation; optional inputs add populated evidence blocks |
| `per_isoform.json` | Same records with list-valued fields preserved |
| `report.html` | Self-contained exploratory report |
| `annotated.h5ad` | Annotation in `var` and, when supplied, per-cell counts in `X` |
| `per_celltype_consequence.tsv` | Molecule-weighted group summaries with `--counts --group-by` |
| `per_celltype_as_nmd.tsv` | Evidence-supported predicted AS–NMD candidates by group |

The authoritative field definitions are in [`docs/features.md`](docs/features.md).

## Interpretation

The main v2 fields are:

- `isoform_evidence_tier`: independent structure evidence; `insufficient_evidence`
  means missing inputs, not a failed transcript.
- `parent_ambiguous`, `parent_match_margin`: whether reference-parent assignment is
  unique enough to support propagation.
- `orf_censoring`: `none`, `left`, `right`, or `both`; a missing transcript end is
  treated as censoring, not as a discovered CDS boundary.
- `resolved_orf_status`: includes `left_censored` and `right_censored` rather than
  inventing a start or calling a terminal truncation a complete ORF.
- `nmd_susceptibility`, `nmd_evidence_tier`, `surveillance_limitations`: structural
  prediction, strength of supporting annotation, and reasons not to over-interpret it.
- `nonstop_decay_candidate`: a separate flag for right-censored/stopless molecules.

Example downstream selection:

```python
import pandas as pd

df = pd.read_csv("craft_out/per_isoform.tsv", sep="\t")

# Structure-supported isoforms with a non-ambiguous propagated parent.
supported = df[
    df["isoform_evidence_tier"].isin(["strong", "supported"])
    & ~df["parent_ambiguous"].fillna(True)
]

# Predicted NMD candidates with complete ORF boundaries and stronger evidence.
nmd_candidates = supported[
    (supported["nmd_susceptibility"] == "likely_sensitive")
    & supported["nmd_evidence_tier"].isin(["high", "moderate"])
    & (supported["orf_censoring"] == "none")
]
```

Do not copy these cuts blindly. The evidence score is not calibrated to a target
false-discovery rate; tune thresholds with spike-ins, technical replicates,
reference holdouts, and plausible artifact decoys from the assay being analysed.

## Algorithm in brief

CRAFT ranks parent candidates using exact/contained intron-chain agreement,
junction precision and recall, overlap, upstream parent/gene hints, and curated
reference metadata. It reports ambiguity instead of hiding close ties. It projects
the chosen reference CDS, reconstructs the observed spliced sequence, and records
whether start and stop codons are observed or censored. Novel/orphan isoforms get a
de novo ORF candidate, clearly labelled as such. Rule-based RNA-surveillance
susceptibility is then reported with its evidence tier and limitations. Coding
potential and Pfam provide supporting features, not ground truth.

See [`docs/v2_design.md`](docs/v2_design.md) for the scientific contract and
literature context, [`docs/migration_v2.md`](docs/migration_v2.md) for changed
semantics, and [`docs/user_guide.md`](docs/user_guide.md) for operations.

## Scope and limitations

- CRAFT annotates and ranks evidence; it is not itself a long-read aligner,
  transcript assembler, cell caller, or universally calibrated isoform filter.
- Reference annotation is evidence, not truth. Reference completeness, phase,
  APPRIS/MANE/CCDS/basic tags, and parent ambiguity are retained in the output.
- Structural rules cannot establish actual NMD. Validation requires perturbation,
  half-life, translation, or other orthogonal data.
- A 5'-truncated long read cannot by itself identify the true translation start.
- The supplied evidence score is transparent and useful for triage, but requires
  assay-specific calibration before being interpreted as a decision boundary.

## Development status

The branch is v2.0.0 development code. The scientific semantics have intentionally
changed from v1; read the migration guide before comparing old and new tables.

## Citation and license

See [`CITATION.cff`](CITATION.cff). CRAFT is MIT licensed; see [`LICENSE`](LICENSE).
