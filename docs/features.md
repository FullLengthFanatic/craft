# CRAFT feature & column reference

This is the single place that documents every output column CRAFT produces and
how each is computed. For the algorithmic rationale and threshold justifications,
see [`methods.md`](methods.md). For operational recipes, see
[`user_guide.md`](user_guide.md).

## The one thing to understand first: geometric vs resolved

CRAFT annotates each isoform's ORF in two ways, and both are written to every
row.

**Geometric (the original v1.4 view).** The parent transcript's CDS coordinates
are projected onto the isoform by interval intersection. It is fast and it is
exact when the isoform's structure matches the parent's, but it never reads the
spliced nucleotide sequence. It uses the *end of the intersected interval* as
the stop. Columns: `orf_outcome`, `propagated_cds_*`, `nmd_status`, `nmd_rule`,
`utr3_length_delta_nt`, and the rest of the v1.4 set.

**Resolved (added in v1.5).** CRAFT reconstructs the isoform's own spliced CDS
from the genome and translates it codon by codon from the parent's start to the
first in-frame stop. This finds the *real* stop, so it catches frameshifts from
alternative splice sites, premature stops from exon skips, and introns retained
inside the CDS. Columns: `resolved_*`, `nmd_status_resolved`, the `*_resolved`
UTR columns, plus the `uorf_*` and `long_utr3_triggers_nmd` advisory flags.

For functional-consequence work, prefer the resolved columns. The geometric
columns are kept unchanged for speed and for continuity with earlier results.
The two agree for structurally intact isoforms and diverge exactly where the
isoform's reading frame departs from the reference, which is where it matters.

## Output files

| File | When | Contents |
| --- | --- | --- |
| `per_isoform.tsv` | always | One row per isoform, 68 columns. List-valued columns are JSON-encoded. |
| `per_isoform.json` | always | Same content as records; list columns stay as lists. |
| `report.html` | always | Self-contained interactive summary. |
| `annotated.h5ad` | always | AnnData: annotations in `var`, per-cell counts in `X` (if `--counts`). |
| `per_celltype_consequence.tsv` | with `--counts` + `--group-by` | Expression-weighted consequence fractions per cell group. |
| `coding_potential_model.json` | unless `--no-coding-potential` | Fitted coding-potential model: feature weights, training counts, 5-fold cross-validated AUC. |

Every `per_isoform` column is listed below, grouped by feature.

## Structural completeness

The isoform's shape relative to its best-matching reference parent. The parent
is the reference transcript with the most exactly shared splice junctions, ties
broken by exon-overlap bp.

| Column | Type | Meaning |
| --- | --- | --- |
| `transcript_id` | str | Isoform ID from the input GTF. |
| `completeness` | categorical | `full_length`, `truncated_5p`, `truncated_3p`, `truncated_both`, `internal_fragment`, `alt_3prime_end`, or `novel_no_match`. `alt_3prime_end` is assigned when a `truncated_3p` isoform's 3' end has poly(A) support. |
| `parent_tx_id` | str | Selected reference parent transcript; empty for `novel_no_match`. |
| `parent_gene_id` | str | Parent gene ID (from the reference). |
| `parent_gene_name` | str | Parent gene name (from the reference). |
| `shared_junctions` | int | Count of splice junctions shared exactly with the parent. |
| `parent_overlap_bp` | int | Total stranded exon-overlap bp with the parent. |
| `has_cds_bearing_parent` | bool | Whether the selected parent has CDS records. Informational; with `--prefer-coding-parent` it also breaks parent-selection ties toward coding transcripts. |

## ORF: geometric propagation

Parent CDS coordinates projected onto the isoform by intersection.

| Column | Type | Meaning |
| --- | --- | --- |
| `orf_outcome` | categorical | `propagated_intact`, `disrupted`, `start_lost`, `stop_not_observed`, `stop_at_alt_polya`, `no_parent`, `no_parent_cds`. |
| `propagated_cds_bp` | int | Parent CDS bp preserved in the isoform. |
| `parent_cds_bp` | int | Total parent CDS bp. |
| `start_codon_covered` | bool | Parent start-codon genomic position observed in the isoform. |
| `stop_codon_covered` | bool | Parent stop-codon genomic position observed in the isoform. |
| `propagated_cds_intervals` | list | Genomic `[chrom, start, end, strand]` of the projected CDS. |

## ORF: sequence resolution (v1.5)

Reconstructed from the isoform's own spliced sequence. Computed for every
isoform whose parent start codon is observed (outcomes other than `start_lost`,
`no_parent`, `no_parent_cds`).

| Column | Type | Meaning |
| --- | --- | --- |
| `resolved_orf_status` | categorical | `intact` (resolved stop matches the parent stop), `ptc_premature` (a premature stop from a frameshift or exon skip), `ptc_intron_retained` (premature stop and a retained CDS intron), `cds_extension` (read-through past the parent stop, a stop-loss), `no_stop_in_read` (translation runs off the 3' end), `resolution_failed` (no usable start anchor). |
| `resolved_stop_pos` | int / null | Genomic position of the last coding base of the resolved CDS; null when no in-frame stop was found or resolution failed. |
| `resolved_cds_bp` | int | Resolved CDS length in bp (0 when no stop found). |
| `resolved_aa_length` | int | Resolved protein length in amino acids. |
| `resolved_cds_intervals` | list | Genomic intervals of the resolved CDS (JSON-encoded in TSV). |
| `ptc_introduced` | bool | The resolved stop is upstream of the parent stop (a premature termination codon). |
| `intron_retained_in_cds` | bool | A parent intron inside the CDS span is carried as exonic sequence by the isoform. Detected by engulfment, so it is not confused with an exon skip. |
| `frame_consistent` | bool | The resolved stop coincides with the parent stop and no intron is retained. |
| `stop_in_transcript` | bool | An in-frame stop was found before the transcript end. |

## ORF: de novo prediction

Used only for isoforms with no usable parent (`no_parent`, `no_parent_cds`,
`start_lost`). Longest ATG-initiated ORF from orfipy on the spliced sequence.

| Column | Type | Meaning |
| --- | --- | --- |
| `denovo_orf_found` | bool | A de novo ORF above the minimum length was found. |
| `denovo_cds_bp` | int | De novo CDS length in bp. |
| `denovo_orf_aa_length` | int | De novo protein length in amino acids. |
| `denovo_start_codon` | str | Start codon reported by orfipy. |
| `denovo_stop_codon` | str | Stop codon reported by orfipy. |
| `denovo_cds_intervals` | list | Genomic intervals of the de novo ORF. |

## ORF confidence

| Column | Type | Meaning |
| --- | --- | --- |
| `orf_confidence` | categorical | `high`, `medium`, `low`, or `none`. Combines the propagation outcome and the completeness penalty. |
| `orf_confidence_score` | float | Numeric score in [0, 1] behind the category. Thresholds: `--orf-high-confidence` (0.85), `--orf-medium-confidence` (0.5). |

## NMD: geometric and resolved

The escape-rule cascade (stop in last exon, within 50 nt of the last junction,
start-proximal, long last exon, else NMD-sensitive) is applied twice: once to
the geometric stop, once to the resolved stop.

### Interpreting the NMD columns

Every NMD call answers one question: *will nonsense-mediated decay degrade this
transcript?* NMD targets mRNAs whose stop codon looks premature, the classic
trigger being a stop more than ~50 nt upstream of the last exon-exon junction
(an exon-junction complex remains downstream of the stop and recruits the decay
machinery). The status takes three values:

- **`sensitive`** = predicted NMD substrate: the stop is >50 nt upstream of the
  last junction and no escape rule fires, so the transcript is predicted to be
  degraded. Read it as "likely unproductive, little-to-no protein" — typically a
  frameshift, exon skip, retained intron, or a regulated AS-NMD isoform.
- **`escaped`** = has a stop but predicted to evade NMD, because one escape rule
  holds: stop in the last exon (the normal case), within ~50 nt of the last
  junction, a very short CDS (re-initiation), or a very long last exon. **Escaped
  does not mean full-length or normal** — a 5'-truncated isoform whose stop lands
  in the last exon is "escaped" too. It only means "not an NMD target."
- **`not_applicable`** = NMD could not be evaluated: no parent CDS, the start
  codon is not observed, or no in-frame stop was found in the read. **It is the
  absence of a call, not "safe."**

`nmd_rule` / `nmd_rule_resolved` name which rule decided the call (one of the four
escapes, or `ptc_50nt_rule` when sensitive).

**Geometric vs resolved vs de novo** is the same question with three ways of
locating the stop:

- **`nmd_status` (geometric)** uses the stop from projecting the parent CDS
  coordinates onto the isoform (no sequence read). Fast and correct when the CDS
  matches the reference, but the projected end is not the real stop for a
  frameshifted / exon-skipped / intron-retaining isoform.
- **`nmd_status_resolved`** uses the real in-frame stop from translating the
  isoform's own spliced CDS. It catches the premature stops the geometric
  projection misses. **Prefer this column** for parent-anchored isoforms.
- **`nmd_status_denovo`** applies the same rules to the de-novo ORF, for orphan
  isoforms with no reference CDS (always `low` confidence).

The geometric and resolved calls agree for structurally intact isoforms and
diverge exactly where the reading frame departs from the reference, which is
where NMD matters most.

| Column | Type | Meaning |
| --- | --- | --- |
| `nmd_status` | categorical | Geometric NMD call: `sensitive`, `escaped`, `not_applicable`. |
| `nmd_rule` | str | Which escape rule fired (or `ptc_50nt_rule` when sensitive). |
| `stop_to_last_junction_nt` | int / null | mRNA distance from the geometric stop to the last exon-exon junction. |
| `last_exon_length_nt` | int / null | Last exon length. |
| `nmd_confidence` | categorical | `high` for intact propagation, `medium` for disrupted, `none` otherwise. |
| `nmd_status_resolved` | categorical | NMD call on the **resolved** stop. Prefer this for parent-anchored isoforms. |
| `nmd_rule_resolved` | str | Escape rule that fired on the resolved stop. |
| `nmd_confidence_resolved` | categorical | `high` when the resolved ORF is intact, `medium` otherwise, `none` when not applicable. |
| `nmd_status_denovo` | categorical | NMD call on the **de novo** ORF, for orphan isoforms (`no_parent` / `no_parent_cds` / `start_lost`) that have no reference-anchored stop. `not_applicable` for parent-anchored isoforms and orphans with no de novo ORF. |
| `nmd_rule_denovo` | str | Escape rule that fired on the de novo stop. |
| `nmd_confidence_denovo` | categorical | Always `low` when a call is made: the stop comes from a predicted ORF, not a curated reference. |

For a single NMD call per isoform, coalesce the resolved call (reference-anchored,
higher confidence) with the de novo call (orphans):

```python
df["nmd_call"] = df["nmd_status_resolved"].where(
    df["nmd_status_resolved"] != "not_applicable", df["nmd_status_denovo"]
)
```

Advisory NMD branches (separate from the primary call, noisier by nature):

| Column | Type | Meaning |
| --- | --- | --- |
| `uorf_count` | int | Upstream ORFs fully contained in the 5'UTR (ATG to in-frame stop before the main start). |
| `uorf_triggers_nmd` | bool | A uORF stop sits more than `--ptc-threshold-nt` upstream of the transcript's last junction (uORF-triggered NMD heuristic). |
| `long_utr3_triggers_nmd` | bool | The resolved 3'UTR is longer than `--long-utr3-nt` (default 1000). |

## 3'UTR and 5'UTR

| Column | Type | Meaning |
| --- | --- | --- |
| `iso_utr3_length_nt` | int / null | Isoform 3'UTR length from the geometric stop. |
| `parent_utr3_length_nt` | int / null | Parent 3'UTR length. |
| `utr3_length_delta_nt` | int / null | Isoform minus parent 3'UTR (geometric). |
| `utr3_length_delta_pct` | float / null | Same delta as a percentage of the parent. |
| `iso_utr3_length_resolved_nt` | int / null | Isoform 3'UTR length from the **resolved** stop. |
| `utr3_length_delta_resolved_nt` | int / null | Resolved 3'UTR delta vs parent. |
| `utr3_length_delta_pct_resolved` | float / null | Resolved 3'UTR delta percentage. |
| `iso_utr5_length_nt` | int / null | Isoform 5'UTR length (upstream of the start codon). Null when the start is not observed. |
| `parent_utr5_length_nt` | int / null | Parent 5'UTR length. |
| `utr5_length_delta_nt` | int / null | Isoform minus parent 5'UTR. |
| `utr5_length_delta_pct` | float / null | Same delta as a percentage of the parent. |

## Poly(A)

| Column | Type | Meaning |
| --- | --- | --- |
| `polya_signal_motif` | str | Strongest canonical poly(A) signal in the isoform 3'UTR; empty if none. |
| `polya_signal_distance_nt` | int / null | nt from the motif to the 3' end. |
| `polya_evidence_source` | str | `polya_db` (atlas hit), `canonical_motif` (motif fallback), or `none`. |
| `polya_db_site_id` | str | Atlas site ID when the source is `polya_db`. |

## Pfam domains

Populated only with `--pfam-hmm`. The isoform protein is taken from the
resolved CDS when available (so frameshift- and intron-retention-truncated
proteins are scored correctly), falling back to the propagated then de novo CDS.

| Column | Type | Meaning |
| --- | --- | --- |
| `iso_pfam_domains` | list | Pfam domains found in the isoform protein. |
| `parent_pfam_domains` | list | Pfam domains in the parent protein. |
| `pfam_preserved` | list | Domains present in both. |
| `pfam_lost` | list | Parent domains absent from the isoform. |
| `pfam_gained` | list | Isoform domains absent from the parent. |

## Coding potential

A coding-potential score self-calibrated to the supplied reference. CRAFT trains
a model from the reference's own transcripts (CDS-bearing as coding, CDS-less as
non-coding): a hexamer coding/non-coding log-likelihood table plus a logistic
regression on four features (hexamer log-likelihood ratio, log10 ORF length, ORF
coverage, and the Fickett TESTCODE statistic). It then scores each isoform's best
ORF (resolved, else propagated, else de novo). No model file is shipped and no
external tool is required; the model fits whatever organism the reference
describes. The fitted model and a 5-fold cross-validated AUC (about 0.86 on
GENCODE v45) are written to `coding_potential_model.json`. Disable with
`--no-coding-potential`; skipped automatically if the reference has no non-coding
transcripts (columns left empty). This is a screening score; confirm borderline
calls with CPC2 or CPAT.

| Column | Type | Meaning |
| --- | --- | --- |
| `coding_potential_score` | float / null | Logistic probability the ORF is coding, in [0, 1]. Null when there is no ORF. |
| `coding_potential_label` | categorical | `coding` if score ≥ 0.5, else `noncoding`. |
| `coding_potential_orf_source` | str | Which ORF was scored: `resolved`, `propagated`, `denovo`, or `none`. |

Use it to gate the orphan tail: a de-novo ORF with `coding_potential_label = coding`
is a credible novel coding isoform (and its `nmd_status_denovo` call is meaningful),
while `noncoding` flags likely lncRNA or spurious ORFs.

```python
# credible novel coding isoforms among orphans:
df[(df["coding_potential_orf_source"] == "denovo") & (df["coding_potential_label"] == "coding")]
# lncRNA candidates: best overlap is non-coding and the ORF scores non-coding:
df[(df["orf_outcome"] == "no_parent_cds") & (df["coding_potential_label"] == "noncoding")]
```

## External classification passthrough

CRAFT does not do structural QC; it assumes the isoform GTF is already curated by
SQANTI3/pigeon. To combine that upstream classification with CRAFT's consequence
calls, pass the classification table with `--classification FILE`. CRAFT joins the
columns named in `--classification-columns` (default `structural_category`) onto
the per-isoform output by transcript id and appends them as new columns.

The table is any TSV/CSV keyed by isoform id (`isoform`, `transcript_id`, `pbid`,
or the first column; SQANTI3's `*_classification.txt` and pigeon both use
`isoform`). Isoforms absent from the table get an empty value; a carried column
whose name collides with a CRAFT column is prefixed `class_`. CRAFT logs the
match rate to stderr.

This is what makes the "novel splice boundary x functional consequence" analysis a
one-liner: CRAFT supplies the consequence half (`resolved_orf_status`,
`ptc_introduced`, `nmd_status_resolved`, `nmd_status_denovo`,
`coding_potential_label`) and the passthrough supplies the SQANTI structural class.

```python
# NNC isoforms that are NMD substrates with a credible ORF:
nnc = df[df["structural_category"] == "novel_not_in_catalog"]
nnc[(nnc["nmd_status_resolved"] == "sensitive") & (nnc["coding_potential_label"] == "coding")]
```

## Per-cell-type consequence aggregation (v1.5)

With `--counts` and `--group-by OBS_COLUMN`, CRAFT writes
`per_celltype_consequence.tsv`: for each cell group, the molecule-weighted
fraction of detected isoform molecules carrying each consequence. The fraction
for a group `g` and class `c` is

```
sum of molecules in g over isoforms with class c   /   sum of all molecules in g
```

so a highly expressed isoform contributes in proportion to its read support.
The same table is stored in `annotated.h5ad` under `uns['celltype_consequences']`.

| Column | Meaning |
| --- | --- |
| `cell_group` | A value of the `--group-by` obs column. |
| `n_cells` | Cells in the group. |
| `total_molecules` | Summed counts over all isoforms in the group (the denominator). |
| `n_isoforms` | Isoforms with non-zero counts in the group. |
| `frac_nmd_sensitive_resolved` | Fraction of molecules from NMD-sensitive (resolved) isoforms. |
| `frac_ptc_introduced` | Fraction from isoforms with a premature stop. |
| `frac_intron_retained_in_cds` | Fraction from isoforms retaining a CDS intron. |
| `frac_truncated_5p` / `frac_truncated_3p` / `frac_truncated_both` | Fraction from truncated isoforms. |
| `frac_internal_fragment` | Fraction from internal fragments. |
| `frac_alt_3prime_end` | Fraction from alternative-3'-end isoforms. |
| `frac_domain_lost` | Fraction from isoforms with at least one lost Pfam domain. |

Isoforms present in the counts but absent from the per-isoform table still count
toward `total_molecules`; they never contribute to a numerator. A group with no
counts yields `NaN` fractions.

## Command-line options

```
craft annotate --isoforms ISO.gtf --reference REF.gtf --genome GENOME.fa --output-dir OUT/
```

| Option | Default | Effect |
| --- | --- | --- |
| `--counts` | none | Per-cell counts (`.h5ad` or 10x MTX dir); populates `annotated.h5ad`. |
| `--pfam-hmm` | none | Enables Pfam domain analysis. |
| `--polya-atlas` | none | Curated poly(A) BED; drives the `alt_3prime_end` reclassification. |
| `--group-by` | none | Obs column to aggregate consequences by; writes `per_celltype_consequence.tsv` (requires `--counts`). |
| `--coding-potential` / `--no-coding-potential` | on | Score each ORF for coding potential against a reference-calibrated model. |
| `--classification` | none | SQANTI3/pigeon (or any) classification TSV; selected columns are joined onto the output by transcript_id. |
| `--classification-columns` | `structural_category` | Comma-separated columns to carry from `--classification`. |
| `--tolerance` | 50 | End slack (bp) before calling a truncation. |
| `--ptc-threshold-nt` | 50 | PTC rule distance to the last junction. |
| `--start-proximal-nt` | 150 | CDS below this (bp) escapes NMD. |
| `--long-last-exon-nt` | 400 | Last exon above this (bp) escapes NMD. |
| `--min-orf-aa` | 50 | Minimum de novo ORF length. |
| `--orf-high-confidence` | 0.85 | High ORF-confidence cutoff. |
| `--orf-medium-confidence` | 0.5 | Medium ORF-confidence cutoff. |
| `--long-utr3-nt` | 1000 | Long-3'UTR NMD flag threshold. |
| `--prefer-coding-parent` | off | Break parent-selection ties toward CDS-bearing transcripts (off keeps results reproducible). |

## Filter recipes

```python
import pandas as pd
df = pd.read_csv("out/per_isoform.tsv", sep="\t")

# High-confidence NMD substrates, using the sequence-resolved call
df[(df["nmd_status_resolved"] == "sensitive") & (df["nmd_confidence_resolved"] == "high")]

# Isoforms where the resolved engine disagrees with the geometric call
df[df["nmd_status"] != df["nmd_status_resolved"]]

# Premature stops introduced by retained CDS introns
df[df["intron_retained_in_cds"]]

# Domain loss driven by a real premature stop (needs --pfam-hmm)
import json
lost = df["pfam_lost"].apply(lambda s: len(json.loads(s)) > 0)
df[lost & df["ptc_introduced"]]

# Alternative 3' ends with a longer UTR (potential APA with regulatory impact)
df[(df["completeness"] == "alt_3prime_end") & (df["utr3_length_delta_resolved_nt"] > 0)]
```
