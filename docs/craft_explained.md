# How CRAFT v2 works

This is the implementation-oriented companion to [`v2_design.md`](v2_design.md).
The source, rather than pinned historical line links, is authoritative during v2
development.

## Pipeline

`craft annotate` calls `pipeline.run_annotate` and performs these stages:

1. Load isoform exons and reference exon/CDS/start/stop features.
2. Build reference CDS-quality and curated-tag metadata.
3. Import unambiguous parent/gene hints from a classification table.
4. Rank candidate parents and classify transcript-end completeness.
5. Project the parent CDS geometrically.
6. Reconstruct spliced sequence and resolve the ORF in the supported frame.
7. Mark unobserved starts/stops as censoring; optionally test an alternative start.
8. Predict a de novo candidate for orphans/no-CDS parents.
9. Compute ORF-confidence features, surveillance susceptibility, and UTR features.
10. Add poly(A), Pfam, coding-potential, count/recurrence, molecule-evidence, and
    independent-ORF comparison blocks when their inputs are available.
11. Write TSV, JSON, HTML, AnnData, and optional group summaries.

CRAFT preserves every input isoform. Missing optional inputs produce empty or
`insufficient_evidence` fields rather than implicit negative evidence.

## Parent ranking

`core/completeness.py` builds overlapping reference candidates on the same strand.
It prioritizes unique upstream hints, gene consistency, exact/contained intron
chains, junction precision/recall/F1, normalized exonic overlap, and curated
reference priority. It reports the candidate count and best-versus-second margin.
Close alternatives set `parent_ambiguous`; downstream evidence is correspondingly
weaker.

End completeness is defined relative to that parent with a configurable genomic
tolerance. Alternative 3' ends can be supported by a curated poly(A) atlas or motif
evidence. Structural completeness is not the same as molecular full-length status.

## ORF projection and resolution

`core/orf/propagation.py` intersects the parent CDS with the isoform and records
start/stop coverage. `core/orf/resolve.py` reconstructs strand-correct spliced
sequence, translates in the supported frame, and maps exact transcript positions
back to the genome.

The resolver distinguishes:

- `intact`: observed stop agrees with the parent;
- `ptc_premature` / `ptc_intron_retained`: observed earlier stop;
- `cds_extension`: observed later stop;
- `left_censored`: the annotated start is outside the molecule;
- `right_censored`: no complete stop is observed before the molecule ends;
- `start_rescued`: explicit opt-in alternative-start hypothesis;
- `resolution_failed`: no defensible sequence placement.

`resolved_stop_pos` remains a legacy last-sense position. New logic should use
`resolved_stop_codon_pos`. The latter anchors the three-base stop and correct UTR
calculation.

## De novo ORFs and coding potential

De novo search is primarily for no-parent/no-CDS cases. Its start and stop are
candidate boundaries, not proof of translation. Coding potential combines hexamer,
ORF length/coverage, and Fickett features. Cross-validation rebuilds all
sequence-derived features and scaling inside each fold. The output remains a
classifier score because no probability calibration is performed.

`--orf-comparator-gtf` adds exact start/stop/CDS-length agreement with another
GTF-native caller. Agreement supports a call; disagreement exposes uncertainty.

## Surveillance susceptibility

`core/nmd.py` applies an ordered structural rule cascade. It records the geometry,
rule, qualitative susceptibility, rule severity, evidence tier, and limitations.
Censored ORFs, incomplete reference CDS, or ambiguous parents cannot receive the
same evidentiary interpretation as an observed complete CDS. Right-censored
stopless molecules can also be labelled as non-stop-decay candidates.

The legacy `nmd_status` remains for compatibility. New analyses should use
`nmd_susceptibility`, `nmd_evidence_tier`, and `surveillance_limitations` together.

## Isoform evidence

`core/evidence.py` accepts common aliases from per-isoform evidence tables. It
normalizes fractions/percentages, retains every component, and combines available
features into a weighted 0–1 score. Strong artifact evidence overrides a favorable
average. Fewer than three measured components gives `insufficient_evidence`.

This model is deliberately transparent and untrained. To replace it with a learned
FLIGHT-seq model, retain sample- and gene-level grouping during nested validation,
calibrate held-out probabilities, and publish the exact labels and error strata.

## Recurrence and groups

With counts, CRAFT reports molecules, cells detected, total cells, detection
fraction, molecules per detected cell, and within-gene fraction. Optional recurrence
null statistics describe dispersion relative to a chosen null. They are exploratory
statistics and are not probabilities that a transcript is genuine.

Group summaries are molecule-weighted descriptions. AS–NMD candidate listings use
strong/supported independent evidence when it is available; expression recurrence
alone is not promoted to structural validation.

## Extension points

The next high-value additions are calibrated FLIGHT-seq labels, technical-replicate
and sample-aware evidence, direct scNoiseMeter/tecap adapters, Ribo-seq/TIS and
protein-homology evidence, splice-site conservation, and benchmark harnesses for
SQANTI3/Pigeon, IsoQuant, Bambu, Isosceles, ORFanage, CPAT/CPC2, and NMDetective.
