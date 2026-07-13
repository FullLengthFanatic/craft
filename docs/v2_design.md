# CRAFT v2 scientific design

This document defines what CRAFT v2 claims, what it computes, and what still needs
experimental calibration. It supersedes the v1 narrative documents when their
terminology conflicts.

## 1. Four separate inference problems

Long-read transcript analysis often collapses four different questions into one
"high-confidence isoform" label. CRAFT keeps them separate:

| Question | Relevant evidence | v2 output |
| --- | --- | --- |
| Is the transcript structure plausible? | unique molecules, splice motifs, junction support, ends, replicates, artifact signals | evidence features, score, tier, warnings |
| Is the molecule complete at each end? | adapters, poly(A) tail/site, reference-coordinate coverage | completeness, observed boundaries, censoring |
| Which ORF is supported? | parent CDS, phase, observed sequence, independent ORF caller | propagated/resolved/de novo fields and caller comparison |
| Could RNA surveillance act? | complete stop, exon-junction geometry, CDS/exon context | susceptibility, rule, evidence tier, limitations |

A strong answer to one question does not imply a strong answer to the others.

## 2. Isoform evidence, not recurrence as truth

Cell recurrence and molecule abundance are useful descriptive summaries, but both
depend on expression, capture efficiency, cell composition, sequencing depth, and
the quantifier. They cannot by themselves distinguish a rare real transcript from
a recurrent technical product.

With `--evidence-table`, CRAFT accepts per-isoform fractions or counts for:

- unique and ambiguous molecule assignment;
- canonical junctions and short-read junction support;
- full-length, 5' adapter, and poly(A)-tail support;
- mapping quality and replicate/sample support;
- internal priming, template switching/strand invasion, and chimera evidence.

The `transparent_uncalibrated_v1` score is a weighted mean of the available
components, with artifact fractions reversed. It reports the number of observed
features and refuses a substantive tier when fewer than three are available. It is
not a posterior probability, false-discovery rate, or trained FLIGHT-seq classifier.

The weights are inspectable in `core/evidence.py`. A future learned model should
replace them only after nested sample-level validation and probability calibration.
Raw features remain in the output so users can audit every decision.

## 3. Parent assignment is an inference

Reference propagation is only as reliable as the parent assignment. CRAFT ranks
candidates using:

1. explicit upstream transcript and gene hints, when unique;
2. exact or contained intron-chain agreement;
3. junction precision, recall, and F1;
4. exonic overlap normalized to isoform length;
5. curated reference priority (MANE Select/Plus Clinical, APPRIS principal, CCDS,
   GENCODE basic), with incomplete-CDS tags penalized;
6. CDS-bearing preference only as an explicit tie behavior.

The number of candidates, best score, margin, selection reason, and ambiguity flag
are emitted. A reference transcript is not "real by definition"; curated metadata
changes the strength of evidence, not ontology.

## 4. Truncation is censoring

For a 5'-truncated transcript, the true translation start may lie outside the
sequenced molecule. CRAFT therefore reports a `left_censored` partial CDS and does
not silently turn the first downstream in-frame ATG into the biological start.
Users can request that hypothesis with `--infer-alternative-start`; it remains
labelled `alternative_start_inferred`.

For a 3'-truncated molecule with no complete stop codon, CRAFT reports
`right_censored`, retains the observed partial CDS, and can flag a non-stop-decay
candidate separately. The last sense codon is not reported as a stop codon.

The output records exact start and stop-codon positions, whether each is observed,
censoring direction, and both complete and partial CDS intervals.

## 5. Reference CDS quality

CRAFT preserves explicit `start_codon`, `stop_codon`, CDS phase, transcript
biotype, transcript-support level, and reference tags. It records whether a CDS is
complete and phase-consistent. Propagating from an incomplete or phase-invalid CDS
is reported with weaker downstream evidence.

This matters especially for reference ISMs, retained introns, read-through loci,
and annotations whose CDS begins or ends at the transcript boundary.

## 6. ORF inference and comparison

The order of evidence is:

1. project an adequately assigned parent CDS;
2. reconstruct and translate the observed spliced sequence in that frame;
3. represent missing boundaries as censoring;
4. use de novo candidates for true orphans/no-CDS parents;
5. optionally compare exact start, stop, and CDS length with an independent GTF
   such as ORFanage output.

Independent agreement is supporting evidence, not proof. Disagreement is retained
for review rather than resolved by majority vote. Ribosome profiling, proteomics,
TIS profiling, conservation, and curated protein evidence are appropriate future
orthogonal inputs.

## 7. RNA-surveillance language

CRAFT implements a transparent structural rule cascade for predicted NMD
susceptibility, including last-junction distance and recognized escape contexts.
It also separates non-stop-decay candidates. The output deliberately uses:

- `nmd_susceptibility`, not measured NMD;
- `nmd_rule_score`, a rule severity index, not a probability;
- `nmd_evidence_tier`, which is reduced for censored ORFs, ambiguous parents, or
  incomplete reference CDS;
- `surveillance_limitations`, an explicit audit string.

Actual degradation depends on translation, tissue/cell state, exon-junction-complex
deposition, alternative initiation, RNA-binding proteins, and other features not
established by transcript structure alone. Perturbation or decay measurements are
needed for biological validation.

## 8. Coding potential

The coding model uses reference-derived hexamer, ORF length/coverage, and Fickett
features. v2 rebuilds hexamer tables and scaling inside each cross-validation fold,
preventing leakage from validation sequences into training features. Its output is
a classifier score. It is not guaranteed to be calibrated, and a CDS-less reference
transcript is an imperfect negative label.

For serious novel-coding claims, compare against established tools such as CPAT or
CPC2 and add conservation, protein homology, Ribo-seq, or proteomics where possible.

## 9. Calibration plan for FLIGHT-seq

The repository does not contain enough labelled FLIGHT-seq data to justify a
trained "real versus fake" probability. A defensible model should use:

1. positive controls: reference spike-ins, concordant high-quality annotated
   transcripts, and orthogonally confirmed novel junctions/ends;
2. negative controls: permuted junctions, internal-priming products, strand-invasion
   signatures, synthetic chimeras, and empty-droplet/technical controls;
3. splits by biological sample or donor, not random isoform rows;
4. gene-held-out analysis to prevent near-identical isoforms leaking across folds;
5. comparison with SQANTI3/Pigeon, IsoQuant, Bambu, and Isosceles outputs;
6. precision-recall curves, calibration plots, error stratification by abundance,
   novelty class, truncation, and gene complexity;
7. a frozen external test set and explicit uncertainty intervals.

Until then, use evidence tiers for transparent triage and keep the component
features in every downstream model.

## 10. Literature and tool context

The design follows the broad consensus that no single signal establishes a novel
isoform. Relevant mature or widely used references include SQANTI/SQANTI3 and the
PacBio Pigeon workflow for structural QC; IsoQuant and Bambu for long-read
transcript discovery/quantification; Isosceles for long-read scRNA-seq artifact
handling; ORFanage and ORFannotate for reference-aware ORF annotation; CPAT and
CPC2 for coding-potential baselines; and NMDetective plus foundational mammalian
NMD-rule studies for surveillance prediction.

These tools answer overlapping but non-identical questions. CRAFT should be
benchmarked against them on the same molecules and truth proxies, not described as
superior from internal concordance alone.

Selected sources:

- Tardaguila et al., [SQANTI](https://genome.cshlp.org/content/28/3/396),
  *Genome Research* (2018), and Pardo-Palacios et al.,
  [SQANTI3](https://www.nature.com/articles/s41592-024-02229-2),
  *Nature Methods* (2024).
- Prjibelski et al., [IsoQuant](https://www.nature.com/articles/s41587-022-01565-y),
  *Nature Biotechnology* (2023).
- Chen et al., [Bambu](https://www.nature.com/articles/s41592-023-01908-w),
  *Nature Methods* (2023).
- Kabza et al., [Isosceles](https://www.nature.com/articles/s41467-024-51584-3),
  *Nature Communications* (2024).
- Varabyou et al., [ORFanage](https://www.nature.com/articles/s43588-023-00496-1),
  *Nature Computational Science* (2023).
- Wang et al., [CPAT](https://academic.oup.com/nar/article/41/6/e74/2902455),
  *Nucleic Acids Research* (2013), and Kang et al.,
  [CPC2](https://academic.oup.com/nar/article/45/W1/W12/3831091),
  *Nucleic Acids Research* (2017).
- Lindeboom et al., [position-dependent NMD rules](https://www.nature.com/articles/ng.3664),
  *Nature Genetics* (2016), and
  [NMDetective](https://www.nature.com/articles/s41588-019-0517-5),
  *Nature Genetics* (2019).
- Kurosaki et al., *Nature Reviews Molecular Cell Biology* (2019), mammalian NMD.
