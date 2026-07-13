# CRAFT v2 output reference

`per_isoform.tsv` and `per_isoform.json` contain the same 124-field v2 schema.
TSV list fields are JSON encoded. Optional analyses leave their fields empty; they
do not remove columns. Additional classification passthrough columns may be appended.

Coordinates follow the internal 0-based, half-open interval convention. Exact
single-base genomic positions are 0-based. Never infer semantic meaning from column
position; select by name.

## Structure and parent assignment

| Field | Meaning |
| --- | --- |
| `transcript_id` | Input isoform identifier. |
| `completeness` | End structure relative to the selected parent: `full_length`, `truncated_5p`, `truncated_3p`, `truncated_both`, `internal_fragment`, `novel_no_match`, or `alt_3prime_end`. |
| `parent_tx_id`, `parent_gene_id`, `parent_gene_name` | Selected reference transcript and gene identifiers/name. |
| `shared_junctions`, `parent_overlap_bp` | Exact shared splice junctions and stranded exon-overlap bases. |
| `has_cds_bearing_parent` | Whether any plausible overlapping reference parent has CDS annotation. |
| `parent_candidate_count` | Number of ranked reference candidates after gene-consistency restriction. |
| `parent_ambiguous` | Best-versus-second score margin is too small for a unique call. |
| `parent_match_score`, `parent_match_margin` | Composite ranking score and difference from runner-up. Scores rank candidates within an isoform; they are not probabilities. |
| `parent_selection_reason` | Semicolon-delimited evidence used for the winner. |
| `junction_precision`, `junction_recall`, `junction_f1` | Shared-junction agreement relative to isoform/reference chains. |
| `exact_intron_chain`, `iso_chain_contained` | Exact chain equality or all isoform junctions contained in the parent. |
| `parent_reference_priority`, `reference_priority_reason` | Curated tag-derived priority and human-readable reasons. |
| `reference_cds_complete` | Reference has no incomplete-CDS tag and has adequate CDS boundary evidence. |
| `reference_has_explicit_start`, `reference_has_explicit_stop` | Explicit GTF codon features are present. |
| `reference_cds_phase_valid` | CDS phases are internally consistent. |

## Geometric CDS propagation

| Field | Meaning |
| --- | --- |
| `orf_outcome` | Structural projection result such as `propagated_intact`, `disrupted`, `start_lost`, `stop_not_observed`, `no_parent`, or `no_parent_cds`. |
| `propagated_cds_bp`, `parent_cds_bp` | Projected and reference CDS bases. |
| `start_codon_covered`, `stop_codon_covered` | Reference boundaries occur on the observed isoform structure. |
| `propagated_cds_intervals` | Projected genomic CDS intervals. |

## Sequence-resolved and censored ORF

| Field | Meaning |
| --- | --- |
| `resolved_orf_status` | `intact`, `ptc_premature`, `ptc_intron_retained`, `cds_extension`, `left_censored`, `right_censored`, opt-in `start_rescued`, legacy `no_stop_in_read`, or `resolution_failed`. |
| `resolved_stop_pos` | Legacy genomic position of the last sense-CDS base. Prefer `resolved_stop_codon_pos` for termination geometry. |
| `resolved_start_pos` | Exact genomic position of the supported/inferred start-codon anchor. |
| `resolved_stop_codon_pos` | Exact genomic anchor for the first base of the complete stop codon in transcript direction. |
| `resolved_cds_bp`, `resolved_aa_length`, `resolved_cds_intervals` | Complete resolved coding length, amino-acid length, and genomic intervals. |
| `ptc_introduced` | Observed stop is premature relative to the parent. |
| `intron_retained_in_cds` | A parent CDS intron is engulfed by an isoform exon. |
| `frame_consistent` | Resolved sequence preserves the expected frame. |
| `stop_in_transcript` | A complete in-frame stop codon was observed. |
| `uorf_count`, `uorf_triggers_nmd` | Upstream ORF candidates and structural uORF-surveillance heuristic. |
| `orf_start_observed`, `orf_stop_observed` | Whether the molecule contains the supported start/stop boundary. |
| `orf_censoring` | `none`, `left`, `right`, or `both`. |
| `partial_cds_bp`, `partial_cds_intervals` | Observed coding-frame segment when a complete ORF cannot be established. |
| `alternative_start_inferred` | A downstream in-frame ATG was selected under explicit opt-in. |

## De novo ORF and confidence

| Field | Meaning |
| --- | --- |
| `denovo_orf_found` | A candidate meeting the configured minimum length was found. |
| `denovo_cds_bp`, `denovo_orf_aa_length` | Candidate length in bases and amino acids. |
| `denovo_start_codon`, `denovo_stop_codon` | Candidate codons. |
| `denovo_cds_intervals` | Candidate genomic intervals. |
| `orf_confidence`, `orf_confidence_score` | Rule-based confidence category and component score; neither is a calibrated correctness probability. |

## RNA-surveillance prediction

| Field | Meaning |
| --- | --- |
| `nmd_status`, `nmd_rule`, `nmd_confidence`, `nmd_basis` | Legacy compatibility fields: structural `sensitive`/`escaped`/`not_applicable`, matched rule, confidence, and resolved/de novo basis. |
| `stop_to_last_junction_nt` | Spliced bases after the complete stop codon to the last exon junction; zero in the last exon. |
| `last_exon_length_nt`, `ptc_exon_length_nt` | Terminal-exon length and length of the stop-containing exon. |
| `nmd_susceptibility` | `likely_sensitive`, `likely_escape`, or `indeterminate`. This is not observed decay. |
| `nmd_rule_score` | Transparent rule severity on 0–1 scale, not a probability. |
| `nmd_evidence_tier` | Strength of structural evidence: `strong`, `moderate`, `limited`, or `none`. |
| `surveillance_status`, `surveillance_mechanism` | Generalized susceptibility status and mechanism (`nmd`, `nonstop_decay`, or `none`). |
| `nonstop_decay_candidate` | Stopless/right-censored structure is compatible with non-stop decay. |
| `surveillance_limitations` | Reason the rule result should be down-weighted or is indeterminate. |
| `long_utr3_triggers_nmd` | Long-3'UTR heuristic relative to the configured threshold. |

## UTR and poly(A)

| Field | Meaning |
| --- | --- |
| `iso_utr3_length_nt`, `parent_utr3_length_nt` | Isoform/reference 3' UTR length after the full three-base stop codon. |
| `utr3_length_delta_nt`, `utr3_length_delta_pct` | Isoform minus parent 3' UTR difference. |
| `iso_utr5_length_nt`, `parent_utr5_length_nt` | Isoform/reference sequence upstream of the supported start. Isoform is empty when start is censored. |
| `utr5_length_delta_nt`, `utr5_length_delta_pct` | Isoform minus parent 5' UTR difference. |
| `polya_signal_motif`, `polya_signal_distance_nt` | Best canonical motif and distance to the observed 3' end. |
| `polya_evidence_source`, `polya_db_site_id` | Curated atlas, motif, or no evidence; matching atlas site ID. |

## Protein domains and coding potential

| Field | Meaning |
| --- | --- |
| `iso_pfam_domains`, `parent_pfam_domains` | Domain calls for isoform and parent proteins. |
| `pfam_preserved`, `pfam_lost`, `pfam_gained` | Domain-set comparison. Empty without `--pfam-hmm`. |
| `coding_potential_score` | Reference-trained logistic classifier score; not probability-calibrated. |
| `coding_potential_label` | Thresholded `coding`/`noncoding` screening label. |
| `coding_potential_orf_source` | `resolved`, `propagated`, `denovo`, or `none`. |

## Counts and recurrence

| Field | Meaning |
| --- | --- |
| `total_count` | UMI-corrected molecules across included cells. |
| `n_cells_detected`, `n_cells_total` | Cells with at least one molecule and all included cells. |
| `detection_fraction` | `n_cells_detected / n_cells_total`. |
| `molecules_per_detected_cell` | `total_count / n_cells_detected`. |
| `isoform_fraction_within_gene` | Isoform count divided by counts for isoforms assigned to the same parent gene. |
| `recurrence_pvalue`, `recurrence_score` | Optional null-model dispersion statistic and `1-p`. Exploratory; not transcript-realness probabilities. |

## Independent molecule/read evidence

Canonical evidence fractions are 0–1. The loader also accepts documented aliases and
0–100 percentages.

| Field | Meaning |
| --- | --- |
| `isoform_evidence_score` | Weighted mean of available favorable evidence and reversed artifact evidence. Uncalibrated. |
| `isoform_evidence_tier` | `strong`, `supported`, `limited`, `artifact_likely`, or `insufficient_evidence`. |
| `evidence_feature_count` | Number of measured components entering the score. |
| `evidence_model` | Model/version identifier. |
| `evidence_warnings` | JSON list, including insufficient features or strong artifact signal. |
| `unique_molecule_fraction`, `ambiguous_molecule_fraction` | Uniquely/ambiguously assigned support among molecules or reads. |
| `canonical_junction_fraction` | Fraction of relevant splice junctions with canonical motifs. |
| `short_read_junction_support_fraction` | Fraction supported by orthogonal short reads. |
| `full_length_fraction` | Molecules classified full length by the upstream assay. |
| `five_prime_adapter_fraction` | 5' adapter/TSO or equivalent completion evidence. |
| `polya_tail_fraction` | Molecules with direct poly(A)-tail evidence. |
| `mapq_fraction` | Mapping-quality summary scaled by 60 and clipped to 0–1. |
| `replicate_support_fraction` | Replicate/sample count scaled to saturation at two. |
| `internal_priming_fraction` | Evidence of oligo(dT) internal priming. |
| `template_switch_fraction` | Reverse-transcription template switching/strand invasion evidence. |
| `chimera_fraction` | Chimeric molecule/alignment evidence. |

## Independent ORF comparison

| Field | Meaning |
| --- | --- |
| `comparator_orf_present` | Comparator GTF contains CDS for this transcript. |
| `comparator_start_pos`, `comparator_stop_codon_pos`, `comparator_cds_bp` | Comparator start, stop anchor, and CDS length. |
| `comparator_start_agrees`, `comparator_stop_agrees` | Exact positional agreement with CRAFT. |
| `comparator_cds_bp_delta` | CRAFT resolved CDS bases minus comparator CDS bases. |

## Interpretation rules

- Empty optional fields mean missing evidence, not evidence of absence.
- `parent_ambiguous` should lower confidence in every reference-propagated conclusion.
- Censored ORFs must not be compared with complete ORFs as if their boundaries were observed.
- Evidence tiers and scores require assay-specific calibration before hard filtering.
- NMD fields report structural susceptibility; actual decay requires orthogonal validation.
