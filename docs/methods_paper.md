# CRAFT: truncation-aware ORF propagation for long-read isoform annotation

Simone Picelli
Institute of Molecular and Clinical Ophthalmology Basel (IOB)
simone.picelli@iob.ch

---

## Abstract

Long-read isoform sequencing (PacBio Iso-Seq, ONT direct RNA, MAS-Seq) recovers full-length transcripts from cells where short-read assembly cannot, but the resulting isoform calls are routinely truncated at one or both ends, especially in single-cell protocols where library bias and read-length variation clip the 5' and 3' termini. De novo ORF prediction on these truncated reads picks the longest plausible ATG, which in practice means a wrong-but-nearby start codon for roughly five percent of transcripts on synthetic truncation and fifteen to twenty-eight percent on real PacBio single-cell data.

CRAFT (Coding Region Annotation From Templates) is a Python toolkit that side-steps this by propagating CDS coordinates from a matched reference parent transcript when shared splice junctions warrant, and emitting de novo predictions only for genuinely novel isoforms. It bundles ORF propagation with a rule-based NMD label, polyA-atlas-aware alternative-3'-end detection, Pfam domain disruption analysis, and an AnnData export for downstream single-cell pipelines.

We validate the central claim through three benchmarks. On simulated truncation of 63,332 GENCODE protein-coding transcripts across a four-rate by three-orientation by three-seed grid, CRAFT recovers the start codon at 0.98 to 1.00 versus 0.94 to 0.95 for orfipy de novo. On 223,976 real PacBio single-cell isoforms whose ORF survives the observed truncation, orfipy alone misses the start at 15 to 28 percent depending on completeness category. CRAFT's rule-based NMD-sensitive label is enriched among transcripts upregulated under UPF1 knockdown in a HeLa bulk RNA-seq dataset (odds ratio 1.46, one-sided Fisher's p = 2.4e-15). CRAFT is MIT-licensed and available at https://github.com/FullLengthFanatic/craft.

---

## 1. Introduction

Long-read RNA sequencing has moved from a specialty tool to a routine method for isoform-level transcript characterization. PacBio Iso-Seq and ONT direct RNA both produce reads long enough to span entire transcripts, and single-cell extensions (MAS-Seq [1], scISOr-Seq [2], scNanoSeq [3]) bring isoform resolution to the per-cell scale. The dominant downstream questions are no longer "which transcripts are present" but "what do those transcripts encode": where the open reading frame starts and stops, whether the predicted protein retains its functional domains, whether the transcript is likely to undergo nonsense-mediated decay.

The functional-annotation tooling has not kept pace. IsoformSwitchAnalyzeR [4] is the most complete option but is R-only and built around bulk switch-detection. IsoAnnotLite [5] is Python-native but only transfers annotations to known isoforms (those that exactly match a reference transcript); genuinely novel isoforms get nothing. tappAS [6] is a broader desktop GUI, last meaningfully updated in 2021. None of these tools were designed for the dominant single-cell long-read failure mode: 5' or 3' truncation of an otherwise-recognizable isoform.

CRAFT addresses the gap. It takes a long-read isoform GTF (from FLAIR, IsoQuant, Bambu, FLAMES, or SQANTI3), a reference annotation GTF with CDS records (GENCODE or Ensembl), and a genome FASTA, and emits per-isoform functional annotations with explicit truncation-aware confidence flags. The novelty is the ORF assignment strategy: for each novel isoform CRAFT identifies its best-matching reference transcript by maximal splice-junction sharing, then propagates the reference CDS coordinates onto the isoform through shared exons. When the propagation succeeds and both start and stop codons fall inside the isoform's exons, the propagated ORF is reported with high confidence. When truncation has clipped past the start or the stop, the outcome is labeled accordingly (`start_lost`, `stop_not_observed`, `stop_at_alt_polya`) and confidence is downgraded. Only isoforms with no usable parent fall through to a de novo `orfipy` call.

Three benchmarks back the central claim that propagation beats de novo prediction on long-read isoform data.

---

## 2. Results

### 2.1 Software overview

CRAFT is a Python ≥ 3.10 package with a Click-based CLI. The pipeline runs in a single `craft annotate` invocation. Required inputs are the isoform GTF, reference GTF, and genome FASTA; the optional `--polya-atlas` flag accepts a BED of curated polyadenylation sites for direct-evidence alternative-3'-end detection (PolyASite v3.0 [7] and PolyA_DB v4 [8] are recommended). Outputs are a per-isoform TSV with one row per input isoform, the equivalent JSON for programmatic consumers, an AnnData (`.h5ad`) export with isoforms as `var` and functional annotations as `var` columns, and a self-contained interactive HTML report.

The internals are organized as composable modules under `src/craft/core/`: `completeness` (structural classification against the reference), `orf/propagation` and `orf/denovo` (the two ORF paths), `orf/confidence` (truncation-aware confidence scoring), `nmd` (rule cascade), `pfam` (domain disruption via pyhmmer [9]), `utr3` (3' UTR features and polyA atlas matching), `report` (plotly HTML), and `export` (AnnData/MuData writers).

Two performance optimizations matter at scale. The polyA atlas match is pre-indexed by (chromosome, strand) with sorted PAS midpoints and resolved via `numpy.searchsorted` in O(log n) per isoform, reducing chr22 runtime with a filtered PolyASite v3.0 atlas from approximately three minutes to under thirty seconds. The pipeline filters isoforms on contigs absent from the genome FASTA up front so that the per-isoform sequence-fetch loop cannot abort mid-run on random or alt contigs that PacBio collapse routinely emits. A full-genome run on a 698k-isoform PacBio bcM0003 sample completes in 71 minutes with peak resident memory of approximately 18 GB.

### 2.2 Bench 1: simulated truncation outperforms de novo prediction

To isolate the propagation versus de novo question from confounders, we truncated GENCODE v45 protein-coding transcripts in silico across a four-by-three-by-three factorial grid: truncation rates {5, 10, 25, 50}%, orientations {5', 3', both}, and three random seeds per cell sampling 3,000 transcripts each from a pool of 63,332 protein-coding transcripts with complete CDS (Methods). For each truncated transcript we computed the GENCODE-truth ORF in the truncated transcript's coordinates and required it to be fully intact (both start and stop codons within the truncated sequence); transcripts where truncation clipped the start or stop were excluded.

CRAFT's propagation was run with the original full-length transcript as the parent reference; orfipy was run de novo on the truncated transcript sequence with no parent information. Across 92,718 scored rows in 8.6 minutes wall time, CRAFT's start-codon exact-match rate sits at 0.98 to 1.00 in every cell; orfipy plateaus at 0.94 to 0.95 (Table 1, Figure 1). CRAFT's mean absolute ORF length error is exactly zero nucleotides for every 3'-truncated cell (3'-end truncation does not move the start codon, and propagation inherits the parent's positions); orfipy sits between 8 and 12 nucleotides across all conditions because the alternative ATGs it selects are typically a few codons downstream of the true start. The one cell where orfipy edges out CRAFT (50% 5'-truncation, n=35 intact-truth isoforms) is too small to read into: a single CRAFT failure shifts the rate by 2.8 percentage points.

**Table 1. Bench 1 selected cells.**

| rate | orientation | n     | CRAFT start | orfipy start | CRAFT \|len_err\| | orfipy \|len_err\| |
|------|-------------|-------|-------------|--------------|-------------------|--------------------|
| 5%   | 3'          | 2,701 | 1.000       | 0.947        | 0.0 nt            | 11.0 nt            |
| 10%  | 3'          | 2,445 | 1.000       | 0.944        | 0.0 nt            | 10.9 nt            |
| 25%  | 3'          | 1,748 | 1.000       | 0.945        | 0.0 nt            | 11.5 nt            |
| 50%  | 3'          |   842 | 1.000       | 0.945        | 0.0 nt            | 12.2 nt            |
| 5%   | 5'          | 1,796 | 0.999       | 0.948        | 1.2 nt            | 8.1 nt             |
| 25%  | 5'          |   253 | 0.990       | 0.936        | 5.8 nt            | 7.9 nt             |
| 25%  | both        |   604 | 0.994       | 0.938        | 3.5 nt            | 11.5 nt            |
| 50%  | both        |    84 | 0.985       | 0.937        | 5.7 nt            | 9.4 nt             |

The five-percentage-point start-exact gap is consistent across truncation rates and orientations, which is the cleanest statement of the central methods claim: a de novo predictor with no parent information picks a wrong-but-nearby ATG for one in every twenty truncated transcripts, and propagation closes that gap.

### 2.3 Bench 2: real-data ORF concordance on bcM0003 single-cell PacBio Iso-Seq

The simulated benchmark above uses clean rectangular truncation, which does not reflect real long-read 5'/3' bias. To test propagation on biologically realistic truncation patterns, we used a 698,049-isoform PacBio single-cell Iso-Seq sample (bcM0003) that CRAFT had already annotated against GENCODE v45 with a filtered PolyASite v3.0 atlas. For each isoform with a GENCODE protein-coding parent we projected the parent's start_codon and stop_codon genomic positions onto the isoform's exon structure; isoforms where either codon fell outside the isoform's exons were dropped. 223,976 isoforms remained eligible.

For these isoforms we compared CRAFT's propagated ORF (parsed from `propagated_cds_intervals`) and orfipy's de novo call against the projected truth, stratified by CRAFT's `completeness` category (Table 2, Figure 2). CRAFT scores 1.00 in every category. Orfipy ranges from 0.727 in the `full_length` category to 0.851 in `truncated_3p`; the gap is 15 to 28 percentage points.

**Table 2. Bench 2 start-codon exact match by completeness category.**

| completeness      | n      | CRAFT | orfipy | gap (ppt) |
|-------------------|--------|-------|--------|-----------|
| full_length       | 95,395 | 1.000 | 0.727  | +27.3     |
| alt_3prime_end    | 79,511 | 1.000 | 0.765  | +23.5     |
| truncated_5p      | 23,576 | 1.000 | 0.740  | +26.0     |
| truncated_3p      |  4,410 | 1.000 | 0.851  | +14.9     |
| truncated_both    |  9,506 | 1.000 | 0.790  | +21.0     |
| internal_fragment | 11,578 | 1.000 | 0.812  | +18.8     |

CRAFT's perfect score is partially tautological: the truth is defined as parent-projected GENCODE coordinates, and CRAFT propagates from the same parent. The valid claim is the orfipy bar: on real long-read single-cell isoforms, even when the truth ORF is fully recoverable, de novo prediction picks the wrong start codon 15 to 28 percent of the time. The orfipy gap is larger than in Bench 1 (Bench 1 sits around five percentage points across all conditions) because real PacBio data has more diverse 5'/3' end variability than clean rectangular truncation, including alternative transcription start sites, alternative polyadenylation, and intronic noise that creates additional plausible ATG candidates downstream of the true start.

### 2.4 Bench 3: CRAFT NMD-sensitive labels track UPF1-KD response

The third benchmark tests whether CRAFT's rule-based NMD label correlates with the biochemical response when the NMD machinery is disabled. We took the GSE86148 dataset (HeLa cells, Lykke-Andersen lab, SRP083135 [10]) and quantified the three scrambled-control and three UPF1-knockdown bulk RNA-seq samples against the GENCODE v45 transcriptome with salmon [11]. Per-sample mapping rates were 83.9 to 85.9 percent with no batch outliers. Differential expression was tested with pydeseq2 [12] using the Wald test on the `~condition` design.

In parallel, we built a CRAFT NMD universe by running CRAFT on the 80,441 GENCODE v45 transcripts that are protein-coding or nonsense_mediated_decay-annotated and have complete CDS, with the full GENCODE GTF as reference. CRAFT labeled 7,282 transcripts (9.1%) as NMD-sensitive, 72,906 as escaped, and 253 as not_applicable. We merged the DE results against this universe, restricted to transcripts whose CRAFT outcome is `propagated_intact` or `disrupted` with a clean stop codon (47,378 eligible transcripts), and tested for enrichment of CRAFT NMD-sensitive labels among UPF1-KD-upregulated transcripts (log2 fold-change ≥ 1, adjusted p < 0.05).

**Table 3. Bench 3 contingency table.**

|              | upregulated | not upregulated |
|--------------|-------------|-----------------|
| NMD-sensitive |   620      |  3,228          |
| NMD-escaped   | 5,070      | 38,460          |

The one-sided Fisher's exact test gives odds ratio 1.457 and p = 2.4e-15. CRAFT-labeled NMD-sensitive transcripts are upregulated under UPF1 KD at 16.1 percent versus 11.6 percent for NMD-escaped transcripts. The log2 fold-change CDF for NMD-sensitive transcripts is visibly right-shifted versus NMD-escaped (Figure 3), which is the directionality NMD biology predicts when UPF1 is silenced.

The modest effect size is expected. CRAFT's NMD label is structural (rule cascade over 50nt distance to last junction, last exon length, start-proximal escape) and does not capture the full set of biochemical determinants that govern UPF1 substrate specificity, particularly in 48-hour siRNA experiments where secondary effects compound. The benchmark answers the directionality question (do CRAFT's labels track real NMD biology) affirmatively rather than claiming the rule cascade is a perfect biochemical predictor.

---

## 3. Methods

### 3.1 ORF propagation

The propagation algorithm operates per isoform in three steps. First, the isoform's exon set is matched against reference transcript exons via maximal splice-junction sharing using pyranges [13] interval operations; ties are broken by exon-coverage fraction. The matched reference transcript becomes the isoform's `parent_tx_id`. Second, the parent's CDS records are walked in genomic order from the start codon. As long as each exon-exon junction in the parent's CDS region is also a junction in the isoform, CRAFT projects the parent's CDS exon onto the isoform with the same genomic coordinates. At the first junction divergence, CRAFT records the structural change (alternative splice site, exon skip, intron retention, alternative transcription start site) and continues propagation only if the reading frame is preserved. Third, the outcome is labeled. `propagated_intact` means both start and stop codons are inside the isoform's exons. `start_lost` means the start codon is outside the isoform (typical for 5'-truncated isoforms with short 5' UTRs). `stop_not_observed` means the start is present but the stop is past the isoform's 3' end. `stop_at_alt_polya` means the isoform terminates at an alternative polyadenylation site upstream of the parent's stop codon, evidenced either by a canonical poly(A) signal motif within fifty nucleotides of the isoform's 3' end or by a hit against a user-supplied PAS atlas BED.

### 3.2 De novo path

When no parent transcript matches the isoform (`novel_no_match` completeness category), the propagation path is skipped and `orfipy` [14] is invoked on the isoform's transcript-orientation sequence with a default minimum ORF length of 75 nucleotides, ATG start codons, and TAA/TAG/TGA stops. The longest ORF on the forward strand is reported. The de novo path is used for less than twelve percent of isoforms in typical PacBio runs (11.3 percent for bcM0003) because most novel isoforms share enough junction structure with at least one reference transcript to support propagation.

### 3.3 Truncation-aware confidence

Each propagated or de novo ORF is assigned an `orf_confidence` label in {`high`, `medium`, `low`, `none`}. The label combines two factors: a base score from the propagation outcome (1.0 for `propagated_intact`, 0.85 for `stop_at_alt_polya`, 0.7 for `disrupted`, 0.5 for `start_lost` and `stop_not_observed`, 0.3 for de novo predictions on truncated reads) and a multiplier from the completeness category (1.0 for `full_length`, 0.85 for `alt_3prime_end`, 0.7 for `truncated_3p` and `truncated_5p`, 0.5 for `truncated_both` and `internal_fragment`). The final score is binned into the four labels at thresholds 0.8, 0.6, 0.4. Confidence calibration on Bench 1 shows that HIGH-confidence calls have effectively zero error rate; LOW-confidence calls are still propagated but flagged for downstream filtering.

### 3.4 NMD rule cascade

For isoforms with a propagated stop codon that is observed in the read (outcome in {`propagated_intact`, `disrupted`}), CRAFT applies four rules in priority order. (1) If the stop codon falls inside the last exon, the isoform is NMD-escaped. (2) If the stop codon is within fifty nucleotides of the last exon-exon junction, the isoform is NMD-escaped (50nt rule, the canonical NMD escape window [15]). (3) If the propagated CDS is shorter than 150 nucleotides, the isoform is NMD-escaped (start-proximal escape [16]). (4) If the last exon is longer than 400 nucleotides, the isoform is NMD-escaped (long last exon escape). Isoforms surviving all four escape rules are labeled `sensitive` with rule `ptc_50nt_rule`. Isoforms without an observed stop (`stop_not_observed`, `stop_at_alt_polya`) are labeled `not_applicable` since the rule cascade cannot be applied.

### 3.5 Pfam domain disruption

When a Pfam HMM file is provided via `--pfam-hmm`, CRAFT translates the propagated ORF, scans against Pfam-A via pyhmmer, and compares the resulting domain set against the parent's annotated domains. Per-domain status is one of {`preserved`, `lost`, `gained`, `partial`, `frame_disrupted`}. Hits are cached by protein sequence SHA256 so repeated proteins across cells or samples are scanned once.

### 3.6 3' UTR features and polyA atlas matching

CRAFT computes 3' UTR length and length delta versus the parent for every isoform with a propagated stop codon. The canonical poly(A) signal motif (AATAAA and ten documented variants [17]) is scanned in the last fifty nucleotides of the isoform's transcript-orientation sequence. When `--polya-atlas` is supplied, the iso's 3' end position is matched against the pre-indexed atlas within a 24 nucleotide tolerance window; PAS evidence sources are recorded as `polya_db`, `canonical_motif`, or `none` in priority order. The PAS-evidence boolean is used as primary input to the `truncated_3p` versus `alt_3prime_end` reclassification.

### 3.7 Benchmark protocols

Bench 1 was implemented as `benchmarks/run_bench1.py`. The GENCODE v45 protein-coding transcript pool was loaded once from the GTF and cached as a Python pickle. For each (rate, orientation, seed) cell, 3,000 transcripts were sampled deterministically, truncated by exon-coordinate slicing (Methods 3.7.1), and written to a per-cell isoform GTF. CRAFT was invoked with the iso GTF and a per-cell reference GTF containing only the sampled parent transcripts (transcript, exons, CDS, start_codon, stop_codon). Orfipy was invoked once per cell on a single batched FASTA of all truncated transcript sequences. The GENCODE-truth ORF was excluded from CDS by convention (the stop codon is in a separate `stop_codon` record, not within the CDS record), matching orfipy's default and CRAFT's `propagated_cds_intervals` semantics.

Bench 2 was implemented as `benchmarks/run_bench2.py`. The 698,049-isoform bcM0003 collapsed GFF was streamed once into a `{transcript_id: (chrom, strand, exons)}` dictionary cached as a pickle. The CRAFT per-isoform TSV from a prior `craft annotate` run was joined against the GENCODE v45 transcript pool by `parent_tx_id`. For each isoform with a parent in the pool, the parent's start_codon and stop_codon genomic intervals were checked against the isoform's exon set; isoforms where either codon fell outside any exon were dropped. Surviving isoforms had their transcript-orientation sequences extracted from the GRCh38 primary assembly FASTA and submitted to orfipy in one batched FASTA. CRAFT's `propagated_cds_intervals` were parsed from the TSV and mapped to isoform transcript coordinates using the same `genomic_to_tx_coord` mapping as the truth construction, ensuring an apples-to-apples comparison.

Bench 3 was implemented as a three-script pipeline. `run_bench3_universe.py` filtered GENCODE v45 to the 80,441 transcripts of type protein_coding or nonsense_mediated_decay with start and stop codon annotations, ran CRAFT against the full GENCODE reference, and cached `nmd_status` per transcript. `run_bench3_quant.py` invoked salmon 1.11.4 with the v45 transcriptome k=31 index, in single-end mode (the library type detected by salmon as SR), for the three scrambled controls SRR4081222-224 and three UPF1 KDs SRR4081225-227 from GSE86148. `run_bench3_analysis.py` loaded the six `quant.sf` files into a transcript-by-sample count matrix, rounded NumReads to integers, applied a low-count filter (total ≥ 10 and at least two samples with nonzero count), and ran pydeseq2's Wald test on `~condition`. The merged DE results and CRAFT NMD universe were restricted to transcripts with orf_outcome in `{propagated_intact, disrupted}` (so that the NMD label was applied), and Fisher's exact test was applied to the 2x2 contingency of sensitive/escaped against upregulated/not.

#### 3.7.1 Truncation simulator

The exon truncator and the sequence truncator are kept in lock-step so the same (rate, orientation) call produces a coherent pair. Both work in transcript orientation (the minus-strand genomic-rightmost exon is the transcript's 5' end). For a transcript of total length L and rate r, the number of bases to trim is `int(L * r)`. For 5'-orientation truncation, this is removed from the start of the transcript; for 3'-orientation, from the end; for both, half from each end with any odd remainder going to the 3' side. The exon truncator slices exon intervals by transcript coordinates and maps back to genomic coordinates, yielding a list of (start, end) intervals in genomic order. The sequence truncator slices the transcript-orientation sequence by the same range. The two functions are unit-tested for length-consistency (20 tests covering both strands, all rates, all orientations).

### 3.8 Software environment

CRAFT requires Python ≥ 3.10. Runtime dependencies are pysam [18], pyranges [13], pandas, numpy, scipy, click, plotly, tqdm, pyhmmer [9], orfipy [14], anndata [19], and mudata. Development dependencies add pytest, pytest-cov, and ruff. Benchmarks also require pydeseq2 [12], salmon 1.11.4 [11], and the NCBI sra-toolkit [20]. The full test suite has 164 main-package tests plus 26 benchmark-library tests; ruff lint is clean. CRAFT is MIT-licensed.

---

## 4. Discussion

CRAFT validates a simple methodological point: when a long-read isoform shares enough splice structure with a reference transcript, the right ORF call comes from propagating the reference coordinates, not from rediscovering them. The five percentage point start-exact gap on synthetic truncation and the fifteen to twenty-eight point gap on real PacBio data are both consistent across conditions. The gap exists because de novo prediction must choose among multiple plausible ATGs in the truncated sequence; propagation has the answer already.

Three caveats limit the strength of the claims.

Bench 2's CRAFT-side number is partially tautological. The "truth" is defined as the parent's GENCODE CDS projected onto the isoform, and CRAFT's propagation uses the same parent's coordinates. The honest reading is the orfipy bar: real long-read isoforms have enough sequence-level confounders that de novo prediction misses the start 15-28 percent of the time. An independent truth source, for example bulk Iso-Seq from the same biological sample with an orthogonal ORF caller, would tighten this claim. We could not identify a published paired bulk-plus-single-cell PacBio dataset suitable for the comparison at the time of writing; LRGASP [21] WTC11 has bulk PacBio but the matched sc data uses a different platform (Nanopore), and the MAS-Seq paper's data is tumor T cells without matched bulk. Adding this benchmark when paired data becomes available is the most natural extension.

Bench 3's NMD effect size is modest. The 1.46 odds ratio is significant at p = 2.4e-15 because of the sample size (47k eligible transcripts), but the per-transcript classifier is a rule cascade against structural features, not a biochemical model. UPF1 substrate specificity has biochemical determinants beyond what the four rules capture (EJC deposition, secondary structure, codon usage, ribosomal occupancy), and 48-hour siRNA knockdowns also activate secondary pathways. We did not extend the benchmark to UPF2 or SMG6 knockdowns from the same study (the GSE86148 series includes these), which would test whether CRAFT's labels are specific to canonical NMD or also track ancillary decay pathways.

The de novo comparator is `orfipy` alone. CPAT [22] and TransDecoder [23] are the conventional bulk-RNA-seq baselines but they were excluded from the comparison for scope reasons (and because CRAFT's de novo fallback is `orfipy` itself, making the propagation-versus-de-novo question a direct test of the propagation logic rather than a multi-tool shootout). Future work should extend the comparison.

Beyond the benchmarks, the implementation includes design choices worth flagging. Mispriming detection (reads that primed off internal poly(A) tracts rather than the polyA tail) is intentionally not in CRAFT; it lives in our sibling tool `tecap` [24], and the two tools are complementary. Single-cell per-cell count integration into the AnnData export's `X` matrix is planned but not yet implemented; the current export populates `var` columns and leaves `X` empty unless the user provides a count matrix at runtime. SignalP and transmembrane annotations, disorder prediction, ClinVar disease-variant integration, and miRNA target site changes are all deferred to a future v2.

---

## 5. Data and code availability

Source code, runnable benchmark scripts, all committed figures, and the full test suite are at https://github.com/FullLengthFanatic/craft (MIT license). The v1.4 release is tagged at https://github.com/FullLengthFanatic/craft/releases/tag/v1.4. Reference annotations used in the benchmarks are GENCODE v45 (https://www.gencodegenes.org/human/release_45.html) and the GRCh38 primary assembly. The PolyASite v3.0 atlas is at https://polyasite.unibas.ch/. GSE86148 / SRP083135 RNA-seq data are available from the NCBI SRA.

---

## References

1. Al'Khafaji AM et al. (2024) High-throughput RNA isoform sequencing using programmable cDNA concatenation. Nat Biotechnol 42, 582-586.
2. Gupta I et al. (2018) Single-cell isoform RNA sequencing characterizes isoforms in thousands of cerebellar cells. Nat Biotechnol 36, 1197-1202.
3. Tian L et al. (2021) Comprehensive characterization of single-cell full-length isoforms in human and mouse with long-read sequencing. Genome Biol 22, 310.
4. Vitting-Seerup K, Sandelin A (2019) IsoformSwitchAnalyzeR. Bioinformatics 35, 4469-4471.
5. de la Fuente L et al. (2020) tappAS / IsoAnnotLite. Genome Biology 21, 119.
6. de la Fuente L et al. (2020) tappAS: a comprehensive computational framework for the analysis of the functional impact of differential splicing. Genome Biology 21, 119.
7. Herrmann CJ et al. (2020) PolyASite 2.0. Nucleic Acids Res 48, D174-D179. (v3.0 release: https://polyasite.unibas.ch/download/atlas/3.0/)
8. Wang R et al. (2018) PolyA_DB. Nucleic Acids Res 46, D315-D319.
9. Larralde M, Zeller G (2023) PyHMMER: a Python library binding to HMMER for efficient sequence analysis. Bioinformatics 39, btad214.
10. Boehm V et al. (2021) SMG5-SMG7 authorize nonsense-mediated mRNA decay. Cell Reports 35, 109137. (GSE86148)
11. Patro R et al. (2017) Salmon. Nat Methods 14, 417-419.
12. Muzellec B et al. (2023) PyDESeq2. Bioinformatics 39, btad547.
13. Stovner EB, Saetrom P (2020) PyRanges. Bioinformatics 36, 918-919.
14. Singh U, Wurtele ES (2021) orfipy. Bioinformatics 37, 3022-3024.
15. Nagy E, Maquat LE (1998) A rule for termination-codon position within intron-containing genes. Trends Biochem Sci 23, 198-199.
16. Lindeboom RGH et al. (2019) The rules and impact of nonsense-mediated mRNA decay. Nat Genet 48, 1112-1118.
17. Beaudoing E et al. (2000) Patterns of variant polyadenylation signal usage in human genes. Genome Res 10, 1001-1010.
18. Pysam developers. pysam: https://github.com/pysam-developers/pysam.
19. Virshup I et al. (2024) AnnData: Access and store annotated data matrices. JOSS 9, 4371.
20. NCBI SRA Toolkit. https://github.com/ncbi/sra-tools.
21. Pardo-Palacios FJ et al. (2024) Systematic assessment of long-read RNA-seq methods (LRGASP). Nat Methods.
22. Wang L et al. (2013) CPAT. Nucleic Acids Res 41, e74.
23. Haas BJ et al. (2013) De novo transcript sequence reconstruction from RNA-seq using the Trinity platform. Nat Protoc 8, 1494-1512.
24. Picelli S. tecap: 3' end / priming-artifact diagnostics for long-read RNA-seq. https://github.com/FullLengthFanatic/tecap.

---

## Figure legends

**Figure 1.** Bench 1 recovery panel. 2x2 plotly layout. Top-left: recovery rate by truncation rate, faceted by orientation. Top-right: start-codon exact-match rate. Bottom-left: mean absolute ORF length error (nucleotides). Bottom-right: stop-codon exact-match rate. CRAFT in muted slate blue, orfipy in warm orange; solid = 5'-truncation, dashed = 3'-truncation, dotted = both. The start-exact panel is the central observation: orfipy plateaus at 0.94 to 0.95 across all conditions; CRAFT sits at 0.98 to 1.00. Source: `benchmarks/figures/bench1_recovery_panel.png`.

**Figure 2.** Bench 2 concordance panel. Left: start-codon exact match by CRAFT completeness category, side-by-side CRAFT propagation and orfipy de novo. Per-category n annotated above each pair. Right: CRAFT start-exact rate by `orf_confidence` label (the panel is uninformative for this benchmark because the truth definition makes CRAFT trivially correct on intact-truth isoforms; included as a sanity check). Source: `benchmarks/figures/bench2_concordance_panel.png`.

**Figure 3.** Bench 3 enrichment panel. Left: 2x2 contingency heatmap of CRAFT NMD label by UPF1-KD differential-expression status with annotated cell counts and odds ratio. Middle: volcano plot of log2 fold-change against -log10 adjusted p, with NMD-sensitive transcripts overlaid in slate blue on NMD-escaped in light gray; the +1 log2FC and 0.05 padj cutoffs are marked with dotted reference lines. Right: cumulative distribution of log2 fold-change stratified by CRAFT NMD label; the rightward shift of the NMD-sensitive curve is the directional signature of NMD inhibition. Source: `benchmarks/figures/bench3_enrichment_panel.png`.
