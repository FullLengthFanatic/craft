# CRAFT methods and design rationale

This document explains every category and metric CRAFT emits, how each is
calculated, and the design decision behind it. It's the reference for
interpreting the per-isoform output table and for understanding why the
pipeline is shaped the way it is.

The companion `user_guide.md` covers how to run the tool. This file is the
methods reference.

## Contents

1. [What CRAFT does and what it doesn't](#what-craft-does-and-what-it-doesnt)
2. [Inputs](#inputs)
3. [Pipeline at a glance](#pipeline-at-a-glance)
4. [Module 1: splice-junction extraction](#module-1-splice-junction-extraction)
5. [Module 2: completeness classification](#module-2-completeness-classification)
6. [Module 3: ORF propagation](#module-3-orf-propagation)
7. [Module 4: de novo ORF prediction](#module-4-de-novo-orf-prediction)
8. [Module 5: ORF confidence scoring](#module-5-orf-confidence-scoring)
9. [Module 6: NMD susceptibility](#module-6-nmd-susceptibility)
10. [Module 7: 3' UTR features](#module-7-3-utr-features)
11. [Module 8: Pfam domain disruption (optional)](#module-8-pfam-domain-disruption-optional)
12. [Output schema](#output-schema)
13. [Interpreting CRAFT calls (FAQ)](#interpreting-craft-calls-faq)
14. [Design rationale and trade-offs](#design-rationale-and-trade-offs)
15. [Known limitations (v1)](#known-limitations-v1)

---

## What CRAFT does and what it doesn't

CRAFT (Coding Region Annotation From Templates) takes the output of any
long-read isoform caller (FLAIR, IsoQuant, Bambu, FLAMES, SQANTI3,
isoseq+pigeon) and emits per-isoform functional consequences:

- structural completeness vs a reference parent transcript
- coding region: where the ORF is, whether propagation from the parent
  worked, and how much we trust the call
- NMD susceptibility (50nt rule + escapes)
- 3' UTR length delta and poly(A) signal motif
- Pfam domain preservation/loss/gain (when `--pfam-hmm` is supplied)

**The methods novelty is reference-isoform ORF propagation with explicit
truncation-aware confidence.** Most isoforms produced by long-read pipelines
are partially truncated. De-novo ORF predictors do badly on truncated
reads. CRAFT projects the parent's CDS coordinates onto the isoform where
structure is preserved, flags exactly where the call becomes uncertain
(start codon not in the read, stop codon not in the read, structural
divergence), and emits a confidence score callers can filter on.

**CRAFT does not do structural QC.** It assumes the iso GTF is already
post-QC. It does not filter "junk" isoforms; it describes what's there.
The `orf_confidence` column gives downstream code a principled way to drop
low-quality calls.

---

## Inputs

| Input              | Required | Format                                           | Used by                                            |
| ------------------ | -------- | ------------------------------------------------ | -------------------------------------------------- |
| Isoform GTF        | yes      | GTF with `exon` rows and `transcript_id`         | every module                                       |
| Reference GTF      | yes      | GTF with `exon` AND `CDS` rows                   | completeness, propagation, NMD, UTR, Pfam          |
| Genome FASTA       | yes      | indexed FASTA (`.fai` built if missing)          | denovo, UTR poly(A) scan, Pfam translation         |
| Per-cell counts    | no       | `.h5ad` or 10x-style MTX directory               | export (per-cell X in `annotated.h5ad`)            |
| Pfam HMM           | no       | `.hmm` file (e.g. Pfam-A.hmm)                    | Pfam domain scan                                   |

All three required inputs must use the same chromosome naming (`chr1` vs
`1`). CRAFT does not do chromosome-name harmonisation; pyranges raises if
the FASTA is missing a chromosome the GTF references.

---

## Pipeline at a glance

```
                       isoform GTF
                            |
                            v
                   io.gtf.load_isoforms
                            |
                            v
                  core.completeness.classify    <-- reference GTF (exon-only)
                            |  parent_tx_id, completeness, shared_junctions
                            v
                  core.orf.propagation.propagate <-- reference GTF (full)
                            |  orf_outcome, propagated_cds_intervals,
                            v  start/stop coverage flags
            orphans (no parent / no CDS / start lost)?
              |                                  \
              v                                   v
   core.orf.denovo.predict                  (pass-through)
   (orfipy ATG search                            |
    on transcript sequence)                      |
              |                                  |
              v                                  v
                       merge into per-isoform
                            |
                            v
                  core.orf.confidence.score
                            |
                            v
                  core.nmd.predict
                            |
                            v
                  core.utr3.annotate <-- genome FASTA (poly(A) scan)
                            |
                            v
                  core.pfam.scan (optional)
                            |
                            v
              per-isoform DataFrame (30 cols)
                            |
              +-------------+-------------+
              v             v             v
     per_isoform.tsv  per_isoform.json  report.html
                                              |
                                              v
                                       annotated.h5ad
                                       (optional counts in X)
```

---

## Module 1: splice-junction extraction

**Source:** `src/craft/core/intervals.py::splice_junctions`

**What it does.** Computes the set of splice junctions for each transcript
(iso or reference). A junction is the intron between two consecutive exons
in transcript order, represented in genomic coordinates as the half-open
interval `[upstream_exon.End, downstream_exon.Start)`.

**Algorithm.** Sort exons by `(transcript_id, Start)`, then use pandas
`groupby().shift(-1)` to pair each exon with its successor within the same
transcript. Single-exon transcripts contribute zero junctions.

**Why this representation.** Two transcripts share a junction if and only
if their `(Chromosome, Start, End, Strand)` tuples match exactly. That
makes "junction overlap" a cheap pandas merge instead of an interval-tree
intersection, which is the main hot path of `completeness.classify` on
large reference GTFs (GENCODE has ~3M junctions).

---

## Module 2: completeness classification

**Source:** `src/craft/core/completeness.py::classify`

**What it does.** For each iso, identifies the best-matching reference
transcript (the "parent") and classifies the iso's structural completeness
relative to that parent.

**Categories.** The enum has six values:

| Category            | Meaning                                                                                   |
| ------------------- | ----------------------------------------------------------------------------------------- |
| `full_length`       | iso spans the parent end-to-end (within tolerance)                                        |
| `truncated_5p`      | iso 3' end matches parent; 5' end is interior                                             |
| `truncated_3p`      | iso 5' end matches parent; 3' end is interior; NO canonical poly(A) signal in last 50 nt  |
| `alt_3prime_end`    | iso 3' end is interior to parent, but a canonical poly(A) signal sits in the last 50 nt — biological APA, not technical truncation |
| `truncated_both`    | both ends interior to parent, but near at least one boundary                              |
| `internal_fragment` | both ends well-interior to parent (>2× tolerance from both)                               |
| `novel_no_match`    | no shared junctions AND no exon overlap with any reference tx                             |

**The poly(A) split between `truncated_3p` and `alt_3prime_end`** runs at the pipeline level (not inside `classify`) because it requires the genome FASTA. The split logic: when an iso would be classified `truncated_3p`, scan its last 50 nt (in transcript orientation; reverse-complemented for `-` strand) for any of the 11 canonical poly(A) signal motifs from `core/utr3.py::POLYA_SIGNALS`. If a motif is found, relabel as `alt_3prime_end`. Otherwise stay `truncated_3p`. For oligo-dT primed long-read cDNA, the vast majority of "shorter than parent annotation" cases land in `alt_3prime_end` because the iso's 3' end *is* the polyadenylation site by construction of the library prep.

**Parent selection.** For each iso, score every candidate reference
transcript by (a) number of exactly shared splice junctions, then (b)
total stranded exon-overlap bp as tiebreaker. The candidate with the
highest (a, b) wins. Ties below the inclusion threshold (zero shared
junctions AND zero exon overlap) are classified `novel_no_match`.

**End-comparison logic.** Once a parent is picked, completeness is decided
by whether each end of the iso reaches the parent's boundary within
`tolerance` (default 50 bp). The check is strand-aware:

- For `+` strand: 5'-complete if `iso.Start <= parent.Start + tolerance`;
  3'-complete if `iso.End >= parent.End - tolerance`.
- For `-` strand: 5'-complete if `iso.End >= parent.End - tolerance`;
  3'-complete if `iso.Start <= parent.Start + tolerance`.

`INTERNAL_FRAGMENT` requires both ends to be deeper than `2 * tolerance`
inside the parent (or at least 100 bp interior, whichever is larger).
Below that depth, `TRUNCATED_BOTH` is used.

**Tolerance default.** 50 bp. Long-read TSS/TES uncertainty is typically
10-30 bp from the true site; 50 bp absorbs that without misclassifying
truly truncated isoforms.

**Rationale.** Junction sharing dominates exon overlap as a parent-match
signal because it's a structural assertion ("these two transcripts use
the same splice sites"), whereas exon overlap can be inflated by long
constitutive exons that many isoforms share. The tiebreak by exon overlap
covers cases where junction count ties or where one transcript is
single-exon (no junctions to share).

**What this is NOT.** This is not SQANTI3 classification. SQANTI3
categories (FSM, ISM, NIC, NNC) describe junction-set relationships at a
finer grain. CRAFT's completeness is end-position-focused because the
downstream propagation step needs to know "does this iso cover the
parent's start/stop codon positions?", and that's a TSS/TES question, not
a junction-set question.

---

## Module 3: ORF propagation

**Source:** `src/craft/core/orf/propagation.py::propagate`

**The motivation.** In single-cell long-read data, the modal isoform is
not a full-length read. A 5'-truncated read of an annotated protein-coding
transcript may still carry the entire CDS region (truncation in the 5'
UTR). De-novo ORF prediction on such a read will either give up (no ATG
in frame) or find a downstream Met that's not the real start. The parent
reference transcript, on the other hand, has a precisely annotated start
and stop. **If the iso's exon structure preserves the parent's CDS
intervals, we can project the parent's CDS onto the iso and call the ORF
with high confidence.** That's propagation.

**Algorithm.**

1. Look up the parent's CDS records (Feature == "CDS") for the
   iso's `parent_tx_id` (assigned by completeness).
2. Identify the parent's start codon genomic position
   (`min(CDS.Start)` on `+`, `max(CDS.End) - 1` on `-`) and stop codon
   genomic position (the opposite end).
3. Check whether the iso has any exon containing the start codon position
   (`start_codon_covered`) and the stop codon position
   (`stop_codon_covered`).
4. Intersect each parent CDS interval with each iso exon (within the same
   chromosome and strand) to compute `propagated_cds_intervals`. Sum the
   intersected lengths to get `propagated_cds_bp`; the parent's CDS length
   gives `parent_cds_bp`.
5. Classify the outcome (see table below).

**`ORFOutcome` values.**

| Outcome              | When                                                                                                                 |
| -------------------- | -------------------------------------------------------------------------------------------------------------------- |
| `propagated_intact`  | start and stop both covered AND propagated_cds_bp == parent_cds_bp                                                   |
| `stop_at_alt_polya`  | start covered, parent's stop not covered, but a canonical poly(A) signal sits in the iso's last 50 nt — alt-polyA upstream of canonical stop (biological APA, not truncation) |
| `stop_not_observed`  | start covered, parent's stop not covered, AND no poly(A) signal nearby — likely real technical 3' truncation         |
| `disrupted`          | start AND stop covered, but propagated_cds_bp < parent_cds_bp (structural change)                                    |
| `start_lost`         | start codon position not in any iso exon                                                                             |
| `no_parent_cds`      | parent was identified but has no CDS records (e.g. lncRNA, pseudogene)                                               |
| `no_parent`          | iso has no usable parent (NOVEL_NO_MATCH)                                                                            |

**Same poly(A) split as completeness** (see above): when an iso would be `stop_not_observed`, the pipeline checks for poly(A) evidence near the iso's 3' end and reclassifies to `stop_at_alt_polya` if found. The split is mechanistically meaningful: oligo-dT primed reads almost always have a polyA tail (and therefore a polyA signal) at their 3' end by construction, so "iso doesn't reach parent's annotated stop" is almost always APA, not truncation.

**Two evidence sources** (added in v1.2):

1. **Atlas match (primary when `--polya-atlas` is provided).** A user-supplied BED file of curated polyadenylation sites (PolyASite v3.0, PolyA_DB v4, or any 6+ column BED) is loaded and indexed by chromosome+strand. The iso's 3' end position is compared against the atlas: a PAS midpoint within ±24 nt on the same chromosome+strand counts as a hit. When a match is found, the iso is reclassified and `polya_evidence_source` is set to `polya_db` with the matched PAS name in `polya_db_site_id`. **The atlas should be pre-filtered by usage score** before passing to CRAFT — the unfiltered PolyASite v3.0 atlas has ~one PAS every 200 bp and produces uninformative 98% hit rates at the default ±24 nt tolerance; `awk '$5 >= 0.01'` cuts that to ~88% atlas hits while keeping the HIGH-confidence ORF fraction effectively unchanged (32.5% vs 33.7% unfiltered) and dropping runtime ~3-fold. See `docs/user_guide.md` for the full recipe.
2. **Canonical motif scan (fallback).** When the atlas is not provided OR the iso's 3' end gets no atlas hit, the v1.1 motif scanner runs: the last 50 nt of the iso's transcript-orientation 3' end is searched for the 11 canonical poly(A) signal motifs (`AATAAA` and known variants). A hit sets `polya_evidence_source = canonical_motif`. No hit → `none`, and the original `TRUNCATED_3P` / `STOP_NOT_OBSERVED` label stands.

Both pathways feed the same boolean ("did we find evidence?") into the reclassification step. The `polya_evidence_source` column lets downstream filters distinguish DB-supported APA from motif-only APA when the user cares.

**Why this priority order.** `START_LOST` is evaluated first because
without the start codon there's no propagation regardless of what else
the iso covers. `STOP_NOT_OBSERVED` is next because if we can't see the
stop, the CDS length is undefined; even if `propagated_cds_bp ==
parent_cds_bp` happens to match (a fluke), we don't trust it.
`PROPAGATED_INTACT` requires both endpoints AND length parity.
`DISRUPTED` is the catch-all for "structurally different but the
endpoints we care about are observable".

**What v1 does NOT yet detect.**

- **Intron retention in the CDS region.** If the iso retains an intron
  inside the CDS but still preserves every parent CDS bp, v1 calls it
  `PROPAGATED_INTACT`. In reality this almost certainly introduces a
  premature stop in the retained intron and disrupts the protein. v1.5
  detects intron retention by comparing iso junctions to parent
  junctions within the CDS region; if iso has fewer junctions, flag
  `intron_retained`.
- **Frame tracking through alt splice sites.** If an iso uses an alt 5'ss
  that shifts the frame by a non-multiple of 3, v1 marks the iso
  `DISRUPTED` but doesn't predict the resulting premature stop position.
  v1.5 walks the propagated CDS in 3-nt steps and finds the first
  in-frame TAA/TAG/TGA after the frame shift.

---

## Module 4: de novo ORF prediction

**Source:** `src/craft/core/orf/denovo.py::predict`

**When it runs.** Only for isoforms whose propagation outcome is
`NO_PARENT`, `NO_PARENT_CDS`, or `START_LOST`. Other outcomes already
have a usable ORF from propagation; running de-novo on top would just
add noise.

**Algorithm.**

1. Build the iso's transcript sequence: concatenate exon sequences in
   genomic order, reverse-complement for `-` strand so the result is
   transcript-orientation (5' → 3').
2. Run `orfipy_core.orfs(seq, minlen=min_bp, starts=["ATG"], strand="f")`.
   `strand="f"` because the sequence is already in transcript orientation.
3. Pick the longest ORF among the returned candidates.
4. Project the ORF's transcript coordinates back to genomic intervals
   (strand-aware: `+` walks exons in genomic order, `-` walks reversed
   and inverts the position math).

**Why ATG-only starts.** orfipy's defaults include `TTG` and `CTG` as
alternative starts. For v1 we stick to the most stringent biologically
common case (ATG). Non-ATG starts exist but are rare and add false
positives. v1.5 can add an option flag.

**Default minimum ORF length: 50 aa.** A protein under 50 aa is
biologically uncommon (most short ORFs are smORFs, which are an
interesting but distinct class). 50 aa = 150 bp of CDS, which is also
above the noise threshold for de-novo ATG-Stop pairs in random sequence
(a random sequence has a TAA/TAG/TGA every ~21 codons by chance).

**Confidence consequence.** Because de-novo ORFs lack a reference
anchor, the pipeline downgrades their `orf_confidence` to `LOW` even
when the de-novo call itself looks clean. See module 5.

---

## Module 5: ORF confidence scoring

**Source:** `src/craft/core/orf/confidence.py::score`

**What it returns.** A categorical confidence (`HIGH`, `MEDIUM`, `LOW`,
`NONE`) and a numeric score in `[0, 1]`. Both go into the per-isoform
output as `orf_confidence` and `orf_confidence_score`.

**Score formula.** `score = base(outcome) * factor(completeness)`,
clipped to `[0, 1]`.

`base(outcome)`:

| outcome             | base |
| ------------------- | ---- |
| propagated_intact   | 1.0  |
| stop_at_alt_polya   | 0.85 |
| stop_not_observed   | 0.55 |
| disrupted           | 0.45 |
| start_lost          | 0.2  |
| no_parent           | (NONE) |
| no_parent_cds       | (NONE) |

`factor(completeness)`:

| completeness        | factor |
| ------------------- | ------ |
| full_length         | 1.0    |
| alt_3prime_end      | 1.0    |
| truncated_5p        | 0.9    |
| truncated_3p        | 0.9    |
| truncated_both      | 0.65   |
| internal_fragment   | 0.5    |
| novel_no_match      | (NONE) |

The `alt_3prime_end` factor is 1.0 (no penalty) because the iso has a clean biological alternative end, not a technical truncation. Likewise, `stop_at_alt_polya` carries the second-highest base score (0.85) because the iso has its own valid stop site via the alt-polyA, just one that's not the parent's annotated stop. Both new categories were added in v1.1 to fix a systematic mislabeling of APA isoforms as truncations.

Categorical thresholds: `>= 0.85` → HIGH, `>= 0.5` → MEDIUM, else LOW.

**Denovo override.** If propagation gave `NO_PARENT` or `NO_PARENT_CDS`
(so `score()` returns `NONE`, 0.0), but the denovo step found an ORF,
the pipeline overrides the confidence to `LOW` (numeric 0.25). De-novo
calls without a reference anchor are by construction less trustworthy
than even a disrupted propagation, so we don't promote them above LOW.

**Why multiplicative.** It enforces a natural ordering: at the same
completeness, PROPAGATED_INTACT >= STOP_NOT_OBSERVED >= DISRUPTED >=
START_LOST; at the same outcome, full-length scores at least as high as
any truncated variant. Both invariants are tested in
`tests/test_orf_confidence.py`.

**Calibration of the constants.** The base values and factors were
chosen so that:

- A full-length intact propagation lands at 1.0 (the strongest signal we
  have).
- A 5'-truncated iso with intact propagation (very common in single-cell
  data) stays HIGH (0.9), because the start codon being observed in the
  read means the ORF is unambiguous regardless of how much 5' UTR is
  missing.
- A full-length DISRUPTED call drops to LOW (0.45), because structural
  divergence inside the CDS region means we can't say what protein the
  iso actually makes.
- An internal fragment with intact propagation lands at MEDIUM (0.5),
  because while we know the propagated region is preserved, we have no
  evidence of the endpoints.

These are not learned from data; they're a deliberate ordering. They can
be tuned by passing different threshold constants if you have a labelled
set to calibrate against.

---

## Module 6: NMD susceptibility

**Source:** `src/craft/core/nmd.py::predict`

**What it predicts.** A single NMD call per isoform from the **resolved** ORF
stop (the real in-frame stop from the sequence-resolution step below), falling
back to the de-novo ORF stop for orphan isoforms. `nmd_basis` records which ORF
was used. As of v1.7 the older geometric NMD (applied to the projected interval
end) has been removed; the rule cascade below now operates on the resolved stop.
It applies a rule cascade to decide whether the transcript would trigger
nonsense-mediated decay.

**Rule cascade** (each is a sufficient condition for ESCAPE, evaluated in
priority order):

| Order | Rule                                              | Default threshold |
| ----- | ------------------------------------------------- | ----------------- |
| 1     | Stop in last exon (or single-exon transcript)     | n/a               |
| 2     | Stop within N nt of last exon-exon junction       | 50 nt (mRNA)      |
| 3     | CDS shorter than start-proximal threshold         | 150 nt            |
| 4     | Last exon longer than long-last-exon threshold    | 400 nt            |

If none of the escape rules fire, the iso is `SENSITIVE`.

**Distance calculation.** Distance from the stop codon to the last
exon-exon junction is measured in mRNA (transcript) nucleotides, not
genomic distance. The algorithm walks the iso exons in transcript order
(strand-aware) from the stop-codon-containing exon toward the
transcript's last exon, summing the remaining length of the stop exon
plus the full lengths of any intermediate exons.

**Stop codon position.** Derived from `propagated_cds_intervals`:
`max(End) - 1` on `+` strand, `min(Start)` on `-` strand. This is the
genomic position of the last CDS base in transcript order.

**Confidence.** `nmd_confidence` is HIGH when the propagation outcome
was `PROPAGATED_INTACT`, MEDIUM when it was `DISRUPTED` (rules still
apply but structural changes mean the stop might be in a different
location than the propagation placed it), and NONE when the rule wasn't
applicable.

**When NMD is `not_applicable`.** Whenever the stop codon is not
observed in the read (`STOP_NOT_OBSERVED`), the ORF was not propagated
(`START_LOST`, `NO_PARENT`, `NO_PARENT_CDS`), or the iso has no
propagated CDS at all. The pipeline does not apply NMD rules in these
cases because the stop-to-junction distance is undefined.

**Rule defaults: where they come from.**

- 50 nt PTC rule: Lindeboom et al. 2019 and the IsoformSwitchAnalyzeR
  default. Some labs use 55 nt. Configurable.
- 150 nt start-proximal escape: covers the re-initiation regime where
  short uORFs evade NMD because the ribosome stays loaded.
- 400 nt long-last-exon escape: Lindeboom et al. observed that
  unusually long last exons correlate with NMD escape even when the
  stop is upstream of the last junction. Mechanism less clear than the
  50nt rule; the value is empirical.

**Important interpretation note.** NMD calls in a healthy non-disease
sample are not artefactual. NMD is also a regulatory mechanism: cells
use alternative splicing to insert premature stops in some isoforms to
control parent-gene expression (regulated unproductive splicing and
translation, RUST; alternatively, AS-NMD). Splicing factors auto-regulate
this way. So a baseline rate of NMD-sensitive isoforms (single-digit
percent in our test data) is expected. See the FAQ section below for a
fuller treatment, including how to distinguish biological NMD from
truncation-artefact apparent NMD.

---

## Module 7: 3' UTR features

**Source:** `src/craft/core/utr3.py::annotate`

**Two outputs.**

1. **3' UTR length delta vs the parent.** For each iso with a propagated
   ORF, sum the iso's exonic bp downstream of the stop codon (`+`
   strand: positions strictly greater than the stop position; `-`
   strand: positions strictly less). Do the same for the parent's CDS +
   exons. Report `iso_utr3_length_nt`, `parent_utr3_length_nt`,
   `utr3_length_delta_nt`, and `utr3_length_delta_pct`.

2. **Poly(A) signal motif scan.** When `--genome` is provided, extract
   the iso's 3' UTR sequence (reverse-complement for `-` strand so the
   scan runs in transcript orientation) and search for canonical poly(A)
   signal motifs.

**Motif priority list** (canonical first, then variants):

```
AATAAA, ATTAAA, AGTAAA, TATAAA, CATAAA, GATAAA,
AATATA, AATACA, AATAGA, AAAAAG, ACTAAA
```

For each motif in priority order, take the *rightmost* (most 3'-proximal)
occurrence. First motif with any hit wins. Priority dominates distance:
if both AATAAA and ATTAAA are present, AATAAA is reported regardless of
which is closer to the cleavage site.

**Distance reported.** Nucleotides from the end of the motif to the
end of the 3' UTR sequence (i.e. roughly nt to the inferred cleavage
site). Canonical signals usually sit 10-30 nt upstream of the cleavage
site.

**When UTR analysis is NOT applicable.** Same gate as NMD: the iso must
have a propagated ORF (`PROPAGATED_INTACT` or `DISRUPTED`) with stop
codon observed. Otherwise iso/parent UTR lengths are `None` and the
poly(A) motif is empty.

**Why not detect internal priming here.** That's tecap's job. CRAFT
intentionally does not duplicate tecap's mispriming classification. If
you've run tecap, its per-read classifications can be joined to CRAFT's
per-isoform output downstream by transcript_id.

---

## Module 8: Pfam domain disruption (optional)

**Source:** `src/craft/core/pfam.py::scan`

**When it runs.** Only when `--pfam-hmm /path/to/Pfam-A.hmm` is supplied.
Skipped otherwise (the five `pfam_*` columns are present in the output
but populated with empty lists).

**Workflow.**

1. Translate each iso's CDS to protein. CDS source: propagated CDS if
   present, otherwise the de-novo CDS (so novel-with-de-novo isoforms
   still get a domain set, but with no parent comparison they end up in
   `pfam_gained`).
2. Translate each parent's CDS the same way.
3. Run `pyhmmer.hmmsearch` on each unique protein sequence against the
   HMM database. Hits are filtered by `hit.included` (HMMer's default
   inclusion threshold, evalue ≤ 0.01).
4. Cache hits by SHA256 of the protein sequence. Repeated proteins
   (common across cells in single-cell data) scan only once.
5. For each iso, compute `pfam_preserved = iso ∩ parent`,
   `pfam_lost = parent - iso`, `pfam_gained = iso - parent`.

**Performance.** v1 uses `hmmsearch` (HMMs as the outer loop). On full
Pfam-A (~20k HMMs) this is slow: roughly an hour per few thousand
isoforms. v1.5 will switch to `hmmscan` against a pressed HMM database
(`hmmpress Pfam-A.hmm` first) which is asymptotically much faster
because it uses a k-mer prefilter.

**Codon table.** Standard nuclear codon table, frame-0, stops at the
first stop codon, unknown codons (any base outside ACGT) → `X`. Not
configurable in v1; the table is in `_CODON_TABLE` if you need to
override.

---

## Module 9: Per-cell recurrence signals (v1.8)

**Source:** `src/craft/core/recurrence.py` (`compute_recurrence`,
`within_gene_fraction`, `load_cell_whitelist`), wired into `src/craft/pipeline.py`.

**When it runs.** Only when `--counts` is supplied. Populates three new columns
in the per-isoform output: `total_count`, `n_cells_detected`, and
`isoform_fraction_within_gene`. All remain null/empty when `--counts` is absent
or the isoform is absent from the count matrix.

**Rationale.** Raw molecule counts are depth-dependent: a cell with 10x higher
sequencing depth will show 10x higher counts for every isoform, inflating total
counts and drowning out rare isoforms. Cell recurrence (the count of independent
cells with at least one molecule) is depth-stable; an isoform seen in many cells
is supported regardless of per-cell sequencing depth. For filtering, cell
recurrence is more robust than read/count thresholds.

**Computation.**

1. Load per-cell counts from `--counts` (either `.h5ad` or 10x MTX). Optionally
   filter to called cells if `--cell-whitelist` is provided (a text file with
   one barcode per line). Without the whitelist, every barcode in the matrix is
   used, including ambient droplets.
2. For each isoform, sum UMI-corrected molecules across the selected cell set
   (`total_count`).
3. Count cells with at least one molecule of the isoform (`n_cells_detected`).
4. For each isoform, divide its `total_count` by the summed `total_count` of all
   isoforms with the same `parent_gene_id` (`isoform_fraction_within_gene`).

The same computation is also exported to `annotated.h5ad` as `var` columns.

---

## Output schema

CRAFT writes four files to `output_dir/`:

### `per_isoform.tsv` (one row per isoform, 63 columns)

The full annotation table. List-valued columns (CDS intervals, Pfam
domain lists) are JSON-encoded so the TSV is grep-friendly while still
being round-trippable.

| Column                          | Source                | Type      | Meaning                                                                  |
| ------------------------------- | --------------------- | --------- | ------------------------------------------------------------------------ |
| `transcript_id`                 | iso GTF               | string    | iso identifier                                                           |
| `completeness`                  | completeness          | category  | full_length / truncated_5p / 3p / alt_3prime_end / truncated_both / internal_fragment / novel_no_match |
| `parent_tx_id`                  | completeness          | string    | best-matching reference transcript_id; "" if none                       |
| `parent_gene_id`                | pipeline lookup       | string    | reference `gene_id` for the parent transcript (looked up from the reference GTF) |
| `parent_gene_name`              | pipeline lookup       | string    | reference `gene_name` for the parent (empty when the reference GTF lacks it)    |
| `shared_junctions`              | completeness          | int       | exactly-shared splice junctions with parent                              |
| `parent_overlap_bp`             | completeness          | int       | total stranded exon overlap with parent                                  |
| `orf_outcome`                   | propagation           | category  | propagated_intact / stop_at_alt_polya / stop_not_observed / disrupted / start_lost / no_parent / no_parent_cds |
| `propagated_cds_bp`             | propagation           | int       | bp of parent CDS preserved in iso                                        |
| `parent_cds_bp`                 | propagation           | int       | parent's total CDS bp                                                    |
| `start_codon_covered`           | propagation           | bool      | iso exons span the parent's start codon position                         |
| `stop_codon_covered`            | propagation           | bool      | iso exons span the parent's stop codon position                          |
| `propagated_cds_intervals`      | propagation           | list[tup] | propagated CDS as genomic intervals (chr, start, end, strand)            |
| `denovo_orf_found`              | denovo                | bool      | de-novo orfipy search returned a candidate ≥ min_orf_aa                  |
| `denovo_cds_bp`                 | denovo                | int       | de-novo ORF length (nt)                                                  |
| `denovo_orf_aa_length`          | denovo                | int       | de-novo ORF length (aa)                                                  |
| `denovo_start_codon`            | denovo                | string    | start codon used (always "ATG" in v1)                                    |
| `denovo_stop_codon`             | denovo                | string    | stop codon ("TAA"/"TAG"/"TGA"/"")                                        |
| `denovo_cds_intervals`          | denovo                | list[tup] | de-novo CDS as genomic intervals                                         |
| `orf_confidence`                | confidence + denovo   | category  | high / medium / low / none                                               |
| `orf_confidence_score`          | confidence            | float     | numeric in [0, 1]                                                        |
| `nmd_status`                    | nmd                   | category  | sensitive / escaped / not_applicable                                     |
| `nmd_rule`                      | nmd                   | string    | which escape rule fired, or `ptc_50nt_rule` for sensitive                |
| `stop_to_last_junction_nt`      | nmd                   | int       | mRNA distance from stop to last junction (0 if stop in last exon)        |
| `last_exon_length_nt`           | nmd                   | int       | iso's last exon length (transcript order)                                |
| `nmd_confidence`                | nmd                   | category  | high / medium / none                                                     |
| `iso_utr3_length_nt`            | utr3                  | float     | iso's 3' UTR length; null if NMD-non-applicable                          |
| `parent_utr3_length_nt`         | utr3                  | float     | parent's 3' UTR length                                                   |
| `utr3_length_delta_nt`          | utr3                  | float     | iso - parent UTR length                                                  |
| `utr3_length_delta_pct`         | utr3                  | float     | percent change relative to parent                                        |
| `polya_signal_motif`            | utr3                  | string    | strongest poly(A) motif found in the UTR; "" if none                     |
| `polya_signal_distance_nt`      | utr3                  | float     | nt from motif end to UTR 3' end                                          |
| `polya_evidence_source`         | polya_atlas + utr3    | string    | `polya_db` (atlas hit), `canonical_motif` (motif fallback), or `none`    |
| `polya_db_site_id`              | polya_atlas           | string    | PAS name from the atlas BED's column 4 when matched; empty otherwise     |
| `iso_pfam_domains` *            | pfam                  | list[str] | Pfam HMM names hitting the iso's protein                                 |
| `parent_pfam_domains` *         | pfam                  | list[str] | Pfam HMM names hitting the parent's protein                              |
| `pfam_preserved` *              | pfam                  | list[str] | iso_domains ∩ parent_domains                                             |
| `pfam_lost` *                   | pfam                  | list[str] | parent_domains - iso_domains                                             |
| `pfam_gained` *                 | pfam                  | list[str] | iso_domains - parent_domains                                             |
| `total_count`                   | recurrence (v1.8)     | int / null | UMI-corrected molecules summed across cells. Computed over called cells when `--cell-whitelist` is given; over every barcode otherwise. Null without `--counts` or if isoform not in count matrix. |
| `n_cells_detected`              | recurrence (v1.8)     | int / null | Cells with at least one molecule; depth-stable recurrence signal. Null without `--counts`. |
| `isoform_fraction_within_gene`  | recurrence (v1.8)     | float / null | `total_count` / summed `total_count` of isoforms sharing `parent_gene_id`. Null for isoforms without parent gene or without measured counts. |

*Pfam columns are empty lists unless `--pfam-hmm` was supplied.

### `per_isoform.json`

Same content as the TSV but as a JSON array of records. List columns
stay as lists (not JSON-encoded), so it's the more convenient format for
programmatic consumers.

### `report.html`

Self-contained interactive HTML. Three sections:

1. **Summary cards** — per-category counts and percentages for
   completeness, ORF outcome, NMD status, and ORF confidence.
2. **Distributions** — plotly bar charts of the same four categorical
   fields. Plotly.js is inlined in the first figure block (single
   `#5b7a9d` slate fill across all charts); no CDN required.
3. **Notable findings** — three small focused tables instead of a
   per-isoform dump (the full data lives in `per_isoform.tsv`):

   - **Top NMD-sensitive isoforms** (max 10 rows): isoforms where
     `nmd_status == "sensitive" AND nmd_confidence == "high"`, sorted
     by `orf_confidence_score` descending. Filters to the
     biologically-trustworthy NMD substrates.
   - **Top ORF-disrupted isoforms** (max 10 rows): isoforms where
     `orf_outcome == "disrupted" AND orf_confidence == "high"`, sorted
     by `(parent_cds_bp - propagated_cds_bp)` descending. Surfaces the
     iso/parent pairs where the most CDS bp were lost in the
     structural change.
   - **Genes with most functional isoform diversity** (max 10 rows):
     for each parent gene, the number of *distinct*
     `(parent_tx_id, orf_outcome)` pairs among isoforms with
     `orf_confidence in {"high", "medium"}`. This deliberately collapses
     the PacBio-collapse over-fragmentation noise: a gene with 50 PB.X.Y
     entries that all map to the same parent transcript with the same
     ORF outcome counts as 1 functional variant, not 50. Raw isoform
     row counts per gene routinely hit 200+ in oligo-dT primed
     long-read data because of bp-level TSS/TES variability; the
     distinct-functional-variant count is in the low double digits even
     for the most diverse genes on chr22 (top: RABL2B with 16).

Each section is skipped if its filtered set is empty (e.g. no NMD-sensitive
isoforms, no `parent_gene_id` column, no high/medium-confidence rows). If
all three are empty the section shows a "no notable findings" message
pointing readers at the TSV.

### `annotated.h5ad`

AnnData with isoforms in `var` (indexed by `transcript_id`), all
per-isoform annotations as `var` columns (list columns are JSON-encoded
strings for h5ad compatibility), per-cell counts in `X` when `--counts`
was supplied (otherwise `obs` is empty).

---

## Interpreting CRAFT calls (FAQ)

### "I'm working on a healthy sample. Why are any isoforms NMD-sensitive?"

Three reasons, in decreasing order of biological interest:

1. **NMD is a normal regulatory mechanism, not just a disease pathway.**
   Cells use alternative splicing to insert premature stops in some
   isoforms of normal genes, then degrade those isoforms via NMD. The
   purpose is to titrate the parent gene's protein output. This is called
   RUST (regulated unproductive splicing and translation) or AS-NMD.
   Famously, many splicing factors auto-regulate this way (SR proteins,
   hnRNPs). So a small fraction of NMD-sensitive isoforms in any healthy
   sample is biologically expected. Our chr22 smoke test showed 3.4% of
   isoforms NMD-sensitive, which is consistent with published estimates
   of 5-15% of all human transcripts being potential NMD targets.

2. **Long-read sequencing captures NMD-sensitive transcripts BEFORE
   degradation more readily than short-read RNA-seq does.** Short-read
   RNA-seq sees a degradation-modulated steady-state. Long-read protocols
   that target full-length transcripts pull in low-abundance
   NMD-sensitive species that would otherwise be hard to detect. So
   long-read data systematically over-represents NMD substrates vs the
   protein-coding "productive" pool. This is a feature, not a bug; it
   lets you measure AS-NMD regulation directly.

3. **Truncation artefacts can mimic NMD substrates.** A 3'-truncated
   read of a normal transcript can look like a transcript with a premature
   stop codon. To guard against this, CRAFT routes 3'-truncated reads to
   `STOP_NOT_OBSERVED` (the stop codon position is past the iso's last
   exon) and short-circuits the NMD module to `NOT_APPLICABLE`. So if
   your NMD-sensitive call has `nmd_confidence == "high"`, it's a
   real biological call: the stop codon was observed in the read AND
   the propagation was intact AND a downstream EJC (last junction) was
   observed AND the 50nt distance rule was violated. If
   `nmd_confidence == "medium"` (DISRUPTED outcome), the stop codon is
   observed but the iso has some structural divergence in the CDS, so the
   call is correct under the propagated stop but the propagated stop might
   not be the iso's actual in-vivo stop.

**Practical filter for biological NMD only:**
`nmd_status == "sensitive" AND nmd_confidence == "high"`. In our chr22
test this is 446 isoforms (3.4%); high-confidence biological NMD
substrates.

### "Why does CRAFT fail to find an ORF in so many isoforms? Pigeon already QC'd them."

Pigeon does *structural* QC; CRAFT does *ORF* assignment. They're
complementary, not redundant. Pigeon classifies isoforms into FSM
(Full-Length Splice Match), ISM (Incomplete Splice Match), NIC (Novel In
Catalogue), NNC (Novel Not in Catalogue), Antisense, Genic-Genomic, etc.
An iso can be pigeon-valid and still have an ORF problem:

- **ISMs are 5'-truncated by definition.** ISM means the iso's splice
  junctions are a *subset* of an annotated transcript's junctions, but
  the iso doesn't reach the annotated 5' end. Pigeon keeps these; CRAFT
  classifies them as `truncated_5p` and, depending on whether the start
  codon falls in the truncated region, `START_LOST` or
  `PROPAGATED_INTACT`. In our chr22 sample, the 17% `start_lost` rate
  is roughly the ISM rate where the truncation happens to cross the
  start codon.

- **lncRNA parents have no CDS in GENCODE.** GENCODE annotates ~17k
  long noncoding RNAs alongside ~20k protein-coding genes. A PacBio
  isoform matching a lncRNA gets classified as `propagated_intact` at
  the structural level but `NO_PARENT_CDS` at the ORF level: there's no
  CDS to propagate from. Our 18% `no_parent_cds` rate is consistent
  with the lncRNA fraction of GENCODE plus pseudogenes.

- **Novels (no_parent, 10%) include real novels and read-level
  artefacts.** Pigeon's "Novel" categories permit isoforms with
  non-canonical junctions, antisense reads, intergenic reads, etc.
  Some are real (especially in less-annotated regions); some are
  artefacts that survived pigeon's filters. CRAFT's de-novo ORF step
  finds a Met...Stop window in some of these (boosting them to LOW
  confidence) and gives up on others (NONE). Use `denovo_orf_found ==
  True` to keep the de-novo-supported orphans and drop the rest.

- **Stop-not-observed is mostly NOT technical truncation in oligo-dT data.**
  For oligo-dT primed cDNA (PacBio Iso-Seq, ONT cDNA), the read's 3' end
  *is* the polyadenylation site (the polyA tail was the priming
  substrate). True 3' truncation past the stop codon is mechanistically
  rare. CRAFT v1.1 splits this category by poly(A) signal evidence:
  isoforms with a canonical poly(A) signal in their last 50 nt are
  reclassified `stop_at_alt_polya` (alternative polyadenylation upstream
  of the canonical stop, biology); isoforms without such a signal stay
  `stop_not_observed` (rare technical artefact). Earlier versions of CRAFT
  conflated these.

**The 23% `propagated_intact` rate is the success rate, not a failure
rate.** It's the fraction of isoforms where CRAFT can hand you a
confident, parent-anchored ORF call. The other 77% break down into
descriptive categories you can filter on:

| Category               | Action                                                            |
| ---------------------- | ----------------------------------------------------------------- |
| `disrupted`            | iso has real structural changes; treat as "altered ORF"            |
| `start_lost`           | 5' truncation; could be a real ISM, decide based on sc context     |
| `stop_not_observed`    | 3' truncation; the iso may be biologically normal, just truncated  |
| `no_parent_cds`        | parent is lncRNA; iso is probably noncoding                        |
| `no_parent` + denovo   | possible novel coding iso; verify with downstream evidence         |
| `no_parent` + no denovo| likely noncoding or read artefact                                  |

### "How should I filter for high-confidence isoforms downstream?"

Depends on the question. Some recipes:

- **"Show me ORFs I can trust":** `orf_confidence in ("high", "medium")`.
- **"Show me biological NMD substrates":** `nmd_status == "sensitive"
  AND nmd_confidence == "high"`.
- **"Show me alternative 3' UTRs":** `utr3_length_delta_nt != 0` AND
  `polya_signal_motif != ""`.
- **"Drop probable artefacts":** `orf_confidence != "none" OR
  denovo_orf_found == True`.

### "Why is my isoform `propagated_intact` AND `truncated_5p` AND `start_codon_covered = True`?"

The most common single state in long-read data. The iso is missing 5'
UTR (so completeness = `truncated_5p`) but the truncation didn't cross
the start codon (so the start is observed and the full CDS got
propagated). This is exactly what the truncation-aware propagation is
designed to handle. The HIGH confidence (0.9) is justified: we're
confident in the ORF, just not in the 5' UTR.

---

## Design rationale and trade-offs

### Why DataFrame outputs instead of PyRanges everywhere

The plan's original stubs returned `pr.PyRanges` from most modules. In
practice, most CRAFT outputs are *per-isoform* (one row per transcript),
not *per-interval* (one row per exon). PyRanges with one row per
isoform loses the spatial-indexing benefits and adds friction (must
maintain dummy Chromosome/Start/End columns). All per-isoform modules
return pandas DataFrames; only inputs and intermediate splice-junction
representations stay in PyRanges.

### Why we don't filter "junk" isoforms upstream

A common temptation is to drop isoforms with no parent or no ORF. CRAFT
deliberately reports them. Reasons:

1. **Description is more useful than deletion.** A user can always
   filter on a column; they can't unfilter a row that was dropped
   silently.
2. **The boundary between "junk" and "novel" is the user's call, not
   ours.** A `no_parent` iso with a de-novo ORF in a poorly-annotated
   region might be a real novel transcript. Aggressively filtering would
   throw away discoveries.
3. **The confidence score does the filtering job better.** Downstream
   code that wants only trustworthy ORFs can filter on
   `orf_confidence`; the underlying data is still there for users who
   want it.

### Why ATG-only starts in de novo

orfipy supports TTG, CTG, GTG as alternative starts. We disabled them
in v1 because: (a) they're rare (<5% of vertebrate CDSs), (b) random
sequence has them at higher rates than ATG so the false-positive rate
goes up, (c) any real non-ATG start in a long-read iso likely has a
known ATG start in some isoform of the same gene, in which case
propagation handles it. v1.5 can add an option.

### Why minimum ORF length 50 aa

Above the noise threshold for random ATG-Stop windows. Below the typical
smORF range (smORFs are mostly 11-100 aa but the well-characterised
ones cluster around 30-80 aa). Configurable per-run.

### Why JSON-encode list columns in TSV

A column like `propagated_cds_intervals = [("chr1", 100, 200, "+"),
("chr1", 300, 400, "+")]` doesn't have a natural TSV representation.
Options: (a) explode to multiple rows (loses one-row-per-iso shape),
(b) use a custom separator like `chr1:100-200,chr1:300-400` (custom
parser needed), (c) JSON-encode (a real string parser exists, and the
data round-trips). We chose (c). For interactive grep use, JSON-encoded
lists are still legible: `grep '"propagated_cds_intervals":\[\["chr1"'`
just works.

### Why we don't use scanpy

Two reasons: (a) the only thing CRAFT needs from the scverse stack is
AnnData, which is its own much smaller package, (b) keeping scanpy out
of the dependency tree makes CRAFT installable in production pipelines
where the user has tight version constraints on scanpy itself.

### Why hmmsearch instead of hmmscan

hmmsearch loops HMMs as the outer iteration; hmmscan needs a pressed
HMM database (`hmmpress` step). For v1 we accept the slower hmmsearch
path so users don't have to press Pfam-A first. v1.5 will detect the
pressed database (look for `.h3i`/`.h3f`/`.h3m`/`.h3p` siblings) and
auto-switch to hmmscan when available.

### Why 50 bp default tolerance in completeness

PacBio TSS/TES uncertainty for full-length capture is typically 10-30 bp
from the true site (alt promoters and APA sites notwithstanding). 50 bp
absorbs that without absorbing real biology. ONT's TSS uncertainty is
larger (often >50 bp); for ONT-only datasets, consider passing
`tolerance=100` or `200`.

### Why we trust the parent for start/stop positions

The alternative is to re-find the start in the iso's sequence (look for
the most upstream ATG in frame with the parent's start). v1 doesn't do
this because the parent's start in GENCODE is curated; finding "the
most upstream ATG" in the iso opens us up to short uORFs and Kozak-poor
starts that aren't real. v1.5 can offer this as an option for users who
care about uORFs explicitly.

---

## Sequence-level ORF resolution (v1.5)

The geometric propagation above never reads the spliced sequence, so it cannot
see a frameshift, an exon-skip premature stop, or an intron retained inside the
CDS. v1.5 adds a sequence-level pass (`src/craft/core/orf/resolve.py`) that runs
alongside the geometric propagation. As of v1.7, NMD and the UTR metrics are
computed once from this resolved stop (the older geometric `nmd_status` /
geometric 3'UTR columns were removed); the geometric propagation columns
(`orf_outcome`, `propagated_cds_*`) remain as the structural classification and
the anchor the resolution builds on.

For every isoform whose parent start codon is observed, CRAFT projects the
parent start into the isoform's transcript coordinates, builds the isoform's own
spliced sequence, and walks it in 3-nt codons to the first in-frame stop. The
resulting stop drives a resolved NMD call (`nmd_status_resolved`) and a resolved
3'UTR length. Intron retention inside the CDS is detected by engulfment: a parent
intron in the CDS span that the isoform carries as continuous exonic sequence
(this does not misfire on exon skips, which also lower the junction count). The
resolved status is one of `intact`, `ptc_premature`, `ptc_intron_retained`,
`cds_extension`, `no_stop_in_read`, `resolution_failed`. uORF detection and a
long-3'UTR flag are emitted as advisory NMD branches.

5'UTR length deltas are computed symmetrically to the 3'UTR. With `--counts` and
`--group-by`, CRAFT aggregates per-isoform consequences into molecule-weighted
fractions per cell group (`per_celltype_consequence.tsv`). Every column is
documented in [`features.md`](features.md).

## Known limitations

- **Pfam scan uses hmmsearch (slow).** Switch to hmmscan against a
  pressed database is still planned.
- **No per-gene HTML track view.** Side-by-side exon/CDS/UTR views for all
  isoforms of a gene are still planned.
- **uORF and long-3'UTR NMD are advisory.** They are reported as separate flags,
  not folded into `nmd_status`, because 5' ends are frequently truncated
  in long-read data and these branches are noisier than the EJC rule.
- **Non-canonical poly(A) signals not scanned.** Eleven canonical
  variants are listed in `POLYA_SIGNALS`. Cell-type-specific or
  organism-specific extensions can be added by editing that tuple.
- **Chromosome name harmonisation not done.** All three inputs must agree
  on naming (`chr1` vs `1`). pyranges raises if the FASTA is missing a
  chromosome the GTF references.

Now implemented (were deferred in v1): intron retention inside the CDS
(`intron_retained_in_cds`), frame tracking to the real premature stop
(`resolved_orf_status`, `resolved_stop_pos`), 5'UTR analysis
(`iso_utr5_length_nt` and deltas), and cell-type-specific consequence
aggregation (`per_celltype_consequence.tsv`).
