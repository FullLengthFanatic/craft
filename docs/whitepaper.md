# CRAFT: functional consequences for long-read isoforms, and how we filter them

*A primer for collaborators. Tells you the logic so you can judge it and push back.*

## Plain-language summary

A gene can be read out as many different RNA versions, built by joining its pieces in different combinations. Long-read sequencing now reads each version from end to end, and a single experiment can turn up hundreds of thousands of them. Knowing a version exists is the easy part. The useful question is what it does: does it still make a working protein, a shortened one, or none at all, and is the cell flagging it as faulty and degrading it. CRAFT works that out for every version, by using what is already known about the gene instead of guessing each one from scratch.

The second half of this note is a quieter problem: deciding which versions are real enough to keep. The common approach keeps the ones supported by enough sequencing reads. That sounds safe but it is biased, because the read count depends on how deeply you happened to sequence, and it discards rare versions that are genuinely present. A fairer test is how many separate cells a version shows up in. One seen loudly in a single cell is weaker evidence than one seen quietly across many. Switching to that test recovers tens of thousands of real, full-length transcripts that a read-count cutoff throws away, which is the change this work makes.

---

Long-read sequencing gives you isoform structures. It does not tell you what those structures do to the protein. CRAFT closes that gap: for every isoform in a long-read catalog it works out the reading frame, whether the protein stays intact or is truncated, whether the transcript is an NMD target, what happens to the UTRs, and which Pfam domains survive. This note covers three things, in order of how much recent work they took: how CRAFT assigns an ORF, how it calls NMD, and how we filter a single-cell isoform set.

The worked examples use **BD70**, a FLIGHT-seq single-cell sample from a human organoid, quantified isoform by cell.

---

## 1. The problem

A long-read caller (pigeon/isoseq, FLAIR, IsoQuant, Bambu) hands you two things: a GTF of isoform structures and a cell-by-isoform count matrix. Between them they say "this isoform has these exons and shows up in these cells." They do not say whether it still codes, whether it shifts frame, whether it is degraded by NMD, or whether it drops a functional domain. Those are the questions that tie an isoform to biology, and most single-cell long-read pipelines stop just short of them.

CRAFT answers them per isoform, with one design commitment: when an isoform overlaps a known gene, trust the reference instead of guessing. A novel isoform is rarely novel along its whole length. It usually shares most of its structure with an annotated transcript, and that annotation carries a curated CDS. CRAFT uses it.

---

## 2. Assigning an ORF: propagate, then resolve

Two steps, in that order.

**Propagation.** CRAFT picks the best-matching reference transcript for an isoform by maximal splice-junction sharing, classifies the isoform's structural completeness from its ends (`full_length`, `truncated_5p`, `alt_3prime_end`, `internal_fragment`, and so on), then projects the parent's CDS coordinates onto the isoform's exons. It records whether the start and stop codons are actually inside the read or fall off the end. This step is geometric and fast, but it never re-reads the sequence, so on its own it cannot see a frame change.

**Resolution.** This is the part that reads sequence. From the projected start, CRAFT reconstructs the isoform's *own* spliced CDS and walks it codon by codon to the first in-frame stop. Frameshifts, exon-skip premature stops, and retained introns all surface here, because they move the real stop. The result is `resolved_orf_status`, one of `intact`, `ptc_premature`, `ptc_intron_retained`, `cds_extension`, `no_stop_in_read`, `resolution_failed`.

**De novo fallback.** Genuinely orphan isoforms, with no reference parent, get a de-novo ORF and a coding-potential score (a hexamer model plus ORF-length and Fickett features, self-calibrated on the reference's own coding and non-coding transcripts). This is the only place CRAFT guesses, and it labels the guess as such.

**Confidence.** Every call gets `high` / `medium` / `low`, tied to how much of the ORF was observed in the read rather than inferred. A start codon you can see is worth more than one you assume.

Why propagate rather than predict from scratch: on simulated truncations, propagation hits the true start codon 98 to 100 percent of the time, against about 95 percent for de-novo prediction, and the ORF-length error is 0 nt for 3'-truncated transcripts versus 8 to 12 nt for de-novo. Reference structure is information; throwing it away to re-predict is a step backward.

---

## 3. Calling NMD

NMD prediction is a rule cascade on the resolved stop. The core is the 50-nucleotide rule: if translation ends more than 50 nt upstream of the last exon-exon junction, the ribosome stalls with a downstream exon-junction complex still in place, and the transcript is a target. Then come the escapes, each a known way a premature-looking stop avoids decay:

- the stop sits in the last exon,
- the stop is within 50 nt of the last junction,
- the CDS is start-proximal (under ~150 bp, so reinitiation rescues it),
- the last exon is long (over ~400 bp).

The output is `nmd_status`: `sensitive`, `escaped`, or `not_applicable`. A companion column, `nmd_basis`, says where the call came from: the resolved ORF, the de-novo ORF (for orphans), or neither.

Two traps, stated plainly because they catch people:

- **`escaped` is not `full-length`.** An isoform can dodge NMD and still encode a wrecked protein. Escape is about degradation, not function.
- **`not_applicable` is not `safe`.** It means CRAFT could not place a reliable ORF, usually a 5'-truncated read, not that the transcript is fine.

Does the structural rule track real decay? On a UPF1-knockdown bulk RNA-seq dataset, transcripts CRAFT calls NMD-sensitive are enriched among the up-regulated ones (odds ratio 1.46, one-sided p = 3e-15). The cascade is structural rather than biochemical, so the effect is modest, but it points the right way.

---

## 4. Filtering single-cell isoforms

This is the part worth scrutinising, because it is where a read-count threshold quietly does the wrong thing.

### Why a read count is the wrong currency

A raw read count is, roughly, abundance times depth times capture. How deeply a sample is sequenced and how many cells it carries vary from run to run, so the depth term is a property of the experiment, not of the isoform. A fixed cutoff, say "keep isoforms with at least 20 reads," is therefore harsher on a shallow sample and softer on a deep one, and it preferentially deletes rare-but-real isoforms of lowly-expressed genes. Concretely, on BD70, lifting the read cutoff from 20 to 50 would discard 49,000 isoforms that are each detected in 5 or more independent cells.

### What to use instead: recurrence

Count how many cells an isoform appears in. An isoform seen in 4 separate cells is stronger evidence than one confined to a single cell, however high its molecule count there: four independent observations beat one. Recurrence carries no depth term, so a threshold on it transfers from one sample to the next. It is the single-cell analogue of requiring detection across replicates.

One fact made this simpler than expected. After UMI correction, abundance and recurrence are essentially the same quantity: on BD70 the median is 1.07 molecules per detected cell, and the rank correlation between total molecules and number of cells is 0.998. So there is really one number to filter on, and it is depth-stable. The full-length read count that a min-reads filter (pigeon's `min20`) thresholds on is depth-dependent in a way recurrence is not: it follows how deeply the sample was sequenced and how efficiently each molecule was captured, not the number of cells the isoform was actually seen in.

CRAFT now writes three columns from the count matrix:

- `total_count`: UMI-corrected molecules summed across cells.
- `n_cells_detected`: the recurrence signal.
- `isoform_fraction_within_gene`: an isoform's share of its gene's molecules, which cancels depth as well.

It computes them over a called-cell whitelist (`--cell-whitelist`). This matters: the raw matrices carry every barcode, and on BD70 roughly 140,000 of the 144,000 barcodes are ambient droplets, not cells. Counting recurrence over all of them would inflate every isoform.

### One cut does not fit every class

Structural classes earn different amounts of trust. Reference-matching isoforms (FSM) are corroborated by the annotation, so they can be kept at low recurrence. Novel isoforms (NIC/NNC) have to prove themselves. The catch-all classes, intergenic, antisense, genic, fusion, are the suspect corner. BD70 makes the pattern obvious:

| class | n | median cells | coding | NMD-sensitive |
|---|---|---|---|---|
| FSM | 24,433 | 22 | 62% | 2% |
| ISM | 36,427 | 9 | 51% | 4% |
| novel (NIC/NNC) | 27,184 | 26 | 61% | 17% |
| intergenic/antisense/genic/fusion | 26,151 | 6 | 4% | 1% |

Two things to read off it. The novel isoforms are the *most* recurrent class and carry most of the NMD signal (17 percent sensitive against 2 percent for FSM), so they are real biology, not noise. The suspect classes sit at a median of 6 cells with 4 percent coding potential, which is the corner to filter hardest.

The rule we settled on: keep FSM, ISM, and novel at `n_cells_detected >= 3`; require the suspect classes to clear `n_cells_detected >= 5` and to be free of intra-priming and non-canonical junctions. On BD70 that retains 82,710 of the 114,202 quantified isoforms, dropping the single-cell and artefact-looking tail without touching abundance.

---

## 5. Recovering the isoforms a read filter drops

A read-filtered de-novo catalog cannot show you what it deleted, so we measured it against ground truth. GENCODE-annotated transcripts are real by definition, so we quantified all of them in the same cells and asked which reproducibly-detected ones the read-filtered catalog is missing.

Of 71,000 annotated transcripts detected in 3 or more BD70 cells, only about 28,000 appear in the read-filtered (min20) catalog. The rest are real, annotated, full-length, and absent.

Stated with the honest bounds:

- **Floor: 13,722 recovered transcripts (4,069 protein-coding).** These sit in genes the read catalog missed entirely, so no sibling isoform can be confused for them. This is the number you can stand behind without caveat.
- **Ceiling: 43,052 (17,354 protein-coding).** This adds transcripts in genes where a sibling *was* captured, where the matrix's best-match read assignment could credit some molecules to the wrong isoform of the gene. Real transcripts, but their per-isoform support is softer.

Either way, the recovered transcripts are genuinely low-abundance (a median of 12 molecules, against 55 for captured transcripts at the same 3-or-more-cell recurrence), which is precisely why a read threshold never built a model for them. Low abundance is not low confidence once you have recurrence.

The operational answer is the **hybrid catalog**: keep every annotated transcript by construction (reference-guided quantification, no read threshold on known isoforms), add the novel isoforms from the de-novo collapse, and apply the recurrence filter only to the novel fraction. Reference for what is known, recurrence for what is new.

---

## 6. What CRAFT does not do

It does not do structural QC. It assumes the isoform GTF is already curated (pigeon, SQANTI3) and describes what is there rather than deleting rows, because a row dropped inside the tool cannot be recovered downstream. It does not call cell types. It does not harmonise chromosome naming between inputs. Filtering is a decision it hands you the signals for, not one it makes for you.

---

## 7. The columns you will actually read

| column | meaning |
|---|---|
| `completeness` | structural class vs the parent (full_length, truncated_5p, alt_3prime_end, ...) |
| `resolved_orf_status` | sequence-level ORF outcome (intact, ptc_premature, ptc_intron_retained, ...) |
| `orf_confidence` | high / medium / low, by how much of the ORF was observed |
| `nmd_status` | sensitive / escaped / not_applicable |
| `nmd_basis` | resolved / denovo / none (which ORF the NMD call used) |
| `coding_potential_label` | coding / noncoding, from the self-calibrated model |
| `total_count` | UMI molecules across called cells (needs `--counts`) |
| `n_cells_detected` | recurrence: cells the isoform is seen in (the depth-stable filter) |
| `isoform_fraction_within_gene` | the isoform's share of its gene's molecules |

A typical filter for trustworthy, expressed, functionally-called isoforms:

```python
df[(df.orf_confidence.isin(["high", "medium"])) & (df.n_cells_detected >= 3)]
```

---

## 8. Open questions, where your input would help

- **A calibrated recurrence threshold.** `n_cells_detected >= 3` is a sensible default, not a principled one. Bambu sets its novelty threshold to hit a target precision against a reference; an equivalent for recurrence would beat a fixed number, especially across samples of different depth.
- **Unique versus ambiguous support.** The recovery ceiling (43k) is soft because best-match assignment can credit a low-abundance isoform with a sibling's reads. A unique-read or discriminating-junction count per isoform would let us report a single honest recovery number instead of a range.
- **Checking the intron-retention NMD calls.** Intron-retained premature stops are the least recurrent class (a median of 5 cells). The matched alternative-splicing event tables (intron-retention PSI) are an independent handle to validate or reject them.
- **NMD beyond the rules.** The 50-nt cascade is structural. uORFs and long 3'UTRs are flagged but not modelled, and sequence context (GC, secondary structure) is ignored. If that matters for your biology, it is the obvious next layer.

---

## Appendix: the BD70 run

Reference GENCODE v44, genome GRCh38 primary assembly. Isoforms from the trimmed (poly(A)-clipped) pigeon collapse, filtered to >=20 full-length reads. Recurrence computed over the top 3,000 barcodes from the cell-calling knee. CRAFT v1.8.0. The per-isoform table, the recurrence and recovery analyses, and the scripts that produced every number above live under `test_dataset/` and `test_dataset/input_files/analysis/`; all figures are reproducible from the saved tables.
