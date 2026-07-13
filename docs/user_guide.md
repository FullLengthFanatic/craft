# CRAFT v2 user guide

## Prepare inputs

The isoform GTF needs exon rows and `transcript_id`. The reference should include
exon, CDS, `start_codon`, and `stop_codon` rows where available. Do not strip phase,
tags, transcript support level, biotype, MANE/APPRIS, or CCDS attributes. GTF and
FASTA must use the same assembly and contig names.

```bash
awk '$3 == "exon" {print $1; exit}' isoforms.gtf
awk '$3 == "CDS" {print $1; exit}' reference.gtf
head -1 genome.fa
```

The FASTA may be indexed with `samtools faidx`; CRAFT attempts to create a missing
`.fai` where possible.

## Basic run

```bash
craft annotate \
  --isoforms isoforms.gtf \
  --reference reference.gtf \
  --genome genome.fa \
  --output-dir out
```

Add a Pigeon or SQANTI3 table to retain structural category and improve parent/gene
assignment:

```bash
craft annotate \
  --isoforms isoforms.gtf \
  --reference reference.gtf \
  --genome genome.fa \
  --classification classification.txt \
  --classification-columns structural_category,associated_transcript,associated_gene \
  --output-dir out
```

Multi-valued upstream parent hints are not silently resolved; CRAFT falls back to
its candidate ranking and reports ambiguity.

## Molecule/read evidence

Pass a TSV/CSV keyed by `transcript_id`, `isoform`, `pbid`, or `id`:

```bash
craft annotate ... --evidence-table isoform_evidence.tsv
```

Recognized fields include fractions (0–1 or percentages) or common aliases for
unique/ambiguous molecules, canonical and short-read-supported junctions,
full-length/5'-adapter/poly(A) evidence, mapping quality, replicate count, internal
priming, template switching/strand invasion, and chimeras. See `features.md` for
canonical names.

The resulting evidence score is uncalibrated. A missing field is ignored rather
than treated as zero. Fewer than three measured components yields
`insufficient_evidence`.

## Counts and cell groups

`--counts` accepts `.h5ad` (cells in `obs`, isoforms in `var`) or a 10x MTX
directory. Feature names must match GTF transcript IDs.

```bash
craft annotate ... \
  --counts counts.h5ad \
  --cell-whitelist called_cells.txt \
  --group-by cell_type
```

Use a called-cell whitelist to avoid treating ambient/empty barcodes as biological
replicates. Counts add abundance and detection summaries. `--recurrence-null
occupancy` or `betabinom` adds exploratory dispersion statistics; do not interpret
them as isoform-validity probabilities.

Group output is molecule weighted. If independent evidence tiers are supplied,
AS–NMD candidate tables retain `strong` and `supported` structures; otherwise the
table is expression-only and should be described that way.

## Alternative ORF comparison

Run an independent reference-aware caller on the same transcript models and pass
its CDS-annotated GTF:

```bash
craft annotate ... --orf-comparator-gtf orfanage.gtf
```

CRAFT reports whether a comparator ORF exists and whether exact genomic start,
stop-codon anchor, and CDS length agree. Caller agreement is supporting evidence,
not a truth label.

## Truncated transcripts

By default, a lost 5' CDS boundary is `left_censored`. If a hypothesis about the
first downstream in-frame ATG is useful, opt in:

```bash
craft annotate ... --infer-alternative-start
```

Always retain `alternative_start_inferred` in exported results. A right-truncated
molecule is `right_censored`; CRAFT preserves the partial CDS and may flag a
non-stop-decay candidate.

## Poly(A), Pfam, and coding potential

```bash
craft annotate ... \
  --polya-atlas filtered_polya_sites.bed \
  --pfam-hmm Pfam-A.hmm
```

Use a tissue-appropriate, quality-filtered poly(A) atlas; overly dense unfiltered
atlases make matches uninformative. Full Pfam scanning can be slow and is often
best run on a biologically selected subset. Coding potential is on by default and
can be disabled with `--no-coding-potential`.

## Recommended review workflow

1. Check input/reference match rates and parent ambiguity.
2. Separate complete ORFs from left/right-censored partial CDS.
3. Inspect raw evidence components and artifact warnings before selecting a tier.
4. Treat recurrence as expression/sampling context.
5. Review independent ORF disagreement and weak reference CDS metadata.
6. Report surveillance results as predicted susceptibility and include limitations.
7. Calibrate final filtering thresholds on controls from the same assay and pipeline.

## Common problems

| Symptom | Likely cause | Action |
| --- | --- | --- |
| Almost no parents | assembly/contig mismatch | compare GTF and FASTA contigs |
| All parents lack CDS | incomplete reference GTF | use the full annotation with CDS rows |
| Many ambiguous parents | repetitive/short structures or weak hints | import upstream gene/parent fields and inspect margins |
| Many left-censored ORFs | 5' truncation in the assay | do not reinterpret internal ATGs as known starts |
| Evidence fields empty | transcript IDs or aliases do not match | inspect table key and canonical names |
| Implausibly many poly(A) matches | atlas too dense | filter atlas by usage/quality and tissue |
| Counts do not join | count feature IDs differ from GTF IDs | normalize IDs before running |

Run `craft annotate --help` for thresholds and [`features.md`](features.md) for the
complete output contract.
