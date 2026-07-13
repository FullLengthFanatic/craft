# CRAFT v2: the quick scientific view

Long-read single-cell RNA sequencing exposes transcript structures that short reads
cannot assemble reliably, but it also exposes assay artifacts and incomplete
molecules. A single “minimum reads” or “seen in N cells” cutoff cannot solve that
problem: expression, capture, mapping ambiguity, splice support, molecular ends,
and artifact mechanisms are different dimensions.

CRAFT v2 therefore annotates four layers independently:

1. **Structure evidence.** Unique versus ambiguous molecules, canonical and
   short-read-supported junctions, 5' and poly(A) evidence, replicate support,
   mapping quality, and artifact signals are retained as auditable features.
2. **Completeness.** Missing transcript ends are represented as left or right
   censoring. CRAFT does not pretend that the first observed ATG or last observed
   codon is necessarily a biological CDS boundary.
3. **ORF consequence.** A carefully ranked reference parent supplies a CDS frame
   when justified; observed sequence then reveals premature stops, extensions,
   exon skipping, and retained introns. True orphans receive a labelled de novo
   candidate. Independent ORF callers can be compared without declaring a winner.
4. **RNA surveillance.** Transparent structural rules report predicted NMD or
   non-stop-decay susceptibility, an evidence tier, and limitations. They do not
   claim that degradation was observed.

Reference annotation helps, but it is not ground truth. CRAFT retains CDS
completeness, phase, explicit start/stop features, MANE/APPRIS/CCDS/basic tags, and
parent ambiguity so users can see when propagation is weak.

The v2 evidence score is deliberately simple and uncalibrated. It is useful for
ranking and triage, not as a universal probability that an isoform is real. A
publishable FLIGHT-seq filter needs assay-specific labelled controls, sample- and
gene-held-out validation, comparison with respected callers/QC tools, calibration
plots, and a frozen external test set. The code preserves the raw evidence needed
to build that model later.

Read [`v2_design.md`](v2_design.md) for the full scientific contract,
[`features.md`](features.md) for fields, and [`user_guide.md`](user_guide.md) to run it.
