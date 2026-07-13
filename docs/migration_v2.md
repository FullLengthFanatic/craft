# Migrating from CRAFT v1 to v2

v2 changes scientific semantics. Pipelines should select fields by name rather
than assuming a fixed table width.

| v1 behavior or term | v2 behavior | Why |
| --- | --- | --- |
| First in-frame ATG rescued a lost start | 5' loss is `left_censored`; inference is opt-in | Internal ATG presence does not establish initiation |
| `no_stop_in_read` mixed absence and truncation | `right_censored` plus partial CDS and non-stop flag | A missing 3' end is censoring |
| `resolved_stop_pos` | retained as legacy last-sense position; use `resolved_stop_codon_pos` | Stop geometry must use the complete three-base stop |
| UTR began after one stop base | 3' UTR begins after all three stop-codon bases | Correct transcript-coordinate definition |
| Best parent chosen mainly by shared junction count | multi-signal ranking with margin and ambiguity | Parent errors propagate into every ORF consequence |
| GTF reference reduced to exon/CDS | start/stop, phase, TSL, tags, biotype retained | Reference CDS quality affects inference strength |
| Recurrence presented as depth-stable truth proxy | descriptive detection/dispersion statistics | Expression and sampling confound recurrence |
| `nmd_status` read as biological state | legacy field retained; use susceptibility, tier, rule, limitations | Structural rules do not observe degradation |
| stopless transcripts only became not-applicable | non-stop-decay candidate can be reported separately | Distinct surveillance mechanism |
| coding CV used globally built sequence features | hexamers/scaling rebuilt in each fold | Prevent validation leakage |
| output documented as a fixed number of columns | versioned named schema | Optional evidence and comparison blocks evolve |

## New inputs

`--evidence-table` joins measured support/artifact features and builds an explicitly
uncalibrated evidence tier. `--orf-comparator-gtf` compares CRAFT with an independent
CDS caller. `--classification` now also uses unambiguous parent/gene fields as hints.

## Compatibility

Legacy fields remain where practical, but their old interpretation is deprecated.
Update report code and filters to use:

```python
df["nmd_susceptibility"]
df["nmd_evidence_tier"]
df["orf_censoring"]
df["isoform_evidence_tier"]
df["parent_ambiguous"]
```

When an evidence table is supplied, unmatched rows are `insufficient_evidence`;
they are not negative controls and should not be discarded automatically. When no
evidence table is supplied, the evidence block is empty.
