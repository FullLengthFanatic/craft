# CRAFT user guide

How to run CRAFT in your workflow, what to watch out for, and how to use the
outputs downstream. This is the operational reference. For algorithm details
see [`methods.md`](methods.md); for a high-level overview see
[`../README.md`](../README.md).

## Contents

1. [Preparing inputs](#preparing-inputs)
2. [Running CRAFT](#running-craft)
3. [Reading the outputs](#reading-the-outputs)
4. [Common errors and fixes](#common-errors-and-fixes)
5. [Tuning parameters](#tuning-parameters)
6. [Integration with upstream tools](#integration-with-upstream-tools)
7. [Single-cell workflow](#single-cell-workflow)
8. [Performance and resources](#performance-and-resources)

---

## Preparing inputs

The three required inputs (iso GTF, reference GTF, genome FASTA) must all use
the same chromosome naming convention. The most common cause of "no isoforms
got annotated" is a `chr1` vs `1` mismatch across files. Check first:

```bash
awk '$3 == "exon" {print $1; exit}' iso.gtf
awk '$3 == "exon" {print $1; exit}' reference.gtf
head -1 genome.fa | head -c 30   # FASTA contig name on the > line
```

All three should agree.

### Iso GTF

CRAFT needs `exon` rows with a `transcript_id` attribute. That's it. The file
can be plain or gzipped (pyranges handles both). It can also be a `.gff` file
provided the attribute syntax is GTF-style (`key "value"; key2 "value2";`).

Verify before running:

```bash
awk '$3 == "exon"' iso.gtf | head -3
```

Each line should have a `transcript_id "..."` attribute. If your file is
actually GFF3 (`key=value;key2=value2`), pre-convert with `gffread`.

### Reference GTF

GENCODE basic or comprehensive both work. Ensembl works too. The file must
contain both `exon` AND `CDS` rows; CRAFT uses `exon` rows for completeness
classification and `CDS` rows for ORF propagation. If your reference only has
`exon` rows, every isoform will be classified `NO_PARENT_CDS`.

Check:

```bash
zcat gencode.v45.annotation.gtf.gz | awk '$3=="CDS"' | head -2
```

If that returns nothing, you have an annotation-only file. Download a
complete GENCODE/Ensembl release instead.

### Genome FASTA

Plain `.fa` or BGZF-compressed `.fa.gz` are both supported. For `.fa.gz` you
also need a `.gzi` index alongside the `.fai`. CRAFT auto-builds the `.fai`
if missing; for BGZF you may need to pre-index with `samtools faidx`.

The FASTA must contain every chromosome name that appears in the iso GTF and
reference GTF. Missing contigs raise during sequence extraction (de novo ORF
or poly(A) scan).

### Subsetting for development

For a fast smoke test (~1 minute), subset everything to chr22:

```bash
awk '$1 == "chr22"' iso.gtf > iso.chr22.gtf
zcat gencode.v45.annotation.gtf.gz | awk '$1 == "chr22"' > gencode.chr22.gtf
samtools faidx genome.fa chr22 > chr22.fa
samtools faidx chr22.fa
```

Avoid `head -N file.gtf.gz | zcat`. `head` reads binary bytes, not lines; the
output of that pipeline is mangled. Always `zcat first, awk second`.

---

## Running CRAFT

### Basic invocation

```bash
craft annotate \
    --isoforms  iso.gtf \
    --reference reference.gtf \
    --genome    genome.fa \
    --output-dir out/
```

Produces `out/per_isoform.tsv`, `out/per_isoform.json`, `out/report.html`,
`out/annotated.h5ad`.

### With per-cell counts

```bash
craft annotate \
    --isoforms  iso.gtf \
    --reference reference.gtf \
    --genome    genome.fa \
    --counts    counts.h5ad \
    --output-dir out/
```

The `--counts` argument accepts either a `.h5ad` file (cells in `obs`,
isoforms in `var`) or a 10x-style MTX directory (`matrix.mtx[.gz]` +
`barcodes.tsv[.gz]` + `features.tsv[.gz]` or `genes.tsv[.gz]`). The output
`annotated.h5ad` then has the cell-level counts in `X`, with the per-isoform
annotations as `var` columns.

The counts file's `var_names` must match the iso GTF's `transcript_id`s.
Isoforms in the counts that aren't in the iso GTF are dropped; isoforms in
the iso GTF that aren't in the counts get zero-filled count columns.

### With Pfam domain analysis

```bash
craft annotate \
    --isoforms  iso.gtf \
    --reference reference.gtf \
    --genome    genome.fa \
    --pfam-hmm  /path/to/Pfam-A.hmm \
    --output-dir out/
```

Pfam scanning in v1 uses `hmmsearch` (HMMs as the outer loop). On full
Pfam-A (~20k HMMs), this is slow: budget roughly an hour per few thousand
isoforms. For a typical full-genome run, plan a separate Pfam pass on the
filtered set of high-confidence ORFs:

```bash
# First run without Pfam, filter to high-confidence ORFs, then Pfam-scan only those.
craft annotate --isoforms iso.gtf --reference ref.gtf --genome genome.fa --output-dir step1/
python -c "
import pandas as pd
df = pd.read_csv('step1/per_isoform.tsv', sep='\t')
df = df[df['orf_confidence'].isin(['high', 'medium'])]
df['transcript_id'].to_csv('high_conf_tx.txt', index=False, header=False)
"
# Filter iso.gtf to the high-confidence set, then re-run with --pfam-hmm.
grep -F -f high_conf_tx.txt iso.gtf > iso.filtered.gtf
craft annotate --isoforms iso.filtered.gtf --reference ref.gtf --genome genome.fa \
    --pfam-hmm Pfam-A.hmm --output-dir step2/
```

v1.5 will detect a pressed Pfam database (`hmmpress Pfam-A.hmm`) and switch
to `hmmscan` for a substantial speedup.

---

## Reading the outputs

### TSV in pandas

```python
import json
import pandas as pd

df = pd.read_csv("out/per_isoform.tsv", sep="\t")

# List columns are JSON-encoded strings; decode when you need to use them.
for col in [
    "propagated_cds_intervals",
    "denovo_cds_intervals",
    "iso_pfam_domains",
    "parent_pfam_domains",
    "pfam_preserved",
    "pfam_lost",
    "pfam_gained",
]:
    if col in df.columns:
        df[col] = df[col].apply(lambda v: json.loads(v) if isinstance(v, str) and v else [])
```

### JSON

```python
import json
with open("out/per_isoform.json") as fh:
    records = json.load(fh)
# records is a list of dicts; list columns are real lists already.
```

### AnnData

```python
import anndata as ad
adata = ad.read_h5ad("out/annotated.h5ad")
# adata.var has all per-iso annotations (list columns are JSON strings)
# adata.X has per-cell counts (when --counts was provided)
# adata.obs has cell metadata (passed through from the input counts AnnData)
```

If you provided `--counts`, the result drops into scanpy directly:

```python
import scanpy as sc
adata = ad.read_h5ad("out/annotated.h5ad")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# Filter to high-confidence ORFs before downstream analysis:
adata = adata[:, adata.var["orf_confidence"].isin(["high", "medium"])]
```

### HTML report

Open `out/report.html` in any browser. Fully self-contained; works offline,
no CDN. Sections: summary cards, distribution bar charts, per-isoform table
(first 1000 rows; full data is in the TSV).

---

## Common errors and fixes

### `pandas.errors.ParserError: Error tokenizing data. C error: Expected 9 fields in line N, saw M`

The reference (or iso) GTF has a malformed line. Most common cause: trying
to preserve a comment header by piping a gzipped file through `head` then
`zcat`. `head` reads binary bytes from gzip streams, mangling the output.
**Fix:** always decompress first, then filter:

```bash
# wrong:
(head -5 file.gtf.gz | zcat; ...) > out.gtf
# right:
zcat file.gtf.gz | awk '$1 == "chr22"' > out.gtf
```

### Every isoform is `no_parent` or `novel_no_match`

Almost always a chromosome-name mismatch. The iso GTF says `chr22` but the
reference says `22`, or vice versa. Verify with `awk '$3=="exon" {print $1;
exit}'` on both files. If they disagree, normalise one side (typically
`sed -i 's/^chr//' file.gtf` to strip prefixes, or the reverse).

### Every isoform is `no_parent_cds`

The reference GTF has no `CDS` rows. Either you have an exon-only annotation
file, or you accidentally filtered them out during subsetting. Check with:

```bash
awk '$3=="CDS"' reference.gtf | wc -l
```

If 0, download GENCODE basic or Ensembl with CDS records.

### `KeyError: 'transcript_id'`

The iso GTF has `exon` rows but no `transcript_id` attribute on them. Some
tools (rare, but FLAIR-deprecated outputs) emit `transcript_name` instead.
Pre-process to rename the attribute, or run the tool with default attribute
naming.

### `IndexError` or `OSError` in pysam.faidx

The genome FASTA doesn't have a `.fai` next to it and CRAFT couldn't build
one. For plain `.fa`, this is usually a permissions issue (you need write
access to the FASTA's directory). For BGZF `.fa.gz`, the file must be
BGZF-compressed (not regular gzip); use `samtools faidx` once with write
permission.

### `MemoryError` or OOM on large datasets

Most likely the `_exon_overlap_bp` step in completeness classification.
PyRanges' interval join is memory-proportional to the number of candidate
parent-iso pairs. For a 50k-isoform iso GTF against full GENCODE, expect
~5-10 GB peak RAM. If you have less, run per-chromosome:

```bash
for chrom in chr{1..22} chrX chrY chrM; do
    awk -v c="$chrom" '$1 == c' iso.gtf > iso.${chrom}.gtf
    awk -v c="$chrom" '$1 == c' reference.gtf > ref.${chrom}.gtf
    samtools faidx genome.fa ${chrom} > ${chrom}.fa && samtools faidx ${chrom}.fa
    craft annotate \
        --isoforms iso.${chrom}.gtf \
        --reference ref.${chrom}.gtf \
        --genome   ${chrom}.fa \
        --output-dir out/${chrom}/
done
```

Then concatenate the TSVs (drop duplicate headers):

```bash
head -1 out/chr1/per_isoform.tsv > all.tsv
for c in chr{1..22} chrX chrY chrM; do
    tail -n +2 out/${c}/per_isoform.tsv >> all.tsv
done
```

### Pipeline runs but `report.html` looks wrong in a browser

The HTML uses plotly.js inlined in the first figure block. Some corporate
browsers strip inline scripts. Test in a different browser or save the
file locally before opening.

---

## Tuning parameters

The CLI doesn't yet expose every tunable, but the public Python API does.
Drop into a notebook and call `run_annotate` directly for parameter sweeps.

### Completeness tolerance (default 50 bp)

The slack on each end before a position is called "truncated". Increase for
ONT data (TSS uncertainty is larger), decrease if you trust your iso 5'/3'
ends precisely (e.g., post-CAGE).

```python
from craft.core.completeness import classify
classified = classify(isoforms, reference_exons, tolerance=100)
```

### Minimum de-novo ORF length (default 50 aa)

Above the noise threshold for random ATG-Stop windows. Lower (e.g., 20 aa)
if you care about smORFs.

```python
from craft.core.orf.denovo import predict
denovo = predict(orphan_isoforms, genome_fasta, min_orf_aa=20)
```

### NMD rule thresholds (defaults 50 nt / 150 nt / 400 nt)

```python
from craft.core.nmd import predict
nmd = predict(
    classified, propagated,
    ptc_threshold_nt=55,         # some labs use 55 instead of 50
    start_proximal_nt=150,
    long_last_exon_nt=400,
)
```

If you want to expose these through the CLI for a particular project, add
them in `src/craft/cli.py` as `@click.option` decorators and pass through
to `run_annotate`.

---

## Integration with upstream tools

The iso GTF file CRAFT consumes is named differently by every long-read
caller. None of them care; the file just needs `exon` rows with
`transcript_id`. Pointers per tool:

### isoseq3 + pigeon (SQANTI3-style)

Use the `<prefix>_corrected.gtf` from `pigeon` (or
`sqanti3_qc.py` if you're on the older standalone). Skip the
`IsoAnnotLite` GFF3 (that's tappAS annotations, not structure) and skip
the classification TXTs (CRAFT does its own classification).

### FLAIR

Use `<prefix>.isoforms.gtf` from `flair collapse`. FLAIR's transcript IDs
are stable across `flair quantify` runs, so you can match them to your
per-cell counts.

### IsoQuant

Use `<prefix>.transcript_models.gtf` (default mode) or
`<prefix>.extended_annotation.gtf` (extended mode that includes reference
transcripts too). CRAFT works fine with either; the extended mode will
produce `propagated_intact` for every reference transcript that's also
present, which inflates your "novel" rate as a fraction.

### Bambu

Use `extended_annotations.gtf` from `bambu::writeBambuOutput`. Bambu's
default behaviour annotates both novel and reference transcripts; if you
only want CRAFT to analyse novels, filter the bambu output first.

### FLAMES

FLAMES emits GFF3 (`isoform_annotated.gff3` or
`transcript_assembly.gff3`). `pyranges.read_gtf` handles GFF3 when the
attributes use GTF-style quoting; if you see parse errors, run the file
through `gffread -F -T file.gff3 -o file.gtf` first.

---

## Single-cell workflow

For PacBio MAS-Seq / Kinnex or ONT scNanoSeq / Curio data:

1. **Get per-cell isoform counts** from your upstream pipeline. Conventions
   vary; what CRAFT needs is an AnnData (or 10x MTX dir) where `var_names`
   are the same `transcript_id`s as in the iso GTF.

2. **Run CRAFT with `--counts`:**

   ```bash
   craft annotate \
       --isoforms  iso.gtf \
       --reference gencode.v45.annotation.gtf \
       --genome    GRCh38.fa \
       --counts    cells_x_isoforms.h5ad \
       --output-dir out/
   ```

3. **Load into scanpy and filter:**

   ```python
   import anndata as ad
   import scanpy as sc

   adata = ad.read_h5ad("out/annotated.h5ad")

   # Drop isoforms with no trustworthy ORF.
   adata = adata[:, adata.var["orf_confidence"].isin(["high", "medium"])]

   # Drop isoforms with very low total counts.
   sc.pp.filter_genes(adata, min_counts=10)

   # Normalise and log.
   sc.pp.normalize_total(adata, target_sum=1e4)
   sc.pp.log1p(adata)

   # Now you have a cells-x-isoforms AnnData with per-isoform functional
   # annotations available as adata.var columns. Standard scverse from here.
   ```

4. **Cell-type-aware isoform calls** are the next step beyond CRAFT v1.
   They're planned as a separate tool that consumes CRAFT's
   `annotated.h5ad`. Until that exists, you can do this manually in scanpy:

   ```python
   # After clustering, find isoforms whose expression differs across clusters.
   sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
   # Top hits per cluster are isoforms (with their CRAFT annotations
   # available as adata.var columns) that mark each cell type.
   ```

---

## Performance and resources

### Runtime

Rough scaling on a single core:

| Iso row count | No --pfam-hmm | With --pfam-hmm full Pfam-A |
| ------------- | ------------- | ---------------------------- |
| 13,000 (chr22)| ~1 min        | ~1 hour                      |
| 100,000       | ~10 min       | ~8 hours                     |
| 600,000 (full genome) | ~1 hour | not practical with hmmsearch |

The pipeline is single-threaded except for pyhmmer (which uses all
available cores by default during hmmsearch). If your dataset is large,
prefer per-chromosome parallel jobs.

### Memory

Peak RAM is dominated by the completeness step (PyRanges interval join).
Rough estimates:

| Iso row count | Peak RAM |
| ------------- | -------- |
| 13,000 (chr22)| ~1 GB    |
| 100,000       | ~5-8 GB  |
| 600,000       | ~15-25 GB |

If you're memory-constrained, run per-chromosome (see the OOM fix above).

### Cluster / SLURM

A simple per-chromosome SLURM array:

```bash
#!/bin/bash
#SBATCH --job-name=craft
#SBATCH --array=1-24
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00

CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
CHROM=${CHROMS[$SLURM_ARRAY_TASK_ID-1]}

awk -v c="$CHROM" '$1 == c' iso.gtf > iso.$CHROM.gtf
awk -v c="$CHROM" '$1 == c' reference.gtf > ref.$CHROM.gtf
samtools faidx genome.fa $CHROM > $CHROM.fa && samtools faidx $CHROM.fa

craft annotate \
    --isoforms iso.$CHROM.gtf \
    --reference ref.$CHROM.gtf \
    --genome $CHROM.fa \
    --output-dir out/$CHROM/

rm iso.$CHROM.gtf ref.$CHROM.gtf $CHROM.fa $CHROM.fa.fai
```

### Disk

Output sizes per 10k isoforms (approximate):

- `per_isoform.tsv`: 5 MB
- `per_isoform.json`: 20 MB
- `report.html`: 5 MB (plotly.js inlined)
- `annotated.h5ad` without counts: 4 MB
- `annotated.h5ad` with 10k cells of counts: 50-200 MB depending on sparsity

For a full-genome run on a sample with 600k isoforms and 10k cells, plan
for ~5 GB of output disk.
