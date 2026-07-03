# CRAFT explained: the annotated design

*What CRAFT does, why each choice was made, how it compares to existing tools, and exactly where in the code every rule and parameter lives so you can change it.*

This is the long companion to [`whitepaper.md`](whitepaper.md) (the plain-language primer) and [`features.md`](features.md) (the per-column output reference). It is written for someone who wants to read the logic, judge it, and tune it. Every stage below names the code, pastes the lines that implement the rule, and tells you whether the parameter is a command-line flag or a source-only constant.

**Code links are pinned.** Every link points at the `v1.8.0` release commit `bdca4e6`, so the line numbers do not drift as the code changes. Permalink base: `https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/`. The source tree has since advanced to v1.9.0; the v1.9.0 additions (frame-aware start rescue, `--recurrence-null` calibration, the cell-type AS-NMD output, and the `ptc_exon_length_nt` column) are described in text here and are not present at the pinned commit. Re-pin on the next tagged release.

**Numbers are traced.** The benchmark numbers (Section 2, Section 6 items 1-2) come from the committed outputs in `benchmarks/figures/`. The BD70 single-cell numbers (Section 6 items 3-4) were reproduced by running the committed analysis scripts under `test_dataset/input_files/analysis/` (the scripts are in the repo; the large count matrices they read and their multi-MB output tables are not). Anything that cannot be reproduced from those sources is marked as reported rather than verified.

---

## Contents

1. [What CRAFT is, and where it sits](#1-what-craft-is-and-where-it-sits)
2. [The central commitment: propagate, then resolve](#2-the-central-commitment-propagate-then-resolve)
3. [The pipeline, stage by stage](#3-the-pipeline-stage-by-stage)
   - [3.1 Splice junctions](#31-splice-junctions)
   - [3.2 Completeness](#32-completeness)
   - [3.3 ORF propagation](#33-orf-propagation)
   - [3.4 The poly(A) split (pipeline level)](#34-the-polya-split-pipeline-level)
   - [3.5 De-novo ORF](#35-de-novo-orf)
   - [3.6 ORF resolution](#36-orf-resolution)
   - [3.7 ORF confidence](#37-orf-confidence)
   - [3.8 NMD](#38-nmd)
   - [3.9 UTR features](#39-utr-features)
   - [3.10 Coding potential](#310-coding-potential)
   - [3.11 Pfam domains](#311-pfam-domains)
   - [3.12 Per-cell recurrence](#312-per-cell-recurrence)
   - [3.13 Export and report](#313-export-and-report)
4. [Every knob: CLI flags and source-only parameters](#4-every-knob-cli-flags-and-source-only-parameters)
5. [The 66-column output](#5-the-66-column-output)
6. [Strengths (with the evidence)](#6-strengths-with-the-evidence)
7. [Limitations, and how CRAFT compares to existing tools](#7-limitations-and-how-craft-compares-to-existing-tools)
8. [Interpreting the calls](#8-interpreting-the-calls)
9. [Implementation notes](#9-implementation-notes)
10. [Open questions](#10-open-questions)

---

## 1. What CRAFT is, and where it sits

A long-read isoform caller (isoseq+pigeon, FLAIR, IsoQuant, Bambu, FLAMES, SQANTI3) hands you two things: a GTF of isoform structures and, for single-cell protocols, a cell-by-isoform count matrix. They tell you an isoform has these exons and shows up in these cells. They do not tell you whether it still codes, whether it shifts frame, whether it is degraded by NMD, or whether it drops a functional domain. CRAFT answers those per isoform. It sits one step downstream of the caller and one step upstream of biological interpretation.

CRAFT takes three required inputs and several optional ones, and emits a 66-column per-isoform table plus an HTML report and an AnnData. The entry point is the `craft annotate` command ([`src/craft/cli.py:16`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L16-L22)); the orchestration is `run_annotate` in [`src/craft/pipeline.py`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/pipeline.py).

Two design commitments shape everything else.

**Trust the reference over de-novo guessing.** A novel isoform is rarely novel along its whole length. It usually shares most of its structure with an annotated transcript that carries a curated CDS. CRAFT projects that CDS onto the isoform rather than re-predicting an ORF from scratch. Section 2 gives the evidence for why this matters.

**Describe, don't filter.** CRAFT does not delete rows. It assumes the isoform GTF is already structurally QC'd (pigeon, SQANTI3) and reports what is there, attaching a confidence score and recurrence signals that downstream code can filter on. A row dropped inside the tool cannot be recovered downstream; a column can always be filtered. This is stated as policy in the methods reference and is the reason CRAFT keeps no-parent and no-ORF isoforms in the output.

---

## 2. The central commitment: propagate, then resolve

ORF assignment happens in two passes. **Propagation** is geometric: it projects the parent's CDS coordinates onto the isoform's exons and never reads the spliced sequence (Section 3.3). **Resolution** reads the sequence: it reconstructs the isoform's own spliced CDS and walks it codon by codon to the first in-frame stop, which is where frameshifts, exon-skip premature stops, and retained introns surface (Section 3.6).

Why propagate at all, instead of running a de-novo ORF finder on every isoform? Because de-novo prediction does badly on truncated reads, which are the modal isoform in single-cell long-read data. Two committed benchmarks measure the gap.

**Simulated truncation (`benchmarks/run_bench1.py`).** GENCODE protein-coding transcripts were truncated in silico across a grid of rates {5, 10, 25, 50}%, orientations {5', 3', both}, and three seeds, then scored only where the truncated transcript still contains the true ORF. CRAFT's propagation hits the exact start codon at 0.98 to 1.00 across every cell; orfipy de-novo plateaus at 0.94 to 0.95. CRAFT's mean absolute ORF-length error is 0 nt for every 3'-truncated cell (3' truncation does not move the start, and propagation inherits the parent's positions); orfipy sits at 8 to 12 nt because the alternative ATG it picks is usually a few codons off the true start. (Verified against `benchmarks/figures/bench1_recovery_panel.json`.)

**Real data (`benchmarks/run_bench2.py`).** On the 698,049-isoform bcM0003 FLIGHT-seq sample, 223,976 isoforms had a recoverable protein-coding-parent ORF. Stratified by completeness, CRAFT scores 1.00 in every category, which is partly tautological here (the truth is the parent projection, and CRAFT propagates from the same parent). The honest read is the orfipy bar: its start-codon exact-match rate runs 0.727 (full_length) to 0.851 (truncated_3p). Even when the true ORF is fully present in the read, de-novo prediction misses the start 15 to 28 percent of the time, which is also the per-category gap to CRAFT (14.9 to 27.3 points, since CRAFT sits at 1.00). (Verified against `benchmarks/figures/bench2_concordance_panel.json`.)

The conclusion is narrow and well-supported: reference structure is information, and re-predicting an ORF you could have inherited throws it away. De-novo prediction is therefore used only as a fallback for genuinely orphan isoforms (Section 3.5).

---

## 3. The pipeline, stage by stage

Each subsection: what the stage does, why, the code (snippet + pinned link), and how to change it.

### 3.1 Splice junctions

**Code:** [`src/craft/core/intervals.py:7-46`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/intervals.py#L7-L46)

A junction is the intron between two consecutive exons, represented as the half-open interval `[upstream_exon.End, downstream_exon.Start)`. The computation is a `groupby().shift(-1)` over exons sorted within each transcript:

```python
df = exons.df.sort_values(["transcript_id", "Start"], kind="stable").reset_index(drop=True)
df["next_start"] = df.groupby("transcript_id")["Start"].shift(-1)
df["junction_index"] = df.groupby("transcript_id").cumcount()
```

**Why this representation.** Two transcripts share a junction if and only if their `(Chromosome, Start, End, Strand)` tuples match exactly, which turns "junction overlap" into a pandas merge instead of an interval-tree intersection. That is the hot path of completeness classification against a GENCODE-scale reference (~3M junctions). Single-exon transcripts contribute zero junctions.

**To change:** structural, no parameters.

### 3.2 Completeness

**Code:** [`src/craft/core/completeness.py`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/completeness.py). Enum [`:11-20`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/completeness.py#L11-L20), parent selection [`:72-97`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/completeness.py#L72-L97), end comparison [`:100-128`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/completeness.py#L100-L128).

For each isoform CRAFT picks a parent reference transcript and classifies the isoform's structural completeness relative to it. There are seven categories: `full_length`, `truncated_5p`, `truncated_3p`, `truncated_both`, `internal_fragment`, `novel_no_match`, and `alt_3prime_end` (the last is assigned later, in Section 3.4).

**Parent selection** scores every candidate by shared junctions first, then by exon-overlap bp as a tiebreaker:

```python
sort_cols = ["iso_tx", "shared_jx", "overlap_bp"]
ascending = [True, False, False]
if prefer_coding_parent and cds_tx_ids:
    scored["ref_has_cds"] = scored["ref_tx"].isin(cds_tx_ids).astype("int64")
    sort_cols.append("ref_has_cds")
    ascending.append(False)
```

Junction sharing dominates because it is a structural assertion (same splice sites), whereas exon overlap is inflated by long constitutive exons that many isoforms share. The `--prefer-coding-parent` flag adds a lowest-priority tiebreak toward a CDS-bearing parent; it is off by default so selection stays byte-identical to v1.4.

**End comparison** is strand-aware and uses a tolerance:

```python
if strand == "+":
    five_prime_complete = iso_start <= parent_start + tolerance
    three_prime_complete = iso_end >= parent_end - tolerance
elif strand == "-":
    five_prime_complete = iso_end >= parent_end - tolerance
    three_prime_complete = iso_start <= parent_start + tolerance
...
interior_pad = max(tolerance * 2, 100)
if iso_start > parent_start + interior_pad and iso_end < parent_end - interior_pad:
    return Completeness.INTERNAL_FRAGMENT
return Completeness.TRUNCATED_BOTH
```

**Parameters.** `tolerance` defaults to 50 bp and is the `--tolerance` flag ([`cli.py:79-85`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L79-L85)). Long-read TSS/TES uncertainty for full-length capture is typically 10 to 30 bp, so 50 absorbs it without misclassifying real truncations; raise it (100 to 200) for protocols whose ends are noisier. The `interior_pad = max(tolerance * 2, 100)` that separates `internal_fragment` from `truncated_both` is **source-only** ([`completeness.py:125`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/completeness.py#L125)).

This is not SQANTI3 classification. SQANTI3's FSM/ISM/NIC/NNC describe junction-set relationships; CRAFT's completeness is end-position-focused because the next stage needs to know whether the isoform covers the parent's start and stop codon positions, which is a TSS/TES question.

### 3.3 ORF propagation

**Code:** [`src/craft/core/orf/propagation.py`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/orf/propagation.py). Enum [`:15-24`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/orf/propagation.py#L15-L24), codon positions [`:27-42`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/orf/propagation.py#L27-L42), outcome logic [`:85-92`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/orf/propagation.py#L85-L92).

CRAFT finds the parent's start and stop codon genomic positions, checks whether the isoform's exons cover each, intersects the parent's CDS intervals with the isoform's exons, and classifies the outcome by codon coverage and preserved CDS length:

```python
if not start_covered:
    outcome = ORFOutcome.START_LOST
elif not stop_covered:
    outcome = ORFOutcome.STOP_NOT_OBSERVED
elif propagated_bp == parent_bp:
    outcome = ORFOutcome.PROPAGATED_INTACT
else:
    outcome = ORFOutcome.DISRUPTED
```

The priority order matters. `start_lost` is first because without the start there is no propagation. `stop_not_observed` is next because without the stop the CDS length is undefined, so even an accidental length match is not trusted. `propagated_intact` requires both endpoints and length parity. `disrupted` is the catch-all for "endpoints observed, but the preserved CDS length changed." Isoforms with no parent get `no_parent`; isoforms whose parent has no CDS (lncRNA, pseudogene) get `no_parent_cds`.

This pass is geometric and fast but blind to frame: it never reads the sequence, so it cannot see a frameshift or a retained intron that preserves every CDS bp. That is what resolution (Section 3.6) is for.

**To change:** structural, no tunable parameters.

### 3.4 The poly(A) split (pipeline level)

**Code:** evidence gathering [`src/craft/pipeline.py:258-286`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/pipeline.py#L258-L286), reclassification [`:289-322`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/pipeline.py#L289-L322); atlas matcher [`src/craft/core/polya_atlas.py:117-173`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/polya_atlas.py#L117-L173); motif list [`src/craft/core/utr3.py:16-28`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/utr3.py#L16-L28).

This split runs at the pipeline level, not inside `classify`, because it needs the genome FASTA. For oligo-dT primed long-read cDNA, the read's 3' end is the polyadenylation site by construction, so an isoform that ends short of the parent's annotated stop is almost always alternative polyadenylation, not technical truncation. CRAFT splits the two:

```python
if row["completeness"] == Completeness.TRUNCATED_3P.value and _found(str(row["transcript_id"])):
    return Completeness.ALT_3PRIME_END.value
...
if row["orf_outcome"] == ORFOutcome.STOP_NOT_OBSERVED.value and _found(str(row["transcript_id"])):
    return ORFOutcome.STOP_AT_ALT_POLYA.value
```

Evidence comes from one of two sources, atlas first when available, motif as fallback:

```python
if atlas_index:
    chrom, iso_3p = _iso_3prime_pos(exons, strand)
    hit = match_iso_end(iso_3p, chrom, strand, atlas_index)
    if hit["matched"]:
        evidence[str(tx_id)] = {"found": True, "source": "polya_db", "pas_id": str(hit["pas_id"])}
        continue
sig = polya_near_3prime_end(exons, strand, genome)
```

**Parameters.** The atlas match window is `tolerance: int = 24` nt in `match_iso_end` ([`polya_atlas.py:117-122`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/polya_atlas.py#L117-L122)), **source-only**. The motif scan window is `window: int = 50` nt in `polya_near_3prime_end`, also **source-only**. The 11 canonical poly(A) motifs are the `POLYA_SIGNALS` tuple, **source-only**: edit that tuple to add organism- or cell-type-specific signals. Providing an atlas is the `--polya-atlas` flag ([`cli.py:57`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L57-L65)). Pre-filter the atlas by usage score (`awk '$5 >= 0.01'` for PolyASite v3.0) before passing it: the unfiltered atlas has roughly one site every 200 bp and produces uninformative ~98% hit rates. The matched site id is reported in `polya_db_site_id`, and `polya_evidence_source` records `polya_db` / `canonical_motif` / `none`.

### 3.5 De-novo ORF

**Code:** [`src/craft/core/orf/denovo.py`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/orf/denovo.py). orfipy call [`:141-148`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/orf/denovo.py#L141-L148); trigger set [`src/craft/pipeline.py:148-154`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/pipeline.py#L148-L154).

De-novo prediction runs only for orphan isoforms, defined as the three outcomes where propagation produced no usable ORF:

```python
_DENOVO_TRIGGER_OUTCOMES = frozenset(
    {
        ORFOutcome.NO_PARENT.value,
        ORFOutcome.NO_PARENT_CDS.value,
        ORFOutcome.START_LOST.value,
    }
)
```

For each, CRAFT builds the transcript sequence (reverse-complemented on the minus strand), runs orfipy, and keeps the longest ORF:

```python
candidates = list(
    orfipy_core.orfs(
        seq,
        minlen=min_bp,
        starts=["ATG"],
        strand="f",
    )
)
```

**Parameters.** Minimum ORF length is `min_orf_aa = 50` aa, i.e. `min_bp = 150` ([`denovo.py:111`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/orf/denovo.py#L108-L112) and [`:132`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/orf/denovo.py#L132)), exposed as `--min-orf-aa` ([`cli.py:108-114`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L108-L114)). 50 aa is above the noise floor for random ATG-stop windows (a random sequence hits a stop every ~21 codons) and below the typical smORF range. Start codons are **ATG-only** and **source-only** (the `starts=["ATG"]` literal): orfipy also supports TTG/CTG/GTG, but those are rare and inflate false positives, and a real non-ATG start usually has a known ATG in a sibling isoform that propagation already handles.

### 3.6 ORF resolution

**Code:** [`src/craft/core/orf/resolve.py`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/orf/resolve.py). Status enum [`:52-60`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/orf/resolve.py#L52-L60), codon walk [`:141-153`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/orf/resolve.py#L141-L153), intron-retention test [`:165-187`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/orf/resolve.py#L165-L187), classification [`:350-374`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/orf/resolve.py#L350-L374).

For every isoform whose parent start codon is observed, CRAFT projects that start into the isoform's transcript coordinates, builds the isoform's own spliced sequence, and walks it in 3-nt codons to the first in-frame stop:

```python
def _walk_to_stop(seq: str, start_tx: int) -> tuple[int, bool]:
    n = len(seq)
    i = start_tx
    while i + 3 <= n:
        if seq[i : i + 3] in _STOP_CODONS:
            return i, True
        i += 3
    return n, False
```

The resulting stop drives the `resolved_orf_status` column, one of `intact`, `ptc_premature`, `ptc_intron_retained`, `cds_extension`, `start_rescued`, `no_stop_in_read`, `resolution_failed`. Classification compares the resolved stop to the parent's:

```python
if parent_stop_tx is not None and resolved_last_coding_tx == parent_stop_tx:
    return ResolvedORFStatus.INTACT, False, True
if parent_stop_tx is not None and resolved_last_coding_tx > parent_stop_tx:
    return ResolvedORFStatus.CDS_EXTENSION, False, False
if ir:
    return ResolvedORFStatus.PTC_INTRON_RETAINED, True, False
return ResolvedORFStatus.PTC_PREMATURE, True, False
```

**Start rescue.** A `start_lost` isoform (5' truncation past the annotated start) is not surrendered to de-novo blindly. CRAFT keeps the parent CDS reading frame, anchored on the 5'-most parent-CDS base still present in the isoform, and takes the first in-frame ATG from the isoform's 5' end, translating to the first in-frame stop. The call is labelled `start_rescued`, its NMD confidence is capped at `low` (the start is inferred, not observed), and only when no in-frame ATG-to-stop exists does the isoform fall back to the de-novo longest ORF. This extends "trust the reference frame" to the truncated-start case instead of guessing a start with orfipy's longest-ORF heuristic.

**Intron retention** is detected by engulfment, not by a junction-count difference: a parent CDS-region intron that the isoform carries as one continuous exon. This is the precise test, and it does not misfire on exon skips (which also lower the junction count but are not retention):

```python
for js, je in parent_introns:
    if js < cds_lo or je > cds_hi:
        continue
    for s, e in zip(starts, ends, strict=True):
        if int(s) <= js and je <= int(e):
            return True
```

Resolution also scans the 5'UTR for upstream ORFs and emits `uorf_count` and `uorf_triggers_nmd` as advisory flags (the uORF NMD window reuses `UORF_PTC_THRESHOLD_NT = 50`, source-only at [`resolve.py:39`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/orf/resolve.py#L39)). These are advisory because 5' ends are frequently truncated in long-read data, so the branch is noisier than the EJC rule.

**Parameters.** The PTC threshold passed in is the same `--ptc-threshold-nt` used by NMD (Section 3.8). The status vocabulary and the engulfment test are source-only.

### 3.7 ORF confidence

**Code:** [`src/craft/core/orf/confidence.py:8-27`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/orf/confidence.py#L8-L27) (constants), [`:68-76`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/orf/confidence.py#L68-L76) (scoring); de-novo override [`src/craft/pipeline.py:210-222`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/pipeline.py#L210-L222).

Every call gets `high` / `medium` / `low` / `none` and a numeric score, computed as `base(outcome) * factor(completeness)`:

```python
HIGH_THRESHOLD = 0.85
MEDIUM_THRESHOLD = 0.5

_BASE_BY_OUTCOME: dict[ORFOutcome, float] = {
    ORFOutcome.PROPAGATED_INTACT: 1.0,
    ORFOutcome.STOP_AT_ALT_POLYA: 0.85,
    ORFOutcome.STOP_NOT_OBSERVED: 0.55,
    ORFOutcome.DISRUPTED: 0.45,
    ORFOutcome.START_LOST: 0.2,
}

_COMPLETENESS_FACTOR: dict[Completeness, float] = {
    Completeness.FULL_LENGTH: 1.0,
    Completeness.ALT_3PRIME_END: 1.0,
    Completeness.TRUNCATED_5P: 0.9,
    Completeness.TRUNCATED_3P: 0.9,
    Completeness.TRUNCATED_BOTH: 0.65,
    Completeness.INTERNAL_FRAGMENT: 0.5,
    Completeness.NOVEL_NO_MATCH: 0.0,
}
```

The multiplicative form enforces a natural ordering: at the same completeness, `propagated_intact` ≥ `stop_not_observed` ≥ `disrupted` ≥ `start_lost`; at the same outcome, full-length scores at least as high as any truncated variant. A 5'-truncated isoform with intact propagation stays HIGH (0.9), because an observed start codon makes the ORF unambiguous regardless of missing 5'UTR. A full-length but disrupted call drops to LOW (0.45). De-novo calls are capped at LOW by an override, since they lack a reference anchor:

```python
if category == ORFConfidence.NONE and denovo_found:
    return ORFConfidence.LOW.value, 0.25
```

**Parameters.** The two category thresholds are `--orf-high-confidence` (0.85, [`cli.py:115-121`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L115-L121)) and `--orf-medium-confidence` (0.5, [`cli.py:122-128`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L122-L128)). The base and factor tables themselves are **source-only** and are a deliberate ordering, not learned weights; edit the two dicts if you have a labelled set to calibrate against. The de-novo override score (0.25) is **source-only**.

### 3.8 NMD

**Code:** [`src/craft/core/nmd.py`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/nmd.py). Thresholds [`:28-30`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/nmd.py#L28-L30), cascade [`:117-134`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/nmd.py#L117-L134), ORF selection [`:194-220`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/nmd.py#L194-L220).

One NMD call per isoform, from the resolved stop, with a de-novo fallback for orphans. `nmd_basis` records which ORF the call used:

```python
res = resolved_by_tx.get(tx_id)
if (res is not None and str(res["resolved_orf_status"]) in _RESOLVED_WITH_STOP
        and bool(res["stop_in_transcript"]) and res["resolved_cds_intervals"]):
    ...
    basis = "resolved"
    confidence = HIGH if status == "intact" else (LOW if status == "start_rescued" else MEDIUM)
else:
    dn = denovo_by_tx.get(tx_id)
    if (dn is not None and bool(dn["denovo_orf_found"]) ...):
        ...
        basis = "denovo"
        confidence = ORFConfidence.LOW
```

The decision is a four-rule escape cascade on the resolved stop; if none fire, the isoform is `sensitive`:

```python
if in_last:
    return NMDStatus.ESCAPED, "stop_in_last_exon"
if distance <= ptc_threshold_nt:
    return NMDStatus.ESCAPED, "within_50nt_of_last_junction"
if cds_bp < start_proximal_nt:
    return NMDStatus.ESCAPED, "start_proximal"
if ptc_exon_len > long_last_exon_nt:
    return NMDStatus.ESCAPED, "long_exon"
return NMDStatus.SENSITIVE, "ptc_50nt_rule"
```

`nmd_status` is `sensitive` / `escaped` / `not_applicable`; `not_applicable` is used when there is no resolved or de-novo stop to measure from. The distance to the last junction is measured in mRNA nucleotides, not genomic distance.

**Parameters.** Three flags: `--ptc-threshold-nt` (50, the EJC 50-nt rule; Lindeboom et al. 2019, the IsoformSwitchAnalyzeR default), `--start-proximal-nt` (150, the re-initiation escape), `--long-last-exon-nt` (400, the long-exon rule). See [`cli.py:86-107`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L86-L107). The `nmd_rule` column tells you which rule fired. The long-exon rule (Lindeboom et al. 2016) is evaluated on the exon that **contains the PTC**, not the transcript's terminal exon: a long exon deposits fewer EJCs per unit length, so a PTC inside it escapes. That exon's length is reported as `ptc_exon_length_nt`; `last_exon_length_nt` is kept alongside it as descriptive context.

Two interpretation traps worth stating: `escaped` is not the same as full-length (an isoform can dodge NMD and still encode a wrecked protein), and `not_applicable` is not the same as safe (it usually means a 5'-truncated read where no reliable ORF could be placed). A baseline rate of NMD-sensitive isoforms is expected even in healthy samples, because regulated unproductive splicing (AS-NMD) is a normal mechanism.

### 3.9 UTR features

**Code:** [`src/craft/core/utr3.py`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/utr3.py). Motif scan [`:44-65`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/utr3.py#L44-L65), long-UTR constant [`:37`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/utr3.py#L37).

CRAFT measures 3' and 5' UTR lengths against the parent from the resolved ORF, and scans the 3'UTR for the strongest poly(A) signal in priority order, taking the rightmost occurrence of the first motif that hits:

```python
for motif in POLYA_SIGNALS:
    idx = upper.rfind(motif)
    if idx >= 0:
        distance = len(upper) - (idx + len(motif))
        return {"motif": motif, "distance_from_3p_end": distance}
return {"motif": "", "distance_from_3p_end": -1}
```

Priority dominates distance: if both AATAAA and ATTAAA are present, AATAAA wins regardless of which is closer to the cleavage site. The columns are `iso_utr3_length_nt`, `parent_utr3_length_nt`, `utr3_length_delta_nt`, `utr3_length_delta_pct` (and the symmetric 5'UTR set), plus `polya_signal_motif` and `polya_signal_distance_nt`.

**Parameters.** `--long-utr3-nt` (1000, [`cli.py:129-135`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L129-L135)) sets the advisory `long_utr3_triggers_nmd` flag. The `POLYA_SIGNALS` motif list is source-only. Internal-priming detection is intentionally out of scope (it lives in the sibling tool `tecap`).

### 3.10 Coding potential

**Code:** [`src/craft/core/coding_potential.py`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/coding_potential.py). Constants [`:33-48`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/coding_potential.py#L33-L48), labeling [`:383`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/coding_potential.py#L379-L385).

CRAFT builds a coding-potential model from the supplied reference itself: coding transcripts (CDS rows) versus non-coding ones (exons, no CDS) train a hexamer log-likelihood table plus three features (hexamer score, log10 ORF length, ORF coverage) and a Fickett TESTCODE statistic, fit by logistic regression. No model file ships; it self-calibrates to whatever organism the reference describes. The label is a threshold on the calibrated probability:

```python
"coding_potential_label": "coding" if score >= threshold else "noncoding",
```

The scored ORF is taken from the resolved CDS, else the propagated CDS, else the de-novo CDS, recorded in `coding_potential_orf_source`. When the reference has no non-coding transcripts to train on, the model is `None` and the label is left empty.

**Parameters.** The whole stage is on by default and toggled with `--coding-potential` / `--no-coding-potential` ([`cli.py:143-150`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L143-L150)). The 0.5 score cutoff (`threshold`, [`coding_potential.py:332`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/coding_potential.py#L329-L333)) is **source-only**, as are `MAX_TRAIN_PER_CLASS = 4000` (a cap sweep on GENCODE v45 showed AUC plateaus by 4000), the 5-fold CV, the pseudocount, and the L2 penalty. This is a screening score; confirm borderline calls with CPC2 or CPAT.

### 3.11 Pfam domains

**Code:** [`src/craft/core/pfam.py`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/pfam.py). Hit filter [`:106-110`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/pfam.py#L106-L110), set comparison [`:190-192`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/pfam.py#L190-L192).

Optional, only when `--pfam-hmm Pfam-A.hmm` is supplied ([`cli.py:50-56`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L50-L56)). CRAFT translates the isoform and parent CDS, scans both against the HMM database, and reports `pfam_preserved` / `pfam_lost` / `pfam_gained` as set differences:

```python
"pfam_preserved": sorted(iso_domains & parent_domains),
"pfam_lost": sorted(parent_domains - iso_domains),
"pfam_gained": sorted(iso_domains - parent_domains),
```

Hits are filtered by HMMer's default inclusion threshold (`hit.included`, evalue ≤ 0.01) and cached by SHA256 of the protein sequence, so a protein repeated across many cells scans once. The scan uses `hmmsearch`, which is slow on full Pfam-A (~20k HMMs); a switch to `hmmscan` against a pressed database is planned. The codon table and the inclusion threshold are source-only.

### 3.12 Per-cell recurrence

**Code:** [`src/craft/core/recurrence.py`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/recurrence.py). Recurrence [`:45-93`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/recurrence.py#L45-L93), within-gene fraction [`:96-108`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/recurrence.py#L96-L108).

Optional, only with `--counts`. CRAFT emits three columns: `total_count` (UMI-corrected molecules summed across cells), `n_cells_detected` (the count of independent cells with at least one molecule), and `isoform_fraction_within_gene`. Recurrence is computed over a called-cell whitelist when one is given:

```python
total = np.asarray(x.sum(axis=0)).ravel()
n_cells = np.asarray((x > 0).sum(axis=0)).ravel()
```

The point of `n_cells_detected` is depth-stability. A per-cell count scales with how deeply each cell was sequenced; the number of cells an isoform appears in does not. The `--cell-whitelist` flag ([`cli.py:41-49`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L41-L49)) matters because the raw matrices carry every barcode, most of which are ambient droplets; without it, recurrence is computed over every barcode and inflated. If the whitelist matches zero barcodes, CRAFT logs a note and falls back to all cells. Section 6 quantifies what this recovers.

**Calibrated recurrence (optional).** `--recurrence-null` replaces the fixed `n_cells_detected >= 3` cut with a dataset-calibrated score (`recurrence_pvalue`, `recurrence_score = 1 - pvalue`). `occupancy` scatters each isoform's molecules across cells in proportion to per-cell depth (multinomial) and takes the upper tail of the resulting Poisson-binomial occupied-cell count via its normal approximation, so it conditions on isoform total, cell count and depth; a broadly-dispersed isoform scores high, a single-cell burst scores low. `betabinom` fits an empirical beta-binomial (by moments) to the observed cells-detected counts, stratified by `structural_category` when `--classification` is supplied. Both default off (`none`), leaving the prior output unchanged. This is the first cut at the open question in Section 10, not a final calibration.

### 3.13 Export and report

**Code:** [`src/craft/export/anndata.py`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/export/anndata.py), [`src/craft/export/celltype.py`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/export/celltype.py), [`src/craft/report/html.py`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/report/html.py).

`annotated.h5ad` carries the per-isoform table as `var` (list columns JSON-encoded) and the per-cell counts in `X` when `--counts` is given. With `--counts --group-by <obs_column>`, `aggregate_consequences` writes molecule-weighted consequence fractions per cell group to `per_celltype_consequence.tsv`. The HTML report has summary cards, inlined plotly distributions (no CDN), and three small "notable findings" tables: top NMD-sensitive isoforms (filtered to `nmd_status == "sensitive"` and `nmd_confidence == "high"`), intron-retained-in-CDS isoforms, and genes with the most distinct functional variants (collapsing the PacBio over-fragmentation by counting distinct `(parent_tx_id, resolved_orf_status)` pairs).

---

## 4. Every knob: CLI flags and source-only parameters

The full option set is in [`src/craft/cli.py:16-165`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L16-L165). There is no other CLI surface; anything not listed here is changed by editing the source.

**Required inputs**

| Flag | cli.py | Meaning |
|---|---|---|
| `--isoforms` | [:17](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L17-L22) | isoform GTF |
| `--reference` | [:23](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L23-L28) | reference GTF with CDS rows |
| `--genome` | [:29](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L29-L34) | indexed genome FASTA |
| `--output-dir` | [:66](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L66-L71) | output directory |

**Optional inputs / behavior**

| Flag | cli.py | Default | Feeds |
|---|---|---|---|
| `--counts` | [:35](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L35-L40) | none | recurrence, AnnData X |
| `--cell-whitelist` | [:41](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L41-L49) | none | recurrence cell set |
| `--pfam-hmm` | [:50](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L50-L56) | none | Pfam scan |
| `--polya-atlas` | [:57](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L57-L65) | none | poly(A) split |
| `--group-by` | [:72](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L72-L78) | none | per-celltype aggregation |
| `--classification` | [:151](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L151-L159) | none | joins SQANTI columns |
| `--classification-columns` | [:160](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L160-L165) | `structural_category` | which columns to carry |
| `--prefer-coding-parent` | [:136](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L136-L142) | off | parent-selection tiebreak |
| `--coding-potential / --no-coding-potential` | [:143](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L143-L150) | on | coding-potential stage |
| `--recurrence-null` | `cli.py` (post-`bdca4e6`) | none | recurrence calibration (occupancy / betabinom) |

**Tunable numeric thresholds**

| Flag | cli.py | Default | Stage |
|---|---|---|---|
| `--tolerance` | [:79](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L79-L85) | 50 bp | completeness (3.2) |
| `--ptc-threshold-nt` | [:86](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L86-L93) | 50 nt | NMD + resolution uORF (3.6, 3.8) |
| `--start-proximal-nt` | [:94](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L94-L100) | 150 bp | NMD (3.8) |
| `--long-last-exon-nt` | [:101](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L101-L107) | 400 bp | NMD (3.8) |
| `--min-orf-aa` | [:108](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L108-L114) | 50 aa | de-novo (3.5) |
| `--orf-high-confidence` | [:115](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L115-L121) | 0.85 | confidence (3.7) |
| `--orf-medium-confidence` | [:122](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L122-L128) | 0.5 | confidence (3.7) |
| `--long-utr3-nt` | [:129](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/cli.py#L129-L135) | 1000 nt | UTR advisory (3.9) |

**Source-only parameters** (no flag; edit the source line and re-run)

| Parameter | Value | Where |
|---|---|---|
| `internal_fragment` interior pad | `max(2*tolerance, 100)` | [`completeness.py:125`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/completeness.py#L125) |
| de-novo start codons | `["ATG"]` | [`denovo.py:145`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/orf/denovo.py#L141-L148) |
| confidence base scores | dict | [`confidence.py:11-17`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/orf/confidence.py#L11-L17) |
| confidence completeness factors | dict | [`confidence.py:19-27`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/orf/confidence.py#L19-L27) |
| de-novo confidence override | 0.25 | [`pipeline.py:221`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/pipeline.py#L220-L221) |
| poly(A) atlas match window | 24 nt | [`polya_atlas.py:122`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/polya_atlas.py#L117-L122) |
| poly(A) motif scan window | 50 nt | [`utr3.py:119-124`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/utr3.py#L119-L124) |
| canonical poly(A) motifs (11) | tuple | [`utr3.py:16-28`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/utr3.py#L16-L28) |
| coding-potential cutoff | 0.5 | [`coding_potential.py:332`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/coding_potential.py#L329-L333) |
| coding-potential training cap | 4000/class | [`coding_potential.py:39`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/coding_potential.py#L39) |
| Pfam inclusion threshold | evalue ≤ 0.01 | [`pfam.py:108`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/core/pfam.py#L106-L110) |

---

## 5. The 66-column output

The canonical order is `_OUTPUT_COLUMNS` in [`src/craft/pipeline.py:57-132`](https://github.com/FullLengthFanatic/craft/blob/bdca4e6f3dcb8e6669f60fdc2033e0df5f1d664b/src/craft/pipeline.py#L57-L132); every column is defined in [`features.md`](features.md). The 66 columns group into eleven blocks:

| Block | Count | Columns |
|---|---|---|
| Completeness + parent | 8 | `transcript_id`, `completeness`, `parent_tx_id`, `parent_gene_id`, `parent_gene_name`, `shared_junctions`, `parent_overlap_bp`, `has_cds_bearing_parent` |
| ORF: propagation | 6 | `orf_outcome`, `propagated_cds_bp`, `parent_cds_bp`, `start_codon_covered`, `stop_codon_covered`, `propagated_cds_intervals` |
| ORF: resolution | 11 | `resolved_orf_status`, `resolved_stop_pos`, `resolved_cds_bp`, `resolved_aa_length`, `resolved_cds_intervals`, `ptc_introduced`, `intron_retained_in_cds`, `frame_consistent`, `stop_in_transcript`, `uorf_count`, `uorf_triggers_nmd` |
| ORF: de novo | 6 | `denovo_orf_found`, `denovo_cds_bp`, `denovo_orf_aa_length`, `denovo_start_codon`, `denovo_stop_codon`, `denovo_cds_intervals` |
| ORF confidence | 2 | `orf_confidence`, `orf_confidence_score` |
| NMD | 8 | `nmd_status`, `nmd_rule`, `nmd_confidence`, `nmd_basis`, `stop_to_last_junction_nt`, `last_exon_length_nt`, `ptc_exon_length_nt`, `long_utr3_triggers_nmd` |
| UTRs | 8 | `iso_utr3_length_nt`, `parent_utr3_length_nt`, `utr3_length_delta_nt`, `utr3_length_delta_pct`, `iso_utr5_length_nt`, `parent_utr5_length_nt`, `utr5_length_delta_nt`, `utr5_length_delta_pct` |
| Poly(A) | 4 | `polya_signal_motif`, `polya_signal_distance_nt`, `polya_evidence_source`, `polya_db_site_id` |
| Pfam | 5 | `iso_pfam_domains`, `parent_pfam_domains`, `pfam_preserved`, `pfam_lost`, `pfam_gained` |
| Coding potential | 3 | `coding_potential_score`, `coding_potential_label`, `coding_potential_orf_source` |
| Per-cell recurrence | 5 | `total_count`, `n_cells_detected`, `isoform_fraction_within_gene`, `recurrence_pvalue`, `recurrence_score` |

A typical filter for trustworthy, expressed, functionally-called isoforms:

```python
df[(df.orf_confidence.isin(["high", "medium"])) & (df.n_cells_detected >= 3)]
```

---

## 6. Strengths (with the evidence)

Every number here was reproduced from its source: the benchmark figures committed under `benchmarks/figures/`, or the committed BD70 analysis scripts under `test_dataset/input_files/analysis/` (the scripts are in the repo; the large input matrices and the output tables they produce are not). The dataset for the single-cell numbers is BD70 (FLIGHT-seq, single-cell, human organoid), GENCODE v44, recurrence over the top 3,000 called barcodes.

**1. Truncation-aware propagation is measurably better than de-novo on truncated reads.** Simulated truncation: CRAFT start-codon exact match 0.98 to 1.00 vs orfipy 0.94 to 0.95; CRAFT ORF-length error 0 nt on 3'-truncated transcripts vs orfipy 8 to 12 nt (`bench1`). Real data: orfipy misses the start 15 to 28 percent of the time on bcM0003 isoforms whose true ORF is fully present (`bench2`). See Section 2.

**2. The structural NMD labels track real decay.** On the GSE86148 UPF1-knockdown HeLa dataset, CRAFT NMD-sensitive transcripts are enriched among UPF1-KD-upregulated ones: odds ratio 1.455, one-sided Fisher p = 3.0e-15, over 47,378 eligible transcripts (CRAFT-sensitive upregulated at 16.1% vs 11.6% for escaped). For scale, GENCODE's own curated `nonsense_mediated_decay` annotation gives OR 1.415 on the same test, so a four-rule structural cascade captures essentially the same biology as hand curation. (Verified against `benchmarks/figures/bench3_enrichment.tsv`.) The effect size is modest because NMD has biochemical determinants no rule cascade captures; the benchmark answers the directionality question, not a magnitude one.

**3. Depth-stable recurrence recovers real isoforms a read filter discards.** On BD70, after UMI correction abundance and recurrence are nearly the same quantity (median 1.07 molecules per detected cell; Spearman correlation between molecules and cells 0.998), so cell recurrence is a depth-stable handle on abundance that a full-length read threshold is not. Switching the filter from read count to cell recurrence matters: lifting a read cutoff from 20 to 50 would discard ~49,000 isoforms each seen in 5 or more independent cells. Measured against GENCODE ground truth, of 71,000 annotated transcripts detected in 3+ cells only ~28,000 appear in the read-filtered (min20) catalog; the rest are real, annotated, and absent. The recovered set is between a floor of 13,722 transcripts (4,069 protein-coding, in genes the read catalog missed entirely) and a ceiling of 43,052 (17,354 protein-coding). Recovered transcripts are genuinely low-abundance (median 12 molecules vs 55 for captured transcripts at the same 3-or-more-cell recurrence), which is exactly why a read threshold never modelled them. (All reproduced from `test_dataset/input_files/analysis/`: `recurrence_v18.py`, `recurrence_pigeon_min20.py`, `recovery_gencode.py`. The floor is reproducible from the same inputs via the gene-level "no captured sibling" split, though `recovery_gencode.py` itself emits only the transcript-level ceiling.)

**4. Class-aware filtering separates signal from artefact.** On BD70, novel isoforms (NIC/NNC) are the most recurrent class (median 26 cells) and carry most of the NMD signal (17% sensitive vs 2% for FSM), so they are biology, not noise; the catch-all classes (intergenic/antisense/genic/fusion) sit at a median of 6 cells with 4% coding potential, the corner to filter hardest. Keeping FSM/ISM/novel at `n_cells_detected >= 3` and the suspect classes at `>= 5` retains 82,710 of 114,202 quantified isoforms. (Class table and retention reproduced from `recurrence_v18.py`.)

**5. Operational properties.** Pure Python, no R, no Java, no scanpy (only AnnData, a much smaller dependency). The coding-potential model self-calibrates to the reference organism with no model file. Runs at full-genome scale: the bcM0003 sample (698,049 isoforms, with `--polya-atlas` and coding potential, no Pfam) is reported in the README at roughly 1h45m and ~19 GB on a 32-core VM. (Runtime is a reported operational figure, not reproducible from a committed timing file.)

---

## 7. Limitations, and how CRAFT compares to existing tools

CRAFT is a functional-consequence annotator, not an isoform caller and not a structural-QC tool. The honest framing is where it adds something the existing tools do not, and where it deliberately does less.

| Tool | Language | Scope | ORF on truncated isoforms | Single-cell |
|---|---|---|---|---|
| **CRAFT** | Python | per-isoform ORF / NMD / UTR / coding-potential / Pfam annotation (no contrast needed) | reads start/stop off the parent when the ORF survives truncation; flags `start_lost` when it doesn't | depth-stable recurrence + per-celltype aggregation |
| IsoformSwitchAnalyzeR v2 | R | the most complete suite: 37 consequences (domains, IDR, topology, localization, signal peptide, ORF/PTC, coding potential) on top of a statistical isoform-switch (DTU) framework | ORF/PTC annotation; not specialised for the truncated-read regime | yes (v2): pseudo-bulk per cell type, then DTU |
| IsoAnnotLite / tappAS | Python / GUI | transfers annotation to isoforms | only to isoforms matching a reference transcript; genuine novels get nothing | tappAS GUI, last meaningfully updated ~2021 |
| SQANTI3 + pigeon | Python | structural QC + classification (FSM/ISM/NIC/NNC), some ORF/NMD | upstream of CRAFT; classifies, does not propagate | classification only |
| orfipy / CPAT / CPC2 / TransDecoder | Python | de-novo ORF / coding potential | de-novo, no parent anchor (the failure mode in Section 2) | n/a |

**Overlap with IsoformSwitchAnalyzeR (ISAR).** As of v2 (late 2025), ISAR explicitly handles long-read and single-cell data and annotates ORF/PTC, Pfam domain changes, and coding potential, so the per-isoform functional-annotation layer now overlaps substantially with CRAFT. ISAR is the broader and more established tool: 37 functional consequences (including IDR, protein topology, sub-cellular localization, and signal peptide), a mature differential-transcript-usage (DTU) statistical framework, and a large R/Bioconductor ecosystem. CRAFT does not try to match that breadth. Where it still differs:

- **Annotation, not switch detection.** ISAR's unit is an isoform *switch* between conditions; it needs a contrast and a DTU test. CRAFT annotates every isoform standalone, no contrast required, which is what you want for a long-read catalog where the question is "what does each isoform do."
- **Truncation-aware ORF.** CRAFT reads the parent's start/stop off the reference when the ORF survives truncation, resolves the real stop from the isoform's own sequence (catching frameshift, exon-skip, and intron-retention PTCs), and flags `start_lost` when the start is genuinely gone. ISAR's ORF/PTC step is not specialised for this regime. When the start is truncated away, CRAFT attempts a frame-aware rescue (first in-frame ATG in the parent frame, labelled `start_rescued` at low confidence) before falling back to de novo.
- **Single-cell filtering.** ISAR pseudo-bulks per cell type and runs DTU. CRAFT instead scores depth-stable per-cell recurrence and recovers isoforms a read-count filter drops, a different question (which isoforms are real) from differential usage.
- **Self-contained Python.** No R and no external annotation services (SignalP / NetSurfP / DeepTMHMM / DeepLoc); the trade-off is that CRAFT does not offer those consequence types.

The two are largely complementary: CRAFT annotates and filters a catalog per isoform, and its table can feed an ISAR switch analysis. What CRAFT deliberately does not do:

- **No structural QC.** It assumes the isoform GTF is post-QC (pigeon, SQANTI3). It describes, it does not delete.
- **No cell typing.** `--group-by` summarises consequences over an existing cell grouping; clustering is upstream.
- **No chromosome-name harmonisation.** All three inputs must agree on `chr1` vs `1`.
- **NMD is structural, not biochemical.** Four rules (Section 3.8), not EJC deposition, secondary structure, codon usage, or ribosome occupancy. The OR ~1.4 against UPF1-KD is the ceiling of what a rule cascade reaches.
- **uORF and long-3'UTR NMD are advisory only**, reported as separate flags, because 5' ends are frequently truncated in long-read data.
- **Pfam uses `hmmsearch`** (slow on full Pfam-A); `hmmscan` against a pressed database is planned.
- **Only 11 canonical poly(A) motifs** are scanned; non-canonical signals are not.
- **Coding potential is a screening score**, not a curated classifier; confirm borderline calls with CPC2/CPAT.
- **The de-novo comparator in the benchmarks is orfipy alone**, not CPAT or TransDecoder, and the bench2 CRAFT-side number is partly tautological (the truth is the parent projection). The valid claim is the orfipy gap.

---

## 8. Interpreting the calls

**Why a healthy sample still has NMD-sensitive isoforms.** Three reasons, in
decreasing order of biological interest. NMD is a normal regulatory mechanism, not
only a disease pathway: cells insert premature stops into some isoforms by
alternative splicing and degrade them to titrate the parent gene's output
(regulated unproductive splicing, RUST, or AS-NMD; splicing factors such as SR
proteins and hnRNPs auto-regulate this way), so a baseline NMD-sensitive fraction
is expected. Long-read protocols that target full-length transcripts also pull in
low-abundance NMD substrates before degradation that short-read steady-state
RNA-seq misses, so long-read data over-represents NMD substrates, which lets you
measure AS-NMD directly. And a 3'-truncated read of a normal transcript can mimic
a premature stop; CRAFT guards against this by routing reads whose stop is not
observed to `stop_not_observed` / `not_applicable` rather than calling them NMD. So
for biological NMD only, filter `nmd_status == "sensitive" and nmd_confidence ==
"high"`: the stop was observed, the ORF resolved intact, a downstream junction was
present, and the 50-nt rule was violated.

**Why many isoforms lack a confident ORF, even after pigeon/SQANTI3 QC.**
Structural QC and ORF assignment are different jobs. A pigeon-valid isoform can
still have an ORF problem:

- ISMs are 5'-truncated by definition (their junctions are a subset of an
  annotated transcript's, not reaching the 5' end). CRAFT calls them
  `truncated_5p`, and `start_lost` when the truncation crosses the start codon.
- lncRNA and pseudogene parents have no CDS in GENCODE, so a structurally clean
  match gets `no_parent_cds`: there is nothing to propagate.
- Novel isoforms (`no_parent`) mix real novels and read-level artefacts; the
  de-novo step rescues the ones with a Met-to-stop window (`denovo_orf_found ==
  True`) and gives up on the rest.
- `stop_not_observed` is mostly not technical truncation in oligo-dT data: the
  read's 3' end is the poly(A) site by construction, so a short 3' end is usually
  alternative polyadenylation, which is why the poly(A) split (Section 3.4)
  reclassifies it to `stop_at_alt_polya` when there is poly(A) evidence.

`propagated_intact` is the success rate, not a failure rate: it is the fraction
where CRAFT hands you a confident, parent-anchored ORF. The rest are descriptive
categories to filter on:

| Outcome | What to do with it |
|---|---|
| `disrupted` | real structural change; treat as an altered ORF |
| `start_lost` | 5' truncation; a frame-aware `start_rescued` ORF may be reconstructed (low confidence), else de-novo. Could be a real ISM; decide on single-cell context |
| `stop_not_observed` | 3' truncation; the isoform may be biologically normal, just clipped |
| `no_parent_cds` | parent is lncRNA/pseudogene; isoform is probably noncoding |
| `no_parent` + de-novo ORF | possible novel coding isoform; verify with other evidence |
| `no_parent`, no de-novo | likely noncoding or read artefact |

**Common filter recipes** (the same set lives in `features.md`):

```python
df[df.orf_confidence.isin(["high", "medium"])]                                   # ORFs you can trust
df[(df.nmd_status == "sensitive") & (df.nmd_confidence == "high")]               # biological NMD substrates
df[(df.coding_potential_orf_source == "denovo") & (df.coding_potential_label == "coding")]  # novel coding orphans
df[df.n_cells_detected >= 3]                                                     # recurrent isoforms (with --counts)
```

One frequent state worth naming: an isoform that is `propagated_intact` and
`truncated_5p` with `start_codon_covered == True` at HIGH confidence is the modal
long-read case. It is missing 5'UTR but not the start codon, so the ORF is
unambiguous and only the 5'UTR length is uncertain; the HIGH confidence is earned.

---

## 9. Implementation notes

A few internal choices that are not parameters but explain the shape of the code
and the output:

- **Per-isoform modules return pandas DataFrames, not PyRanges.** Most outputs are
  one row per transcript, not per interval; PyRanges would lose its
  spatial-indexing benefit and force dummy coordinate columns. Only the inputs and
  the splice-junction representation stay in PyRanges.
- **List columns are JSON-encoded in the TSV** (`propagated_cds_intervals`, Pfam
  lists, and so on) so the table stays one-row-per-isoform and round-trips;
  `per_isoform.json` keeps them as real lists.
- **No scanpy dependency.** CRAFT needs only AnnData from the scverse stack (a much
  smaller package), which keeps it installable in pipelines with tight scanpy
  version constraints.
- Per-parameter rationale (why 50 bp tolerance, 50 aa, ATG-only starts, hmmsearch,
  trusting the parent's start/stop) is inline with each stage in Section 3.

---

## 10. Open questions

These are the places where input would change the tool, drawn from the whitepaper and the design notes above.

- **A calibrated recurrence threshold.** `n_cells_detected >= 3` is a sensible default, not a principled one. Bambu sets its novelty threshold to hit a target precision against a reference; an equivalent for recurrence would beat a fixed number across samples of different depth. `--recurrence-null` (Section 3.12) is a first cut: it scores each isoform against a depth-aware occupancy null or an empirical beta-binomial. What is still missing is a target-precision calibration against a truth set (the Bambu-style step) and unique-read support (below), which would turn the score into a defensible cutoff.
- **Unique versus ambiguous support.** The recovery ceiling (43,052) is soft because best-match read assignment can credit a low-abundance isoform with a sibling's reads. A unique-read or discriminating-junction count per isoform would collapse the floor/ceiling range to a single honest number.
- **Validating the intron-retention NMD calls.** Intron-retained premature stops are the least recurrent class; the matched alternative-splicing event tables (intron-retention PSI) are an independent handle to confirm or reject them.
- **NMD beyond the rules.** The cascade is structural. uORFs and long 3'UTRs are flagged but not modelled, and sequence context (GC, secondary structure) is ignored. That is the obvious next layer if it matters for your biology.

---

*CRAFT v1.9.0, MIT-licensed. Source: https://github.com/FullLengthFanatic/craft. This document pins all code links to the v1.8.0 commit `bdca4e6` for stable line numbers; see the note at the top for the v1.9.0 additions described here but not present at that commit.*
