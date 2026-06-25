#!/usr/bin/env python3
"""Reconstruct the hybrid_v3 isoform GTF: PB exons from the trimmed.min20count
collapse + ENST exons from GENCODE v44, with transcript_ids matching the hybrid
count-matrix ids so --counts recurrence joins."""
import re
from pathlib import Path

BASE = Path("/data/long_read_annotation_tool/test_dataset/input_files")
HYB_IDS = BASE / "hybrid_isoform_v3/BD70_hybridv3_isoform_matrix_isoforms.txt"
PB_GTF = BASE / "20260608_FLIGHT-seq_horg_BD70_BD70.trimmed.min20count.gtf"
ENST_GTF = Path("/data/tmp/gencode.v44.annotation.gtf")
OUT = BASE / "hybrid_isoform_v3/hybrid_v3.gtf"

ids = [x.strip() for x in open(HYB_IDS) if x.strip()]
pb = {i for i in ids if i.startswith("PB")}
enst = {i for i in ids if i.startswith("ENST")}
print(f"hybrid ids: {len(ids)}  PB={len(pb)}  ENST={len(enst)}")

tid_re = re.compile(r'transcript_id "([^"]+)"')


def extract(gtf, keep, out_fh, label):
    written, seen = 0, set()
    with open(gtf) as fh:
        for line in fh:
            if "\texon\t" not in line:
                continue
            m = tid_re.search(line)
            if m and m.group(1) in keep:
                out_fh.write(line)
                written += 1
                seen.add(m.group(1))
    print(f"  {label}: {len(seen)}/{len(keep)} transcripts, {written} exon rows")
    return seen


with open(OUT, "w") as out:
    seen_pb = extract(PB_GTF, pb, out, "PB (trimmed.min20count)")
    seen_enst = extract(ENST_GTF, enst, out, "ENST (gencode v44)")

miss_pb = pb - seen_pb
miss_enst = enst - seen_enst
print(f"missing PB: {len(miss_pb)}  missing ENST: {len(miss_enst)}")
if miss_enst:
    print("  sample missing ENST:", list(miss_enst)[:3])
print(f"wrote {OUT}  ({(seen_pb and seen_enst) and len(seen_pb)+len(seen_enst)} transcripts)")
