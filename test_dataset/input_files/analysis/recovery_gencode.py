#!/usr/bin/env python3
"""Recovery demo: does reference-guided quant + recurrence rescue real isoforms
that the read-filtered (min20) de-novo catalog drops?

GENCODE-annotated transcripts are ground-truth real and full-length. We quantify
their recurrence over the same 3000 called cells, then ask how many reproducibly
detected (n_cells>=3) annotated transcripts the pigeon min20 catalog missed.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread

BASE = Path("/data/long_read_annotation_tool/test_dataset/input_files")
G = BASE / "isoform_gencodev44"
MTX = G / "BD70_isoforms_gencode_umi_corrected_isoform_matrix.mtx"
ROWS = G / "BD70_isoforms_gencode_isoforms.txt"
COLS = G / "BD70_isoforms_gencode_cells.txt"
META = G / "BD70_isoforms_gencode_isoform_metadata.txt"
WL = BASE / "pigeon_min20_top3000_cells.txt"
SQ = BASE / "classification/BD70.trimmed.min20count_classification.txt"
HYB = BASE / "hybrid_isoform_v3/BD70_hybridv3_isoform_matrix_isoform_metadata.txt"
OUT = BASE / "analysis"


def hr(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


def base_enst(s):
    return s.astype(str).str.split(".").str[0]


# ----------------------------------------------------- gencode recurrence over called cells
hr("1. GENCODE v44 recurrence over the 3000 called cells")
ensts = pd.read_csv(ROWS, header=None)[0].to_numpy()
barcodes = pd.read_csv(COLS, header=None)[0].to_numpy()
whitelist = set(pd.read_csv(WL, header=None)[0])
mask = np.fromiter((b in whitelist for b in barcodes), bool, len(barcodes))
print(f"annotated transcripts={len(ensts)}  cells in matrix={len(barcodes)}  whitelisted={int(mask.sum())}")
m = mmread(str(MTX)).tocsc()[:, mask]
total = np.asarray(m.sum(axis=1)).ravel()
ncells = np.asarray((m > 0).sum(axis=1)).ravel()
df = pd.DataFrame({"enst": ensts, "total_count": total.astype(int), "n_cells": ncells.astype(int)})
df["base"] = base_enst(df["enst"])

meta = pd.read_csv(META, sep="\t")[["transcript_id", "transcript_type"]]
meta["base"] = base_enst(meta["transcript_id"])
df = df.merge(meta[["base", "transcript_type"]].drop_duplicates("base"), on="base", how="left")
df["is_pc"] = df["transcript_type"].eq("protein_coding")

detected = df["n_cells"] >= 1
recurrent = df["n_cells"] >= 3
print(f"detected (>=1 cell): {int(detected.sum())}  recurrent (>=3 cells): {int(recurrent.sum())}")
print(f"  protein_coding detected: {int((detected & df['is_pc']).sum())}  "
      f"recurrent: {int((recurrent & df['is_pc']).sum())}")
print(f"molecules per cell (median over detected): "
      f"{(df.loc[detected,'total_count']/df.loc[detected,'n_cells']).median():.2f}")

# ----------------------------------------------------- rare-but-real population
hr("2. Rare-but-real: real (annotated) isoforms a count cut would drop")
rare_real = (df["total_count"] < 20) & (df["n_cells"] >= 3)
print(f"annotated transcripts with total_count<20 BUT detected in >=3 cells: {int(rare_real.sum())}")
print(f"  of which protein_coding (full-length real): {int((rare_real & df['is_pc']).sum())}")
print(f"  their median n_cells={int(df.loc[rare_real,'n_cells'].median())}, "
      f"median total_count={int(df.loc[rare_real,'total_count'].median())}")

# ----------------------------------------------------- vs pigeon min20 catalog
hr("3. RECOVERY vs the pigeon min20 read-filtered catalog")
sq = pd.read_csv(SQ, sep="\t", usecols=["associated_transcript"], low_memory=False)
captured = set(base_enst(sq.loc[sq["associated_transcript"].astype(str).str.startswith("ENST"),
                                "associated_transcript"]))
print(f"distinct ENSTs captured by the pigeon min20 catalog (associated_transcript): {len(captured)}")
df["captured"] = df["base"].isin(captured)
recovered = recurrent & ~df["captured"]
print(f"\nreproducibly-detected (>=3 cells) annotated transcripts: {int(recurrent.sum())}")
print(f"  captured by pigeon min20:        {int((recurrent & df['captured']).sum())}")
print(f"  RECOVERED (missed by min20):     {int(recovered.sum())}  "
      f"({100*recovered.sum()/recurrent.sum():.1f}% of recurrent annotated)")
print(f"  of recovered, protein_coding (full-length): {int((recovered & df['is_pc']).sum())}")
print(f"  recovered abundance: median n_cells={int(df.loc[recovered,'n_cells'].median())}, "
      f"median total_count={int(df.loc[recovered,'total_count'].median())}  "
      f"(low: that is why a read filter missed them)")
print("\nrecovered by transcript_type (top):")
print(df.loc[recovered, "transcript_type"].value_counts().head(6).to_string())

# ----------------------------------------------------- hybrid framing
hr("4. Hybrid set = the practical embodiment (keep all reference, add novel)")
hyb = pd.read_csv(HYB, sep="\t")
nref = int((hyb["structural_category"] == "reference").sum())
nnov = int((hyb["structural_category"] != "reference").sum())
print(f"hybrid_v3: {len(hyb)} isoforms = {nref} reference (annotated) + {nnov} novel/other")
print("By construction the hybrid keeps every annotated transcript, so the recovered set above")
print("is included; the novel fraction is what you filter on recurrence (the v18 analysis).")

# ----------------------------------------------------- save
hr("5. SAVE")
df.to_csv(OUT / "gencode_recovery.tsv", sep="\t", index=False)
print(f"wrote {OUT/'gencode_recovery.tsv'}  ({len(df)} transcripts)")
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    ax[0].scatter(df.loc[detected, "total_count"].clip(1, 1e4),
                  df.loc[detected, "n_cells"].clip(1, 3000),
                  s=3, alpha=0.05,
                  c=np.where(df.loc[detected, "captured"], "#9aa0a6", "#d1495b"))
    ax[0].axhline(3, color="k", ls=":", lw=0.8); ax[0].axvline(20, color="k", ls=":", lw=0.8)
    ax[0].set(xscale="log", yscale="log", xlabel="total molecules", ylabel="n_cells",
              title="annotated tx (red=missed by min20)")
    rec = df[recovered]
    ax[1].hist(rec["n_cells"].clip(3, 200), bins=np.logspace(np.log10(3), np.log10(200), 30),
               color="#d1495b")
    ax[1].set(xscale="log", xlabel="n_cells (recovered)", ylabel="transcripts",
              title=f"recovered annotated tx (n={int(recovered.sum())})")
    fig.tight_layout(); fig.savefig(OUT / "gencode_recovery.png", dpi=110)
    print(f"wrote {OUT/'gencode_recovery.png'}")
except Exception as e:
    print(f"(plots skipped: {e})")
print("\nDONE")
