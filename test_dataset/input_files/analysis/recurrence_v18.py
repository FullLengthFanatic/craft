#!/usr/bin/env python3
"""Structure x recurrence analysis on the shipped CRAFT v1.8 BD70 output.

Reads craft_out_v18_counts/per_isoform.tsv (recurrence columns now produced by the
pipeline over the 3000-cell whitelist) and joins SQANTI structural_category +
artifact signals + FL reads. Adds CRAFT's own functional consequence (resolved ORF,
NMD, coding potential) to the recurrence picture, which the standalone v17 analysis
could not do in one table.
"""
from pathlib import Path

import numpy as np
import pandas as pd

V18 = Path(
    "/data/long_read_annotation_tool/test_dataset/BD70_min20count/"
    "craft_out_v18_counts/per_isoform.tsv"
)
SQANTI = Path(
    "/data/long_read_annotation_tool/test_dataset/input_files/classification/"
    "BD70.trimmed.min20count_classification.txt"
)
OUT = Path("/data/long_read_annotation_tool/test_dataset/input_files/analysis")


def hr(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


craft_cols = [
    "transcript_id", "completeness", "resolved_orf_status", "orf_confidence",
    "coding_potential_label", "nmd_status", "nmd_confidence", "parent_gene_id",
    "total_count", "n_cells_detected", "isoform_fraction_within_gene",
]
df = pd.read_csv(V18, sep="\t", usecols=craft_cols, low_memory=False)

sq = pd.read_csv(SQANTI, sep="\t", low_memory=False).rename(columns={"isoform": "transcript_id"})
sq["perc_A"] = pd.to_numeric(sq["perc_A_downstream_TTS"], errors="coerce")
sq["FL_reads"] = pd.to_numeric(sq["FL.sample"], errors="coerce")
df = df.merge(
    sq[["transcript_id", "structural_category", "perc_A", "RTS_stage", "all_canonical", "FL_reads"]],
    on="transcript_id", how="left",
)
df["intrapriming"] = df["perc_A"] >= 60
df["noncanonical"] = df["all_canonical"].astype(str).str.lower().eq("noncanonical")

fsm = df["structural_category"].eq("full-splice_match")
novel = df["structural_category"].isin(["novel_not_in_catalog", "novel_in_catalog"])
junk = df["structural_category"].isin(
    ["intergenic", "antisense", "genic", "fusion", "genic_intron"]
)
measured = df["total_count"].notna()

hr("0. SET + COVERAGE")
print(f"CRAFT isoforms: {len(df)}")
print(f"measured (in UMI matrix, recurrence present): {int(measured.sum())} "
      f"({100*measured.mean():.1f}%)")
print(f"molecule-less (FL reads >=20 but 0 UMI molecules in matrix): {int((~measured).sum())} "
      f"({100*(~measured).mean():.1f}%)")

hr("1. molecules ~= n_cells (depth currency check, from shipped columns)")
r = df.loc[measured, "total_count"] / df.loc[measured, "n_cells_detected"].replace(0, np.nan)
print(f"molecules per detected cell: median {r.median():.2f}  p95 {r.quantile(.95):.2f}")
print(f"spearman(total_count, n_cells): "
      f"{df.loc[measured, ['total_count','n_cells_detected']].corr(method='spearman').iloc[0,1]:.3f}")

hr("2. the molecule-less set is read-supported but suspect")
ml = df[~measured]
print(f"FL reads of molecule-less: median {ml['FL_reads'].median():.0f}  "
      f"(all passed the >=20 FL-read filter)")
print("structural_category of molecule-less:")
print((ml["structural_category"].value_counts(normalize=True) * 100).round(1).head(6).to_string())
print(f"intra-primed among molecule-less: {100*ml['intrapriming'].mean():.1f}%  "
      f"vs measured: {100*df.loc[measured,'intrapriming'].mean():.1f}%")

hr("3. STRUCTURE x RECURRENCE grid (measured isoforms; CRAFT consequence joined)")
g = df[measured].assign(
    klass=np.select(
        [fsm[measured], df.loc[measured, "structural_category"].eq("incomplete-splice_match"),
         novel[measured], junk[measured]],
        ["FSM", "ISM", "novel", "junk(interg/anti/genic/fusion)"], default="other",
    )
)
tab = g.groupby("klass").agg(
    n=("transcript_id", "size"),
    median_ncells=("n_cells_detected", "median"),
    median_FLreads=("FL_reads", "median"),
    frac_intraprimed=("intrapriming", "mean"),
    frac_ncells_ge3=("n_cells_detected", lambda s: (s >= 3).mean()),
    frac_craft_full_length=("completeness", lambda s: s.eq("full_length").mean()),
    frac_resolved_intact=("resolved_orf_status", lambda s: s.eq("intact").mean()),
    frac_nmd_sensitive=("nmd_status", lambda s: s.eq("sensitive").mean()),
    frac_coding=("coding_potential_label", lambda s: s.eq("coding").mean()),
).reindex(["FSM", "ISM", "novel", "junk(interg/anti/genic/fusion)"])
print(tab.to_string(float_format=lambda v: f"{v:.2f}"))

hr("4. CRAFT consequence x recurrence (new on v18: are PTC/NMD calls recurrent?)")
for col in ("resolved_orf_status", "nmd_status"):
    print(f"\nmedian n_cells by {col} (measured):")
    print(df[measured].groupby(col)["n_cells_detected"].agg(["size", "median"]).to_string())

hr("5. FL-read filter vs recurrence (confirm the depth-currency argument)")
m = df[measured]
base = m["n_cells_detected"] >= 5
print(f"{'FL_cut':>8}{'kept':>10}{'kept n_cells>=5':>18}{'lost n_cells>=5':>18}")
for c in (20, 50, 100, 200):
    kept = m["FL_reads"] >= c
    print(f"{c:>8}{int(kept.sum()):>10}{int((kept & base).sum()):>18}{int((~kept & base).sum()):>18}")

hr("6. RECOMMENDED class-conditional recurrence filter, applied")
keep = measured & (
    ((fsm | novel) & (df["n_cells_detected"] >= 3))
    | (junk & (df["n_cells_detected"] >= 5) & ~df["intrapriming"] & ~df["noncanonical"])
    | (df["structural_category"].eq("incomplete-splice_match") & (df["n_cells_detected"] >= 3))
)
print(f"kept: {int(keep.sum())} of {int(measured.sum())} measured "
      f"({100*keep.sum()/measured.sum():.1f}%); of all {len(df)} CRAFT isoforms "
      f"({100*keep.sum()/len(df):.1f}%)")
print("kept composition by class:")
print(df[keep].assign(
    klass=np.select([fsm[keep], novel[keep], junk[keep],
                     df.loc[keep, "structural_category"].eq("incomplete-splice_match")],
                    ["FSM", "novel", "junk", "ISM"], default="other")
)["klass"].value_counts().to_string())
print(f"\nfull_length kept: {int((df[keep]['completeness'].eq('full_length')).sum())}; "
      f"resolved intact kept: {int((df[keep]['resolved_orf_status'].eq('intact')).sum())}")

hr("7. SAVE enriched table")
df.to_csv(OUT / "v18_structure_recurrence.tsv", sep="\t", index=False)
print(f"wrote {OUT/'v18_structure_recurrence.tsv'}  ({len(df)} rows)")
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
    ax[0].scatter(df.loc[measured, "FL_reads"].clip(20, 1e5),
                  df.loc[measured, "n_cells_detected"].clip(1, 3000),
                  s=3, alpha=0.05,
                  c=np.where(fsm[measured], "#2e8b57", np.where(junk[measured], "#d1495b", "#9aa0a6")))
    ax[0].set(xscale="log", yscale="log", xlabel="FL reads", ylabel="n_cells_detected",
              title="FL reads vs recurrence (green FSM, red junk)")
    for kl, c in [("FSM", "#2e8b57"), ("novel", "#3b6fb6"), ("junk(interg/anti/genic/fusion)", "#d1495b")]:
        s = g[g["klass"] == kl]["n_cells_detected"].clip(1, 3000)
        ax[1].hist(s, bins=np.logspace(0, np.log10(3000), 40), histtype="step", label=kl, color=c)
    ax[1].set(xscale="log", xlabel="n_cells_detected", ylabel="isoforms",
              title="recurrence by class"); ax[1].legend(fontsize=7)
    order = ["intact", "ptc_premature", "ptc_intron_retained", "cds_extension",
             "no_stop_in_read", "resolution_failed"]
    med = df[measured].groupby("resolved_orf_status")["n_cells_detected"].median().reindex(order).dropna()
    ax[2].barh(range(len(med)), med.values, color="#4c72b0")
    ax[2].set_yticks(range(len(med))); ax[2].set_yticklabels(med.index, fontsize=8)
    ax[2].set(xlabel="median n_cells", title="recurrence by resolved ORF status")
    fig.tight_layout(); fig.savefig(OUT / "v18_structure_recurrence.png", dpi=110)
    print(f"wrote {OUT/'v18_structure_recurrence.png'}")
except Exception as e:
    print(f"(plots skipped: {e})")
print("\nDONE")
