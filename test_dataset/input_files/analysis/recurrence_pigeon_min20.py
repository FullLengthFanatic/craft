#!/usr/bin/env python3
"""Depth-stable filtering analysis on the BD70 pigeon_min20 lineage.

Currencies for keeping an isoform:
  - FL reads (SQANTI FL.sample)  : what the >=20 filter used; PCR-amplified, depth-dependent
  - molecules (UMI-corrected)    : ~= n_cells here (1 molecule/cell), so collapses onto recurrence
  - recurrence (n_cells_detected, over CALLED cells)  <- depth-stable
joined to SQANTI structural class + artifact signals and CRAFT functional consequence.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import mmread

BASE = Path("/data/long_read_annotation_tool/test_dataset/input_files")
PG = BASE / "pigeon_isoforms_min20counts"
MTX = PG / "BD70_min20counts_isoforms_trimmed_umi_corrected_isoform_matrix.mtx"
ROWS = PG / "BD70_min20counts_isoforms_trimmed_isoforms.txt"
COLS = PG / "BD70_min20counts_isoforms_trimmed_cells.txt"
META = PG / "BD70_min20counts_isoforms_trimmed_isoform_metadata.txt"
KNEE = BASE / "BD70_strict_barcode_knee_knee_data.tsv"
SQANTI = BASE / "classification/BD70.trimmed.min20count_classification.txt"
CRAFT = Path("/data/long_read_annotation_tool/test_dataset/BD70_min20count/craft_out_v17/per_isoform.tsv")
OUT = BASE / "analysis"; OUT.mkdir(exist_ok=True)
N_CELLS = 3000          # primary called-cell count (top of the rank/read cliff)

def hr(t): print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)

# ---------------------------------------------------------------- knee / cells
hr("1. CELL CALLING")
knee = pd.read_csv(KNEE, sep="\t")
knee = knee[knee["read_count"] > 0].sort_values("rank").reset_index(drop=True)
# transparency: distance-to-chord knee on the curve truncated above the singleton tail
sig = knee[knee["read_count"] >= 3]
lr, lc = np.log10(sig["rank"].values), np.log10(sig["read_count"].values)
xn = (lr - lr.min()) / (lr.max() - lr.min()); yn = (lc - lc.min()) / (lc.max() - lc.min())
chord = yn[0] + (yn[-1] - yn[0]) * (xn - xn[0]) / (xn[-1] - xn[0])
auto = int(sig.iloc[int(np.argmax(chord - yn))]["rank"])
print(f"auto knee (truncated chord): rank {auto}, read_count "
      f"{int(knee.loc[knee['rank']==auto,'read_count'].iloc[0])}")
print(f"read_count at candidate cutoffs:")
for r in (2000, 2500, 3000, 3500, 4000):
    rc = int(knee.loc[knee["rank"] == r, "read_count"].iloc[0])
    cum = 100 * knee.loc[knee["rank"] <= r, "read_count"].sum() / knee["read_count"].sum()
    print(f"  rank<= {r}: read_count {rc:>6}  cum_reads {cum:4.1f}%")
print(f"--> using top {N_CELLS} barcodes as called cells (sensitivity at 2500 shown later)")

def whitelist(n): return set(knee.loc[knee["rank"] <= n, "cell_barcode"])

# ---------------------------------------------------------------- matrix
hr("2. MATRIX -> recurrence over called cells")
iso_ids = pd.read_csv(ROWS, header=None)[0].to_numpy()
barcodes = pd.read_csv(COLS, header=None)[0].to_numpy()
m = mmread(str(MTX)).tocsc()
print(f"matrix {m.shape}, nnz={m.nnz:,}")

def per_iso(n):
    mask = np.fromiter((b in whitelist(n) for b in barcodes), bool, len(barcodes))
    sub = m[:, mask]
    return (np.asarray(sub.sum(1)).ravel().astype(int),
            np.asarray((sub > 0).sum(1)).ravel().astype(int), int(mask.sum()))

tot, nc, ncol = per_iso(N_CELLS)
print(f"called barcodes present in matrix: {ncol} (top {N_CELLS})")
df = pd.DataFrame({"transcript_id": iso_ids, "total_count": tot, "n_cells_detected": nc})
mpc = (df["total_count"] / df["n_cells_detected"].replace(0, np.nan))
print(f"molecules per detected cell: median {mpc.median():.2f}  p95 {mpc.quantile(.95):.2f}  "
      f"-> molecules ~= n_cells (use n_cells as the currency)")
# sensitivity: n_cells at 2500 vs 3000
_, nc25, _ = per_iso(2500)
df["n_cells_2500"] = nc25

# ---------------------------------------------------------------- class + gene fraction
hr("3. STRUCTURAL CLASS + within-gene fraction")
meta = pd.read_csv(META, sep="\t").rename(columns={"isoform_id": "transcript_id"})
df = df.merge(meta[["transcript_id", "associated_gene", "structural_category", "exons"]],
              on="transcript_id", how="left")
g = df.groupby("associated_gene")["total_count"].transform("sum")
df["isoform_fraction_within_gene"] = np.where(g > 0, df["total_count"] / g, np.nan)
print(df["structural_category"].value_counts().to_string())

# ---------------------------------------------------------------- SQANTI artifact + FL reads
hr("4. SQANTI artifact axis + FL reads")
sq = pd.read_csv(SQANTI, sep="\t", low_memory=False).rename(columns={"isoform": "transcript_id"})
sq["perc_A"] = pd.to_numeric(sq["perc_A_downstream_TTS"], errors="coerce")
sq["FL_reads"] = pd.to_numeric(sq["FL.sample"], errors="coerce")
df = df.merge(sq[["transcript_id", "perc_A", "RTS_stage", "all_canonical", "FL_reads"]],
              on="transcript_id", how="left")
df["intrapriming"] = df["perc_A"] >= 60
df["rt_switch"] = df["RTS_stage"].astype(str).str.upper().eq("TRUE")
df["noncanonical"] = df["all_canonical"].astype(str).str.lower().eq("noncanonical")
print(f"intra-priming (perc_A>=60): {df['intrapriming'].sum()}   "
      f"RT-switch: {df['rt_switch'].sum()}   noncanonical-jxn: {df['noncanonical'].sum()}")
df["pcr_factor"] = df["FL_reads"] / df["total_count"].replace(0, np.nan)
print(f"FL_reads: median {df['FL_reads'].median():.0f}  max {df['FL_reads'].max():.0f}; "
      f"PCR factor (FL/molec): median {df['pcr_factor'].median():.2f}  p99 {df['pcr_factor'].quantile(.99):.1f}")

# ---------------------------------------------------------------- CRAFT
hr("5. CRAFT functional consequence")
cr = pd.read_csv(CRAFT, sep="\t", low_memory=False,
                 usecols=lambda c: c in {"transcript_id", "completeness", "resolved_orf_status",
                                         "orf_confidence", "coding_potential_label", "nmd_status"})
df = df.merge(cr, on="transcript_id", how="left")
print(f"CRAFT joined for {df['completeness'].notna().sum()} / {len(df)}")

# ---------------------------------------------------------------- THE ARGUMENT
hr("6. FL-READ FILTER vs RECURRENCE")
fsm = df["structural_category"].eq("full-splice_match")
novel = df["structural_category"].isin(["novel_not_in_catalog", "novel_in_catalog"])
junk = df["structural_category"].isin(["intergenic", "antisense", "genic", "fusion",
                                       "genic_intron"])

print("\n(a) FL reads and recurrence rank isoforms differently (PCR/depth):")
hi_fl_low_cell = (df["FL_reads"] >= 100) & (df["n_cells_detected"] <= 2)
print(f"  FL_reads>=100 but n_cells<=2 (FL-filter keeps, recurrence rejects): {hi_fl_low_cell.sum()}")
print(f"     of those: {junk[hi_fl_low_cell].sum()} junk-class, "
      f"{df['intrapriming'][hi_fl_low_cell].sum()} intra-primed, {fsm[hi_fl_low_cell].sum()} FSM")
lo_fl_hi_cell = (df["FL_reads"] < 50) & (df["n_cells_detected"] >= 10)
print(f"  FL_reads<50 but n_cells>=10 (recurrence keeps, a tighter FL cut would drop): "
      f"{lo_fl_hi_cell.sum()}  ({fsm[lo_fl_hi_cell].sum()} FSM)")

print("\n(b) tightening the FL-read filter sheds recurrent isoforms:")
print(f"{'FL_cut':>8}{'kept':>10}{'kept n_cells>=5':>18}{'lost n_cells>=5':>18}")
base = df["n_cells_detected"] >= 5
for c in (20, 50, 100, 200):
    kept = df["FL_reads"] >= c
    print(f"{c:>8}{kept.sum():>10}{(kept & base).sum():>18}{((~kept) & base).sum():>18}")

print("\n(c) recurrence cut, by class (depth-stable, class-aware):")
print(f"{'ncells_cut':>10}{'FSM':>8}{'novel':>8}{'junk':>8}")
for k in (1, 2, 3, 5, 10):
    kept = df["n_cells_detected"] >= k
    print(f"{k:>10}{(kept & fsm).sum():>8}{(kept & novel).sum():>8}{(kept & junk).sum():>8}")

# ---------------------------------------------------------------- grid
hr("7. STRUCTURE x RECURRENCE grid")
df["klass"] = np.select([fsm, df["structural_category"].eq("incomplete-splice_match"), novel, junk],
                        ["FSM", "ISM", "novel", "junk(interg/anti/genic/fusion)"], default="other")
tab = df.groupby("klass").agg(
    n=("transcript_id", "size"),
    median_ncells=("n_cells_detected", "median"),
    median_FLreads=("FL_reads", "median"),
    frac_intraprimed=("intrapriming", "mean"),
    frac_ncells_ge3=("n_cells_detected", lambda s: (s >= 3).mean()),
    frac_craft_fl_intact=("transcript_id", lambda idx: (
        df.loc[idx.index, "completeness"].eq("full_length") &
        df.loc[idx.index, "resolved_orf_status"].eq("intact")).mean()),
).reindex(["FSM", "ISM", "novel", "junk(interg/anti/genic/fusion)"])
print(tab.to_string(float_format=lambda v: f"{v:.2f}"))

# sensitivity of n_cells>=3 retention to the cell cutoff
hr("8. SENSITIVITY to cell cutoff (3000 vs 2500) + SAVE")
for col, lab in [("n_cells_detected", "3000"), ("n_cells_2500", "2500")]:
    print(f"  cells={lab}: FSM n_cells>=3 = {(df[col].ge(3) & fsm).sum()}, "
          f"novel = {(df[col].ge(3) & novel).sum()}, junk = {(df[col].ge(3) & junk).sum()}")
keep_cols = ["transcript_id", "structural_category", "associated_gene", "FL_reads",
             "total_count", "n_cells_detected", "isoform_fraction_within_gene",
             "intrapriming", "rt_switch", "noncanonical", "completeness",
             "resolved_orf_status", "orf_confidence", "coding_potential_label", "nmd_status"]
df[keep_cols].to_csv(OUT / "pigeon_min20_recurrence.tsv", sep="\t", index=False)
print(f"wrote {OUT/'pigeon_min20_recurrence.tsv'}")
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
    ax[0].loglog(knee["rank"], knee["read_count"], lw=1)
    ax[0].axvline(N_CELLS, color="r", ls="--", label=f"cutoff={N_CELLS}")
    ax[0].set(xlabel="barcode rank", ylabel="reads", title="cell-calling knee"); ax[0].legend()
    ax[1].scatter(df["FL_reads"].clip(20, 1e5), df["n_cells_detected"].clip(1, N_CELLS),
                  s=3, alpha=0.06, c=np.where(fsm, "#2e8b57", np.where(junk, "#d1495b", "#9aa0a6")))
    ax[1].set(xscale="log", yscale="log", xlabel="FL reads", ylabel="n_cells_detected",
              title="FL reads vs recurrence (green=FSM, red=junk)")
    for kl, c in [("FSM", "#2e8b57"), ("novel", "#3b6fb6"), ("junk(interg/anti/genic/fusion)", "#d1495b")]:
        s = df[df["klass"] == kl]["n_cells_detected"].clip(1, N_CELLS)
        ax[2].hist(s, bins=np.logspace(0, np.log10(N_CELLS), 40), histtype="step", label=kl, color=c)
    ax[2].set(xscale="log", xlabel="n_cells_detected", ylabel="isoforms",
              title="recurrence by class"); ax[2].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(OUT / "pigeon_min20_recurrence.png", dpi=110)
    print(f"wrote {OUT/'pigeon_min20_recurrence.png'}")
except Exception as e:
    print(f"(plots skipped: {e})")
print("\nDONE")
