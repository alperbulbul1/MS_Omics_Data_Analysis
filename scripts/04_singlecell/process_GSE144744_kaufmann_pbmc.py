"""
process_GSE144744_kaufmann_pbmc.py
-------------------------
MS PBMC + CSF scRNA-seq — Ramesh et al. 2020 (PNAS, PMID 32907943, GSE144744).

Cohort: 62 donors, 497,706 cells across three sub-cohorts:
    - RRMS_HI   : RRMS patients + matched healthy controls
    - NAT_HI    : MS patients on natalizumab + matched controls
    - PPMS_HI   : PPMS patients + matched controls
Disease groups in column `group`: HI1/HI2/HI3 (controls), MS1, MS2, MS1_nat, PPMS.

The GEO deposit ships rich per-cell metadata with the paper's exact clusters
(`cluster_names`: T01..T11, M01..M07, NK01..NK02, B01, CDC01, PDC01, PT01,
PLC01) and broad types (`basictype`: t_cells / monocytes / nk_cells / b_cells
/ cdc / pdc / plasma_cells / platelets).

The full SaverX-imputed expression tarball is only a 425-gene immune panel
that does NOT include our candidate genes (LXN, SH3BP4, THRB, CHL1, RPAP2,
PCNP), and the full normalized RNA matrix is 1.5–1.8 GB which doesn't fit
in this sandbox's per-call time budget. So:

  IN-SANDBOX run ➜ paper-cluster composition + stratification figures
  LOCAL-LAPTOP run ➜ also ingests RNA_normalised matrix to compute
                     candidate-gene expression per paper cluster

Outputs in results/figures/blood_Ramesh2020/:
  fig1_celltype_composition.png      stacked bar (basictype, MS vs HC, per cohort)
  fig1b_cluster_composition.png      paper cluster_names by group
  fig2_donors_per_group.png
  fig3_cohort_summary.png
  composition_basictype.csv
  composition_cluster.csv
  cohort_donor_breakdown.csv
"""
from __future__ import annotations
import gzip, time, warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "blood_Ramesh2020_GSE144744"
FIG  = ROOT / "results" / "figures" / "blood_Ramesh2020"
FIG.mkdir(parents=True, exist_ok=True)

GROUP_PALETTE = {
    "HI1":"#9ecae1", "HI2":"#6baed6", "HI3":"#3182bd",
    "MS1":"#fb6a4a", "MS2":"#cb181d", "MS1_nat":"#fcae91", "PPMS":"#a50f15",
}
BASIC_PALETTE = {
    "t_cells":"#1f77b4", "monocytes":"#d62728", "nk_cells":"#e377c2",
    "b_cells":"#ff7f0e", "cdc":"#9467bd", "pdc":"#c5b0d5",
    "plasma_cells":"#bcbd22", "platelets":"#7f7f7f",
}

# Coarse paper cluster name -> color family
def cluster_palette(cluster_names):
    base = {
        "T":"#1f77b4", "M":"#d62728", "N":"#e377c2", "B":"#ff7f0e",
        "C":"#9467bd", "P":"#7f7f7f",
    }
    return {c: base.get(c[0], "#444") for c in cluster_names}


def main():
    t0 = time.time()
    print("Reading cell metadata ...")
    meta = pd.read_csv(DATA/"cell_meta.csv.gz", low_memory=False)
    print(f"  {len(meta):,} cells, {meta['donor'].nunique()} donors")
    meta["disease_status"] = np.where(
        meta["group"].astype(str).str.startswith("HI"), "HC", "MS")

    # ---------- composition: basictype × group ------------------------------
    bt = (meta.groupby(["cohort","disease_status","basictype"]).size()
              .reset_index(name="n"))
    bt["pct"] = bt.groupby(["cohort","disease_status"])["n"].transform(
        lambda s: 100*s/s.sum())
    bt.to_csv(FIG/"composition_basictype.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
    cohorts = ["RRMS_HI","NAT_HI","PPMS_HI"]
    for ax, coh in zip(axes, cohorts):
        sub = bt[bt["cohort"]==coh]
        pivot = sub.pivot(index="basictype", columns="disease_status", values="pct").fillna(0)
        pivot = pivot.reindex(["t_cells","nk_cells","b_cells","monocytes",
                                "cdc","pdc","plasma_cells","platelets"])
        pivot.plot(kind="bar", ax=ax,
                    color={"HC":"#377eb8","MS":"#d62728"},
                    edgecolor="black", linewidth=0.4)
        ax.set_title(f"{coh}", fontsize=12, fontweight="bold")
        ax.set_xlabel(""); ax.set_ylabel("% of cells")
        for t in ax.get_xticklabels(): t.set_rotation(35); t.set_ha("right")
        if ax is not axes[0]: ax.set_ylabel("")
    fig.suptitle("Ramesh 2020 — basic cell type composition, MS vs HC",
                  fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(FIG/"fig1_celltype_composition.png", dpi=300, bbox_inches="tight")
    plt.close()

    # ---------- composition: paper cluster_names × group --------------------
    cl = (meta.groupby(["cohort","disease_status","cluster_names"]).size()
              .reset_index(name="n"))
    cl["pct"] = cl.groupby(["cohort","disease_status"])["n"].transform(
        lambda s: 100*s/s.sum())
    cl.to_csv(FIG/"composition_cluster.csv", index=False)

    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    for ax, coh in zip(axes, cohorts):
        sub = cl[cl["cohort"]==coh]
        pivot = sub.pivot(index="cluster_names", columns="disease_status",
                          values="pct").fillna(0)
        # order: T01..T11, then M, NK, B, etc.
        pivot["sort"] = pivot.index.str.replace(r"\d+","",regex=True)
        pivot = pivot.sort_index().sort_values("sort", kind="mergesort")
        pivot.drop(columns="sort").plot(kind="bar", ax=ax,
            color={"HC":"#377eb8","MS":"#d62728"},
            edgecolor="black", linewidth=0.3)
        ax.set_title(f"{coh}", fontsize=11, fontweight="bold")
        ax.set_ylabel("% cells")
        if ax is axes[-1]:
            for t in ax.get_xticklabels(): t.set_rotation(40); t.set_ha("right")
            ax.set_xlabel("paper cluster")
    fig.suptitle("Ramesh 2020 — paper cluster_names composition, MS vs HC",
                  fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(FIG/"fig1b_cluster_composition.png", dpi=300,
                bbox_inches="tight")
    plt.close()

    # ---------- donor breakdown --------------------------------------------
    donors = meta.groupby(["cohort","group"])["donor"].nunique().reset_index(
        name="n_donors")
    cells_per = meta.groupby(["cohort","group"]).size().reset_index(name="n_cells")
    summary = donors.merge(cells_per, on=["cohort","group"])
    summary.to_csv(FIG/"cohort_donor_breakdown.csv", index=False)

    fig, ax = plt.subplots(figsize=(10, 5))
    pivot = summary.pivot(index="group", columns="cohort", values="n_donors").fillna(0)
    pivot.plot(kind="bar", ax=ax, edgecolor="black", linewidth=0.4)
    ax.set_ylabel("Number of donors")
    ax.set_title("Ramesh 2020 — donors per group × cohort",
                  fontsize=12, fontweight="bold")
    for t in ax.get_xticklabels(): t.set_rotation(0); t.set_ha("center")
    plt.tight_layout()
    plt.savefig(FIG/"fig2_donors_per_group.png", dpi=300, bbox_inches="tight")
    plt.close()

    # ---------- overall summary --------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    overall_bt = (meta.groupby(["disease_status","basictype"]).size()
                       .reset_index(name="n"))
    overall_bt["pct"] = overall_bt.groupby("disease_status")["n"].transform(
        lambda s: 100*s/s.sum())
    pivot = overall_bt.pivot(index="basictype", columns="disease_status",
                              values="pct").fillna(0)
    pivot = pivot.reindex(["t_cells","nk_cells","b_cells","monocytes",
                            "cdc","pdc","plasma_cells","platelets"])
    pivot.plot(kind="bar", ax=axes[0],
                color={"HC":"#377eb8","MS":"#d62728"},
                edgecolor="black", linewidth=0.4)
    axes[0].set_title("Pooled — basic types", fontweight="bold")
    axes[0].set_ylabel("% cells")
    for t in axes[0].get_xticklabels(): t.set_rotation(35); t.set_ha("right")

    overall_cl = (meta.groupby(["disease_status","cluster_names"]).size()
                       .reset_index(name="n"))
    overall_cl["pct"] = overall_cl.groupby("disease_status")["n"].transform(
        lambda s: 100*s/s.sum())
    pivot2 = overall_cl.pivot(index="cluster_names", columns="disease_status",
                                values="pct").fillna(0)
    pivot2["delta"] = pivot2["MS"] - pivot2["HC"]
    pivot2 = pivot2.sort_values("delta", ascending=False)
    pivot2["delta"].plot(kind="bar", ax=axes[1],
                          color=np.where(pivot2["delta"]>0, "#d62728","#377eb8"),
                          edgecolor="black", linewidth=0.3)
    axes[1].set_title("Δ paper cluster (MS − HC)", fontweight="bold")
    axes[1].set_ylabel("% point difference")
    axes[1].axhline(0, color="black", lw=0.5)
    for t in axes[1].get_xticklabels(): t.set_rotation(60); t.set_ha("right"); t.set_fontsize(7)
    plt.tight_layout()
    plt.savefig(FIG/"fig3_cohort_summary.png", dpi=300, bbox_inches="tight")
    plt.close()

    print(f"done {time.time()-t0:.1f}s -> {FIG}")
    print("\nNote: candidate gene expression analysis requires the full RNA "
          "matrix (GSE144744_RNA_normalised.tar.gz, 1.8 GB). Run on your "
          "laptop with scripts/process_GSE144744_kaufmann_genes.py once the file "
          "is downloaded into data/blood_Ramesh2020_GSE144744/.")


if __name__ == "__main__":
    main()
