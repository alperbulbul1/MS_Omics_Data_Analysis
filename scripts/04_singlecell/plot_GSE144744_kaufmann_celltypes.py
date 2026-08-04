"""
plot_GSE144744_kaufmann_celltypes.py
--------------------
Re-render Ramesh 2020 figures with a clearer palette using the cached
adata_ramesh_umap.h5ad. Avoids re-running the slow UMAP step.
"""
import warnings
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=150, frameon=False, fontsize=9)

# In the authors' tree this file sat in SingleCell_CELLxGENE/scripts/, so parent.parent was the
# single-cell data root. In the release it resolves to <repo>/scripts/, which would stream the
# ~18 GB single-cell tree into the git checkout and put it where no downstream script looks.
# Use the same placeholder the sibling pseudobulk scripts use.
ROOT = Path("__MS_GEO_ROOT__") / "SingleCell_CELLxGENE"
FIG  = ROOT / "results" / "figures" / "blood_Ramesh2020_UMAP"
H5   = FIG / "adata_ramesh_umap.h5ad"

GROUP_PALETTE = {"MS":"#d62728","HC":"#377eb8"}
COHORT_PALETTE = {"PPMS_HI":"#a50f15","NAT_HI":"#fd8d3c","RRMS_HI":"#fb6a4a"}
BASIC_PALETTE = {
    "t_cells":"#1f77b4","monocytes":"#d62728","nk_cells":"#e377c2",
    "b_cells":"#ff7f0e","cdc":"#9467bd","pdc":"#c5b0d5",
    "plasma_cells":"#bcbd22","platelets":"#7f7f7f",
}
# 25 clearly distinct colors for paper clusters
DISTINCT25 = [
    "#1f77b4","#aec7e8","#5b8cb1","#94c0db",  # T cells (blue family)
    "#2ca02c","#98df8a","#5fb360",
    "#d62728","#ff9896","#a31515","#e7211d",  # monocytes (red)
    "#7f0000","#ff5757","#fb9a99",
    "#e377c2","#f7b6d2","#cd5790",            # NK (pink)
    "#ff7f0e","#ffbb78",                       # B (orange)
    "#9467bd","#c5b0d5",                       # DC (purple)
    "#bcbd22","#dbdb8d",                       # Plasma
    "#7f7f7f","#000000",
]


def main():
    a = sc.read_h5ad(H5)
    cluster_names = sorted(a.obs["cluster_names"].astype(str).unique())
    PAL = {n: DISTINCT25[i % len(DISTINCT25)] for i, n in enumerate(cluster_names)}

    sc.pl.umap(a, color="cluster_names", palette=PAL, size=8,
                legend_loc="right margin", frameon=False, legend_fontsize=8,
                title="Ramesh 2020 — paper cluster_names",
                show=False)
    plt.savefig(FIG/"fig1_umap_paperclusters.png", dpi=300,
                bbox_inches="tight"); plt.close()

    sc.pl.umap(a, color="basictype", palette=BASIC_PALETTE, size=8,
                legend_loc="right margin", frameon=False,
                title="Ramesh 2020 — basic cell type", show=False)
    plt.savefig(FIG/"fig1b_umap_basictype.png", dpi=300,
                bbox_inches="tight"); plt.close()

    sc.pl.umap(a, color="cohort", palette=COHORT_PALETTE, size=8,
                legend_loc="right margin", frameon=False,
                title="Ramesh 2020 — cohort", show=False)
    plt.savefig(FIG/"fig2_umap_cohort.png", dpi=300, bbox_inches="tight")
    plt.close()

    sc.pl.umap(a, color="disease_status", palette=GROUP_PALETTE, size=8,
                legend_loc="right margin", frameon=False,
                title="Ramesh 2020 — MS vs HC", show=False)
    plt.savefig(FIG/"fig3_umap_MS_vs_HC.png", dpi=300, bbox_inches="tight")
    plt.close()

    sc.pl.umap(a, color="donor", size=6, legend_loc=None,
                frameon=False, title="Ramesh 2020 — donor (62 donors)",
                show=False)
    plt.savefig(FIG/"fig4_umap_donor.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Side-by-side MS vs HC by basictype (immune lineage colors)
    fig, axes = plt.subplots(1, 2, figsize=(14, 7), sharex=True, sharey=True)
    for ax, ds in zip(axes, ["MS","HC"]):
        sub = a[a.obs["disease_status"] == ds]
        for ct in sorted(sub.obs["basictype"].unique()):
            mask = (sub.obs["basictype"] == ct).values
            if not mask.any(): continue
            ax.scatter(sub.obsm["X_umap"][mask,0], sub.obsm["X_umap"][mask,1],
                        s=4, c=[BASIC_PALETTE.get(ct,"#000")],
                        label=ct, linewidths=0, alpha=0.85, rasterized=True)
        ax.set_title(f"{ds} (n={sub.n_obs:,})", fontsize=12, fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(False)
    axes[0].legend(loc="upper left", bbox_to_anchor=(-0.05,-0.05),
                    ncol=4, frameon=False, fontsize=9)
    fig.suptitle("Ramesh 2020 PBMC — MS vs HC by basic cell type",
                  fontsize=13, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(FIG/"fig3b_umap_split_MS_HC.png", dpi=300,
                bbox_inches="tight"); plt.close()

    # Composition delta plot
    comp = (a.obs.groupby(["disease_status","cluster_names"]).size()
                 .reset_index(name="n"))
    comp["pct"] = comp.groupby("disease_status")["n"].transform(
        lambda s: 100*s/s.sum())
    pivot = comp.pivot(index="cluster_names", columns="disease_status",
                        values="pct").fillna(0)
    pivot["delta"] = pivot["MS"] - pivot["HC"]
    pivot = pivot.sort_values("delta", ascending=False)
    fig, ax = plt.subplots(figsize=(11, 4.5))
    pivot["delta"].plot(kind="bar", ax=ax,
                        color=np.where(pivot["delta"]>0, "#d62728", "#377eb8"),
                        edgecolor="black", linewidth=0.4)
    ax.axhline(0, color="black", lw=0.5)
    ax.set_ylabel("MS − HC (pp)")
    ax.set_xlabel("paper cluster_names")
    ax.set_title("Δ paper cluster (MS − HC) — Ramesh 2020 PBMC subsample",
                  fontsize=12, fontweight="bold")
    for t in ax.get_xticklabels(): t.set_rotation(60); t.set_ha("right"); t.set_fontsize(8)
    plt.tight_layout()
    plt.savefig(FIG/"fig5_composition_delta.png", dpi=300, bbox_inches="tight")
    plt.close()

    print("Redraw done.")


if __name__ == "__main__":
    main()
