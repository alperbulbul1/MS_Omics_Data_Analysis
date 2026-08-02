"""
process_GSE118257_jakel_brain.py
------------------------
MS brain snRNA-seq — Jäkel et al. 2019 (Nature, PMID 30747918, GSE118257).

The GEO deposit ships the paper's FINAL cluster labels + cell-type names:
  Sample, Condition (Ctrl vs MS), Lesion (Ctrl/A/CA/CI/NAWM/RM),
  Clusters_res08 (paper Seurat clusters, 0-12),
  Celltypes (named: Oligo1..6, OPCs, COPs, ImOlGs, Astrocytes, Astrocytes2,
             Microglia_Macrophages, Macrophages, Immune_cells,
             Endothelial_cells1/2, Pericytes, Vasc_smooth_muscle,
             Neuron1..5)

We use those labels as-is — they are the paper's clusters.

UMAP topology is anchored on the paper Celltypes via PAGA initialisation
(paper's UMAP coords are not deposited, but topology will mirror Fig 1).

Outputs in results/figures/brain_Jakel2019/:
  fig1_umap_celltypes.png       (paper celltypes)
  fig1b_umap_clusters.png       (paper Clusters_res08)
  fig1c_smallmultiples.png
  fig2_umap_lesion.png          (Ctrl / NAWM / A / CA / CI / RM)
  fig3_umap_MS_vs_HC.png
  fig4_candidate_genes.png
  fig6_violin_paperclusters.png
  stats_per_celltype_MSvsHC.csv
"""
from __future__ import annotations
import gzip, time, warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests

import anndata as ad
import scanpy as sc
import bbknn
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=150, frameon=False, fontsize=9)

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "brain_Jakel2019_GSE118257"
FIG  = ROOT / "results" / "figures" / "brain_Jakel2019"
FIG.mkdir(parents=True, exist_ok=True)
H5   = FIG / "adata_jakel.h5ad"

CAND = ["LXN", "SH3BP4", "THRB", "CHL1", "RPAP2", "PCNP"]

# Group cell types by lineage colour family for cleaner palette
CELLTYPE_PALETTE = {
    # Oligodendrocyte lineage --------------------
    "Oligo1": "#1f77b4", "Oligo2": "#3a8ec0", "Oligo3": "#5fa7d0",
    "Oligo4": "#88bee0", "Oligo5": "#aed3ec", "Oligo6": "#cae3f4",
    "OPCs":  "#9467bd",  "COPs":  "#ba9bd6", "ImOlGs": "#7e57c2",
    # Microglia / immune -------------------------
    "Microglia_Macrophages": "#d62728", "Macrophages": "#a31515",
    "Immune_cells":          "#fc9272",
    # Astrocytes ---------------------------------
    "Astrocytes":  "#2ca02c", "Astrocytes2": "#5dbf65",
    # Vascular ----------------------------------
    "Endothelial_cells1": "#7f7f7f", "Endothelial_cells2": "#bcbcbc",
    "Pericytes":          "#404040", "Vasc_smooth_muscle": "#000000",
    # Neurons ----------------------------------
    "Neuron1": "#ff7f0e", "Neuron2": "#fda863", "Neuron3": "#fdc485",
    "Neuron4": "#fdd9a8", "Neuron5": "#fdebd0",
}
GROUP_PALETTE = {"MS":"#d62728","HC":"#377eb8"}
LESION_PALETTE = {
    "Ctrl": "#377eb8", "NAWM": "#74add1",
    "A": "#fc9272", "CA": "#de2d26", "CI": "#fdae6b", "RM": "#fee0b6",
}


def cohens_d(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if len(a) < 2 or len(b) < 2: return np.nan
    s = np.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))
                 /(len(a)+len(b)-2))
    return (a.mean()-b.mean())/s if s else 0.0


def bh_fdr(df, pcol="wilcoxon_p", by="gene"):
    out = df.copy(); out["fdr"] = np.nan
    for k, sub in out.groupby(by):
        _, q, _, _ = multipletests(sub[pcol].fillna(1.0), method="fdr_bh")
        out.loc[sub.index, "fdr"] = q
    return out


def build_anndata():
    if H5.exists():
        try:
            print(f"[cache] {H5}")
            return sc.read_h5ad(H5)
        except Exception as e:
            print(f"  cache unreadable, rebuilding: {e}")
    t0 = time.time()
    expr_path = DATA / "expr.txt.gz"
    anno_path = DATA / "anno.txt.gz"
    print(f"Reading {expr_path}")
    # gene x cell tab matrix (21582 genes x 17799 cells, ~24 MB compressed)
    chunks, genes, cells_ref = [], [], None
    for ch in pd.read_csv(expr_path, sep="\t", chunksize=2500,
                           index_col=0, compression="gzip"):
        if cells_ref is None: cells_ref = ch.columns.tolist()
        genes.extend(ch.index.astype(str).tolist())
        chunks.append(sp.csr_matrix(ch.to_numpy(dtype=np.float32, copy=False)))
    X = sp.vstack(chunks, format="csr").T.tocsr()  # cells x genes
    a = ad.AnnData(X=X, obs=pd.DataFrame(index=cells_ref),
                    var=pd.DataFrame(index=genes))
    a.var_names_make_unique()

    anno = pd.read_csv(anno_path, sep="\t", index_col=0)
    common = a.obs_names.intersection(anno.index)
    a = a[common].copy()
    a.obs = anno.loc[common].copy()
    a.obs.rename(columns={"Clusters_res08":"paper_cluster",
                           "Celltypes":"celltype",
                           "Sample":"sample",
                           "Condition":"condition",
                           "Lesion":"lesion"}, inplace=True)
    a.obs["paper_cluster"] = a.obs["paper_cluster"].astype(str)
    a.obs["celltype"]      = a.obs["celltype"].astype(str)
    a.obs["group"] = np.where(a.obs["condition"].astype(str)=="MS", "MS", "HC")
    print(f"AnnData {a.shape} built in {time.time()-t0:.1f}s")
    print("Group:", a.obs["group"].value_counts().to_dict())

    a.write_h5ad(H5, compression="gzip")
    return a


def main():
    t0 = time.time()
    a = build_anndata()

    # Normalize -> log -> HVG -> PCA -> bbknn -> PAGA-init UMAP
    if "X_umap" not in a.obsm:
        sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
        sc.pp.highly_variable_genes(a, n_top_genes=3000, flavor="cell_ranger")
        a.raw = a
        hv = a[:, a.var["highly_variable"] | a.var_names.isin(CAND)].copy()
        sc.pp.scale(hv, max_value=10)
        sc.tl.pca(hv, n_comps=50)
        bbknn.bbknn(hv, batch_key="sample", n_pcs=40, neighbors_within_batch=6)
        sc.tl.paga(hv, groups="celltype")
        sc.pl.paga(hv, show=False, plot=False)
        sc.tl.umap(hv, init_pos="paga", min_dist=0.35, spread=1.4)
        a.obsm["X_umap"] = hv.obsm["X_umap"]
        a.obsm["X_pca"]  = hv.obsm["X_pca"]
        a.write_h5ad(H5, compression="gzip")

    # Figures ------------------------------------------------------------
    sc.pl.umap(a, color="celltype", palette=CELLTYPE_PALETTE, size=8,
                legend_loc="right margin", frameon=False,
                title="Jäkel 2019 — paper celltypes", show=False)
    plt.savefig(FIG/"fig1_umap_celltypes.png", dpi=300, bbox_inches="tight")
    plt.close()

    sc.pl.umap(a, color="paper_cluster", palette="tab20", size=8,
                legend_loc="on data", frameon=False, legend_fontsize=8,
                title="Jäkel 2019 — paper Clusters_res08", show=False)
    plt.savefig(FIG/"fig1b_umap_clusters.png", dpi=300, bbox_inches="tight")
    plt.close()

    sc.pl.umap(a, color="lesion", palette=LESION_PALETTE, size=8,
                legend_loc="right margin", frameon=False,
                title="Jäkel 2019 — lesion category", show=False)
    plt.savefig(FIG/"fig2_umap_lesion.png", dpi=300, bbox_inches="tight")
    plt.close()

    sc.pl.umap(a, color="group", palette=GROUP_PALETTE, size=8,
                legend_loc="right margin", frameon=False,
                title="Jäkel 2019 — MS vs HC", show=False)
    plt.savefig(FIG/"fig3_umap_MS_vs_HC.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Small multiples per celltype
    cts = sorted(a.obs["celltype"].unique())
    cols = 5; rows = int(np.ceil(len(cts)/cols))
    fig, axes = plt.subplots(rows, cols, figsize=(3.6*cols, 3.2*rows),
                              sharex=True, sharey=True)
    for ax, ct in zip(axes.flat, cts):
        mask = (a.obs["celltype"].astype(str) == ct).values
        ax.scatter(a.obsm["X_umap"][~mask,0], a.obsm["X_umap"][~mask,1],
                    s=1.5, c="#dddddd", linewidths=0, alpha=0.7,
                    rasterized=True)
        ax.scatter(a.obsm["X_umap"][mask,0], a.obsm["X_umap"][mask,1],
                    s=4, c=CELLTYPE_PALETTE.get(ct,"#333"), linewidths=0,
                    alpha=0.95, rasterized=True)
        ax.set_title(f"{ct} ({mask.sum()})", fontsize=8, fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(False)
    for ax in axes.flat[len(cts):]: ax.axis("off")
    fig.suptitle("Jäkel 2019 — paper celltype small multiples",
                  fontsize=13, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(FIG/"fig1c_smallmultiples.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Candidate genes
    present = [g for g in CAND if g in a.raw.var_names]
    fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharex=True, sharey=True)
    for ax, g in zip(axes.flat, present):
        x = a.raw[:, g].X
        x = x.toarray().ravel() if sp.issparse(x) else np.asarray(x).ravel()
        order = np.argsort(x)
        sc_ = ax.scatter(a.obsm["X_umap"][order,0], a.obsm["X_umap"][order,1],
                          s=4, c=x[order], cmap="magma_r", vmin=0,
                          vmax=max(float(np.quantile(x,0.99)), 0.2),
                          linewidths=0, alpha=0.9, rasterized=True)
        ax.set_title(g, fontsize=13, fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([])
        for sp_ in ax.spines.values(): sp_.set_visible(False)
        plt.colorbar(sc_, ax=ax, shrink=0.7, pad=0.02)
    fig.suptitle("Candidate MS genes (Jäkel 2019 brain snRNA-seq)",
                  fontsize=14, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(FIG/"fig4_candidate_genes.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Split violins
    raw = a.raw.to_adata()
    longs = []
    for g in present:
        x = raw[:, g].X
        x = x.toarray().ravel() if sp.issparse(x) else np.asarray(x).ravel()
        longs.append(pd.DataFrame({"gene":g,"expr":x,
                                    "celltype":a.obs["celltype"].astype(str).values,
                                    "group":a.obs["group"].astype(str).values}))
    long = pd.concat(longs, ignore_index=True)
    fig, axes = plt.subplots(2, 3, figsize=(20, 10))
    for ax, g in zip(axes.flat, present):
        sub = long[long["gene"]==g]
        keep = sub.groupby("celltype")["group"].apply(
            lambda s: (s=="MS").sum()>=20 and (s=="HC").sum()>=20)
        sub = sub[sub["celltype"].isin(keep[keep].index)]
        if sub.empty: continue
        order = sorted(sub["celltype"].unique())
        sns.violinplot(data=sub, x="celltype", y="expr", hue="group",
                        palette=GROUP_PALETTE, split=True, inner="quartile",
                        order=order, ax=ax, linewidth=0.4, cut=0)
        ax.set_title(g, fontsize=12, fontweight="bold")
        ax.set_xlabel(""); ax.set_ylabel("log-norm expr")
        for lbl in ax.get_xticklabels():
            lbl.set_rotation(60); lbl.set_ha("right"); lbl.set_fontsize(7)
        if ax is not axes.flat[0] and ax.get_legend() is not None:
            ax.get_legend().remove()
    fig.suptitle("Candidate gene expression — Jäkel 2019, MS vs HC by celltype",
                  fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(FIG/"fig6_violin_paperclusters.png", dpi=300,
                bbox_inches="tight"); plt.close()

    # Stats per celltype
    rows = []
    for g in present:
        x = raw[:, g].X
        x = x.toarray().ravel() if sp.issparse(x) else np.asarray(x).ravel()
        for ct, idx in raw.obs.groupby("celltype").groups.items():
            sel = raw.obs_names.get_indexer(idx)
            grp = raw.obs.loc[idx,"group"].values
            sub = x[sel]; ms=sub[grp=="MS"]; hc=sub[grp=="HC"]
            if len(ms)<10 or len(hc)<10: continue
            stat, p = ranksums(ms, hc)
            rows.append({"gene":g,"celltype":ct,
                          "n_MS":int(len(ms)),"n_HC":int(len(hc)),
                          "mean_MS":float(ms.mean()),"mean_HC":float(hc.mean()),
                          "logfc":float(ms.mean()-hc.mean()),
                          "cohens_d":float(cohens_d(ms,hc)),
                          "pct_MS":float((ms>0).mean()*100),
                          "pct_HC":float((hc>0).mean()*100),
                          "wilcoxon_p":float(p)})
    df = pd.DataFrame(rows)
    df = bh_fdr(df, pcol="wilcoxon_p", by="gene").sort_values(["gene","fdr"])
    df.to_csv(FIG/"stats_per_celltype_MSvsHC.csv", index=False)
    print(f"done {time.time()-t0:.1f}s -> {FIG}")


if __name__ == "__main__":
    main()
