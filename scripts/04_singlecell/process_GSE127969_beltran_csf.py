"""
process_GSE127969_beltran_csf.py
--------------------------
MS CSF + PBMC scRNA-seq — Beltrán et al. 2019 (Brain).

Single-cell RNA-seq in monozygotic twins discordant for multiple sclerosis
plus auto-immune encephalitis (Anti-LGI1, Anti-NMDA) controls.

GEO: GSE127969 ships:
   - GSE127969_counts_TPM_ALL.csv.gz     gene × cell TPM matrix (TAB-sep, 58k × 3730)
   - GSE127969_sc_cell_info.txt.gz       paper per-cell metadata
       columns: Cell, Twin, Case, Sample, index.sort, Clones
       Case   ∈ {HD, MS, HD_PBMCs, MS_PBMCs, Enc, SCNI}
       Sample ∈ {CSF, CD4_PBMCs, CD8_PBMCs, CD19CD27_PBMCs, Others_PBMCs}
       index.sort = paper's FACS-sorted cell type
                  ∈ {CD4, CD8, CD4Tcell, CD8Tcell, CD3, Bcell, Plasmablast,
                     NK, Other, CD19CD27, Rest}

We use Case for MS/HC, index.sort for cell-type, Twin for batch correction.
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
DATA = ROOT / "data" / "blood_Beltran2019_GSE127969"
FIG  = ROOT / "results" / "figures" / "blood_Beltran2019"
FIG.mkdir(parents=True, exist_ok=True)
H5   = FIG / "adata_beltran.h5ad"

CAND = ["LXN", "SH3BP4", "THRB", "CHL1", "RPAP2", "PCNP"]

CELLTYPE_PALETTE = {
    "CD4":"#1f77b4", "CD4Tcell":"#aec7e8",
    "CD8":"#2ca02c", "CD8Tcell":"#98df8a",
    "CD3":"#17becf", "NK":"#e377c2",
    "Bcell":"#ff7f0e", "Plasmablast":"#bcbd22",
    "CD19CD27":"#fdbf6b", "Other":"#7f7f7f", "Rest":"#cccccc",
}
GROUP_PALETTE = {"MS":"#d62728","HC":"#377eb8"}


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
            return sc.read_h5ad(H5)
        except Exception as e:
            print(f"cache unreadable: {e} — rebuilding")
    t0 = time.time()
    print("Reading TPM matrix (tab-separated despite .csv extension)...")
    chunks, genes, cells_ref = [], [], None
    for ch in pd.read_csv(DATA/"GSE127969_counts_TPM_ALL.csv.gz",
                           sep="\t", index_col=0, chunksize=4000,
                           compression="gzip", encoding="latin-1"):
        if cells_ref is None: cells_ref = ch.columns.tolist()
        genes.extend(ch.index.astype(str).tolist())
        chunks.append(sp.csr_matrix(ch.to_numpy(dtype=np.float32, copy=False)))
    X = sp.vstack(chunks, format="csr").T.tocsr()
    a = ad.AnnData(X=X, obs=pd.DataFrame(index=cells_ref),
                    var=pd.DataFrame(index=genes))
    a.var_names_make_unique()
    print(f"matrix {a.shape} in {time.time()-t0:.1f}s")

    # cell info
    info = pd.read_csv(DATA/"GSE127969_sc_cell_info.txt.gz", sep="\t",
                       index_col=0)
    common = a.obs_names.intersection(info.index)
    a = a[common].copy()
    a.obs = info.loc[common].copy()
    # MS / HC label: only HD = healthy, MS = patient (drop encephalitis controls
    # to keep contrast clean)
    a.obs["group"] = a.obs["Case"].astype(str).map(
        {"HD":"HC", "MS":"MS",
         "HD_PBMCs":"HC", "MS_PBMCs":"MS"})
    a.obs.rename(columns={"index.sort":"celltype",
                            "Sample":"compartment",
                            "Twin":"twin"}, inplace=True)
    a.obs["celltype"] = a.obs["celltype"].astype(str).str.strip()
    print("Compartment counts:", a.obs["compartment"].value_counts().to_dict())
    print("Group counts:", a.obs["group"].value_counts(dropna=False).to_dict())
    print("Celltype counts:", a.obs["celltype"].value_counts().to_dict())

    a.write_h5ad(H5, compression="gzip")
    return a


def main():
    t0 = time.time()
    a = build_anndata()
    a.obs["celltype"] = a.obs["celltype"].astype(str).str.strip()

    # Drop encephalitis / SCNI controls to keep MS vs HC contrast clean
    keep = a.obs["group"].notna()
    a = a[keep].copy()
    print(f"After group filter: {a.shape}, MS={int((a.obs['group']=='MS').sum())} "
          f"HC={int((a.obs['group']=='HC').sum())}")

    # TPM is already library-size normalized
    if "X_umap" not in a.obsm:
        sc.pp.log1p(a)
        sc.pp.filter_genes(a, min_cells=5)
        sc.pp.highly_variable_genes(a, n_top_genes=3000, flavor="seurat")
        a.raw = a
        hv = a[:, a.var["highly_variable"] | a.var_names.isin(CAND)].copy()
        sc.pp.scale(hv, max_value=10)
        sc.tl.pca(hv, n_comps=40)
        bbknn.bbknn(hv, batch_key="twin", n_pcs=30, neighbors_within_batch=4)
        sc.tl.paga(hv, groups="celltype")
        sc.pl.paga(hv, show=False, plot=False)
        sc.tl.umap(hv, init_pos="paga", min_dist=0.4, spread=1.6)
        a.obsm["X_umap"] = hv.obsm["X_umap"]
        a.obsm["X_pca"]  = hv.obsm["X_pca"]
        a.write_h5ad(H5, compression="gzip")

    # ---------------- Figures ----------------
    sc.pl.umap(a, color="celltype", palette=CELLTYPE_PALETTE, size=18,
                legend_loc="right margin", frameon=False,
                title="Beltrán 2019 — paper celltype (FACS index.sort)",
                show=False)
    plt.savefig(FIG/"fig1_umap_celltype.png", dpi=300, bbox_inches="tight")
    plt.close()

    sc.pl.umap(a, color="compartment", size=18,
                legend_loc="right margin", frameon=False,
                title="Beltrán 2019 — compartment (CSF / PBMC sub-sort)",
                show=False)
    plt.savefig(FIG/"fig2_umap_compartment.png", dpi=300, bbox_inches="tight")
    plt.close()

    sc.pl.umap(a, color="group", palette=GROUP_PALETTE, size=18,
                legend_loc="right margin", frameon=False,
                title="Beltrán 2019 — MS vs HC", show=False)
    plt.savefig(FIG/"fig3_umap_MS_vs_HC.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Small multiples
    cts = sorted(a.obs["celltype"].unique())
    cols = 4; rows = int(np.ceil(len(cts)/cols))
    fig, axes = plt.subplots(rows, cols, figsize=(4.2*cols, 3.6*rows),
                              sharex=True, sharey=True)
    for ax, ct in zip(axes.flat, cts):
        mask = (a.obs["celltype"].astype(str) == ct).values
        ax.scatter(a.obsm["X_umap"][~mask,0], a.obsm["X_umap"][~mask,1],
                    s=3, c="#dddddd", linewidths=0, alpha=0.7,
                    rasterized=True)
        ax.scatter(a.obsm["X_umap"][mask,0], a.obsm["X_umap"][mask,1],
                    s=10, c=CELLTYPE_PALETTE.get(ct,"#333"), linewidths=0,
                    alpha=0.95, rasterized=True)
        ax.set_title(f"{ct} ({mask.sum()})", fontsize=10, fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(False)
    for ax in axes.flat[len(cts):]: ax.axis("off")
    fig.suptitle("Beltrán 2019 — paper celltype small multiples",
                  fontsize=13, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(FIG/"fig1c_smallmultiples.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Candidate gene UMAP
    present = [g for g in CAND if g in a.raw.var_names]
    fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharex=True, sharey=True)
    for ax, g in zip(axes.flat, present):
        x = a.raw[:, g].X
        x = x.toarray().ravel() if sp.issparse(x) else np.asarray(x).ravel()
        order = np.argsort(x)
        sc_ = ax.scatter(a.obsm["X_umap"][order,0], a.obsm["X_umap"][order,1],
                          s=8, c=x[order], cmap="magma_r", vmin=0,
                          vmax=max(float(np.quantile(x,0.99)), 0.2),
                          linewidths=0, alpha=0.95, rasterized=True)
        ax.set_title(g, fontsize=13, fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([])
        for sp_ in ax.spines.values(): sp_.set_visible(False)
        plt.colorbar(sc_, ax=ax, shrink=0.7, pad=0.02)
    fig.suptitle("Candidate MS genes — Beltrán 2019 CSF+PBMC scRNA",
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
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    for ax, g in zip(axes.flat, present):
        sub = long[long["gene"]==g]
        keep = sub.groupby("celltype")["group"].apply(
            lambda s: (s=="MS").sum()>=10 and (s=="HC").sum()>=10)
        sub = sub[sub["celltype"].isin(keep[keep].index)]
        if sub.empty: continue
        order = sorted(sub["celltype"].unique())
        sns.violinplot(data=sub, x="celltype", y="expr", hue="group",
                        palette=GROUP_PALETTE, split=True, inner="quartile",
                        order=order, ax=ax, linewidth=0.4, cut=0)
        ax.set_title(g, fontsize=12, fontweight="bold")
        ax.set_xlabel(""); ax.set_ylabel("log1p(TPM)")
        for lbl in ax.get_xticklabels():
            lbl.set_rotation(40); lbl.set_ha("right"); lbl.set_fontsize(8)
        if ax is not axes.flat[0] and ax.get_legend() is not None:
            ax.get_legend().remove()
    fig.suptitle("Candidate gene expression — Beltrán 2019, MS vs HC by celltype",
                  fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(FIG/"fig6_violin_paperclusters.png", dpi=300,
                bbox_inches="tight"); plt.close()

    # Stats
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
