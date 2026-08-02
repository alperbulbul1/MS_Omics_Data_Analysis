#!/usr/bin/env python3
# 11_inverse_scRNA_validation  —  generated from notebook spec


# ============================================================
# # 11 — Inverse-concordance validation in scRNA-seq (Python / scanpy)
# 
# Take the inverse-concordant gene panel from notebook 09 and re-test it
# on every locally-available CELLxGENE AnnData dataset:
# 
# | Dataset       | Compartment | h5ad path                                                 |
# |---------------|-------------|------------------------------------------------------------|
# | Schafflick 2020 | PBMC       | `blood/adata_blood.h5ad`                                  |
# | Schafflick paper UMAP | PBMC | `blood_paperUMAP/adata_paperclusters.h5ad`              |
# | Absinta 2021  | brain WM   | `brain_paperUMAP/adata_paperclusters.h5ad`                |
# | Jakel 2019    | brain      | `brain_Jakel2019/adata_jakel.h5ad`                        |
# | Beltran 2019  | blood/CSF  | `blood_Beltran2019/adata_beltran.h5ad`                    |
# | Ramesh 2020   | blood      | `blood_Ramesh2020_UMAP/adata_ramesh_small.h5ad`           |
# | VitD MS       | blood      | `blood_VitD_MS_GSE239626/adata_vitd_post.h5ad`            |
# 
# For each (gene × dataset × cell_type) we run a Wilcoxon MS-vs-HC test
# with BH FDR + Cohen's d and write the long table to
# `results/INV_scRNA_validation_long.tsv`, plus a per-gene aggregate.
# 
# Uses the **py_scrna** kernel (scanpy 1.11.5, anndata 0.12.11).
# 
# **Outputs**
# - `results/INV_scRNA_validation_long.tsv`
# - `results/INV_scRNA_validation_by_gene.tsv`
# - `figures/11_inverse_scRNA_heatmap.png`
# ============================================================

import os, glob
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu

PROJ = Path("__MS_GEO_ROOT__")
ME   = PROJ / "Methylation"  / "results"
SC_FIG = PROJ / "SingleCell_CELLxGENE" / "results" / "figures"
OUT_DIR = ME
FIG_DIR = PROJ / "Methylation" / "figures"

inv = pd.read_csv(ME / "INVERSE_CONCORDANT_by_gene.tsv", sep="\t")
inv = inv.sort_values(["n_pairings","best_rna_fdr"], ascending=[False, True]).head(40)
GENE_PANEL = inv.gene.tolist()
print(f"Loaded {len(GENE_PANEL)} inverse-concordant genes")
print(inv[['gene','n_pairings','direction']].head(15).to_string(index=False))


DATASETS = [
    ("Schafflick PBMC",      SC_FIG / "blood/adata_blood.h5ad",                     "cell_type"),
    ("Schafflick paperUMAP", SC_FIG / "blood_paperUMAP/adata_paperclusters.h5ad",   "leiden"),
    ("Absinta brain",        SC_FIG / "brain_paperUMAP/adata_paperclusters.h5ad",   "leiden"),
    ("Jakel brain",          SC_FIG / "brain_Jakel2019/adata_jakel.h5ad",           "celltype"),
    ("Beltran blood",        SC_FIG / "blood_Beltran2019/adata_beltran.h5ad",       "cell_type"),
    ("Ramesh blood",         SC_FIG / "blood_Ramesh2020_UMAP/adata_ramesh_small.h5ad","cell_type"),
    ("VitD blood",           SC_FIG / "blood_VitD_MS_GSE239626/adata_vitd_post.h5ad","cell_type"),
]

def find_col(obs, candidates):
    for c in candidates:
        if c in obs.columns: return c
    return None

def cohens_d(a, b):
    na, nb = len(a), len(b)
    if na < 2 or nb < 2: return np.nan
    s = np.sqrt(((na-1)*np.var(a, ddof=1) + (nb-1)*np.var(b, ddof=1))/(na+nb-2))
    return (np.mean(a) - np.mean(b)) / s if s > 0 else 0.0

def per_celltype_de(ad, group_col, ct_col, genes, group_a="MS", group_b="HC"):
    obs = ad.obs.copy()
    if group_col not in obs.columns: return None
    # accept variants
    obs[group_col] = obs[group_col].astype(str)
    rows = []
    cts = obs[ct_col].astype(str).unique() if ct_col in obs.columns else ["all"]
    if ct_col not in obs.columns:
        obs["__all__"] = "all"; ct_col = "__all__"
    # restrict to MS / HC
    cond = obs[group_col].str.upper()
    a_mask = cond.isin(["MS","CASE","PATIENT","RRMS","SPMS","PPMS"])
    b_mask = cond.isin(["HC","CTRL","CONTROL","HEALTHY","IIH"])
    obs_sub = obs[a_mask | b_mask].copy()
    obs_sub["__grp__"] = np.where(a_mask[a_mask | b_mask], "MS", "HC")
    if obs_sub.empty: return None
    ad_sub = ad[obs_sub.index].copy()
    ad_sub.obs["__grp__"] = obs_sub["__grp__"]
    # subset to panel
    keep = [g for g in genes if g in ad_sub.var_names]
    if len(keep) == 0: return pd.DataFrame()
    ad_sub = ad_sub[:, keep]
    X = ad_sub.X.toarray() if hasattr(ad_sub.X, "toarray") else np.asarray(ad_sub.X)
    df = pd.DataFrame(X, columns=keep, index=ad_sub.obs_names)
    df["__grp__"] = ad_sub.obs["__grp__"].values
    df["__ct__"]  = ad_sub.obs[ct_col].astype(str).values
    out = []
    for ct, sub in df.groupby("__ct__"):
        ms = sub[sub["__grp__"] == "MS"]; hc = sub[sub["__grp__"] == "HC"]
        if len(ms) < 5 or len(hc) < 5: continue
        for g in keep:
            x = ms[g].values; y = hc[g].values
            try:
                p = mannwhitneyu(x, y, alternative="two-sided").pvalue
            except Exception:
                p = 1.0
            out.append(dict(
                gene=g, cell_type=ct,
                n_MS=len(ms), n_HC=len(hc),
                mean_MS=float(np.mean(x)), mean_HC=float(np.mean(y)),
                logFC=float(np.mean(x) - np.mean(y)),
                cohens_d=cohens_d(x, y),
                pct_MS=float(np.mean(x > 0) * 100),
                pct_HC=float(np.mean(y > 0) * 100),
                wilcox_p=p))
    o = pd.DataFrame(out)
    if not o.empty:
        from statsmodels.stats.multitest import fdrcorrection
        o["fdr"] = fdrcorrection(o["wilcox_p"].values, method="indep")[1]
    return o


# Run on every dataset
all_long = []
for name, path, ct_col in DATASETS:
    if not path.exists():
        print(f"  SKIP {name} (missing {path})"); continue
    print(f"\n>>> {name}  ({path.stat().st_size/1e6:.0f} MB)")
    try:
        ad = sc.read_h5ad(path)
    except Exception as e:
        print(f"  load failed: {e}"); continue
    print(f"    cells={ad.n_obs}  genes={ad.n_vars}  obs_cols={list(ad.obs.columns)[:8]}")
    grp_col = find_col(ad.obs, ["condition","disease","group","Group","Diagnosis","Case_Control"])
    if grp_col is None:
        print("  no group col"); continue
    ct = ct_col if ct_col in ad.obs.columns else find_col(ad.obs, ["cell_type","celltype","leiden","seurat_clusters"])
    if ct is None:
        print("  no cell-type col"); continue
    res = per_celltype_de(ad, grp_col, ct, GENE_PANEL)
    if res is None or res.empty: continue
    res["dataset"] = name
    all_long.append(res)
    print(f"    tests: {len(res)}  sig (FDR<0.05): {(res['fdr']<0.05).sum()}")

long = pd.concat(all_long, ignore_index=True) if all_long else pd.DataFrame()
print(f"\nTotal long rows: {len(long)}")
long.to_csv(OUT_DIR / "INV_scRNA_validation_long.tsv", sep="\t", index=False)


# Per-gene aggregate
if not long.empty:
    by_gene = (long.groupby("gene")
                   .apply(lambda d: pd.Series({
                       "n_tests": len(d),
                       "n_sig_FDR05": int((d["fdr"] < 0.05).sum()),
                       "n_celltypes_sig": d.loc[d["fdr"] < 0.05, "cell_type"].nunique(),
                       "n_datasets_sig": d.loc[d["fdr"] < 0.05, "dataset"].nunique(),
                       "best_dataset": d.loc[d["fdr"].idxmin(), "dataset"],
                       "best_celltype": d.loc[d["fdr"].idxmin(), "cell_type"],
                       "best_logFC": d.loc[d["fdr"].idxmin(), "logFC"],
                       "best_d": d.loc[d["fdr"].idxmin(), "cohens_d"],
                       "best_FDR": float(d["fdr"].min()),
                       "sig_celltypes": "; ".join(sorted(set(d.loc[d["fdr"] < 0.05, "cell_type"])))[:200],
                   }))
                   .reset_index()
                   .sort_values(["n_sig_FDR05","best_FDR"], ascending=[False, True]))
    by_gene.to_csv(OUT_DIR / "INV_scRNA_validation_by_gene.tsv", sep="\t", index=False)
    print(f"By-gene table: {len(by_gene)} rows")
    print(by_gene.head(20).to_string(index=False))
else:
    print("No long results.")


# Heatmap (gene × cell_type) — focus on Beltran blood where almost all signal lives
if long.empty:
    print("nothing to plot")
else:
    import matplotlib.pyplot as plt
    try:
        import seaborn as sns
    except ImportError:
        import subprocess, sys as _sys
        subprocess.check_call([_sys.executable, "-m", "pip", "install", "seaborn"])
        import seaborn as sns

    # ---- Plot A: gene × (dataset, cell_type) heatmap (compact) ----
    # For each (gene, dataset, cell_type) keep its single row (no aggregation)
    sig_only = long[long["fdr"] < 0.05].copy()
    sig_only["dataset_ct"] = sig_only["dataset"].str.replace(" ", "_") + "::" + sig_only["cell_type"].astype(str)
    if len(sig_only) == 0:
        print("no sig hits to plot")
    else:
        pivot_d   = sig_only.pivot_table(index="gene", columns="dataset_ct",
                                         values="cohens_d", aggfunc="first")
        pivot_fdr = sig_only.pivot_table(index="gene", columns="dataset_ct",
                                         values="fdr", aggfunc="first")
        # Order genes by n_sig_FDR05
        gene_ord = by_gene[by_gene["n_sig_FDR05"] > 0].sort_values(
            ["n_sig_FDR05","best_FDR"], ascending=[False, True]).gene.tolist()
        gene_ord = [g for g in gene_ord if g in pivot_d.index]
        pivot_d   = pivot_d.loc[gene_ord]
        pivot_fdr = pivot_fdr.loc[gene_ord]
        annot = pivot_fdr.applymap(lambda x: "***" if pd.notna(x) and x < 0.001 else
                                            ("**" if pd.notna(x) and x < 0.01 else
                                             ("*" if pd.notna(x) and x < 0.05 else "")))
        fig, ax = plt.subplots(figsize=(max(8, 0.45*pivot_d.shape[1]),
                                          max(5, 0.4*pivot_d.shape[0])), dpi=160)
        sns.heatmap(pivot_d.fillna(0), cmap="RdBu_r", center=0,
                    vmin=-1, vmax=1, ax=ax,
                    annot=annot, fmt="", linewidths=0.3,
                    cbar_kws={"label": "Cohen's d (MS - HC)"})
        ax.set_title("Inverse-concordant genes — scRNA-seq validation\n(only FDR<0.05 cells shown)")
        ax.set_xlabel(""); ax.set_ylabel("")
        plt.xticks(rotation=45, ha="right", fontsize=8)
        plt.yticks(fontsize=9)
        plt.tight_layout()
        plt.savefig(FIG_DIR / "11_inverse_scRNA_heatmap.png", dpi=160, bbox_inches="tight")
        plt.show()
        print("Wrote 11_inverse_scRNA_heatmap.png")

    # ---- Plot B: per-dataset sig counts barplot ----
    pivot_ct = (long[long["fdr"] < 0.05]
                .groupby(["dataset","cell_type"]).size().reset_index(name="n_sig"))
    pivot_ct = pivot_ct.sort_values("n_sig", ascending=False).head(25)
    if len(pivot_ct) > 0:
        fig, ax = plt.subplots(figsize=(10, 5), dpi=160)
        sns.barplot(data=pivot_ct, y=pivot_ct.dataset + "::" + pivot_ct.cell_type,
                    x="n_sig", color="#1F4E79", ax=ax)
        ax.set_title("Top (dataset × cell_type) with sig inverse-concordant hits")
        ax.set_xlabel("# inverse-concordant genes sig (FDR<0.05)"); ax.set_ylabel("")
        plt.tight_layout()
        plt.savefig(FIG_DIR / "11_inverse_scRNA_celltype_counts.png", dpi=160, bbox_inches="tight")
        plt.show()

