#!/usr/bin/env python3
# 13_perstudy_scRNA_validation_py  —  generated from notebook spec


# ============================================================
# # 13 — Per-study scRNA validation (tissue + cell-type aware)
# 
# Re-runs the inverse-concordant + CO7 panel **inside each study's own
# proper context**:
# 
# - **Brain studies** (Absinta, Jakel) — test in brain cell types
#   (astrocytes, OPCs, oligodendrocytes, neurons, microglia, vascular)
#   stratified by lesion pathology where available.
# - **Blood studies** (Beltran, Ramesh) — test in blood cell types,
#   Beltran additionally split into **CSF vs PBMC compartment**.
# - **VitD MS** — *NOT* tested for MS-vs-HC (no HC; all subjects are MS);
#   excluded from this notebook.
# - **Schafflick PBMC** — h5ad is corrupt; reuse pre-computed CSV stats
#   from prior runs.
# 
# This corrects three issues in nb 11:
# 1. VitD was wrongly treated as MS-vs-HC (actually Placebo vs VitaminD).
# 2. Beltran lumped PBMC and CSF; separating them lets us see compartment-
#    specific effects.
# 3. Brain studies are now tested ONLY in brain cell types (and stratified
#    by lesion pathology), not pooled across mixed tissues.
# 
# **Outputs**
# - `results/INV_scRNA_per_study_long.tsv`
# - `results/INV_scRNA_per_study_by_gene.tsv`
# - `figures/13_perstudy_brain_heatmap.png`     (brain studies)
# - `figures/13_perstudy_blood_heatmap.png`     (blood studies)
# - `figures/13_perstudy_csf_heatmap.png`       (Beltran CSF)
# ============================================================

import os
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection

PROJ   = Path("__MS_GEO_ROOT__")
ME     = PROJ / "Methylation" / "results"
SC_FIG = PROJ / "SingleCell_CELLxGENE" / "results" / "figures"
SC_RAW = PROJ / "SingleCell_CELLxGENE" / "data"
FIG    = PROJ / "Methylation" / "figures"

CO7 = ["LXN","SH3BP4","CHL1","CTSZ","RPAP2","PCNP","THRB"]
# Tier-2 non-concordant proteomic anchors. They are displayed in the Figure 6 intersection
# matrix, which carries an scRNA block, so they must be tested there rather than left blank.
# The legacy CO7 list predates the current tiering and happened to cover only CTSZ and CHL1,
# which is why ICAM1, ITGAL and FOXP3 were silently missing from the scRNA columns.
PROT_ANCHORS = ["ICAM1","ITGAL","FOXP3"]
# Every gene displayed in Figure 6 must be tested here, or it shows an empty scRNA block that reads
# as "no effect" when it actually means "never tested". CO7 + the top-30 inverse-concordant cut
# covered 17 of the 18 by luck: CASP8 ranks 41st on that ordering (one pairing) and SH3BP4 31st,
# and only SH3BP4 happened to sit in the legacy CO7 list.
DISPLAY_PANEL = ["ITGB2","IKZF1","HLA-E","CD79B","LXN","SH3BP4","CASP6","CASP8","DGKQ","MX1",
                 "IFIT1","NUP210","RUNX3","CTSZ","CHL1","ICAM1","FOXP3","ITGAL"]
inv = pd.read_csv(ME / "INVERSE_CONCORDANT_by_gene.tsv", sep="\t")
INV_TOP = inv.sort_values(["n_pairings","best_rna_fdr"], ascending=[False, True]).gene.tolist()[:30]
GENE_PANEL = list(dict.fromkeys(CO7 + PROT_ANCHORS + DISPLAY_PANEL + INV_TOP))
print(f"Gene panel: {len(GENE_PANEL)} genes  (CO7 + top 30 inverse)")


# ----- Per-study config (the heart of the fix) -----
DATASETS = [
    dict(
        name="Absinta brain (snRNA WM)",
        tissue="brain",
        path=SC_FIG / "brain_paperUMAP" / "adata_paperclusters.h5ad",
        group_col="group", case_labels={"MS"}, ctrl_labels={"HC"},
        ct_col="cell_type",
        stratify_col="pathology",
        notes="Brain WM lesion snRNA. Stratify by pathology (lesion core/edge/NAWM/HC)."
    ),
    dict(
        name="Jakel brain (snRNA)",
        tissue="brain",
        path=SC_FIG / "brain_Jakel2019" / "adata_jakel.h5ad",
        group_col="condition", case_labels={"MS"}, ctrl_labels={"Ctrl"},
        ct_col="celltype",
        stratify_col="lesion",
        notes="Adult brain snRNA. 23 fine cell types, lesion-stratified."
    ),
    dict(
        name="Beltran PBMC (10x)",
        tissue="blood",
        path=SC_FIG / "blood_Beltran2019" / "adata_beltran.h5ad",
        group_col="Case", case_labels={"MS","MS_PBMCs"}, ctrl_labels={"HD","HD_PBMCs"},
        ct_col="celltype",
        stratify_col="compartment",
        compartment_filter={"CD4_PBMCs","CD8_PBMCs","CD19CD27_PBMCs","Others_PBMCs"},
        notes="PBMC compartment only (filtered from CSF)."
    ),
    dict(
        name="Beltran CSF (10x)",
        tissue="csf",
        path=SC_FIG / "blood_Beltran2019" / "adata_beltran.h5ad",
        group_col="Case", case_labels={"MS"}, ctrl_labels={"HD"},
        ct_col="celltype",
        stratify_col="compartment",
        compartment_filter={"CSF"},
        notes="CSF compartment only (filtered from PBMC)."
    ),
    dict(
        name="Ramesh PBMC",
        tissue="blood",
        path=SC_FIG / "blood_Ramesh2020_UMAP" / "adata_ramesh_small.h5ad",
        group_col="disease_status", case_labels={"MS"}, ctrl_labels={"HC"},
        ct_col="basictype",
        stratify_col=None,
        cohort_filter=lambda obs: ~obs["cohort"].astype(str).isin(["NAT_HI"]),  # drop natalizumab
        notes="PBMC, natalizumab subgroup filtered out."
    ),
]

for d in DATASETS:
    if not d["path"].exists():
        print(f"  [MISSING] {d['name']}  {d['path']}")
print(f"\nConfigured {len(DATASETS)} datasets.\n")


# ----- Helper: per-cell-type MS-vs-HC test with proper config -----
def cohens_d(a, b):
    na, nb = len(a), len(b)
    if na < 2 or nb < 2: return np.nan
    s = np.sqrt(((na-1)*np.var(a, ddof=1) + (nb-1)*np.var(b, ddof=1))/(na+nb-2))
    return (np.mean(a) - np.mean(b)) / s if s > 0 else 0.0

def run_study(cfg, genes):
    if not cfg["path"].exists():
        return pd.DataFrame()
    print(f"\n>>> {cfg['name']}  ({cfg['path'].name})")
    ad = sc.read_h5ad(cfg["path"])
    obs = ad.obs.copy()

    # Apply optional compartment filter (Beltran)
    if "compartment_filter" in cfg:
        mask = obs[cfg["stratify_col"]].astype(str).isin(cfg["compartment_filter"])
        print(f"  Compartment filter: {sum(mask)}/{len(mask)} cells kept ({cfg['compartment_filter']})")
        obs = obs[mask]
    # Apply optional cohort filter (Ramesh)
    if "cohort_filter" in cfg:
        mask = cfg["cohort_filter"](obs)
        print(f"  Cohort filter: {sum(mask)}/{len(mask)} cells kept")
        obs = obs[mask]

    # Assign MS/HC labels
    cond = obs[cfg["group_col"]].astype(str)
    grp = pd.Series(index=obs.index, dtype=str)
    grp[cond.isin(cfg["case_labels"])] = "MS"
    grp[cond.isin(cfg["ctrl_labels"])] = "HC"
    grp = grp.dropna()
    obs = obs.loc[grp.index].copy()
    obs["__grp__"] = grp
    print(f"  After group filter: {len(obs)} cells (MS={(obs['__grp__']=='MS').sum()} HC={(obs['__grp__']=='HC').sum()})")

    if obs.empty: return pd.DataFrame()
    ad_sub = ad[obs.index].copy()
    ad_sub.obs["__grp__"] = obs["__grp__"]
    if cfg.get("stratify_col"):
        ad_sub.obs["__strat__"] = obs[cfg["stratify_col"]].astype(str).values
    ad_sub.obs["__ct__"] = obs[cfg["ct_col"]].astype(str).values

    keep_genes = [g for g in genes if g in ad_sub.var_names]
    if not keep_genes: return pd.DataFrame()
    ad_sub = ad_sub[:, keep_genes]
    X = ad_sub.X.toarray() if hasattr(ad_sub.X, "toarray") else np.asarray(ad_sub.X)
    df = pd.DataFrame(X, columns=keep_genes, index=ad_sub.obs_names)
    df["__grp__"]   = ad_sub.obs["__grp__"].values
    df["__ct__"]    = ad_sub.obs["__ct__"].values
    if cfg.get("stratify_col"):
        df["__strat__"] = ad_sub.obs["__strat__"].values

    out = []
    # Pass 1: per cell-type (pooled over strata)
    for ct, sub in df.groupby("__ct__"):
        ms = sub[sub["__grp__"] == "MS"]; hc = sub[sub["__grp__"] == "HC"]
        if len(ms) < 5 or len(hc) < 5: continue
        for g in keep_genes:
            x = ms[g].values; y = hc[g].values
            try: p = mannwhitneyu(x, y, alternative="two-sided").pvalue
            except Exception: p = 1.0
            out.append(dict(
                dataset=cfg["name"], tissue=cfg["tissue"], stratum="ALL",
                cell_type=ct, gene=g, n_MS=len(ms), n_HC=len(hc),
                logFC=float(np.mean(x)-np.mean(y)),
                cohens_d=cohens_d(x, y),
                pct_MS=float(np.mean(x > 0)*100), pct_HC=float(np.mean(y > 0)*100),
                wilcox_p=p))
    # Pass 2: per cell-type × stratum (e.g. pathology / lesion)
    if cfg.get("stratify_col"):
        for (ct, st), sub in df.groupby(["__ct__","__strat__"]):
            ms = sub[sub["__grp__"] == "MS"]; hc = sub[sub["__grp__"] == "HC"]
            if len(ms) < 5 or len(hc) < 5: continue
            for g in keep_genes:
                x = ms[g].values; y = hc[g].values
                try: p = mannwhitneyu(x, y, alternative="two-sided").pvalue
                except Exception: p = 1.0
                out.append(dict(
                    dataset=cfg["name"], tissue=cfg["tissue"], stratum=str(st),
                    cell_type=ct, gene=g, n_MS=len(ms), n_HC=len(hc),
                    logFC=float(np.mean(x)-np.mean(y)),
                    cohens_d=cohens_d(x, y),
                    pct_MS=float(np.mean(x > 0)*100), pct_HC=float(np.mean(y > 0)*100),
                    wilcox_p=p))
    o = pd.DataFrame(out)
    if not o.empty:
        # BH-FDR within each (cell_type × stratum) group across genes
        o["fdr"] = np.nan
        for (ct, st), idx in o.groupby(["cell_type","stratum"]).groups.items():
            ps = o.loc[idx, "wilcox_p"].values
            if len(ps) > 0:
                o.loc[idx, "fdr"] = fdrcorrection(ps, method="indep")[1]
        print(f"  Tests: {len(o)}  sig (FDR<0.05): {(o['fdr']<0.05).sum()}")
    return o

# Run all studies
all_long = []
for cfg in DATASETS:
    o = run_study(cfg, GENE_PANEL)
    if not o.empty: all_long.append(o)

long = pd.concat(all_long, ignore_index=True) if all_long else pd.DataFrame()
print(f"\nTotal rows: {len(long)}  ·  unique studies: {long.dataset.nunique() if not long.empty else 0}")
long.to_csv(ME / "INV_scRNA_per_study_long.tsv", sep="\t", index=False)


# Per-gene per-study aggregate
if long.empty:
    print("no results")
else:
    by_gene_study = (long.dropna(subset=["fdr"])
                         .groupby(["gene","dataset","tissue"])
                         .apply(lambda d: pd.Series({
                             "n_tests": len(d),
                             "n_sig_FDR05": int((d["fdr"] < 0.05).sum()),
                             "best_ct": d.loc[d["fdr"].idxmin(), "cell_type"],
                             "best_stratum": d.loc[d["fdr"].idxmin(), "stratum"],
                             "best_d": d.loc[d["fdr"].idxmin(), "cohens_d"],
                             "best_logFC": d.loc[d["fdr"].idxmin(), "logFC"],
                             "best_FDR": float(d["fdr"].min()),
                             "sig_cts": "; ".join(sorted(set(
                                 d.loc[d["fdr"] < 0.05, "cell_type"] + "[" +
                                 d.loc[d["fdr"] < 0.05, "stratum"] + "]"
                             )))[:300],
                         })).reset_index())
    by_gene_study = by_gene_study.sort_values(["gene","tissue","best_FDR"])
    by_gene_study.to_csv(ME / "INV_scRNA_per_study_by_gene.tsv", sep="\t", index=False)
    print("By (gene,dataset) summary (top 30 sig):")
    print(by_gene_study[by_gene_study["n_sig_FDR05"] > 0].sort_values("n_sig_FDR05", ascending=False).head(30).to_string(index=False))


# ----- Tissue-stratified heatmaps -----
def tissue_heatmap(tissue, title, out_png, max_genes=30):
    sub = long[long["tissue"] == tissue].copy()
    if sub.empty: print(f"no data for {tissue}"); return
    sub = sub.dropna(subset=["fdr"])
    # Build "dataset::cell_type[stratum]" composite label
    sub["x_label"] = (sub["dataset"].str.replace(" ", "_") + "::" +
                      sub["cell_type"].astype(str) + "[" +
                      sub["stratum"].astype(str) + "]")
    # Per (gene, x_label) take row with smallest FDR
    sub = sub.loc[sub.groupby(["gene","x_label"])["fdr"].idxmin()].copy()
    # Keep cell-type columns where at least one gene is sig
    sig_x = sub.loc[sub["fdr"] < 0.05, "x_label"].unique()
    if len(sig_x) == 0:
        print(f"no sig hits in {tissue}"); return
    sub = sub[sub["x_label"].isin(sig_x)]
    pivot_d   = sub.pivot_table(index="gene", columns="x_label", values="cohens_d", aggfunc="first")
    pivot_fdr = sub.pivot_table(index="gene", columns="x_label", values="fdr", aggfunc="first")
    # Order rows by number of sig hits across columns
    order = (sub[sub["fdr"] < 0.05].groupby("gene").size()
             .sort_values(ascending=False).index.tolist())
    extra = [g for g in pivot_d.index if g not in order]
    order = (order + extra)[:max_genes]
    pivot_d = pivot_d.loc[order]; pivot_fdr = pivot_fdr.loc[order]
    annot = pivot_fdr.applymap(lambda x: "***" if pd.notna(x) and x < 0.001 else
                                       ("**" if pd.notna(x) and x < 0.01 else
                                        ("*" if pd.notna(x) and x < 0.05 else "")))
    w = max(8, 0.5 * pivot_d.shape[1])
    h = max(5, 0.4 * pivot_d.shape[0])
    fig, ax = plt.subplots(figsize=(w, h), dpi=150)
    sns.heatmap(pivot_d.fillna(0), cmap="RdBu_r", center=0, vmin=-1, vmax=1,
                ax=ax, annot=annot, fmt="", linewidths=0.3,
                cbar_kws={"label": "Cohen's d (MS - HC)"})
    ax.set_title(title)
    ax.set_xlabel(""); ax.set_ylabel("")
    plt.xticks(rotation=45, ha="right", fontsize=7)
    plt.yticks(fontsize=9)
    plt.tight_layout()
    plt.savefig(out_png, dpi=160, bbox_inches="tight")
    plt.show()
    print(f"Wrote {out_png}")

print(">>> BRAIN tissue heatmap")
tissue_heatmap("brain", "scRNA validation in BRAIN studies (Absinta + Jakel)\nwith pathology / lesion stratification",
                FIG / "13_perstudy_brain_heatmap.png")
print("\n>>> BLOOD tissue heatmap (PBMC)")
tissue_heatmap("blood", "scRNA validation in BLOOD studies (Beltran PBMC + Ramesh PBMC)",
                FIG / "13_perstudy_blood_heatmap.png")
print("\n>>> CSF tissue heatmap (Beltran)")
tissue_heatmap("csf", "scRNA validation in CSF (Beltran)",
                FIG / "13_perstudy_csf_heatmap.png")


# Per-tissue × per-gene summary table
if long.empty:
    print("no results")
else:
    tab = (long.dropna(subset=["fdr"])
              .groupby(["gene","tissue"])
              .apply(lambda d: pd.Series({
                  "n_tests": len(d),
                  "n_sig_FDR05": int((d["fdr"] < 0.05).sum()),
                  "best_dataset_ct": (d.loc[d["fdr"].idxmin(), "dataset"] + "::" +
                                       d.loc[d["fdr"].idxmin(), "cell_type"] + "[" +
                                       d.loc[d["fdr"].idxmin(), "stratum"] + "]"),
                  "best_d": d.loc[d["fdr"].idxmin(), "cohens_d"],
                  "best_FDR": float(d["fdr"].min()),
                  "sig_locations": "; ".join(sorted(set(
                      d.loc[d["fdr"] < 0.05, "dataset"] + "::" +
                      d.loc[d["fdr"] < 0.05, "cell_type"] + "[" +
                      d.loc[d["fdr"] < 0.05, "stratum"] + "]"
                  )))[:400]
              })).reset_index())
    # Pivot to wide: gene × tissue with sig counts
    wide = tab.pivot_table(index="gene", columns="tissue",
                           values="n_sig_FDR05", aggfunc="first", fill_value=0)
    wide["TOTAL_sig"] = wide.sum(axis=1)
    wide = wide.sort_values("TOTAL_sig", ascending=False)
    print("Per-gene × tissue sig count (top 30):")
    print(wide.head(30).to_string())
    tab.to_csv(ME / "INV_scRNA_per_gene_per_tissue.tsv", sep="\t", index=False)
    print(f"\nSaved INV_scRNA_per_gene_per_tissue.tsv")

