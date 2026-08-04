#!/usr/bin/env python3
# 12_celltype_4layer_master_py  —  generated from notebook spec


# ============================================================
# # 12 — Cell-type-specific scRNA validation + 4-layer master summary
# 
# For every gene of interest (CO7 panel + top inverse-concordant), produce:
# 
# 1. **Cell-type-specific scRNA breakdown** across all 7 scRNA datasets
#    (Schafflick PBMC, Schafflick paperUMAP, Absinta brain, Jakel brain,
#    Beltran blood, Ramesh blood, VitD blood) showing per-cell-type
#    Cohen's d + FDR — not just the "best" cell type as in nb 11.
# 
# 2. **4-layer integrated master table** per gene with:
#    - RNA: best stratum (PBMC / T cells / Brain WM / IFN-β / Pan-tissue / etc.)
#    - Methylation: best gene-level stratum + mCSEA promoter NES
#    - Proteomics: best assay (CSF Astral/timsTOF/combined, T-lineage, Magliozzi 4)
#    - scRNA: significant cell types per dataset
# 
# **Outputs**
# - `results/MASTER_4layer_validation.tsv`
# - `figures/12_celltype_scRNA_heatmap_CO7.png`
# - `figures/12_celltype_scRNA_heatmap_INVERSE.png`
# - `figures/12_master_4layer_overview.png`
# 
# Uses py_scrna kernel (uses already-saved TSVs from nb 01-11).
# ============================================================

import os, glob
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

PROJ = Path("__MS_GEO_ROOT__")
TX   = PROJ / "Transcriptome" / "results"
ME   = PROJ / "Methylation"   / "results"
PT   = PROJ / "Proteomics"    / "processed" / "META"
FIG  = PROJ / "Methylation"   / "figures"
SC_FIG = PROJ / "SingleCell_CELLxGENE" / "results" / "figures"

CO7 = ["LXN","SH3BP4","CHL1","CTSZ","RPAP2","PCNP","THRB"]

# ----- Top inverse-concordant from nb09 -----
inv = pd.read_csv(ME / "INVERSE_CONCORDANT_by_gene.tsv", sep="\t")
INV_TOP = inv.sort_values(["n_pairings","best_rna_fdr"], ascending=[False, True]).gene.tolist()[:25]
print(f"CO7: {CO7}")
print(f"INV_TOP (top 25): {INV_TOP[:15]}...")


# ----- Pull per-cell-type scRNA stats for all genes -----
# Sources: (A) nb11's INV_scRNA_validation_long.tsv (40 inverse genes)
#          (B) prior SingleCell_CELLxGENE stats CSVs (6 CO7 genes only, but
#              cover more datasets/cell types)
sc_rows = []

# (A) nb11 long table — for inverse genes + LXN
if (ME / "INV_scRNA_validation_long.tsv").exists():
    d = pd.read_csv(ME / "INV_scRNA_validation_long.tsv", sep="\t")
    d["source"] = "nb11_scan"
    d = d.rename(columns={"wilcox_p": "pval"})
    sc_rows.append(d[["gene","dataset","cell_type","n_MS","n_HC",
                      "logFC","cohens_d","pval","fdr","source"]])

# (B) prior CSVs — only 6 CO7 genes
prior_files = [
    ("Absinta brain (Schafflick paperUMAP)",
     SC_FIG / "brain/stats_per_celltype_MSvsHC.csv"),
    ("Schafflick blood",
     SC_FIG / "blood/stats_per_celltype_MSvsHC.csv"),
    ("Jakel brain",
     SC_FIG / "brain_Jakel2019/stats_per_celltype_MSvsHC.csv"),
    ("Beltran blood",
     SC_FIG / "blood_Beltran2019/stats_per_celltype_MSvsHC.csv"),
    ("Ramesh blood",
     SC_FIG / "blood_Ramesh2020/stats_per_celltype_MSvsHC.csv"
     if (SC_FIG / "blood_Ramesh2020/stats_per_celltype_MSvsHC.csv").exists() else None),
    ("VitD blood",
     SC_FIG / "blood_VitD_MS_GSE239626/stats_per_celltype_VitDvsPlac.csv"),
]
for tag, fp in prior_files:
    if fp is None or not fp.exists(): continue
    d = pd.read_csv(fp)
    ct_col = "cell_type" if "cell_type" in d.columns else ("celltype" if "celltype" in d.columns else None)
    if ct_col is None: continue
    d["dataset"] = tag
    d["source"]  = "prior_scRNA"
    d = d.rename(columns={
        "logfc": "logFC", "cohens_d": "cohens_d",
        "wilcoxon_p": "pval", "fdr": "fdr",
        ct_col: "cell_type"})
    if "n_MS" not in d.columns: d["n_MS"] = np.nan
    if "n_HC" not in d.columns: d["n_HC"] = np.nan
    sc_rows.append(d[["gene","dataset","cell_type","n_MS","n_HC",
                      "logFC","cohens_d","pval","fdr","source"]])

sc_long = pd.concat(sc_rows, ignore_index=True)
print(f"Combined scRNA stats: {len(sc_long)} rows across {sc_long.gene.nunique()} genes, "
      f"{sc_long.dataset.nunique()} datasets, "
      f"{sc_long['cell_type'].astype(str).nunique()} cell types")
sc_long.to_csv(ME / "scRNA_all_stats_combined.tsv", sep="\t", index=False)


# ----- Cell-type heatmap function -----
def celltype_heatmap(genes, title, out_png, max_cells=30):
    sub = sc_long[sc_long["gene"].isin(genes)].copy()
    sub = sub.dropna(subset=["fdr"])  # idxmin can't handle NaN
    if sub.empty:
        print(f"  no rows for genes={genes[:5]}..."); return None, None
    sub["dataset_ct"] = sub["dataset"].str.replace(" ", "_") + "::" + sub["cell_type"].astype(str)
    # Per (gene, dataset_ct) take best (smallest FDR) row
    sub = sub.loc[sub.groupby(["gene","dataset_ct"])["fdr"].idxmin()].copy()
    # Filter to cell types where at least one gene is sig (FDR<0.05) — keeps figure compact
    sig_cts = sub.loc[sub["fdr"] < 0.05, "dataset_ct"].unique()
    keep_ct = [c for c in sub["dataset_ct"].unique() if c in sig_cts]
    if len(keep_ct) == 0:
        # fallback: keep all cell types
        keep_ct = sub["dataset_ct"].unique().tolist()
    sub = sub[sub["dataset_ct"].isin(keep_ct)]
    # Heatmap values: cohens_d (only show cells with at least p<0.05 in lighter colour)
    pivot_d   = sub.pivot_table(index="gene", columns="dataset_ct",
                                values="cohens_d", aggfunc="first")
    pivot_fdr = sub.pivot_table(index="gene", columns="dataset_ct",
                                values="fdr", aggfunc="first")
    # order genes by sig count
    order = (sub[sub["fdr"] < 0.05]
             .groupby("gene").size().sort_values(ascending=False).index.tolist())
    rest = [g for g in genes if g not in order and g in pivot_d.index]
    order = order + rest
    pivot_d   = pivot_d.loc[order]
    pivot_fdr = pivot_fdr.loc[order]
    annot = pivot_fdr.applymap(lambda x: "***" if pd.notna(x) and x < 0.001 else
                                       ("**" if pd.notna(x) and x < 0.01 else
                                        ("*" if pd.notna(x) and x < 0.05 else "")))
    h = max(3, 0.35 * len(order))
    w = max(8, 0.4 * pivot_d.shape[1])
    fig, ax = plt.subplots(figsize=(w, h), dpi=160)
    sns.heatmap(pivot_d.fillna(0), cmap="RdBu_r", center=0,
                vmin=-1, vmax=1, ax=ax,
                annot=annot, fmt="", linewidths=0.3,
                cbar_kws={"label": "Cohen's d (MS - HC)"})
    ax.set_title(title)
    ax.set_xlabel(""); ax.set_ylabel("")
    plt.xticks(rotation=45, ha="right", fontsize=7)
    plt.yticks(fontsize=9)
    plt.tight_layout()
    plt.savefig(out_png, dpi=160, bbox_inches="tight")
    plt.show()
    print(f"Wrote {out_png}")
    return pivot_d, pivot_fdr

print(">>> CO7 panel cell-type heatmap")
_, _ = celltype_heatmap(CO7, "CO7 panel — cell-type-specific scRNA validation",
                         FIG / "12_celltype_scRNA_heatmap_CO7.png")

print("\n>>> Top inverse-concordant cell-type heatmap")
_, _ = celltype_heatmap(INV_TOP, "Top inverse-concordant genes — cell-type scRNA validation",
                         FIG / "12_celltype_scRNA_heatmap_INVERSE.png")


# ----- 4-layer master summary -----
def load_tx(fp, name, gene_col="gene", fc_col="logFC", p_col="adj.P.Val"):
    if not Path(fp).exists(): return None
    d = pd.read_csv(fp, sep="\t")
    if gene_col not in d.columns and "Gene" in d.columns: d["gene"] = d["Gene"]
    return d[[gene_col, fc_col, p_col]].rename(columns={fc_col:"logFC", p_col:"fdr"}).assign(assay=name)

# ---- RNA layer ----
rna_assays = [
    (TX / "01_pbmc_DE.tsv", "PBMC"),
    (TX / "02_tcells_DE.tsv", "T cells"),
    (TX / "03_bcells_DE.tsv", "B cells"),
    (TX / "04_brainwm_DE.tsv", "Brain WM"),
    (TX / "05_whole_blood_DE.tsv", "Whole blood"),
    (TX / "06_pbmc_ifnb_DE.tsv", "IFN-b PBMC"),
    (TX / "07_pan_tissue_DE.tsv", "Pan-tissue"),
]
rna_long = pd.concat([load_tx(fp, n) for fp, n in rna_assays if fp.exists()],
                     ignore_index=True)
rna_long["layer"] = "RNA"

# ---- Methylation gene-level ----
meth_assays = [
    (ME / "01_tcells_meth_gene.tsv", "Meth T cells"),
    (ME / "02_wb_dmf_meth_gene.tsv", "Meth WB DMF"),
    (ME / "03_wb_ocrelizumab_meth_gene.tsv", "Meth WB Ocrelizumab"),
    (ME / "04_tcells_remission_meth_gene.tsv", "Meth T cells remission"),
    (ME / "05_combined_meth_gene.tsv", "Meth combined cohort"),
]
meth_rows = []
for fp, n in meth_assays:
    if not fp.exists(): continue
    d = pd.read_csv(fp, sep="\t")
    meth_rows.append(d[["gene","mean_logFC","adj.P.Val"]].rename(
        columns={"mean_logFC":"logFC","adj.P.Val":"fdr"}).assign(assay=n))
meth_long = pd.concat(meth_rows, ignore_index=True)
meth_long["layer"] = "Methylation"

# mCSEA promoter
mc = pd.read_csv(ME / "06_mCSEA_promoter.tsv", sep="\t")
mc_long = mc[["gene","NES","padj"]].rename(columns={"NES":"logFC","padj":"fdr"}).assign(
    assay="mCSEA promoter (combined)", layer="Methylation")
meth_long = pd.concat([meth_long, mc_long], ignore_index=True)

# ---- Proteomics ----
prot_assays = [
    (PT / "CSF_Astral_CC_results.tsv", "Prot CSF Astral"),
    (PT / "CSF_timsTOF_CC_results.tsv", "Prot CSF timsTOF"),
    # EXCLUDED under the no-imputation decision: CSF_combined_R_ComBat_DE.tsv is the only
    # remaining table built on MinProb-imputed input (03_csf_cross_platform_meta.R).
#     (PT / "CSF_combined_R_ComBat_DE.tsv", "Prot CSF combined"),
    (PT / "T_lineage_R_combined_DE.tsv", "Prot T-lineage meta"),
    (PT / "Pegram_R_DE_gene.tsv", "Prot Pegram NK8"),
    (PT / "Magliozzi_CC_MS_CTX_vs_ODC_CTX.tsv", "Prot Brain CTX"),
    (PT / "Magliozzi_CC_MS_NAWM_vs_ODC_WM.tsv", "Prot Brain NAWM"),
    (PT / "Magliozzi_CC_MS_WML_vs_ODC_WM.tsv", "Prot Brain WML-vs-WM"),
    (PT / "Magliozzi_CC_MS_WML_vs_MS_NAWM.tsv", "Prot Brain WML-vs-NAWM"),
]
prot_long = pd.concat([load_tx(fp, n) for fp, n in prot_assays if fp.exists()],
                      ignore_index=True)
prot_long["layer"] = "Proteomics"

all_layers = pd.concat([rna_long, meth_long, prot_long], ignore_index=True)
print(f"All-layer long table: {len(all_layers)} rows  · "
      f"{all_layers.gene.nunique()} genes  · "
      f"{all_layers.assay.nunique()} assays  · "
      f"layers: {sorted(all_layers.layer.unique())}")


# ----- Build per-gene 4-layer master summary -----
PANEL = list(dict.fromkeys(CO7 + INV_TOP))

def best(d, group_col, only_sig=False):
    if only_sig: d = d[d["fdr"] < 0.05]
    if d.empty: return None
    r = d.loc[d["fdr"].idxmin()]
    return f"{r['assay']} | logFC={r['logFC']:+.2f} | FDR={r['fdr']:.1e}"

def sig_summary(d, label_layer):
    sig = d[d["fdr"] < 0.05]
    n = len(sig)
    return n, "; ".join(sig["assay"].tolist())

master = []
for g in PANEL:
    row = {"gene": g, "in_CO7": g in CO7}
    # RNA
    rna = all_layers[(all_layers["gene"] == g) & (all_layers["layer"] == "RNA")]
    sig_rna = rna[rna["fdr"] < 0.05]
    row["RNA_n_sig"] = len(sig_rna)
    row["RNA_sig_strata"] = "; ".join(sig_rna["assay"].tolist())
    row["RNA_best"] = best(rna, "assay")
    # Methylation
    me = all_layers[(all_layers["gene"] == g) & (all_layers["layer"] == "Methylation")]
    sig_me = me[me["fdr"] < 0.05]
    row["Meth_n_sig"] = len(sig_me)
    row["Meth_sig_strata"] = "; ".join(sig_me["assay"].tolist())
    row["Meth_best"] = best(me, "assay")
    # Proteomics
    pr = all_layers[(all_layers["gene"] == g) & (all_layers["layer"] == "Proteomics")]
    sig_pr = pr[pr["fdr"] < 0.05]
    row["Prot_n_sig"] = len(sig_pr)
    row["Prot_sig_assays"] = "; ".join(sig_pr["assay"].tolist())
    row["Prot_best"] = best(pr, "assay")
    # scRNA
    sc = sc_long[sc_long["gene"] == g]
    sc_sig = sc[sc["fdr"] < 0.05]
    row["scRNA_n_sig"] = len(sc_sig)
    row["scRNA_sig_celltypes"] = "; ".join(
        (sc_sig["dataset"].str.replace(" ","_") + "::" + sc_sig["cell_type"].astype(str)).tolist()
    )
    sc_valid = sc.dropna(subset=["fdr"])
    sc_best = sc_valid.loc[sc_valid["fdr"].idxmin()] if not sc_valid.empty else None
    if sc_best is not None and not pd.isna(sc_best.get("fdr")):
        row["scRNA_best"] = f"{sc_best['dataset']}::{sc_best['cell_type']} | d={sc_best['cohens_d']:+.2f} | FDR={sc_best['fdr']:.1e}"
    else:
        row["scRNA_best"] = ""
    # Total layers with at least 1 sig
    row["layers_with_sig"] = sum(row[k] > 0 for k in
                                  ["RNA_n_sig","Meth_n_sig","Prot_n_sig","scRNA_n_sig"])
    master.append(row)

master_df = pd.DataFrame(master)
master_df = master_df.sort_values(["in_CO7","layers_with_sig","RNA_n_sig","Meth_n_sig","Prot_n_sig","scRNA_n_sig"],
                                   ascending=[False, False, False, False, False, False])
master_df.to_csv(ME / "MASTER_4layer_validation.tsv", sep="\t", index=False)
print(f"Wrote MASTER_4layer_validation.tsv  ({len(master_df)} genes × 4 layers)")
print()
cols_compact = ["gene","in_CO7","RNA_n_sig","Meth_n_sig","Prot_n_sig","scRNA_n_sig","layers_with_sig"]
print(master_df[cols_compact].to_string(index=False))


# ----- 4-layer overview heatmap -----
panel_ord = master_df["gene"].tolist()
M = master_df.set_index("gene")[["RNA_n_sig","Meth_n_sig","Prot_n_sig","scRNA_n_sig"]]

fig, ax = plt.subplots(figsize=(7, max(5, 0.3*len(panel_ord))), dpi=160)
sns.heatmap(M, cmap="YlOrRd", annot=True, fmt="d",
            linewidths=0.4, ax=ax,
            cbar_kws={"label": "# sig assays (FDR<0.05)"})
ax.set_title("4-layer validation matrix — CO7 + top inverse-concordant\n(rows ordered: CO7 first, then by total layer coverage)")
ax.set_xlabel(""); ax.set_ylabel("")
plt.yticks(fontsize=9)
plt.tight_layout()
plt.savefig(FIG / "12_master_4layer_overview.png", dpi=160, bbox_inches="tight")
plt.show()


# ----- Pretty-print per-gene fact sheets for the top 15 -----
print("="*80)
print("PER-GENE 4-LAYER FACT SHEETS (top 15 by layers_with_sig + CO7)")
print("="*80)
for _, r in master_df.head(15).iterrows():
    print(f"\n[{r['gene']}]  CO7={r['in_CO7']}  layers={r['layers_with_sig']}/4")
    if r['RNA_n_sig']  > 0: print(f"  RNA  ({r['RNA_n_sig']}): {r['RNA_sig_strata']}")
    if r['Meth_n_sig'] > 0: print(f"  Meth ({r['Meth_n_sig']}): {r['Meth_sig_strata']}")
    if r['Prot_n_sig'] > 0: print(f"  Prot ({r['Prot_n_sig']}): {r['Prot_sig_assays']}")
    if r['scRNA_n_sig']> 0: print(f"  scRNA({r['scRNA_n_sig']}): {r['scRNA_sig_celltypes']}")
    if r['RNA_best']:   print(f"    ► best RNA:   {r['RNA_best']}")
    if r['Meth_best']:  print(f"    ► best Meth:  {r['Meth_best']}")
    if r['Prot_best']:  print(f"    ► best Prot:  {r['Prot_best']}")
    if r['scRNA_best']: print(f"    ► best scRNA: {r['scRNA_best']}")

