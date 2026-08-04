#!/usr/bin/env python3
"""figure2_rna_volcanoes.py — per-stratum RNA volcanoes with INV-only Tier-1 highlighting.

Updates the prior per_celltype_volcanoes_v3.png by:
  • Re-colouring labels into the current evidence groups
  • Inverse-concordant Tier-1 (2 genes: ITGB2, IKZF1) get RED bold * + thick ring
  • Tier-2 auxiliary inverse-concordant (11 genes, including CD79B and HLA-E) get ORANGE ring
  • Tier-2 non-concordant proteomic anchors are PURPLE
  • Force-label both inverse-concordant Tier-1 genes in every volcano where they have signal
  • Stratum panels: PBMC, T cells, IFN-β PBMC, B cells, Whole blood, Brain WM

Output: figures/per_celltype_volcanoes_INV.png
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from adjustText import adjust_text
import warnings; warnings.filterwarnings("ignore")

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.titleweight": "bold",
    "axes.titlesize": 10.5,
    # Axis labels and ticks were 8.5/7.5, far too small next to the 28.8 pt panel letters and the
    # 11-14 pt gene labels, and illegible once the figure is scaled to column width. Raised to the
    # same scale Figure 3 uses for its axis labels so the two figures match.
    "axes.labelsize": 15.2,
    "xtick.labelsize": 13.0,
    "ytick.labelsize": 13.0,
    "axes.spines.top": False, "axes.spines.right": False,
})

DGE_ROOT = Path("__MS_GEO_ROOT__/Stratified_Analyses/Expression")
OUT = Path("__MS_GEO_ROOT__/Poster_v2/figures/per_celltype_volcanoes_INV.png")
SUMMARY = pd.read_csv("__MS_GEO_ROOT__/Stratified_Analyses/Expression_Subgroup_Results_Summary.csv")

# Map out_dir-derived subdir → human stratum label used in SUMMARY['subgroup']
SUBDIR2LABEL = {
    "cell_tissue_case_control_b_cells":      "B cells",
    "cell_tissue_case_control_brain_wm":     "Brain / WM",
    "cell_tissue_case_control_pbmc":         "PBMC",
    "cell_tissue_case_control_t_cells":      "T cells",
    "cell_tissue_case_control_whole_blood":  "Whole blood",
    "label_context_case_control_pbmc_ifnb":  "IFNb in PBMC",
}
SUMMARY_BY_DIR = SUMMARY.set_index("subgroup")

NAVY="#0D3B66"; TEAL="#3E92CC"; RED_HOT="#B71C1C"; SUG_TEAL="#00796B"; ORANGE="#E65100"; PURPLE="#6A1B9A"
GREY_DARK="#424242"; MUTE="#9E9E9E"

# Gene tier definitions
INV_TIER1 = {"ITGB2","IKZF1"}
SUGGESTIVE = set()
TIER2_AUX_INV = {"CD79B","CASP6","CASP8","DGKQ","MX1","IFIT1","NUP210","RUNX3","SH3BP4","LXN","HLA-E"}
# FOXP3 was removed from this group on revision: it is not quantified in ANY of the seven
# proteomic compartments (both CSF instruments, all four brain-region contrasts, UK Biobank-PPP),
# so it could not be a "strong proteomic candidate", which is what defines this group. It is
# retained in the STRING display as a canonical MS immune context gene, which is the role it
# actually plays (the IKZF1-RUNX3/FOXP3-STAT1-STAT3 axis).
TIER2_PROT = {"CTSZ","CHL1","ICAM1","ITGAL"}
ALL_TRACKED = INV_TIER1 | SUGGESTIVE | TIER2_AUX_INV | TIER2_PROT

def gene_style(g):
    if g in INV_TIER1:
        return dict(facecolor="#ffcdd2", edgecolor=RED_HOT, lw=2.6,
                     star="* ", fontsize=13.8, fontweight="bold",
                     marker_size=160, color=RED_HOT)
    if g in SUGGESTIVE:
        return dict(facecolor="#e0f2f1", edgecolor=SUG_TEAL, lw=2.0,
                     star="", fontsize=12.6, fontweight="bold",
                     marker_size=120, color=SUG_TEAL)
    if g in TIER2_AUX_INV:
        return dict(facecolor="#fff3e0", edgecolor=ORANGE, lw=1.6,
                     star="◆ ", fontsize=12.0, fontweight="bold",
                     marker_size=95, color=ORANGE)
    return dict(facecolor="#f3e5f5", edgecolor=PURPLE, lw=1.2,
                 star="", fontsize=11.2, fontweight="semibold",
                 marker_size=78, color=PURPLE)

FIGSIZE = (18.485, 12.7)  # 3 rows: A–F per-stratum + G combined pan-tissue

def load_dge(subdir):
    fp = DGE_ROOT / subdir / "DGE_Results_MS_vs_HC.csv"
    df = pd.read_csv(fp)
    df["adj.P.Val"] = pd.to_numeric(df["adj.P.Val"], errors="coerce")
    df["logFC"]    = pd.to_numeric(df["logFC"], errors="coerce")
    return df.dropna(subset=["logFC","adj.P.Val"])

def ns(subdir):
    label = SUBDIR2LABEL.get(subdir)
    if label is None or label not in SUMMARY_BY_DIR.index: return ""
    row = SUMMARY_BY_DIR.loc[label]
    return f"{int(row['ms_samples'])} MS / {int(row['hc_samples'])} HC"

def volcano(ax, df, title, sample_n, fdr_thr=0.05,
             color_up="#D62828", color_dn="#3E92CC"):
    df = df.copy()
    df["logp"] = -np.log10(df["adj.P.Val"].clip(lower=1e-300))
    sig = (df["adj.P.Val"] < fdr_thr)  # FDR-only for DEG count (manuscript-consistent)
    # background non-sig
    ax.scatter(df.loc[~sig, "logFC"], df.loc[~sig, "logp"],
                s=2.5, c="#dadada", alpha=0.55, rasterized=True, linewidths=0)
    up = sig & (df["logFC"] > 0); dn = sig & (df["logFC"] < 0)
    ax.scatter(df.loc[up, "logFC"], df.loc[up, "logp"], s=7,
                c=color_up, alpha=0.85, rasterized=True, linewidths=0)
    ax.scatter(df.loc[dn, "logFC"], df.loc[dn, "logp"], s=7,
                c=color_dn, alpha=0.85, rasterized=True, linewidths=0)

    # Label every tracked gene that passes FDR<0.05 in this stratum. Earlier versions capped the
    # Tier-2 lists at the four/three most significant, which silently dropped genes that were
    # significant here - LXN and RUNX3 disappeared from panel G that way. adjustText resolves
    # the resulting collisions.
    df_sig_only = df[df["adj.P.Val"] < fdr_thr]
    inv_rows = df_sig_only[df_sig_only.Gene.isin(INV_TIER1)].sort_values("adj.P.Val")
    sug_rows = df_sig_only[df_sig_only.Gene.isin(SUGGESTIVE)].sort_values("adj.P.Val")
    aux_rows = df_sig_only[df_sig_only.Gene.isin(TIER2_AUX_INV)].sort_values("adj.P.Val")
    prot_rows = df_sig_only[df_sig_only.Gene.isin(TIER2_PROT)].sort_values("adj.P.Val")
    cands = pd.concat([inv_rows, sug_rows, aux_rows, prot_rows]).drop_duplicates(subset=["Gene"]).sort_values("adj.P.Val")
    # no overall cap: every significant tracked gene is shown

    # Draw rings + create text objects (adjustText will reposition them)
    texts = []
    for r in cands.itertuples():
        g = r.Gene; x = r.logFC; ypos = r.logp
        st = gene_style(g)
        ax.scatter(x, ypos, s=st["marker_size"], edgecolor=st["edgecolor"],
                    facecolor="none", linewidth=st["lw"], zorder=5)
        t = ax.text(x, ypos, f"{st['star']}{g}",
                     fontsize=st["fontsize"], fontweight=st["fontweight"],
                     color=st["color"],
                     bbox=dict(boxstyle="round,pad=0.18",
                                 facecolor=st["facecolor"],
                                 edgecolor=st["edgecolor"], linewidth=0.6),
                     zorder=6)
        texts.append(t)
    if texts:
        adjust_text(texts, ax=ax,
                     expand=(1.6, 1.8),
                     arrowprops=dict(arrowstyle="-", color="#555", lw=0.7,
                                       shrinkA=3, shrinkB=4),
                     force_text=(0.5, 0.7), force_points=(0.3, 0.5),
                     time_lim=2)

    ax.axhline(-np.log10(fdr_thr), ls=":", color="#777", lw=0.5)
    # FC-magnitude guide lines removed: significance is BH-FDR < fdr_thr alone.
    nsig = int(sig.sum())
    ax.set_title(title, fontsize=28.8, fontweight="bold", loc="left", pad=4)
    ax.set_xlabel("log₂ fold-change (MS vs HC)")
    ax.set_ylabel("−log₁₀(BH-FDR)")

# ─── Figure layout: A–F per-stratum (2×3) + G combined pan-tissue (row 3) ───
import matplotlib.gridspec as gridspec
fig = plt.figure(figsize=FIGSIZE, dpi=300)
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.42, wspace=0.27)

panels = [
    ("cell_tissue_case_control_pbmc",    "PBMC"),
    ("cell_tissue_case_control_t_cells", "T cells"),
    ("label_context_case_control_pbmc_ifnb", "IFN-β PBMC"),
    ("cell_tissue_case_control_b_cells", "B cells"),
    ("cell_tissue_case_control_whole_blood", "Whole blood"),
    ("cell_tissue_case_control_brain_wm", "Brain WM"),
]
positions = [gs[0,0],gs[0,1],gs[0,2],gs[1,0],gs[1,1],gs[1,2]]
for _k,(pos,(subdir,title)) in enumerate(zip(positions, panels)):
    ax = fig.add_subplot(pos)
    df = load_dge(subdir)
    volcano(ax, df, chr(65+_k), ns(subdir))

# Panel G — combined pan-tissue DEG (13 datasets, 462 samples, tissue-adjusted).
# Source is 07_pan_tissue_DE.tsv, the output of 01_transcriptome/07_total_combined_de.R. This panel
# previously read Poster_v2/figures/COMBINED_pantissue_proper_DEG.csv, an artefact with no producer
# script that disagreed with the pipeline: under it IKZF1 was pan-tissue-significant (FDR 0.036) and
# HLA-E was not (0.295), whereas the inverse-concordance scan, the four-layer master, the unified
# assay table and Figure 6 all read 07_pan_tissue_DE.tsv, where IKZF1 is not significant (0.573) and
# HLA-E is (0.012). The pairing counts reported in the manuscript derive from 07, so 07 is canonical.
axG = fig.add_subplot(gs[2,1])
comb = pd.read_csv("__MS_GEO_ROOT__/Transcriptome/results/07_pan_tissue_DE.tsv", sep="\t")
comb = comb.rename(columns={"gene": "Gene"})
comb["adj.P.Val"] = pd.to_numeric(comb["adj.P.Val"], errors="coerce")
comb["logFC"]     = pd.to_numeric(comb["logFC"], errors="coerce")
comb = comb.dropna(subset=["logFC","adj.P.Val"])
volcano(axG, comb, "G", "")

# Add tier legend in bottom-left empty cell
# suptitle removed (bare A–G panel letters only)

from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor="#ffcdd2", edgecolor=RED_HOT, linewidth=2.0,
           label=f"Inverse-concordant Tier-1 ({len(INV_TIER1)} genes)"),
    Patch(facecolor="#fff3e0", edgecolor=ORANGE, linewidth=1.4,
           label=f"Tier-2 auxiliary inverse-concordant ({len(TIER2_AUX_INV)} genes)"),
    Patch(facecolor="#f3e5f5", edgecolor=PURPLE, linewidth=1.0,
           label=f"Tier-2 non-concordant proteomic anchors ({len(TIER2_PROT)} genes)"),
]
fig.legend(handles=legend_elements, loc="center", ncol=1,
            bbox_to_anchor=(0.170, 0.135), frameon=True, fontsize=13.6)
plt.subplots_adjust(left=0.05, right=0.97, top=0.96, bottom=0.05)
plt.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
print(f"✓ saved → {OUT}")
from PIL import Image
img = Image.open(OUT)
print(f"  → {img.size[0]} × {img.size[1]}")
