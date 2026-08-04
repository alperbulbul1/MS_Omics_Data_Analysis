#!/usr/bin/env python3
"""figure6_intersection_heatmap.py — POSTER-STYLE 4-layer × tissue
intersection matrix with Inverse-concordant Tier-1 highlighting.

Same data as v5_poster but with:
  • Rows fixed in four evidence blocks (Tier-1 / suggestive / Tier-2 auxiliary / proteomic anchor)
  • Horizontal separator lines between tiers
  • Gene-label coloring per tier (red bold / orange / grey)
  • Tier badge text on left edge

Columns (unchanged): 4 layer groups × 22 cohort-strata readouts
  RNA:    PBMC · T cells · B cells · BrainWM · WholeBlood · IFN-β PBMC · Pan-tissue
  Meth:   T cells · WB DMF · WB Ocrelizumab · Brain WM · Combined  (Figure-3 strata; mCSEA removed)
  Prot:   CSF Astral · CSF timsTOF · Brain CTX/NAWM/WML×2 · UKB plasma  (T-lineage microarray meta removed)
  scRNA:  Brain · Blood · CSF
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

PROJ = Path("__MS_GEO_ROOT__")
OUT  = PROJ / "Poster_v2" / "figures" / "intersection_heatmap_v_INV.png"

# Evidence groups.  HLA-E is a Tier-2 auxiliary gene: it satisfies both discovery arms
# (RNA FDR 0.031 in whole blood, methylation FDR 8.1e-8) but carries no proteomic or
# donor-level single-cell anchor, which is exactly the Tier-2 auxiliary definition.  It was
# previously displayed as a separate 'suggestive' category, which contradicted that rule.
INV_TIER1 = ["ITGB2","IKZF1"]
SUGGESTIVE = []
TIER2_AUX_INV = ["CD79B","LXN","HLA-E","SH3BP4","CASP6","CASP8","DGKQ","MX1","IFIT1","NUP210","RUNX3"]
# FOXP3 was removed from this group on revision: it is not quantified in ANY of the seven
# proteomic compartments (both CSF instruments, all four brain-region contrasts, UK Biobank-PPP),
# so it could not be a "strong proteomic candidate", which is what defines this group. It is
# retained in the STRING display as a canonical MS immune context gene, which is the role it
# actually plays (the IKZF1-RUNX3/FOXP3-STAT1-STAT3 axis).
TIER2_PROT_ANCHOR = ["CTSZ","CHL1","ICAM1","ITGAL"]

RED_HOT="#B71C1C"; SUG_TEAL="#00796B"; ORANGE="#E65100"; PURPLE="#6A1B9A"; GREY_DARK="#424242"

COLUMNS = [
    ("RNA",         "PBMC",                       "RNA-PBMC"),
    ("RNA",         "T cells",                    "RNA-T cells"),
    ("RNA",         "B cells",                    "RNA-B cells"),
    ("RNA",         "Brain WM",                   "RNA-BrainWM"),
    ("RNA",         "Whole blood",                "RNA-WholeBlood"),
    ("RNA",         "IFN-b PBMC",                 "RNA-IFN-β PBMC †"),
    ("RNA",         "Pan-tissue",                 "RNA-Pan-tissue"),
    ("Methylation", "Meth T cells",               "Meth-T cells"),
    ("Methylation", "Meth WB DMF",                "Meth-WB DMF"),
    ("Methylation", "Meth WB Ocrelizumab",        "Meth-WB Ocre"),
    ("Methylation", "Meth Brain WM",              "Meth-Brain WM"),
    ("Methylation", "Meth combined cohort",       "Meth-Combined"),
    # Region-level mCSEA promoter test. Added because it is the ONLY assay in which CASP8 reaches
    # methylation significance (FDR 0.019); without this column CASP8's methylation block was blank,
    # contradicting both Figure 3 panel F and the Results text, which names CASP8 among the six
    # mCSEA-promoter-retained genes. NOTE: this statistic is a directional region-level enrichment
    # score (range about +/-2.3), not a per-CpG beta log fold-change (range about +/-0.6), so under
    # the shared +/-1 clip its colour saturates and is NOT magnitude-comparable with the columns to
    # its left; read its sign and its asterisks only. Stated in the figure caption.
    ("Methylation", "mCSEA promoter (combined)",  "Meth-mCSEA prom"),
    ("Proteomics",  "raw:astral",                 "Prot-CSF Astral"),
    ("Proteomics",  "raw:timstof",                "Prot-CSF timsTOF"),
    ("Proteomics",  "mag:CTX",                    "Prot-Brain CTX"),
    ("Proteomics",  "mag:NAWM",                   "Prot-Brain NAWM"),
    ("Proteomics",  "mag:WMLWM",                  "Prot-Brain WML/WM"),
    ("Proteomics",  "mag:WMLNAWM",                "Prot-Brain WML/NAWM"),
    ("Proteomics",  "raw:ukb",                    "Prot-Blood UKB"),
    ("scRNA",       "brain",                      "scRNA-Brain"),
    ("scRNA",       "blood",                      "scRNA-Blood"),
    ("scRNA",       "csf",                        "scRNA-CSF"),   # data uses lowercase 'csf'
]
LAYER_BOUNDARIES = [7, 13, 20]  # RNA 0-6 | Meth 7-12 | Prot 13-19 | scRNA 20-22

# ── Load data ──────────────────────────────────────────────────────────────
ua = pd.read_csv(PROJ/"Methylation/results/Unified_All_Assays_Long.tsv", sep="\t")
sc = pd.read_csv(PROJ/"Methylation/results/INV_scRNA_per_gene_per_tissue.tsv", sep="\t")
inv = pd.read_csv(PROJ/"Methylation/results/INVERSE_CONCORDANT_by_gene.tsv", sep="\t")
inv["combined"] = -np.log10(inv.best_rna_fdr.clip(1e-300)) + \
                   -np.log10(inv.best_meth_fdr.clip(1e-300))

# Brain WM methylation (Huynh 2014) — external load from concordance file
brain_meth = pd.read_csv(PROJ/"Methylation/results/07_BrainWM_RNA_vs_Meth_concordance.tsv",
                          sep="\t")
brain_meth_idx = brain_meth.set_index("gene")

# ── RAW proteomic sources (authoritative; match Figure 5) ──────────────────
def _dedup(df, fc, fdr, gene):
    df = df.rename(columns={fc:"_fc", fdr:"_fdr", gene:"_g"}).dropna(subset=["_fc","_fdr","_g"])
    df["_g"]=df["_g"].astype(str).str.strip()
    return df.sort_values("_fdr").drop_duplicates("_g", keep="first").set_index("_g")
astral_idx  = _dedup(pd.read_csv(PROJ/"Proteomics/processed/RDEP_CC/Astral_RDEP.tsv",sep="\t"),
                      "log2FC_MSvsCtrl","FDR","Genes")
timstof_idx = _dedup(pd.read_csv(PROJ/"Proteomics/processed/RDEP_CC/timsTOF_RDEP.tsv",sep="\t"),
                      "log2FC","FDR","Genes")
ukb_idx     = _dedup(pd.read_csv(PROJ/"Proteomics/blood_raw/Jacobs2024_UKB_primary_DE.tsv",sep="\t"),
                      "beta","fdr","gene")
MAGFILE = {"CTX":"Magliozzi2026_MS_CTX_vs_ODC_CTX.tsv","NAWM":"Magliozzi2026_MS_NAWM_vs_ODC_WM.tsv",
           "WMLWM":"Magliozzi2026_MS_WML_vs_ODC_WM.tsv","WMLNAWM":"Magliozzi2026_MS_WML_vs_MS_NAWM.tsv"}
def _dedup_mag(df):  # brain n=8: retain nominal p (matches the Figure 4 significance criterion)
    df = df.rename(columns={"log2FC":"_fc","FDR":"_fdr","p":"_p","Gene":"_g"}).dropna(subset=["_fc","_p","_g"])
    df["_g"]=df["_g"].astype(str).str.strip()
    return df.sort_values("_p").drop_duplicates("_g", keep="first").set_index("_g")
mag_idx = {k:_dedup_mag(pd.read_csv(PROJ/f"Proteomics/processed/RDEP_CC/{v}",sep="\t")) for k,v in MAGFILE.items()}
PROT_RAW = {"raw:astral":astral_idx, "raw:timstof":timstof_idx, "raw:ukb":ukb_idx,
            "mag:CTX":mag_idx["CTX"], "mag:NAWM":mag_idx["NAWM"],
            "mag:WMLWM":mag_idx["WMLWM"], "mag:WMLNAWM":mag_idx["WMLNAWM"]}

# Every displayed row is explicitly discussed in the manuscript.  The strict tiered
# panel contains 17 tiered genes; there is no separate suggestive row.
genes = INV_TIER1 + SUGGESTIVE + TIER2_AUX_INV + TIER2_PROT_ANCHOR
print(f"Selected {len(genes)} displayed candidates (all tiered)")

# ── Build matrices ─────────────────────────────────────────────────────────
val_mat = np.full((len(genes), len(COLUMNS)), np.nan)
sig_mat = np.full((len(genes), len(COLUMNS)), "", dtype=object)
for i, g in enumerate(genes):
    for j, (layer, assay, _label) in enumerate(COLUMNS):
        if layer == "scRNA":
            sub = sc[(sc.gene == g) & (sc.tissue == assay)]
            if len(sub):
                v = float(sub.best_d.iloc[0])
                fdr = float(sub.best_FDR.iloc[0])
                val_mat[i, j] = v
                if fdr < 0.001: sig_mat[i, j] = "***"
                elif fdr < 0.01: sig_mat[i, j] = "**"
                elif fdr < 0.05: sig_mat[i, j] = "*"
        elif assay in PROT_RAW:
            # RAW proteomic source (Astral / timsTOF / UKB / Magliozzi regions)
            idx = PROT_RAW[assay]
            gkey = g if g in idx.index else (g.upper() if g.upper() in idx.index else None)
            if gkey is not None:
                row = idx.loc[gkey]
                v = float(row["_fc"]); val_mat[i, j] = v
                # brain (mag) → nominal p (n=8 underpowered, matches Figure 4); CSF/UKB → BH-FDR
                pv = float(row["_p"]) if assay.startswith("mag:") else float(row["_fdr"])
                if pv < 0.001: sig_mat[i, j] = "***"
                elif pv < 0.01: sig_mat[i, j] = "**"
                elif pv < 0.05: sig_mat[i, j] = "*"
        elif assay == "Meth Brain WM":
            # External: load from Brain WM concordance file
            if g in brain_meth_idx.index:
                row = brain_meth_idx.loc[g]
                v = float(row["meth"]); fdr = float(row["meth_fdr"])
                val_mat[i, j] = v
                if fdr < 0.001: sig_mat[i, j] = "***"
                elif fdr < 0.01: sig_mat[i, j] = "**"
                elif fdr < 0.05: sig_mat[i, j] = "*"
        else:
            sub = ua[(ua.gene == g) & (ua.layer == layer) & (ua.assay == assay)]
            if len(sub):
                row = sub.iloc[sub.logFC.abs().argmax()]
                v = float(row.logFC); fdr = float(row.fdr)
                val_mat[i, j] = v
                if fdr < 0.001: sig_mat[i, j] = "***"
                elif fdr < 0.01: sig_mat[i, j] = "**"
                elif fdr < 0.05: sig_mat[i, j] = "*"

# ── Plot ───────────────────────────────────────────────────────────────────
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':16.0})
fig, ax = plt.subplots(figsize=(20, 11.5), dpi=200)

cmap = LinearSegmentedColormap.from_list("RdBu_r_clip",
        ["#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0",
         "#FFFFFF",
         "#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"])
plot_mat = np.clip(val_mat, -1, 1)
from matplotlib.colors import Normalize
norm = Normalize(vmin=-1, vmax=1)
rgba = cmap(norm(np.nan_to_num(plot_mat, nan=0.0)))
GREY = np.array([0.93, 0.93, 0.93])
for i in range(len(genes)):
    for j in range(len(COLUMNS)):
        if np.isnan(val_mat[i, j]):
            rgba[i, j, :3] = 1.0                                    # no data -> white
        elif not sig_mat[i, j]:
            rgba[i, j, :3] = 0.32 * rgba[i, j, :3] + 0.68 * GREY    # non-significant -> faded (only significant cells keep full colour)
im = ax.imshow(rgba, aspect='auto')
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])

# Significance stars overlay
for i in range(len(genes)):
    for j in range(len(COLUMNS)):
        if sig_mat[i, j]:
            col = 'white' if abs(plot_mat[i, j]) > 0.55 else '#212121'
            ax.text(j, i, sig_mat[i, j], ha='center', va='center',
                     fontsize=14.4, fontweight='bold', color=col)

# Vertical layer separators
for sep in LAYER_BOUNDARIES:
    ax.axvline(sep - 0.5, color='black', linewidth=2.5)

# Horizontal tier separators
n_inv = sum(1 for g in genes if g in INV_TIER1)
n_sug = sum(1 for g in genes if g in SUGGESTIVE)
n_aux = sum(1 for g in genes if g in TIER2_AUX_INV)
n_prot = sum(1 for g in genes if g in TIER2_PROT_ANCHOR)

if n_inv > 0:
    ax.axhline(n_inv - 0.5, color="#000", linewidth=3.0, alpha=0.9)
if n_sug > 0:
    ax.axhline(n_inv + n_sug - 0.5, color="#000", linewidth=1.4,
                linestyle=":", alpha=0.8)
if n_aux > 0:
    ax.axhline(n_inv + n_sug + n_aux - 0.5, color="#000", linewidth=2.0,
                linestyle="--", alpha=0.8)

# Evidence-group badges on the far left.
ax.text(-2.7, (n_inv - 1) / 2, f"Inverse-concordant Tier-1\n({n_inv} genes)",
         ha='right', va='center', fontsize=18.4, fontweight='bold',
         color=RED_HOT, transform=ax.transData)
if n_sug > 0:
    ax.text(-2.7, n_inv + (n_sug-1)/2, "Suggestive\n(non-tier-1)",
             ha='right', va='center', fontsize=16.8, fontweight='bold',
             color=SUG_TEAL, transform=ax.transData)
ax.text(-2.7, n_inv + n_sug + n_aux/2 - 0.5, f"Tier-2\nauxiliary inverse-concordant\n({n_aux} genes)",
         ha='right', va='center', fontsize=16.8, fontweight='semibold',
         color=ORANGE, transform=ax.transData)
ax.text(-2.7, n_inv + n_sug + n_aux + n_prot/2 - 0.5,
         f"Tier-2\nnon-concordant\nproteomic\nanchors\n({n_prot} genes)",
         ha='right', va='center', fontsize=16.0,
         fontweight='semibold', color=PURPLE, transform=ax.transData)

# Column labels
ax.set_xticks(range(len(COLUMNS)))
ax.set_xticklabels([c[2] for c in COLUMNS], rotation=45, ha='right',
                    fontsize=16.0)

# Row labels (gene names, italic) with tier-specific coloring
ax.set_yticks(range(len(genes)))
ax.set_yticklabels(genes, fontsize=17.6, fontweight='bold', fontstyle='italic')
for ytl, g in zip(ax.get_yticklabels(), genes):
    if g in INV_TIER1:
        ytl.set_color(RED_HOT); ytl.set_fontweight('bold')
    elif g in SUGGESTIVE:
        ytl.set_color(SUG_TEAL); ytl.set_fontweight('bold')
    elif g in TIER2_AUX_INV:
        ytl.set_color(ORANGE); ytl.set_fontweight('semibold')
    elif g in TIER2_PROT_ANCHOR:
        ytl.set_color(PURPLE); ytl.set_fontweight('semibold')
    else:
        ytl.set_color(GREY_DARK)

# (title removed — descriptive text moved to the figure caption)

# Colorbar
cbar = plt.colorbar(sm, ax=ax, fraction=0.025, pad=0.01)
cbar.set_label("log₂ fold change  (clipped −1..+1)", fontsize=17.6)
cbar.ax.tick_params(labelsize=14.4)

# Light cell grid
for y in range(len(genes) + 1):
    ax.axhline(y - 0.5, color='#EEEEEE', linewidth=0.3, zorder=0)
for x in range(len(COLUMNS) + 1):
    if x - 0.5 not in [s - 0.5 for s in LAYER_BOUNDARIES]:
        ax.axvline(x - 0.5, color='#EEEEEE', linewidth=0.3, zorder=0)

# Top-margin column-group labels (above tick labels, below title)
group_labels = [("RNA bulk", 0, 6), ("Methylation", 7, 11),
                  ("Proteomics", 12, 18), ("scRNA", 19, 21)]
for lab, lo, hi in group_labels:
    ax.text((lo + hi) / 2, -1.1, lab, ha='center', va='center',
             fontsize=17.6, fontweight='bold', color='#333',
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#F5F5F5",
                         edgecolor="#888", linewidth=0.5),
             transform=ax.transData,
             clip_on=False)

plt.subplots_adjust(left=0.13)
plt.tight_layout()
plt.savefig(OUT, dpi=300, bbox_inches='tight', facecolor='white')
print(f"✓ saved → {OUT}")
from PIL import Image
img = Image.open(OUT)
print(f"  → {img.size[0]} × {img.size[1]} ({img.size[0]/img.size[1]:.2f}:1)")

print("\nGenes per tier:")
print(f"  Inverse-concordant Tier-1 ({n_inv}): {[g for g in genes if g in INV_TIER1]}")
print(f"  Suggestive non-tier-1 ({n_sug}): {[g for g in genes if g in SUGGESTIVE]}")
print(f"  Tier-2 auxiliary inverse-concordant ({n_aux}): {[g for g in genes if g in TIER2_AUX_INV]}")
print(f"  Tier-2 non-concordant prot anchors ({n_prot}): {[g for g in genes if g in TIER2_PROT_ANCHOR]}")
