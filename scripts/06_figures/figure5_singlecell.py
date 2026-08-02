#!/usr/bin/env python3
"""figure5_singlecell.py — composite Figure 6 v3 with
Inverse-concordant Tier-1 priority highlighting.

Re-organises the heatmap with:
  - Inverse-concordant Tier-1 (6 genes) at TOP, label colour = red bold
  - Tier-2 auxiliary inverse-concordant (7 genes) middle, label colour = orange
  - Other aux candidates (non-concordant) at bottom, label colour = grey
  - Thick horizontal separator lines between tiers
  - Panel C gene UMAP overlays focus on Inverse-concordant Tier-1 (HLA-E, ITGB2, LXN, CD79B,
    IKZF1, SH3BP4) plus the auxiliary inverse-concordant gene CASP6
"""
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpecFromSubplotSpec
import numpy as np, pandas as pd
from scipy.sparse import issparse
from pathlib import Path
import warnings; warnings.filterwarnings('ignore')

FIG_DIR = Path("__MS_GEO_ROOT__/SingleCell_CELLxGENE/results/figures")
OUT     = Path("__MS_GEO_ROOT__/Poster_v2/figures/singlecell_combined_figure_v3_INV.png")
PB      = Path("__MS_GEO_ROOT__/Poster_v2/figures/scrna_WILCOXON_v1.tsv")

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':14.4})

# Tiered gene lists ─────────────────────────────────────────────────────────
INV_TIER1 = ["ITGB2","IKZF1"]
SUGGESTIVE = ["HLA-E"]   # demoted: suggestive, non-tier-1 (single-cell only)
TIER2_INV_AUX = ["LXN", "CD79B", "SH3BP4", "CASP6", "CASP8", "DGKQ", "MX1", "IFIT1", "NUP210", "RUNX3"]
TIER2_OTHER = ["CTSZ", "CHL1", "ITGAL", "IFI44L", "ICAM1", "FOXP3", "TYK2", "STAT3"]
ALL_GENES = INV_TIER1 + SUGGESTIVE + TIER2_INV_AUX + TIER2_OTHER

COL_INV   = "#C62828"   # red bold
COL_AUX   = "#E65100"   # orange
COL_OTHER = "#616161"   # grey
COL_SUG   = "#9C27B0"   # purple (suggestive, non-tier-1)

def label_color(gene):
    if gene in INV_TIER1: return COL_INV, 'bold'
    if gene in SUGGESTIVE: return COL_SUG, 'bold'
    if gene in TIER2_INV_AUX: return COL_AUX, 'semibold'
    return COL_OTHER, 'normal'

COHORTS = [
    {"name":"Jäkel 2019\nBrain snRNA-seq",
     "info":"17,799 nuclei · GSE118257", "short":"Jakel_2019_brain",
     "path":FIG_DIR/"brain_Jakel2019/adata_jakel.h5ad",
     "ct":"celltype", "dx":"group", "ms":"MS", "hc":"HC",
     "genes":["ITGB2","LXN","CD79B","IKZF1","SH3BP4"]},
    {"name":"Kaufmann 2021\nPBMC 10x Chromium 3′",
     "info":"497,705 cells · GSE144744", "short":"Ramesh_2020_PBMC",
     "path":FIG_DIR/"blood_Ramesh2020_UMAP/adata_ramesh_umap.h5ad",
     "ct":"basictype", "dx":"disease_status", "ms":"MS", "hc":"HC",
     "genes":["CD79B","IKZF1","SH3BP4","HLA-E","ITGB2"]},
    {"name":"Beltrán 2019\nCSF + PBMC twins",
     "info":"2,029 cells · 5 twin pairs · GSE127969", "short":"Beltran_2019_CSFPBMC",
     "path":FIG_DIR/"blood_Beltran2019/adata_beltran.h5ad",
     "ct":"celltype", "dx":"group", "ms":"MS", "hc":"HC",
     "genes":["HLA-E","ITGB2","CD79B","IKZF1","LXN"]},
]

print("Loading h5ad files...")
ads = []
for c in COHORTS:
    ad = sc.read_h5ad(c["path"])
    obs = ad.obs.copy()
    obs[c["ct"]] = obs[c["ct"]].astype(str)
    obs[c["dx"]] = obs[c["dx"]].astype(str)
    ad.obs = obs
    ads.append(ad)

pb = pd.read_csv(PB, sep="\t")
print(f"Wilcoxon tests: {len(pb)}, sig BH-FDR<0.05: {(pb.fdr<0.05).sum()}")

# Build heatmap row order = INV_TIER1 then TIER2_INV_AUX then TIER2_OTHER
pb["co_ct"] = pb["cohort"] + "  ·  " + pb["cell_type"]
heat_d   = pb.pivot_table(index="gene", columns="co_ct", values="logFC", aggfunc="mean").clip(-2.5,2.5)
heat_fdr = pb.pivot_table(index="gene", columns="co_ct", values="fdr",  aggfunc="min")  # Benjamini–Hochberg FDR over all 629 candidate gene × cell-type tests (multiple-testing corrected)
heat_p   = pb.pivot_table(index="gene", columns="co_ct", values="pval", aggfunc="min")  # nominal tie-corrected Wilcoxon p (directional)
heat_d   = heat_d.loc[heat_d.index.isin(ALL_GENES)]
heat_fdr = heat_fdr.loc[heat_fdr.index.isin(ALL_GENES)]
heat_p   = heat_p.loc[heat_p.index.isin(ALL_GENES)]
_sig=(heat_p<0.05).sum(); _abs=heat_d.clip(-2.5,2.5).abs().sum()
ranked = pd.DataFrame({"s":_sig,"a":_abs}).sort_values(["s","a"],ascending=False).index
# guarantee each gene's strongest nominally-significant cell-type column is displayed (every gene incl. MX1, e.g. LXN·Oligo3)
must=set()
for _g in heat_p.index:
    _sr=heat_p.loc[_g]; _sr=_sr[_sr<0.05]
    if len(_sr): must.add(_sr.idxmin())
_keep=set(list(ranked[:8])+list(must))   # 8 busiest columns ∪ each gene's strongest nominal column → ~15 cols (keeps figure readable, every gene shown)
col_order=[c for c in ranked if c in _keep][:28]
heat_d = heat_d[col_order]; heat_fdr = heat_fdr[col_order]; heat_p = heat_p[col_order]
row_order = [g for g in ALL_GENES if g in heat_d.index]
heat_d = heat_d.loc[row_order]
heat_fdr = heat_fdr.loc[row_order]
heat_p = heat_p.loc[row_order]
print(f"Heatmap shape: {heat_d.shape}")

fig = plt.figure(figsize=(22, 23.5), dpi=200)
gs = gridspec.GridSpec(7, 6, figure=fig,
                        height_ratios=[0.80, 0.80, 0.05, 3.0, 0.62, 0.85, 0.85],
                        hspace=0.32, wspace=0.30)

PALETTE = (plt.cm.tab20.colors + plt.cm.tab20b.colors)[:28]

# ─── Panel A ───────────────────────────────────────────────────────────────
for col, (cohort, ad) in enumerate(zip(COHORTS, ads)):
    U = ad.obsm["X_umap"]; obs = ad.obs
    ax = fig.add_subplot(gs[0, col*2:col*2+2])
    cts = sorted(obs[cohort["ct"]].unique(),
                  key=lambda x: -(obs[cohort["ct"]]==x).sum())
    cmap_ct = {ct: PALETTE[i % len(PALETTE)] for i, ct in enumerate(cts)}
    # plot largest class FIRST so small populations render on top (not overplotted)
    for ct in cts:
        m = (obs[cohort["ct"]] == ct).values
        ax.scatter(U[m,0], U[m,1], s=2.2, c=[cmap_ct[ct]], alpha=0.75,
                    edgecolors='none', rasterized=True)
    ax.set_title(f"{cohort['name']}\n{cohort['info']}",
                  fontsize=16.0, fontweight='bold')
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    _CTLBL = {"t_cells":"T cells","monocytes":"Monocytes","nk_cells":"NK cells",
              "b_cells":"B cells","cdc":"cDC","pdc":"pDC","platelets":"Platelets",
              "plasma_cells":"Plasma cells","Microglia_Macrophages":"Microglia/Macro",
              "Endothelial_cells1":"Endothelial 1","Endothelial_cells2":"Endothelial 2"}
    def _clean(ct): return _CTLBL.get(ct, ct.replace("_"," "))
    handles = [mpatches.Patch(color=cmap_ct[ct],
                                label=f"{_clean(ct)} (n={(obs[cohort['ct']]==ct).sum():,})")
                for ct in cts[:8]]
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.005, 0.42),
                fontsize=7.6, frameon=False, ncol=1, handlelength=1.0,
                handletextpad=0.4, labelspacing=0.28, borderaxespad=0.1)

    ax2 = fig.add_subplot(gs[1, col*2:col*2+2])
    dx = obs[cohort["dx"]].values
    ms_mask = np.isin(dx, [cohort["ms"], f"{cohort['ms']}_PBMCs", "MS"])
    hc_mask = np.isin(dx, [cohort["hc"], f"{cohort['hc']}_PBMCs", "HD",
                              "HD_PBMCs", "Ctrl"])
    if hc_mask.sum() == 0: hc_mask = ~ms_mask
    ax2.scatter(U[:,0], U[:,1], s=1.8, c="#E0E0E0", alpha=0.4, edgecolors='none',
                 rasterized=True)
    ax2.scatter(U[hc_mask,0], U[hc_mask,1], s=2.5, c="#1976D2", alpha=0.55,
                 edgecolors='none', rasterized=True,
                 label=f"HC (n={int(hc_mask.sum()):,})")
    ax2.scatter(U[ms_mask,0], U[ms_mask,1], s=2.5, c="#D32F2F", alpha=0.55,
                 edgecolors='none', rasterized=True,
                 label=f"MS (n={int(ms_mask.sum()):,})")
    ax2.set_title("MS vs HC overlay", fontsize=16.0, fontweight='bold')
    ax2.set_xticks([]); ax2.set_yticks([])
    for sp in ax2.spines.values(): sp.set_visible(False)
    ax2.legend(loc='lower right', fontsize=12.8, frameon=True,
                facecolor='white', edgecolor='#888')

fig.text(0.015, 0.985, "A", fontsize=38.4, fontweight='bold', va='top', ha='left')
plt.subplots_adjust(top=0.965, bottom=0.02)

# ─── Panel B: Tier-stratified log₂ fold-change heatmap ─────────────────────
# (Panel B label "B" is drawn on the heatmap axes; descriptive text moved to the caption)

ax_heat = fig.add_subplot(gs[3, :])
cmap = LinearSegmentedColormap.from_list("RdBu_r_clip",
        ["#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0",
         "#FFFFFF",
         "#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"])
plot_mat = np.clip(heat_d.values, -2.5, 2.5)
from matplotlib.colors import Normalize
norm = Normalize(vmin=-2.5, vmax=2.5)
rgba = cmap(norm(np.nan_to_num(plot_mat, nan=0.0)))
for i in range(heat_d.shape[0]):
    for j in range(heat_d.shape[1]):
        if pd.isna(heat_d.values[i, j]):
            rgba[i, j, :3] = 1.0                                    # WHITE = gene not testable here (no data / <20 cells / absent)
        # every TESTED cell keeps its log₂FC colourmap colour — significant AND non-significant (no grey-out)
im = ax_heat.imshow(rgba, aspect='auto')
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])

# Every tested cell is coloured by log₂FC and shows its fold-change value; stars mark BH-FDR significance (FDR<0.05)
for i in range(heat_d.shape[0]):
    for j in range(heat_d.shape[1]):
        v = heat_d.values[i, j]; p = heat_p.values[i, j]; fdr = heat_fdr.values[i, j]
        if pd.isna(v):
            continue                                               # not testable → leave white/blank
        tc = 'white' if abs(plot_mat[i,j]) > 1.5 else '#212121'
        sig = pd.notna(fdr) and fdr < 0.05                         # BH-FDR significance (consistent with Methods/Results/caption)
        if sig:                                                    # FDR-significant: fold-change VALUE + graded star (both)
            ax_heat.text(j, i-0.18, f"{v:+.2f}", ha='center', va='center',
                          fontsize=12.8, fontweight='bold', color=tc)
            stars = "***" if fdr < 0.001 else ("**" if fdr < 0.01 else "*")  # graded star by BH-FDR
            ax_heat.text(j, i+0.26, stars, ha='center', va='center',
                          fontsize=18.0, fontweight='bold', color=tc)
        else:                                                      # tested but NS: show fold-change value (colour-mapped), no star
            ax_heat.text(j, i, f"{v:+.2f}", ha='center', va='center',
                          fontsize=11.0, fontweight='normal', color=tc, alpha=0.92)

# Horizontal separator lines between gene tiers
n_t1  = sum(1 for g in row_order if g in INV_TIER1)
n_sug = sum(1 for g in row_order if g in SUGGESTIVE)
n_aux = sum(1 for g in row_order if g in TIER2_INV_AUX)
N_rows = len(row_order)
# inverse-concordant Tier-1 block
ax_heat.axhline(n_t1 - 0.5, color='#000', linewidth=2.0, alpha=0.7)
ax_heat.text(-2.3, (n_t1-1)/2, f"Inverse-concordant\nTier-1 ({n_t1} genes)",
              ha='right', va='center', fontsize=14.4, fontweight='bold',
              color=COL_INV, transform=ax_heat.transData)
# suggestive (HLA-E) row — demoted, non-tier-1
if n_sug:
    ax_heat.axhline(n_t1 + n_sug - 0.5, color='#000', linewidth=1.1, alpha=0.6, linestyle=':')
    ax_heat.text(-2.3, n_t1 + (n_sug-1)/2, "suggestive\n(HLA-E)",
                  ha='right', va='center', fontsize=13.6, fontweight='bold',
                  color=COL_SUG, transform=ax_heat.transData)
# Tier-2 auxiliary inverse-concordant block
if n_aux:
    ax_heat.axhline(n_t1 + n_sug + n_aux - 0.5, color='#000', linewidth=1.5, alpha=0.5, linestyle='--')
    ax_heat.text(-2.3, n_t1 + n_sug + n_aux/2 - 0.5, "Tier-2 auxiliary\ninverse-concordant",
                  ha='right', va='center', fontsize=14.4, fontweight='semibold',
                  color=COL_AUX, transform=ax_heat.transData)
# non-concordant block
if n_t1 + n_sug + n_aux < N_rows:
    ax_heat.text(-2.3, n_t1 + n_sug + n_aux + (N_rows - n_t1 - n_sug - n_aux)/2 - 0.5,
                  "non-concordant\naux", ha='right', va='center',
                  fontsize=14.4, color=COL_OTHER, transform=ax_heat.transData)

ax_heat.set_xticks(range(heat_d.shape[1]))
short_cols=[c.replace("Jakel_2019_brain","Jäkel").replace("Ramesh_2020_PBMC","Kaufmann").replace("Beltran_2019_CSFPBMC","Beltrán").replace("Microglia_Macrophages","Microglia/Mac").replace("Endothelial_cells","Endo").replace("  ·  "," · ") for c in heat_d.columns]
ax_heat.set_xticklabels(short_cols, rotation=42, ha='right', fontsize=12)
ax_heat.set_yticks(range(heat_d.shape[0]))
yt_labels = []
for g in heat_d.index:
    yt_labels.append(g)
ax_heat.set_yticklabels(yt_labels, fontsize=16.0, fontstyle='italic')
# Recolour y-tick labels per tier
for tick, g in zip(ax_heat.get_yticklabels(), heat_d.index):
    col, fw = label_color(g)
    tick.set_color(col)
    tick.set_fontweight(fw)
cbar = plt.colorbar(sm, ax=ax_heat, fraction=0.018, pad=0.012)
cbar.set_label("log₂ fold change (MS vs HC, capped ±2.5)", fontsize=14.4)
cbar.ax.tick_params(labelsize=12.8)

# (long explanatory note removed — all detail lives in the Figure 5 caption)
ax_heat.text(-0.045, 1.025, "B", transform=ax_heat.transAxes, fontsize=38.4, fontweight='bold', ha='right', va='bottom')

# ─── Panel C: gene expression UMAPs (focus on Inverse-concordant Tier-1) ──────────────────
inner = GridSpecFromSubplotSpec(3, 5, subplot_spec=gs[5:7, :],
                                  hspace=0.40, wspace=0.18)
for row, (cohort, ad) in enumerate(zip(COHORTS, ads)):
    U = ad.obsm["X_umap"]
    X = ad.X
    if issparse(X): X = X.toarray()
    for c_idx, gene in enumerate(cohort["genes"]):
        ax_g = fig.add_subplot(inner[row, c_idx])
        if row == 0 and c_idx == 0:
            ax_g.text(-0.32, 1.28, "C", transform=ax_g.transAxes, fontsize=38.4, fontweight='bold', ha='left', va='top')
        if gene not in ad.var_names:
            ax_g.text(0.5, 0.5, f"{gene}\n(n.q.)", ha='center', va='center',
                       transform=ax_g.transAxes, fontsize=14.4, fontstyle='italic',
                       color='#888')
            ax_g.set_xticks([]); ax_g.set_yticks([])
            for sp in ax_g.spines.values(): sp.set_visible(False)
            continue
        gi = list(ad.var_names).index(gene)
        expr = X[:, gi]
        vmax = np.quantile(expr[expr > 0], 0.99) if (expr > 0).sum() else 1.0
        order = np.argsort(expr)
        sc_plot = ax_g.scatter(U[order,0], U[order,1], c=expr[order],
                                 s=1.8, cmap='viridis', vmin=0, vmax=vmax,
                                 alpha=0.85, edgecolors='none', rasterized=True)
        pct = 100 * (expr > 0).sum() / len(expr)
        # Title color marks gene tier
        col, fw = label_color(gene)
        ax_g.set_title(f"${gene}$  ({pct:.0f}%>0)", fontsize=14.4,
                        fontweight='bold', color=col)
        if c_idx == 0:
            ax_g.set_ylabel(cohort["name"].split('\n')[0], fontsize=14.4,
                             fontweight='bold', rotation=90)
        ax_g.set_xticks([]); ax_g.set_yticks([])
        for sp in ax_g.spines.values(): sp.set_visible(False)
        cb = plt.colorbar(sc_plot, ax=ax_g, shrink=0.65, pad=0.02)
        cb.ax.tick_params(labelsize=9.6)

# (suptitle removed — descriptive panel text moved to the figure caption; bare A/B/C letters only)

plt.savefig(OUT, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\n✓ saved → {OUT}")
from PIL import Image
img = Image.open(OUT)
print(f"  → {img.size[0]} × {img.size[1]} ({img.size[0]/img.size[1]:.2f}:1)")
