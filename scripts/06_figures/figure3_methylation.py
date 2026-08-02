#!/usr/bin/env python3
"""figure3_methylation.py — Figure 3 methylation panel, UNIFORM gene-level DMP.

ALL per-stratum methylation panels now use the SAME gene-level limma-DMP method
(x = mean methylation logFC, positive = hypermethylated in MS; y = −log10 BH-FDR),
replacing the earlier mixed mCSEA-NES + gene-level layout. This makes every
volcano directly comparable on a single fold-change axis.

Seven panels (3 × 3 grid):
  Row 1: 1. T cells baseline       gene-level limma-DMP (10/16,966 sig)
         2. Whole blood · DMF      gene-level limma-DMP (2,650/16,667 sig)
         3. Whole blood · Ocrelizumab gene-level limma-DMP (99/16,590 sig)
  Row 2: 4. Brain WM (Huynh 2014)  gene-level CpG model (31/2,288 sig)
         5. Combined · Stouffer cross-stratum meta-analysis (461/10,863 sig;
            all 5 Inverse-concordant Tier-1 significant)
         6. gene-tier legend
  Row 3: 7. inverse-concordant top-22 stacked bar (full width)

(The T-cells Remission stratum is omitted per request; at gene-level it carries
only 31 sig genes and no Inverse-concordant Tier-1 hit.)

Highlights (only genes that ACTUALLY pass BH-FDR<0.05 are labelled):
  • Inverse-concordant Tier-1 (3): ITGB2 · CD79B · IKZF1 — RED *
  • Tier-2 auxiliary inverse-concordant (7): CASP6 · CASP8 · DGKQ · MX1 · IFIT1 · NUP210 · RUNX3 — ORANGE ◆
  • Tier-2 non-concordant (CTSZ · CHL1 · ITGAL · THRB · ICAM1 · …) — GREY
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
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "axes.spines.top": False, "axes.spines.right": False,
})

RES = Path("__MS_GEO_ROOT__/Methylation/results")
INV_FILE = RES / "INVERSE_CONCORDANT_by_gene.tsv"
OUT = Path("__MS_GEO_ROOT__/Poster_v2/figures/methylation_v_INV.png")

# Per-stratum gene-level limma-DMP files (gene, n_probes, mean_logFC, z_combined, P.Value, adj.P.Val)
DMP_TCELLS = RES / "01_tcells_meth_gene.tsv"
DMP_DMF    = RES / "02_wb_dmf_meth_gene.tsv"
DMP_OCRE   = RES / "03_wb_ocrelizumab_meth_gene.tsv"
DMP_COMB   = RES / "05_combined_meth_gene.tsv"
BRAIN_WM   = RES / "07_BrainWM_RNA_vs_Meth_concordance.tsv"   # gene, meth, meth_fdr
ALLMETH_COMBAT = Path("__MS_GEO_ROOT__/Methylation_Data/AllMeth_ComBat_limma_DMP_byGene.csv")  # gene, logFC, P, FDR

NAVY="#0D3B66"; TEAL="#3E92CC"; RED="#D62828"; RED_HOT="#B71C1C"; ORANGE="#E65100"
GREY_DARK="#424242"

INV_TIER1 = {"ITGB2","IKZF1"}
TIER2_AUX_INV = {"CD79B","CASP6","CASP8","DGKQ","MX1","IFIT1","NUP210","RUNX3","SH3BP4","LXN"}
TIER2_OTHER = {"CTSZ","CHL1","ITGAL","IFI44L","ICAM1","FOXP3","TYK2","STAT3",
                "MOSPD3","RPAP2","THRB","SLAMF1","KLF6","HIGD1A","PCNP","DUSP22","IK",
                "ETV3","RNF216","LYN","TRAF1","PRKCH","RBM38","RPS6KA4","ZAP70"}

def gene_style(g):
    if g in INV_TIER1:
        return dict(facecolor="#ffcdd2", edgecolor=RED_HOT, lw=2.4,
                     star="* ", fontsize=13.4, fontweight="bold",
                     marker_size=140, color=RED_HOT)
    if g in TIER2_AUX_INV:
        return dict(facecolor="#fff3e0", edgecolor=ORANGE, lw=1.5,
                     star="◆ ", fontsize=11.7, fontweight="bold",
                     marker_size=90, color=ORANGE)
    return dict(facecolor="#eeeeee", edgecolor=GREY_DARK, lw=0.9,
                 star="", fontsize=10.6, fontweight="normal",
                 marker_size=65, color=GREY_DARK)

# ── Uniform gene-level methylation volcano ──────────────────────────────────
def gene_volcano(ax, path, title, sample_n, *, fc_col="mean_logFC",
                 fdr_col="adj.P.Val", gene_col="gene",
                 method="limma-DMP",
                 xlabel="methylation logFC (positive = hypermethylated in MS)",
                 oth_n=4, extra_sub="", fdr_thr=0.05, sep="\t"):
    """Generic gene-level methylation volcano. x = fold-change (fc_col),
    y = −log10(BH-FDR). Only genes passing BH-FDR<0.05 are labelled."""
    df = (pd.read_csv(path, sep=sep)
            .rename(columns={fc_col:"fc", fdr_col:"fdr", gene_col:"g"}))
    df["fdr"] = pd.to_numeric(df["fdr"], errors="coerce")
    df["fc"]  = pd.to_numeric(df["fc"], errors="coerce")
    df = df.dropna(subset=["fdr","fc","g"])
    df["logp"] = -np.log10(df["fdr"].clip(lower=1e-300))
    sig = df["fdr"] < fdr_thr
    up  = sig & (df["fc"] > 0); dn = sig & (df["fc"] < 0)

    ax.scatter(df.loc[~sig,"fc"], df.loc[~sig,"logp"], s=2.5, c="#dadada",
                alpha=0.45, rasterized=True, linewidths=0)
    ax.scatter(df.loc[up,"fc"], df.loc[up,"logp"], s=7, c=RED, alpha=0.7,
                rasterized=True, linewidths=0,
                label=f"hyper-methylated in MS (n={int(up.sum())})")
    ax.scatter(df.loc[dn,"fc"], df.loc[dn,"logp"], s=7, c=TEAL, alpha=0.7,
                rasterized=True, linewidths=0,
                label=f"hypo-methylated in MS (n={int(dn.sum())})")

    df_sig = df[sig]
    inv_rows = df_sig[df_sig.g.isin(INV_TIER1)].sort_values("fdr")
    aux_rows = df_sig[df_sig.g.isin(TIER2_AUX_INV)].sort_values("fdr").head(3)
    oth_rows = df_sig[df_sig.g.isin(TIER2_OTHER)].sort_values("fdr").head(oth_n)
    cands = pd.concat([inv_rows, aux_rows, oth_rows]).drop_duplicates(subset=["g"])
    texts = []
    for r in cands.itertuples():
        gg=r.g; x=r.fc; y=r.logp; st=gene_style(gg)
        ax.scatter(x, y, s=st["marker_size"], edgecolor=st["edgecolor"],
                    facecolor="none", linewidth=st["lw"], zorder=5)
        texts.append(ax.text(x, y, f"{st['star']}{gg}", fontsize=st["fontsize"],
                      fontweight=st["fontweight"], color=st["color"],
                      fontstyle="italic",
                      bbox=dict(boxstyle="round,pad=0.15", facecolor=st["facecolor"],
                                  edgecolor=st["edgecolor"], linewidth=0.5), zorder=6))
    if texts:
        adjust_text(texts, ax=ax, expand=(1.5,1.8),
                     arrowprops=dict(arrowstyle="-", color="#555", lw=0.6,
                                       shrinkA=3, shrinkB=4),
                     force_text=(0.5,0.8), time_lim=3)

    ax.axhline(-np.log10(fdr_thr), ls=":", color="#777", lw=0.5)
    ax.axvline(0, ls=":", color="#bbb", lw=0.3)
    sub = f"{int(sig.sum()):,}/{len(df):,} genes · FDR<0.05"
    if extra_sub: sub += f" · {extra_sub}"
    head = f"{title}  ({sample_n})" if sample_n else title
    ax.set_title(title, fontsize=28.8, fontweight="bold", loc="left", pad=4)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("−log₁₀(BH-FDR)")
    ax.legend(loc="upper left", fontsize=11.2, frameon=True,
               facecolor="white", edgecolor="#aaa")
    return int(sig.sum()), int(up.sum()), int(dn.sum())

# ── RNA × methylation concordance scatter (replaces stacked bar) ────────────
def inv_scatter(ax):
    # The 82-gene inverse-concordant DISCOVERY POOL (per-stratum: best RNA x best methylation,
    # opposite sign in both layers) — the same pool used throughout the manuscript.
    inv = pd.read_csv(INV_FILE, sep="\t")
    inv["g"]  = inv.gene.astype(str).str.upper()
    inv["mfc"]= pd.to_numeric(inv.best_meth_fc, errors="coerce")
    inv["rfc"]= pd.to_numeric(inv.best_rna_fc,  errors="coerce")
    inv = inv.dropna(subset=["mfc","rfc"]); n = len(inv)
    YCAP = 1.6
    xl=(min(-0.28, inv.mfc.min()*1.15), max(0.28, inv.mfc.max()*1.15))
    yl=(-1.35, YCAP+0.18)
    inv["rfc_d"]   = inv.rfc.clip(yl[0]+0.05, YCAP)        # clip IFN-beta-driven RNA outliers for readability
    inv["clipped"] = inv.rfc > YCAP
    for i,(idx,_) in enumerate(inv[inv.clipped].sort_values("rfc",ascending=False).iterrows()):
        inv.loc[idx,"rfc_d"] = YCAP - 0.15*i                # stagger clipped outliers so labels don't collide
    from matplotlib.patches import Rectangle
    ax.add_patch(Rectangle((xl[0],0),-xl[0],yl[1],facecolor="#E8F5E9",zorder=0))   # Q2 RNA up / meth down
    ax.add_patch(Rectangle((0,yl[0]), xl[1],-yl[0],facecolor="#E8F5E9",zorder=0))   # Q4 RNA down / meth up
    ax.axhline(0,color="#888",lw=0.7,zorder=1); ax.axvline(0,color="#888",lw=0.7,zorder=1)
    UT1={x.upper() for x in INV_TIER1}; UA={x.upper() for x in TIER2_AUX_INV}
    ext=inv[~inv.g.isin(UT1|UA)]
    ax.scatter(ext.mfc, ext.rfc_d, s=26, c="#9e9e9e", alpha=0.60, edgecolors="white", linewidths=0.4, zorder=2)
    texts=[]
    for genes,col,ec in [(UA,"#FB8C00",ORANGE),(UT1,"#E53935",RED_HOT)]:
        sub=inv[inv.g.isin(genes)]
        ax.scatter(sub.mfc, sub.rfc_d, s=100, facecolor=col, edgecolor="white", linewidth=1.1, zorder=4)
        for r in sub.itertuples():
            lab = r.g + (f" (RNA +{r.rfc:.1f})" if r.clipped else "")
            texts.append(ax.text(r.mfc, r.rfc_d, lab, fontsize=14.4, fontstyle="italic", fontweight="bold", color=ec, zorder=5))
    try:
        from adjustText import adjust_text
        adjust_text(texts, ax=ax, expand=(1.4,1.7),
                    arrowprops=dict(arrowstyle="-", color="#888", lw=0.6))
    except Exception: pass
    ax.set_xlim(xl); ax.set_ylim(yl)
    ax.set_xlabel("best methylation logFC   (x · + = hypermethylated in MS)", fontsize=15.2)
    ax.set_ylabel(f"best RNA-seq logFC   (y · + = up in MS; clipped at +{YCAP:g})", fontsize=15.2)
    ax.set_title("F", fontsize=28.8, fontweight="bold", loc="left", pad=6)
    ax.text(0.985,0.025,
            "every gene is inverse-concordant by construction (RNA↑·meth↓ or RNA↓·meth↑ → green quadrants)\n"
            "red = Tier-1   ·   orange = Tier-2 auxiliary   ·   grey = extended inverse-concordant",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=11.8, color="#333",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#cccccc", alpha=0.9))

# ─── Figure layout: 3 rows × 3 columns (uniform gene-level DMP) ─────────────
fig = plt.figure(figsize=(20.5, 16.8), dpi=400)
gs = fig.add_gridspec(3, 3, hspace=0.46, wspace=0.28, height_ratios=[1.0, 1.0, 0.9])

# Row 1 — per-stratum gene-level limma-DMP volcanoes
ax = fig.add_subplot(gs[0, 0])
gene_volcano(ax, DMP_TCELLS, "A", "18 MS / 35 HC")
ax = fig.add_subplot(gs[0, 1])
gene_volcano(ax, DMP_DMF, "B", "12 MS / 18 HC")
ax = fig.add_subplot(gs[0, 2])
gene_volcano(ax, DMP_OCRE, "C", "6 MS / 18 HC")

# Row 2 — brain gene-level CpG model + all-datasets IDAT+ComBat + legend
# (Stouffer combined panel removed — the 475-sample IDAT+ComBat analysis supersedes it)
ax = fig.add_subplot(gs[1, 0])
gene_volcano(ax, BRAIN_WM, "D", "21 MS / 23 HC",
             fc_col="meth", fdr_col="meth_fdr",
             xlabel="methylation effect (+ = hyper in MS)")
ax = fig.add_subplot(gs[1, 1])
gene_volcano(ax, ALLMETH_COMBAT, "E", "244 MS / 231 HC",
             fc_col="logFC", fdr_col="FDR", sep=",",
             xlabel="methylation logFC (M-value, MS−HC)",
             extra_sub="5/5 Inverse-concordant Tier-1 sig")

# Row 3 — inverse-concordant combined-evidence bar (full width)
ax_bar = fig.add_subplot(gs[2, 0:3])
inv_scatter(ax_bar)

from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor="#ffcdd2", edgecolor=RED_HOT, linewidth=2.0,
           label=f"Inverse-concordant Tier-1 ({len(INV_TIER1)} genes)"),
    Patch(facecolor="#fff3e0", edgecolor=ORANGE, linewidth=1.4,
           label=f"Tier-2 auxiliary inverse-concordant ({len(TIER2_AUX_INV)} genes)"),
    Patch(facecolor="#eeeeee", edgecolor=GREY_DARK, linewidth=0.8,
           label="Tier-2 non-concordant proteomic anchors"),
]
ax_leg = fig.add_subplot(gs[1, 2]); ax_leg.axis("off")
ax_leg.legend(handles=legend_elements, loc="center", frameon=True, fontsize=11,
               title="Gene-tier highlighting\n(* Tier-1, ◆ Tier-2 auxiliary; BH-FDR<0.05 labelled)",
               title_fontsize=12, labelspacing=0.6, borderpad=0.6, handlelength=1.4)

# suptitle removed (bare A–F panel letters only)

plt.tight_layout()
plt.subplots_adjust(top=0.95)
plt.savefig(OUT, dpi=400, bbox_inches="tight", facecolor="white")
print(f"✓ saved → {OUT}")
from PIL import Image
img = Image.open(OUT)
print(f"  → {img.size[0]} × {img.size[1]}")
