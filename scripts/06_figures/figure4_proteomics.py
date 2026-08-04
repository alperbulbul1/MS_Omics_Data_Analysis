#!/usr/bin/env python3
"""figure4_proteomics.py — Figure 5 (proteomic validation), CORRECTED.

KEY CORRECTION: the "timsTOF" dataset is CSF (Bader & Mann 2024, 2nd platform of
the Astral CSF cohort; n=1,536 MS / 2,363 HC), NOT a separate "Wang & Julien
brain" cohort. The genuine brain proteomics is Wang & Julien 2025 (region-resolved
CTX / NAWM / WML, n=8 per group).

Layout (3 rows):
  Row 1 — three large-cohort volcanoes: CSF Astral · CSF timsTOF · Blood UKB
  Row 2 — four brain-region volcanoes (Wang & Julien 2025): CTX · NAWM · WML-vs-WM · WML-vs-NAWM
  Row 3 — directional-consistency heatmap across all 7 proteomic compartments

All panels re-analysed from primary or summary-statistic data.
"""
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from adjustText import adjust_text
from pathlib import Path
import warnings; warnings.filterwarnings('ignore')

PROJ = Path("__MS_GEO_ROOT__")
OUT  = PROJ / "Poster_v2" / "figures" / "proteomics_v_INV.png"

INV_TIER1  = ["ITGB2","IKZF1"]
# Tier-2 auxiliary inverse-concordant genes shown in the proteomic panels. Only LXN is carried
# here: it is the auxiliary gene whose (nominal) proteomic support the text cites, and it is
# quantified in the brain contrasts and in UK Biobank plasma. It falls below the completeness
# threshold on both CSF instruments, so it is legitimately absent from panels A and B.
# Coloured ORANGE, the Tier-2 auxiliary inverse-concordant colour used in Figures 2, 3 and 6.
# NOT purple: purple denotes the Tier-2 non-concordant PROTEOMIC ANCHORS (CTSZ/CHL1/ICAM1),
# a different subcategory. NOT teal: teal is reserved for the suggestive candidate HLA-E.
TIER2_AUX  = ["LXN"]
# FOXP3 was removed from this group on revision: it is not quantified in ANY of the seven
# proteomic compartments (both CSF instruments, all four brain-region contrasts, UK Biobank-PPP),
# so it could not be a "strong proteomic candidate", which is what defines this group. It is
# retained in the STRING display as a canonical MS immune context gene, which is the role it
# actually plays (the IKZF1-RUNX3/FOXP3-STAT1-STAT3 axis).
TIER2_PROT = ["CTSZ","CHL1","ICAM1","ITGAL"]
RED_HOT="#B71C1C"; PURPLE="#6A1B9A"; ORANGE="#E65100"   # ORANGE = Tier-2 auxiliary inverse-concordant, as in Figures 2, 3 and 6

def gene_style(g):
    if g in INV_TIER1:
        return dict(color=RED_HOT, fc="#ffcdd2", marker_c="#D32F2F", fs=9, fw="bold", ms=170)
    if g in TIER2_AUX:
        return dict(color=ORANGE, fc="#FFE0B2", marker_c="#F57C00", fs=8.7, fw="bold", ms=150)
    if g in TIER2_PROT:
        return dict(color=PURPLE, fc="#E1BEE7", marker_c="#8E24AA", fs=8.4, fw="bold", ms=135)
    return None

# ── Load CSF (Astral + timsTOF, both Bader/Mann) ───────────────────────────
def load_dedup(path, rename):
    df = pd.read_csv(path, sep="\t").rename(columns=rename).dropna(subset=["log2FC","FDR","gene"])
    df = df.sort_values("FDR").drop_duplicates("gene", keep="first").reset_index(drop=True)
    df["neg10"] = -np.log10(df.FDR.clip(1e-300))
    return df
astral = load_dedup(PROJ/"Proteomics/processed/RDEP_CC/Astral_RDEP.tsv",
                     {"log2FC_MSvsCtrl":"log2FC","Genes":"gene"})
ctimstof = load_dedup(PROJ/"Proteomics/processed/RDEP_CC/timsTOF_RDEP.tsv", {"Genes":"gene"})

# ── Blood UKB (Jacobs 2024 Olink) ──────────────────────────────────────────
ukb = pd.read_csv(PROJ/"Proteomics/blood_raw/Jacobs2024_UKB_primary_DE.tsv", sep="\t")
ukb = ukb.dropna(subset=["gene","beta","fdr"])
ukb["gene"]=ukb.gene.astype(str).str.strip()
ukb = ukb.sort_values("fdr").drop_duplicates("gene",keep="first").reset_index(drop=True)
ukb["neg10"]=-np.log10(ukb.fdr.clip(1e-300))

# ── Brain Wang & Julien 2025 (4 region comparisons, n=8/group) ─────────────────
MAG = {
    "Brain CTX (MS vs ctrl)":      "Magliozzi2026_MS_CTX_vs_ODC_CTX.tsv",
    "Brain NAWM (MS vs ctrl)":     "Magliozzi2026_MS_NAWM_vs_ODC_WM.tsv",
    "Brain WML vs ctrl-WM":        "Magliozzi2026_MS_WML_vs_ODC_WM.tsv",
    "Brain WML vs MS-NAWM":        "Magliozzi2026_MS_WML_vs_MS_NAWM.tsv",
}
mags = {}
for k,v in MAG.items():
    df = pd.read_csv(PROJ/f"Proteomics/processed/RDEP_CC/{v}", sep="\t")
    df = df.rename(columns={"Gene":"gene"}).dropna(subset=["log2FC","p","gene"])
    df = df.sort_values("p").drop_duplicates("gene",keep="first").reset_index(drop=True)
    df["neg10p"] = -np.log10(df.p.clip(1e-300))
    mags[k]=df

# ── T-lineage (Unified) ────────────────────────────────────────────────────
ua = pd.read_csv(PROJ/"Methylation/results/Unified_All_Assays_Long.tsv", sep="\t")
prot_u = ua[ua.layer=="Proteomics"].copy()

print(f"CSF Astral: {len(astral)} genes, {int((astral.FDR<0.05).sum())} sig (978/306)")
print(f"CSF timsTOF: {len(ctimstof)} genes, {int((ctimstof.FDR<0.05).sum())} sig (1536/2363)")
print(f"Blood UKB: {len(ukb)} genes, {int((ukb.fdr<0.05).sum())} sig (407/39979)")
for k,df in mags.items(): print(f"  {k}: {len(df)} genes, {int((df.FDR<0.05).sum())} FDR-sig, {int((df.p<0.05).sum())} nominal-p<0.05 (n=8/8)")

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':17.6})
fig = plt.figure(figsize=(18, 18), dpi=190)
gs = gridspec.GridSpec(3, 12, figure=fig, height_ratios=[1.0, 0.9, 1.25],
                        hspace=0.55, wspace=2.2)

# ── generic volcano ────────────────────────────────────────────────────────
def volcano(ax, df, xcol, ycol, sigmask, title, subtitle, xlabel,
             extra_context=None, use_nominal=False):
    up = sigmask & (df[xcol]>0); dn = sigmask & (df[xcol]<0); ns=~(up|dn)
    ax.scatter(df.loc[ns,xcol], df.loc[ns,ycol], s=5, c="#C8C8C8", alpha=0.4,
                edgecolors='none', rasterized=True)
    ax.scatter(df.loc[up,xcol], df.loc[up,ycol], s=8, c="#D32F2F", alpha=0.6,
                edgecolors='none', rasterized=True, label=f"MS↑ (n={int(up.sum())})")
    ax.scatter(df.loc[dn,xcol], df.loc[dn,ycol], s=8, c="#1976D2", alpha=0.6,
                edgecolors='none', rasterized=True, label=f"MS↓ (n={int(dn.sum())})")
    texts=[]
    for tier in (INV_TIER1, TIER2_AUX, TIER2_PROT):
        for g in tier:
            r=df[df.gene.str.upper()==g.upper()]
            if not len(r): continue
            r=r.iloc[0]; st=gene_style(g)
            ax.scatter(r[xcol], r[ycol], s=st["ms"], c=st["marker_c"],
                        marker="*" if g in INV_TIER1 else ("o" if g in TIER2_AUX else "D"),
                        edgecolors='#000', linewidth=0.6, zorder=10)
            texts.append(ax.text(r[xcol], r[ycol], g, fontsize=st["fs"],
                          fontweight=st["fw"], fontstyle='italic', color=st["color"],
                          bbox=dict(boxstyle='round,pad=0.18', fc=st["fc"],
                                      ec=st["color"], lw=0.5, alpha=0.9), zorder=11))
    if extra_context:
        for g in extra_context:
            r=df[df.gene.str.upper()==g]
            if len(r):
                r=r.iloc[0]
                ax.scatter(r[xcol], r[ycol], s=45, facecolor="none",
                            edgecolor="#00897B", linewidth=1.1, zorder=9)
                texts.append(ax.text(r[xcol], r[ycol], g, fontsize=11.2,
                              color="#00695C", fontstyle='italic', zorder=11))
    if texts:
        adjust_text(texts, ax=ax, expand=(1.4,1.7),
                     arrowprops=dict(arrowstyle="-", color="#777", lw=0.5, shrinkA=3, shrinkB=4),
                     force_text=(0.5,0.8), time_lim=2)
    thr = -np.log10(0.05)
    ax.axhline(thr, color="#666", linestyle=":", lw=0.7)
    ax.axvline(0, color="#666", lw=0.4)
    ax.set_title(title, fontsize=28.8, fontweight="bold", loc="left", pad=4)
    ax.set_xlabel(xlabel, fontsize=15.2)
    ax.set_ylabel("−log₁₀(nominal p)" if use_nominal else "−log₁₀(BH-FDR)", fontsize=15.2)
    ax.tick_params(labelsize=13.6); ax.grid(alpha=0.1)
    ax.legend(loc='upper left', fontsize=12.0, frameon=True, facecolor='white', edgecolor='#aaa')
    for sp in ['top','right']: ax.spines[sp].set_visible(False)

# ── ROW 1: 3 large-cohort volcanoes ────────────────────────────────────────
ax = fig.add_subplot(gs[0, 0:4])
volcano(ax, astral, "log2FC", "neg10", astral.FDR<0.05,
        "A",
        f"978 MS / 306 HC · {int((astral.FDR<0.05).sum())}/{len(astral):,} at FDR<0.05",
        "log₂ fold-change (MS vs HC)")
ax = fig.add_subplot(gs[0, 4:8])
volcano(ax, ctimstof, "log2FC", "neg10", ctimstof.FDR<0.05,
        "B",
        f"1,536 MS / 2,363 HC · {int((ctimstof.FDR<0.05).sum())}/{len(ctimstof):,} at FDR<0.05",
        "log₂ fold-change (MS vs HC)")
ax = fig.add_subplot(gs[0, 8:12])
volcano(ax, ukb, "beta", "neg10", ukb.fdr<0.05,
        "C",
        f"407 MS / 39,979 HC · {int((ukb.fdr<0.05).sum())}/{len(ukb):,} of 2,911 at FDR<0.05",
        "association β (− = lower in MS plasma)",
        extra_context=["ITGAV","ITGAM","ITGA11","NEFL"])

# ── ROW 2: 4 Magliozzi brain-region volcanoes (n=8/group) ──────────────────
brain_axes_pos = [(1,0,3),(1,3,6),(1,6,9),(1,9,12)]
for _bk,((k,df),(r,c0,c1)) in enumerate(zip(mags.items(), brain_axes_pos)):
    ax = fig.add_subplot(gs[r, c0:c1])
    # n=8/group → no FDR<0.05; use nominal p for structure, honest subtitle
    volcano(ax, df, "log2FC", "neg10p", df.p<0.05,
            chr(68+_bk),
            f"n=8/8 · 0/{len(df):,} FDR<0.05 · {int((df.p<0.05).sum())} nominal p<0.05",
            "log₂FC (MS vs ctrl)", use_nominal=True)

# ── ROW 3: directional-consistency heatmap (7 compartments) ────────────────
ax = fig.add_subplot(gs[2, 1:11])
GENES_D = INV_TIER1 + TIER2_AUX + TIER2_PROT
COLS_D = [
    ("CSF Astral","astral"), ("CSF timsTOF","ctimstof"),
    ("Brain CTX","CTX"), ("Brain NAWM","NAWM"),
    ("Brain WML/WM","WMLWM"), ("Brain WML/NAWM","WMLNAWM"),
    ("Blood UKB (Jacobs)","ukb"),
]
mag_key = {"CTX":"Brain CTX (MS vs ctrl)","NAWM":"Brain NAWM (MS vs ctrl)",
           "WMLWM":"Brain WML vs ctrl-WM","WMLNAWM":"Brain WML vs MS-NAWM"}
val=np.full((len(GENES_D),len(COLS_D)),np.nan); sig=np.zeros_like(val,dtype=bool); sigp=np.full_like(val,np.nan)
for i,g in enumerate(GENES_D):
    for j,(lab,key) in enumerate(COLS_D):
        if key=="astral":
            r=astral[astral.gene.str.upper()==g.upper()]
            if len(r): val[i,j]=r.log2FC.iloc[0]; sigp[i,j]=r.FDR.iloc[0]; sig[i,j]=r.FDR.iloc[0]<0.05
        elif key=="ctimstof":
            r=ctimstof[ctimstof.gene.str.upper()==g.upper()]
            if len(r): val[i,j]=r.log2FC.iloc[0]; sigp[i,j]=r.FDR.iloc[0]; sig[i,j]=r.FDR.iloc[0]<0.05
        elif key=="ukb":
            r=ukb[ukb.gene.str.upper()==g.upper()]
            if len(r): val[i,j]=r.beta.iloc[0]; sigp[i,j]=r.fdr.iloc[0]; sig[i,j]=r.fdr.iloc[0]<0.05
        elif key in mag_key:
            df=mags[mag_key[key]]; r=df[df.gene.str.upper()==g.upper()]
            if len(r): val[i,j]=r.log2FC.iloc[0]; sigp[i,j]=r.p.iloc[0]; sig[i,j]=r.p.iloc[0]<0.05   # brain n=8: nominal p (matches volcano)
        else:
            r=prot_u[(prot_u.gene==g)&(prot_u.assay==key)]
            if len(r):
                rr=r.iloc[r.logFC.abs().argmax()]; val[i,j]=rr.logFC; sigp[i,j]=rr.fdr; sig[i,j]=rr.fdr<0.05
cmap=LinearSegmentedColormap.from_list("BuRd",
        ["#1565C0","#5E9BD4","#BBDEFB","#FFFFFF","#FFCDD2","#E57373","#C62828"])
norm=TwoSlopeNorm(vmin=-0.5,vcenter=0,vmax=0.5)
# Colour cells ONLY where significant (matching the volcano panels: FDR<0.05 for CSF/blood/T-lineage,
# nominal p<0.05 for the underpowered n=8 brain); non-significant measured cells fade toward grey.
rgba=cmap(norm(np.nan_to_num(val,nan=0.0)))
_GREY=np.array([0.93,0.93,0.93])
for _i in range(len(GENES_D)):
    for _j in range(len(COLS_D)):
        if np.isnan(val[_i,_j]): rgba[_i,_j,:3]=np.array([0.965,0.965,0.965])
        elif not sig[_i,_j]:     rgba[_i,_j,:3]=0.28*rgba[_i,_j,:3]+0.72*_GREY
rgba[...,3]=1.0
im=ax.imshow(rgba, aspect="auto")
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm); sm.set_array([])
ax.set_facecolor("#F5F5F5")
for i in range(len(GENES_D)):
    for j in range(len(COLS_D)):
        if not np.isnan(val[i,j]):
            if sig[i,j]:
                _pv=sigp[i,j]; _st="***" if _pv<0.001 else ("**" if _pv<0.01 else "*")
                ax.text(j,i,_st,ha='center',va='center',fontsize=15.0,
                         color='white' if abs(val[i,j])>0.3 else '#000', fontweight='bold')
        else:
            ax.text(j,i,"n.q.",ha='center',va='center',fontsize=9.6,color='#999')
ax.axhline(len(INV_TIER1)-0.5, color='#000', linewidth=2.0)
ax.axhline(len(INV_TIER1)+len(TIER2_AUX)-0.5, color='#000', linewidth=1.1)
ax.axvline(1.5, color='#000', lw=1.3); ax.axvline(5.5, color='#000', lw=1.3); ax.axvline(6.5, color='#000', lw=1.3)
ax.text(-1.5,(len(INV_TIER1)-1)/2,"Inverse-concordant Tier-1",ha='right',va='center',fontsize=14.4,fontweight='bold',color=RED_HOT)
ax.text(-1.5,len(INV_TIER1)+(len(TIER2_AUX)-1)/2,"Tier-2 aux\n(LXN)",ha='right',va='center',fontsize=12.4,fontweight='bold',color=ORANGE)
ax.text(-1.5,len(INV_TIER1)+len(TIER2_AUX)+(len(TIER2_PROT)-1)/2,"Tier-2\nprot anchor",ha='right',va='center',fontsize=13.6,fontweight='bold',color=PURPLE)
ax.set_xticks(range(len(COLS_D))); ax.set_xticklabels([c[0] for c in COLS_D],rotation=35,ha='right',fontsize=14.4)
ax.set_yticks(range(len(GENES_D))); ax.set_yticklabels(GENES_D,fontsize=15.2,fontstyle='italic',fontweight='bold')
for tick,g in zip(ax.get_yticklabels(),GENES_D): tick.set_color(RED_HOT if g in INV_TIER1 else (ORANGE if g in TIER2_AUX else PURPLE))
ax.set_title("H", fontsize=28.8, fontweight="bold", loc="left", pad=4)
cbar=plt.colorbar(sm,ax=ax,fraction=0.018,pad=0.01); cbar.set_label("log₂FC / β  (coloured only where significant)",fontsize=12.8); cbar.ax.tick_params(labelsize=12.0)

# suptitle removed (bare A–H panel letters only)

plt.savefig(OUT, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\n✓ saved → {OUT}")
from PIL import Image
img=Image.open(OUT); print(f"  → {img.size[0]} × {img.size[1]} ({img.size[0]/img.size[1]:.2f}:1)")
