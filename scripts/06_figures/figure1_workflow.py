#!/usr/bin/env python3
"""figure1_workflow.py — workflow with the current two-gene Tier-1."""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(7.5, 10.5), dpi=300)
ax.set_xlim(0, 10); ax.set_ylim(0, 14); ax.axis('off')

CI="#E3F2FD"; CP="#FFF3E0"; CT="#E8F5E9"; CIN="#FCE4EC"; CV="#F3E5F5"; CO="#FFEBEE"
ED="#37474F"; TX="#212121"

def box(x,y,w,h,t,fc,fs=9,fw='normal'):
    ax.add_patch(FancyBboxPatch((x,y),w,h, boxstyle="round,pad=0.04,rounding_size=0.10",
                                  facecolor=fc, edgecolor=ED, linewidth=1.4))
    ax.text(x+w/2, y+h/2, t, ha='center', va='center', fontsize=fs, fontweight=fw, color=TX, wrap=True)
def arr(x1,y1,x2,y2):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2), arrowstyle='->', mutation_scale=18,
                                    color=ED, linewidth=1.6))

ax.text(5, 13.55, "Integrated MS Multi-Omics Pipeline", ha='center', va='center',
         fontsize=14, fontweight='bold', color=TX)
ax.text(5, 13.15, "30 datasets · 4 omic layers · Tier-1 + Tier-2-aux candidates",
         ha='center', va='center', fontsize=10, color="#555", style='italic')

# Tier 1: Inputs
y=11.3
box(0.4,y,2.15,1.40,"BULK RNA\n15 datasets\n552 samples",CI,8.5,'bold')
box(2.7,y,2.15,1.40,"METHYLATION\n8 arrays · 475\n+ 1 WGBS series",CI,8.5,'bold')
box(5.0,y,2.15,1.40,"scRNA-seq\n3 cohorts\n81 donors",CI,8.5,'bold')
box(7.3,y,2.25,1.40,"PROTEOMICS\nCSF · brain\n· blood",CI,8.5,'bold')
for x in [1.475,3.775,6.075,8.425]: arr(x,y, x,y-0.5)

y=9.1
box(0.4,y,2.15,1.40,"Per-dataset norm\n(quantile / log₂-CPM)\nComBat batch\n+ limma-trend",CP,8.5)
box(2.7,y,2.15,1.40,"minfi β→M\nComBat\nlimma DMP\n+ mCSEA",CP,8.5)
box(5.0,y,2.15,1.40,"scanpy QC\nlog-norm · 3k HVG\nBBKNN\nLeiden",CP,8.5)
box(7.3,y,2.25,1.40,"log₂ LFQ/DIA\nlimma::removeBatchEffect\nDEP ≥50% filter\n+ raw-matrix rescue",CP,6.7)
for x in [1.475,3.775,6.075,8.425]: arr(x,y, x,y-0.4)

y=7.0
box(0.4,y,4.45,1.45,"Per-stratum DGE (BH-FDR<0.05)\nT cells: 2,177 DEGs · PBMC: 1,668\nIFN-β: 1,400 · B cells: 614\nWhole blood: 107 · Brain WM: 185",CT,8.5)
box(5.0,y,4.55,1.45,"Per-stratum mCSEA promoter GSEA\nDMF whole blood: 115 / 3,259 regions\nOcrelizumab WB: 12 regions · T cells: 1\nComBat matrix: both tier-1 genes recovered",CT,8.5)
arr(2.625,y, 2.625,y-0.4); arr(7.275,y, 7.275,y-0.4)

y=5.0
box(1.0,y,8.0,1.4,
     "Inverse-Concordant RNA × Methylation Filter (STRICT)\n"
     "RNA↓ × meth↑   (or   RNA↑ × meth↓)\n"
     "BH-FDR<0.05 in BOTH layers · same gene · matched tissue",
     CIN,9.5,'bold')
arr(5,y, 5,y-0.4)

y=3.1
box(0.4,y,2.85,1.30,"STRING physical PPI\n38-gene panel · 53 edges\nadhesion / Th17 hubs\n(Figure 7)",CV,8.5)
box(3.45,y,2.85,1.30,"g:Profiler (g:SCS)\n9 GO:BP immune terms\npan-lymphocyte\nactivation/adhesion",CV,8.5)
box(6.5,y,3.05,1.30,"GWAS + drug overlap\nIMSGC 2019 (200 loci)\nITGB2 + IKZF1 axes\nLXN biomarker",CV,8.5)
for x in [1.825,4.875,8.025]: arr(x,y, x,y-0.4)

# Output: 2 Tier-1 candidates + 10 Tier-2-auxiliary candidates
y=0.55
box(0.6,y,8.8,2.0,
     "Cross-modal candidate genes\n"
     "Tier-1 (2): ITGB2 · IKZF1\n"
     "Tier-2 aux (11): CD79B · LXN · HLA-E · CASP6 · CASP8\n"
     "DGKQ · MX1 · IFIT1 · NUP210 · RUNX3 · SH3BP4",
     CO,9.2,'bold')

out = "__MS_GEO_ROOT__/Poster_v2/figures/workflow_v5_6tier1.png"
plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
print(f"✓ saved → {out}")
