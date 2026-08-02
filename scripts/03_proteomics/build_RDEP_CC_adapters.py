#!/usr/bin/env python3
"""Build figure-adapter tables from the complete-case (no-imputation) proteomic re-analysis.

Rationale: MinProb imputation assumes missingness is left-censored (MNAR) and fills absent
values with low intensities. Where missingness differs between groups this manufactures an
effect. ITGB2 in the Astral CSF cohort is measured in 71.6% of MS but only 63.1% of control
samples, so MinProb injects proportionally more low values into controls and returns a
spurious MS-up call (log2FC=+0.166, FDR=0.023) that the measured values do not support
(complete-case log2FC=+0.046, FDR=0.55; raw median difference +0.015). The same imputation
suppresses a real effect in the opposite direction for LXN in white-matter lesions, where
4 of 5 lesion samples are measured and the single absent value is filled low.

Complete-case (MAR) is therefore used for all three re-analysed cohorts uniformly. Proteins
with 0% missingness (CTSZ, ICAM1) are unaffected, confirming the change is confined to the
incomplete rows.

Output mirrors the legacy column layout consumed by figure4_proteomics.py and
figure6_intersection_heatmap.py, so only the input directory changes.
"""
import pandas as pd
from pathlib import Path

META = Path("__MS_GEO_ROOT__/Proteomics/processed/META")
OUT = Path("__MS_GEO_ROOT__/Proteomics/processed/RDEP_CC")
OUT.mkdir(exist_ok=True)

# Astral and timsTOF use different fold-change column names in the legacy layout; keep both.
CSF = {"Astral_RDEP.tsv": ("CSF_Astral_CC_results.tsv", "log2FC_MSvsCtrl"),
       "timsTOF_RDEP.tsv": ("CSF_timsTOF_CC_results.tsv", "log2FC")}
MAG = ["MS_CTX_vs_ODC_CTX", "MS_NAWM_vs_ODC_WM", "MS_WML_vs_ODC_WM", "MS_WML_vs_MS_NAWM"]

for dst, (src, fccol) in CSF.items():
    d = pd.read_csv(META / src, sep="\t")
    out = pd.DataFrame({"Genes": d["gene"], fccol: d["logFC"],
                        "t": d["t"], "pval": d["P.Value"], "FDR": d["adj.P.Val"]})
    out.to_csv(OUT / dst, sep="\t", index=False)
    print(f"{dst:22s} {len(out):5d} genes  FDR<0.05={int((out.FDR<0.05).sum())}")

for cn in MAG:
    d = pd.read_csv(META / f"Magliozzi_CC_{cn}.tsv", sep="\t")
    out = pd.DataFrame({"Gene": d["gene"], "log2FC": d["logFC"],
                        "t": d["t"], "p": d["P.Value"], "FDR": d["adj.P.Val"]})
    out.to_csv(OUT / f"Magliozzi2026_{cn}.tsv", sep="\t", index=False)
    print(f"{cn:22s} {len(out):5d} genes  nominal p<0.05={int((out.p<0.05).sum())}")

print(f"\nwrote {OUT}")
