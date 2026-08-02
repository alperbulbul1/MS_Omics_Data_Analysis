#!/usr/bin/env python3
"""build_data_sources_table.py — comprehensive data-sources table for the manuscript.

Generates Supplementary Table S1 listing every contributing study with:
  - GSE / PRIDE / cohort identifier
  - Tissue + platform + MS / HC sample sizes
  - PMID and citation anchor
  - Stratum used in analysis

Output:
  Poster_v2/figures/SupplementaryTableS1_DataSources.tsv
  Poster_v2/figures/SupplementaryTableS1_DataSources.xlsx (formatted, color-coded)
"""
import json, pandas as pd, numpy as np
from pathlib import Path

PROJ = Path("__MS_GEO_ROOT__")
P = PROJ / "Poster_v2"
OUT = P / "figures"
OUT.mkdir(exist_ok=True)

with open(P / "datasets_used.json") as f:
    reg = json.load(f)

# ════════════════════════════════════════════════════════════════════════════════
# Bulk transcriptomics (14 series)
# ════════════════════════════════════════════════════════════════════════════════
expr_rows = []
for e in reg["expr"]:
    expr_rows.append({
        "Layer": "Bulk transcriptomics",
        "Accession": e["GSE"],
        "Cohort_label": "",                  # will fill from author lookup below
        "Tissue": e["Tissue"],
        "Platform": e["Platform"] + " (" + e["Type"] + ")",
        "n_MS": e["MS"], "n_HC": e["HC"], "n_total": e["MS"] + e["HC"],
        "PMID": e["PMID"],
        "Stratum_used": "",                   # filled below
    })
# Author lookups for bulk
expr_authors = {
    "GSE103005": ("Yang 2022",   "Whole-blood discovery"),
    "GSE137143": ("Mexhitaj 2019","CD4/CD8 T discovery"),
    "GSE138064": ("Hagan 2020", "PBMC + WB"),
    "GSE172009": ("Unpublished",  "CD4 T"),
    "GSE173789": ("Aktas 2022",  "B-cell"),
    "GSE190847": ("Kular 2023",  "B-cell array"),
    "GSE207680": ("Schirmer 2023","Brain WM (small n)"),
    "GSE209596": ("Unpublished",  "T cells (large n)"),
    "GSE211358": ("Glanzman 2022","B cells"),
    "GSE214334": ("Unpublished",  "Brain WM"),
    "GSE21942":  ("Kemppinen 2011","PBMC + WB legacy"),
    "GSE288904": ("2025 release", "Whole blood DMF/Ocrelizumab"),
    "GSE38010":  ("Han 2012",    "Brain WM legacy"),
    "GSE43591":  ("Ottoboni 2013","Whole-blood IFNβ"),
}
for r in expr_rows:
    a = expr_authors.get(r["Accession"], ("—",""))
    r["Cohort_label"] = a[0]; r["Stratum_used"] = a[1]

# ════════════════════════════════════════════════════════════════════════════════
# Methylation (9 series)
# ════════════════════════════════════════════════════════════════════════════════
meth_rows = []
for e in reg["meth"]:
    meth_rows.append({
        "Layer": "DNA methylation",
        "Accession": e["GSE"],
        "Cohort_label": "",
        "Tissue": e["Tissue"],
        "Platform": e["Platform"],
        "n_MS": e["MS"], "n_HC": e["HC"], "n_total": e["MS"] + e["HC"],
        "PMID": e["PMID"],
        "Stratum_used": "",
    })
meth_authors = {
    "GSE106648": ("Kular 2018",       "Whole-blood discovery (largest n)"),
    "GSE130029": ("Maltby 2020",      "CD4 T cells"),
    "GSE130030": ("Maltby 2020",      "CD4/CD8 T cells"),
    "GSE151017": ("Roostaei 2021",    "BAL (relapse triggers)"),
    "GSE189255": ("Souren 2022",      "CD4 T cells — remission"),
    "GSE189256": ("Souren 2022",      "CD14 monocytes"),
    "GSE219293": ("Bos 2023",         "WB DMF + Ocrelizumab"),
    "GSE40360":  ("Huynh 2014",       "Brain WM (chronic active)"),
    "GSE88824":  ("Kular 2018 (sub)", "T cells / WB"),
}
for r in meth_rows:
    a = meth_authors.get(r["Accession"], ("—",""))
    r["Cohort_label"] = a[0]; r["Stratum_used"] = a[1]

# ════════════════════════════════════════════════════════════════════════════════
# Single-cell RNA-seq (3 cohorts: Jäkel, Ramesh, Beltran)
# ════════════════════════════════════════════════════════════════════════════════
scrna_rows = []
for e in reg["scrna"]:
    scrna_rows.append({
        "Layer": "Single-cell RNA-seq",
        "Accession": e["GSE"],
        "Cohort_label": e["Cohort"],
        "Tissue": e["Tissue"],
        "Platform": e["Platform"],
        "n_MS": "—", "n_HC": "—",
        "n_total": f"{e['Cells']} cells / {e['Donors']} donors",
        "PMID": e["PMID"],
        "Stratum_used": "After per-cohort QC + Tier-1 candidate evaluation",
    })
# Add Beltran 2019 (MS-discordant monozygotic twins, CSF + PBMC scRNA)
scrna_rows.append({
    "Layer": "Single-cell RNA-seq",
    "Accession": "GSE138266",
    "Cohort_label": "Beltran 2019",
    "Tissue": "CSF + PBMC (twin pairs)",
    "Platform": "10x Chromium (5′)",
    "n_MS": 5, "n_HC": 5,
    "n_total": "2,029 cells / 10 twins (5 MS-HC pairs)",
    "PMID": "31566580",
    "Stratum_used": "CSF CD8 T cells — HLA-E pseudobulk re-analysis",
})

# ════════════════════════════════════════════════════════════════════════════════
# Proteomics (2 cohorts: Bader/Mann CSF Astral, Wang/Julien brain timsTOF)
# ════════════════════════════════════════════════════════════════════════════════
prot_rows = [
    {"Layer": "Proteomics (CSF Astral)",
     "Accession": "PXD046288",  # Bader/Mann CSF Astral
     "Cohort_label": "Bader, Mann 2024",
     "Tissue": "CSF",
     "Platform": "Orbitrap Astral DIA",
     "n_MS": 978, "n_HC": 306,
     "n_total": 1284,
     "PMID": "38684892",
     "Stratum_used": "DEP::filter_proteins ≥50% + raw Mann-Whitney rescue"},
    {"Layer": "Proteomics (brain timsTOF)",
     "Accession": "PXD060064",
     "Cohort_label": "Wang, Julien 2025",
     "Tissue": "Brain (cortex, NAWM, lesion, lesion-edge)",
     "Platform": "Bruker timsTOF dia-PASEF",
     "n_MS": "20-32", "n_HC": "12-20",
     "n_total": "≈52 (n=5-8 per contrast)",
     "PMID": "40594720",
     "Stratum_used": "Magliozzi-style four-contrast tracking"},
]

# ════════════════════════════════════════════════════════════════════════════════
# Assemble + save
# ════════════════════════════════════════════════════════════════════════════════
all_rows = expr_rows + meth_rows + scrna_rows + prot_rows
df = pd.DataFrame(all_rows)

# Compute totals where numeric
expr_total_ms = sum(r["n_MS"] for r in expr_rows)
expr_total_hc = sum(r["n_HC"] for r in expr_rows)
meth_total_ms = sum(r["n_MS"] for r in meth_rows)
meth_total_hc = sum(r["n_HC"] for r in meth_rows)

print(f"BULK TRANSCRIPTOMICS: {len(expr_rows)} series, "
       f"{expr_total_ms} MS + {expr_total_hc} HC = {expr_total_ms+expr_total_hc} total samples")
print(f"DNA METHYLATION:      {len(meth_rows)} series, "
       f"{meth_total_ms} MS + {meth_total_hc} HC = {meth_total_ms+meth_total_hc} total samples")
print(f"SINGLE-CELL:          {len(scrna_rows)} cohorts (Jäkel, Ramesh, Beltran)")
print(f"PROTEOMICS:           {len(prot_rows)} cohorts (Bader/Mann CSF, Wang/Julien brain)")
print(f"\nTOTAL: {len(df)} datasets")

# Save TSV
tsv_fp = OUT / "SupplementaryTableS1_DataSources.tsv"
df.to_csv(tsv_fp, sep="\t", index=False)
print(f"\n✓ TSV saved → {tsv_fp}")

# Save XLSX with formatting
xlsx_fp = OUT / "SupplementaryTableS1_DataSources.xlsx"
with pd.ExcelWriter(xlsx_fp, engine='openpyxl') as wr:
    df.to_excel(wr, sheet_name='DataSources', index=False)
    # Summary sheet
    summary = pd.DataFrame([
        {"Layer": "Bulk transcriptomics", "Series": 14, "MS samples": expr_total_ms,
          "HC samples": expr_total_hc, "Total samples": expr_total_ms + expr_total_hc},
        {"Layer": "DNA methylation", "Series": 9, "MS samples": meth_total_ms,
          "HC samples": meth_total_hc, "Total samples": meth_total_ms + meth_total_hc},
        {"Layer": "Single-cell RNA-seq", "Series": 3, "MS samples": "≈515k cells / 76 donors",
          "HC samples": "—", "Total samples": "517,533 cells"},
        {"Layer": "Proteomics", "Series": 2, "MS samples": "978 + 20-32",
          "HC samples": "306 + 12-20", "Total samples": "~1,336"},
    ])
    summary.to_excel(wr, sheet_name='Summary', index=False)
print(f"✓ XLSX saved → {xlsx_fp}")

print("\n" + "="*100)
print("SUMMARY FOR MANUSCRIPT")
print("="*100)
print(summary.to_string(index=False))
