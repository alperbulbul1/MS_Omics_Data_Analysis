#!/usr/bin/env python3
"""run_pseudobulk_reanalysis.py — donor-level pseudobulk MS-vs-HC for Tier 1 candidates.

Addresses verdict critical issue #1: cell-level Wilcoxon FDR magnitudes (10⁻¹⁵⁶)
inflated by pseudo-replication. This script:

  1. Aggregates per-cell log-normalised expression to DONOR-level means
     per (gene, cell_type) combination
  2. Runs t-test (Welch) MS vs HC at the donor level
  3. BH-FDR per gene across cell types within cohort
  4. Reports deflation factor: -log10(cell_FDR) - (-log10(donor_FDR))

Cohorts: Beltran (PBMC + CSF) and Ramesh PBMC where Tier 1 gene
significance was previously reported at cell-level resolution.

Tier 1 candidates evaluated: ITGB2, CD79B, SLAMF1, ITGAL, HLA-E,
CTSZ, IFI44L, IKZF1, RPAP2, SH3BP4, CHL1, LXN, PCNP, THRB.
"""
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from scipy.stats import ttest_ind, mannwhitneyu
from statsmodels.stats.multitest import multipletests

PROJ = Path("__MS_GEO_ROOT__")
OUT_DIR = PROJ / "Poster_v2" / "figures"
OUT_DIR.mkdir(exist_ok=True)

TIER1 = ["ITGB2","CD79B","SLAMF1","ITGAL","HLA-E","CTSZ","IFI44L","IKZF1",
         "RPAP2","SH3BP4","CHL1","LXN","PCNP","THRB"]

# ────────────────────────────────────────────────────────────────
# Pseudo-bulk function
# ────────────────────────────────────────────────────────────────
def pseudobulk_test(ad, donor_col, dx_col, ct_col,
                     dx_ms_label="MS", dx_hc_label="HC",
                     min_cells_per_donor=5, min_donors=3,
                     genes=TIER1, cohort_name=""):
    """Donor-mean per (cell_type × gene), then Welch t-test MS vs HC across
    donors. Returns long-format DataFrame with cell-level + donor-level FDRs
    and the deflation factor."""
    X = ad.X
    if issparse(X): X = X.toarray()
    obs = ad.obs.copy()
    var_names = list(ad.var_names)
    gene_idx = {g: i for i, g in enumerate(var_names) if g in genes}
    if not gene_idx:
        print(f"  no Tier 1 genes in {cohort_name}"); return pd.DataFrame()

    cell_types = obs[ct_col].dropna().unique()
    rows = []
    for ct in cell_types:
        mask_ct = (obs[ct_col] == ct).values
        if mask_ct.sum() < 30:  # too few cells in this type
            continue
        ad_ct_X = X[mask_ct]; obs_ct = obs.loc[mask_ct]

        # Aggregate per donor
        for g, gi in gene_idx.items():
            vals = ad_ct_X[:, gi]
            ms_cells = vals[(obs_ct[dx_col] == dx_ms_label).values]
            hc_cells = vals[(obs_ct[dx_col] == dx_hc_label).values]
            n_ms_cells = len(ms_cells); n_hc_cells = len(hc_cells)
            if n_ms_cells < 5 or n_hc_cells < 5: continue

            # Cell-level reference (Mann–Whitney, what we previously reported)
            try:
                if np.var(ms_cells) == 0 and np.var(hc_cells) == 0:
                    cell_p = 1.0
                else:
                    _, cell_p = mannwhitneyu(ms_cells, hc_cells, alternative="two-sided")
            except Exception:
                cell_p = 1.0

            # Donor-level pseudobulk: mean expression per donor
            df_donor = pd.DataFrame({
                'expr': vals,
                'donor': obs_ct[donor_col].values,
                'dx': obs_ct[dx_col].values
            })
            # cell count per donor (drop donors with too few cells)
            counts = df_donor.groupby('donor').size()
            keep_donors = counts[counts >= min_cells_per_donor].index
            df_donor = df_donor[df_donor.donor.isin(keep_donors)]

            donor_means = df_donor.groupby(['donor','dx']).expr.mean().reset_index()
            ms_donors = donor_means[donor_means.dx == dx_ms_label].expr.values
            hc_donors = donor_means[donor_means.dx == dx_hc_label].expr.values
            n_ms_donors = len(ms_donors); n_hc_donors = len(hc_donors)

            if n_ms_donors < min_donors or n_hc_donors < min_donors:
                donor_p = np.nan; donor_t = np.nan; donor_d = np.nan
                logfc = np.nan
            else:
                try:
                    t_stat, donor_p = ttest_ind(ms_donors, hc_donors, equal_var=False)
                    donor_t = float(t_stat)
                except Exception:
                    donor_p = 1.0; donor_t = 0.0
                # Cohen's d on donor-means
                pooled_sd = np.sqrt((np.var(ms_donors, ddof=1) + np.var(hc_donors, ddof=1)) / 2)
                donor_d = float((np.mean(ms_donors) - np.mean(hc_donors)) / pooled_sd) if pooled_sd > 0 else 0.0
                logfc = float(np.mean(ms_donors) - np.mean(hc_donors))

            rows.append({
                'cohort': cohort_name, 'cell_type': ct, 'gene': g,
                'n_ms_cells': int(n_ms_cells), 'n_hc_cells': int(n_hc_cells),
                'n_ms_donors': int(n_ms_donors), 'n_hc_donors': int(n_hc_donors),
                'cell_wilcox_p': float(cell_p),
                'donor_logfc': logfc, 'donor_cohens_d': donor_d,
                'donor_t': donor_t, 'donor_ttest_p': float(donor_p) if not pd.isna(donor_p) else np.nan,
            })

    df = pd.DataFrame(rows)
    if df.empty: return df

    # BH-FDR per gene (across cell types within cohort)
    df['cell_wilcox_fdr'] = np.nan
    df['donor_ttest_fdr'] = np.nan
    for g in df.gene.unique():
        m = df.gene == g
        # cell-level
        _, q_cell, _, _ = multipletests(df.loc[m, 'cell_wilcox_p'].fillna(1.0).values,
                                          method='fdr_bh')
        df.loc[m, 'cell_wilcox_fdr'] = q_cell
        # donor-level (skip NaN — too few donors)
        donor_p = df.loc[m, 'donor_ttest_p'].values
        valid = ~np.isnan(donor_p)
        if valid.sum() > 0:
            q_donor = np.full(len(donor_p), np.nan)
            _, q_donor[valid], _, _ = multipletests(donor_p[valid], method='fdr_bh')
            df.loc[m, 'donor_ttest_fdr'] = q_donor
    # Deflation: how many orders of magnitude did FDR shrink?
    df['log10_deflation'] = -np.log10(df.cell_wilcox_fdr.clip(1e-300)) - \
                             (-np.log10(df.donor_ttest_fdr.clip(1e-300)))
    return df

# ════════════════════════════════════════════════════════════════
# Run for Beltran (PBMC + CSF — KEY: HLA-E CSF CD8 signal)
# ════════════════════════════════════════════════════════════════
print("=" * 75)
print("PSEUDO-BULK REANALYSIS — Tier 1 candidates × cell types × donors")
print("=" * 75)

print("\n[1/2] Beltran (PBMC + CSF) ...")
ad = sc.read_h5ad(PROJ / "SingleCell_CELLxGENE" / "results" / "figures" /
                   "blood_Beltran2019" / "adata_beltran.h5ad")
# Note: 'group' column gives MS/HC consolidated; 'twin' is donor; 'celltype' or 'compartment' as cell context
# We test by compartment × celltype combinations
ad.obs['ct_full'] = ad.obs['compartment'].astype(str) + '_' + ad.obs['celltype'].astype(str)
beltran = pseudobulk_test(ad, donor_col='twin', dx_col='group', ct_col='ct_full',
                            dx_ms_label='MS', dx_hc_label='HC',
                            cohort_name="Beltran_2019")
print(f"  → {len(beltran)} rows")

# ════════════════════════════════════════════════════════════════
# Run for Ramesh PBMC (42 donors)
# ════════════════════════════════════════════════════════════════
print("\n[2/2] Ramesh PBMC ...")
ad = sc.read_h5ad(PROJ / "SingleCell_CELLxGENE" / "results" / "figures" /
                   "blood_Ramesh2020_UMAP" / "adata_ramesh_umap.h5ad")
ramesh = pseudobulk_test(ad, donor_col='donor', dx_col='disease_status', ct_col='basictype',
                          dx_ms_label='MS', dx_hc_label='HC',
                          cohort_name="Ramesh_2020")
print(f"  → {len(ramesh)} rows")

# ════════════════════════════════════════════════════════════════
# Combine + save
# ════════════════════════════════════════════════════════════════
combined = pd.concat([beltran, ramesh], ignore_index=True)
out_fp = OUT_DIR / "scrna_DONOR_PSEUDOBULK_reanalysis.tsv"
combined.to_csv(out_fp, sep='\t', index=False)
print(f"\n✓ saved {out_fp}  ({len(combined)} total rows)")

# ════════════════════════════════════════════════════════════════
# Report: most-deflated signals
# ════════════════════════════════════════════════════════════════
print("\n" + "=" * 75)
print("DEFLATION REPORT — cell-level vs donor-level FDR for key claims")
print("=" * 75)

# Sort by cell_fdr (most-significant first) and show what happens at donor level
key_signals = combined.sort_values('cell_wilcox_fdr').head(25).copy()
key_signals['cell_log10'] = -np.log10(key_signals.cell_wilcox_fdr.clip(1e-300))
key_signals['donor_log10'] = -np.log10(key_signals.donor_ttest_fdr.clip(1e-300))
print(key_signals[['cohort','cell_type','gene','n_ms_donors','n_hc_donors',
                    'cell_wilcox_fdr','donor_ttest_fdr','log10_deflation']].to_string(index=False))

# HLA-E spotlight
print("\n" + "─" * 75)
print("HLA-E CSF spotlight (key manuscript claim — was FDR=2.89e-18)")
print("─" * 75)
hla = combined[(combined.gene == 'HLA-E') &
                (combined.cell_type.str.contains('CSF'))]
if len(hla):
    print(hla[['cell_type','n_ms_cells','n_hc_cells','n_ms_donors','n_hc_donors',
                'donor_logfc','donor_cohens_d','cell_wilcox_fdr','donor_ttest_p',
                'donor_ttest_fdr']].to_string(index=False))

# ITGAL spotlight
print("\n" + "─" * 75)
print("ITGAL t_cells spotlight (was FDR=4.5e-156 in mark6)")
print("─" * 75)
itgal = combined[(combined.gene == 'ITGAL') & (combined.cohort == 'Ramesh_2020')]
if len(itgal):
    print(itgal[['cell_type','n_ms_cells','n_hc_cells','n_ms_donors','n_hc_donors',
                  'donor_logfc','donor_cohens_d','cell_wilcox_fdr','donor_ttest_p',
                  'donor_ttest_fdr']].to_string(index=False))

# ITGB2 spotlight
print("\n" + "─" * 75)
print("ITGB2 t_cells spotlight (was FDR=1.1e-37 in mark6)")
print("─" * 75)
itgb = combined[(combined.gene == 'ITGB2') & (combined.cohort == 'Ramesh_2020')]
if len(itgb):
    print(itgb[['cell_type','n_ms_cells','n_hc_cells','n_ms_donors','n_hc_donors',
                 'donor_logfc','donor_cohens_d','cell_wilcox_fdr','donor_ttest_p',
                 'donor_ttest_fdr']].to_string(index=False))
