#!/usr/bin/env python3
"""Cell-level Wilcoxon scan over the 25-gene candidate panel (Figure 5B, Methods 4.7 and 4.10).

WHY THIS EXISTS. Results reports 793 gene x cohort x cell-type comparisons of which 111 are
significant, and Figure 5B is drawn from them, but no script in the release produced the table:
Poster_v2/figures/scrna_WILCOXON_v1.tsv was written by an ad-hoc session and nothing regenerated
it. This reproduces it from the three cohort artefacts.

WHAT IT DOES. For each cohort, and within it each cell type retaining at least MIN_CELLS cells in
BOTH arms, every panel gene is tested MS versus control with a tie-corrected Wilcoxon rank-sum test
over the individual cells, via scanpy's rank_genes_groups. Benjamini-Hochberg correction is applied
ONCE across the whole scan, not per cohort or per cell type - the joint correction is what gives
the reported 111 of 793, and correcting within cohort (115) or within cohort x cell type (133)
would not.

WHAT IT IS NOT. These are cell-level tests: the unit is the cell, so the p-values are inflated by
the number of cells rather than by the number of donors, and the manuscript treats them as
directional support only, never as a Tier-1 anchor. The donor-level pseudobulk scripts in this
directory are the honest unit and are what the tier rule uses.

GENE COVERAGE. The Kaufmann cohort's cached UMAP object carries only the 21 panel genes that were
needed for plotting. The remaining four (LYN, MOSPD3, RPAP2, THRB) are read from the deposited
full normalised matrix, which is why that path is required rather than optional.
"""
import os
import sys

import numpy as np
import pandas as pd

ROOT = "__MS_GEO_ROOT__"
SC = os.path.join(ROOT, "SingleCell_CELLxGENE")
OUT = os.path.join(ROOT, "Poster_v2", "figures", "scrna_WILCOXON_rebuilt.tsv")

PANEL = ["CASP6", "CASP8", "CD79B", "CHL1", "CTSZ", "DGKQ", "FOXP3", "HLA-E", "ICAM1", "IFI44L",
         "IFIT1", "IKZF1", "ITGAL", "ITGB2", "LXN", "LYN", "MOSPD3", "MX1", "NUP210", "RPAP2",
         "RUNX3", "SH3BP4", "STAT3", "THRB", "TYK2"]

MIN_CELLS = 20      # per arm; excludes e.g. Beltran NK (19 control cells) and Jakel Macrophages (1)

COHORTS = {
    "Jakel_2019_brain": dict(
        h5ad=os.path.join(SC, "results/figures/brain_Jakel2019/adata_jakel.h5ad"),
        cell_type="celltype", group="group"),
    "Beltran_2019_CSFPBMC": dict(
        h5ad=os.path.join(SC, "results/figures/blood_Beltran2019/adata_beltran.h5ad"),
        cell_type="celltype", group="group"),
    "Ramesh_2020_PBMC": dict(
        h5ad=os.path.join(SC, "results/figures/blood_Ramesh2020_UMAP/adata_ramesh_umap.h5ad"),
        cell_type="basictype", group="disease_status",
        mtx=os.path.join(SC, "data/blood_Ramesh2020_GSE144744/RNA_normalised"),
        meta=os.path.join(SC, "data/blood_Ramesh2020_GSE144744/cell_meta.csv.gz")),
}


def bh(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    q_sorted = np.minimum.accumulate(
        (p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    q = np.empty_like(p)
    q[order] = q_sorted
    return np.clip(q, 0, 1)


def norm_group(s):
    s = pd.Series(s).astype(str).str.upper()
    return s.replace({"CTRL": "HC", "CONTROL": "HC", "HEALTHY": "HC"})


def genes_from_mtx(mtx_dir, meta_path, wanted, barcodes_index):
    """Pull specific genes out of the deposited 10x triplet without loading the whole matrix."""
    import scipy.io as sio
    import scipy.sparse as sp

    genes = pd.read_csv(os.path.join(mtx_dir, "genes.tsv"), sep="\t", header=None)
    sym = genes.iloc[:, -1].astype(str).values
    want = [g for g in wanted if g in set(sym)]
    if not want:
        return None
    rows = {g: int(np.where(sym == g)[0][0]) for g in want}
    bc = pd.read_csv(os.path.join(mtx_dir, "barcodes.tsv"), sep="\t", header=None).iloc[:, 0].astype(str).values
    print(f"    reading {len(want)} genes from the full matrix ({len(bc):,} cells)")
    M = sio.mmread(os.path.join(mtx_dir, "matrix.mtx")).tocsr()
    sub = pd.DataFrame({g: np.asarray(M[rows[g], :].todense()).ravel() for g in want}, index=bc)
    return sub.reindex(barcodes_index)


def scan_cohort(name, cfg):
    import anndata as ad
    import scanpy as sc

    a = ad.read_h5ad(cfg["h5ad"])
    a.obs["_grp"] = norm_group(a.obs[cfg["group"]]).values
    a.obs["_ct"] = a.obs[cfg["cell_type"]].astype(str).values
    have = [g for g in PANEL if g in set(a.var_names)]
    missing = [g for g in PANEL if g not in set(a.var_names)]
    print(f"  {name}: {a.shape[0]:,} cells, {len(have)}/{len(PANEL)} panel genes in the object")

    extra = None
    if missing and cfg.get("mtx"):
        extra = genes_from_mtx(cfg["mtx"], cfg.get("meta"), missing, a.obs_names)
        if extra is not None:
            print(f"    recovered {list(extra.columns)} from the full matrix")
    elif missing:
        print(f"    not measured in this cohort: {missing}")

    a = a[:, have].copy()
    if extra is not None:
        for g in extra.columns:
            a.obs[f"__{g}"] = extra[g].values

    rows = []
    for ct, idx in a.obs.groupby("_ct").groups.items():
        sub = a[idx]
        n_ms = int((sub.obs._grp == "MS").sum())
        n_hc = int((sub.obs._grp == "HC").sum())
        if n_ms < MIN_CELLS or n_hc < MIN_CELLS:
            continue
        s = sub.copy()
        s.obs["_grp"] = s.obs["_grp"].astype("category")
        # use_raw=False is essential: these objects carry a .raw with the full transcriptome, and
        # scanpy prefers it by default, which would silently test every gene instead of the panel.
        sc.tl.rank_genes_groups(s, "_grp", groups=["MS"], reference="HC",
                                method="wilcoxon", tie_correct=True, use_raw=False)
        r = sc.get.rank_genes_groups_df(s, group="MS")
        X = s.X.toarray() if hasattr(s.X, "toarray") else np.asarray(s.X)
        expr = {g: int((X[:, i] > 0).sum()) for i, g in enumerate(s.var_names)}
        for _, x in r.iterrows():
            rows.append(dict(cohort=name, cell_type=ct, gene=x["names"], n_MS=n_ms, n_HC=n_hc,
                             n_expr=expr.get(x["names"], 0),
                             logFC=x["logfoldchanges"], pval=x["pvals"]))
        # genes recovered from the full matrix, tested the same way on the same cells
        from scipy.stats import mannwhitneyu
        for col in [c for c in sub.obs.columns if c.startswith("__")]:
            g = col[2:]
            v = sub.obs[col].astype(float).values
            m, h = v[sub.obs._grp == "MS"], v[sub.obs._grp == "HC"]
            if np.all(np.isnan(m)) or np.all(np.isnan(h)):
                continue
            try:
                p = mannwhitneyu(m, h, alternative="two-sided", use_continuity=True).pvalue
            except ValueError:
                continue
            lfc = float(np.log2((np.nanmean(m) + 1e-9) / (np.nanmean(h) + 1e-9)))
            rows.append(dict(cohort=name, cell_type=ct, gene=g, n_MS=n_ms, n_HC=n_hc,
                             n_expr=int(np.nansum(v > 0)), logFC=lfc, pval=p))
    print(f"    {len(set(x['cell_type'] for x in rows))} cell types passed the "
          f">= {MIN_CELLS}-cells-per-arm filter -> {len(rows)} comparisons")
    return rows


def main():
    allrows = []
    for name, cfg in COHORTS.items():
        if not os.path.exists(cfg["h5ad"]):
            print(f"  {name}: {cfg['h5ad']} missing, skipped")
            continue
        allrows += scan_cohort(name, cfg)

    d = pd.DataFrame(allrows)
    # A gene detected in ZERO cells of a cell type carries no information and its p-value is an
    # artefact of the tie correction. The published table dropped these silently; the rule is made
    # explicit here so the comparison count is reproducible.
    n_all = len(d)
    d = d[d.n_expr > 0].reset_index(drop=True)
    print(f"\n  dropped {n_all - len(d)} comparisons in which the gene was detected in no cell")
    d["fdr"] = bh(d.pval.values)          # ONE correction across the whole scan
    d = d.sort_values(["cohort", "cell_type", "gene"])
    d.to_csv(OUT, sep="\t", index=False)
    print(f"\n  {len(d)} comparisons, {int((d.fdr < 0.05).sum())} significant at BH-FDR < 0.05")
    print(f"  wrote {OUT}")


if __name__ == "__main__":
    sys.exit(main())
