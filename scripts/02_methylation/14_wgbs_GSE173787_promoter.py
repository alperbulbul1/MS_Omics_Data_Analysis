#!/usr/bin/env python3
"""Sorted-cell WGBS promoter assessment (GSE173787), Methods 4.4 and Results 2.4.

WHY THIS EXISTS. Methods describes this analysis and Results reports its outcome, but no script in
the release performed it, so the reported numbers could not be re-derived from the deposit.

WHAT IT DOES. GSE173787 provides whole-genome bisulfite sequencing of four sorted immune
populations (CD19 B cells, CD4 T, CD8 T, CD14 monocytes) from MS patients and controls. Promoter
methylation was aggregated over TSS +/- 2 kb windows (GENCODE v49) into
Methylation_Data/GSE173787_promoter_meth_long.tsv, a long table of
sample / gene / mean-promoter-methylation / n-CpG. This script performs the test on that table:
per cell type, Welch's t-test of MS versus control for each candidate gene, with Benjamini-Hochberg
correction WITHIN cell type. That is the "tested per cell type with Welch tests and BH correction"
of Methods 4.4.

WHAT IT ESTABLISHES. Neither inverse-concordant Tier-1 promoter reaches FDR < 0.05 in any of the
four sorted populations - 0 of 8 gene x cell-type combinations - which is why the manuscript calls
this an orthogonal ASSESSMENT rather than a validation. The promoters that do pass are the
auxiliary and context genes RUNX3, CTSZ and ICAM1 in B cells, plus ICAM1 in CD4 T cells; CD79B is
nominal in B cells (delta-beta -0.014, p = 0.025) but does not survive correction.

The BH denominator is the number of genes actually testable in that cell type (>= 3 samples per
arm), which is 18 for all four populations in the deposited table. This matters: under a 25-gene
denominator only RUNX3 survives. The denominator is printed so the choice is visible.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

ROOT = "__MS_GEO_ROOT__"
IN = os.path.join(ROOT, "Methylation_Data", "GSE173787_promoter_meth_long.tsv")
OUT = os.path.join(ROOT, "Methylation", "results", "14_wgbs_GSE173787_promoter.tsv")

TIER1 = ["ITGB2", "IKZF1"]
MIN_PER_ARM = 3


def bh(p):
    """Benjamini-Hochberg. The rank must be the rank of each p-value in SORTED order; applying
    np.arange in the original element order silently inflates or deflates individual q-values."""
    p = np.asarray(p, float)
    order = np.argsort(p)
    ranked = p[order]
    q_sorted = np.minimum.accumulate(
        (ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    q = np.empty_like(p)
    q[order] = q_sorted
    return np.clip(q, 0, 1)


def main():
    d = pd.read_csv(IN, sep="\t", header=None,
                    names=["sample", "gene", "meth", "n_cpg"])
    # sample ids are GSM<id>_<group>-<donor>-<celltype>
    parts = d["sample"].str.extract(r"^(GSM\d+)_(MS|HC)-(\d+)-(CD\d+)$")
    d["group"], d["donor"], d["cell"] = parts[1], parts[2], parts[3]
    if d["group"].isna().any():
        raise ValueError("unparsed sample identifiers: "
                         f"{d.loc[d.group.isna(), 'sample'].unique()[:5]}")

    rows = []
    for cell, g1 in d.groupby("cell"):
        genes, pvals, deltas, ns = [], [], [], []
        for gene, g2 in g1.groupby("gene"):
            a = g2.loc[g2.group == "MS", "meth"].dropna()
            b = g2.loc[g2.group == "HC", "meth"].dropna()
            if len(a) < MIN_PER_ARM or len(b) < MIN_PER_ARM:
                continue
            genes.append(gene)
            pvals.append(stats.ttest_ind(a, b, equal_var=False)[1])
            deltas.append(a.mean() - b.mean())
            ns.append((len(a), len(b)))
        q = bh(pvals)
        print(f"  {cell}: {len(genes)} genes testable (BH denominator), "
              f"{g1['sample'].nunique()} samples")
        for gene, p, qq, dl, (na, nb) in zip(genes, pvals, q, deltas, ns):
            rows.append(dict(cell_type=cell, gene=gene, n_MS=na, n_HC=nb,
                             delta_beta=dl, p_value=p, BH_FDR=qq))

    r = pd.DataFrame(rows).sort_values(["cell_type", "BH_FDR"])
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    r.to_csv(OUT, sep="\t", index=False)

    sig = r[r.BH_FDR < 0.05]
    print(f"\n  FDR < 0.05: {len(sig)} gene x cell-type combinations")
    for _, x in sig.iterrows():
        print(f"    {x.gene:<8} {x.cell_type:<6} delta-beta {x.delta_beta:+.4f}  "
              f"p {x.p_value:.4g}  FDR {x.BH_FDR:.3f}")

    t1 = r[r.gene.isin(TIER1)]
    n_sig_t1 = int((t1.BH_FDR < 0.05).sum())
    print(f"\n  Tier-1 promoters ({', '.join(TIER1)}): {n_sig_t1} of {len(t1)} "
          "gene x cell-type combinations reach FDR < 0.05")
    for _, x in t1.sort_values(["gene", "BH_FDR"]).iterrows():
        print(f"    {x.gene:<8} {x.cell_type:<6} delta-beta {x.delta_beta:+.4f}  "
              f"p {x.p_value:.4g}  FDR {x.BH_FDR:.3f}")
    if n_sig_t1:
        print("\n  NOTE: the manuscript states that no Tier-1 promoter is significant in any sorted "
              "population. That statement no longer holds for this input.")
    print(f"\n  wrote {OUT}")


if __name__ == "__main__":
    sys.exit(main())
