"""
process_GSE144744_kaufmann_genes.py
-----------------------------
Companion to process_GSE144744_kaufmann_pbmc.py — runs on a local laptop where the
full ~1.5–1.8 GB count matrix can be downloaded.

Steps:
  1. If data/blood_Ramesh2020_GSE144744/RNA_normalised.tar.gz is missing,
     stream it from NCBI FTP. (~10 min on a 5 MB/s link.)
  2. Untar into rna_normalised/ (mtx + barcodes + genes).
  3. Build an AnnData, attach paper metadata, compute candidate-gene
     expression per paper cluster name, and write per-cluster Wilcoxon
     stats and the same publication-style figures the other 4 datasets
     produced.

Outputs land in results/figures/blood_Ramesh2020/genes/.
"""
from __future__ import annotations
import gzip, os, time, urllib.request, tarfile, warnings
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests

import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=150, frameon=False, fontsize=9)

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "blood_Ramesh2020_GSE144744"
FIG  = ROOT / "results" / "figures" / "blood_Ramesh2020" / "genes"
FIG.mkdir(parents=True, exist_ok=True)

URL  = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144744/suppl/GSE144744_RNA_normalised.tar.gz"
TAR  = DATA / "RNA_normalised.tar.gz"
EXTR = DATA / "rna_normalised"

CAND = ["LXN", "SH3BP4", "THRB", "CHL1", "RPAP2", "PCNP"]
GROUP_PALETTE = {"MS":"#d62728","HC":"#377eb8"}


def cohens_d(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if len(a) < 2 or len(b) < 2: return np.nan
    s = np.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))
                 /(len(a)+len(b)-2))
    return (a.mean()-b.mean())/s if s else 0.0


def bh_fdr(df, pcol="wilcoxon_p", by="gene"):
    out = df.copy(); out["fdr"] = np.nan
    for k, sub in out.groupby(by):
        _, q, _, _ = multipletests(sub[pcol].fillna(1.0), method="fdr_bh")
        out.loc[sub.index, "fdr"] = q
    return out


def fetch_and_extract():
    if EXTR.exists() and any(EXTR.iterdir()):
        print(f"[exists] {EXTR}")
        return
    if not TAR.exists():
        print(f"Downloading {URL}\n  -> {TAR}  (1.8 GB; this can take ~10 min)")
        t0 = time.time()
        urllib.request.urlretrieve(URL, TAR)
        print(f"  done in {time.time()-t0:.0f}s")
    print(f"Extracting {TAR}")
    EXTR.mkdir(exist_ok=True)
    with tarfile.open(TAR, "r:gz") as tf:
        tf.extractall(EXTR)
    print(f"Files: {[p.name for p in EXTR.rglob('*') if p.is_file()][:5]}")


def main():
    t0 = time.time()
    fetch_and_extract()

    # Find the standard 10x triplet (matrix.mtx, barcodes.tsv, genes.tsv)
    base = next((p for p in EXTR.rglob("matrix.mtx*")), None)
    assert base is not None, "matrix.mtx not found in extraction"
    base = base.parent
    print(f"Using 10x folder: {base}")

    a = sc.read_mtx(base/"matrix.mtx").T  # cells x genes
    barcodes = pd.read_csv(base/"barcodes.tsv", sep="\t", header=None)[0].astype(str).values
    genes = pd.read_csv(base/"genes.tsv", sep="\t", header=None)
    a.obs_names = barcodes
    a.var_names = genes.iloc[:, -1].astype(str).values
    a.var_names_make_unique()

    meta = pd.read_csv(DATA/"cell_meta.csv.gz", low_memory=False)
    meta.set_index("cell_names", inplace=True)
    common = a.obs_names.intersection(meta.index)
    print(f"Cells: matrix={a.n_obs:,}  meta={len(meta):,}  common={len(common):,}")
    a = a[common].copy()
    a.obs = meta.loc[common, ["donor","sample","cohort","group",
                                "basictype","cluster_names"]].copy()
    a.obs["disease_status"] = np.where(
        a.obs["group"].astype(str).str.startswith("HI"), "HC", "MS")
    a.obs.rename(columns={"cluster_names":"paper_cluster",
                            "basictype":"celltype"}, inplace=True)

    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)

    # candidate gene presence
    present = [g for g in CAND if g in a.var_names]
    print(f"Candidate genes present: {present}")

    # stats per paper cluster (Wilcoxon MS vs HC)
    rows = []
    for g in present:
        x = a[:, g].X
        x = x.toarray().ravel() if sp.issparse(x) else np.asarray(x).ravel()
        for ct, idx in a.obs.groupby("paper_cluster").groups.items():
            sel = a.obs_names.get_indexer(idx)
            grp = a.obs.loc[idx,"disease_status"].values
            sub = x[sel]; ms=sub[grp=="MS"]; hc=sub[grp=="HC"]
            if len(ms)<30 or len(hc)<30: continue
            stat, p = ranksums(ms, hc)
            rows.append({"gene":g,"paper_cluster":ct,
                          "n_MS":int(len(ms)),"n_HC":int(len(hc)),
                          "mean_MS":float(ms.mean()),"mean_HC":float(hc.mean()),
                          "logfc":float(ms.mean()-hc.mean()),
                          "cohens_d":float(cohens_d(ms,hc)),
                          "pct_MS":float((ms>0).mean()*100),
                          "pct_HC":float((hc>0).mean()*100),
                          "wilcoxon_p":float(p)})
    df = pd.DataFrame(rows)
    df = bh_fdr(df, pcol="wilcoxon_p", by="gene").sort_values(["gene","fdr"])
    df.to_csv(FIG/"stats_per_paperCluster_MSvsHC.csv", index=False)
    print(f"wrote {FIG/'stats_per_paperCluster_MSvsHC.csv'} ({len(df)} rows)")

    # Quick split-violin per paper cluster
    rawX = a.X
    longs = []
    for g in present:
        x = a[:, g].X
        x = x.toarray().ravel() if sp.issparse(x) else np.asarray(x).ravel()
        longs.append(pd.DataFrame({"gene":g,"expr":x,
                                    "paper_cluster":a.obs["paper_cluster"].astype(str).values,
                                    "group":a.obs["disease_status"].astype(str).values}))
    long = pd.concat(longs, ignore_index=True)
    fig, axes = plt.subplots(2, 3, figsize=(20, 10))
    for ax, g in zip(axes.flat, present):
        sub = long[long["gene"]==g]
        keep = sub.groupby("paper_cluster")["group"].apply(
            lambda s: (s=="MS").sum()>=30 and (s=="HC").sum()>=30)
        sub = sub[sub["paper_cluster"].isin(keep[keep].index)]
        if sub.empty: continue
        order = sorted(sub["paper_cluster"].unique())
        sns.violinplot(data=sub, x="paper_cluster", y="expr", hue="group",
                        palette=GROUP_PALETTE, split=True, inner="quartile",
                        order=order, ax=ax, linewidth=0.4, cut=0)
        ax.set_title(g, fontsize=12, fontweight="bold")
        ax.set_xlabel(""); ax.set_ylabel("log-norm expr")
        for lbl in ax.get_xticklabels():
            lbl.set_rotation(50); lbl.set_ha("right"); lbl.set_fontsize(7)
        if ax is not axes.flat[0] and ax.get_legend() is not None:
            ax.get_legend().remove()
    fig.suptitle("Ramesh 2020 — candidate gene expression by paper cluster, MS vs HC",
                  fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(FIG/"fig6_violin_paperclusters.png", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"done {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
