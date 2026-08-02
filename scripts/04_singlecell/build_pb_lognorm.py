#!/usr/bin/env python3
"""Build donor x cell-type pseudobulk from the DEPOSITED LOG-NORMALISED Kaufmann matrix
(GSE144744 RNA_normalised: MatrixMarket 'real', 15354 genes x 497705 cells).

Because the deposited values are log-normalised expression (NOT counts), the statistically
appropriate pseudobulk summary is the MEAN of the log-normalised values per donor x cell type,
which is then modelled with limma-trend (Law et al.; Squair et al. 2021) -- voom is NOT
applicable to already-normalised data since it needs counts + library sizes.

Streams the 529,680,471 non-zero entries in chunks and accumulates group sums.
Outputs: PB_lognorm_matrix.csv (genes x donor__celltype) + PB_lognorm_coldata.csv
"""
import numpy as np, pandas as pd, os, sys
BASE="__MS_GEO_ROOT__/SingleCell_CELLxGENE/data/blood_Ramesh2020_GSE144744"
MTX=f"{BASE}/RNA_normalised/matrix.mtx"
OUT="__MS_GEO_ROOT__/Poster_v2/figures/pseudobulk_proper"
MIN_CELLS=10

genes=pd.read_csv(f"{BASE}/RNA_normalised/genes.tsv",sep="\t",header=None)[0].astype(str).values
barc =pd.read_csv(f"{BASE}/RNA_normalised/barcodes.tsv",header=None)[0].astype(str).values
print(f"genes={len(genes)} cells={len(barc)}",flush=True)

meta=pd.read_csv(f"{BASE}/cell_meta.csv.gz",
                 usecols=["cell_names","donor","basictype","group","cohort","batch_pair"])
meta=meta.drop_duplicates("cell_names").set_index("cell_names")
meta=meta.reindex(barc)                                  # align to matrix column order
print(f"metadata matched: {meta['donor'].notna().sum()}/{len(barc)}",flush=True)

def cond_of(g):
    g=str(g)
    if g.startswith("HI"): return "HC"
    if g.startswith(("MS","PPMS","RRMS","SPMS")): return "MS"
    return None
cond=np.array([cond_of(g) for g in meta["group"].values],dtype=object)

key=pd.Series(meta["donor"].astype(str)+"__"+meta["basictype"].astype(str))
ok=pd.notna(meta["donor"]).values & pd.notna(meta["basictype"]).values & np.isin(cond,["MS","HC"])
key=key.where(ok, other=None)
groups=pd.Index(sorted(set(key.dropna())))
gidx=pd.Series(np.arange(len(groups)),index=groups)
col_group=np.full(len(barc),-1,dtype=np.int64)
valid=key.notna().values
col_group[valid]=gidx.reindex(key[valid]).values
ncell=np.bincount(col_group[col_group>=0],minlength=len(groups))
print(f"groups (donor x celltype) = {len(groups)} | cells assigned = {(col_group>=0).sum()}",flush=True)

acc=np.zeros((len(genes),len(groups)),dtype=np.float64)
rows_done=0
reader=pd.read_csv(MTX,sep=r"\s+",skiprows=3,header=None,names=["g","c","v"],
                   dtype={"g":np.int32,"c":np.int32,"v":np.float32},chunksize=25_000_000)
for ch in reader:
    gi=ch["g"].values-1; ci=ch["c"].values-1; vv=ch["v"].values.astype(np.float64)
    grp=col_group[ci]
    m=grp>=0
    np.add.at(acc,(gi[m],grp[m]),vv[m])
    rows_done+=len(ch)
    print(f"  {rows_done:,}/529,680,471 entries",flush=True)

mean=acc/np.maximum(ncell,1)[None,:]                     # MEAN log-normalised expression
sel=ncell>=MIN_CELLS
print(f"groups kept (>={MIN_CELLS} cells): {sel.sum()}/{len(groups)}",flush=True)
mean=mean[:,sel]; gkept=groups[sel]; nk=ncell[sel]

md=pd.DataFrame({"sample":gkept})
md["donor"]=[s.split("__")[0] for s in gkept]
md["cell_type"]=[s.split("__")[1] for s in gkept]
dm=meta.reset_index().drop_duplicates("donor").set_index("donor")
md["condition"]=[cond_of(dm.loc[d,"group"]) for d in md.donor]
md["group"]=[dm.loc[d,"group"] for d in md.donor]
md["cohort"]=[dm.loc[d,"cohort"] for d in md.donor]
md["batch_pair"]=[dm.loc[d,"batch_pair"] for d in md.donor]
md["n_cells"]=nk

os.makedirs(OUT,exist_ok=True)
pd.DataFrame(mean,index=genes,columns=gkept).to_csv(f"{OUT}/PB_lognorm_matrix.csv")
md.to_csv(f"{OUT}/PB_lognorm_coldata.csv",index=False)
print(f"\nwrote {OUT}/PB_lognorm_matrix.csv  ({mean.shape[0]} genes x {mean.shape[1]} pseudobulk samples)")
print("donors:",md.groupby('condition').donor.nunique().to_dict())
print("cell types:",md.cell_type.nunique(),"| pairs:",md.batch_pair.nunique())
