#!/usr/bin/env python3
"""Normalisation-based pseudobulk aggregation, as a sensitivity comparison against the
summed-raw-count (muscat standard) aggregation.

The deposited Kaufmann matrix is log-normalised: lognorm = log1p(count / total * 1e4).
Therefore expm1(lognorm) = count/total*1e4 = the cell's CP10K (counts per 10,000) value, i.e.
the depth-normalised expression. For 10x 3' tag counting, TPM reduces to CPM because effective
transcript length is not identifiable from 3'-biased reads, so CP10K/CPM is the appropriate
normalised unit here (CPM = CP10K * 100).

Aggregation schemes produced (each per donor x cell type):
  MEANCP10K : mean over cells of the per-cell CP10K  -> "normalise each cell, then average"
  SUMCP10K  : sum  over cells of the per-cell CP10K  -> "normalise each cell, then sum"
                (differs from mean only by the per-group cell count)
Both weight every cell equally regardless of sequencing depth, unlike summing raw counts which
weights deep cells more. Downstream these go to limma-trend on log2, since they are not counts.
"""
import numpy as np, pandas as pd, os
BASE="__MS_GEO_ROOT__/SingleCell_CELLxGENE/data/blood_Ramesh2020_GSE144744"
OUT="__MS_GEO_ROOT__/Poster_v2/figures/pseudobulk_proper"
MIN_CELLS=10

genes=pd.read_csv(f"{BASE}/RNA_normalised/genes.tsv",sep="\t",header=None)[0].astype(str).values
barc =pd.read_csv(f"{BASE}/RNA_normalised/barcodes.tsv",header=None)[0].astype(str).values
meta=pd.read_csv(f"{BASE}/cell_meta.csv.gz",
    usecols=["cell_names","donor","basictype","cluster_names","group","cohort","batch_pair","sex","age_sampling"])
meta=meta.drop_duplicates("cell_names").set_index("cell_names").reindex(barc)

def cond_of(g):
    g=str(g)
    if g.startswith("HI"): return "HC"
    if g.startswith(("MS","PPMS","RRMS","SPMS")): return "MS"
    return None
cond=np.array([cond_of(g) for g in meta["group"].values],dtype=object)
base_ok=meta["donor"].notna().values & np.isin(cond,["MS","HC"])

specs={}
for name,col in [("coarse","basictype"),("fine","cluster_names")]:
    k=(meta["donor"].astype(str)+"__"+meta[col].astype(str))
    specs[name]=k.where(base_ok & meta[col].notna().values, other=None)
specs["donor"]=meta["donor"].astype(str).where(base_ok, other=None)

acc={}; cmap={}; ncell={}
for name,key in specs.items():
    groups=pd.Index(sorted(set(key.dropna())))
    gi=pd.Series(np.arange(len(groups)),index=groups)
    cg=np.full(len(barc),-1,dtype=np.int64); v=key.notna().values
    cg[v]=gi.reindex(key[v]).values
    acc[name]=np.zeros((len(genes),len(groups)),dtype=np.float64)
    cmap[name]=(cg,groups); ncell[name]=np.bincount(cg[cg>=0],minlength=len(groups))
    print(f"{name}: {len(groups)} groups",flush=True)

done=0
for ch in pd.read_csv(f"{BASE}/RNA_normalised/matrix.mtx",sep=r"\s+",skiprows=3,header=None,
                      names=["g","c","v"],dtype={"g":np.int32,"c":np.int32,"v":np.float64},
                      chunksize=25_000_000):
    gi=ch["g"].values-1; ci=ch["c"].values-1
    cp=np.expm1(ch["v"].values)                 # per-cell CP10K (depth-normalised)
    for name in specs:
        cg,_=cmap[name]; grp=cg[ci]; m=grp>=0
        np.add.at(acc[name],(gi[m],grp[m]),cp[m])
    done+=len(ch); print(f"  {done:,}/529,680,471",flush=True)

dm=meta.reset_index().drop_duplicates("donor").set_index("donor")
for name in specs:
    _,groups=cmap[name]; nk=ncell[name]; sel=nk>=MIN_CELLS
    S=acc[name][:,sel]; gk=groups[sel]; nn=nk[sel]
    md=pd.DataFrame({"sample":gk})
    md["donor"]=[s.split("__")[0] for s in gk] if name!="donor" else list(gk)
    md["cell_type"]=[s.split("__")[1] for s in gk] if name!="donor" else "all_PBMC"
    md["condition"]=[cond_of(dm.loc[d,"group"]) for d in md.donor]
    md["group"]=[dm.loc[d,"group"] for d in md.donor]
    md["batch_pair"]=[dm.loc[d,"batch_pair"] for d in md.donor]
    md["n_cells"]=nn
    # MEAN CP10K (equal weight per cell); CPM = *100
    pd.DataFrame(S/nn[None,:],index=genes,columns=gk).to_csv(f"{OUT}/PBN_{name}_meanCP10K.csv")
    # SUM CP10K
    pd.DataFrame(S,index=genes,columns=gk).to_csv(f"{OUT}/PBN_{name}_sumCP10K.csv")
    md.to_csv(f"{OUT}/PBN_{name}_coldata.csv",index=False)
    print(f"wrote PBN_{name}: {S.shape[0]} x {S.shape[1]} | meanCP10K col-sums med={np.median((S/nn[None,:]).sum(0)):,.0f} (expect ~1e4)")
print("\nsanity: mean-CP10K column sums should be ~10,000 by construction")
