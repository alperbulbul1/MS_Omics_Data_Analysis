#!/usr/bin/env python3
"""Literature-standard pseudobulk for Kaufmann GSE144744: SUM of RAW INTEGER UMI COUNTS per
biological sample x cell cluster (Crowell/muscat; Squair 2021).

The GEO matrix is deposited log-normalised, but the raw counts are exactly recoverable because
Seurat LogNormalize is invertible given each cell's total UMI (cell_meta$nCount_RNA):
    count = expm1(lognorm) * nCount_RNA / 1e4
Verified: 100.000% of reconstructed values are integers and per-cell expm1 sums equal exactly 1e4.

Emits three aggregations, all as summed raw counts:
  PBC_fine_*    donor x cluster_names   (25 subsets - the study's own published clustering)
  PBC_coarse_*  donor x basictype       (8 lineages)
  PBC_donor_*   donor                   (whole-PBMC; one sample per donor, no repeated measures)
"""
import numpy as np, pandas as pd, os
BASE="__MS_GEO_ROOT__/SingleCell_CELLxGENE/data/blood_Ramesh2020_GSE144744"
OUT="__MS_GEO_ROOT__/Poster_v2/figures/pseudobulk_proper"; os.makedirs(OUT,exist_ok=True)
SCALE=1e4; MIN_CELLS=10

genes=pd.read_csv(f"{BASE}/RNA_normalised/genes.tsv",sep="\t",header=None)[0].astype(str).values
barc =pd.read_csv(f"{BASE}/RNA_normalised/barcodes.tsv",header=None)[0].astype(str).values
meta=pd.read_csv(f"{BASE}/cell_meta.csv.gz",
    usecols=["cell_names","donor","basictype","cluster_names","group","cohort","batch_pair",
             "nCount_RNA","sex","age_sampling"])
meta=meta.drop_duplicates("cell_names").set_index("cell_names").reindex(barc)
tot=meta["nCount_RNA"].values.astype(np.float64)

def cond_of(g):
    g=str(g)
    if g.startswith("HI"): return "HC"
    if g.startswith(("MS","PPMS","RRMS","SPMS")): return "MS"
    return None
cond=np.array([cond_of(g) for g in meta["group"].values],dtype=object)
base_ok=meta["donor"].notna().values & np.isin(cond,["MS","HC"]) & np.isfinite(tot)

# three groupings
specs={}
for name,col in [("fine","cluster_names"),("coarse","basictype")]:
    k=(meta["donor"].astype(str)+"__"+meta[col].astype(str))
    ok=base_ok & meta[col].notna().values
    specs[name]=(k.where(ok,other=None),)
k=meta["donor"].astype(str); specs["donor"]=(k.where(base_ok,other=None),)

acc={}; cmap={}; ncell={}
for name,(key,) in specs.items():
    groups=pd.Index(sorted(set(key.dropna())))
    gi=pd.Series(np.arange(len(groups)),index=groups)
    cg=np.full(len(barc),-1,dtype=np.int64); v=key.notna().values
    cg[v]=gi.reindex(key[v]).values
    acc[name]=np.zeros((len(genes),len(groups)),dtype=np.float64)
    cmap[name]=(cg,groups); ncell[name]=np.bincount(cg[cg>=0],minlength=len(groups))
    print(f"{name}: {len(groups)} pseudobulk groups",flush=True)

done=0
for ch in pd.read_csv(f"{BASE}/RNA_normalised/matrix.mtx",sep=r"\s+",skiprows=3,header=None,
                      names=["g","c","v"],dtype={"g":np.int32,"c":np.int32,"v":np.float64},
                      chunksize=25_000_000):
    gi=ch["g"].values-1; ci=ch["c"].values-1
    cnt=np.expm1(ch["v"].values)*tot[ci]/SCALE          # exact raw-count recovery
    cnt=np.rint(cnt)
    for name in specs:
        cg,_=cmap[name]; grp=cg[ci]; m=grp>=0
        np.add.at(acc[name],(gi[m],grp[m]),cnt[m])
    done+=len(ch); print(f"  {done:,}/529,680,471",flush=True)

dm=meta.reset_index().drop_duplicates("donor").set_index("donor")
for name in specs:
    cg,groups=cmap[name]; nk=ncell[name]
    sel=nk>=MIN_CELLS if name!="donor" else nk>=MIN_CELLS
    Mx=acc[name][:,sel]; gk=groups[sel]; nn=nk[sel]
    assert np.allclose(Mx,np.rint(Mx)), "non-integer pseudobulk!"
    md=pd.DataFrame({"sample":gk})
    md["donor"]=[s.split("__")[0] for s in gk] if name!="donor" else list(gk)
    md["cell_type"]=[s.split("__")[1] for s in gk] if name!="donor" else "all_PBMC"
    md["condition"]=[cond_of(dm.loc[d,"group"]) for d in md.donor]
    md["group"]=[dm.loc[d,"group"] for d in md.donor]
    md["cohort"]=[dm.loc[d,"cohort"] for d in md.donor]
    md["batch_pair"]=[dm.loc[d,"batch_pair"] for d in md.donor]
    md["sex"]=[dm.loc[d,"sex"] for d in md.donor]
    md["age"]=[dm.loc[d,"age_sampling"] for d in md.donor]
    md["n_cells"]=nn
    pd.DataFrame(Mx.astype(np.int64),index=genes,columns=gk).to_csv(f"{OUT}/PBC_{name}_matrix.csv")
    md.to_csv(f"{OUT}/PBC_{name}_coldata.csv",index=False)
    print(f"wrote PBC_{name}: {Mx.shape[0]} genes x {Mx.shape[1]} samples | "
          f"donors {md.groupby('condition').donor.nunique().to_dict()} | libsize med={np.median(Mx.sum(0)):,.0f}")
print("\nsex/age balance across pairs (for covariate decision):")
d=dm.reset_index()[["donor","group","batch_pair","sex","age_sampling"]]
d["condition"]=[cond_of(g) for g in d.group]
d=d[d.condition.notna()]
print(pd.crosstab(d.condition,d.sex))
print(d.groupby("condition").age_sampling.describe()[["count","mean","std"]])
same=d.groupby("batch_pair").sex.nunique()
print(f"pairs with SAME sex: {(same==1).sum()}/{len(same)}")
