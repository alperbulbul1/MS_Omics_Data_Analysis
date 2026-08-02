#!/usr/bin/env python3
"""Literature-standard pseudobulk inputs for the brain (Jaekel GSE118257) and CSF/blood
(Beltran GSE127969) cohorts, so all three single-cell cohorts are processed the same way.

Jaekel  : expr.txt.gz holds RAW integer UMI counts -> SUM per patient x cell type (muscat standard).
          Region samples are collapsed to the patient (MS176_NAWM/_RM/_CA -> MS176) because the
          biological replicate is the donor, not the lesion block.
Beltran : only TPM-like normalised values were deposited (values are non-integer and no per-cell
          library size is recorded), so raw counts are NOT recoverable. Aggregated as the MEAN of the
          deposited normalised values per twin x cell type and analysed with limma-trend; this is
          reported as a normalisation-based aggregation, not a count model.
"""
import numpy as np, pandas as pd, re, os, anndata as ad
OUT="__MS_GEO_ROOT__/Poster_v2/figures/pseudobulk_proper"
MIN_CELLS=10

# ---------------- Jaekel: raw counts ----------------
J="__MS_GEO_ROOT__/SingleCell_CELLxGENE/data/brain_Jakel2019_GSE118257"
print("Jaekel: reading raw count matrix...",flush=True)
expr=pd.read_csv(f"{J}/expr.txt.gz",sep="\t",index_col=0)
anno=pd.read_csv(f"{J}/anno.txt.gz",sep="\t",index_col=0)
common=[c for c in expr.columns if c in anno.index]
expr=expr[common]; anno=anno.loc[common]
assert np.allclose(expr.values[:500],np.rint(expr.values[:500])), "Jaekel not integer!"
def patient(s): return re.sub(r"[_/](A\d?|CA\d?|CI|NAWM|WM|GM|RM|RIM|LR|\d+)$","",str(s))
anno["patient"]=anno["Sample"].map(patient)
anno["condition"]=np.where(anno["Condition"].astype(str).str.upper().str.startswith("MS"),"MS","HC")
print(f"  {expr.shape[0]} genes x {expr.shape[1]} cells | patients: "
      f"{anno.groupby('condition').patient.nunique().to_dict()} | cell types: {anno.Celltypes.nunique()}",flush=True)

key=anno["patient"].astype(str)+"__"+anno["Celltypes"].astype(str)
groups=pd.Index(sorted(key.unique()))
gi=pd.Series(np.arange(len(groups)),index=groups)
idx=gi.reindex(key.values).values
M=np.zeros((expr.shape[0],len(groups)))
X=expr.values
for j in range(len(groups)):
    sel=idx==j
    if sel.sum(): M[:,j]=X[:,sel].sum(1)
nc=np.bincount(idx,minlength=len(groups))
sel=nc>=MIN_CELLS
md=pd.DataFrame({"sample":groups[sel]})
md["donor"]=[s.split("__")[0] for s in md["sample"]]
md["cell_type"]=[s.split("__")[1] for s in md["sample"]]
pc=anno.drop_duplicates("patient").set_index("patient")["condition"]
md["condition"]=[pc[d] for d in md.donor]
md["batch_pair"]="none"; md["group"]=md.condition; md["n_cells"]=nc[sel]
pd.DataFrame(M[:,sel].astype(np.int64),index=expr.index,columns=md["sample"]).to_csv(f"{OUT}/PBC_jakel_matrix.csv")
md.to_csv(f"{OUT}/PBC_jakel_coldata.csv",index=False)
print(f"  wrote PBC_jakel: {M.shape[0]} x {sel.sum()} | donors {md.groupby('condition').donor.nunique().to_dict()}",flush=True)

# ---------------- Beltran: TPM-like normalised ----------------
print("\nBeltran: reading h5ad (normalised, counts NOT recoverable)...",flush=True)
a=ad.read_h5ad("__MS_GEO_ROOT__/SingleCell_CELLxGENE/results/figures/blood_Beltran2019/adata_beltran.h5ad")
Xb=a.raw.X if a.raw is not None else a.X
Xb=Xb.toarray() if hasattr(Xb,"toarray") else np.asarray(Xb)
vn=list(a.raw.var_names) if a.raw is not None else list(a.var_names)
ob=a.obs.copy()
ob["donor"]=ob["twin"].astype(str).str.replace("_PBMCs","",regex=False)
ob["condition"]=np.where(ob["group"].astype(str).str.upper().str.startswith("MS"),"MS","HC")
print(f"  {Xb.shape[0]} cells x {Xb.shape[1]} genes | donors {ob.groupby('condition').donor.nunique().to_dict()}",flush=True)
key=ob["donor"].astype(str)+"__"+ob["celltype"].astype(str)
groups=pd.Index(sorted(key.unique())); gi=pd.Series(np.arange(len(groups)),index=groups)
idx=gi.reindex(key.values).values
Mb=np.zeros((Xb.shape[1],len(groups)))
for j in range(len(groups)):
    sel=idx==j
    if sel.sum(): Mb[:,j]=Xb[sel].mean(0)          # MEAN of normalised values
nc=np.bincount(idx,minlength=len(groups)); sel=nc>=MIN_CELLS
mdb=pd.DataFrame({"sample":groups[sel]})
mdb["donor"]=[s.split("__")[0] for s in mdb["sample"]]
mdb["cell_type"]=[s.split("__")[1] for s in mdb["sample"]]
dc=ob.drop_duplicates("donor").set_index("donor")["condition"]
mdb["condition"]=[dc[d] for d in mdb.donor]
mdb["batch_pair"]="none"; mdb["group"]=mdb.condition; mdb["n_cells"]=nc[sel]
pd.DataFrame(Mb[:,sel],index=vn,columns=mdb["sample"]).to_csv(f"{OUT}/PBN_beltran_meanNorm.csv")
mdb.to_csv(f"{OUT}/PBN_beltran_coldata.csv",index=False)
print(f"  wrote PBN_beltran: {Mb.shape[0]} x {sel.sum()} | donors {mdb.groupby('condition').donor.nunique().to_dict()}")
print(f"  cell types kept: {sorted(mdb.cell_type.unique())}")
