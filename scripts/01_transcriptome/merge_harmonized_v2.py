#!/usr/bin/env python3
"""Merge the two harmonization tracks (RNA-seq v3 + microarray v2) into one symbol x sample
matrix on gene symbols, ready for per-dataset normalisation + ComBat (downstream)."""
import os, glob
import pandas as pd, numpy as np
ED="__MS_GEO_ROOT__/Expression_Data"; H=f"{ED}/harmonized_v2"

mats=sorted(glob.glob(f"{H}/*_symbol_matrix.csv"))
frames={};
for f in mats:
    g=os.path.basename(f).replace("_symbol_matrix.csv","")
    df=pd.read_csv(f,index_col=0)
    df=df[~df.index.duplicated(keep='first')]
    df.columns=[f"{g}__{c}" for c in df.columns]
    frames[g]=df
print(f"loaded {len(frames)} datasets: "+", ".join(f"{g}({frames[g].shape[1]})" for g in frames))

# gene density: keep symbols present in >= 50% of datasets
from collections import Counter
cnt=Counter()
for df in frames.values(): cnt.update(set(df.index))
ndat=len(frames); keep=[g for g,c in cnt.items() if c>=0.5*ndat]
print(f"genes in >=50% of {ndat} datasets: {len(keep)}")

merged=pd.concat([df.reindex(keep) for df in frames.values()], axis=1)
merged.index.name="gene"
# metadata
meta=pd.concat([pd.read_csv(f"{H}/rnaseq_v2_metadata.csv"), pd.read_csv(f"{H}/microarray_v2_metadata.csv")], ignore_index=True)
meta=meta[meta.sample_id.isin(merged.columns)]
merged=merged[meta.sample_id.tolist()]
merged.to_csv(f"{H}/Global_Harmonized_v2_Expression.csv")
meta.to_csv(f"{H}/Global_Harmonized_v2_Metadata.csv",index=False)

rs=set(["GSE209596","GSE211739","GSE207680","GSE211358","GSE288904","GSE214334","GSE172009","GSE137143","GSE173789","GSE66573"])
meta["platform"]=meta.dataset.map(lambda d:"RNA-seq" if d in rs else "microarray")
print(f"\n=== GLOBAL HARMONIZED v2 ===")
print(f"  matrix: {merged.shape[0]} genes x {merged.shape[1]} samples")
print(f"  datasets: {meta.dataset.nunique()} | MS={ (meta.condition=='MS').sum() } / HC={ (meta.condition=='HC').sum() }")
print(f"  RNA-seq: {(meta.platform=='RNA-seq').sum()} samples | microarray: {(meta.platform=='microarray').sum()} samples")
print(f"  saved -> Global_Harmonized_v2_Expression.csv + _Metadata.csv")
