#!/usr/bin/env python3
"""RNA-seq harmonization v3 — robust per-dataset gene-ID -> symbol + metadata-based MS/HC.
Loads the 9 RNA-seq series from DEPOSITED supplementary files; assigns condition from
series_matrix metadata where sample codes don't encode it."""
import os, re, gzip, io, tarfile
import pandas as pd, numpy as np
ED="__MS_GEO_ROOT__/Expression_Data"
OUT=f"{ED}/harmonized_v2"; os.makedirs(OUT,exist_ok=True)
# clean re-downloads replacing corrupt RAW.tar members (GSM -> external clean file)
PATCH={"GSE137143":{"GSM4071601":f"{ED}/GSM4071601_56514a-CD14.genes.results.txt.gz"}}

def read_tab(path, sep=None):
    op=gzip.open(path,'rt',errors='replace')
    if sep is None:
        h=op.readline(); sep='\t' if h.count('\t')>=h.count(',') else ','; op.seek(0)
    df=pd.read_csv(op, sep=sep, index_col=0); df.columns=[str(c).strip() for c in df.columns]; return df

def sm_meta(gse):
    """return dict GSM->disease('MS'/'HC'/None) and ordered list of disease by GSM order."""
    p=f"{ED}/{gse}_series_matrix.txt.gz"
    if not os.path.exists(p): return {},[]
    L=[l.rstrip() for l in gzip.open(p,'rt',errors='replace')]
    def row(k):
        r=[l for l in L if l.startswith(k)]; return [v.strip('"') for v in r[0].split('\t')[1:]] if r else []
    gsm=row('!Sample_geo_accession'); title=row('!Sample_title')
    dis=None
    for l in L:
        if l.startswith('!Sample_characteristics_ch1'):
            vals=[v.strip('"') for v in l.split('\t')[1:]]
            key=vals[0].split(':')[0].lower() if ':' in vals[0] else ''
            if any(k in key for k in ['disease state','patient status','disease status','group','diagnosis','condition']):
                dis=vals; break
    src = dis if dis else title
    def lab(v):
        v=str(v).lower()
        if any(k in v for k in ['healthy','control','hc','non-ms','normal']) and 'ms' not in v.replace('healthy',''): return 'HC'
        if any(k in v for k in ['multiple sclerosis','ms','rrms','ppms','spms','cis','progressive','relapsing']): return 'MS'
        return None
    g2d={g:lab(s) for g,s in zip(gsm,src)}
    order=[lab(s) for s in src]
    return g2d, order

def cond_name(s):
    s=str(s).strip().upper()
    if re.search(r'(RRMS|PPMS|SPMS|PATIENT|MULTIPLE|(?:^|_)MS(?:[_\d]|$))', s): return 'MS'
    if re.search(r'((?:^|_)(HC|HD)(?:[_\d]|$)|CONTROL|HEALTHY|CTRL|CTL|NORMAL|(?:^|_)C\d)', s): return 'HC'
    return None

# ---- load gene matrices (idtype) ----
def load_matrix():
    D={}
    D["GSE209596"]=(read_tab(f"{ED}/GSE209596_rsem.filtered.counts.exvivo.memory.v2.txt.gz"),"ensg","name")
    D["GSE211739"]=(read_tab(f"{ED}/GSE211739_all_FPKMs.txt.gz"),"ensg","name")
    d=read_tab(f"{ED}/GSE207680_Human_RNAseq_Data_Processed.csv.gz",sep=","); d.index=[str(i).split('.')[0] for i in d.index]
    D["GSE207680"]=(d,"ensg","positional")
    D["GSE211358"]=(read_tab(f"{ED}/GSE211358_MSHD_mB.txt.gz"),"ensg","gse211358")
    D["GSE288904"]=(read_tab(f"{ED}/GSE288904_gene_TPM_raw.txt.gz"),"symbol","hm_prefix")
    g=read_tab(f"{ED}/GSE214334_countmatrix_MJ.txt.gz")
    for c in ["Chr","Start","End","Strand","Length"]:
        if c in g.columns: g=g.drop(columns=c)
    D["GSE214334"]=(g,"symbol","name")
    D["GSE172009"]=(read_tab(f"{ED}/GSE172009_Raw_gene_counts_matrix.txt.gz"),"refseq","name")
    g66=read_tab(f"{ED}/GSE66573_FPKM_table.txt.gz")
    if "gene_name" in g66.columns: g66=g66.drop(columns="gene_name")  # keep ENSG index
    D["GSE66573"]=(g66,"ensg","gse66573")   # whole blood RNA-seq FPKM (6 RRMS/8 HC), _EST=MS _CTRL=HC
    return D

def load_tar(gse, valuecol, idsplit, key_gsm):
    import zlib
    t=tarfile.open(f"{ED}/{gse}_RAW.tar"); cols={}; bad=0; patched=0
    for m in t.getmembers():
        if not m.name.endswith(('.gz','.txt')): continue
        gsm=m.name.split('_')[0]
        clean=re.sub(r'^GSM\d+_','',m.name).replace('.genes.results.txt.gz','').replace('_cpm.txt.gz','').replace('.txt.gz','')
        key = gsm if key_gsm else clean
        txt=None
        try:
            raw=t.extractfile(m).read()
            txt = gzip.decompress(raw).decode('utf-8','replace') if m.name.endswith('.gz') else raw.decode('utf-8','replace')
        except Exception:
            pf=PATCH.get(gse,{}).get(gsm)            # clean re-download for known-corrupt member
            if pf and os.path.exists(pf):
                txt=gzip.decompress(open(pf,'rb').read()).decode('utf-8','replace'); patched+=1
            else:
                bad+=1; continue
        try:
            d=pd.read_csv(io.StringIO(txt),sep='\t',index_col=0); d=d[~d.index.duplicated(keep='first')]
        except Exception:
            bad+=1; continue
        val=d[valuecol] if (valuecol and valuecol in d.columns) else d.iloc[:,0]
        cols[key]=val
    if bad: print(f"   ({gse}: skipped {bad} unreadable)")
    if patched: print(f"   ({gse}: used {patched} clean re-download(s) for corrupt members)")
    df=pd.DataFrame(cols)
    if idsplit: df.index=[str(i).split('_',1)[1] if '_' in str(i) else str(i) for i in df.index]
    return df

print("loading RNA-seq matrices + RAW.tar ...")
D=load_matrix()
D["GSE137143"]=(load_tar("GSE137143","expected_count",True,True),"symbol","gsm137143")
D["GSE173789"]=(load_tar("GSE173789",None,False,False),"ensg","name")
for g,(df,t,_) in D.items():
    df.index=[str(i).split('.')[0] for i in df.index] if t=="ensg" else [str(i) for i in df.index]

# ---- map ENSG/RefSeq -> symbol ----
import mygene; mg=mygene.MyGeneInfo()
ensg=set(); refs=set()
for g,(df,t,_) in D.items():
    if t=="ensg": ensg|=set(df.index)
    elif t=="refseq": refs|=set(str(i).split('.')[0] for i in df.index)
print(f"mapping {len(ensg)} ENSG + {len(refs)} RefSeq...")
e2s={r['query']:r['symbol'] for r in mg.querymany(list(ensg),scopes='ensembl.gene',fields='symbol',species='human',verbose=False) if 'symbol' in r}
r2s={r['query']:r['symbol'] for r in mg.querymany(list(refs),scopes='refseq',fields='symbol',species='human',verbose=False) if 'symbol' in r} if refs else {}

# ---- per-dataset condition + collapse + save ----
meta=[]
for g,(df,t,mode) in D.items():
    if t=="ensg": df.index=[e2s.get(i) for i in df.index]
    elif t=="refseq": df.index=[r2s.get(str(i).split('.')[0]) for i in df.index]
    df=df[[x is not None for x in df.index]]
    df=df.apply(pd.to_numeric,errors='coerce').groupby(level=0).max().dropna(how='all')
    # condition per sample-column
    cond={}
    if mode=="gsm137143" or mode=="positional":
        g2d,order=sm_meta(g)
        if mode=="gsm137143":
            cond={c:g2d.get(c) for c in df.columns}          # columns are GSM
        else:  # positional: assume column order == GSM order
            cond={c:(order[i] if i<len(order) else None) for i,c in enumerate(df.columns)}
    elif mode=="hm_prefix":
        cond={c:('MS' if str(c).upper().startswith('M') else 'HC' if str(c).upper().startswith('H') else None) for c in df.columns}
    elif mode=="gse211358":
        cond={c:('MS' if str(c).upper().startswith('MS') else 'HC' if str(c).upper().startswith('HD') else None) for c in df.columns}  # HS excluded
    elif mode=="gse66573":
        cond={c:('MS' if str(c).upper().endswith('_EST') else 'HC' if str(c).upper().endswith('_CTRL') else None) for c in df.columns}  # _EST=RRMS _CTRL=HC
    else:
        cond={c:cond_name(c) for c in df.columns}
    keep=[c for c in df.columns if cond[c] in ('MS','HC')]
    df=df[keep]
    # rename GSE137143 GSM columns to readable
    df.to_csv(f"{OUT}/{g}_symbol_matrix.csv")
    for c in df.columns: meta.append({"sample_id":f"{g}__{c}","dataset":g,"condition":cond[c]})
    ms=sum(1 for c in df.columns if cond[c]=='MS'); print(f"  {g}: {df.shape[0]} symbols x {df.shape[1]} ({ms} MS / {df.shape[1]-ms} HC)")
md=pd.DataFrame(meta); md.to_csv(f"{OUT}/rnaseq_v2_metadata.csv",index=False)
print(f"\nRNA-seq v3 TOTAL: {len(md)} samples ({(md.condition=='MS').sum()} MS / {(md.condition=='HC').sum()} HC) / {md.dataset.nunique()} datasets")
