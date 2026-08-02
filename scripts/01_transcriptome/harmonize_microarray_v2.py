#!/usr/bin/env python3
"""Microarray harmonization (track 2) — read the 6 array series_matrix expression tables,
map platform probe IDs -> HGNC symbol (GPL570 via mygene reporter; GPL10558/GPL17586/GPL23126
via the GEO platform SOFT annotation), assign MS/HC from series_matrix metadata."""
import os, re, gzip, io, requests
import pandas as pd, numpy as np
ED="__MS_GEO_ROOT__/Expression_Data"
OUT=f"{ED}/harmonized_v2"; os.makedirs(OUT,exist_ok=True)
MICRO={"GSE21942":"GPL570","GSE43591":"GPL570","GSE38010":"GPL570",
       "GSE103005":"GPL10558","GSE138064":"GPL17586","GSE190847":"GPL23126"}

def sm_table_and_meta(gse):
    """return (expr df probe x GSM, {GSM:'MS'/'HC'})"""
    L=[l.rstrip("\n") for l in gzip.open(f"{ED}/{gse}_series_matrix.txt.gz",'rt',errors='replace')]
    def row(k):
        r=[l for l in L if l.startswith(k)]; return [v.strip('"') for v in r[0].split('\t')[1:]] if r else []
    gsm=row('!Sample_geo_accession'); title=row('!Sample_title')
    dis=None
    for l in L:
        if l.startswith('!Sample_characteristics_ch1'):
            vals=[v.strip('"') for v in l.split('\t')[1:]]
            key=vals[0].split(':')[0].lower() if ':' in vals[0] else ''
            if any(k in key for k in ['disease','status','group','diagnosis','condition','phenotype','health']):
                dis=vals; break
    src=dis if dis else title
    def lab(v):
        v=str(v).lower()
        if any(k in v for k in ['healthy','control','hc','normal','non-ms']) and 'sclerosis' not in v: return 'HC'
        if any(k in v for k in ['multiple sclerosis',' ms','rrms','ppms','spms','cis','progressive','relapsing','patient','ms ']) or v.strip()=='ms': return 'MS'
        return None
    g2d={g:lab(s) for g,s in zip(gsm,src)}
    # expression table
    rows=[]; intab=False
    for l in L:
        if l.startswith('!series_matrix_table_begin'): intab=True; continue
        if l.startswith('!series_matrix_table_end'): break
        if intab: rows.append(l)
    df=pd.read_csv(io.StringIO("\n".join(rows)),sep='\t',index_col=0)
    df.columns=[str(c).strip().strip('"') for c in df.columns]
    return df, g2d

def gpl_soft_map(gpl):
    """probe -> symbol from GEO platform SOFT (full)."""
    url=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gpl}&targ=self&form=text&view=full"
    txt=requests.get(url,timeout=180).text
    lines=txt.splitlines(); intab=False; hdr=None; symi=-1; gene_asg=False; mp={}
    for l in lines:
        if l.startswith('#') or l.startswith('^') or l.startswith('!'):
            if l.startswith('!platform_table_begin'): intab=True
            elif l.startswith('!platform_table_end'): break
            continue
        if not intab: continue
        parts=l.split('\t')
        if hdr is None:
            hdr=[h.strip().lower() for h in parts]
            for cand in ['gene symbol','symbol','ilmn_gene','gene_symbol']:
                if cand in hdr: symi=hdr.index(cand); break
            if symi<0 and 'gene_assignment' in hdr: symi=hdr.index('gene_assignment'); gene_asg=True
            continue
        if symi>=0 and len(parts)>symi:
            pid=parts[0].strip(); raw=parts[symi].strip()
            sym = (raw.split('//')[1].strip() if (gene_asg and '//' in raw) else raw)
            sym = sym.split('///')[0].split('//')[0].strip()
            if sym and sym not in ('---','','NA'): mp[pid]=sym
    return mp

# ---- mygene reporter for GPL570 ----
import mygene; mg=mygene.MyGeneInfo()
print("loading microarray series_matrix tables...")
DAT={}
for g,gpl in MICRO.items():
    df,g2d=sm_table_and_meta(g)
    DAT[g]=(df,gpl,g2d)
    print(f"  {g} ({gpl}): {df.shape[0]} probes x {df.shape[1]} samples")

# build probe->symbol maps per platform
soft_cache={}
def get_map(gpl, probes):
    if gpl=="GPL570":
        r=mg.querymany(list(probes),scopes="reporter",fields="symbol",species="human",verbose=False)
        return {x['query']:x['symbol'] for x in r if 'symbol' in x}
    if gpl not in soft_cache:
        print(f"  fetching GPL SOFT annotation: {gpl} ...")
        soft_cache[gpl]=gpl_soft_map(gpl)
        print(f"    {gpl}: {len(soft_cache[gpl])} probe->symbol")
    return soft_cache[gpl]

meta=[]
for g,(df,gpl,g2d) in DAT.items():
    pm=get_map(gpl, df.index.tolist())
    df.index=[pm.get(str(i)) for i in df.index]
    df=df[[x is not None for x in df.index]]
    df=df.apply(pd.to_numeric,errors='coerce').groupby(level=0).max().dropna(how='all')
    cond={c:g2d.get(c) for c in df.columns}
    keep=[c for c in df.columns if cond[c] in ('MS','HC')]
    df=df[keep]
    df.to_csv(f"{OUT}/{g}_symbol_matrix.csv")
    for c in df.columns: meta.append({"sample_id":f"{g}__{c}","dataset":g,"condition":cond[c]})
    ms=sum(1 for c in df.columns if cond[c]=='MS')
    print(f"  {g} ({gpl}): {df.shape[0]} symbols x {df.shape[1]} ({ms} MS / {df.shape[1]-ms} HC)")
md=pd.DataFrame(meta); md.to_csv(f"{OUT}/microarray_v2_metadata.csv",index=False)
print(f"\nMICROARRAY v2 TOTAL: {len(md)} samples ({(md.condition=='MS').sum()} MS / {(md.condition=='HC').sum()} HC) / {md.dataset.nunique()} datasets")
