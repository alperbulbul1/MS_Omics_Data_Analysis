import os
import pandas as pd
import gzip
import GEOparse
import numpy as np
import mygene

dest_dir = "__MS_GEO_ROOT__/Expression_Data"

def get_ensg_to_symbol_mapping(ensg_list):
    mg = mygene.MyGeneInfo()
    print("Querying mygene for symbols...")
    results = mg.querymany(ensg_list, scopes='ensembl.gene', fields='symbol', species='human', as_dataframe=True)
    mapping = {}
    for ensg, row in results.iterrows():
        if isinstance(row, pd.Series):
            symbol = row.get('symbol')
        else: # dataframe with multiple matching
            symbol = row.iloc[0].get('symbol') if 'symbol' in row.columns else None
            
        if pd.notna(symbol):
            mapping[ensg] = symbol
    return mapping

def load_gse173789():
    print("Parsing GSE173789 (RNA-Seq)...")
    dir_path = os.path.join(dest_dir, "GSE173789")
    files = [f for f in os.listdir(dir_path) if f.endswith(".txt.gz")]
    dfs = []
    meta = []
    
    for f in files:
        sample_id = f.split("_")[0]
        condition = "MS" if "_MS_" in f else "HC"
        
        df = pd.read_csv(os.path.join(dir_path, f), sep="\t")
        df = df.rename(columns={df.columns[1]: sample_id})
        df = df.set_index("Gene_id")
        dfs.append(df)
        
        meta.append({"sample_id": sample_id, "dataset": "GSE173789", "condition": condition})
        
    expr = pd.concat(dfs, axis=1)
    # Log2 transform to emulate normal distribution for ComBat (CPMs)
    expr = np.log2(expr + 1)
    
    # Map ENSG to SYMBOL
    mapping = get_ensg_to_symbol_mapping(expr.index.tolist())
    expr = expr.rename(index=mapping)
    # Filter out anything that wasn't mapped (still starts with ENSG)
    expr = expr[~expr.index.str.startswith('ENSG')]
    # Aggregate by symbol
    expr = expr.groupby(expr.index).mean()
    
    return expr, pd.DataFrame(meta)

def load_gse235357():
    print("Parsing GSE235357 (RNA-Seq)...")
    csv_path = os.path.join(dest_dir, "GSE235357_normalized_annotated.csv.gz")
    df = pd.read_csv(csv_path)
    df = df.dropna(subset=['SYMBOL'])
    df = df.groupby('SYMBOL').mean(numeric_only=True).reset_index()
    df = df.set_index('SYMBOL')
    if "Unnamed: 0" in df.columns:
        df = df.drop(columns=["Unnamed: 0"])
        
    expr = np.log2(df + 1)
    
    # Needs metadata mapping
    gse = GEOparse.get_GEO(filepath=os.path.join(dest_dir, "GSE235357_family.soft.gz"), silent=True)
    meta = []
    
    # We will map based on numerical order since SM002604_1 corresponds to the first GSM.
    # GEO samples are dicts. The order in .gsms.values() is usually the intended order.
    gsms = list(gse.gsms.values())
    
    for i, col in enumerate(expr.columns):
        if i < len(gsms):
            title = gsms[i].metadata.get("title", [""])[0]
            cond = "HC" if "HC" in title.upper() or "HEALTHY" in title.upper() else "MS"
            meta.append({"sample_id": col, "dataset": "GSE235357", "condition": cond})
        else:
            meta.append({"sample_id": col, "dataset": "GSE235357", "condition": "Unknown"})
            
    return expr, pd.DataFrame(meta)

def load_gse137143():
    print("Parsing GSE137143 (RNA-Seq)...")
    dir_path = os.path.join(dest_dir, "GSE137143")
    if not os.path.exists(dir_path):
        return None, None
        
    gse = GEOparse.get_GEO(filepath=os.path.join(dest_dir, "GSE137143_family.soft.gz"), silent=True)
    files = [f for f in os.listdir(dir_path) if f.endswith(".txt.gz")]
    dfs = []
    meta = []
    
    for f in files:
        gsm_id = f.split("_")[0]
        if gsm_id not in gse.gsms:
            continue
            
        gsm = gse.gsms[gsm_id]
        chars = " ".join(gsm.metadata.get("characteristics_ch1", [])).lower()
        
        cond = "Unknown"
        if "healthy control" in chars:
            cond = "HC"
        elif "multiple sclerosis" in chars:
            cond = "MS"
            
        if cond == "Unknown":
            continue
            
        df = pd.read_csv(os.path.join(dir_path, f), sep="\t")
        df['SYMBOL'] = df['gene_id'].apply(lambda x: x.split('_')[1] if isinstance(x, str) and '_' in x else None)
        df = df.dropna(subset=['SYMBOL'])
        
        df = df[['SYMBOL', 'FPKM']].rename(columns={'FPKM': gsm_id})
        df = df.groupby('SYMBOL').mean(numeric_only=True)
        dfs.append(df)
        
        meta.append({"sample_id": gsm_id, "dataset": "GSE137143", "condition": cond})
        
    if not dfs:
        return None, None
        
    expr = pd.concat(dfs, axis=1)
    expr = np.log2(expr + 1)
    return expr, pd.DataFrame(meta)

def load_microarray(gse_id, platform_id, symbol_col):
    print(f"Parsing {gse_id} (Microarray)...")
    # Due to size, we fetch the series matrix bypassing massive soft if possible, or just use the local soft
    soft_path = os.path.join(dest_dir, f"{gse_id}_family.soft.gz")
    if not os.path.exists(soft_path):
        return None, None
        
    try:
        gse = GEOparse.get_GEO(filepath=soft_path, silent=True)
    except Exception as e:
        print(e)
        return None, None
        
        
    matrix_path = os.path.join(dest_dir, f"{gse_id}_series_matrix.txt.gz")
    if os.path.exists(matrix_path):
        expr = pd.read_csv(matrix_path, sep="\t", comment="!", index_col=0)
    else:
        print(f"Matrix not found: {matrix_path}")
        return None, None
    if platform_id not in gse.gpls:
        platform_id = list(gse.gpls.keys())[0] # Just use first if mismatch
        
    gpl = gse.gpls[platform_id]
    
    # Try mapping
    # Symbol column can be 'Gene Symbol', 'Symbol', 'gene_assignment' (clariom D)
    mapping = {}
    if symbol_col in gpl.table.columns:
        for idx, row in gpl.table.iterrows():
            sym = str(row[symbol_col])
            # clariom d format: "ENST0000.. // symbol // ..."
            if '//' in sym:
                parts = sym.split('//')
                if len(parts) > 1:
                    sym = parts[1].strip()
            if sym and sym != 'nan' and sym != '---':
                mapping[row['ID']] = sym
                
    expr = expr.rename(index=mapping)
    # Keep only mapped
    expr = expr[expr.index.isin(mapping.values())]
    # Aggregate
    expr = expr.groupby(expr.index).mean()
    
    meta = []
    for gsm_name, gsm in gse.gsms.items():
        # Condition identification
        chars = " ".join(gsm.metadata.get("characteristics_ch1", [])).upper()
        title = gsm.metadata.get("title", [""])[0].upper()
        
        cond = "Unknown"
        if "HEALTHY" in chars or "HEALTHY" in title or "NORMAL" in chars or " HC" in title or "CONTROL" in chars:
            cond = "HC"
        elif "MULTIPLE SCLEROSIS" in chars or " MS " in title or " MS" in title or "RELAPSING" in chars or "DISEASE: MS" in chars:
            cond = "MS"
            
        meta.append({"sample_id": gsm_name, "dataset": gse_id, "condition": cond})
        
    return expr, pd.DataFrame(meta)

def main():
    print("Starting integration pipeline...")
    expr1, meta1 = load_gse173789()
    expr2, meta2 = load_gse235357()
    expr3, meta3 = load_gse137143()
    
    # GSE255952: GPL23126 (Clariom D)
    expr4, meta4 = load_microarray("GSE255952", "GPL23126", "gene_assignment")
    
    # Identify intersection of genes
    genes = set(expr1.index) & set(expr2.index)
    if expr3 is not None:
        genes = genes & set(expr3.index)
    if expr4 is not None:
        genes = genes & set(expr4.index)
        
    genes = list(genes)
    print(f"Intersected Genes across datasets: {len(genes)}")
    
    if len(genes) < 100:
        print("Warning: Very few intersecting genes.")
        
    # Slice expressions
    expr1 = expr1.loc[genes]
    expr2 = expr2.loc[genes]
    exprs = [expr1, expr2]
    metas = [meta1, meta2]
    
    if expr3 is not None:
        exprs.append(expr3.loc[genes])
        metas.append(meta3)
        
    if expr4 is not None:
        exprs.append(expr4.loc[genes])
        metas.append(meta4)
        
        
    # Concat
    final_expr = pd.concat(exprs, axis=1)
    final_meta = pd.concat(metas, ignore_index=True)
    
    # Save raw harmonized
    final_expr.to_csv(os.path.join(dest_dir, "Harmonized_Expression_Matrix.csv"))
    final_meta.to_csv(os.path.join(dest_dir, "Harmonized_Metadata.csv"), index=False)
    
    print("Saved Harmonized matrices safely.")

if __name__ == '__main__':
    main()
