import os
import gzip
import pandas as pd
import numpy as np
import GEOparse
import json
import mygene

dest_dir = "__MS_GEO_ROOT__/Expression_Data"

def fast_parse_metadata(soft_path):
    # Instead of reading the entire multi-GB array table, we just parse the soft file headers to get MS vs HC metadata
    metadata = {}
    current_gsm = None
    try:
        with gzip.open(soft_path, 'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('^SAMPLE = '):
                    current_gsm = line.split('=')[1].strip()
                    metadata[current_gsm] = []
                elif current_gsm and (line.startswith('!Sample_characteristics_ch1 = ') or line.startswith('!Sample_title = ')):
                    metadata[current_gsm].append(line.split('=', 1)[1].strip().lower())
                elif line.startswith('^SERIES = '):
                    # Done with samples
                    pass
    except Exception as e:
        print(f"Error parsing soft headers: {e}")
    
    return metadata

def fast_parse_gpl(soft_path):
    # Extract only the ^PLATFORM table to map probes quickly
    mapping = {}
    in_platform_table = False
    columns = []
    symbol_idx = -1
    id_idx = -1
    try:
        with gzip.open(soft_path, 'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('!platform_table_begin'):
                    in_platform_table = True
                    continue
                elif line.startswith('!platform_table_end'):
                    break
                    
                if in_platform_table:
                    parts = line.split('\t')
                    if not columns:
                        columns = [c.lower() for c in parts]
                        # Find ID and Symbol
                        if 'id' in columns:
                            id_idx = columns.index('id')
                        for candidate in ['gene symbol', 'symbol', 'gene_assignment', 'ilmn_gene']:
                            if candidate in columns:
                                symbol_idx = columns.index(candidate)
                                break
                        if id_idx == -1 or symbol_idx == -1:
                            return {} # Cannot map this platform
                        continue
                        
                    if len(parts) > max(id_idx, symbol_idx):
                        probe_id = parts[id_idx]
                        sym = parts[symbol_idx]
                        if '//' in sym:
                            sym = sym.split('//')[1].strip()
                        if sym and sym != 'nan' and sym != '---':
                            mapping[probe_id] = sym
    except Exception as e:
        pass
    return mapping

def global_harmonize():
    df = pd.read_csv('__MS_GEO_ROOT__/Refined_MS_HC_Global_Datasets.csv')
    gse_ids = df['gse_id'].tolist()
    
    # We will exclude massive single cell repositories if they slipped through
    drop_list = ["GSE118257", "GSE123496"]
    gse_ids = [g for g in gse_ids if g not in drop_list]
    print(f"Loaded {len(gse_ids)} valid targets for fast global harmonization.")
    
    matrices = {}
    master_metadata = []
    
    for gse in gse_ids:
        print(f"\n--- {gse} ---")
        gse_dir = os.path.join(dest_dir, gse)
        soft_path = os.path.join(gse_dir, f"{gse}_family.soft.gz")
        matrix_path = os.path.join(gse_dir, f"{gse}_series_matrix.txt.gz")
        
        if not os.path.exists(soft_path) or not os.path.exists(matrix_path):
            print(f"Skipping {gse}: Missing files.")
            continue
            
        try:
            # Fast parse metadata
            raw_meta = fast_parse_metadata(soft_path)
            if not raw_meta:
                print(f"No metadata found in soft headers.")
                continue
                
            meta_clean = []
            for gsm_id, lines in raw_meta.items():
                text = " ".join(lines)
                cond = "Unknown"
                if "healthy" in text or "control" in text or " hc " in text:
                    cond = "HC"
                elif "multiple sclerosis" in text or " ms " in text or "relapsing" in text:
                    cond = "MS"
                if cond != "Unknown":
                    meta_clean.append({"sample_id": gsm_id, "dataset": gse, "condition": cond})
            
            if not meta_clean:
                print(f"Skipping {gse}: Cohort not exclusively MS/HC.")
                continue
                
            meta_df = pd.DataFrame(meta_clean)
            
            # Fast load series matrix
            expr = pd.read_csv(matrix_path, sep='\t', comment='!', index_col=0, low_memory=False)
            # Filter samples
            valid_samples = [s for s in meta_df['sample_id'] if s in expr.columns]
            meta_df = meta_df[meta_df['sample_id'].isin(valid_samples)]
            expr = expr[valid_samples]
            
            if expr.empty:
                print("Skipping: Expression matrix empty post-filtering.")
                continue
                
            mapped = False
            # 1. Try GPL Probe Mapping
            probe_map = fast_parse_gpl(soft_path)
            if probe_map:
                expr = expr.rename(index=probe_map)
                expr = expr[expr.index.isin(probe_map.values())]
                expr = expr.groupby(expr.index).mean()
                mapped = True
            
            if not mapped:
                # 2. Heuristic for RNA-Seq ENSG
                if any(str(idx).startswith('ENSG') for idx in expr.index[:20]):
                    mg = mygene.MyGeneInfo()
                    results = mg.querymany(expr.index.tolist(), scopes='ensembl.gene', fields='symbol', species='human', as_dataframe=True, verbose=False)
                    mapping = {}
                    for ensg, row in results.iterrows():
                        if isinstance(row, pd.Series):
                            symbol = row.get('symbol')
                        else:
                            symbol = row.iloc[0].get('symbol') if 'symbol' in row.columns else None
                        if pd.notna(symbol):
                            mapping[ensg] = symbol
                            
                    expr = expr.rename(index=mapping)
                    expr = expr[~expr.index.str.startswith('ENSG')]
                    expr = expr.groupby(expr.index).mean()
                    
                    mapped = True
                    
            if mapped and not expr.empty:
                matrices[gse] = expr
                master_metadata.extend(meta_df.to_dict('records'))
                print(f"Success! Integrated {expr.shape[1]} samples over {expr.shape[0]} genes.")
            else:
                print(f"Skipping {gse}: ID Mapping failed.")
                
        except Exception as e:
            print(f"Error parsing {gse}: {e}")
            
    print(f"\n============================")
    print(f"Harmonization complete. Loaded {len(matrices)} valid datasets.")
    
    if matrices:
        common_genes = list(matrices.values())[0].index
        for gse, expr in matrices.items():
            common_genes = common_genes.intersection(expr.index)
            
        print(f"Strict absolute consensus: {len(common_genes)} universal genes")
        
        # Dynamic 75% scaling
        if len(common_genes) < 8000:
            print("Consensus strictly too low! Applying dynamic 75% dataset-density thresholding imputator...")
            gene_counts = pd.Series(dtype=int)
            for expr in matrices.values():
                gene_counts = gene_counts.add(pd.Series(1, index=expr.index), fill_value=0)
                
            threshold = int(len(matrices) * 0.75)
            common_genes = gene_counts[gene_counts >= threshold].index
            print(f"Recovered robust {len(common_genes)} genes present in >= 75% of clinical trials.")
            
        final_exprs = []
        for gse, expr in matrices.items():
            expr = expr.reindex(common_genes).fillna(0)
            final_exprs.append(expr)
            
        final_matrix = pd.concat(final_exprs, axis=1)
        final_meta_df = pd.DataFrame(master_metadata)
        
        final_matrix.to_csv(os.path.join(dest_dir, "Global_Harmonized_Expression.csv"))
        final_meta_df.to_csv(os.path.join(dest_dir, "Global_Harmonized_Metadata.csv"), index=False)
        print(f"\nDone! Exported Global Meta-Matrix: {final_matrix.shape[1]} patients tracking {final_matrix.shape[0]} genes.")

if __name__ == '__main__':
    global_harmonize()
