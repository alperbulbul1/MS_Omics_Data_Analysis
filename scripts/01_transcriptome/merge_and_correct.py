"""
Merge the recovered supplementary datasets with the original 6-dataset global matrix,
then re-run neuroCombat and Limma DGE on the combined super-matrix.
"""
import os
import pandas as pd
import numpy as np
from neuroCombat import neuroCombat
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import subprocess

dest_dir = "__MS_GEO_ROOT__/Expression_Data"

def merge_and_run():
    # Load original global matrix
    global_expr = pd.read_csv(os.path.join(dest_dir, "Global_Harmonized_Expression.csv"), index_col=0)
    global_meta = pd.read_csv(os.path.join(dest_dir, "Global_Harmonized_Metadata.csv"))
    
    # Load newly recovered matrix
    recovered_expr = pd.read_csv(os.path.join(dest_dir, "Recovered_Expression.csv"), index_col=0)
    recovered_meta = pd.read_csv(os.path.join(dest_dir, "Recovered_Metadata.csv"))
    
    print(f"Original global: {global_expr.shape[1]} samples, {global_expr.shape[0]} genes")
    print(f"Recovered:       {recovered_expr.shape[1]} samples, {recovered_expr.shape[0]} genes")
    
    # Find common genes
    common_genes = global_expr.index.intersection(recovered_expr.index)
    print(f"Common genes:    {len(common_genes)}")
    
    global_expr = global_expr.reindex(common_genes)
    recovered_expr = recovered_expr.reindex(common_genes)
    
    # Combine
    combined_expr = pd.concat([global_expr, recovered_expr], axis=1)
    combined_meta = pd.concat([global_meta, recovered_meta], ignore_index=True)
    
    # Deduplicate: keep only the first occurrence of each sample ID
    combined_expr = combined_expr.loc[:, ~combined_expr.columns.duplicated(keep='first')]
    combined_meta = combined_meta.drop_duplicates(subset='sample_id', keep='first')
    
    # Filter to MS/HC only and match
    combined_meta = combined_meta[combined_meta['condition'].isin(["MS", "HC"])]
    valid_samples = [s for s in combined_meta['sample_id'] if s in combined_expr.columns]
    combined_meta = combined_meta[combined_meta['sample_id'].isin(valid_samples)]
    combined_expr = combined_expr[valid_samples]
    
    # Drop zero-variance
    combined_expr = combined_expr.loc[combined_expr.var(axis=1) > 0]
    
    print(f"\nCombined: {combined_expr.shape[1]} samples x {combined_expr.shape[0]} genes across {combined_meta['dataset'].nunique()} datasets")
    
    # Save combined pre-correction
    combined_expr.to_csv(os.path.join(dest_dir, "Combined_Expression_Pre_ComBat.csv"))
    combined_meta.to_csv(os.path.join(dest_dir, "Combined_Metadata.csv"), index=False)
    
    # neuroCombat
    print("\nRunning neuroCombat batch correction...")
    combat_data = neuroCombat(
        dat=combined_expr,
        covars=combined_meta,
        batch_col="dataset",
        categorical_cols=["condition"]
    )
    corrected = pd.DataFrame(combat_data["data"], index=combined_expr.index, columns=combined_expr.columns)
    corrected.to_csv(os.path.join(dest_dir, "Combined_Batch_Corrected.csv"))
    print("Batch correction done.")
    
    # Save metadata for R
    combined_meta.to_csv(os.path.join(dest_dir, "Combined_Metadata.csv"), index=False)
    
    print(f"\nFinal matrix: {corrected.shape}")
    return corrected, combined_meta

if __name__ == '__main__':
    corrected, meta = merge_and_run()
    print("\nReady for Limma!")
