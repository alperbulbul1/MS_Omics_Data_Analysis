"""
correct_and_normalize.py
========================
Corrected combined MS transcriptomics normalization pipeline.

Root causes of inflated DGE results:
  1. Raw count datasets (max ~ 50,000–75,000) NOT properly log2-transformed
  2. neuroCombat applied to mixed-scale data → LFC ±3000
  3. No within-dataset quantile normalization before cross-platform merge

This script:
  1. Loads each already-integrated dataset from the Global_Harmonized_Expression + Recovered_Expression matrices
  2. Detects scale per dataset using robust percentile statistics
  3. For raw-count datasets: applies CPM → log2(CPM+1)
  4. For already-log2 microarray datasets: validates range (0–25)
  5. Applies within-dataset quantile normalization to put all datasets on same scale
  6. Merges into a combined matrix
  7. Runs neuroCombat (only on properly validated log2 data)
  8. Saves corrected matrix + PCA plots for validation
"""

import os
import sys
import numpy as np
import pandas as pd
from neuroCombat import neuroCombat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

dest_dir = "__MS_GEO_ROOT__/Expression_Data"
out_dir = "__MS_GEO_ROOT__"

# ─────────────────────────────────────────────────────────────
# 1.  Scale detection
# ─────────────────────────────────────────────────────────────

def detect_scale(df_sub):
    """
    Returns 'raw_counts', 'log2', or 'microarray_log2'.

    Logic:
    - If >10% of values exceed 100  →  raw (counts or unnormalized intensities)
    - If values are tightly bounded in [0, 25] or even negative  →  log2
    - We use the 99th percentile to avoid outlier sensitivity
    """
    vals = df_sub.values.astype(float).flatten()
    vals = vals[np.isfinite(vals)]
    p99 = np.percentile(vals, 99)
    p50 = np.percentile(vals, 50)
    frac_above_100 = np.mean(vals > 100)

    if frac_above_100 > 0.05 or p99 > 200:
        return "raw_counts"
    elif p50 < 0 or (p50 < 2 and p99 < 15):
        return "log2"  # likely quantile-normalized microarray
    else:
        return "log2"  # normal log2 intensity or logCPM


# ─────────────────────────────────────────────────────────────
# 2.  Per-dataset normalization
# ─────────────────────────────────────────────────────────────

def normalize_dataset(df_sub, dataset_id):
    """
    Bring a single dataset's expression sub-matrix to log2 scale,
    then quantile-normalize across samples within the dataset.

    df_sub: genes × samples DataFrame (already subset to this dataset's columns)
    Returns: normalized genes × samples DataFrame
    """
    # Ensure all numeric
    df_sub = df_sub.apply(pd.to_numeric, errors='coerce')
    df_sub = df_sub.loc[df_sub.isnull().mean(axis=1) < 0.5]  # drop > 50% NA rows

    scale = detect_scale(df_sub)
    print(f"  [{dataset_id}] Detected scale: {scale}  (99pct={np.percentile(df_sub.values[np.isfinite(df_sub.values)], 99):.1f})")

    if scale == "raw_counts":
        # 2a. Column-wise CPM (counts per million) then log2
        col_sums = df_sub.sum(axis=0)
        # Guard against zero-sum columns
        col_sums = col_sums.replace(0, 1)
        cpm = df_sub.divide(col_sums, axis=1) * 1e6
        df_sub = np.log2(cpm + 1)
        print(f"       → Applied CPM + log2(CPM+1). New 99pct: {np.percentile(df_sub.values[np.isfinite(df_sub.values)], 99):.1f}")
    else:
        # Already log2 — clip extreme negatives to 0 (some datasets have small negatives)
        df_sub = df_sub.clip(lower=0)

    # 2b. Within-dataset quantile normalization across samples
    # (normalizes each sample to the same distribution — standard for microarray integration)
    df_sub = df_sub.fillna(0.0)

    if df_sub.shape[1] >= 2:  # need at least 2 samples
        # Sort each column independently
        sorted_vals = np.sort(df_sub.values, axis=0)
        # Compute row means (the average distribution)
        mean_dist = np.mean(sorted_vals, axis=1)
        
        # Rank the original dataframe
        ranks = df_sub.rank(method='min').astype(int) - 1
        
        # Map the ranks to the mean distribution
        def map_ranks(col):
            return mean_dist[col]
            
        df_norm = pd.DataFrame({col: map_ranks(ranks[col]) for col in ranks.columns}, index=df_sub.index)
    else:
        df_norm = df_sub

    return df_norm


# ─────────────────────────────────────────────────────────────
# 3.  Main pipeline
# ─────────────────────────────────────────────────────────────

def run():
    print("=" * 60)
    print("STEP 1: Loading combined expression data")
    print("=" * 60)

    # Load the pre-combat combined matrix and metadata
    combined_expr = pd.read_csv(os.path.join(dest_dir, "Combined_Expression_Pre_ComBat.csv"), index_col=0)
    combined_meta = pd.read_csv(os.path.join(dest_dir, "Combined_Metadata.csv"))

    # Filter only MS / HC
    combined_meta = combined_meta[combined_meta['condition'].isin(["MS", "HC"])]
    valid_samples = [s for s in combined_meta['sample_id'] if s in combined_expr.columns]
    combined_meta = combined_meta[combined_meta['sample_id'].isin(valid_samples)].copy()
    combined_meta = combined_meta.drop_duplicates(subset='sample_id', keep='first')
    combined_expr = combined_expr[combined_meta['sample_id'].tolist()]

    print(f"Loaded: {combined_expr.shape[1]} samples × {combined_expr.shape[0]} genes")
    print(f"Datasets: {combined_meta['dataset'].unique().tolist()}")
    print(f"Conditions: {combined_meta['condition'].value_counts().to_dict()}")

    # ─── Step 2: Per-dataset normalization ───────────────────
    print("\n" + "=" * 60)
    print("STEP 2: Per-dataset scale detection & normalization")
    print("=" * 60)

    normalized_datasets = {}
    for ds in combined_meta['dataset'].unique():
        ds_samples = combined_meta[combined_meta['dataset'] == ds]['sample_id'].tolist()
        ds_expr = combined_expr[ds_samples].copy()
        print(f"\n  Processing {ds}: {len(ds_samples)} samples × {ds_expr.shape[0]} genes")
        ds_norm = normalize_dataset(ds_expr, ds)
        normalized_datasets[ds] = ds_norm
        print(f"  → Normalized range: [{ds_norm.values.min():.2f}, {ds_norm.values.max():.2f}]  median={np.nanmedian(ds_norm.values):.2f}")

    # ─── Step 3: Merge normalized datasets ───────────────────
    print("\n" + "=" * 60)
    print("STEP 3: Merging normalized datasets")
    print("=" * 60)

    # Find genes present in at least 75% of datasets
    all_genes = pd.Index([])
    for df in normalized_datasets.values():
        all_genes = all_genes.union(df.index)

    gene_counts = pd.Series(0, index=all_genes)
    for df in normalized_datasets.values():
        gene_counts[df.index] += 1

    n_ds = len(normalized_datasets)
    threshold = max(1, int(n_ds * 0.75))
    common_genes = gene_counts[gene_counts >= threshold].index
    print(f"Genes in ≥75% of datasets: {len(common_genes)}")

    # Reindex each dataset to common genes (fill missing with 0)
    dfs_aligned = []
    for ds, df in normalized_datasets.items():
        df_aligned = df.reindex(common_genes, fill_value=0)
        dfs_aligned.append(df_aligned)

    merged_expr = pd.concat(dfs_aligned, axis=1)
    # Align metadata to merged columns
    merged_expr = merged_expr[combined_meta['sample_id'].tolist()]

    print(f"Merged matrix: {merged_expr.shape[1]} samples × {merged_expr.shape[0]} genes")

    # Drop zero-variance genes (can't batch-correct them)
    var_mask = merged_expr.var(axis=1) > 1e-6
    merged_expr = merged_expr.loc[var_mask]
    print(f"After removing zero-variance genes: {merged_expr.shape[0]} genes remaining")

    # Final range check — must be in reasonable log2 range
    global_max = merged_expr.values.max()
    global_min = merged_expr.values.min()
    print(f"Pre-ComBat range: [{global_min:.2f}, {global_max:.2f}]")
    if global_max > 50:
        print("WARNING: Values still > 50 detected. Clipping extreme outliers at 99.9th percentile.")
        p999 = np.percentile(merged_expr.values, 99.9)
        merged_expr = merged_expr.clip(upper=p999)

    # ─── Step 4: PCA before ComBat ───────────────────────────
    print("\n" + "=" * 60)
    print("STEP 4: PCA before ComBat")
    print("=" * 60)
    _plot_pca(merged_expr, combined_meta, "Corrected_PCA_Before_ComBat.png", "Before Batch Correction (Per-Dataset QN)")

    # Save pre-combat matrix
    merged_expr.to_csv(os.path.join(dest_dir, "Corrected_Expression_Pre_ComBat.csv"))
    combined_meta.to_csv(os.path.join(dest_dir, "Corrected_Metadata.csv"), index=False)

    # ─── Step 5: neuroCombat batch correction ────────────────
    print("\n" + "=" * 60)
    print("STEP 5: neuroCombat batch correction")
    print("=" * 60)

    single_cond = []  # track excluded datasets

    # neuroCombat requires metadata rows aligned to expression columns
    # Use explicit column-order alignment (avoids set_index/reset_index pandas bugs)
    expr_cols = merged_expr.columns.tolist()
    meta_map = combined_meta[combined_meta['sample_id'].isin(expr_cols)].copy()
    meta_map = meta_map.drop_duplicates(subset='sample_id', keep='first')
    # Sort meta to match column order of merged_expr
    col_order = {s: i for i, s in enumerate(expr_cols)}
    meta_aligned = meta_map.iloc[meta_map['sample_id'].map(col_order).argsort()].reset_index(drop=True)

    # Check each batch has both conditions — neuroCombat can fail otherwise
    batch_check = meta_aligned.groupby('dataset')['condition'].nunique()
    single_cond = batch_check[batch_check < 2].index.tolist()
    if single_cond:
        print(f"WARNING: Datasets with only one condition: {single_cond}")
        print("  These datasets will be EXCLUDED from ComBat to avoid rank-deficient design.")
        meta_aligned = meta_aligned[~meta_aligned['dataset'].isin(single_cond)].reset_index(drop=True)
        keep_samples = meta_aligned['sample_id'].tolist()
        merged_expr = merged_expr[[s for s in expr_cols if s in keep_samples]]
        var_mask2 = merged_expr.var(axis=1) > 1e-6
        merged_expr = merged_expr.loc[var_mask2]
        # Re-sort meta to match new column order
        new_cols = merged_expr.columns.tolist()
        col_order2 = {s: i for i, s in enumerate(new_cols)}
        meta_aligned = meta_aligned[meta_aligned['sample_id'].isin(new_cols)].copy()
        meta_aligned = meta_aligned.iloc[meta_aligned['sample_id'].map(col_order2).argsort()].reset_index(drop=True)
        print(f"  After exclusion: {merged_expr.shape[1]} samples × {merged_expr.shape[0]} genes")

    # Also check batches with too few samples
    batch_sizes = meta_aligned['dataset'].value_counts()
    small_batches = batch_sizes[batch_sizes < 2].index.tolist()
    if small_batches:
        print(f"WARNING: Datasets with <2 samples: {small_batches} — excluding from ComBat")
        meta_aligned = meta_aligned[~meta_aligned['dataset'].isin(small_batches)].reset_index(drop=True)
        keep_samples2 = meta_aligned['sample_id'].tolist()
        merged_expr = merged_expr[[s for s in merged_expr.columns if s in keep_samples2]]
        var_mask3 = merged_expr.var(axis=1) > 1e-6
        merged_expr = merged_expr.loc[var_mask3]
        new_cols2 = merged_expr.columns.tolist()
        col_order3 = {s: i for i, s in enumerate(new_cols2)}
        meta_aligned = meta_aligned[meta_aligned['sample_id'].isin(new_cols2)].copy()
        meta_aligned = meta_aligned.iloc[meta_aligned['sample_id'].map(col_order3).argsort()].reset_index(drop=True)

    print(f"Running neuroCombat on: {merged_expr.shape[1]} samples × {merged_expr.shape[0]} genes")
    assert list(meta_aligned['sample_id']) == list(merged_expr.columns), "Metadata/expression column mismatch!"

    # neuroCombat covars must be a DataFrame with sample_id matching column order
    covars_df = meta_aligned[['sample_id', 'dataset', 'condition']].copy()

    combat_result = neuroCombat(
        dat=merged_expr,
        covars=covars_df,
        batch_col="dataset",
        categorical_cols=["condition"]
    )

    corrected_expr = pd.DataFrame(
        combat_result["data"],
        index=merged_expr.index,
        columns=merged_expr.columns
    )

    print(f"Post-ComBat range: [{corrected_expr.values.min():.2f}, {corrected_expr.values.max():.2f}]")

    # ─── Step 6: Save + PCA after ComBat ────────────────────
    corrected_path = os.path.join(dest_dir, "Corrected_Batch_Corrected_Expression.csv")
    corrected_expr.to_csv(corrected_path)

    # Save metadata (restore original combined_meta for the records, tag excluded datasets)
    combined_meta['combat_excluded'] = combined_meta['dataset'].isin(single_cond if single_cond else [])
    combined_meta.to_csv(os.path.join(dest_dir, "Corrected_Metadata_Final.csv"), index=False)
    # Also save the combat-aligned metadata (only samples that went through ComBat)
    meta_aligned.to_csv(os.path.join(dest_dir, "Corrected_Metadata_ComBat.csv"), index=False)

    _plot_pca(corrected_expr, meta_aligned, "Corrected_PCA_After_ComBat.png", "After Batch Correction (ComBat)")

    print("\n" + "=" * 60)
    print("DONE!")
    print(f"  Corrected matrix  → {corrected_path}")
    print(f"  Metadata          → Expression_Data/Corrected_Metadata_Final.csv")
    print(f"  PCA plots         → {out_dir}/Corrected_PCA_*.png")
    print("  Ready for Limma DGE with corrected_limma_dge.R")
    print("=" * 60)


# ─────────────────────────────────────────────────────────────
# Helper: PCA plot
# ─────────────────────────────────────────────────────────────

def _plot_pca(expr, meta, filename, title):
    try:
        # Samples × genes
        X = expr.T.values
        X = np.nan_to_num(X, nan=0.0)
        pca = PCA(n_components=2)
        coords = pca.fit_transform(X)

        # Align metadata  
        meta_idx = meta.set_index('sample_id')
        samples = expr.columns.tolist()

        colors_cond = {'MS': '#E05252', 'HC': '#5278E0'}
        dataset_list = meta['dataset'].unique().tolist()
        cmap = plt.cm.get_cmap('tab20', len(dataset_list))
        colors_ds = {ds: cmap(i) for i, ds in enumerate(dataset_list)}

        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        fig.suptitle(title, fontsize=14, fontweight='bold')

        for ax, (color_map, label, legend_title) in zip(axes, [
            (colors_cond, 'condition', 'Condition'),
            (colors_ds,   'dataset',   'Dataset')
        ]):
            for category, color in color_map.items():
                mask = [meta_idx.loc[s, label] == category if s in meta_idx.index else False for s in samples]
                ax.scatter(
                    coords[mask, 0], coords[mask, 1],
                    c=[color], label=category, alpha=0.7, s=40, edgecolors='none'
                )
            ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
            ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
            ax.legend(title=legend_title, fontsize=8, loc='upper right')
            ax.grid(True, alpha=0.3)

        plt.tight_layout()
        path = os.path.join(out_dir, filename)
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  PCA saved: {path}")
    except Exception as e:
        print(f"  PCA plot failed: {e}")


if __name__ == '__main__':
    run()
