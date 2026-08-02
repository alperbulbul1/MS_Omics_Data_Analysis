#!/usr/bin/env Rscript
## 07_brainwm_rna_meth_rerun.R  —  generated from notebook spec
## Run: Rscript 07_brainwm_rna_meth_rerun.R


## ============================================================
## # 07 — Brain WM RNA + methylation rerun (R/limma)
## 
## Re-run the brain white-matter stratum used in the bulk MS_GEO analysis,
## this time in R. RNA matrix comes from
## `Stratified_Analyses/Expression/cell_tissue_case_control_brain_wm/`
## (batch-corrected via removeBatchEffect in the original Python pipeline);
## methylation gene-level results are loaded from the existing rerun TSV
## (`processed/rerun/BrainWM_meth_genelevel_rerun.tsv`) for the R-side
## inverse-concordance scatter.
## 
## **Outputs**
## - `processed/META/BrainWM_R_RNA_DE.tsv`
## - `figures/BrainWM_R_RNA_volcano.png`
## - `figures/BrainWM_R_meth_volcano.png`
## - `figures/BrainWM_R_inverse_concordance.png`
## ============================================================

suppressPackageStartupMessages({
  library(limma); library(data.table); library(ggplot2); library(ggrepel)
  library(dplyr); library(gridExtra)
})
source("helpers.R")


rna_dir <- file.path(PROJ_ROOT, "Stratified_Analyses", "Expression",
                      "cell_tissue_case_control_brain_wm")
expr <- fread(file.path(rna_dir, "Batch_Corrected_Expression.csv"))
meta <- fread(file.path(rna_dir, "metadata.csv"))
cat(sprintf("Expression: %d × %d  |  metadata: %d rows\n",
            nrow(expr), ncol(expr), nrow(meta)))

gene_col <- "Gene"
genes <- expr[[gene_col]]
sample_cols <- setdiff(colnames(expr), gene_col)
mat <- as.matrix(expr[, ..sample_cols])
storage.mode(mat) <- "numeric"
mat <- gene_dedup(mat, genes)
cat(sprintf("Gene-dedup'd: %d × %d\n", nrow(mat), ncol(mat)))

# Match metadata to columns
keep <- intersect(colnames(mat), meta$sample_id)
mat <- mat[, keep, drop = FALSE]
groups <- meta$condition[match(keep, meta$sample_id)]
groups <- ifelse(groups == "MS", "MS", "HC")
cat(sprintf("After matching: %d samples (%s)\n",
            ncol(mat), paste(table(groups), collapse = "/")))


mat <- mat[apply(!is.na(mat), 1, all), ]   # complete cases only
cat(sprintf("After complete-case filter: %d genes × %d samples\n",
            nrow(mat), ncol(mat)))
res_rna <- moderated_t_safe(mat, groups, "MS", "HC")
res_rna <- res_rna[order(res_rna$P.Value), ]
res_rna$is_cross_omics <- res_rna$gene %in% CROSS_OMICS
res_rna$is_recurring   <- res_rna$gene %in% RECURRING

cat(sprintf("Brain WM RNA: %d  |  P<0.05: %d  |  FDR<0.05: %d\n",
            nrow(res_rna), sum(res_rna$P.Value < 0.05),
            sum(res_rna$adj.P.Val < 0.05)))
print(subset(res_rna, is_cross_omics)[, c("gene","logFC","P.Value","adj.P.Val")])

write.table(res_rna, file.path(OUT_DIR, "BrainWM_R_RNA_DE.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


p_rna <- dep_volcano_gg(res_rna,
            title = "Brain WM transcriptomics — R/limma rerun",
            subtitle = sprintf("%d genes · %d MS vs %d HC · eBayes(trend, robust)",
                               nrow(res_rna), sum(groups=="MS"), sum(groups=="HC")))
ggsave(file.path(FIG_DIR, "BrainWM_R_RNA_volcano.png"), p_rna,
       width = 11, height = 8, dpi = 200)
print(p_rna)


# ---- Methylation rerun loading + volcano ----
meth_fp <- file.path(PROT_ROOT, "processed", "rerun", "BrainWM_meth_genelevel_rerun.tsv")
res_meth <- fread(meth_fp)
setnames(res_meth, c("mean_logFC","P.Value","adj.P.Val"),
                   c("logFC","P.Value","adj.P.Val"), skip_absent = TRUE)
res_meth <- as.data.frame(res_meth)
res_meth$gene <- res_meth$Gene
res_meth$is_cross_omics <- res_meth$gene %in% CROSS_OMICS
res_meth$is_recurring   <- res_meth$gene %in% RECURRING
cat(sprintf("Brain WM meth (gene-level): %d  |  FDR<0.05: %d\n",
            nrow(res_meth), sum(res_meth$adj.P.Val < 0.05, na.rm = TRUE)))
print(subset(res_meth, is_cross_omics)[, c("gene","logFC","P.Value","adj.P.Val")])

p_meth <- dep_volcano_gg(res_meth,
              title = "Brain WM methylation (gene-level, Stouffer)",
              subtitle = "From processed/rerun/BrainWM_meth_genelevel_rerun.tsv")
ggsave(file.path(FIG_DIR, "BrainWM_R_meth_volcano.png"), p_meth,
       width = 11, height = 8, dpi = 200)
print(p_meth)


# ---- Inverse concordance: RNA logFC vs Meth logFC ----
common <- intersect(res_rna$gene, res_meth$gene)
df <- data.frame(
  gene = common,
  rna = res_rna$logFC[match(common, res_rna$gene)],
  meth = res_meth$logFC[match(common, res_meth$gene)],
  rna_sig = res_rna$adj.P.Val[match(common, res_rna$gene)] < 0.05,
  meth_sig = res_meth$adj.P.Val[match(common, res_meth$gene)] < 0.05,
  is_co = common %in% CROSS_OMICS,
  stringsAsFactors = FALSE
)
df$both_sig <- df$rna_sig & df$meth_sig
df$inverse  <- df$both_sig & (sign(df$rna) != sign(df$meth))
cat(sprintf("Common: %d  |  both sig: %d  |  inverse-concordant: %d\n",
            nrow(df), sum(df$both_sig, na.rm = TRUE),
            sum(df$inverse, na.rm = TRUE)))

p_inv <- ggplot(df, aes(rna, meth)) +
  geom_point(data = subset(df, !both_sig), colour = "grey80", size = 0.6) +
  geom_point(data = subset(df, both_sig & !inverse),
             colour = "#FFC107", size = 1.4, alpha = 0.7) +
  geom_point(data = subset(df, inverse),
             colour = "#1F4E79", size = 2.0) +
  geom_point(data = subset(df, is_co),
             shape = 8, colour = "#D62828", size = 5, stroke = 1.4) +
  geom_text_repel(data = subset(df, is_co),
                  aes(label = gene), size = 3.5, fontface = "bold",
                  colour = "#D62828") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed",
              colour = "#D62828") +
  labs(title = "Brain WM RNA vs methylation (inverse-concordance)",
       subtitle = sprintf("%d common genes · cross-omics highlighted (red star)",
                          nrow(df)),
       x = "RNA logFC (MS / HC)", y = "Methylation mean logFC (MS / HC)") +
  theme_classic(base_size = 11)
ggsave(file.path(FIG_DIR, "BrainWM_R_inverse_concordance.png"),
       p_inv, width = 9, height = 8, dpi = 200)
print(p_inv)

