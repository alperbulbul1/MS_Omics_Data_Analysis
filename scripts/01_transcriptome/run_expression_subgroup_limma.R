#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_expression_subgroup_limma.R <meta_csv> <matrix_csv> <out_dir> [precorrected]")
}

meta_path <- args[[1]]
matrix_path <- args[[2]]
out_dir <- args[[3]]
precorrected <- length(args) >= 4 && args[[4]] == "precorrected"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

meta <- read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
mat_df <- read.csv(matrix_path, stringsAsFactors = FALSE, check.names = FALSE)

gene_ids <- mat_df[[1]]
mat <- as.matrix(mat_df[, -1, drop = FALSE])
rownames(mat) <- gene_ids
storage.mode(mat) <- "numeric"

meta <- meta[meta$condition %in% c("MS", "HC"), , drop = FALSE]
samples <- meta$sample_id[meta$sample_id %in% colnames(mat)]
meta <- meta[match(samples, meta$sample_id), , drop = FALSE]
mat <- mat[, samples, drop = FALSE]

keep <- apply(mat, 1, function(x) var(x, na.rm = TRUE) > 0)
mat <- mat[keep, , drop = FALSE]

if (nrow(meta) < 8 || sum(meta$condition == "MS") < 4 || sum(meta$condition == "HC") < 4) {
  stop("Not enough samples for subgroup analysis")
}

if (nrow(mat) == 0) {
  write.csv(
    data.frame(Gene = character(), logFC = numeric(), AveExpr = numeric(), t = numeric(), P.Value = numeric(), adj.P.Val = numeric(), B = numeric()),
    file.path(out_dir, "DGE_Results_MS_vs_HC.csv"),
    row.names = FALSE
  )
  summary_df <- data.frame(
    metric = c("samples", "datasets", "ms_samples", "hc_samples", "genes_tested", "significant_fdr_0_05", "batch_method"),
    value = c(nrow(meta), length(unique(meta$dataset)), sum(meta$condition == "MS"), sum(meta$condition == "HC"), 0, 0, "none")
  )
  write.csv(summary_df, file.path(out_dir, "Summary.csv"), row.names = FALSE)
  cat("Expression subgroup done: no variable genes remained after filtering\n")
  quit(save = "no", status = 0)
}

batch_corrected <- mat
batch_method <- if (precorrected) "precomputed_batch_corrected" else "none"
if (!precorrected && length(unique(meta$dataset)) > 1) {
  batch_corrected <- removeBatchEffect(mat, batch = meta$dataset, design = model.matrix(~ meta$condition))
  batch_method <- "removeBatchEffect"
}

write.csv(
  data.frame(Gene = rownames(batch_corrected), batch_corrected, check.names = FALSE),
  file.path(out_dir, "Batch_Corrected_Expression.csv"),
  row.names = FALSE
)

condition <- factor(meta$condition, levels = c("HC", "MS"))

if (!precorrected && length(unique(meta$dataset)) > 1) {
  # Batch correction via design covariate (correct limma method)
  batch <- factor(meta$dataset)
  design <- model.matrix(~ condition + batch)
  # Intercept will be HC baseline, conditionMS will be the MS vs HC effect
  colnames(design)[1:2] <- c("Intercept", "MS_vs_HC")
  # Remaining columns are batch effects, which we want to adjust for but not test
} else {
  design <- model.matrix(~ condition)
  colnames(design) <- c("Intercept", "MS_vs_HC")
}

# IMPORTANT: use uncorrected `mat` for lmFit, not `batch_corrected`
fit <- lmFit(mat, design)
fit <- eBayes(fit, robust = TRUE, trend = TRUE)

res <- topTable(fit, coef = "MS_vs_HC", number = Inf, adjust.method = "BH", sort.by = "P")
res$Gene <- rownames(res)
res <- res[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
write.csv(res, file.path(out_dir, "DGE_Results_MS_vs_HC.csv"), row.names = FALSE)

summary_df <- data.frame(
  metric = c("samples", "datasets", "ms_samples", "hc_samples", "genes_tested", "significant_fdr_0_05", "batch_method"),
  value = c(
    nrow(meta),
    length(unique(meta$dataset)),
    sum(meta$condition == "MS"),
    sum(meta$condition == "HC"),
    nrow(res),
    sum(res$adj.P.Val < 0.05, na.rm = TRUE),
    batch_method
  )
)
write.csv(summary_df, file.path(out_dir, "Summary.csv"), row.names = FALSE)

cat(sprintf(
  "Expression subgroup done: %d samples, %d datasets, %d significant genes\n",
  nrow(meta), length(unique(meta$dataset)), sum(res$adj.P.Val < 0.05, na.rm = TRUE)
))
