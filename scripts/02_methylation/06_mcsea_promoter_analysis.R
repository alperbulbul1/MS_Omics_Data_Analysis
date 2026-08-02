#!/usr/bin/env Rscript
## 06_mcsea_promoter_analysis.R  —  generated from notebook spec
## Run: Rscript 06_mcsea_promoter_analysis.R


## ============================================================
## # 06 — mCSEA promoter + gene-body methylation enrichment
## 
## Uses the `mCSEA` Bioconductor package on the combined cohort DMP results
## (notebook 05) to identify genes with concentrated promoter or gene-body
## methylation changes. mCSEA does a GSEA-like enrichment on probe ranks
## within each gene region (promoter / body / 3'UTR / etc.).
## 
## **Outputs**
## - `results/06_mCSEA_promoter.tsv`
## - `results/06_mCSEA_gene_body.tsv`
## ============================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr)
})
source("helpers.R")
if (!requireNamespace("mCSEA", quietly = TRUE)) {
  cat("mCSEA missing — installing...\n")
  BiocManager::install("mCSEA", update = FALSE, ask = FALSE)
}
library(mCSEA)


# mCSEA requires the M-value matrix AND a phenotype frame; reload
# the combined cohort from disk (subset to PBMC + T cells to fit RAM
# while still being representative of the dominant strata).
meta <- fread(file.path(COMBINED_DIR, "Combined_Methylation_Metadata.csv"))
keep <- meta$sample_id[meta$condition %in% c("MS","HC")]
mat_fp <- file.path(COMBINED_DIR, "Combined_Methylation_Batch_Corrected.csv")
cat("Reading combined matrix...\n")
hdr <- fread(mat_fp, nrows = 0); cn <- colnames(hdr)
use_cols <- c(cn[1], intersect(keep, cn))
expr <- fread(mat_fp, select = use_cols, showProgress = TRUE)
probes <- expr[[cn[1]]]; sample_cols <- setdiff(colnames(expr), cn[1])
mat <- as.matrix(expr[, ..sample_cols]); storage.mode(mat) <- "numeric"
rownames(mat) <- probes
pheno <- data.frame(group = meta$condition[match(sample_cols, meta$sample_id)],
                    row.names = sample_cols)
cat(sprintf("Loaded: %d probes × %d samples\n", nrow(mat), ncol(mat)))
rm(expr, hdr); invisible(gc())


# rankProbes uses limma DE internally + returns ranking
rank_vec <- rankProbes(mat, pheno, refGroup = "HC", caseGroup = "MS",
                        typeInput = "M", typeAnalysis = "M")
cat(sprintf("Ranked %d probes\n", length(rank_vec)))

# Promoter analysis
prom_res <- mCSEATest(rank = rank_vec, methData = mat, pheno = pheno,
                       regionsTypes = "promoters", platform = "450k",
                       minCpGs = 5, nproc = 1)
prom_df <- as.data.frame(prom_res$promoters)
prom_df$gene <- rownames(prom_df)
prom_df <- prom_df[order(prom_df$pval), ]
cat(sprintf("Promoter regions tested: %d  ·  FDR<0.05: %d  ·  FDR<0.001: %d\n",
            nrow(prom_df),
            sum(prom_df$padj < 0.05, na.rm = TRUE),
            sum(prom_df$padj < 0.001, na.rm = TRUE)))
cat("\n=== Top 15 promoter hits ===\n")
print(head(prom_df[, c("gene","size","NES","pval","padj","leadingEdge")], 15))
cat("\n=== Cross-omics 7-gene promoters ===\n")
co_prom <- prom_df[prom_df$gene %in% CROSS_OMICS, ]
print(co_prom[, c("gene","size","NES","pval","padj")])
write.table(prom_df, file.path(OUT_DIR, "06_mCSEA_promoter.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# Gene-body (regionsTypes = "genes")
body_res <- try(mCSEATest(rank = rank_vec, methData = mat, pheno = pheno,
                           regionsTypes = "genes", platform = "450k",
                           minCpGs = 5, nproc = 1), silent = TRUE)
if (!inherits(body_res, "try-error")) {
  body_df <- as.data.frame(body_res$genes)
  body_df$gene <- rownames(body_df)
  body_df <- body_df[order(body_df$pval), ]
  cat(sprintf("Gene-body regions tested: %d  ·  FDR<0.05: %d\n",
              nrow(body_df), sum(body_df$padj < 0.05, na.rm = TRUE)))
  print(head(body_df[, c("gene","size","NES","pval","padj")], 10))
  write.table(body_df, file.path(OUT_DIR, "06_mCSEA_gene_body.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  cat("Gene-body analysis failed:\n",
      conditionMessage(attr(body_res, "condition")), "\n")
}

