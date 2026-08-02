#!/usr/bin/env Rscript
## 05_combined_meth_dmp.R  —  generated from notebook spec
## Run: Rscript 05_combined_meth_dmp.R


## ============================================================
## # 05 — Combined cohort methylation DMP (all 549 samples, R/limma)
## 
## Re-runs limma on the full combined Methylation_Data cohort
## (~549 samples across 17 GSE series; M-value matrix is 317 MB so we
## load it carefully via `data.table::fread` with `select=` to limit RAM).
## 
## Uses **study as batch covariate** in the design matrix.
## 
## **Outputs**
## - `results/05_combined_meth_DMP.tsv` (probe level)
## - `results/05_combined_meth_gene.tsv` (gene level)
## - `figures/05_combined_meth_volcano_{probe,gene}.png`
## ============================================================

suppressPackageStartupMessages({
  library(limma); library(data.table); library(ggplot2); library(dplyr)
})
source("helpers.R")


meta <- fread(file.path(COMBINED_DIR, "Combined_Methylation_Metadata.csv"))
cat(sprintf("Total samples in metadata: %d\n", nrow(meta)))
print(table(meta$condition, meta$cell_type))


# Load combined matrix (selecting only annotated samples)
keep <- meta$sample_id[meta$condition %in% c("MS","HC")]
mat_fp <- file.path(COMBINED_DIR, "Combined_Methylation_Batch_Corrected.csv")
cat("Reading combined matrix (this may take ~30 s, 317 MB)...\n")
hdr <- fread(mat_fp, nrows = 0); cn <- colnames(hdr)
probe_col <- cn[1]
use_cols <- c(probe_col, intersect(keep, cn))
expr <- fread(mat_fp, select = use_cols, showProgress = TRUE)
probes <- expr[[probe_col]]; sample_cols <- setdiff(colnames(expr), probe_col)
mat <- as.matrix(expr[, ..sample_cols]); storage.mode(mat) <- "numeric"
rownames(mat) <- probes
meta_sub <- meta[match(sample_cols, meta$sample_id)]
cat(sprintf("Combined matrix: %d probes × %d samples\n", nrow(mat), ncol(mat)))
print(table(meta_sub$condition, meta_sub$cell_type))
rm(expr, hdr); invisible(gc())


# Run limma with study covariate
dmp <- run_limma_meth(mat, meta_sub$condition, batch = meta_sub$dataset,
                       group_a = "MS", group_b = "HC")
cat(sprintf("Probes: %d  ·  FDR<0.05: %d  ·  FDR<0.001: %d\n",
            nrow(dmp), sum(dmp$adj.P.Val < 0.05),
            sum(dmp$adj.P.Val < 0.001)))
print(head(dmp[, c("Probe","logFC","P.Value","adj.P.Val")], 15))


# Probe-to-gene aggregation
if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE))
  BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19",
                       update = FALSE, ask = FALSE)
p2g <- annotate_probes_to_genes(dmp$Probe)
gene_level <- probe_to_gene_stouffer(dmp, p2g)
gene_level$is_cross_omics <- gene_level$gene %in% CROSS_OMICS
cat(sprintf("Gene-level: %d  ·  FDR<0.05: %d  ·  FDR<0.001: %d\n",
            nrow(gene_level),
            sum(gene_level$adj.P.Val < 0.05),
            sum(gene_level$adj.P.Val < 0.001)))

cat("\n=== Cross-omics 7-gene panel (combined cohort) ===\n")
print(subset(gene_level, is_cross_omics)[, c("gene","n_probes","mean_logFC","P.Value","adj.P.Val")])

write.table(dmp, file.path(OUT_DIR, "05_combined_meth_DMP.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(gene_level, file.path(OUT_DIR, "05_combined_meth_gene.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


dmp_lbl <- merge(dmp, p2g, by = "Probe", all.x = TRUE)
dmp_lbl$gene[is.na(dmp_lbl$gene)] <- dmp_lbl$Probe[is.na(dmp_lbl$gene)]
p_probe <- meth_volcano_gg(dmp_lbl,
              title = "Combined cohort methylation DMP (R/limma)",
              subtitle = sprintf("%d probes · %d samples · study covariate",
                                 nrow(dmp_lbl), ncol(mat)))
ggsave(file.path(FIG_DIR, "05_combined_meth_volcano_probe.png"), p_probe,
       width = 11, height = 8, dpi = 200)
gene_level$logFC <- gene_level$mean_logFC
p_gene <- meth_volcano_gg(gene_level,
              title = "Combined cohort methylation — gene-level (Stouffer)",
              subtitle = sprintf("%d genes · %d samples", nrow(gene_level), ncol(mat)))
ggsave(file.path(FIG_DIR, "05_combined_meth_volcano_gene.png"), p_gene,
       width = 11, height = 8, dpi = 200)
print(p_gene)

