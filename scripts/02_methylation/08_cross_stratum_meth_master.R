#!/usr/bin/env Rscript
## 08_cross_stratum_meth_master.R  —  generated from notebook spec
## Run: Rscript 08_cross_stratum_meth_master.R


## ============================================================
## # 08 — Cross-stratum methylation master heatmap + cross-omics panel
## 
## Same as transcriptome notebook 08 but for methylation. Pulls the
## gene-level TSVs from notebooks 01-05 and builds a 7 × N heatmap of the
## cross-omics panel across all methylation strata + the combined cohort.
## 
## **Outputs**
## - `results/CrossStratum_Meth_R_Summary.tsv`
## - `figures/CrossStratum_Meth_R_Heatmap.png`
## ============================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(pheatmap); library(reshape2)
})
source("helpers.R")


assays <- list(
  list(name = "T cells (case-ctrl)",   fp = file.path(OUT_DIR, "01_tcells_meth_gene.tsv")),
  list(name = "WB DMF",                fp = file.path(OUT_DIR, "02_wb_dmf_meth_gene.tsv")),
  list(name = "WB Ocrelizumab",        fp = file.path(OUT_DIR, "03_wb_ocrelizumab_meth_gene.tsv")),
  list(name = "T cells remission",     fp = file.path(OUT_DIR, "04_tcells_remission_meth_gene.tsv")),
  list(name = "Combined cohort",       fp = file.path(OUT_DIR, "05_combined_meth_gene.tsv"))
)
panel_genes <- unique(c(CROSS_OMICS, METH_TOP))
rows <- list()
for (a in assays) {
  if (!file.exists(a$fp)) { message("SKIP missing: ", a$fp); next }
  d <- fread(a$fp); d <- as.data.frame(d)
  for (g in panel_genes) {
    r <- d[d$gene == g, , drop = FALSE]
    if (nrow(r) == 0) next
    r <- r[order(r$P.Value), ][1, ]
    rows[[length(rows) + 1]] <- data.frame(
      stratum = a$name, gene = g,
      mean_logFC = r$mean_logFC,
      P.Value = r$P.Value, adj.P.Val = r$adj.P.Val,
      is_cross_omics = g %in% CROSS_OMICS,
      stringsAsFactors = FALSE)
  }
}
summ <- do.call(rbind, rows)
summ$is_sig_FDR05 <- summ$adj.P.Val < 0.05
fwrite(summ, file.path(OUT_DIR, "CrossStratum_Meth_R_Summary.tsv"), sep = "\t")
cat(sprintf("Summary: %d rows (%d strata × %d genes)\n",
            nrow(summ), uniqueN(summ$stratum), uniqueN(summ$gene)))


heat <- dcast(summ, gene ~ stratum, value.var = "mean_logFC", fun.aggregate = mean)
sig  <- dcast(summ, gene ~ stratum, value.var = "adj.P.Val", fun.aggregate = min)
heat_mat <- as.matrix(heat[, -1]); rownames(heat_mat) <- heat$gene
sig_mat  <- as.matrix(sig [, -1]); rownames(sig_mat)  <- sig$gene
star <- ifelse(sig_mat < 0.001, "***",
        ifelse(sig_mat < 0.01,  "**",
        ifelse(sig_mat < 0.05,  "*", "")))
star[is.na(star)] <- ""

row_order <- c(intersect(CROSS_OMICS, rownames(heat_mat)),
               intersect(METH_TOP, rownames(heat_mat)))
row_order <- unique(row_order[row_order %in% rownames(heat_mat)])

pheatmap(heat_mat[row_order, ],
         display_numbers = star[row_order, ], fontsize_number = 11,
         color = colorRampPalette(c("#1F4E79","white","#D62828"))(50),
         na_col = "grey90",
         cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Cross-stratum methylation × cross-omics + meth-top panel",
         filename = file.path(FIG_DIR, "CrossStratum_Meth_R_Heatmap.png"),
         width = 9, height = 9)
cat("Wrote heatmap.\n")

