#!/usr/bin/env Rscript
## 08_cross_stratum_master.R  —  generated from notebook spec
## Run: Rscript 08_cross_stratum_master.R


## ============================================================
## # 08 — Cross-stratum master heatmap + cross-omics 7-gene panel
## 
## Pulls per-stratum DE TSVs from notebooks 01–07 and builds:
## 
## 1. **Master 16-gene × 7-stratum significance heatmap** — for the cross-omics
##    7-gene panel + 9 prior top-recurring genes (STAT3 / TYK2 / CXCL13 etc.).
## 2. **Effect-size matrix** with colour = log2FC, text = FDR stars.
## 3. **Triple-validated set** — genes with FDR<0.05 in ≥3 strata simultaneously.
## 
## **Outputs**
## - `results/CrossStratum_R_Summary.tsv`
## - `figures/CrossStratum_R_Heatmap.png`
## - `results/Triple_Stratum_Validated_R.tsv`
## ============================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(pheatmap)
})
source("helpers.R")


assays <- list(
  list(name = "PBMC",                 fp = file.path(OUT_DIR, "01_pbmc_DE.tsv")),
  list(name = "T cells",              fp = file.path(OUT_DIR, "02_tcells_DE.tsv")),
  list(name = "B cells",              fp = file.path(OUT_DIR, "03_bcells_DE.tsv")),
  list(name = "Brain WM",             fp = file.path(OUT_DIR, "04_brainwm_DE.tsv")),
  list(name = "Whole blood",          fp = file.path(OUT_DIR, "05_whole_blood_DE.tsv")),
  list(name = "IFN-beta PBMC",        fp = file.path(OUT_DIR, "06_pbmc_ifnb_DE.tsv")),
  list(name = "Combined (pan-tissue)",fp = file.path(OUT_DIR, "07_total_DE.tsv"))
)
panel_genes <- unique(c(CROSS_OMICS, RECURRING))
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
      logFC = r$logFC, P.Value = r$P.Value, adj.P.Val = r$adj.P.Val,
      is_cross_omics = g %in% CROSS_OMICS,
      stringsAsFactors = FALSE)
  }
}
summ <- do.call(rbind, rows)
summ$is_sig_FDR05 <- summ$adj.P.Val < 0.05
fwrite(summ, file.path(OUT_DIR, "CrossStratum_R_Summary.tsv"), sep = "\t")
cat(sprintf("Built summary: %d rows (%d strata × %d genes)\n",
            nrow(summ), uniqueN(summ$stratum), uniqueN(summ$gene)))


# Pivot to heatmap
library(reshape2)
heat <- dcast(summ, gene ~ stratum, value.var = "logFC", fun.aggregate = mean)
sig  <- dcast(summ, gene ~ stratum, value.var = "adj.P.Val", fun.aggregate = min)
heat_mat <- as.matrix(heat[, -1]); rownames(heat_mat) <- heat$gene
sig_mat  <- as.matrix(sig [, -1]); rownames(sig_mat)  <- sig$gene

star <- ifelse(sig_mat < 0.001, "***",
        ifelse(sig_mat < 0.01,  "**",
        ifelse(sig_mat < 0.05,  "*", "")))
star[is.na(star)] <- ""

# Reorder rows: cross-omics first then recurring
row_order <- c(intersect(CROSS_OMICS, rownames(heat_mat)),
               intersect(RECURRING, rownames(heat_mat)))
row_order <- unique(row_order[row_order %in% rownames(heat_mat)])

pheatmap(heat_mat[row_order, ],
         display_numbers = star[row_order, ], fontsize_number = 11,
         color = colorRampPalette(c("#1F4E79","white","#D62828"))(50),
         na_col = "grey90",
         cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Cross-stratum cross-omics + recurring × bulk transcriptome strata",
         filename = file.path(FIG_DIR, "CrossStratum_R_Heatmap.png"),
         width = 11, height = 7.5)
cat("Wrote heatmap.\n")


# Multi-stratum validated: gene sig (FDR<0.05) in >=2 strata
gene_sig_count <- aggregate(is_sig_FDR05 ~ gene, data = summ, FUN = sum)
gene_sig_count <- gene_sig_count[order(-gene_sig_count$is_sig_FDR05), ]
multi <- gene_sig_count[gene_sig_count$is_sig_FDR05 >= 2, ]
cat(sprintf("Genes sig (FDR<0.05) in >=2 strata: %d\n", nrow(multi)))
print(multi)
cat("\nFull per-gene sig counts:\n"); print(gene_sig_count)
fwrite(multi, file.path(OUT_DIR, "MultiStratum_Validated_R.tsv"), sep = "\t")

