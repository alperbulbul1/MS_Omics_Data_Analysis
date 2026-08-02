#!/usr/bin/env Rscript
## 05_whole_blood_de.R  —  generated from notebook spec
## Run: Rscript 05_whole_blood_de.R


## ============================================================
## # 05 — Whole blood stratum DE (R/limma)
## 
## R/limma rerun of stratum `cell_tissue_case_control_whole_blood` from
## `Stratified_Analyses/Expression/`. Uses the already batch-corrected
## expression matrix (Python ComBat) as input and runs limma eBayes
## with `trend=TRUE, robust=TRUE`.
## 
## **Outputs**
## - `results/05_whole_blood_DE.tsv`
## - `figures/05_whole_blood_volcano.png` + `.pdf`
## 
## **Cohort**: 26 MS / 47 HC across 3 GSE series. ~107 sig in Python pipeline.
## ============================================================

suppressPackageStartupMessages({
  library(limma); library(data.table); library(ggplot2); library(ggrepel); library(dplyr)
})
source("helpers.R")


s <- load_stratum("cell_tissue_case_control_whole_blood")
cat(sprintf("Stratum: %s  ·  %d genes × %d samples\n",
            s$stratum, nrow(s$mat), ncol(s$mat)))
print(table(s$groups))


res <- run_limma_stratum(s$mat, s$groups, "MS", "HC")
cat(sprintf("Total: %d  ·  P<0.05: %d  ·  FDR<0.05: %d  ·  FDR<0.001: %d\n",
            nrow(res), sum(res$P.Value < 0.05),
            sum(res$adj.P.Val < 0.05),
            sum(res$adj.P.Val < 0.001)))
cat("\n=== Top 15 by FDR ===\n")
print(head(res[, c("gene","logFC","P.Value","adj.P.Val")], 15))
cat("\n=== Cross-omics 7-gene panel ===\n")
print(subset(res, is_cross_omics)[, c("gene","logFC","P.Value","adj.P.Val")])


out_fp <- file.path(OUT_DIR, "05_whole_blood_DE.tsv")
write.table(res, out_fp, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Wrote %s\n", out_fp))

p <- tx_volcano_gg(res,
       title = sprintf("%s — R/limma DE  (eBayes trend+robust)", s$stratum),
       subtitle = sprintf("%d genes · %s vs %s · n=%d",
                          nrow(res), "MS", "HC", ncol(s$mat)))
ggsave(file.path(FIG_DIR, "05_whole_blood_volcano.png"), p,
       width = 11, height = 8, dpi = 200)
ggsave(file.path(FIG_DIR, "05_whole_blood_volcano.pdf"), p, width = 11, height = 8)
print(p)

