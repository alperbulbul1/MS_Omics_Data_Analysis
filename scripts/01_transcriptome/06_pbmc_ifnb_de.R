#!/usr/bin/env Rscript
## 06_pbmc_ifnb_de.R  —  generated from notebook spec
## Run: Rscript 06_pbmc_ifnb_de.R


## ============================================================
## # 06 — IFN-β treatment context (PBMC)
## 
## R/limma rerun of stratum `label_context_case_control_pbmc_ifnb` from
## `Stratified_Analyses/Expression/`. Uses the already batch-corrected
## expression matrix (Python ComBat) as input and runs limma eBayes
## with `trend=TRUE, robust=TRUE`.
## 
## **Outputs**
## - `results/06_pbmc_ifnb_DE.tsv`
## - `figures/06_pbmc_ifnb_volcano.png` + `.pdf`
## 
## **Treatment-context stratum**: 45 MS-on-IFN-β / 8 HC (matched on
## PBMC tissue). Useful for separating disease vs treatment-induced effects.
## ============================================================

suppressPackageStartupMessages({
  library(limma); library(data.table); library(ggplot2); library(ggrepel); library(dplyr)
})
source("helpers.R")


s <- load_stratum("label_context_case_control_pbmc_ifnb")
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


out_fp <- file.path(OUT_DIR, "06_pbmc_ifnb_DE.tsv")
write.table(res, out_fp, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Wrote %s\n", out_fp))

p <- tx_volcano_gg(res,
       title = sprintf("%s — R/limma DE  (eBayes trend+robust)", s$stratum),
       subtitle = sprintf("%d genes · %s vs %s · n=%d",
                          nrow(res), "MS", "HC", ncol(s$mat)))
ggsave(file.path(FIG_DIR, "06_pbmc_ifnb_volcano.png"), p,
       width = 11, height = 8, dpi = 200)
ggsave(file.path(FIG_DIR, "06_pbmc_ifnb_volcano.pdf"), p, width = 11, height = 8)
print(p)

