#!/usr/bin/env Rscript
## 07_total_combined_de.R  —  generated from notebook spec
## Run: Rscript 07_total_combined_de.R


## ============================================================
## # 07 — Pan-tissue combined cohort DE (R/limma + tissue covariate)
## 
## Combines the 5 case-control strata (PBMC, T cells, B cells, Brain WM,
## Whole blood) into a single multi-tissue cohort, then runs limma with
## **tissue + condition** in the design matrix. The condition coefficient
## is the pan-tissue MS-vs-HC effect adjusted for tissue.
## 
## **Outputs**
## - `results/07_pan_tissue_DE.tsv`
## - `figures/07_pan_tissue_volcano.png`
## ============================================================

suppressPackageStartupMessages({
  library(limma); library(data.table); library(ggplot2); library(dplyr)
})
source("helpers.R")


strata <- c("cell_tissue_case_control_pbmc",
             "cell_tissue_case_control_t_cells",
             "cell_tissue_case_control_b_cells",
             "cell_tissue_case_control_brain_wm",
             "cell_tissue_case_control_whole_blood")
tissues <- c("PBMC","T cells","B cells","Brain WM","Whole blood")

mats <- list(); meta_all <- list()
for (i in seq_along(strata)) {
  s <- load_stratum(strata[i])
  meta_i <- data.frame(
    sample_id = colnames(s$mat),
    condition = s$groups,
    tissue    = tissues[i],
    stratum   = s$stratum,
    stringsAsFactors = FALSE)
  mats[[i]]     <- s$mat
  meta_all[[i]] <- meta_i
}
meta <- do.call(rbind, meta_all)
common_genes <- Reduce(intersect, lapply(mats, rownames))
cat(sprintf("Common genes across 5 strata: %d\n", length(common_genes)))
mat_all <- do.call(cbind, lapply(mats, function(m) m[common_genes, , drop = FALSE]))
cat(sprintf("Combined matrix: %d × %d  (samples)\n", nrow(mat_all), ncol(mat_all)))
print(table(meta$tissue, meta$condition))


# limma with tissue covariate
tissue_f    <- factor(meta$tissue)
condition_f <- factor(meta$condition, levels = c("HC","MS"))
design <- model.matrix(~ tissue_f + condition_f)
coef_name <- tail(colnames(design), 1)
cat("Design coef of interest:", coef_name, "\n")
fit <- lmFit(mat_all, design)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)
res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
res$gene <- rownames(res)
res$is_cross_omics <- res$gene %in% CROSS_OMICS
res$is_recurring   <- res$gene %in% RECURRING

cat(sprintf("\nPan-tissue DE: %d genes  ·  P<0.05: %d  ·  FDR<0.05: %d  ·  FDR<0.001: %d\n",
            nrow(res), sum(res$P.Value < 0.05),
            sum(res$adj.P.Val < 0.05),
            sum(res$adj.P.Val < 0.001)))
cat("\n=== Top 15 by FDR ===\n")
print(head(res[, c("gene","logFC","P.Value","adj.P.Val")], 15))
cat("\n=== Cross-omics in pan-tissue ===\n")
print(subset(res, is_cross_omics)[, c("gene","logFC","P.Value","adj.P.Val")])

write.table(res, file.path(OUT_DIR, "07_pan_tissue_DE.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


p <- tx_volcano_gg(res,
       title = "Pan-tissue combined MS-vs-HC DE  (R/limma, tissue covariate)",
       subtitle = sprintf("%d common genes · %d samples  (PBMC+Tcells+Bcells+BrainWM+WB)",
                          nrow(res), ncol(mat_all)))
ggsave(file.path(FIG_DIR, "07_pan_tissue_volcano.png"), p,
       width = 11, height = 8, dpi = 200)
ggsave(file.path(FIG_DIR, "07_pan_tissue_volcano.pdf"), p, width = 11, height = 8)
print(p)

