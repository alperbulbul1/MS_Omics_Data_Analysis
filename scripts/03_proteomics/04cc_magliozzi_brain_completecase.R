#!/usr/bin/env Rscript
## 04_magliozzi_brain_dep.R  —  generated from notebook spec
## Run: Rscript 04_magliozzi_brain_dep.R


## ============================================================
## # 04 — Magliozzi 2026 brain proteomics: DEP/limma, 4 contrasts
## 
## DIA-MS on **post-mortem brain tissue** (Magliozzi et al. 2026, Nat Commun,
## 10.1038/s41467-025-68118-0). 3,575 detection-filtered proteins × 37 samples:
## 8 ODC CTX · 8 MS CTX · 8 ODC WM · 8 MS NAWM · 5 MS WML.
## 
## Four contrasts (each via separate limma fit):
## - MS CTX vs ODC CTX
## - MS NAWM vs ODC WM
## - MS WML vs ODC WM
## - **MS WML vs MS NAWM ← LXN's home contrast**
## 
## **Outputs**
## - `processed/META/Magliozzi_CC_<contrast>.tsv` × 4
## - `figures/Magliozzi_CC_<contrast>_volcano.png` × 4
## - `figures/Magliozzi_CC_4panel.png`
## ============================================================

suppressPackageStartupMessages({
  library(limma); library(ggplot2); library(dplyr)
  library(data.table); library(gridExtra)
})
source("helpers.R")


xl_fp <- file.path(PROT_ROOT, "processed", "Magliozzi_S1_sheet3.tsv")
sh3 <- as.data.frame(data.table::fread(xl_fp, sep="	", header=TRUE))
cat(sprintf("Loaded: %d proteins × %d cols\n", nrow(sh3), ncol(sh3)))

# openxlsx replaces spaces with dots in column names — match both forms
sample_cols <- grep("^(ODC|MS)[. ]", colnames(sh3), value = TRUE)
cat("Sample columns:", length(sample_cols), "\n")
cond_raw <- gsub("[.]\\d+$", "", sample_cols)
cond_raw <- gsub(" \\d+$", "", cond_raw)
print(table(cond_raw))

# Build numeric matrix
expr_raw <- as.matrix(sh3[, sample_cols])
expr_raw <- apply(expr_raw, 2, function(x) suppressWarnings(as.numeric(x)))
expr_raw <- log2(expr_raw + 1)  # original is linear intensity
rownames(expr_raw) <- sh3$Gene
cat(sprintf("Expression matrix: %d × %d  (log2 transformed)\n",
            nrow(expr_raw), ncol(expr_raw)))

cond <- gsub("[. ]", "_", cond_raw)  # 'ODC.CTX' -> 'ODC_CTX'
print(table(cond))


# Drop proteins without a gene symbol
keep <- !is.na(rownames(expr_raw)) & rownames(expr_raw) != ""
cat(sprintf("Drop NA/empty gene rows: %d -> %d\n", nrow(expr_raw), sum(keep)))
expr_raw <- expr_raw[keep, , drop = FALSE]

# ---- DEP-equivalent filter+normalise+impute once, contrast per contrast ----
expr <- filter_missval_R(expr_raw, cond, thr = 0.5)
expr <- vsn_with_fallback(expr)
# COMPLETE-CASE: no imputation (MAR assumption)
cat(sprintf("After filter+norm+impute: %d × %d\n", nrow(expr), ncol(expr)))

# Collapse duplicate gene symbols by mean (only meaningful when >1 row/gene)
if (anyDuplicated(rownames(expr)) > 0) {
  ids <- rownames(expr)
  expr <- limma::avereps(expr, ID = ids)
  cat(sprintf("After avereps: %d genes\n", nrow(expr)))
}


# ---- 4 contrasts via limma ----
contrasts_list <- list(
  MS_CTX_vs_ODC_CTX = c(grp_a = "MS_CTX",  grp_b = "ODC_CTX"),
  MS_NAWM_vs_ODC_WM = c(grp_a = "MS_NAWM", grp_b = "ODC_WM"),
  MS_WML_vs_ODC_WM  = c(grp_a = "MS_WML",  grp_b = "ODC_WM"),
  MS_WML_vs_MS_NAWM = c(grp_a = "MS_WML",  grp_b = "MS_NAWM")
)

all_results <- list()
all_plots <- list()
for (cn in names(contrasts_list)) {
  ga <- contrasts_list[[cn]]["grp_a"]
  gb <- contrasts_list[[cn]]["grp_b"]
  keep <- cond %in% c(ga, gb)
  res <- moderated_t_safe(expr[, keep, drop = FALSE], cond[keep], ga, gb)
  res <- res[order(res$P.Value), ]
  res$contrast <- cn
  res$is_cross_omics <- res$gene %in% CROSS_OMICS
  res$is_ECM        <- res$gene %in% ECM_FAMILY
  res$is_recurring  <- res$gene %in% RECURRING
  cat(sprintf("\n%-22s  n=%d/%d  |  p<0.05: %d  |  FDR<0.05: %d\n",
              cn, sum(cond == ga), sum(cond == gb),
              sum(res$P.Value < 0.05),
              sum(res$adj.P.Val < 0.05)))
  cat("  cross-omics:\n"); print(subset(res, is_cross_omics)[, c("gene","logFC","P.Value","adj.P.Val")])

  tsv <- file.path(OUT_DIR, sprintf("Magliozzi_CC_%s.tsv", cn))
  write.table(res, tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  p <- dep_volcano_gg(res, title = sprintf("Magliozzi 2026 — %s", cn),
                      subtitle = sprintf("n=%d vs %d  ·  limma eBayes(trend, robust)",
                                         sum(cond == ga), sum(cond == gb)))
  ggsave(file.path(FIG_DIR, sprintf("Magliozzi_CC_%s_volcano.png", cn)),
         p, width = 10, height = 7, dpi = 200)
  all_results[[cn]] <- res
  all_plots[[cn]]   <- p
}


# ---- 2x2 grid ----
ga <- gridExtra::arrangeGrob(grobs = all_plots, ncol = 2)
ggsave(file.path(FIG_DIR, "Magliozzi_CC_4panel.png"), ga,
       width = 18, height = 14, dpi = 180)
grid::grid.draw(ga)

