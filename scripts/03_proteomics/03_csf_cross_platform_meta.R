#!/usr/bin/env Rscript
## SUPERSEDED / STILL IMPUTED. This cross-platform ComBat meta calls the MinProb helper on
## both platforms before concatenation, which is the procedure the revision withdrew: it is
## what manufactured the spurious MS-up ITGB2 CSF call. Its output
## CSF_combined_R_ComBat_DE.tsv is no longer read by any validation script and no reported
## number depends on it. Retained for provenance only.
## 03_csf_cross_platform_meta.R  —  generated from notebook spec
## Run: Rscript 03_csf_cross_platform_meta.R


## ============================================================
## # 03 — CSF cross-platform meta (Astral + timsTOF, ComBat + limma)
## 
## Properly integrate the two Cell-2026 CSF platforms with **sva::ComBat**
## (batch = platform, mod = ~group) rather than the Stouffer-Z meta used
## in the Python pipeline. ComBat preserves biological signal while
## removing the systematic platform shift.
## 
## Workflow
## 1. Read both DEP-pipeline results (notebooks 01 + 02) for the per-platform
##    reference numbers.
## 2. Re-load the gene-keyed matrices, inner-join on common genes.
## 3. ComBat(batch = platform, mod = model.matrix(~ group)).
## 4. PCA before vs after (sanity check).
## 5. limma::lmFit → makeContrasts(MS - Control) → eBayes.
## 6. EnhancedVolcano-style output.
## 
## **Outputs**
## - `processed/META/CSF_combined_R_ComBat_DE.tsv`
## - `figures/CSF_combined_R_PCA.png`
## - `figures/CSF_combined_R_volcano.png`
## ============================================================

suppressPackageStartupMessages({
  library(limma); library(sva); library(data.table); library(ggplot2)
  library(dplyr); library(gridExtra)
})
source("helpers.R")


# ---- Load both gene-keyed matrices ----
astral_fp  <- file.path(PROT_ROOT, "processed", "astral_discovery_gene_keyed.tsv")
timstof_fp <- file.path(PROT_ROOT, "processed", "timsTOF_gene_mapped.tsv")
cat("Reading Astral...\n"); A <- fread(astral_fp, sep="\t", showProgress=FALSE)
cat("Reading timsTOF...\n"); T <- fread(timstof_fp, sep="\t", showProgress=FALSE,
                                          nThread=4)

ann_fp <- file.path(PROT_ROOT, "osfstorage-archive",
  "processed proteomic data", "0_sample_annotations",
  "annotations_v42_49_2_10_4_10_interimSky17_PL01-PL56_PepResCustv01_resubmission.tsv")
ann <- fread(ann_fp, sep = "\t", header = TRUE)

# Match each platform's sample columns
meta <- c("protein","Genes","Gene","data_completeness","completeness")

ann_a <- ann[Run_Astral_Measurement != "" & !is.na(Run_Astral_Measurement),
             .(run = Run_Astral_Measurement, Diagnosis_group, MSgroup)]
ann_a[, group := fifelse(MSgroup == "MS", "MS",
                  fifelse(Diagnosis_group %in% c("Other","Neurological Control"),
                          "Control", NA_character_))]
ann_a <- ann_a[!is.na(group)]

run_col_t <- grep("timsTOF|tTOF|Run_timstof|Run_tims", colnames(ann),
                  ignore.case=TRUE, value=TRUE)[1]
ann_t <- ann[get(run_col_t) != "" & !is.na(get(run_col_t)),
             .(run = get(run_col_t), Diagnosis_group, MSgroup)]
ann_t[, group := fifelse(MSgroup == "MS", "MS",
                  fifelse(Diagnosis_group %in% c("Other","Neurological Control"),
                          "Control", NA_character_))]
ann_t <- ann_t[!is.na(group)]


# ---- Build per-platform matrices keyed on Genes ----
gene_dedupe_dt <- function(dt, sample_cols, gene_col = "Genes") {
  dt <- as.data.frame(dt)
  m <- as.matrix(dt[, sample_cols, drop = FALSE])
  storage.mode(m) <- "numeric"
  rownames(m) <- dt[[gene_col]]
  v <- apply(m, 1, var, na.rm = TRUE)
  ord <- order(-v, na.last = TRUE)
  m   <- m[ord, , drop = FALSE]
  g   <- rownames(m)
  keep <- !duplicated(g) & !is.na(g) & g != ""
  m[keep, , drop = FALSE]
}

# Sample columns are .raw (Astral) or .d (timsTOF) run names
A_runs <- grep("\\.raw$|^[0-9]{8}_", colnames(A), value = TRUE)
T_runs <- grep("\\.d$|\\.raw$|^[0-9]{8}_", colnames(T), value = TRUE)
A_cols <- intersect(A_runs, ann_a$run)
T_cols <- intersect(T_runs, ann_t$run)
A_g <- ann_a$group[match(A_cols, ann_a$run)]
T_g <- ann_t$group[match(T_cols, ann_t$run)]

Am <- gene_dedupe_dt(A[, c("Genes", A_cols), with = FALSE], A_cols)
Tm <- gene_dedupe_dt(T[, c("Genes", T_cols), with = FALSE], T_cols)
cat(sprintf("Astral matrix:  %d genes × %d samples (MS=%d Ctrl=%d)\n",
            nrow(Am), ncol(Am), sum(A_g=="MS"), sum(A_g=="Control")))
cat(sprintf("timsTOF matrix: %d genes × %d samples (MS=%d Ctrl=%d)\n",
            nrow(Tm), ncol(Tm), sum(T_g=="MS"), sum(T_g=="Control")))

common <- intersect(rownames(Am), rownames(Tm))
cat(sprintf("Common genes: %d\n", length(common)))


# ---- Filter per-platform to >=50% valid per condition, then median-norm ----
Am1 <- filter_missval_R(Am[common, ], A_g, thr = 0.5)
Tm1 <- filter_missval_R(Tm[common, ], T_g, thr = 0.5)

common2 <- intersect(rownames(Am1), rownames(Tm1))
cat(sprintf("Common after filtering: %d genes\n", length(common2)))

Am2 <- vsn_with_fallback(Am1[common2, ])
Tm2 <- vsn_with_fallback(Tm1[common2, ])
Am3 <- impute_minprob_R(Am2)
Tm3 <- impute_minprob_R(Tm2)

# ---- Concatenate + ComBat ----
combined <- cbind(Am3, Tm3)
platform <- c(rep("Astral", ncol(Am3)), rep("timsTOF", ncol(Tm3)))
groups   <- c(A_g, T_g)
cat(sprintf("Combined: %d genes × %d samples  (Astral=%d, timsTOF=%d)\n",
            nrow(combined), ncol(combined), ncol(Am3), ncol(Tm3)))

mod <- model.matrix(~ groups)
combat_mat <- ComBat(dat = combined, batch = platform, mod = mod,
                     par.prior = TRUE)


# ---- PCA before vs after ----
make_pca <- function(mat, ttl) {
  pc <- prcomp(t(mat), scale. = TRUE)
  data.frame(PC1 = pc$x[,1], PC2 = pc$x[,2],
             platform = platform, group = groups,
             stringsAsFactors = FALSE) |>
    ggplot(aes(PC1, PC2, colour = platform, shape = group)) +
    geom_point(size = 1.8, alpha = 0.75) +
    scale_colour_manual(values = c(Astral = "#1F4E79", timsTOF = "#D62828")) +
    theme_classic(base_size = 11) +
    ggtitle(ttl) +
    theme(legend.position = "bottom")
}
g1 <- make_pca(combined,   "Before ComBat")
g2 <- make_pca(combat_mat, "After ComBat")
ga <- gridExtra::arrangeGrob(g1, g2, ncol = 2)
ggsave(file.path(FIG_DIR, "CSF_combined_R_PCA.png"), ga,
       width = 13, height = 5.5, dpi = 200)
print(grid::grid.draw(ga))


# ---- limma DE on ComBat-corrected matrix ----
res <- moderated_t_safe(combat_mat, groups, "MS", "Control")
res <- res[order(res$P.Value), ]
res$is_cross_omics <- res$gene %in% CROSS_OMICS
res$is_paper_top   <- res$gene %in% PAPER_TOP
res$is_recurring   <- res$gene %in% RECURRING

cat(sprintf("Genes: %d  |  FDR<0.05: %d  |  FDR<0.001: %d\n",
            nrow(res), sum(res$adj.P.Val < 0.05),
            sum(res$adj.P.Val < 0.001)))

cat("\n=== Cross-omics (combined) ===\n")
print(subset(res, is_cross_omics)[, c("gene","logFC","P.Value","adj.P.Val")])

tsv_fp <- file.path(OUT_DIR, "CSF_combined_R_ComBat_DE.tsv")
write.table(res, tsv_fp, sep = "\t", row.names = FALSE, quote = FALSE)


p <- dep_volcano_gg(res,
       title = "CSF Astral + timsTOF — R/ComBat + limma meta",
       subtitle = sprintf("%d genes · %d MS / %d Ctrl · sva::ComBat(platform) + eBayes",
                          nrow(res), sum(groups=="MS"), sum(groups=="Control")))
ggsave(file.path(FIG_DIR, "CSF_combined_R_volcano.png"), p,
       width = 11, height = 8, dpi = 200)
print(p)

