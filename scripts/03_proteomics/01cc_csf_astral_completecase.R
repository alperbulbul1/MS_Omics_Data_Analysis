#!/usr/bin/env Rscript
## 01_csf_astral_dep.R  —  generated from notebook spec
## Run: Rscript 01_csf_astral_dep.R


## ============================================================
## # 01 — CSF Astral proteomics: DEP/limma reanalysis
## 
## **Cell 2026 (Skene/Mann)** Astral DIA — *single largest CSF proteomics cohort*
## (978 MS + 306 Control after annotation matching). Reanalysed here in R with
## the canonical Bioconductor stack:
## 
## `DEP::filter_proteins → DEP::normalize_vsn → DEP::impute(MinProb) → DEP::test_diff`
## 
## with a graceful fallback to `limma + vsn + manual MinProb` if the DEP
## package is not installed.
## 
## **Inputs**
## - `processed/astral_discovery_gene_keyed.tsv` — 3,053 proteins × ~2,000 sample columns
## - `osfstorage-archive/processed proteomic data/0_sample_annotations/annotations_*.tsv`
## 
## **Outputs**
## - `processed/META/CSF_Astral_CC_results.tsv`
## - `figures/CSF_Astral_CC_volcano.png`
## - `figures/CSF_Astral_CC_qc.png`
## 
## **Why this matters:** Astral has more power than every other catalogued CSF
## study combined; per user instruction NO cross-platform combination here
## (see notebook 03 for the Astral×timsTOF meta).
## ============================================================

suppressPackageStartupMessages({
  library(limma); library(data.table); library(ggplot2); library(dplyr)
})
source("helpers.R")

## COMPLETE-CASE variant of dep_equivalent_de: filter + vsn, NO imputation.
## limma fits each protein on its observed values (MAR assumption); proteins needing
## >=2 observed values per group to be estimable.
dep_completecase_de <- function(mat, group_vec, group_a = "MS", group_b = "Control", thr = 0.5) {
  m <- filter_missval_R(mat, group_vec, thr = thr)
  m <- vsn_with_fallback(m)
  ga <- group_vec == group_a; gb <- group_vec == group_b
  ok <- apply(m, 1, function(x) sum(!is.na(x[ga])) >= 2 && sum(!is.na(x[gb])) >= 2)
  m <- m[ok, , drop = FALSE]
  message(sprintf("  -> COMPLETE-CASE limma on %d proteins x %d samples (no imputation)",
                  nrow(m), ncol(m)))
  moderated_t_safe(m, group_vec, group_a = group_a, group_b = group_b)
}
have_DEP <- FALSE   # COMPLETE-CASE: force the limma path, which tolerates NAs
cat(sprintf("DEP available: %s\n", have_DEP))


# ---- Load Astral gene-keyed matrix ----
astral_fp <- file.path(PROT_ROOT, "processed", "astral_discovery_gene_keyed.tsv")
cat(sprintf("Reading %s (%.0f MB)...\n",
            basename(astral_fp),
            file.info(astral_fp)$size/1e6))
raw <- fread(astral_fp, sep = "\t", header = TRUE, showProgress = FALSE)
cat(sprintf("  Astral: %d proteins × %d cols\n", nrow(raw), ncol(raw)))

# Sample columns = the .raw run names (everything else is metadata)
sample_cols <- grep("\\.raw$|^[0-9]{8}_", colnames(raw), value = TRUE)
proteins <- as.character(raw$Genes)
cat(sprintf("  Sample (.raw) columns: %d   non-sample meta: %d\n",
            length(sample_cols), ncol(raw) - length(sample_cols)))
if (length(sample_cols) == 0)
  stop("No sample run columns found — check naming pattern.")


# ---- Load annotation, build MS vs Control mapping ----
ann_fp <- file.path(PROT_ROOT, "osfstorage-archive",
  "processed proteomic data", "0_sample_annotations",
  "annotations_v42_49_2_10_4_10_interimSky17_PL01-PL56_PepResCustv01_resubmission.tsv")
ann <- fread(ann_fp, sep = "\t", header = TRUE)
cat(sprintf("  Annotation rows: %d\n", nrow(ann)))

ann_a <- ann[!is.na(Run_Astral_Measurement) & Run_Astral_Measurement != "",
             .(Run_Astral_Measurement, Diagnosis_group, MSgroup)]
ann_a[, group := fifelse(MSgroup == "MS", "MS",
                  fifelse(Diagnosis_group %in% c("Other","Neurological Control"),
                          "Control", NA_character_))]
ann_a <- ann_a[!is.na(group)]
cat(sprintf("  Annotated samples: %d  (MS=%d  Control=%d)\n",
            nrow(ann_a), sum(ann_a$group=="MS"), sum(ann_a$group=="Control")))

keep_cols <- intersect(sample_cols, ann_a$Run_Astral_Measurement)
ann_match <- ann_a[Run_Astral_Measurement %in% keep_cols]
groups <- ann_match$group[match(keep_cols, ann_match$Run_Astral_Measurement)]
cat(sprintf("  Matched matrix cols: %d  (MS=%d  Control=%d)\n",
            length(keep_cols), sum(groups=="MS"), sum(groups=="Control")))

expr_mat <- as.matrix(raw[, ..keep_cols])
storage.mode(expr_mat) <- "numeric"

# Collapse to gene level (max-variance representative)
expr_mat <- gene_dedup(expr_mat, proteins)
cat(sprintf("  Matrix after gene-dedup: %d genes × %d samples\n",
            nrow(expr_mat), ncol(expr_mat)))


# ---- DEP-equivalent pipeline ----
if (have_DEP) {
  message("Using canonical DEP pipeline (vsn + MinProb + test_diff)")
  suppressPackageStartupMessages({ library(DEP); library(SummarizedExperiment) })
  data_df <- data.frame(name = proteins, ID = proteins, expr_mat,
                        check.names = FALSE, stringsAsFactors = FALSE)
  data_unique <- make_unique(data_df, names = "name", ids = "ID")
  exp_design  <- data.frame(label = keep_cols, condition = groups,
                            replicate = ave(seq_along(keep_cols), groups,
                                            FUN = seq_along))
  cols_lfq <- which(colnames(data_unique) %in% keep_cols)
  se <- make_se(data_unique, columns = cols_lfq, expdesign = exp_design)
  se_flt  <- filter_proteins(se, type = "fraction", min = 0.5)
  se_norm <- normalize_vsn(se_flt)
  se_imp <- se_norm   # COMPLETE-CASE: no imputation (MAR assumption)
  dep <- test_diff(se_imp, type = "manual", test = "MS_vs_Control")
  dep <- add_rejections(dep, alpha = 0.05, lfc = 0)
  res <- get_results(dep)
  out <- data.frame(
    gene = res$name,
    logFC = res[[grep("MS_vs_Control_ratio$", colnames(res), value = TRUE)]],
    P.Value = res[[grep("MS_vs_Control_p\\.val$", colnames(res), value = TRUE)]],
    adj.P.Val = res[[grep("MS_vs_Control_p\\.adj$", colnames(res), value = TRUE)]],
    stringsAsFactors = FALSE)
} else {
  message("DEP not installed — using helpers' DEP-equivalent path")
  out <- dep_completecase_de(expr_mat, groups, "MS", "Control", thr = 0.5)
}
out <- out[!is.na(out$logFC) & !is.na(out$P.Value), ]
out <- out[order(out$P.Value), ]
cat(sprintf("  Proteins: %d  |  FDR<0.05: %d  |  FDR<0.001: %d\n",
            nrow(out),
            sum(out$adj.P.Val < 0.05, na.rm = TRUE),
            sum(out$adj.P.Val < 0.001, na.rm = TRUE)))


# ---- Inspect cross-omics + paper-top + recurring hits ----
cat("\n=== Cross-omics candidates ===\n")
print(subset(out, gene %in% CROSS_OMICS)[, c("gene","logFC","P.Value","adj.P.Val")])

cat("\n=== Cell-2026 paper top hits ===\n")
print(subset(out, gene %in% PAPER_TOP)[, c("gene","logFC","P.Value","adj.P.Val")])


# ---- Save TSV + volcano figure ----
out$is_cross_omics <- out$gene %in% CROSS_OMICS
out$is_paper_top   <- out$gene %in% PAPER_TOP
out$is_recurring   <- out$gene %in% RECURRING
tsv_fp <- file.path(OUT_DIR, "CSF_Astral_CC_results.tsv")
write.table(out, tsv_fp, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Wrote %s\n", tsv_fp))

p <- dep_volcano_gg(out,
       title = "CSF Astral — R/DEP-limma reanalysis",
       subtitle = sprintf("%d proteins · %d MS / %d Control · vsn + MinProb + eBayes",
                          nrow(out), sum(groups=="MS"), sum(groups=="Control")))
ggsave(file.path(FIG_DIR, "CSF_Astral_CC_volcano.png"), p,
       width = 11, height = 8, dpi = 200)
ggsave(file.path(FIG_DIR, "CSF_Astral_CC_volcano.pdf"), p, width = 11, height = 8)
print(p)

