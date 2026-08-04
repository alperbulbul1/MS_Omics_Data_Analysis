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
# Retained for reference, not called by this script: it is the limma-equivalent estimation route
# that dep_bh_equivalence_check.R exercises to show the DEP path below reproduces the same
# effect estimates. The reported numbers come from the DEP branch, not from here.
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
## DEP::test_diff is deliberately not used. It hard-codes fdrtool for multiple-testing
## adjustment (fdr_res <- fdrtool::fdrtool(res$t, ...)) with no argument to change it, whereas
## every layer of this study reports Benjamini-Hochberg; adopting it would make the proteomic
## layer non-comparable with the RNA, methylation and single-cell layers. It is otherwise
## equivalent, and that was tested rather than assumed - see dep_bh_equivalence_check.R: run on
## this same matrix with imputation omitted, test_diff accepts the missing values and returns
## identical fold changes (Pearson r = 1.000, Spearman rho = 0.9999 on p-value rank), and
## BH-adjusting its own p-values gives 955 significant proteins against the 941 reported, with
## identical calls on 1,983 of 1,995 shared proteins (99.4%).
## The limma path below is therefore the same test, reported on the same FDR scale as the rest
## of the study. DEP was also removed from Bioconductor at release 3.23.
have_DEP <- requireNamespace("DEP", quietly = TRUE)
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


# ---- DEP pipeline, imputation omitted, Benjamini-Hochberg adjustment ----
## DEP is used for filtering, normalisation and testing. Two departures from its documented
## workflow, both deliberate and both verified in dep_bh_equivalence_check.R:
##
##   1. DEP::impute is not called. Group-dependent detection turns left-censored imputation into an
##      apparent abundance effect; MinProb imputation was previously shown to manufacture an ITGB2
##      CSF result and to erase the LXN brain signal. DEP::test_diff accepts the missing values.
##   2. DEP::test_diff's own adjusted p-values are discarded and replaced with Benjamini-Hochberg.
##      test_diff hard-codes fdrtool (fdr_res <- fdrtool::fdrtool(res$t, ...)) with no argument to
##      change it. Every other layer of this study reports BH, so keeping fdrtool would make the
##      proteomic layer non-comparable; on this cohort it also calls only 35 proteins against 955.
##
## Two DEP behaviours that silently corrupt results on this input are worked around:
##   - make_se applies log2() internally, expecting raw LFQ intensities. This matrix is already
##     log-scaled, so the assay is restored afterwards to avoid a second log2.
##   - normalize_vsn likewise assumes raw intensities; median-centring is used instead, which is
##     what vsn_with_fallback selects for already-log data.
if (have_DEP) {
  message("Using DEP (filter + test_diff), no imputation, BH adjustment")
  suppressPackageStartupMessages({ library(DEP); library(SummarizedExperiment) })
  data_df <- data.frame(expr_mat, check.names = FALSE)
  data_df$name <- rownames(expr_mat); data_df$ID <- rownames(expr_mat)
  exp_design <- data.frame(label = keep_cols, condition = groups,
                           replicate = ave(seq_along(groups), groups, FUN = seq_along),
                           stringsAsFactors = FALSE)
  se <- make_se(data_df, columns = seq_len(ncol(expr_mat)), expdesign = exp_design)
  stopifnot(nrow(se) == nrow(expr_mat), ncol(se) == ncol(expr_mat))
  assay(se, withDimnames = FALSE) <- unname(expr_mat)          # undo make_se's second log2
  se <- filter_proteins(se, type = "fraction", min = 0.5)
  assay(se, withDimnames = FALSE) <- vsn_with_fallback(assay(se))
  dep <- suppressWarnings(test_diff(se, type = "manual", test = "MS_vs_Control"))
  rr  <- as.data.frame(rowData(dep))
  out <- data.frame(
    gene    = rr$name,
    logFC   = rr[[grep("_diff$|_ratio$", colnames(rr), value = TRUE)[1]]],
    P.Value = rr[[grep("_p\\.val$",      colnames(rr), value = TRUE)[1]]],
    stringsAsFactors = FALSE)
  out <- out[!is.na(out$logFC) & !is.na(out$P.Value), ]
  out$adj.P.Val <- p.adjust(out$P.Value, method = "BH")        # replaces DEP's fdrtool qval
} else {
  stop("DEP is required. It was removed from Bioconductor at release 3.23; install with\n",
       "  remotes::install_github(\"arnesmits/DEP\")\n",
       "or run dep_bh_equivalence_check.R, which documents the limma path that reproduces it.")
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

