#!/usr/bin/env Rscript
## dep_bh_equivalence_check.R
## =========================
## Shows that the complete-case limma path used for the reported CSF proteomics is equivalent to
## running DEP itself, without imputation and with Benjamini-Hochberg adjustment.
##
## WHY THIS EXISTS. Methods states that protein intensities were tested with limma rather than with
## DEP, and a reader is entitled to ask whether that choice changed anything. It did not. This
## script runs DEP on the same matrix, skipping DEP::impute, and compares the two results directly.
##
## WHAT IT ESTABLISHES (CSF Astral, 1,284 samples, 978 MS / 306 control)
##
##   DEP::test_diff runs on non-imputed data. It emits "Missing values in 'se_n'" and completes.
##   Fold changes are identical to the reported analysis:      Pearson r = 1.000
##   The p-value ranking is identical:                         Spearman rho = 0.9999
##   Significance calls agree on 1,983 of 1,995 shared proteins (99.4%).
##
##   The only material difference is the multiple-testing procedure. DEP::test_diff has no argument
##   for it; the method is hard-coded in the function body as
##
##       fdr_res <- fdrtool::fdrtool(res$t, plot = FALSE, verbose = FALSE)
##
##   fdrtool estimates an empirical null from the t-statistics. In this cohort - unbalanced at
##   978 MS versus 306 controls, with widespread real signal - that null comes out wide and only
##   35 proteins survive. Re-adjusting DEP's own p-values with BH gives 955, against the 941
##   reported. Every layer of this study reports BH-FDR, so adopting DEP's default would have made
##   the proteomic layer non-comparable with the RNA, methylation and single-cell layers.
##
## TWO PITFALLS, both of which silently corrupt the comparison if not handled:
##
##   1. DEP::make_se applies log2() internally because it expects raw LFQ intensities. This matrix
##      is already log-scaled (median 19.9), so make_se log2-transforms it a second time and
##      compresses every fold change by roughly a factor of ten. The assay is restored after
##      make_se.
##   2. DEP::normalize_vsn likewise assumes raw intensities. Applied to already-log data it
##      compresses fold changes by about 4,000x. The same median-centring the reported pipeline
##      uses is applied instead.
##
## DEP was removed from Bioconductor at release 3.23. Install from source:
##   remotes::install_github("arnesmits/DEP")

suppressPackageStartupMessages({
  library(data.table); library(limma)
})

ROOT <- "__MS_GEO_ROOT__"
PROT <- file.path(ROOT, "Proteomics")
source("helpers.R")   # layer convention: run from scripts/03_proteomics/

if (!requireNamespace("DEP", quietly = TRUE)) {
  cat("DEP is not installed. It was removed from Bioconductor at release 3.23.\n",
      "Install with: remotes::install_github(\"arnesmits/DEP\")\n", sep = "")
  quit(save = "no", status = 0)
}
suppressPackageStartupMessages({ library(DEP); library(SummarizedExperiment) })

## ---- data prep: identical to 01cc_csf_astral_completecase.R ----
raw <- fread(file.path(PROT, "processed", "astral_discovery_gene_keyed.tsv"),
             sep = "\t", header = TRUE, showProgress = FALSE)
sample_cols <- grep("\\.raw$|^[0-9]{8}_", colnames(raw), value = TRUE)
proteins <- as.character(raw$Genes)

ann <- fread(file.path(PROT, "osfstorage-archive", "processed proteomic data",
                       "0_sample_annotations",
  "annotations_v42_49_2_10_4_10_interimSky17_PL01-PL56_PepResCustv01_resubmission.tsv"),
  sep = "\t", header = TRUE)
ann_a <- ann[!is.na(Run_Astral_Measurement) & Run_Astral_Measurement != "",
             .(Run_Astral_Measurement, Diagnosis_group, MSgroup)]
ann_a[, group := fifelse(MSgroup == "MS", "MS",
                  fifelse(Diagnosis_group %in% c("Other", "Neurological Control"),
                          "Control", NA_character_))]
ann_a <- ann_a[!is.na(group)]
keep_cols <- intersect(sample_cols, ann_a$Run_Astral_Measurement)
ann_match <- ann_a[Run_Astral_Measurement %in% keep_cols]
groups <- ann_match$group[match(keep_cols, ann_match$Run_Astral_Measurement)]

expr_mat <- as.matrix(raw[, ..keep_cols]); storage.mode(expr_mat) <- "numeric"
expr_mat <- gene_dedup(expr_mat, proteins)
cat(sprintf("matrix: %d genes x %d samples (MS=%d Control=%d)\n",
            nrow(expr_mat), ncol(expr_mat),
            sum(groups == "MS"), sum(groups == "Control")))

## ---- DEP, with impute() omitted ----
df <- data.frame(expr_mat, check.names = FALSE)
df$name <- rownames(expr_mat); df$ID <- rownames(expr_mat)
exp_design <- data.frame(label = keep_cols, condition = groups,
                         replicate = ave(seq_along(groups), groups, FUN = seq_along),
                         stringsAsFactors = FALSE)
se <- make_se(df, columns = seq_len(ncol(expr_mat)), expdesign = exp_design)
stopifnot(nrow(se) == nrow(expr_mat), ncol(se) == ncol(expr_mat))
assay(se, withDimnames = FALSE) <- unname(expr_mat)     # undo make_se's second log2
se <- filter_proteins(se, type = "fraction", min = 0.5)
assay(se, withDimnames = FALSE) <- vsn_with_fallback(assay(se))   # median-centre, not vsn

res <- suppressWarnings(test_diff(se, type = "manual", test = "MS_vs_Control"))
rr <- as.data.frame(rowData(res))
out <- data.table(
  gene    = rr$name,
  logFC   = rr[[grep("_diff$|_ratio$", colnames(rr), value = TRUE)[1]]],
  P.Value = rr[[grep("_p\\.val$",      colnames(rr), value = TRUE)[1]]],
  FDR_dep = rr[[grep("_p\\.adj$",      colnames(rr), value = TRUE)[1]]])
out <- out[!is.na(logFC) & !is.na(P.Value)]
out[, FDR_BH := p.adjust(P.Value, method = "BH")]

## ---- compare with the reported complete-case limma result ----
cc <- fread(file.path(PROT, "processed", "META", "CSF_Astral_CC_results.tsv"))
m  <- merge(out, cc, by = "gene")

cat(sprintf("\n%-34s %8s\n", "significant at FDR < 0.05", "n"))
cat(paste(rep("-", 44), collapse = ""), "\n")
cat(sprintf("%-34s %8d\n", "DEP, its own fdrtool adjustment", sum(out$FDR_dep < 0.05, na.rm = TRUE)))
cat(sprintf("%-34s %8d\n", "DEP, BH-adjusted",               sum(out$FDR_BH  < 0.05, na.rm = TRUE)))
cat(sprintf("%-34s %8d\n", "reported (limma + BH)",          sum(cc$adj.P.Val < 0.05, na.rm = TRUE)))

cat(sprintf("\nshared proteins: %d\n", nrow(m)))
cat(sprintf("  logFC Pearson r          = %.5f\n", cor(m$logFC.x, m$logFC.y, use = "complete.obs")))
cat(sprintf("  BH-FDR Spearman rho      = %.5f\n",
            cor(m$FDR_BH, m$adj.P.Val, method = "spearman", use = "complete.obs")))
cat(sprintf("  identical calls at 0.05  = %d / %d (%.1f%%)\n",
            sum((m$FDR_BH < 0.05) == (m$adj.P.Val < 0.05), na.rm = TRUE), nrow(m),
            100 * mean((m$FDR_BH < 0.05) == (m$adj.P.Val < 0.05), na.rm = TRUE)))

G <- c("ITGB2", "IKZF1", "CD79B", "LXN", "HLA-E", "CTSZ", "CHL1", "ICAM1", "ITGAL", "SH3BP4")
cat(sprintf("\n  %-9s %21s %25s\n", "gene", "DEP + BH", "reported (limma + BH)"))
for (g in G) {
  r <- m[gene == g]
  if (nrow(r)) cat(sprintf("  %-9s %+8.3f  FDR %-9.2g %+8.3f  FDR %.2g\n",
                            g, r$logFC.x, r$FDR_BH, r$logFC.y, r$adj.P.Val))
}

fwrite(out, file.path(PROT, "processed", "META", "CSF_Astral_DEP_BH_check.tsv"), sep = "\t")
cat(sprintf("\nWrote %s\n",
            file.path(PROT, "processed", "META", "CSF_Astral_DEP_BH_check.tsv")))
