#!/usr/bin/env Rscript
## 00_run_all.R — execute every methylation notebook in order
##
## Outputs:
##   ../results/   per-stratum DMP + gene-level + mCSEA + RNA-vs-meth concordance
##   ../figures/   per-stratum volcanos + cross-stratum heatmap
##
## Run:
##   cd Methylation/r_notebooks
##   Rscript 00_run_all.R

HERE <- dirname(normalizePath(sys.frame(1)$ofile))
setwd(HERE)

SCRIPTS <- c(
  "01_tcells_meth_dmp.R",
  "02_wb_dmf_meth_dmp.R",
  "03_wb_ocrelizumab_meth_dmp.R",
  "04_tcells_remission_meth_dmp.R",
  "05_combined_meth_dmp.R",        # large — 317 MB matrix
  "06_mcsea_promoter_analysis.R",  # depends on 05 output
  "07_brainwm_rna_vs_meth.R",
  "08_cross_stratum_meth_master.R"
)

for (s in SCRIPTS) {
  cat("\n", strrep("=", 70), "\nRUN: ", s, "\n", strrep("=", 70), "\n\n", sep = "")
  t0 <- Sys.time()
  res <- try(source(s, echo = FALSE, max.deparse.length = 1e6))
  t1 <- Sys.time()
  cat(sprintf("\n>>> %s finished in %.1f s\n", s,
              as.numeric(difftime(t1, t0, units = "secs"))))
  if (inherits(res, "try-error")) {
    cat("!!! stopping pipeline.\n"); break
  }
}
cat("\nAll done.\n")
