#!/usr/bin/env Rscript
## 00_run_all.R — execute every transcriptome notebook in order
##
## Outputs:
##   ../results/   per-stratum DE TSVs + cross-stratum summary
##   ../figures/   per-stratum volcanos + cross-stratum heatmap
##
## Run:
##   cd Transcriptome/r_notebooks
##   Rscript 00_run_all.R
##
## Or execute the .ipynb form with output capture:
##   PATH=$RENV/bin:$PATH jupyter nbconvert --to notebook --execute --inplace \
##       --ExecutePreprocessor.kernel_name=ir_methylation \
##       --ExecutePreprocessor.timeout=900 0?_*.ipynb

HERE <- dirname(normalizePath(sys.frame(1)$ofile))
setwd(HERE)

SCRIPTS <- c(
  "01_pbmc_de.R",
  "02_tcells_de.R",
  "03_bcells_de.R",
  "04_brainwm_de.R",
  "05_whole_blood_de.R",
  "06_pbmc_ifnb_de.R",
  "07_total_combined_de.R",
  "08_cross_stratum_master.R"
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
