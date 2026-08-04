#!/usr/bin/env Rscript
## 00_run_all.R — execute every notebook in order (R-script form)
##
## Run each .R produced by _build_ipynb.py end-to-end. Outputs go to
##   ../processed/META/      *.tsv
##   ../figures/             *.png + *.pdf
##
## Re-running is idempotent (overwrites).
##
## Notebook execution order matters because nb 09 + 10 consume
## results from 01-07.
##
## Usage:
##   cd Proteomics/r_notebooks
##   Rscript 00_run_all.R
##
## Or to execute the .ipynb forms with all stdout / figures captured:
##   PATH=$RENV/bin:$PATH jupyter nbconvert --to notebook --execute --inplace \
##       --ExecutePreprocessor.kernel_name=ir_methylation \
##       --ExecutePreprocessor.timeout=3600 \
##       0?_*.ipynb 10_*.ipynb

HERE <- local({
  # sys.frame(1)$ofile only exists under source(); this file's own header documents
  # `Rscript 00_run_all.R`, which would abort here. Resolve from --file= when run that way.
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else dirname(normalizePath(sys.frame(1)$ofile))
})
setwd(HERE)

SCRIPTS <- c(
  "08_per_group_consistency.R",   # no I/O dep, runs first
  "04cc_magliozzi_brain_completecase.R",  # local xlsx only; complete-case, no imputation

  "06_pegram_gse32915_de.R",      # GEOquery download GSE32915
  "05_t_lineage_meta.R",          # GEOquery GSE32915+GSE78244 (cache hits 06)
  "01cc_csf_astral_completecase.R",       # 79 MB matrix; complete-case, no imputation
  "02cc_csf_timstof_completecase.R",      # 203 MB matrix; complete-case, no imputation
  "03_csf_cross_platform_meta.R", # depends on Astral + timsTOF raw
  "07_brainwm_rna_meth_rerun.R",  # local CSV only
  "09_cross_assay_lxn.R",         # depends on 01-07 outputs
  "10_master_validation.R",       # depends on 09 summary
  "11_itgb2_csf_pleocytosis.R"    # CSF leukocyte-count model (Section 2.5)
)

for (s in SCRIPTS) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("RUN: ", s, "\n", sep = "")
  cat(strrep("=", 70), "\n\n", sep = "")
  t0 <- Sys.time()
  res <- try(source(s, echo = FALSE, max.deparse.length = 1e6), silent = FALSE)
  t1 <- Sys.time()
  cat(sprintf("\n>>> %s finished in %.1f s\n", s,
              as.numeric(difftime(t1, t0, units = "secs"))))
  if (inherits(res, "try-error")) {
    cat("!!! ERROR — stopping pipeline.\n")
    break
  }
}
cat("\nAll done.\n")

# The adapter is the only step that converts the complete-case tables into the RDEP_CC layout
# that Figures 4 and 6 read; without it a clean run leaves both figures without inputs.
adapter <- file.path(HERE, "build_RDEP_CC_adapters.py")
if (file.exists(adapter)) {
  message("Running build_RDEP_CC_adapters.py")
  st <- system2("__PYTHON_BIN__", adapter)
  if (st != 0) stop("build_RDEP_CC_adapters.py failed with status ", st)
} else stop("missing: ", adapter)
