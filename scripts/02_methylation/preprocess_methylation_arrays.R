#!/usr/bin/env Rscript
# =============================================================================
# preprocess_methylation_arrays.R
# =============================================================================
# IDAT-based methylation preprocessing using minfi, following NBIS tutorial:
#   https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/
#         tutorials/methylationArray/Array_Tutorial.html
#
# Steps (per dataset):
#  1. Read IDAT files → RGChannelSet
#  2. QC: detection p-values (sample + probe level), plotQC, controlStripPlot
#  3. Normalize with preprocessIllumina (bg subtract + control normalization)
#  4. Filter: bead count, detection p-value, SNP probes, cross-reactive, sex chr
#  5. Save beta + M-value matrices and RDS per dataset
#  6. Merge across datasets → combined matrices (saved for analysis pipeline)
# =============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})

# ── Cross-reactive probe list ─────────────────────────────────────────────────
load_xreactive_probes <- function() {
  if (requireNamespace("maxprobes", quietly = TRUE)) {
    message("Using maxprobes package for cross-reactive probe list")
    probes <- tryCatch(
      maxprobes::xreactive_probes(array_type = "450K"),
      error = function(e) character(0)
    )
    if (length(probes)) return(probes)
  }
  # Fallback: SNP-in-probe column from 450K annotation
  ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  xr  <- rownames(ann)[!is.na(ann$Probe_rs) & nchar(ann$Probe_rs) > 0]
  message(sprintf("Loaded %d SNP-in-probe probes as cross-reactive proxy", length(xr)))
  return(xr)
}

# ── Args ──────────────────────────────────────────────────────────────────────
args      <- commandArgs(trailingOnly = TRUE)
meta_path <- if (length(args) >= 1) args[[1]] else "__MS_GEO_ROOT__/Methylation_Data/Combined_Methylation_Metadata.csv"
base_dir  <- if (length(args) >= 2) args[[2]] else "__MS_GEO_ROOT__/Methylation_Data"
out_dir   <- if (length(args) >= 3) args[[3]] else file.path(base_dir, "Strict_Array_Preprocessed")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Thresholds ────────────────────────────────────────────────────────────────
sample_fail_threshold <- 0.05   # remove samples where > 5 % probes fail
probe_detp_threshold  <- 0.01   # p-value cutoff per probe
bead_count_threshold  <- 3      # probes with < 3 beads in ANY sample → dropped

xreactive_probes <- load_xreactive_probes()

# ── Metadata ──────────────────────────────────────────────────────────────────
meta <- read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
meta$dataset      <- as.character(meta$dataset)
meta$base_dataset <- sub("__.*$", "", meta$dataset)

# ── Helper: find IDAT basenames ───────────────────────────────────────────────
find_idat_basenames <- function(dataset_dir, sample_ids) {
  files <- list.files(dataset_dir, pattern = "\\.idat(\\.gz)?$", full.names = TRUE)
  if (!length(files)) {
    return(data.frame(sample_id = character(), Basename = character(),
                      stringsAsFactors = FALSE))
  }
  basenames <- list()
  for (sample_id in sample_ids) {
    sample_files <- files[grepl(
      sprintf("^%s_.*_(Grn|Red)\\.idat(\\.gz)?$", sample_id),
      basename(files)
    )]
    if (!length(sample_files)) next
    root <- unique(sub("_(Grn|Red)\\.idat(\\.gz)?$", "", sample_files))
    if (length(root) != 1) next
    grn_exists <- any(grepl("_Grn\\.idat(\\.gz)?$", sample_files))
    red_exists <- any(grepl("_Red\\.idat(\\.gz)?$", sample_files))
    if (!grn_exists || !red_exists) next
    basenames[[sample_id]] <- root[[1]]
  }
  if (!length(basenames)) {
    return(data.frame(sample_id = character(), Basename = character(),
                      stringsAsFactors = FALSE))
  }
  data.frame(
    sample_id = names(basenames),
    Basename  = unname(unlist(basenames)),
    stringsAsFactors = FALSE
  )
}

# ── Helper: array type string ─────────────────────────────────────────────────
array_to_dmrcate <- function(annotation_array) {
  text <- tolower(paste(annotation_array, collapse = " "))
  if (grepl("epic", text, fixed = TRUE)) return("EPIC")
  return("450K")
}

# =============================================================================
# Per-dataset loop
# =============================================================================
dataset_results  <- list()
dataset_summaries <- list()

for (base_dataset in sort(unique(meta$base_dataset))) {
  dataset_meta <- meta[meta$base_dataset == base_dataset, , drop = FALSE]
  dataset_dir  <- file.path(base_dir, base_dataset)

  if (!dir.exists(dataset_dir)) {
    dataset_summaries[[base_dataset]] <- data.frame(
      dataset = base_dataset, status = "skipped",
      reason = "dataset directory missing",
      samples_requested = nrow(dataset_meta), samples_with_idat = 0,
      samples_after_qc = 0, probes_after_filter = 0, arraytype = NA_character_,
      stringsAsFactors = FALSE
    )
    next
  }

  targets <- merge(
    dataset_meta,
    find_idat_basenames(dataset_dir, dataset_meta$sample_id),
    by = "sample_id", all.x = FALSE, all.y = FALSE
  )
  targets <- targets[!duplicated(targets$sample_id), , drop = FALSE]

  if (!nrow(targets)) {
    dataset_summaries[[base_dataset]] <- data.frame(
      dataset = base_dataset, status = "skipped",
      reason = "no matched idat basenames",
      samples_requested = nrow(dataset_meta), samples_with_idat = 0,
      samples_after_qc = 0, probes_after_filter = 0, arraytype = NA_character_,
      stringsAsFactors = FALSE
    )
    next
  }

  message(sprintf("[%s] reading %d IDAT samples", base_dataset, nrow(targets)))
  rgset <- tryCatch(
    read.metharray(as.character(targets$Basename), force = TRUE, extended = TRUE, verbose = FALSE),
    error = function(e) e
  )
  if (inherits(rgset, "error")) {
    dataset_summaries[[base_dataset]] <- data.frame(
      dataset = base_dataset, status = "failed",
      reason = conditionMessage(rgset),
      samples_requested = nrow(dataset_meta), samples_with_idat = nrow(targets),
      samples_after_qc = 0, probes_after_filter = 0, arraytype = NA_character_,
      stringsAsFactors = FALSE
    )
    next
  }
  sampleNames(rgset) <- targets$sample_id

  # ── QC 1: Detection p-values (NBIS tutorial approach) ─────────────────────
  det_p <- detectionP(rgset)

  # Bar chart of mean detection p-values per sample
  dataset_out <- file.path(out_dir, base_dataset)
  dir.create(dataset_out, recursive = TRUE, showWarnings = FALSE)

  png(file.path(dataset_out, paste0(base_dataset, "_sample_detection_pval.png")),
      width = 900, height = 600, res = 120)
  tryCatch({
    barplot(colMeans(det_p), las = 2, cex.names = 0.7,
            ylab = "Mean detection p-value",
            main = paste0(base_dataset, " – sample detection p-values"))
    abline(h = 0.05, col = "red", lty = 2)
  }, error = function(e) NULL)
  dev.off()

  # ── QC 2: plotQC (MethylSet log-median M vs U intensities) ────────────────
  mset_raw <- preprocessRaw(rgset)
  png(file.path(dataset_out, paste0(base_dataset, "_plotQC.png")),
      width = 800, height = 700, res = 120)
  tryCatch({
    qc <- getQC(mset_raw)
    plotQC(qc)
    title(main = paste0(base_dataset, " – sample QC"))
  }, error = function(e) NULL)
  dev.off()

  # ── QC 3: Bisulfite conversion control strip ───────────────────────────────
  png(file.path(dataset_out, paste0(base_dataset, "_bisulfite_control.png")),
      width = 1200, height = 600, res = 120)
  tryCatch({
    controlStripPlot(rgset, controls = "BISULFITE CONVERSION II")
  }, error = function(e) NULL)
  dev.off()

  # ── QC 4: Density plot BEFORE normalization ────────────────────────────────
  png(file.path(dataset_out, paste0(base_dataset, "_density_raw.png")),
      width = 900, height = 600, res = 120)
  tryCatch({
    raw_beta_pre <- getBeta(mset_raw)
    colnames(raw_beta_pre) <- targets$sample_id
    densityPlot(raw_beta_pre,
                sampGroups = targets$condition,
                main = paste0(base_dataset, " – raw beta"),
                legend = FALSE)
    legend("topright",
           legend = levels(factor(targets$condition)),
           col    = 1:2, lwd = 2, cex = 0.8)
    rm(raw_beta_pre)
  }, error = function(e) NULL)
  dev.off()
  rm(mset_raw)

  # ── Remove failed samples (> sample_fail_threshold probes fail) ────────────
  sample_fail_fraction <- colMeans(det_p > probe_detp_threshold, na.rm = TRUE)
  keep_samples <- sample_fail_fraction <= sample_fail_threshold
  if (!any(keep_samples)) {
    dataset_summaries[[base_dataset]] <- data.frame(
      dataset = base_dataset, status = "failed",
      reason = "all samples failed detection p-value QC",
      samples_requested = nrow(dataset_meta), samples_with_idat = nrow(targets),
      samples_after_qc = 0, probes_after_filter = 0, arraytype = NA_character_,
      stringsAsFactors = FALSE
    )
    next
  }

  n_removed_samples <- sum(!keep_samples)
  if (n_removed_samples > 0) {
    message(sprintf("[%s] removing %d failed samples: %s", base_dataset,
                    n_removed_samples,
                    paste(colnames(det_p)[!keep_samples], collapse = ", ")))
  }
  rgset  <- rgset[, keep_samples]
  det_p  <- det_p[, keep_samples, drop = FALSE]
  targets <- targets[keep_samples, , drop = FALSE]

  # ── Normalize with preprocessIllumina (NBIS tutorial recommended) ─────────
  # bg.correct: background subtraction (as in GenomeStudio)
  # normalize="controls": control-probe normalization
  message(sprintf("[%s] preprocessIllumina on %d retained samples", base_dataset, nrow(targets)))
  mset  <- preprocessIllumina(rgset, bg.correct = TRUE, normalize = "controls", reference = 1)
  gmset <- mapToGenome(mset)

  # ── Filter 1: Bead count (probes with < 3 beads in ANY sample) ────────────
  bc <- tryCatch(getNBeads(rgset), error = function(e) NULL)
  if (!is.null(bc)) {
    fail_bc <- rownames(bc)[
      rowSums(bc[rownames(bc) %in% rownames(gmset), , drop = FALSE] < bead_count_threshold,
              na.rm = TRUE) > 0
    ]
    gmset <- gmset[!rownames(gmset) %in% fail_bc, ]
    message(sprintf("[%s] bead count filter: removed %d probes with < %d beads",
                    base_dataset, length(fail_bc), bead_count_threshold))
  }

  # ── Filter 2: Detection p-value probe filter (ALL samples must pass) ───────
  # Align det_p to current gmset probes
  det_p_sub <- det_p[rownames(det_p) %in% rownames(gmset), , drop = FALSE]
  keep_probe_ids <- rownames(det_p_sub)[
    rowSums(det_p_sub < probe_detp_threshold, na.rm = TRUE) == ncol(det_p_sub)
  ]
  n_before <- nrow(gmset)
  gmset <- gmset[rownames(gmset) %in% keep_probe_ids, ]
  message(sprintf("[%s] detection p-value filter: removed %d probes",
                  base_dataset, n_before - nrow(gmset)))

  # ── Filter 3: SNP probes (SBE and CpG sites, MAF > 0.01) ─────────────────
  gmset <- dropLociWithSnps(gmset, snps = c("SBE", "CpG"), maf = 0.01)
  message(sprintf("[%s] after SNP filter: %d probes", base_dataset, nrow(gmset)))

  # ── Filter 4: Cross-reactive probes ───────────────────────────────────────
  if (length(xreactive_probes)) {
    n_before_xr <- nrow(gmset)
    gmset <- gmset[!rownames(gmset) %in% xreactive_probes, ]
    message(sprintf("[%s] cross-reactive filter: removed %d probes",
                    base_dataset, n_before_xr - nrow(gmset)))
  }

  # ── Filter 5: Sex chromosomes ──────────────────────────────────────────────
  annotation_df <- getAnnotation(gmset)
  autosomal <- !(annotation_df$chr %in% c("chrX", "chrY"))
  autosomal[is.na(autosomal)] <- FALSE
  gmset <- gmset[autosomal, ]
  annotation_df <- getAnnotation(gmset)
  message(sprintf("[%s] after sex-chr filter: %d probes", base_dataset, nrow(gmset)))

  # ── Extract values ─────────────────────────────────────────────────────────
  beta <- getBeta(gmset)
  mval <- getM(gmset)
  colnames(beta) <- targets$sample_id
  colnames(mval) <- targets$sample_id

  # ── QC 5: Density plot AFTER normalization ────────────────────────────────
  png(file.path(dataset_out, paste0(base_dataset, "_density_normalized.png")),
      width = 900, height = 600, res = 120)
  tryCatch({
    densityPlot(beta, sampGroups = targets$condition,
                main = paste0(base_dataset, " – normalized beta"),
                legend = FALSE)
    legend("topright",
           legend = levels(factor(targets$condition)),
           col = 1:2, lwd = 2, cex = 0.8)
  }, error = function(e) NULL)
  dev.off()

  # ── QC 6: MDS plot (limma::plotMDS) after filtering ───────────────────────
  png(file.path(dataset_out, paste0(base_dataset, "_MDS_filtered.png")),
      width = 900, height = 800, res = 120)
  tryCatch({
    pal <- c(MS = "#E05252", HC = "#5278E0")
    col_vec <- pal[targets$condition]
    plotMDS(mval, top = 1000, gene.selection = "common",
            col = col_vec,
            main = paste0(base_dataset, " – MDS (M-values, post-filter)"),
            xlab = "MDS1", ylab = "MDS2")
    legend("bottomright", legend = names(pal), fill = pal, cex = 0.8)
  }, error = function(e) NULL)
  dev.off()

  # ── Determine array type ───────────────────────────────────────────────────
  arraytype <- array_to_dmrcate(annotation(mset)[["array"]])

  # ── Save per-dataset outputs ───────────────────────────────────────────────
  write.csv(targets,
            file.path(dataset_out, "Filtered_Metadata.csv"), row.names = FALSE)
  write.csv(data.frame(Probe = rownames(beta),  beta, check.names = FALSE),
            file.path(dataset_out, "Normalized_Beta.csv"), row.names = FALSE)
  write.csv(data.frame(Probe = rownames(mval),  mval, check.names = FALSE),
            file.path(dataset_out, "Normalized_M.csv"),    row.names = FALSE)
  write.csv(
    data.frame(
      sample_id            = colnames(det_p),
      failed_probe_fraction = sample_fail_fraction[match(colnames(det_p), names(sample_fail_fraction))],
      stringsAsFactors = FALSE
    ),
    file.path(dataset_out, "Sample_QC.csv"), row.names = FALSE
  )
  write.csv(
    data.frame(Probe = rownames(annotation_df),
               chr   = annotation_df$chr,
               pos   = annotation_df$pos,
               stringsAsFactors = FALSE),
    file.path(dataset_out, "Probe_Annotation.csv"), row.names = FALSE
  )

  # Save filtered GenomicMethylSet as RDS for downstream use
  saveRDS(gmset, file.path(dataset_out, "Filtered_gmset.rds"))

  dataset_results[[base_dataset]] <- list(
    beta      = beta,
    mval      = mval,
    meta      = targets[, c("sample_id", "condition", "cell_type", "raw_text",
                             "dataset", "base_dataset"), drop = FALSE],
    arraytype = arraytype
  )
  dataset_summaries[[base_dataset]] <- data.frame(
    dataset            = base_dataset,
    status             = "ok",
    reason             = "",
    samples_requested  = nrow(dataset_meta),
    samples_with_idat  = nrow(find_idat_basenames(dataset_dir, dataset_meta$sample_id)),
    samples_after_qc   = nrow(targets),
    probes_after_filter = nrow(beta),
    arraytype          = arraytype,
    stringsAsFactors   = FALSE
  )
}

# =============================================================================
# Merge and save combined matrices
# =============================================================================
summary_df <- do.call(rbind, dataset_summaries)
write.csv(summary_df, file.path(out_dir, "Dataset_QC_Summary.csv"), row.names = FALSE)

successful <- dataset_results[summary_df$dataset[summary_df$status == "ok"]]
if (!length(successful)) stop("No datasets passed IDAT preprocessing")

common_probes <- Reduce(intersect, lapply(successful, function(x) rownames(x$beta)))
if (!length(common_probes)) stop("No common probes remained after preprocessing")

combined_beta <- do.call(cbind, lapply(successful, function(x)
  x$beta[common_probes, x$meta$sample_id, drop = FALSE]))
combined_m    <- do.call(cbind, lapply(successful, function(x)
  x$mval[common_probes, x$meta$sample_id, drop = FALSE]))
combined_meta <- do.call(rbind, lapply(successful, function(x) x$meta))
combined_meta <- combined_meta[match(colnames(combined_m), combined_meta$sample_id), , drop = FALSE]

# Batch correction (for visualisation / downstream merging reference)
batch_method <- "none"
combined_batch_corrected <- combined_m
if (length(unique(combined_meta$base_dataset)) > 1) {
  combined_batch_corrected <- tryCatch(
    removeBatchEffect(
      combined_m,
      batch  = combined_meta$base_dataset,
      design = model.matrix(~ factor(combined_meta$condition, levels = c("HC", "MS")))
    ),
    error = function(e) {
      message("removeBatchEffect failed: ", conditionMessage(e))
      combined_m
    }
  )
  batch_method <- "removeBatchEffect"
}

write.csv(
  data.frame(Probe = common_probes, combined_beta, check.names = FALSE),
  file.path(out_dir, "Combined_Methylation_Strict_Beta.csv"), row.names = FALSE
)
write.csv(
  data.frame(Probe = common_probes, combined_m, check.names = FALSE),
  file.path(out_dir, "Combined_Methylation_Strict_M.csv"), row.names = FALSE
)
write.csv(
  data.frame(Probe = common_probes, combined_batch_corrected, check.names = FALSE),
  file.path(out_dir, "Combined_Methylation_Strict_M_BatchCorrected.csv"), row.names = FALSE
)
write.csv(
  combined_meta[, c("sample_id", "condition", "cell_type", "raw_text",
                     "dataset", "base_dataset"), drop = FALSE],
  file.path(out_dir, "Combined_Methylation_Strict_Metadata.csv"), row.names = FALSE
)
write.csv(
  data.frame(
    metric = c("datasets_included", "samples_total", "ms_samples", "hc_samples",
               "common_probes", "combined_arraytype", "batch_method"),
    value  = c(length(successful), nrow(combined_meta),
               sum(combined_meta$condition == "MS"),
               sum(combined_meta$condition == "HC"),
               length(common_probes), "450K", batch_method),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "Combined_Summary.csv"), row.names = FALSE
)

message(sprintf(
  "IDAT preprocessing complete: %d datasets, %d samples, %d shared probes",
  length(successful), nrow(combined_meta), length(common_probes)
))
