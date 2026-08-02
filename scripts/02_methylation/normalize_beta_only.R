#!/usr/bin/env Rscript
# =============================================================================
# normalize_beta_only.R
# =============================================================================
# Normalizes datasets that only provide beta-value matrices (no IDAT files)
# following the NBIS Array Tutorial approach:
#   https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/
#         tutorials/methylationArray/Array_Tutorial.html
#
# Strategy for beta-only datasets:
#  "preprocessQuantile is recommended for datasets without global differences
#   (e.g. blood). preprocessFunnorm for datasets with large-scale differences."
#  Since we cannot use minfi's preprocessors (they require RGChannelSet),
#  we use:
#   1. Within-dataset BMIQ normalization (wateRmelon::BMIQ) — corrects for
#      Type I vs Type II probe bias in beta values.
#      If wateRmelon is unavailable, falls back to quantile normalization.
#   2. Per-dataset detection p-value filtering is unavailable without IDAT;
#      instead filter probes with > 20% missing values.
#   3. Remove SNP-proximal probes (from annotation package).
#   4. Remove sex chromosome probes.
#   5. Convert normalized beta → M-values.
#
# Input:
#   - Combined_Methylation_Pre_Batch.csv  (probes × samples, beta or M values)
#   - Combined_Methylation_Metadata.csv
#   - Strict_Array_Preprocessed/ (IDAT-derived normalized matrices)
#
# Output (in Normalized_Beta_Only/):
#   - Per-dataset normalized beta + M CSV files
#   - Combined_BetaOnly_Beta.csv / Combined_BetaOnly_M.csv
#   - Combined_BetaOnly_Metadata.csv
#   - QC density plots before / after normalization
#
# After this script, the analysis pipeline merges IDAT and beta-only datasets.
# =============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})

# ── Args ──────────────────────────────────────────────────────────────────────
args       <- commandArgs(trailingOnly = TRUE)
base_dir   <- if (length(args) >= 1) args[[1]] else "__MS_GEO_ROOT__/Methylation_Data"
strict_dir <- if (length(args) >= 2) args[[2]] else file.path(base_dir, "Strict_Array_Preprocessed")
out_dir    <- if (length(args) >= 3) args[[3]] else file.path(base_dir, "Normalized_Beta_Only")

pre_batch_path <- file.path(base_dir, "Combined_Methylation_Pre_Batch.csv")
meta_path      <- file.path(base_dir, "Combined_Methylation_Metadata.csv")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) message(paste0("[normalize_beta_only] ", ...))

# ── Probe annotation (for SNP and sex chromosome filtering) ───────────────────
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Probes with SNPs at SBE or CpG (minor allele freq ANY threshold — be strict)
snp_probes <- rownames(ann450k)[!is.na(ann450k$Probe_rs) & nchar(ann450k$Probe_rs) > 0]
# Also include probes flagged in CpG SNP columns if available
if ("CpG_rs" %in% colnames(ann450k)) {
  snp_probes <- union(snp_probes,
                      rownames(ann450k)[!is.na(ann450k$CpG_rs) & nchar(ann450k$CpG_rs) > 0])
}
# Sex chromosome probes
sex_probes <- rownames(ann450k)[ann450k$chr %in% c("chrX", "chrY") & !is.na(ann450k$chr)]

msg("Annotation loaded: ", nrow(ann450k), " probes")
msg("  SNP probes to remove:        ", length(snp_probes))
msg("  Sex chromosome probes to remove: ", length(sex_probes))

# ── BMIQ normalization ────────────────────────────────────────────────────────
# Requires probe type information (Type I vs Type II) from annotation
probe_type <- rep("II", nrow(ann450k))
names(probe_type) <- rownames(ann450k)
if ("Type" %in% colnames(ann450k)) {
  probe_type[rownames(ann450k)] <- as.character(ann450k$Type)
}
type1_probes <- names(probe_type)[probe_type == "I"]
type2_probes <- names(probe_type)[probe_type == "II"]

has_bmiq <- requireNamespace("wateRmelon", quietly = TRUE)
msg(if (has_bmiq) "BMIQ normalization available (wateRmelon)" else "wateRmelon not found — using quantile normalization fallback")

normalize_beta_dataset <- function(beta_mat, dataset_id) {
  # beta_mat: probes × samples, values in [0,1]
  msg(dataset_id, ": normalizing ", nrow(beta_mat), " probes × ", ncol(beta_mat), " samples")

  if (has_bmiq) {
    # BMIQ: within-array normalization for Type I / Type II probe bias
    norm_mat <- tryCatch({
      # BMIQ operates on a named vector of probe types
      common_ann <- intersect(rownames(beta_mat), names(probe_type))
      if (length(common_ann) < 1000) {
        msg(dataset_id, ": too few annotated probes for BMIQ, falling back to quantile")
        stop("too few probes")
      }
      beta_sub <- beta_mat[common_ann, , drop = FALSE]
      pt       <- probe_type[common_ann]
      design_vec <- ifelse(pt == "I", 1L, 2L)

      result <- matrix(NA_real_, nrow = nrow(beta_sub), ncol = ncol(beta_sub),
                       dimnames = dimnames(beta_sub))
      for (j in seq_len(ncol(beta_sub))) {
        sample_beta <- beta_sub[, j]
        valid <- !is.na(sample_beta) & sample_beta > 0 & sample_beta < 1
        if (sum(valid) < 500) {
          result[, j] <- sample_beta
          next
        }
        norm_j <- tryCatch(
          wateRmelon::BMIQ(sample_beta[valid], design.v = design_vec[valid],
                           plots = FALSE, pri = FALSE)$nbeta,
          error = function(e) sample_beta[valid]
        )
        result[valid, j] <- norm_j
        result[!valid, j] <- NA_real_
      }
      beta_mat[common_ann, ] <- result
      beta_mat
    }, error = function(e) {
      msg(dataset_id, ": BMIQ failed (", conditionMessage(e), ") — using quantile")
      # Quantile normalization fallback
      tryCatch(
        preprocessCore::normalize.quantiles(beta_mat),
        error = function(e2) beta_mat
      )
    })
  } else {
    # Quantile normalization (limma approach: normalizeQuantiles)
    norm_mat <- tryCatch(
      limma::normalizeQuantiles(beta_mat),
      error = function(e) beta_mat
    )
  }

  # Clip to [0,1] after normalization
  norm_mat <- pmin(pmax(norm_mat, 0), 1)
  dimnames(norm_mat) <- dimnames(beta_mat)
  norm_mat
}

beta_to_m <- function(beta, eps = 1e-6) {
  b <- pmin(pmax(beta, eps), 1 - eps)
  log2(b / (1 - b))
}

# =============================================================================
# STEP 1: Identify which datasets have IDAT-preprocessed results
# =============================================================================
msg("STEP 1: Identifying beta-only datasets")

meta_all <- read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
meta_all$base_dataset <- sub("__.*$", "", as.character(meta_all$dataset))
meta_all <- meta_all[meta_all$condition %in% c("MS", "HC"), , drop = FALSE]
meta_all <- meta_all[!duplicated(meta_all$sample_id), , drop = FALSE]

# Datasets with IDAT preprocessing
idat_datasets <- character(0)
if (dir.exists(strict_dir)) {
  idat_datasets <- list.dirs(strict_dir, full.names = FALSE, recursive = FALSE)
  idat_datasets <- idat_datasets[nchar(idat_datasets) > 0]
}
msg("IDAT-preprocessed datasets: ", paste(idat_datasets, collapse = ", "))

# =============================================================================
# STEP 2: Load and parse Combined_Methylation_Pre_Batch.csv (beta or M values)
# =============================================================================
msg("STEP 2: Loading pre-batch matrix from Python pipeline")

if (!file.exists(pre_batch_path)) {
  msg("WARNING: Combined_Methylation_Pre_Batch.csv not found — nothing to process")
  quit(save = "no", status = 0)
}

smx_df <- read.csv(pre_batch_path, check.names = FALSE, stringsAsFactors = FALSE,
                   row.names = 1)
smx_mat <- as.matrix(smx_df)
storage.mode(smx_mat) <- "double"

# Determine whether values are beta [0,1] or M values
sample_min <- min(smx_mat, na.rm = TRUE)
sample_max <- max(smx_mat, na.rm = TRUE)
is_beta <- (sample_min >= -1e-4 && sample_max <= 1.0001)
msg(sprintf("Matrix range: [%.4f, %.4f] → interpreted as %s values",
            sample_min, sample_max, if (is_beta) "beta" else "M"))

# =============================================================================
# STEP 3: Per beta-only dataset processing
# =============================================================================
msg("STEP 3: Processing each beta-only dataset")

beta_only_results <- list()  # dataset_id → list(beta, mval, meta)
beta_only_summaries <- list()

all_datasets <- sort(unique(meta_all$base_dataset))
beta_only_datasets <- setdiff(all_datasets, idat_datasets)

msg("Beta-only datasets to process: ", paste(beta_only_datasets, collapse = ", "))
if (!length(beta_only_datasets)) {
  msg("No beta-only datasets found — all datasets have IDAT preprocessing. Exiting.")
  quit(save = "no", status = 0)
}

for (ds in beta_only_datasets) {
  ds_meta <- meta_all[meta_all$base_dataset == ds, , drop = FALSE]
  ds_samples <- intersect(ds_meta$sample_id, colnames(smx_mat))

  if (length(ds_samples) < 4) {
    msg(ds, ": only ", length(ds_samples), " sample(s) in matrix — skipping")
    beta_only_summaries[[ds]] <- data.frame(
      dataset = ds, status = "skipped", reason = "too few samples",
      n_ms = 0, n_hc = 0, probes = 0, stringsAsFactors = FALSE
    )
    next
  }

  ds_meta_filt <- ds_meta[ds_meta$sample_id %in% ds_samples, , drop = FALSE]
  n_ms <- sum(ds_meta_filt$condition == "MS")
  n_hc <- sum(ds_meta_filt$condition == "HC")
  if (n_ms < 2 || n_hc < 2) {
    msg(ds, ": insufficient MS/HC samples (MS=", n_ms, " HC=", n_hc, ") — skipping")
    beta_only_summaries[[ds]] <- data.frame(
      dataset = ds, status = "skipped", reason = "insufficient case/control",
      n_ms = n_ms, n_hc = n_hc, probes = 0, stringsAsFactors = FALSE
    )
    next
  }

  # Extract this dataset's matrix
  if (is_beta) {
    beta_raw <- smx_mat[, ds_samples, drop = FALSE]
  } else {
    # Convert M → beta for normalization, then back
    m_raw    <- smx_mat[, ds_samples, drop = FALSE]
    beta_raw <- 2^m_raw / (1 + 2^m_raw)
  }

  # ── Per-dataset output directory ──────────────────────────────────────────
  ds_out <- file.path(out_dir, ds)
  dir.create(ds_out, recursive = TRUE, showWarnings = FALSE)

  # ── Filter 1: Remove probes with > 20% missing values ────────────────────
  na_frac <- rowMeans(is.na(beta_raw))
  beta_raw <- beta_raw[na_frac <= 0.20, , drop = FALSE]
  msg(ds, ": after NA filter: ", nrow(beta_raw), " probes")

  # ── Filter 2: Remove SNP probes ───────────────────────────────────────────
  beta_raw <- beta_raw[!rownames(beta_raw) %in% snp_probes, , drop = FALSE]
  msg(ds, ": after SNP filter: ", nrow(beta_raw), " probes")

  # ── Filter 3: Remove sex chromosome probes ────────────────────────────────
  beta_raw <- beta_raw[!rownames(beta_raw) %in% sex_probes, , drop = FALSE]
  msg(ds, ": after sex-chr filter: ", nrow(beta_raw), " probes")

  if (nrow(beta_raw) < 1000) {
    msg(ds, ": too few probes after filtering — skipping")
    beta_only_summaries[[ds]] <- data.frame(
      dataset = ds, status = "skipped", reason = "too few probes after filter",
      n_ms = n_ms, n_hc = n_hc, probes = nrow(beta_raw), stringsAsFactors = FALSE
    )
    next
  }

  # ── Density plot BEFORE normalization ────────────────────────────────────
  png(file.path(ds_out, paste0(ds, "_density_before_norm.png")),
      width = 900, height = 600, res = 120)
  tryCatch({
    pal <- c(MS = "#E05252", HC = "#5278E0")
    col_vec <- pal[ds_meta_filt$condition[match(colnames(beta_raw), ds_meta_filt$sample_id)]]
    densityPlot(beta_raw, sampGroups = ds_meta_filt$condition[match(colnames(beta_raw), ds_meta_filt$sample_id)],
                main = paste0(ds, " – raw beta (before BMIQ)"), legend = FALSE)
    legend("topright", legend = names(pal), fill = pal, cex = 0.8)
  }, error = function(e) NULL)
  dev.off()

  # ── BMIQ / quantile normalization ────────────────────────────────────────
  beta_norm <- normalize_beta_dataset(beta_raw, ds)

  # ── Density plot AFTER normalization ─────────────────────────────────────
  png(file.path(ds_out, paste0(ds, "_density_after_norm.png")),
      width = 900, height = 600, res = 120)
  tryCatch({
    densityPlot(beta_norm, sampGroups = ds_meta_filt$condition[match(colnames(beta_norm), ds_meta_filt$sample_id)],
                main = paste0(ds, " – normalized beta (after BMIQ)"), legend = FALSE)
    legend("topright", legend = names(pal[1:2]), fill = pal[1:2], cex = 0.8)
  }, error = function(e) NULL)
  dev.off()

  # ── Convert to M-values ───────────────────────────────────────────────────
  mval_norm <- beta_to_m(beta_norm)

  # ── MDS plot after normalization ─────────────────────────────────────────
  png(file.path(ds_out, paste0(ds, "_MDS_normalized.png")),
      width = 900, height = 800, res = 120)
  tryCatch({
    pal <- c(MS = "#E05252", HC = "#5278E0")
    col_vec <- pal[ds_meta_filt$condition[match(colnames(mval_norm), ds_meta_filt$sample_id)]]
    plotMDS(mval_norm, top = 1000, gene.selection = "common",
            col = col_vec,
            main = paste0(ds, " – MDS (M-values, post-normalization)"))
    legend("bottomright", legend = names(pal), fill = pal, cex = 0.8)
  }, error = function(e) NULL)
  dev.off()

  # ── Save per-dataset results ──────────────────────────────────────────────
  write.csv(data.frame(Probe = rownames(beta_norm), beta_norm, check.names = FALSE),
            file.path(ds_out, "Normalized_Beta.csv"), row.names = FALSE)
  write.csv(data.frame(Probe = rownames(mval_norm), mval_norm, check.names = FALSE),
            file.path(ds_out, "Normalized_M.csv"), row.names = FALSE)
  write.csv(ds_meta_filt, file.path(ds_out, "Filtered_Metadata.csv"), row.names = FALSE)

  beta_only_results[[ds]] <- list(
    beta = beta_norm,
    mval = mval_norm,
    meta = ds_meta_filt[, intersect(c("sample_id", "condition", "cell_type",
                                       "raw_text", "dataset", "base_dataset"),
                                     colnames(ds_meta_filt)), drop = FALSE]
  )
  beta_only_summaries[[ds]] <- data.frame(
    dataset = ds, status = "ok", reason = "",
    n_ms = n_ms, n_hc = n_hc, probes = nrow(beta_norm), stringsAsFactors = FALSE
  )
  msg(ds, ": done — ", nrow(beta_norm), " probes, MS=", n_ms, " HC=", n_hc)
}

# =============================================================================
# STEP 4: Combine beta-only datasets + write combined files
# =============================================================================
msg("STEP 4: Combining beta-only datasets")

summary_df <- do.call(rbind, beta_only_summaries)
write.csv(summary_df, file.path(out_dir, "BetaOnly_QC_Summary.csv"), row.names = FALSE)

ok_results <- beta_only_results[summary_df$dataset[summary_df$status == "ok"]]
if (!length(ok_results)) {
  msg("No beta-only datasets passed QC — no combined output will be written")
  quit(save = "no", status = 0)
}

common_probes <- Reduce(intersect, lapply(ok_results, function(x) rownames(x$beta)))
msg("Common probes across beta-only datasets: ", length(common_probes))

if (!length(common_probes)) {
  msg("WARNING: no common probes across beta-only datasets")
  quit(save = "no", status = 0)
}

combined_beta <- do.call(cbind, lapply(ok_results, function(x)
  x$beta[common_probes, , drop = FALSE]))
combined_mval <- do.call(cbind, lapply(ok_results, function(x)
  x$mval[common_probes, , drop = FALSE]))
combined_meta <- do.call(rbind, lapply(ok_results, function(x) x$meta))
combined_meta <- combined_meta[match(colnames(combined_beta), combined_meta$sample_id), , drop = FALSE]

write.csv(data.frame(Probe = common_probes, combined_beta, check.names = FALSE),
          file.path(out_dir, "Combined_BetaOnly_Beta.csv"), row.names = FALSE)
write.csv(data.frame(Probe = common_probes, combined_mval, check.names = FALSE),
          file.path(out_dir, "Combined_BetaOnly_M.csv"), row.names = FALSE)
write.csv(combined_meta, file.path(out_dir, "Combined_BetaOnly_Metadata.csv"), row.names = FALSE)

msg("Beta-only normalization complete!")
msg("  Datasets processed: ", length(ok_results))
msg("  Total samples:      ", nrow(combined_meta))
msg("  Common probes:      ", length(common_probes))
msg("")
msg("Next step: run methylation_analysis_pipeline.R which merges IDAT + beta-only datasets")
