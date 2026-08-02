#!/usr/bin/env Rscript
# =============================================================================
# run_methylation_subgroup_limma.R
# =============================================================================
# Full per-subgroup methylation analysis following the NBIS Array Tutorial:
#   https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/
#         tutorials/methylationArray/Array_Tutorial.html
#
# Usage:
#   Rscript run_methylation_subgroup_limma.R \
#       <meta_csv> <m_matrix_csv> <out_dir> \
#       [beta_matrix_csv] [arraytype]
#
# Arguments:
#   meta_csv        CSV with columns: sample_id, condition, base_dataset (or dataset)
#   m_matrix_csv    Probes x samples M-value matrix (first col = probe IDs)
#   out_dir         Output directory
#   beta_matrix_csv (optional) Beta-value matrix for plotCpg / mCSEA
#   arraytype       (optional) "450K" (default) or "EPIC"
#
# Steps (NBIS tutorial approach):
#  1. Load and QC M-values / beta-values
#  2. ComBat / removeBatchEffect for QC visualisation (MDS plots)
#  3. DMP — lmFit + eBayes with batch as design covariate (NOT pre-corrected)
#     — annotated with 450K annotation package
#     — plotCpg for top significant probes
#  4. DMR — cpg.annotate → dmrcate → extractRanges → DMR.plot
#  5. mCSEA — named rank vector from DMP logFC (NBIS approach)
#     — mCSEATest (promoters + genes) → mCSEAPlot for top genes
#  6. GO / KEGG enrichment — missMethyl::gometh
# =============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(minfi)
  library(DMRcate)
  library(mCSEA)
  library(data.table)
  library(matrixStats)
})

# ── Performance: detect cores for mCSEA parallelisation ───────────────────────
nCores <- min(4L, parallel::detectCores(logical = FALSE))
msg_start <- proc.time()

# ── Parse arguments ───────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_methylation_subgroup_limma.R <meta_csv> <m_matrix_csv> <out_dir> [beta_matrix_csv] [arraytype]")
}

meta_path   <- args[[1]]
matrix_path <- args[[2]]
out_dir     <- args[[3]]
beta_path   <- if (length(args) >= 4 && nzchar(args[[4]])) args[[4]] else
               sub("_M(\\.csv)$", "_Beta\\1", matrix_path, perl = TRUE)
arraytype   <- if (length(args) >= 5 && nzchar(args[[5]])) args[[5]] else "450K"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) message(paste0("[methylation_subgroup] ", ...))

# ── Annotation package ────────────────────────────────────────────────────────
ann_pkg <- tryCatch({
  if (toupper(arraytype) == "EPIC") {
    requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)
    getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  } else {
    requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)
    getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19)
  }
}, error = function(e) NULL)

ann450kSub <- NULL
if (!is.null(ann_pkg)) {
  ann450kSub <- as.data.frame(ann_pkg[, c("chr", "pos", "strand",
                                            "Name", "UCSC_RefGene_Name",
                                            "Relation_to_Island"), drop = FALSE],
                               stringsAsFactors = FALSE)
  ann450kSub$Probe <- rownames(ann450kSub)
}

# =============================================================================
# STEP 1: Load data
# =============================================================================
msg("STEP 1: Loading data (optimised with data.table::fread)")
t1 <- proc.time()

meta <- fread(meta_path, data.table = FALSE)
meta <- meta[meta$condition %in% c("MS", "HC"), , drop = FALSE]
dataset_col <- if ("base_dataset" %in% names(meta)) "base_dataset" else "dataset"

# ── Fast load M-value matrix ──────────────────────────────────────────────────
msg("  Reading M-value matrix...")
mat_dt <- fread(matrix_path, data.table = FALSE)
probe_ids <- mat_dt[[1]]
mat <- as.matrix(mat_dt[, -1, drop = FALSE])
rownames(mat) <- probe_ids
storage.mode(mat) <- "double"
rm(mat_dt); gc(verbose = FALSE)
msg("  M-values: ", nrow(mat), " probes x ", ncol(mat), " columns")

# ── Align samples BEFORE loading beta (only load needed columns) ──────────────
samples <- meta$sample_id[meta$sample_id %in% colnames(mat)]
meta    <- meta[match(samples, meta$sample_id), , drop = FALSE]
mat     <- mat[, samples, drop = FALSE]

# ── Vectorised variance filter (matrixStats is 10-100x faster than apply) ─────
msg("  Variance filtering...")
row_var <- rowVars(mat, na.rm = TRUE)
keep <- is.finite(row_var) & row_var > 0
mat <- mat[keep, , drop = FALSE]
row_var <- row_var[keep]
msg("  After variance filter: ", nrow(mat), " probes")

# ── Top-variance probe prefilter for speed (max 450K probes) ──────────────────
MAX_PROBES <- 450000L
if (nrow(mat) > MAX_PROBES) {
  top_idx <- order(row_var, decreasing = TRUE)[seq_len(MAX_PROBES)]
  mat <- mat[top_idx, , drop = FALSE]
  msg("  Top-variance prefilter: kept ", MAX_PROBES, " probes")
}

# ── Sample count check ────────────────────────────────────────────────────────
n_ms <- sum(meta$condition == "MS")
n_hc <- sum(meta$condition == "HC")
if (nrow(meta) < 8 || n_ms < 4 || n_hc < 4) {
  stop(sprintf("Not enough samples (MS=%d, HC=%d, need >= 4 each)", n_ms, n_hc))
}

msg("Samples: ", nrow(meta), " (MS=", n_ms, ", HC=", n_hc, ")")
msg("Probes:  ", nrow(mat))
msg("Datasets: ", paste(sort(unique(meta[[dataset_col]])), collapse = ", "))

# ── Vectorised NA → row median replacement (no for-loop) ──────────────────────
msg("  Replacing NAs with row medians...")
mat_nona <- mat
na_mask <- is.na(mat_nona)
if (any(na_mask)) {
  row_meds <- rowMedians(mat_nona, na.rm = TRUE)
  for (j in seq_len(ncol(mat_nona))) {
    col_nas <- na_mask[, j]
    if (any(col_nas)) mat_nona[col_nas, j] <- row_meds[col_nas]
  }
}

# ── Load beta matrix (only needed columns + probes) ───────────────────────────
beta <- NULL
if (file.exists(beta_path)) {
  msg("  Reading Beta matrix...")
  b_dt <- fread(beta_path, data.table = FALSE)
  b_probes <- b_dt[[1]]
  b_mat <- as.matrix(b_dt[, -1, drop = FALSE])
  rownames(b_mat) <- b_probes
  rm(b_dt); gc(verbose = FALSE)
  # Subset to matching probes and samples
  common_b <- intersect(rownames(mat), rownames(b_mat))
  common_s <- intersect(samples, colnames(b_mat))
  beta <- b_mat[common_b, common_s, drop = FALSE]
  storage.mode(beta) <- "double"
  # Also subset mat to common probes
  mat_nona <- mat_nona[common_b, , drop = FALSE]
  mat <- mat[common_b, , drop = FALSE]
  rm(b_mat); gc(verbose = FALSE)
  msg("  Beta aligned: ", nrow(beta), " probes x ", ncol(beta), " samples")
}

msg("  Data loading took ", round((proc.time() - t1)[3], 1), "s")

# =============================================================================
# STEP 2: Batch correction for QC visualisation ONLY
# =============================================================================
msg("STEP 2: Batch correction (visualisation only)")

condition    <- factor(meta$condition, levels = c("HC", "MS"))
batch_vector <- meta[[dataset_col]]
batch_method <- "none"
mat_corrected <- mat_nona

if (length(unique(batch_vector)) > 1) {
  # Use removeBatchEffect (much faster than ComBat for large matrices)
  mat_corrected <- tryCatch(
    removeBatchEffect(mat_nona, batch = batch_vector,
                      design = model.matrix(~ condition)),
    error = function(e) { msg("Batch correction failed, using raw"); mat_nona }
  )
  batch_method <- "removeBatchEffect"
}

# ── MDS before correction ─────────────────────────────────────────────────────
pal_cond <- c(MS = "#E05252", HC = "#5278E0")
pch_vec  <- ifelse(meta$condition == "MS", 16, 17)

png(file.path(out_dir, "MDS_before_correction.png"), width = 900, height = 800, res = 130)
tryCatch({
  plotMDS(mat_nona, top = 1000, gene.selection = "common",
          col = pal_cond[meta$condition], pch = pch_vec,
          main = "MDS – before batch correction")
  legend("topright", legend = names(pal_cond), fill = pal_cond, cex = 0.8)
}, error = function(e) NULL)
dev.off()

# ── MDS after correction ──────────────────────────────────────────────────────
png(file.path(out_dir, "MDS_after_correction.png"), width = 900, height = 800, res = 130)
tryCatch({
  plotMDS(mat_corrected, top = 1000, gene.selection = "common",
          col = pal_cond[meta$condition], pch = pch_vec,
          main = paste0("MDS – after ", batch_method, " (visualisation only)"))
  legend("topright", legend = names(pal_cond), fill = pal_cond, cex = 0.8)
}, error = function(e) NULL)
dev.off()

# Save batch-corrected matrix
write.csv(data.frame(Probe = rownames(mat_corrected), mat_corrected, check.names = FALSE),
          file.path(out_dir, "Batch_Corrected_M.csv"), row.names = FALSE)

# =============================================================================
# STEP 3: DMP — lmFit + eBayes (UN-corrected, batch as covariate)
# =============================================================================
msg("STEP 3: DMP analysis (lmFit + eBayes, batch as covariate)")

if (length(unique(batch_vector)) > 1) {
  design <- model.matrix(~ 0 + condition + factor(batch_vector))
  colnames(design)[seq_len(nlevels(condition))] <- levels(condition)
  colnames(design) <- make.names(colnames(design))
} else {
  design <- model.matrix(~ 0 + condition)
  colnames(design) <- levels(condition)
}

contrast_mat <- makeContrasts(MS_vs_HC = MS - HC, levels = design)

fit  <- lmFit(mat_nona, design)
fit  <- contrasts.fit(fit, contrast_mat)
fit  <- eBayes(fit, robust = TRUE)

if (!is.null(ann450kSub)) {
  ann_match <- ann450kSub[match(rownames(mat_nona), ann450kSub$Probe), ]
  rownames(ann_match) <- rownames(mat_nona)
  dmp_res <- topTable(fit, coef = "MS_vs_HC", number = Inf,
                      adjust.method = "BH", sort.by = "p",
                      genelist = ann_match)
  if (!"Probe" %in% colnames(dmp_res)) dmp_res$Probe <- rownames(dmp_res)
} else {
  dmp_res <- topTable(fit, coef = "MS_vs_HC", number = Inf,
                      adjust.method = "BH", sort.by = "p")
  dmp_res$Probe <- rownames(dmp_res)
}

write.csv(dmp_res, file.path(out_dir, "DMP_Results_MS_vs_HC.csv"), row.names = FALSE)
sig_dmp <- sum(dmp_res$adj.P.Val < 0.05, na.rm = TRUE)
msg("DMP: ", nrow(dmp_res), " probes tested, ", sig_dmp, " significant (FDR<0.05)")

# ── plotCpg for top 3 significant DMPs ─────────────────────────────────────
if (!is.null(beta) && sig_dmp > 0) {
  top_cpgs <- head(dmp_res$Probe[order(dmp_res$adj.P.Val)], 3)
  for (cpg in top_cpgs) {
    png(file.path(out_dir, paste0("DMP_plotCpg_", cpg, ".png")),
        width = 700, height = 600, res = 120)
    tryCatch({
      beta_cols <- intersect(colnames(beta), meta$sample_id)
      plotCpg(beta[, beta_cols, drop = FALSE],
              cpg   = cpg,
              pheno = meta$condition[match(beta_cols, meta$sample_id)],
              ylab  = "Beta value",
              main  = paste0(cpg, " (adj.P=",
                             format(dmp_res$adj.P.Val[dmp_res$Probe == cpg][1], digits = 3), ")"))
    }, error = function(e) NULL)
    dev.off()
  }
}

# =============================================================================
# STEP 4: DMR — cpg.annotate → dmrcate → extractRanges → DMR.plot
# =============================================================================
msg("STEP 4: DMR analysis (DMRcate) — streamlined")

dmr_df  <- data.frame()
sig_dmr <- 0L

tryCatch({
  annot_obj <- cpg.annotate(
    datatype      = "array",
    object        = mat_nona,
    what          = "M",
    arraytype     = arraytype,
    analysis.type = "differential",
    design        = design,
    contrasts     = TRUE,
    cont.matrix   = contrast_mat,
    coef          = "MS_vs_HC",
    fdr           = 0.05
  )
  n_sig_probes <- sum(annot_obj@ranges$is.sig, na.rm = TRUE)
  msg("cpg.annotate done: ", n_sig_probes, " significant probes")

  dmr_obj    <- dmrcate(annot_obj, lambda = 1000, C = 2, min.cpgs = 3)
  dmr_ranges <- extractRanges(dmr_obj, genome = "hg19")
  dmr_df     <- as.data.frame(dmr_ranges)

  write.csv(dmr_df, file.path(out_dir, "DMR_Results_MS_vs_HC.csv"), row.names = FALSE)
  sig_col <- if ("min_smoothed_fdr" %in% names(dmr_df)) "min_smoothed_fdr" else
             if ("FDR" %in% names(dmr_df)) "FDR" else NA_character_
  sig_dmr <- if (!is.na(sig_col)) sum(dmr_df[[sig_col]] < 0.05, na.rm = TRUE) else nrow(dmr_df)
  msg("DMR: ", nrow(dmr_df), " regions, ", sig_dmr, " with FDR<0.05")

  # ── DMR.plot for top 3 DMRs ─────────────────────────────────────────────
  if (!is.null(beta) && nrow(dmr_ranges) > 0 && !is.null(ann_pkg)) {
    msg("  Generating DMR plots for top 3 regions...")
    beta_cols <- intersect(colnames(beta), meta$sample_id)
    common_b  <- intersect(rownames(beta), rownames(ann_pkg))
    if (length(common_b) > 500 && length(beta_cols) > 0) {
      beta_sub <- beta[common_b, beta_cols, drop = FALSE]
      gr_ann   <- GenomicRanges::GRanges(
        seqnames = ann_pkg[common_b, "chr"],
        ranges   = IRanges::IRanges(start = ann_pkg[common_b, "pos"], width = 1),
        strand   = ann_pkg[common_b, "strand"]
      )
      names(gr_ann) <- common_b
      arr_str <- if (toupper(arraytype) == "EPIC") "IlluminaHumanMethylationEPIC" else "IlluminaHumanMethylation450k"
      ann_str <- if (toupper(arraytype) == "EPIC") "ilm10b4.hg19" else "ilmn12.hg19"
      grs <- tryCatch(
        minfi::GenomicRatioSet(gr_ann, Beta = beta_sub,
                               annotation = c(array = arr_str, annotation = ann_str)),
        error = function(e) NULL
      )
      if (!is.null(grs)) {
        col_groups <- pal_cond[meta$condition[match(beta_cols, meta$sample_id)]]
        for (idx in seq_len(min(3L, nrow(dmr_ranges)))) {
          png(file.path(out_dir, paste0("DMR_plot_", idx, ".png")),
              width = 1600, height = 1000, res = 140)
          tryCatch(
            DMR.plot(ranges = dmr_ranges, dmr = idx, CpGs = grs,
                     phen.col = col_groups, genome = "hg19"),
            error = function(e) msg("DMR.plot failed for DMR ", idx, ": ", conditionMessage(e))
          )
          dev.off()
        }
        rm(grs, beta_sub, gr_ann)
      }
    }
  }

}, error = function(e) {
  msg("WARNING: DMR analysis failed — ", conditionMessage(e))
  write.csv(data.frame(error = conditionMessage(e)),
            file.path(out_dir, "DMR_Results_MS_vs_HC.csv"), row.names = FALSE)
})

# =============================================================================
# STEP 5: mCSEA — logFC rank from DMP (NBIS tutorial approach)
# =============================================================================
msg("STEP 5: mCSEA gene-set enrichment")

promoter_tested  <- 0L; promoter_sig  <- 0L
gene_body_tested <- 0L; gene_body_sig <- 0L

tryCatch({
  # Build named rank vector from logFC (NBIS tutorial: myRank <- DMPs$logFC)
  logfc_col <- if ("logFC" %in% colnames(dmp_res)) "logFC" else "B"
  probe_col <- if ("Probe" %in% colnames(dmp_res)) "Probe" else
               if ("Name"  %in% colnames(dmp_res)) "Name"  else NULL
  probes_vec <- if (!is.null(probe_col)) dmp_res[[probe_col]] else rownames(dmp_res)

  myRank        <- dmp_res[[logfc_col]]
  names(myRank) <- probes_vec
  myRank        <- myRank[is.finite(myRank) & !is.na(names(myRank)) & nchar(names(myRank)) > 0]

  # Phenotype: simple Sample_Group column (NBIS tutorial approach)
  pheno_mcsea <- data.frame(
    Sample_Group = meta$condition,
    row.names    = meta$sample_id,
    stringsAsFactors = FALSE
  )

  # Use beta values (NBIS uses bVals in mCSEATest)
  if (!is.null(beta)) {
    beta_cols  <- intersect(rownames(pheno_mcsea), colnames(beta))
    bVals_keep <- intersect(names(myRank), rownames(beta))
    bVals_mcsea <- beta[bVals_keep, beta_cols, drop = FALSE]
    pheno_use   <- pheno_mcsea[beta_cols, , drop = FALSE]
    myRank      <- myRank[bVals_keep]
  } else {
    # Convert M → beta if no beta available
    bVals_all   <- 2^mat_nona / (1 + 2^mat_nona)
    bVals_keep  <- intersect(names(myRank), rownames(bVals_all))
    bVals_mcsea <- bVals_all[bVals_keep, , drop = FALSE]
    pheno_use   <- pheno_mcsea
    myRank      <- myRank[bVals_keep]
    rm(bVals_all)
  }

  platform_str <- if (toupper(arraytype) == "EPIC") "EPIC" else "450k"
  msg("mCSEA rank vector: ", length(myRank), " probes (logFC-ranked)")

  msg("  Running mCSEA with ", nCores, " cores (promoters only for speed)...")
  myResults <- mCSEATest(
    myRank,
    bVals_mcsea,
    pheno_use,
    column       = "Sample_Group",
    regionsTypes = c("promoters"),
    platform     = platform_str,
    minCpGs      = 5,
    nproc        = nCores
  )

  saveRDS(myResults, file.path(out_dir, "mCSEA_Results.rds"))

  # ── Promoter results ────────────────────────────────────────────────────
  prom_df <- as.data.frame(myResults[["promoters"]])
  if (nrow(prom_df) > 0) {
    prom_df$Gene <- rownames(prom_df)
    prom_df      <- prom_df[order(prom_df$padj), , drop = FALSE]
    write.csv(prom_df, file.path(out_dir, "Promoter_Results_mCSEA.csv"), row.names = FALSE)
    promoter_tested <- nrow(prom_df)
    promoter_sig    <- sum(prom_df$padj < 0.05, na.rm = TRUE)
    msg("mCSEA promoters: ", promoter_tested, " tested, ", promoter_sig, " sig")

    # mCSEAPlot for top 5 significant promoter genes (with gene name headers)
    top_genes <- head(prom_df$Gene[prom_df$padj < 0.05], 5)
    if (length(top_genes) > 0) {
      msg("  Generating mCSEA promoter plots for: ", paste(top_genes, collapse = ", "))
      for (gene in top_genes) {
        safe_name <- gsub("[^A-Za-z0-9_-]", "_", gene)
        png(file.path(out_dir, paste0("Promoter_mCSEA_", safe_name, ".png")),
            width = 1800, height = 1200, res = 160)
        tryCatch({
          mCSEAPlot(myResults, regionType = "promoters", dmrName = gene,
                    transcriptAnnotation = "symbol", makePDF = FALSE)
          title(main = paste0("Promoter: ", gene), line = -1, cex.main = 1.4, font.main = 2)
        }, error = function(e) {
          msg("  mCSEAPlot failed for ", gene, ": ", conditionMessage(e))
          if (dev.cur() > 1) dev.off()
        })
        if (dev.cur() > 1) dev.off()
      }
    }
  }

  # ── Gene body results ───────────────────────────────────────────────────
  if ("genes" %in% names(myResults)) {
    gene_df <- as.data.frame(myResults[["genes"]])
    if (nrow(gene_df) > 0) {
      gene_df$Gene <- rownames(gene_df)
      gene_df      <- gene_df[order(gene_df$padj), , drop = FALSE]
      write.csv(gene_df, file.path(out_dir, "Gene_Results_mCSEA.csv"), row.names = FALSE)
      gene_body_tested <- nrow(gene_df)
      gene_body_sig    <- sum(gene_df$padj < 0.05, na.rm = TRUE)
      msg("mCSEA gene bodies: ", gene_body_tested, " tested, ", gene_body_sig, " sig")
    }
  }

}, error = function(e) {
  msg("WARNING: mCSEA failed — ", conditionMessage(e))
})

# =============================================================================
# STEP 6: GO / KEGG enrichment — missMethyl::gometh
# =============================================================================
msg("STEP 6: GO / KEGG enrichment (missMethyl::gometh)")

tryCatch({
  if (!requireNamespace("missMethyl", quietly = TRUE)) stop("missMethyl not installed")

  probe_col2 <- if ("Probe" %in% colnames(dmp_res)) "Probe" else
                if ("Name"  %in% colnames(dmp_res)) "Name"  else NULL
  if (is.null(probe_col2)) stop("Probe column not found in DMP results")

  sig_cpgs <- dmp_res[[probe_col2]][dmp_res$adj.P.Val < 0.05 & !is.na(dmp_res$adj.P.Val)]
  all_cpgs <- dmp_res[[probe_col2]][!is.na(dmp_res[[probe_col2]])]
  msg("GO: ", length(sig_cpgs), " significant CpGs / ", length(all_cpgs), " total")

  if (length(sig_cpgs) >= 10) {
    gst_go <- missMethyl::gometh(sig.cpg = sig_cpgs, all.cpg = all_cpgs,
                                  collection = "GO", array.type = arraytype,
                                  plot.bias = FALSE)
    write.csv(missMethyl::topGSA(gst_go, number = 100),
              file.path(out_dir, "GO_Enrichment_BP.csv"), row.names = FALSE)

    gst_kegg <- missMethyl::gometh(sig.cpg = sig_cpgs, all.cpg = all_cpgs,
                                    collection = "KEGG", array.type = arraytype,
                                    plot.bias = FALSE)
    write.csv(missMethyl::topGSA(gst_kegg, number = 50),
              file.path(out_dir, "GO_Enrichment_KEGG.csv"), row.names = FALSE)

    msg("GO sig terms: ", sum(gst_go$FDR < 0.05, na.rm = TRUE),
        " | KEGG sig pathways: ", sum(gst_kegg$FDR < 0.05, na.rm = TRUE))
  } else {
    msg("Too few significant CpGs (", length(sig_cpgs), ") — skipping GO enrichment")
  }
}, error = function(e) {
  msg("WARNING: GO enrichment failed — ", conditionMessage(e))
})

# =============================================================================
# Summary
# =============================================================================
summary_df <- data.frame(
  metric = c(
    "samples", "ms_samples", "hc_samples", "datasets",
    "probes_tested", "significant_DMP_fdr005",
    "DMR_total_regions", "significant_DMR_fdr005",
    "promoter_genes_tested", "significant_promoter_fdr005",
    "gene_body_genes_tested", "significant_gene_body_fdr005",
    "batch_method_visualisation", "batch_method_stats",
    "mcsea_platform", "arraytype"
  ),
  value = c(
    nrow(meta), n_ms, n_hc, length(unique(batch_vector)),
    nrow(dmp_res), sig_dmp,
    nrow(dmr_df), sig_dmr,
    promoter_tested, promoter_sig,
    gene_body_tested, gene_body_sig,
    batch_method,
    "dataset as covariate in model.matrix (M-values, UN-corrected)",
    if (toupper(arraytype) == "EPIC") "EPIC" else "450k",
    arraytype
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_df, file.path(out_dir, "Summary.csv"), row.names = FALSE)

msg("══════════════════════════════════════════════")
msg("DONE: ", out_dir)
msg("  DMP:          DMP_Results_MS_vs_HC.csv       (", sig_dmp, " sig)")
msg("  DMR:          DMR_Results_MS_vs_HC.csv        (", sig_dmr, " sig)")
msg("  mCSEA prom:   Promoter_Results_mCSEA.csv     (", promoter_sig, " sig)")
msg("  mCSEA genes:  Gene_Results_mCSEA.csv          (", gene_body_sig, " sig)")
msg("══════════════════════════════════════════════")
