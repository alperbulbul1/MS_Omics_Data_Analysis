## helpers.R — shared utilities for Methylation/r_notebooks
## Loaded with: source("helpers.R")
##
## Provides:
##   - PROJ_ROOT, METH_ROOT, OUT_DIR, FIG_DIR, STRATA_DIR, COMBINED_DIR
##   - CROSS_OMICS / METH_TOP / cross-omics-by-tissue panels
##   - load_meth_stratum(name) -> list(mat, groups, meta) for a meth stratum
##   - load_combined_meth(filter_tissue=NULL) -> combined M matrix + meta
##   - run_limma_meth(mat, groups, batch=NULL) -> topTable
##   - meth_volcano_gg(df, ...) -> ggplot
##   - annotate_probes_to_genes(probes) -> data.frame using missMethyl/minfi
##   - bh_fdr()

suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(dplyr); library(data.table); library(limma)
})

PROJ_ROOT     <- "__MS_GEO_ROOT__"
METH_ROOT     <- file.path(PROJ_ROOT, "Methylation")
STRATA_DIR    <- file.path(PROJ_ROOT, "Stratified_Analyses", "Methylation")
COMBINED_DIR  <- file.path(PROJ_ROOT, "Methylation_Data")
OUT_DIR       <- file.path(METH_ROOT, "results")
FIG_DIR       <- file.path(METH_ROOT, "figures")
NB_DIR        <- file.path(METH_ROOT, "r_notebooks")
CACHE_DIR     <- file.path(NB_DIR, "_cache")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

## Cross-omics 7-gene panel (target of validation)
CROSS_OMICS <- c("LXN", "SH3BP4", "CHL1", "CTSZ", "RPAP2", "PCNP", "THRB")
## Methylation top hits previously discovered in the bulk pipeline
METH_TOP <- c("ZIC4", "HOXA3", "DUSP22", "RAI1", "TMEM140",
              "ZNF471", "NKX6-2", "THRB", "SALL1", "CNTNAP2",
              "DIABLO", "LTA", "SPI1", "DGKZ", "C20orf123", "PRF1",
              "MIR886", "MIR596", "PM20D1", "HKR1", "VARS2",
              "CRISP2", "SLFN12", "CTHRC1", "POFUT2", "CTNNB1",
              "KRTAP12-3", "MYCBPAP", "HOXA2")

bh_fdr <- function(p) p.adjust(p, method = "BH")

## Load a stratum directory: returns list(mat, groups, meta).
load_meth_stratum <- function(name) {
  full_path <- name
  if (!dir.exists(name)) {
    full_path <- file.path(STRATA_DIR, name)
    if (!dir.exists(full_path)) {
      candidates <- list.dirs(STRATA_DIR, recursive = FALSE)
      hit <- candidates[grepl(paste0("_", name, "$"), candidates) |
                         grepl(name, basename(candidates), ignore.case = TRUE)]
      if (length(hit) == 0)
        stop(sprintf("No meth stratum matching '%s' in %s", name, STRATA_DIR))
      full_path <- hit[1]
    }
  }
  message(sprintf("[load_meth_stratum] %s", basename(full_path)))
  m_fp <- file.path(full_path, "Batch_Corrected_M.csv")
  expr <- fread(m_fp, showProgress = FALSE)
  meta <- fread(file.path(full_path, "metadata.csv"), showProgress = FALSE)

  probe_col <- "Probe"
  if (!"Probe" %in% colnames(expr)) probe_col <- colnames(expr)[1]
  probes <- expr[[probe_col]]
  sample_cols <- setdiff(colnames(expr), probe_col)
  mat <- as.matrix(expr[, ..sample_cols])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- probes

  keep <- intersect(colnames(mat), meta$sample_id)
  mat <- mat[, keep, drop = FALSE]
  groups <- meta$condition[match(keep, meta$sample_id)]
  message(sprintf("  %d probes × %d samples  (groups: %s)",
                  nrow(mat), ncol(mat),
                  paste(sprintf("%s=%d", names(table(groups)), table(groups)),
                        collapse = "  ")))
  list(mat = mat, groups = groups, meta = meta,
       stratum = basename(full_path), path = full_path)
}

## Load combined methylation matrix (317 MB — slow). Optionally filter by
## cell_type tissue ("whole_blood", "brain_wm", "T cells", etc.).
load_combined_meth <- function(filter_tissue = NULL, max_probes = NULL) {
  meta_fp <- file.path(COMBINED_DIR, "Combined_Methylation_Metadata.csv")
  mat_fp  <- file.path(COMBINED_DIR, "Combined_Methylation_Batch_Corrected.csv")
  message(sprintf("[load_combined_meth] meta: %s", basename(meta_fp)))
  meta <- fread(meta_fp)
  if (!is.null(filter_tissue)) {
    keep_samp <- meta$sample_id[meta$cell_type %in% filter_tissue]
    message(sprintf("  filtering for tissue=%s -> %d samples",
                    paste(filter_tissue, collapse = ","), length(keep_samp)))
  } else {
    keep_samp <- meta$sample_id
  }
  message(sprintf("[load_combined_meth] M-matrix: %s (317 MB)", basename(mat_fp)))
  # Read only the sample columns we need to save RAM
  hdr <- fread(mat_fp, nrows = 0); cn <- colnames(hdr)
  probe_col <- cn[1]
  use_cols <- c(probe_col, intersect(keep_samp, cn))
  expr <- fread(mat_fp, select = use_cols, showProgress = TRUE,
                nrows = if (is.null(max_probes)) Inf else max_probes)
  probes <- expr[[probe_col]]
  sample_cols <- setdiff(colnames(expr), probe_col)
  mat <- as.matrix(expr[, ..sample_cols])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- probes
  meta_sub <- meta[match(sample_cols, meta$sample_id)]
  message(sprintf("  loaded: %d probes × %d samples", nrow(mat), ncol(mat)))
  list(mat = mat, meta = meta_sub, groups = meta_sub$condition)
}

## limma on M-values (handles optional study covariate)
run_limma_meth <- function(mat, groups, batch = NULL,
                            group_a = "MS", group_b = "HC") {
  group_vec <- factor(groups, levels = c(group_b, group_a))
  if (any(is.na(group_vec))) {
    keep <- !is.na(group_vec)
    mat <- mat[, keep, drop = FALSE]
    group_vec <- droplevels(group_vec[keep])
    if (!is.null(batch)) batch <- batch[keep]
  }
  if (is.null(batch)) {
    design <- model.matrix(~ group_vec)
    coef <- 2
  } else {
    batch <- factor(batch)
    design <- model.matrix(~ batch + group_vec)
    coef  <- tail(colnames(design), 1)
  }
  fit <- lmFit(mat, design)
  fit <- eBayes(fit, robust = TRUE)
  out <- topTable(fit, coef = coef, number = Inf, sort.by = "P")
  out$Probe <- rownames(out)
  out
}

## Annotate probes to nearest gene using IlluminaHumanMethylation450kanno.ilmn12.hg19
.annot_cache <- NULL
annotate_probes_to_genes <- function(probes) {
  pkg <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
  if (!requireNamespace(pkg, quietly = TRUE)) {
    warning("450k annotation package missing — install via BiocManager")
    return(data.frame(Probe = probes, gene = NA_character_, gene_group = NA_character_))
  }
  suppressPackageStartupMessages({
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19); library(minfi)
  })
  if (is.null(.annot_cache)) {
    .annot_cache <<- minfi::getAnnotation(get(pkg))
  }
  ann <- .annot_cache
  hit <- match(probes, ann$Name)
  out <- data.frame(
    Probe = probes,
    gene = sub(";.*$", "", ann$UCSC_RefGene_Name[hit]),
    gene_group = ann$UCSC_RefGene_Group[hit],
    stringsAsFactors = FALSE)
  out$gene[is.na(out$gene) | out$gene == ""] <- NA_character_
  out
}

## Aggregate probe-level DMP to gene-level by Stouffer combined p
## Aggregate probe-level limma results to gene level.
##
## The combination is an UNWEIGHTED signed Stouffer, z_i = sign(logFC_i) * Phi^-1(1 - p_i/2),
## Z = sum(z_i)/sqrt(n). Because the weighted-Z statistic sum(w z)/sqrt(sum w^2) is invariant
## under w -> lambda*w, any constant weight vector gives this same Z; in particular, since every
## probe of a gene is measured on the identical sample set, sample-size weighting w = sqrt(N) is
## constant and Z below IS the sample-size-weighted Stouffer.
##
## Effect size. For a linear summary theta = sum(a_i b_i)/sum(a_i) of probe effects b_i with
## standard errors SE_i, the corresponding Wald statistic is a weighted-Z combination with
## w_i = a_i * SE_i. The Z above uses w_i = 1, so the effect-size estimator consistent with it is
## a_i = 1/SE_i, i.e. the 1/SE-weighted mean `se_logFC`. That is the effect size to report
## alongside these p-values. `mean_logFC` (the unweighted arithmetic mean) is retained for
## backward comparability but is NOT the estimator matching the test, and `ivw_logFC`
## (inverse-variance, a_i = 1/SE_i^2) belongs to a different test (w_i = 1/SE) and is provided
## only for sensitivity analysis. min/max_logFC give the per-gene range.
##
## Caveat, deliberately recorded here: the sqrt(n) denominator assumes the probes are independent.
## CpGs of the same gene are co-methylated (median residual pairwise rho ~0.06-0.10 in these
## strata), so gene-level p-values are anti-conservative. See
## Methylation/r_notebooks/15_genelevel_weighting_corrected.R for the dependence-aware variant
## Z = sum(z)/sqrt(n + 2*sum_{i<j} rho_ij).
probe_to_gene_stouffer <- function(dmp_df, probe_to_gene_df) {
  d <- merge(dmp_df, probe_to_gene_df, by = "Probe")
  d <- d[!is.na(d$gene) & d$gene != "", ]
  d$SE <- abs(ifelse(d$t == 0, NA_real_, d$logFC / d$t))
  res <- d %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      n_probes = dplyr::n(),
      mean_logFC = mean(logFC, na.rm = TRUE),
      se_logFC   = sum(logFC / SE, na.rm = TRUE) / sum(1 / SE, na.rm = TRUE),
      ivw_logFC  = sum(logFC / SE^2, na.rm = TRUE) / sum(1 / SE^2, na.rm = TRUE),
      min_logFC  = min(logFC, na.rm = TRUE),
      max_logFC  = max(logFC, na.rm = TRUE),
      z_combined = sum(qnorm(P.Value / 2, lower.tail = FALSE) * sign(logFC),
                       na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop") %>%
    dplyr::mutate(P.Value = 2 * pnorm(abs(z_combined), lower.tail = FALSE),
                  adj.P.Val = p.adjust(P.Value, method = "BH")) %>%
    dplyr::arrange(P.Value)
  as.data.frame(res)
}

## Volcano for methylation (M-value DE, x-axis is delta-M)
meth_volcano_gg <- function(df, title, subtitle = NULL,
                            fc_col = "logFC", p_col = "P.Value",
                            fdr_col = "adj.P.Val",
                            label_col = "gene",  # or "Probe"
                            fc_cut = 0.3, fdr_cut = 0.05,
                            y_cap = 50, top_n_label = 12) {
  df <- as.data.frame(df)
  df$fc  <- df[[fc_col]]; df$pv <- df[[p_col]]
  df$fdr <- df[[fdr_col]]; df$lbl <- df[[label_col]]
  df$y_raw <- -log10(pmax(df$pv, 1e-300))
  df$y <- pmin(df$y_raw, y_cap)
  df$cat <- "NS"
  df$cat[df$pv < 0.05]                              <- "P"
  df$cat[df$fdr < fdr_cut]                          <- "FDR_only"
  df$cat[abs(df$fc) > fc_cut & df$fdr >= fdr_cut]   <- "FC_only"
  df$cat[df$fdr < fdr_cut & abs(df$fc) > fc_cut]    <- "BOTH"
  df$is_co  <- df$lbl %in% CROSS_OMICS
  df$is_top <- df$lbl %in% METH_TOP

  pal <- c(NS = "#D9D9D9", P = "#A8C7E0", FDR_only = "#3E92CC",
           FC_only = "#C8A464", BOTH = "#1F4E79")
  bg <- df[!df$is_co & !df$is_top, ]
  topm <- df[df$is_top, ]
  co  <- df[df$is_co,  ]
  lbl_top <- df[order(df$pv), ][seq_len(min(top_n_label, nrow(df))), ]
  lbl <- unique(rbind(co, head(lbl_top, top_n_label)))

  ggplot() +
    geom_point(data = bg, aes(fc, y, colour = cat), size = 1.2, alpha = 0.5) +
    scale_colour_manual(values = pal, name = "Category") +
    geom_point(data = topm, aes(fc, y), shape = 23, fill = "#7B3FA0",
               colour = "black", size = 2.8, stroke = 0.3, alpha = 0.9) +
    geom_point(data = co, aes(fc, y), shape = 8, colour = "#D62828",
               size = 5.5, stroke = 1.4) +
    geom_text_repel(data = lbl, aes(fc, y, label = lbl),
                    size = 3.0, fontface = "bold", max.overlaps = 30,
                    box.padding = 0.5, segment.colour = "grey50") +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed",
               colour = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dotted",
               colour = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.3) +
    labs(title = title, subtitle = subtitle,
         x = expression(Delta~"M-value (MS - HC)"),
         y = expression(-log[10]~"p-value (limma eBayes)")) +
    theme_classic(base_size = 11) +
    theme(plot.title.position = "plot", legend.position = "right",
          panel.grid.major = element_line(colour = "grey95"))
}

cat(sprintf("[helpers.R] Methylation  ·  STRATA_DIR=%s\n", STRATA_DIR))
cat(sprintf("[helpers.R] %d cross-omics, %d meth-top\n",
            length(CROSS_OMICS), length(METH_TOP)))
