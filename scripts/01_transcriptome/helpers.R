## helpers.R — shared utilities for Transcriptome/r_notebooks
## Loaded with: source("helpers.R")
##
## Provides:
##   - PROJ_ROOT, TX_ROOT, OUT_DIR, FIG_DIR, STRATA_DIR
##   - CROSS_OMICS / RECURRING / PAPER_TOP / ECM_FAMILY (same panels as proteomics)
##   - load_stratum(name) -> list(mat, groups, meta) for a stratum
##   - run_limma_stratum(mat, groups) -> topTable data.frame
##   - tx_volcano_gg(df, ...) -> ggplot
##   - bh_fdr(), gene_dedup()

suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(dplyr); library(data.table); library(limma)
})

PROJ_ROOT   <- "__MS_GEO_ROOT__"
TX_ROOT     <- file.path(PROJ_ROOT, "Transcriptome")
STRATA_DIR  <- file.path(PROJ_ROOT, "Stratified_Analyses", "Expression")
OUT_DIR     <- file.path(TX_ROOT, "results")
FIG_DIR     <- file.path(TX_ROOT, "figures")
NB_DIR      <- file.path(TX_ROOT, "r_notebooks")
CACHE_DIR   <- file.path(NB_DIR, "_cache")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

CROSS_OMICS <- c("LXN", "SH3BP4", "CHL1", "CTSZ", "RPAP2", "PCNP", "THRB")
RECURRING   <- c("STAT3", "TYK2", "CXCL13", "MBP", "CFI", "STAT1", "APLP1",
                 "CCL20", "CD5", "GZMB", "IFNG", "IL17A", "MAPK14", "A1BG",
                 "CD8A", "ITGAM", "JCHAIN")
PAPER_TOP   <- c("SDC1", "CHIT1", "JCHAIN", "IGHM", "MBP", "CHI3L1", "CHI3L2",
                 "NEFL", "GFAP")
ECM_FAMILY  <- c(paste0("ANXA", 1:13),
                 paste0("S100A", 1:16), "S100B", "S100P",
                 "AHNAK", "AHNAK2")

bh_fdr <- function(p) p.adjust(p, method = "BH")

## Collapse duplicate gene rows to max-variance representative
gene_dedup <- function(mat, gene_vec) {
  v <- apply(mat, 1, var, na.rm = TRUE)
  ord <- order(-v, na.last = TRUE)
  m   <- mat[ord, , drop = FALSE]
  g   <- gene_vec[ord]
  keep <- !duplicated(g) & !is.na(g) & g != ""
  out  <- m[keep, , drop = FALSE]
  rownames(out) <- g[keep]
  out
}

## Load a stratum directory: returns list(mat, groups, meta).
## Accepts both "cell_tissue_case_control_pbmc" and "pbmc" shorthand.
load_stratum <- function(name) {
  full_path <- name
  if (!dir.exists(name)) {
    full_path <- file.path(STRATA_DIR, name)
    if (!dir.exists(full_path)) {
      candidates <- list.dirs(STRATA_DIR, recursive = FALSE)
      hit <- candidates[grepl(paste0("_", name, "$"), candidates) |
                         grepl(name, basename(candidates))]
      if (length(hit) == 0)
        stop(sprintf("No stratum matching '%s' in %s", name, STRATA_DIR))
      full_path <- hit[1]
    }
  }
  message(sprintf("[load_stratum] reading %s", basename(full_path)))
  expr <- fread(file.path(full_path, "Batch_Corrected_Expression.csv"),
                showProgress = FALSE)
  meta <- fread(file.path(full_path, "metadata.csv"), showProgress = FALSE)

  gene_col <- "Gene"
  if (!"Gene" %in% colnames(expr)) gene_col <- colnames(expr)[1]
  genes <- expr[[gene_col]]
  sample_cols <- setdiff(colnames(expr), gene_col)
  mat <- as.matrix(expr[, ..sample_cols])
  storage.mode(mat) <- "numeric"
  mat <- gene_dedup(mat, genes)

  keep <- intersect(colnames(mat), meta$sample_id)
  mat <- mat[, keep, drop = FALSE]
  groups <- meta$condition[match(keep, meta$sample_id)]
  message(sprintf("  %d genes × %d samples  (groups: %s)",
                  nrow(mat), ncol(mat),
                  paste(sprintf("%s=%d", names(table(groups)), table(groups)),
                        collapse = "  ")))
  list(mat = mat, groups = groups, meta = meta, stratum = basename(full_path),
       path = full_path)
}

## limma DE with eBayes(trend, robust)
run_limma_stratum <- function(mat, groups, group_a = "MS", group_b = "HC") {
  group_vec <- factor(groups, levels = c(group_b, group_a))
  if (any(is.na(group_vec))) {
    keep <- !is.na(group_vec)
    mat <- mat[, keep, drop = FALSE]
    group_vec <- droplevels(group_vec[keep])
  }
  design <- model.matrix(~ group_vec)
  colnames(design) <- c("Intercept", paste0(group_a, "_vs_", group_b))
  fit <- lmFit(mat, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  out <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  out$gene <- rownames(out)
  out$is_cross_omics <- out$gene %in% CROSS_OMICS
  out$is_recurring   <- out$gene %in% RECURRING
  out
}

## Pretty volcano with cross-omics overlay
tx_volcano_gg <- function(df, title, subtitle = NULL,
                          fc_col = "logFC", p_col = "P.Value",
                          fdr_col = "adj.P.Val", gene_col = "gene",
                          fc_cut = 0.5, fdr_cut = 0.05,
                          y_cap = 50, top_n_label = 12) {
  df <- as.data.frame(df)
  df$fc  <- df[[fc_col]]; df$pv <- df[[p_col]]
  df$fdr <- df[[fdr_col]]; df$gene <- df[[gene_col]]
  df$y_raw <- -log10(pmax(df$pv, 1e-300))
  df$y <- pmin(df$y_raw, y_cap)
  df$cat <- "NS"
  df$cat[df$pv < 0.05]                                 <- "P"
  df$cat[df$fdr < fdr_cut]                             <- "FDR_only"
  df$cat[abs(df$fc) > fc_cut & df$fdr >= fdr_cut]      <- "FC_only"
  df$cat[df$fdr < fdr_cut & abs(df$fc) > fc_cut]       <- "BOTH"
  df$is_co  <- df$gene %in% CROSS_OMICS
  df$is_rec <- df$gene %in% RECURRING

  pal <- c(NS = "#D9D9D9", P = "#A8C7E0", FDR_only = "#3E92CC",
           FC_only = "#C8A464", BOTH = "#1F4E79")
  bg <- df[!df$is_co & !df$is_rec, ]
  rec <- df[df$is_rec, ]
  co  <- df[df$is_co,  ]
  lbl_top <- df[order(df$pv), ][seq_len(min(top_n_label, nrow(df))), ]
  lbl <- unique(rbind(co, head(lbl_top, top_n_label)))

  ggplot() +
    geom_point(data = bg, aes(fc, y, colour = cat), size = 1.2, alpha = 0.5) +
    scale_colour_manual(values = pal, name = "Category") +
    geom_point(data = rec, aes(fc, y), shape = 23, fill = "#7B3FA0",
               colour = "black", size = 2.8, stroke = 0.3, alpha = 0.9) +
    geom_point(data = co, aes(fc, y), shape = 8, colour = "#D62828",
               size = 5.5, stroke = 1.4) +
    geom_text_repel(data = lbl, aes(fc, y, label = gene),
                    size = 3.0, fontface = "bold", max.overlaps = 30,
                    box.padding = 0.5, segment.colour = "grey50") +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed",
               colour = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dotted",
               colour = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.3) +
    labs(title = title, subtitle = subtitle,
         x = expression(log[2]~"Fold Change  (MS / HC)"),
         y = expression(-log[10]~"p-value (limma eBayes)")) +
    theme_classic(base_size = 11) +
    theme(plot.title.position = "plot", legend.position = "right",
          panel.grid.major = element_line(colour = "grey95"))
}

cat(sprintf("[helpers.R] Transcriptome  ·  STRATA_DIR=%s\n", STRATA_DIR))
cat(sprintf("[helpers.R] %d cross-omics, %d recurring, %d paper-top\n",
            length(CROSS_OMICS), length(RECURRING), length(PAPER_TOP)))
