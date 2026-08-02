## helpers.R — shared utilities for r_notebooks/
## Loaded by every notebook with: source("helpers.R")
## Provides:
##   - PROJ_ROOT, PROT_ROOT, OUT_DIR, FIG_DIR, CACHE_DIR
##   - CROSS_OMICS / RECURRING / PAPER_TOP / ECM_FAMILY gene sets
##   - vsn_with_fallback() — robust VSN that falls back to median-center
##   - dep_volcano_gg() — EnhancedVolcano-style ggplot with cross-omics overlay
##   - moderated_t_safe() — limma::eBayes wrapper that handles edge cases
##   - bh_fdr(), gene_dedup()
##
## Paths are LOCAL (this machine) — NOT the sandbox /sessions/ path.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(data.table)
  library(limma)
})

PROJ_ROOT <- "__MS_GEO_ROOT__"
PROT_ROOT <- file.path(PROJ_ROOT, "Proteomics")
OUT_DIR   <- file.path(PROT_ROOT, "processed", "META")
FIG_DIR   <- file.path(PROT_ROOT, "figures")
NB_DIR    <- file.path(PROT_ROOT, "r_notebooks")
CACHE_DIR <- file.path(NB_DIR, "_cache")
dir.create(OUT_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

CROSS_OMICS <- c("LXN", "SH3BP4", "CHL1", "CTSZ", "RPAP2", "PCNP", "THRB")
RECURRING   <- c("STAT3", "TYK2", "CXCL13", "MBP", "CFI", "STAT1", "APLP1",
                 "CCL20", "CD5", "GZMB", "IFNG", "IL17A", "MAPK14", "A1BG",
                 "CD8A", "ITGAM", "JCHAIN", "MAPK1", "PTEN",
                 "S100A8", "S100A9", "S100A12")
PAPER_TOP   <- c("SDC1", "CHIT1", "JCHAIN", "IGHM", "MBP", "CHI3L1", "CHI3L2",
                 "NEFL", "GFAP")
ECM_FAMILY  <- c(paste0("ANXA",  1:13),
                 paste0("S100A", 1:16), "S100B", "S100P",
                 "AHNAK", "AHNAK2")

CO_COLORS <- c(LXN="#D62828", SH3BP4="#F4A261", CHL1="#2A9D8F",
               CTSZ="#7B3FA0", RPAP2="#1F4E79", PCNP="#264653",
               THRB="#E76F51")

bh_fdr <- function(p) p.adjust(p, method = "BH")

## Robust VSN normalisation — falls back to median-centering when VSN fails.
## VSN expects LINEAR-scale intensities. If the matrix looks log-transformed
## (median ≤ 30 typically), skip VSN and median-center instead.
vsn_with_fallback <- function(mat) {
  med <- median(mat, na.rm = TRUE)
  is_log <- !is.na(med) && med <= 30
  if (!is_log && requireNamespace("vsn", quietly = TRUE)) {
    res <- try(suppressMessages(vsn::justvsn(mat)), silent = TRUE)
    if (!inherits(res, "try-error")) {
      message("  Applied VSN (justvsn)")
      return(res)
    }
    message("  VSN failed — falling back to median-center")
  }
  if (is_log) message("  Matrix already looks log-transformed (median=",
                      sprintf("%.1f", med), ") — using median-center only")
  sm <- apply(mat, 2, median, na.rm = TRUE)
  mat - matrix(sm, nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE) +
    median(sm, na.rm = TRUE)
}

## DEP-equivalent filter: keep proteins with ≥ thr fraction of valid values
## (non-NA) in every condition. Mirrors DEP::filter_proteins(type="fraction").
filter_missval_R <- function(mat, group_vec, thr = 0.5) {
  keep <- rep(TRUE, nrow(mat))
  for (g in unique(group_vec)) {
    cols <- which(group_vec == g)
    valid_frac <- rowMeans(!is.na(mat[, cols, drop = FALSE]))
    keep <- keep & (valid_frac >= thr)
  }
  message(sprintf("  filter_missval (>= %.0f%% valid/cond): %d / %d kept",
                  100*thr, sum(keep), length(keep)))
  mat[keep, , drop = FALSE]
}

## DEP-equivalent MinProb imputation: draw from N(q-quantile, 0.3*sd)
## per protein. Mirrors imputeLCMD::impute.MinProb (DEP default).
impute_minprob_R <- function(mat, q = 0.01, seed = 42) {
  set.seed(seed)
  out <- mat
  for (i in seq_len(nrow(mat))) {
    row <- mat[i, ]
    miss <- is.na(row)
    if (!any(miss)) next
    obs <- row[!miss]
    if (length(obs) == 0) next
    qv <- if (length(obs) >= 5) as.numeric(quantile(obs, q)) else min(obs)
    sdv <- sd(obs, na.rm = TRUE); if (is.na(sdv)) sdv <- 0.01
    out[i, miss] <- rnorm(sum(miss), mean = qv, sd = max(sdv * 0.3, 1e-3))
  }
  out
}

## Full DEP-equivalent pipeline: filter → normalize → impute → limma.
## Returns a topTable-style data.frame.
dep_equivalent_de <- function(mat, group_vec, group_a = "MS",
                               group_b = "Control", thr = 0.5,
                               q_minprob = 0.01) {
  m  <- filter_missval_R(mat, group_vec, thr = thr)
  m  <- vsn_with_fallback(m)
  m  <- impute_minprob_R(m, q = q_minprob)
  message(sprintf("  -> running limma on %d proteins x %d samples", nrow(m), ncol(m)))
  moderated_t_safe(m, group_vec, group_a = group_a, group_b = group_b)
}

## Collapse probes/proteins to gene-level by max-variance representative.
gene_dedup <- function(mat, gene_vec) {
  stopifnot(length(gene_vec) == nrow(mat))
  v <- apply(mat, 1, var, na.rm = TRUE)
  ord <- order(-v)
  keep <- !duplicated(gene_vec[ord]) & !is.na(gene_vec[ord]) & gene_vec[ord] != ""
  out <- mat[ord[keep], , drop = FALSE]
  rownames(out) <- gene_vec[ord[keep]]
  out
}

## Safer eBayes wrapper with trend + robust defaults; handles small-n.
moderated_t_safe <- function(mat, group_vec, group_a = "MS", group_b = "Control",
                              batch = NULL) {
  group_vec <- factor(group_vec, levels = c(group_b, group_a))
  if (is.null(batch)) {
    design <- model.matrix(~ group_vec)
    colnames(design) <- c("Intercept", paste0(group_a, "_vs_", group_b))
    fit <- lmFit(mat, design)
    fit <- eBayes(fit, trend = TRUE, robust = TRUE)
    out <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  } else {
    batch <- factor(batch)
    design <- model.matrix(~ batch + group_vec)
    coef <- tail(colnames(design), 1)
    fit <- lmFit(mat, design)
    fit <- eBayes(fit, trend = TRUE, robust = TRUE)
    out <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  }
  out$gene <- rownames(out)
  out
}

## EnhancedVolcano-style ggplot — works without EnhancedVolcano dependency.
dep_volcano_gg <- function(df, title, subtitle = NULL,
                            fc_col = "logFC", p_col = "P.Value",
                            fdr_col = "adj.P.Val", gene_col = "gene",
                            fc_cut = 0.5, fdr_cut = 0.05,
                            y_cap = 50, top_n_label = 12) {
  df <- as.data.frame(df)
  stopifnot(all(c(fc_col, p_col, fdr_col, gene_col) %in% colnames(df)))
  df$fc <- df[[fc_col]]
  df$pv <- df[[p_col]]
  df$fdr <- df[[fdr_col]]
  df$gene <- df[[gene_col]]
  df$y_raw <- -log10(pmax(df$pv, 1e-300))
  df$y <- pmin(df$y_raw, y_cap)
  df$capped <- df$y_raw > y_cap
  df$cat <- "NS"
  df$cat[df$pv < 0.05] <- "P"
  df$cat[df$fdr < fdr_cut] <- "FDR_only"
  df$cat[abs(df$fc) > fc_cut & df$fdr >= fdr_cut] <- "FC_only"
  df$cat[df$fdr < fdr_cut & abs(df$fc) > fc_cut] <- "BOTH"

  df$is_co  <- df$gene %in% CROSS_OMICS
  df$is_rec <- df$gene %in% RECURRING
  df$is_ecm <- df$gene %in% ECM_FAMILY
  df$is_pap <- df$gene %in% PAPER_TOP

  pal <- c(NS = "#D9D9D9", P = "#A8C7E0", FDR_only = "#3E92CC",
           FC_only = "#C8A464", BOTH = "#1F4E79")
  bg <- df[!df$is_co & !df$is_rec & !df$is_ecm & !df$is_pap, ]
  ecm <- df[df$is_ecm, ]
  rec <- df[df$is_rec, ]
  pap <- df[df$is_pap, ]
  co  <- df[df$is_co,  ]
  lbl_top <- df[order(df$pv), ][seq_len(min(top_n_label, nrow(df))), ]
  lbl <- unique(rbind(co, pap, head(lbl_top, top_n_label)))

  p <- ggplot() +
    geom_point(data = bg,
               aes(x = fc, y = y, colour = cat),
               size = 1.2, alpha = 0.5) +
    scale_colour_manual(values = pal, name = "Category") +
    geom_point(data = ecm,
               aes(x = fc, y = y),
               shape = 22, fill = "#F4A261", colour = "black",
               size = 2.6, stroke = 0.3, alpha = 0.9) +
    geom_point(data = rec,
               aes(x = fc, y = y),
               shape = 23, fill = "#7B3FA0", colour = "black",
               size = 3.0, stroke = 0.3, alpha = 0.9) +
    geom_point(data = pap,
               aes(x = fc, y = y),
               shape = 24, fill = "#0F8B5C", colour = "black",
               size = 3.2, stroke = 0.4, alpha = 0.95) +
    geom_point(data = co,
               aes(x = fc, y = y),
               shape = 8, colour = "#D62828",
               size = 5.5, stroke = 1.4) +
    geom_text_repel(data = lbl,
                    aes(x = fc, y = y, label = gene),
                    size = 3.0, fontface = "bold",
                    max.overlaps = 30,
                    box.padding = 0.5, segment.colour = "grey50") +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed",
               colour = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dotted",
               colour = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.3) +
    labs(title = title, subtitle = subtitle,
         x = expression(log[2]~"Fold Change"),
         y = expression(-log[10]~"p-value")) +
    theme_classic(base_size = 11) +
    theme(plot.title.position = "plot",
          legend.position = "right",
          panel.grid.major = element_line(colour = "grey95"))
  p
}

cat(sprintf("[helpers.R] loaded · PROJ=%s\n", PROJ_ROOT))
cat(sprintf("[helpers.R] %d cross-omics, %d recurring, %d paper-top, %d ECM\n",
            length(CROSS_OMICS), length(RECURRING), length(PAPER_TOP),
            length(ECM_FAMILY)))
