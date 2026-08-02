## 15_genelevel_weighting_corrected.R
##
## Reviewer 1, point 8 — corrected replacement for 14_genelevel_weighting_final.R.
##
## Two errors in script 14 are fixed here.
##
## ERROR 1 — the wrong effect-size comparator.
##   For any linear summary theta_hat = sum(a_i b_i)/sum(a_i) of probe effects b_i with standard
##   errors SE_i, the Wald statistic is
##       Z = sum(a_i b_i)/sqrt(sum(a_i^2 SE_i^2)) = sum((a_i SE_i) z_i)/sqrt(sum((a_i SE_i)^2)),
##   i.e. a Liptak combination with weights w_i = a_i * SE_i. The consistency rule is therefore
##   w_i = a_i * SE_i, not w_i = a_i.
##   The PUBLISHED test is the unweighted Stouffer, w_i = 1, which corresponds to a_i = 1/SE_i.
##   The effect size consistent with the p-values the manuscript already reports is thus the
##   1/SE-weighted mean
##       theta_SE = sum(logFC_i/SE_i) / sum(1/SE_i),
##   NOT the inverse-variance (1/SE^2) weighted mean that script 14 proposed. Reporting theta_SE
##   changes no p-value, no FDR and no candidate call: it is a pure effect-size correction, which
##   is exactly what the reviewer asked for.
##   (The IVW mean is retained below only as a secondary column, because it belongs to a DIFFERENT
##   test — the Liptak w=1/SE combination — and adopting it would silently change the significance
##   calculation as well.)
##
## ERROR 2 — the wrong correlation estimand.
##   The Stouffer denominator needs Corr(z_i, z_j). The z_i are functions of the probe test
##   statistics from the model ~ batch + group, so the required correlation is that of the probe
##   RESIDUALS under that design. Script 14 used the raw Pearson correlation of M-values, which
##   additionally contains the group and study effects the model removes, and is therefore
##   upward-biased precisely for the probes with the strongest case-control signal. That bias is
##   visible in script 14's own output: genes significant under the published scheme had median
##   rho-hat 0.25-0.53 against 0.04-0.10 elsewhere. Here rho is estimated from residuals of the
##   same design limma fitted.
##
## Schemes reported:
##   A  published        Z = sum(z)/sqrt(n)                      effect = mean(logFC)
##   B  consistent       Z = sum(z)/sqrt(n)   [UNCHANGED]        effect = sum(b/SE)/sum(1/SE)
##   C  dependence-aware Z = sum(z)/sqrt(n + 2*sum rho_res)      effect = sum(b/SE)/sum(1/SE)
##   plus per-gene min/max probe logFC (the range the reviewer suggested).
##
## Output: Methylation/results/15_genelevel_weighting_corrected.tsv

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(limma)
})
setwd("__MS_GEO_ROOT__/Methylation/r_notebooks")
source("helpers.R")

RES <- "__MS_GEO_ROOT__/Methylation/results"
PANEL <- c("ITGB2", "CD79B", "IKZF1", "LXN", "SH3BP4",
           "CASP6", "CASP8", "DGKQ", "MX1", "IFIT1", "NUP210", "RUNX3")
MAXPR <- 40

STRATA_DIRS <- list(
  "01_tcells"           = "cell_tissue_case_control_t_cells",
  "02_wb_dmf"           = "label_context_case_control_whole_blood_dmf",
  "03_wb_ocrelizumab"   = "label_context_case_control_whole_blood_ocrelizumab",
  "04_tcells_remission" = "label_context_case_control_t_cells_remission")
STUBS <- list("01_tcells" = "01_tcells_meth", "02_wb_dmf" = "02_wb_dmf_meth",
              "03_wb_ocrelizumab" = "03_wb_ocrelizumab_meth",
              "04_tcells_remission" = "04_tcells_remission_meth",
              "05_combined" = "05_combined_meth")

## Return the residual matrix under the SAME design the DMP fit used, so that the estimated
## correlation is the correlation of the statistics being combined.
residual_mat <- function(nm) {
  if (nm == "05_combined") {
    meta <- fread(file.path(COMBINED_DIR, "Combined_Methylation_Metadata.csv"))
    keep <- meta$sample_id[meta$condition %in% c("MS", "HC")]
    fp <- file.path(COMBINED_DIR, "Combined_Methylation_Batch_Corrected.csv")
    cn <- colnames(fread(fp, nrows = 0))
    expr <- fread(fp, select = c(cn[1], intersect(keep, cn)), showProgress = FALSE)
    pr <- expr[[1]]; sc <- setdiff(colnames(expr), cn[1])
    M <- as.matrix(expr[, ..sc]); storage.mode(M) <- "numeric"; rownames(M) <- pr
    rm(expr); invisible(gc())
    ms <- meta[match(sc, meta$sample_id)]
    grp <- factor(ms$condition, levels = c("HC", "MS")); bat <- factor(ms$dataset)
  } else {
    s <- load_meth_stratum(STRATA_DIRS[[nm]])
    M <- s$mat
    grp <- factor(s$groups, levels = c("HC", "MS"))
    bat <- NULL
    if ("dataset" %in% colnames(s$meta)) {
      b <- s$meta$dataset[match(colnames(M), s$meta$sample_id)]
      if (length(unique(b)) > 1) bat <- factor(b)
    }
  }
  ok <- !is.na(grp)
  M <- M[, ok, drop = FALSE]; grp <- droplevels(grp[ok])
  design <- if (is.null(bat)) model.matrix(~ grp) else {
    bat <- droplevels(bat[ok]); model.matrix(~ bat + grp)
  }
  ## residuals of the fitted linear model, probe by probe
  H <- design %*% solve(crossprod(design)) %*% t(design)
  R <- M - M %*% H
  rownames(R) <- rownames(M)
  R
}

all_out <- list()

for (nm in names(STUBS)) {
  dmp <- fread(file.path(RES, paste0(STUBS[[nm]], "_DMP.tsv")), sep = "\t", header = TRUE)
  p2g <- annotate_probes_to_genes(dmp$Probe)
  d <- merge(as.data.frame(dmp), p2g, by = "Probe")
  d <- d[!is.na(d$gene) & d$gene != "", ]
  d$SE <- abs(ifelse(d$t == 0, NA_real_, d$logFC / d$t))
  d <- d[is.finite(d$SE) & d$SE > 0, ]
  d$z <- qnorm(d$P.Value / 2, lower.tail = FALSE) * sign(d$logFC)

  R <- residual_mat(nm)
  d <- d[d$Probe %in% rownames(R), ]

  spl <- split(d$Probe, d$gene)
  rho_sum <- vapply(spl, function(pr) {
    n_full <- length(pr)
    if (n_full < 2) return(0)
    sub <- if (n_full > MAXPR) pr[seq_len(MAXPR)] else pr
    cm <- suppressWarnings(cor(t(R[sub, , drop = FALSE])))
    o <- cm[upper.tri(cm)]; o <- o[is.finite(o)]
    if (!length(o)) return(0)
    mean(o) * choose(n_full, 2)
  }, numeric(1))
  rm(R); invisible(gc())

  g <- d %>%
    group_by(gene) %>%
    summarise(n_probes  = n(),
              sum_z     = sum(z),
              mean_logFC = mean(logFC),                                  # A: as published
              se_logFC   = sum(logFC / SE) / sum(1 / SE),                # B: matches the published Z
              ivw_logFC  = sum(logFC / SE^2) / sum(1 / SE^2),            # secondary, different test
              min_logFC  = min(logFC),
              max_logFC  = max(logFC),
              .groups = "drop") %>%
    mutate(rho_sum = as.numeric(rho_sum[gene]),
           infl = sqrt(pmax(n_probes + 2 * rho_sum, 1e-6)) / sqrt(n_probes),
           Z_A  = sum_z / sqrt(n_probes),
           Z_C  = sum_z / sqrt(pmax(n_probes + 2 * rho_sum, 1e-6)),
           P_A  = 2 * pnorm(abs(Z_A), lower.tail = FALSE),
           P_C  = 2 * pnorm(abs(Z_C), lower.tail = FALSE),
           FDR_A = p.adjust(P_A, method = "BH"),
           FDR_C = p.adjust(P_C, method = "BH"),
           stratum = nm)
  all_out[[nm]] <- g

  m <- g[g$n_probes > 1, ]
  cat(sprintf("%-20s %6d genes (%5d multi-probe) | median residual rho %+.4f | median inflation %.3f | FDR<0.05: published %d, dependence-aware %d\n",
              nm, nrow(g), nrow(m),
              median(m$rho_sum / choose(m$n_probes, 2), na.rm = TRUE),
              median(m$infl, na.rm = TRUE), sum(g$FDR_A < 0.05), sum(g$FDR_C < 0.05)))
}

out <- bind_rows(all_out)
write.table(out, file.path(RES, "15_genelevel_weighting_corrected.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n", strrep("=", 108), "\n", sep = "")
cat("EFFECT SIZE: does the consistent (1/SE-weighted) estimator change anything?\n")
cat(strrep("=", 108), "\n", sep = "")
for (nm in names(all_out)) {
  m <- all_out[[nm]]; m <- m[m$n_probes > 1, ]
  cat(sprintf("  %-20s r(mean, SE-wtd) = %.4f | sign disagreement %.2f%% | median |diff| %.4f   [IVW: %.2f%% sign disagreement]\n",
              nm, cor(m$mean_logFC, m$se_logFC),
              100 * mean(sign(m$mean_logFC) != sign(m$se_logFC)),
              median(abs(m$mean_logFC - m$se_logFC)),
              100 * mean(sign(m$mean_logFC) != sign(m$ivw_logFC))))
}

cat("\n", strrep("=", 108), "\n", sep = "")
cat("CANDIDATES, pan-tissue combined stratum\n")
cat(strrep("=", 108), "\n", sep = "")
cat(sprintf("  %-8s %4s %9s %9s %9s %16s %5s %11s %11s\n",
            "gene", "nprb", "mean", "SE-wtd", "IVW", "probe range", "infl", "FDR_pub", "FDR_dep"))
for (gn in PANEL) {
  r <- all_out[["05_combined"]][all_out[["05_combined"]]$gene == gn, ]
  if (!nrow(r)) next
  st <- function(x) if (x < 0.05) "*" else " "
  cat(sprintf("  %-8s %4d %+9.4f %+9.4f %+9.4f %7.3f..%+.3f %5.2f %10.2e%s %10.2e%s\n",
              gn, r$n_probes, r$mean_logFC, r$se_logFC, r$ivw_logFC,
              r$min_logFC, r$max_logFC, r$infl,
              r$FDR_A, st(r$FDR_A), r$FDR_C, st(r$FDR_C)))
}

cat("\n-- published candidate calls lost under the dependence-aware test --\n")
lost <- out %>% filter(gene %in% PANEL, FDR_A < 0.05, FDR_C >= 0.05) %>%
  transmute(gene, stratum, FDR_published = signif(FDR_A, 2),
            FDR_dependence_aware = signif(FDR_C, 2), inflation = round(infl, 2))
if (nrow(lost)) print(as.data.frame(lost), row.names = FALSE) else cat("  none\n")

cat("\n-- raw vs residual rho, to document the size of the correction to script 14 --\n")
old <- tryCatch(fread(file.path(RES, "14_genelevel_weighting_final.tsv"), sep = "\t"),
                error = function(e) NULL)
if (!is.null(old)) {
  cmp <- merge(as.data.frame(old)[, c("gene", "stratum", "n_probes", "rho_sum", "infl")],
               as.data.frame(out)[, c("gene", "stratum", "rho_sum", "infl")],
               by = c("gene", "stratum"), suffixes = c("_raw", "_resid"))
  cmp <- cmp[cmp$n_probes > 1, ]
  for (nm in names(all_out)) {
    s <- cmp[cmp$stratum == nm, ]
    if (!nrow(s)) next
    cat(sprintf("  %-20s median mean-rho raw %+.4f -> residual %+.4f | median inflation %.3f -> %.3f\n",
                nm,
                median(s$rho_sum_raw / choose(s$n_probes, 2)),
                median(s$rho_sum_resid / choose(s$n_probes, 2)),
                median(s$infl_raw), median(s$infl_resid)))
  }
}
cat("\nwrote", file.path(RES, "15_genelevel_weighting_corrected.tsv"), "\n")
