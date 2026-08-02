#!/usr/bin/env Rscript
## 09_cross_assay_lxn.R  —  generated from notebook spec
## Run: Rscript 09_cross_assay_lxn.R


## ============================================================
## # 09 — Cross-assay summary: 7 cross-omics genes × N R-rerun assays
## 
## Pulls per-gene log2FC + FDR for `LXN, SH3BP4, CHL1, CTSZ, RPAP2, PCNP, THRB`
## across every R-side analysis (notebooks 01, 02, 03, 04, 05, 06, 07).
## 
## **Outputs**
## - `figures/CrossAssay_R_summary.tsv`
## - `figures/CrossAssay_R_grid.png`
## ============================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(dplyr)
  library(gridExtra)
})
source("helpers.R")


assays <- list(
  list(name = "CSF Astral (R/DEP)",
       fp = file.path(OUT_DIR, "CSF_Astral_R_DEP_results.tsv")),
  list(name = "CSF timsTOF (R/DEP)",
       fp = file.path(OUT_DIR, "CSF_timsTOF_R_DEP_results.tsv")),
  list(name = "CSF combined (R/ComBat)",
       fp = file.path(OUT_DIR, "CSF_combined_R_ComBat_DE.tsv")),
  list(name = "Brain CTX (Magliozzi)",
       fp = file.path(OUT_DIR, "Magliozzi_R_DEP_MS_CTX_vs_ODC_CTX.tsv")),
  list(name = "Brain NAWM (Magliozzi)",
       fp = file.path(OUT_DIR, "Magliozzi_R_DEP_MS_NAWM_vs_ODC_WM.tsv")),
  list(name = "Brain WML vs ctrl (Magliozzi)",
       fp = file.path(OUT_DIR, "Magliozzi_R_DEP_MS_WML_vs_ODC_WM.tsv")),
  list(name = "Brain WML vs NAWM (Magliozzi)",
       fp = file.path(OUT_DIR, "Magliozzi_R_DEP_MS_WML_vs_MS_NAWM.tsv")),
  list(name = "T-lineage meta (R/ComBat)",
       fp = file.path(OUT_DIR, "T_lineage_R_combined_DE.tsv")),
  list(name = "Pegram NK8 (R standalone)",
       fp = file.path(OUT_DIR, "Pegram_R_DE_gene.tsv")),
  list(name = "Brain WM RNA (R)",
       fp = file.path(OUT_DIR, "BrainWM_R_RNA_DE.tsv"))
)

rows <- list()
for (a in assays) {
  if (!file.exists(a$fp)) {
    message(sprintf("  SKIP missing: %s", a$fp))
    next
  }
  d <- fread(a$fp)
  d <- as.data.frame(d)
  if (!"gene" %in% colnames(d) && "Gene" %in% colnames(d)) d$gene <- d$Gene
  for (g in CROSS_OMICS) {
    r <- d[d$gene == g, , drop = FALSE]
    if (nrow(r) == 0) next
    r <- r[order(r$P.Value), ][1, ]
    rows[[length(rows) + 1]] <- data.frame(
      assay = a$name, gene = g,
      log2FC = r$logFC, pval = r$P.Value, FDR = r$adj.P.Val,
      stringsAsFactors = FALSE)
  }
}
summ <- do.call(rbind, rows)
summ$is_sig_FDR05 <- summ$FDR < 0.05
summ$is_sig_p05   <- summ$pval < 0.05
fwrite(summ, file.path(FIG_DIR, "CrossAssay_R_summary.tsv"), sep = "\t")
print(summ)


# ---- Render compact per-assay volcanos with 7 stars overlaid ----
make_assay_panel <- function(a) {
  if (!file.exists(a$fp)) return(NULL)
  d <- fread(a$fp); d <- as.data.frame(d)
  if (!"gene" %in% colnames(d) && "Gene" %in% colnames(d)) d$gene <- d$Gene
  d$y <- pmin(-log10(pmax(d$P.Value, 1e-300)), 50)
  d$is_co <- d$gene %in% CROSS_OMICS
  ggplot() +
    geom_point(data = d[!d$is_co, ],
               aes(logFC, y), colour = "#CCCCCC", size = 0.3, alpha = 0.4) +
    geom_point(data = d[d$is_co, ],
               aes(logFC, y, colour = gene), shape = 8, size = 4, stroke = 1.3) +
    scale_colour_manual(values = CO_COLORS, name = "", guide = "none") +
    geom_text_repel(data = d[d$is_co, ],
                    aes(logFC, y, label = gene, colour = gene),
                    size = 2.5, fontface = "bold", max.overlaps = 30) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               colour = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.3) +
    labs(title = a$name) +
    theme_classic(base_size = 8) +
    theme(plot.title = element_text(size = 8, face = "bold"))
}
panels <- Filter(Negate(is.null), lapply(assays, make_assay_panel))
ga <- gridExtra::arrangeGrob(grobs = panels, ncol = 3)
ggsave(file.path(FIG_DIR, "CrossAssay_R_grid.png"), ga,
       width = 18, height = 4 * ceiling(length(panels)/3), dpi = 160,
       limitsize = FALSE)
cat(sprintf("Wrote CrossAssay_R_grid.png (%d panels)\n", length(panels)))

