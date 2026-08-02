#!/usr/bin/env Rscript
## 10_inverse_proteomics_validation  —  generated from notebook spec


## ============================================================
## # 10 — Inverse-concordance validation in proteomics (R)
## 
## Take the **inverse-concordant gene list** discovered in notebook 09 and
## look each gene up in the R proteomics pipeline outputs:
## 
## - CSF Astral (R/DEP-equivalent) — single largest cohort
## - CSF timsTOF (R/DEP-equivalent)
## - CSF combined (Astral+timsTOF, sva::ComBat)
## - T-lineage meta (GSE32915+GSE78244, ComBat+limma)
## - Magliozzi 2026 brain proteomics: 4 contrasts (CTX, NAWM, WML-vs-WM, WML-vs-NAWM)
## - Pegram NK8 standalone
## 
## For every (gene × proteomics assay) we record log2FC + FDR. Produces:
## 
## **Outputs**
## - `results/INV_proteomics_validation_long.tsv`  (full long table)
## - `results/INV_proteomics_validation_by_gene.tsv` (per-gene aggregate: n_assays_sig)
## - `figures/10_inverse_proteomics_heatmap.png`     (gene × assay heatmap)
## ============================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ggplot2); library(pheatmap); library(reshape2)
})
source("helpers.R")
META <- file.path(PROJ_ROOT, "Proteomics", "processed", "META")


# Load inverse-concordant gene list (top 40 from nb09)
inv <- fread(file.path(OUT_DIR, "INVERSE_CONCORDANT_by_gene.tsv"))
inv <- inv[order(-n_pairings, best_rna_fdr)][1:40]
gene_panel <- inv$gene
cat(sprintf("Loading %d inverse-concordant genes\n", length(gene_panel)))
print(head(inv, 10))


# Load all R proteomics results
proteomics_assays <- list(
  list(name = "CSF Astral (R/DEP)",
       fp = file.path(META, "CSF_Astral_R_DEP_results.tsv")),
  list(name = "CSF timsTOF (R/DEP)",
       fp = file.path(META, "CSF_timsTOF_R_DEP_results.tsv")),
  list(name = "CSF combined (R/ComBat)",
       fp = file.path(META, "CSF_combined_R_ComBat_DE.tsv")),
  list(name = "T-lineage meta (R)",
       fp = file.path(META, "T_lineage_R_combined_DE.tsv")),
  list(name = "Pegram NK8 (R)",
       fp = file.path(META, "Pegram_R_DE_gene.tsv")),
  list(name = "Brain CTX (Magliozzi)",
       fp = file.path(META, "Magliozzi_R_DEP_MS_CTX_vs_ODC_CTX.tsv")),
  list(name = "Brain NAWM (Magliozzi)",
       fp = file.path(META, "Magliozzi_R_DEP_MS_NAWM_vs_ODC_WM.tsv")),
  list(name = "Brain WML vs ctrl (Magliozzi)",
       fp = file.path(META, "Magliozzi_R_DEP_MS_WML_vs_ODC_WM.tsv")),
  list(name = "Brain WML vs NAWM (Magliozzi)",
       fp = file.path(META, "Magliozzi_R_DEP_MS_WML_vs_MS_NAWM.tsv"))
)

rows <- list()
for (a in proteomics_assays) {
  if (!file.exists(a$fp)) { message("SKIP missing: ", a$fp); next }
  d <- as.data.frame(fread(a$fp))
  if (!"gene" %in% colnames(d) && "Gene" %in% colnames(d)) d$gene <- d$Gene
  for (g in gene_panel) {
    r <- d[d$gene == g, , drop = FALSE]
    if (nrow(r) == 0) next
    r <- r[order(r$P.Value), ][1, ]
    rows[[length(rows)+1]] <- data.frame(
      assay = a$name, gene = g,
      log2FC = r$logFC, P.Value = r$P.Value, adj.P.Val = r$adj.P.Val,
      stringsAsFactors = FALSE)
  }
}
long <- do.call(rbind, rows)
long$sig_FDR05  <- long$adj.P.Val < 0.05
long$sig_p05    <- long$P.Value   < 0.05
cat(sprintf("\nLong table: %d rows  ·  %d genes × %d assays\n",
            nrow(long), uniqueN(long$gene), uniqueN(long$assay)))
fwrite(long, file.path(OUT_DIR, "INV_proteomics_validation_long.tsv"), sep="\t")


# Per-gene aggregate
by_gene <- long %>%
  group_by(gene) %>%
  summarise(n_assays_tested  = n(),
            n_assays_sig_FDR = sum(sig_FDR05),
            n_assays_sig_p   = sum(sig_p05),
            best_assay       = assay[which.min(adj.P.Val)],
            best_log2FC      = log2FC[which.min(adj.P.Val)],
            best_FDR         = min(adj.P.Val, na.rm = TRUE),
            sig_assays       = paste(unique(assay[sig_FDR05]), collapse = "; ")) %>%
  arrange(desc(n_assays_sig_FDR), best_FDR)
# Merge with inv panel for n_pairings + direction
by_gene <- merge(by_gene,
                 inv[, c("gene","n_pairings","direction","is_cross_omics",
                         "best_rna_fc","best_rna_fdr","best_meth_fc","best_meth_fdr")],
                 by = "gene", all.x = TRUE)
by_gene <- by_gene[order(-by_gene$n_assays_sig_FDR, by_gene$best_FDR), ]
fwrite(by_gene, file.path(OUT_DIR, "INV_proteomics_validation_by_gene.tsv"), sep="\t")

cat("\n=== Top 25 genes (ranked by n_assays_sig_FDR) ===\n")
print(as.data.frame(head(by_gene[, c("gene","n_pairings","n_assays_tested","n_assays_sig_FDR",
                                      "best_log2FC","best_FDR","best_assay")], 25)),
      row.names = FALSE)


# Heatmap: gene × proteomics assay (log2FC + sig stars)
mat_w <- dcast(as.data.table(long), gene ~ assay, value.var = "log2FC",
               fun.aggregate = mean)
sig_w <- dcast(as.data.table(long), gene ~ assay, value.var = "adj.P.Val",
               fun.aggregate = function(x) suppressWarnings(min(x, na.rm = TRUE)))
mat_w <- as.data.table(mat_w); sig_w <- as.data.table(sig_w)
mat_m <- as.matrix(mat_w[, -1, with = FALSE]); rownames(mat_m) <- mat_w$gene
sig_m <- as.matrix(sig_w[, -1, with = FALSE]); rownames(sig_m) <- sig_w$gene
star <- ifelse(is.na(sig_m), "",
        ifelse(sig_m < 0.001, "***",
        ifelse(sig_m < 0.01,  "**",
        ifelse(sig_m < 0.05,  "*", ""))))

# Order rows by best_FDR
ord <- by_gene$gene[by_gene$gene %in% rownames(mat_m)]
mat_m <- mat_m[ord, , drop = FALSE]
star  <- star [ord, , drop = FALSE]

# Cap extreme values for colour scale
mat_m_capped <- pmin(pmax(mat_m, -2), 2)

row_ann <- data.frame(
  n_RNA_meth_pairings = by_gene$n_pairings[match(ord, by_gene$gene)],
  direction = by_gene$direction[match(ord, by_gene$gene)],
  cross_omics = ifelse(by_gene$is_cross_omics[match(ord, by_gene$gene)], "yes","no"),
  row.names = ord
)

pheatmap(mat_m_capped,
         display_numbers = star, fontsize_number = 9,
         color = colorRampPalette(c("#1F4E79","white","#D62828"))(50),
         na_col = "grey92",
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_row = row_ann,
         annotation_colors = list(
           cross_omics = c(yes = "#D62828", no = "#CCCCCC"),
           direction = c(RNA_UP_meth_DOWN = "#D62828",
                          RNA_DOWN_meth_UP = "#1F4E79")),
         main = sprintf("Inverse-concordant genes (n=%d) × R-proteomics assays\n(log2FC; stars = FDR sig)",
                        nrow(mat_m)),
         filename = file.path(FIG_DIR, "10_inverse_proteomics_heatmap.png"),
         width = 13, height = 12)
cat("Wrote 10_inverse_proteomics_heatmap.png\n")

