#!/usr/bin/env Rscript
## 04_tcells_remission_meth_dmp.R  —  generated from notebook spec
## Run: Rscript 04_tcells_remission_meth_dmp.R


## ============================================================
## # 04 — T cells remission context
## 
## R/limma rerun of methylation stratum `label_context_case_control_t_cells_remission` from
## `Stratified_Analyses/Methylation/`. Uses the already batch-corrected
## M-value matrix (Batch_Corrected_M.csv) as input.
## 
## **Pipeline**: `limma::lmFit + eBayes(robust)` on M-values → BH FDR →
## probe-level volcano → annotate to gene → Stouffer aggregate → gene-level
## volcano with cross-omics overlay.
## 
## **Outputs**
## - `results/04_tcells_remission_meth_DMP.tsv` (probe level)
## - `results/04_tcells_remission_meth_gene.tsv` (gene level via Stouffer)
## - `figures/04_tcells_remission_meth_volcano_probe.png`, `04_tcells_remission_meth_volcano_gene.png`
## 
## **T cells in remission context**: 6 MS-remission / 8 HC.
## Very small n; Python pipeline found promoter-level signals
## (KRTAP12-3/SLFN12/MYCBPAP).
## ============================================================

suppressPackageStartupMessages({
  library(limma); library(data.table); library(ggplot2); library(ggrepel); library(dplyr)
})
source("helpers.R")


s <- load_meth_stratum("label_context_case_control_t_cells_remission")
cat(sprintf("Stratum: %s  ·  %d probes × %d samples\n",
            s$stratum, nrow(s$mat), ncol(s$mat)))
print(table(s$groups))


# Optional batch covariate (study) if multiple datasets
batches <- NULL
if ("dataset" %in% colnames(s$meta) && uniqueN(s$meta$dataset[match(colnames(s$mat), s$meta$sample_id)]) > 1) {
  batches <- s$meta$dataset[match(colnames(s$mat), s$meta$sample_id)]
  cat("Using study as batch covariate. studies:\n"); print(table(batches))
}
dmp <- run_limma_meth(s$mat, s$groups, batch = batches,
                       group_a = "MS", group_b = "HC")
cat(sprintf("\nDMP probes: %d  ·  P<0.05: %d  ·  FDR<0.05: %d\n",
            nrow(dmp), sum(dmp$P.Value < 0.05),
            sum(dmp$adj.P.Val < 0.05)))
cat("\n=== Top 10 probes ===\n")
print(head(dmp[, c("Probe","logFC","P.Value","adj.P.Val")], 10))


# Probe -> gene annotation (lazy install if missing)
have_ann <- requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19",
                              quietly = TRUE)
if (!have_ann) {
  message("Installing 450k annotation package...")
  BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19",
                       update = FALSE, ask = FALSE)
}
p2g <- annotate_probes_to_genes(dmp$Probe)
cat(sprintf("Annotated %d/%d probes to gene\n",
            sum(!is.na(p2g$gene)), nrow(p2g)))

gene_level <- probe_to_gene_stouffer(dmp, p2g)
gene_level$is_cross_omics <- gene_level$gene %in% CROSS_OMICS
cat(sprintf("Gene-level: %d genes  ·  FDR<0.05: %d  ·  FDR<0.001: %d\n",
            nrow(gene_level),
            sum(gene_level$adj.P.Val < 0.05),
            sum(gene_level$adj.P.Val < 0.001)))

cat("\n=== Cross-omics gene-level signals ===\n")
print(subset(gene_level, is_cross_omics)[, c("gene","n_probes","mean_logFC","P.Value","adj.P.Val")])

cat("\n=== Top 10 gene-level hits ===\n")
print(head(gene_level[, c("gene","n_probes","mean_logFC","P.Value","adj.P.Val")], 10))

write.table(dmp, file.path(OUT_DIR, "04_tcells_remission_meth_DMP.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(gene_level, file.path(OUT_DIR, "04_tcells_remission_meth_gene.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# Volcano (probe-level)
dmp_lbl <- merge(dmp, p2g, by = "Probe", all.x = TRUE)
dmp_lbl$gene[is.na(dmp_lbl$gene)] <- dmp_lbl$Probe[is.na(dmp_lbl$gene)]
p_probe <- meth_volcano_gg(dmp_lbl,
              title = sprintf("%s — probe-level DMP  (R/limma)", s$stratum),
              subtitle = sprintf("%d probes · n=%d  %s vs %s",
                                 nrow(dmp_lbl), ncol(s$mat),
                                 "MS", "HC"))
ggsave(file.path(FIG_DIR, "04_tcells_remission_meth_volcano_probe.png"), p_probe,
       width = 11, height = 8, dpi = 200)

# Volcano (gene-level)
gene_level$logFC     <- gene_level$mean_logFC
p_gene <- meth_volcano_gg(gene_level,
              title = sprintf("%s — gene-level (Stouffer)", s$stratum),
              subtitle = sprintf("%d genes · %s vs %s",
                                 nrow(gene_level), "MS", "HC"))
ggsave(file.path(FIG_DIR, "04_tcells_remission_meth_volcano_gene.png"), p_gene,
       width = 11, height = 8, dpi = 200)
print(p_gene)

