#!/usr/bin/env Rscript
## 07_brainwm_rna_vs_meth.R  —  generated from notebook spec
## Run: Rscript 07_brainwm_rna_vs_meth.R


## ============================================================
## # 07 — Brain WM RNA × methylation inverse-concordance scan
## 
## Cross-omics scan: take all genes with both an RNA logFC AND a
## methylation gene-level logFC in the Brain WM stratum, then identify the
## *inverse-concordant* set (RNA↑/meth↓ or RNA↓/meth↑) at joint FDR<0.05.
## 
## Uses the RNA result from `BrainWM_R_RNA_DE.tsv` (produced by
## Proteomics/r_notebooks/07_brainwm_rna_meth_rerun) and the meth
## gene-level rerun TSV.
## 
## **Outputs**
## - `results/07_BrainWM_RNA_vs_Meth_concordance.tsv`
## - `figures/07_BrainWM_RNA_vs_Meth_concordance.png`
## ============================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(dplyr)
})
source("helpers.R")


# Locate prior R-pipeline outputs
rna_fp <- file.path(PROJ_ROOT, "Proteomics", "processed", "META",
                     "BrainWM_R_RNA_DE.tsv")
meth_fp <- file.path(PROJ_ROOT, "Proteomics", "processed", "rerun",
                      "BrainWM_meth_genelevel_rerun.tsv")
stopifnot(file.exists(rna_fp), file.exists(meth_fp))

rna <- fread(rna_fp); meth <- fread(meth_fp)
setnames(meth, c("Gene","mean_logFC","P.Value","adj.P.Val"),
              c("gene","logFC","P.Value","adj.P.Val"), skip_absent = TRUE)
cat(sprintf("RNA:  %d genes  ·  meth: %d genes\n", nrow(rna), nrow(meth)))


common <- intersect(rna$gene, meth$gene)
m <- data.frame(
  gene = common,
  rna  = rna$logFC[match(common, rna$gene)],
  meth = meth$logFC[match(common, meth$gene)],
  rna_fdr  = rna$adj.P.Val[match(common, rna$gene)],
  meth_fdr = meth$adj.P.Val[match(common, meth$gene)],
  is_co = common %in% CROSS_OMICS,
  stringsAsFactors = FALSE)
m$both_sig <- m$rna_fdr < 0.05 & m$meth_fdr < 0.05
m$inverse  <- m$both_sig & (sign(m$rna) != sign(m$meth))
cat(sprintf("Common: %d  ·  both sig: %d  ·  inverse-concordant: %d\n",
            nrow(m), sum(m$both_sig, na.rm = TRUE),
            sum(m$inverse, na.rm = TRUE)))
print(m[m$inverse | m$is_co, ])

fwrite(m, file.path(OUT_DIR, "07_BrainWM_RNA_vs_Meth_concordance.tsv"),
       sep = "\t")


p <- ggplot(m, aes(rna, meth)) +
  geom_point(data = subset(m, !both_sig), colour = "grey80", size = 0.6) +
  geom_point(data = subset(m, both_sig & !inverse),
             colour = "#FFC107", size = 1.4, alpha = 0.7) +
  geom_point(data = subset(m, inverse),
             colour = "#1F4E79", size = 2.0) +
  geom_point(data = subset(m, is_co),
             shape = 8, colour = "#D62828", size = 5, stroke = 1.4) +
  geom_text_repel(data = subset(m, is_co | inverse),
                  aes(label = gene), size = 3.5, fontface = "bold",
                  colour = "#D62828") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed",
              colour = "#D62828") +
  labs(title = "Brain WM RNA × methylation (inverse-concordance scan)",
       subtitle = sprintf("%d common genes  ·  %d inverse-concordant at joint FDR<0.05",
                          nrow(m), sum(m$inverse, na.rm = TRUE)),
       x = "RNA logFC (MS / HC)", y = "Methylation mean logFC (MS / HC)") +
  theme_classic(base_size = 11)
ggsave(file.path(FIG_DIR, "07_BrainWM_RNA_vs_Meth_concordance.png"), p,
       width = 9, height = 8, dpi = 200)
print(p)

