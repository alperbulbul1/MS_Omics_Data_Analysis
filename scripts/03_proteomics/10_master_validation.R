#!/usr/bin/env Rscript
## 10_master_validation.R  —  generated from notebook spec
## Run: Rscript 10_master_validation.R


## ============================================================
## # 10 — Master cross-omics validation panel
## 
## Final integrating figure. For each of the 7 cross-omics candidates,
## build a heatmap showing log2FC × significance across every R-pipeline
## result. Also computes the triple-validated set
## (CSF protein × Brain RNA × Brain methylation, all FDR<0.05).
## 
## **Outputs**
## - `figures/CrossOmics_Master_Heatmap_R.png`
## - `processed/META/CrossOmics_Master_Table_R.tsv`
## - `processed/META/Triple_Validated_R.tsv`
## ============================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(pheatmap)
})
source("helpers.R")
summ <- fread(file.path(FIG_DIR, "CrossAssay_R_summary.tsv"))
cat(sprintf("Loaded summary: %d rows × %d assays × %d genes\n",
            nrow(summ), uniqueN(summ$assay), uniqueN(summ$gene)))


# ---- pivot: rows = gene, cols = assay, values = log2FC * sign(sig) ----
heat <- dcast(summ, gene ~ assay, value.var = "log2FC", fun.aggregate = mean)
heat_mat <- as.matrix(heat[, -1])
rownames(heat_mat) <- heat$gene

sig <- dcast(summ, gene ~ assay, value.var = "FDR", fun.aggregate = function(x) min(x, na.rm=TRUE))
sig_mat <- as.matrix(sig[, -1])
rownames(sig_mat) <- sig$gene

# Significance star matrix (text overlay)
star_mat <- ifelse(sig_mat < 0.001, "***",
            ifelse(sig_mat < 0.01,  "**",
            ifelse(sig_mat < 0.05,  "*", "")))
star_mat[is.na(star_mat)] <- ""

pheatmap(heat_mat,
         display_numbers = star_mat, fontsize_number = 12,
         color = colorRampPalette(c("#1F4E79","white","#D62828"))(50),
         na_col = "grey90",
         cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Cross-omics adaylar × R-pipeline assayleri (log2FC + FDR stars)",
         filename = file.path(FIG_DIR, "CrossOmics_Master_Heatmap_R.png"),
         width = 12, height = 5)
cat("Wrote master heatmap.\n")
fwrite(summ, file.path(OUT_DIR, "CrossOmics_Master_Table_R.tsv"), sep = "\t")


# ---- Triple-validated set ----
csf_fp  <- file.path(OUT_DIR, "CSF_Astral_CC_results.tsv")
rna_fp  <- file.path(OUT_DIR, "BrainWM_R_RNA_DE.tsv")
meth_fp <- file.path(PROT_ROOT, "processed", "rerun",
                      "BrainWM_meth_genelevel_rerun.tsv")

for (.f in c(csf_fp, rna_fp, meth_fp)) if (!file.exists(.f))
  stop("missing input for the triple-validated set: ", .f)   # was a silent skip
if (TRUE) {
  csf  <- fread(csf_fp);  setnames(csf,  c("logFC","adj.P.Val"),
                                     c("CSF_log2FC","CSF_FDR"), skip_absent = TRUE)
  rna  <- fread(rna_fp);  setnames(rna,  c("logFC","adj.P.Val"),
                                     c("BrainRNA_log2FC","BrainRNA_FDR"),
                                     skip_absent = TRUE)
  meth <- fread(meth_fp); setnames(meth, c("mean_logFC","adj.P.Val"),
                                     c("BrainMeth_log2FC","BrainMeth_FDR"),
                                     skip_absent = TRUE)
  csf$gene <- csf$gene; rna$gene <- rna$gene
  m <- merge(merge(csf[, .(gene, CSF_log2FC, CSF_FDR)],
                   rna[, .(gene, BrainRNA_log2FC, BrainRNA_FDR)],
                   by = "gene"),
             meth[, .(gene = Gene, BrainMeth_log2FC, BrainMeth_FDR)],
             by = "gene")
  triple <- m[CSF_FDR < 0.05 & BrainRNA_FDR < 0.05 & BrainMeth_FDR < 0.05]
  cat(sprintf("Triple-validated (3 layers FDR<0.05): %d genes\n", nrow(triple)))
  print(triple[order(CSF_FDR)])
  fwrite(triple, file.path(OUT_DIR, "Triple_Validated_R.tsv"), sep = "\t")
} else {
  cat("Some input result files missing — run notebooks 01 + 07 first.\n")
}

