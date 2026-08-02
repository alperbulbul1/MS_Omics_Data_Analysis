#!/usr/bin/env Rscript
## 06_pegram_gse32915_de.R  —  generated from notebook spec
## Run: Rscript 06_pegram_gse32915_de.R


## ============================================================
## # 06 — GSE32915 (Pegram 2021 NK8+) standalone limma-style DE
## 
## Single-study DE on Pegram et al. 2021 NK8+ MS-vs-Control microarray
## (Agilent GPL6480, 4 MS × 4 Ctrl with technical reps). Demonstrates that
## single-study power for n=4/4 is insufficient for cross-omics validation —
## a deliberately negative result included as a methodological reference.
## 
## **Outputs**
## - `processed/META/Pegram_R_DE_gene.tsv`
## - `figures/Pegram_R_DE_volcano.png`
## ============================================================

suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(Biobase); library(ggplot2)
})
source("helpers.R")
Sys.setenv(GEOQUERY_CACHE = CACHE_DIR)


cache_fp <- file.path(CACHE_DIR, "GSE32915_eset.rds")
eset <- if (file.exists(cache_fp)) readRDS(cache_fp) else {
  e <- getGEO("GSE32915", GSEMatrix = TRUE, destdir = CACHE_DIR, AnnotGPL = TRUE)[[1]]
  saveRDS(e, cache_fp); e
}
cat(sprintf("GSE32915 expression: %d probes × %d samples\n",
            nrow(exprs(eset)), ncol(exprs(eset))))

expr <- exprs(eset)
if (max(expr, na.rm=TRUE) > 30) expr <- log2(expr + 1)
expr <- normalizeBetweenArrays(expr, method = "quantile")

titles <- pData(eset)$title
group <- ifelse(grepl("Multiple Sclerosis", titles), "MS",
          ifelse(grepl("Control", titles), "Control", NA))
patient <- ifelse(group == "MS",
            paste0("MS", sub(".*Sclerosis ([0-9]+)_.*", "\\1", titles)),
            sub(".*(Control [0-9]+).*", "\\1", titles))

# Collapse tech reps
patient_means <- sapply(unique(patient[!is.na(patient)]), function(p)
                          rowMeans(expr[, which(patient == p), drop = FALSE], na.rm = TRUE))
group_pp <- ifelse(grepl("MS", colnames(patient_means)), "MS", "Control")
fd <- fData(eset)
gcol <- intersect(c("Gene symbol","Gene Symbol","Symbol","GENE_SYMBOL"),
                  colnames(fd))[1]
gene <- fd[[gcol]]
ag <- avereps(patient_means, ID = gene)
ag <- ag[!(rownames(ag) %in% c("", "---") | is.na(rownames(ag))), ]
cat(sprintf("Patient-collapsed gene matrix: %d × %d  (%s)\n",
            nrow(ag), ncol(ag),
            paste(table(group_pp), collapse="/")))


res <- moderated_t_safe(ag, group_pp, "MS", "Control")
res <- res[order(res$P.Value), ]
res$is_cross_omics <- res$gene %in% CROSS_OMICS
res$is_recurring   <- res$gene %in% RECURRING

cat(sprintf("\nTotal: %d  |  P<0.05: %d  |  FDR<0.05: %d  |  FDR<0.001: %d\n",
            nrow(res), sum(res$P.Value < 0.05),
            sum(res$adj.P.Val < 0.05),
            sum(res$adj.P.Val < 0.001)))

cat("\n=== Cross-omics (Pegram standalone, expected: weak/ns) ===\n")
print(subset(res, is_cross_omics)[, c("gene","logFC","P.Value","adj.P.Val")])

write.table(res, file.path(OUT_DIR, "Pegram_R_DE_gene.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

p <- dep_volcano_gg(res,
       title = "Pegram 2021 NK8+ — R limma standalone DE",
       subtitle = sprintf("GSE32915 · %d genes · 4 MS / 4 Ctrl · eBayes(trend, robust)",
                          nrow(res)))
ggsave(file.path(FIG_DIR, "Pegram_R_DE_volcano.png"), p,
       width = 11, height = 8, dpi = 200)
print(p)

