#!/usr/bin/env Rscript
## 05_t_lineage_meta.R  —  generated from notebook spec
## Run: Rscript 05_t_lineage_meta.R


## ============================================================
## # 05 — T-lineage microarray meta (GSE32915 + GSE78244)
## 
## Cross-study T-lineage meta-analysis:
## 
## | Series   | Author        | Cell type | Platform   | n            |
## |----------|--------------|-----------|------------|--------------|
## | GSE32915 | Pegram 2021  | NK8+      | GPL6480    | 4 MS / 4 Ctrl (tech reps) |
## | GSE78244 | Hellberg 2016| CD4+ T    | GPL17077   | 14 MS / 14 Ctrl (unstim subset) |
## 
## **Pipeline (R/Bioconductor canonical):**
## 1. `GEOquery::getGEO` → eSet per study (cached in `_cache/`)
## 2. log2 transform if needed, `limma::normalizeBetweenArrays(method="quantile")`
## 3. `limma::avereps` collapse to gene level
## 4. Inner-join on common gene symbols
## 5. `sva::ComBat(batch=study, mod=~group)`
## 6. `limma::lmFit + makeContrasts(MS - Control) + eBayes(trend, robust)`
## 7. `EnhancedVolcano` + `pheatmap` of cross-omics × samples
## 
## **Outputs**
## - `processed/META/T_lineage_R_combined_DE.tsv`
## - `figures/T_lineage_R_pca.png`
## - `figures/T_lineage_R_volcano.png`
## - `figures/T_lineage_R_heatmap.png`
## ============================================================

suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(sva); library(Biobase)
  library(ggplot2); library(dplyr); library(gridExtra); library(pheatmap)
})
source("helpers.R")
geo_cache <- CACHE_DIR
Sys.setenv(GEOQUERY_CACHE = geo_cache)
options(timeout = 600)


load_gse_eset <- function(id) {
  cache_fp <- file.path(geo_cache, paste0(id, "_eset.rds"))
  if (file.exists(cache_fp)) {
    cat(sprintf("  using cache: %s\n", basename(cache_fp)))
    return(readRDS(cache_fp))
  }
  e <- getGEO(id, GSEMatrix = TRUE, destdir = geo_cache, AnnotGPL = TRUE)[[1]]
  saveRDS(e, cache_fp)
  e
}
prep_eset <- function(eset, label) {
  expr <- exprs(eset)
  if (max(expr, na.rm=TRUE) > 30) {
    cat(sprintf("    %s: applying log2(x+1)\n", label))
    expr <- log2(expr + 1)
  }
  expr <- normalizeBetweenArrays(expr, method = "quantile")
  fd <- fData(eset)
  gcol <- intersect(c("Gene symbol","Gene Symbol","Symbol","GENE_SYMBOL"),
                    colnames(fd))[1]
  if (is.na(gcol)) stop(sprintf("no gene-symbol col in %s fData", label))
  list(expr = expr, gene = fd[[gcol]], pdata = pData(eset))
}


cat("== Loading GSE32915 (Pegram NK8+) ==\n")
e1 <- load_gse_eset("GSE32915")
p1 <- prep_eset(e1, "GSE32915")
titles_1 <- p1$pdata$title
group_1 <- ifelse(grepl("Multiple Sclerosis", titles_1), "MS",
            ifelse(grepl("Control", titles_1), "Control", NA))
patient_1 <- ifelse(group_1 == "MS",
            paste0("MS", sub(".*Sclerosis ([0-9]+)_.*", "\\1", titles_1)),
            sub(".*(Control [0-9]+).*", "\\1", titles_1))
# Collapse technical replicates per patient
patient_keep <- patient_1[!is.na(patient_1)]
expr_1 <- p1$expr
patient_means <- sapply(unique(patient_keep), function(pp) {
  cols <- which(patient_1 == pp)
  rowMeans(expr_1[, cols, drop = FALSE], na.rm = TRUE)
})
group_per_patient <- ifelse(grepl("MS", colnames(patient_means)), "MS", "Control")
ag1 <- avereps(patient_means, ID = p1$gene)
ag1 <- ag1[!(rownames(ag1) %in% c("", "---") | is.na(rownames(ag1))), ]
cat(sprintf("  collapsed: %d genes × %d patients (%s)\n",
            nrow(ag1), ncol(ag1),
            paste(table(group_per_patient), collapse="/")))


cat("\n== Loading GSE78244 (Hellberg CD4 T) ==\n")
e2 <- load_gse_eset("GSE78244")
titles_2_all <- pData(e2)$title
keep <- grepl("unstimulated", titles_2_all)
e2s <- e2[, keep]
p2 <- prep_eset(e2s, "GSE78244_unstim")
group_2 <- ifelse(grepl("patient", p2$pdata$title), "MS",
            ifelse(grepl("control", p2$pdata$title), "Control", NA))
ag2 <- avereps(p2$expr, ID = p2$gene)
ag2 <- ag2[!(rownames(ag2) %in% c("", "---") | is.na(rownames(ag2))), ]
cat(sprintf("  unstim: %d genes × %d samples (%s)\n",
            nrow(ag2), ncol(ag2),
            paste(table(group_2), collapse="/")))


common <- intersect(rownames(ag1), rownames(ag2))
cat(sprintf("Common genes: %d\n", length(common)))
combined <- cbind(ag1[common, ], ag2[common, ])
meta <- data.frame(
  sample = colnames(combined),
  study  = c(rep("GSE32915", ncol(ag1)), rep("GSE78244", ncol(ag2))),
  group  = c(group_per_patient, group_2),
  stringsAsFactors = FALSE)

mod <- model.matrix(~ group, data = meta)
combat_mat <- ComBat(dat = combined, batch = meta$study, mod = mod,
                     par.prior = TRUE)
cat(sprintf("Combined: %d × %d\n", nrow(combat_mat), ncol(combat_mat)))


# ---- PCA before vs after ----
plot_pca <- function(mat, ttl) {
  pc <- prcomp(t(mat), scale. = TRUE)
  ggplot(data.frame(PC1 = pc$x[,1], PC2 = pc$x[,2],
                    study = meta$study, group = meta$group),
         aes(PC1, PC2, colour = study, shape = group)) +
    geom_point(size = 3, alpha = 0.85) +
    scale_colour_manual(values = c(GSE32915 = "#1F4E79", GSE78244 = "#D62828")) +
    theme_classic(base_size = 11) + ggtitle(ttl) +
    theme(legend.position = "bottom")
}
gg <- gridExtra::arrangeGrob(plot_pca(combined,   "Before ComBat"),
                              plot_pca(combat_mat, "After ComBat"), ncol = 2)
ggsave(file.path(FIG_DIR, "T_lineage_R_pca.png"), gg,
       width = 13, height = 5.5, dpi = 200)
grid::grid.draw(gg)


# ---- limma DE on ComBat-corrected ----
design <- model.matrix(~ 0 + factor(meta$group, levels = c("Control","MS")))
colnames(design) <- c("Control", "MS")
fit <- lmFit(combat_mat, design)
ct  <- makeContrasts(MS_vs_Ctrl = MS - Control, levels = design)
fit2 <- contrasts.fit(fit, ct)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
res <- topTable(fit2, coef = "MS_vs_Ctrl", number = Inf, sort.by = "P")
res$gene <- rownames(res)
res$is_cross_omics <- res$gene %in% CROSS_OMICS
res$is_recurring   <- res$gene %in% RECURRING

cat(sprintf("Total: %d  |  P<0.05: %d  |  FDR<0.05: %d  |  FDR<0.001: %d\n",
            nrow(res), sum(res$P.Value < 0.05),
            sum(res$adj.P.Val < 0.05),
            sum(res$adj.P.Val < 0.001)))

cat("\n=== Cross-omics in T-lineage meta ===\n")
print(subset(res, is_cross_omics)[, c("gene","logFC","P.Value","adj.P.Val")])

write.table(res, file.path(OUT_DIR, "T_lineage_R_combined_DE.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


p <- dep_volcano_gg(res,
       title = "T-lineage microarray meta (R/ComBat + limma)",
       subtitle = sprintf("GSE32915 (NK8) + GSE78244 (CD4) · %d genes · ComBat batch-corrected",
                          nrow(res)))
ggsave(file.path(FIG_DIR, "T_lineage_R_volcano.png"), p,
       width = 11, height = 8, dpi = 200)
print(p)

# Heatmap of cross-omics + top hits across samples
top_genes <- unique(c(CROSS_OMICS,
                       head(res$gene[res$adj.P.Val < 0.05], 30)))
top_genes <- intersect(top_genes, rownames(combat_mat))
ann_col <- data.frame(study = meta$study, group = meta$group,
                      row.names = meta$sample)
ann_row <- data.frame(category = ifelse(top_genes %in% CROSS_OMICS, "Cross-omics aday",
                              ifelse(top_genes %in% RECURRING, "Recurring", "DE top")),
                      row.names = top_genes)
pheatmap(combat_mat[top_genes, ],
         annotation_col = ann_col, annotation_row = ann_row,
         scale = "row", show_colnames = FALSE,
         color = colorRampPalette(c("#1F4E79","white","#D62828"))(50),
         filename = file.path(FIG_DIR, "T_lineage_R_heatmap.png"),
         width = 9, height = 8)
cat("Saved heatmap.\n")

