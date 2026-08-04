#!/usr/bin/env Rscript
## 09_inverse_concordance_scan.R  —  generated from notebook spec
## Run: Rscript 09_inverse_concordance_scan.R


## ============================================================
## # 09 — Inverse-concordance discovery scan (RNA ↔ methylation)
## 
## Beyond the original cross-omics 7-gene panel, this notebook scans every
## RNA stratum (from `Transcriptome/results/`) against two methylation
## sources — combined-cohort Stouffer gene-level (nb 05) and mCSEA
## promoter NES (nb 06) — to find ALL genes where:
## 
## - **RNA FDR < 0.05** AND **methylation FDR < 0.05** AND
## - **sign(RNA) ≠ sign(methylation)** (reverse / discordant coupling)
## 
## A gene is ranked by `n_pairings` (in how many RNA-stratum × meth-source
## pairings the inverse-concordance holds). Genes recurring across many
## pairings are the strongest candidates for promoter-methylation-driven
## gene regulation in MS — same biological logic that selected the original
## cross-omics 7-gene panel.
## 
## **Outputs**
## - `results/INVERSE_CONCORDANT_full_pairings.tsv` — every hit, every pairing
## - `results/INVERSE_CONCORDANT_by_gene.tsv` — per-gene aggregate (n_pairings, best FDRs)
## - `figures/09_inverse_concordant_top30.png` — bubble/scatter plot
## ============================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ggplot2); library(ggrepel)
})
source("helpers.R")
TX_DIR <- file.path(PROJ_ROOT, "Transcriptome", "results")


# ---- Load methylation sources ----
meth_combined <- fread(file.path(OUT_DIR, "05_combined_meth_gene.tsv"))
mcsea_prom    <- fread(file.path(OUT_DIR, "06_mCSEA_promoter.tsv"))
mcsea_df <- data.frame(gene = mcsea_prom$gene,
                       mean_logFC = mcsea_prom$NES * 0.1,
                       adj.P.Val  = mcsea_prom$padj,
                       stringsAsFactors = FALSE)

# ---- Load all transcriptome strata ----
tx_files <- c("PBMC"       = "01_pbmc_DE.tsv",
              "T cells"    = "02_tcells_DE.tsv",
              "B cells"    = "03_bcells_DE.tsv",
              "Brain WM"   = "04_brainwm_DE.tsv",
              "Whole blood"= "05_whole_blood_DE.tsv",
              ## IFN-b PBMC is deliberately excluded. 06_pbmc_ifnb_DE.tsv is an
              ## IFN-beta-versus-baseline treatment-response contrast, not MS versus
              ## control, so admitting it to the discovery scan would let a drug effect
              ## define a disease candidate. Seven genes entered the pool through that
              ## pairing alone (ATP6V0E2, EPHX1, HDAC4, MCF2L2, MFAP5, SARM1, ZFP36L1);
              ## none of the tiered candidates did. The stratum is still reported as
              ## treatment context in Figure 2C and in the per-assay tables.
              "Pan-tissue" = "07_pan_tissue_DE.tsv")
tx_list <- lapply(tx_files, function(f) fread(file.path(TX_DIR, f)))
cat(sprintf("Loaded %d TX strata, %d meth sources\n",
            length(tx_list), 2))


# ---- Discovery: scan every (meth source × RNA stratum) pairing ----
discover <- function(meth_df, fc_col, p_col, tx_name, tx_df, tag) {
  m <- data.frame(gene = meth_df$gene, meth_fc = meth_df[[fc_col]],
                  meth_fdr = meth_df[[p_col]], stringsAsFactors = FALSE)
  r <- data.frame(gene = tx_df$gene, rna_fc = tx_df$logFC,
                  rna_fdr = tx_df$adj.P.Val, stringsAsFactors = FALSE)
  j <- merge(m, r, by = "gene")
  j$inverse <- (j$rna_fdr < 0.05) & (j$meth_fdr < 0.05) &
               (sign(j$rna_fc) != sign(j$meth_fc))
  inv <- j[j$inverse, ]
  if (nrow(inv) == 0) return(NULL)
  inv$source <- sprintf("%s::%s", tag, tx_name)
  inv[order(inv$rna_fdr), ]
}

hits <- list()
for (n in names(tx_list)) {
  for (src in list(list(df=meth_combined, tag="Combined-meth(Stouffer)",
                         fc="mean_logFC", p="adj.P.Val"),
                    list(df=mcsea_df,       tag="mCSEA-promoter",
                         fc="mean_logFC", p="adj.P.Val"))) {
    h <- discover(src$df, src$fc, src$p, n, tx_list[[n]], src$tag)
    if (!is.null(h)) hits[[length(hits)+1]] <- h
  }
}
all_df <- do.call(rbind, hits)
cat(sprintf("\n=== %d inverse-concordant rows (%d unique genes) ===\n",
            nrow(all_df), length(unique(all_df$gene))))
fwrite(all_df, file.path(OUT_DIR, "INVERSE_CONCORDANT_full_pairings.tsv"), sep="\t")


# ---- Aggregate per-gene ----
by_gene <- all_df %>%
  group_by(gene) %>%
  summarise(n_pairings    = n(),
            best_rna_fdr  = min(rna_fdr),
            best_meth_fdr = min(meth_fdr),
            best_rna_fc   = rna_fc[which.min(rna_fdr)],
            best_meth_fc  = meth_fc[which.min(meth_fdr)],
            direction     = ifelse(best_rna_fc > 0, "RNA_UP_meth_DOWN", "RNA_DOWN_meth_UP"),
            sources       = paste(unique(source), collapse = "; "),
            is_cross_omics = gene[1] %in% CROSS_OMICS) %>%
  arrange(desc(n_pairings), best_rna_fdr)
cat(sprintf("Genes in >=4 pairings: %d  |  >=3: %d  |  ==2: %d\n",
            sum(by_gene$n_pairings >= 4),
            sum(by_gene$n_pairings >= 3),
            sum(by_gene$n_pairings == 2)))
cat("\n=== TOP 40 ===\n")
print(as.data.frame(head(by_gene[, c("gene","n_pairings","best_rna_fc","best_meth_fc",
                                       "best_rna_fdr","best_meth_fdr","direction","is_cross_omics")], 40)),
      row.names = FALSE)
fwrite(by_gene, file.path(OUT_DIR, "INVERSE_CONCORDANT_by_gene.tsv"), sep="\t")


# ---- Per-gene × per-stratum × per-meth-source breakdown ----
# Tag each row with (rna_stratum, meth_source) parsed from `source`
all_df$rna_stratum  <- sub("^[^:]+::", "", all_df$source)
all_df$meth_source  <- sub("::.+$",      "", all_df$source)

# Wide table: gene × RNA stratum, value = "stars (logFC)"
star <- function(p) ifelse(is.na(p), "",
                     ifelse(p < 0.001, "***",
                     ifelse(p < 0.01,  "**",
                     ifelse(p < 0.05,  "*", ""))))
all_df$cell <- sprintf("%s (%+.2f)", star(all_df$rna_fdr), all_df$rna_fc)

wide_rna <- dcast(as.data.table(all_df),
                  gene ~ rna_stratum,
                  value.var = "cell",
                  fun.aggregate = function(x) paste(unique(x), collapse = " / "))
# Order rows by n_pairings desc
wide_rna <- wide_rna[order(-rowSums(wide_rna[, -1, with = FALSE] != ""))]
fwrite(wide_rna, file.path(OUT_DIR, "INVERSE_CONCORDANT_by_gene_by_stratum.tsv"), sep = "\t")

# Per-source presence table (wide): gene × meth_source (Stouffer / mCSEA)
wide_src <- dcast(as.data.table(all_df),
                  gene ~ meth_source,
                  value.var = "rna_stratum",
                  fun.aggregate = function(x) paste(sort(unique(x)), collapse = ", "))
fwrite(wide_src, file.path(OUT_DIR, "INVERSE_CONCORDANT_by_gene_by_meth_source.tsv"), sep = "\t")

cat("\n=== PER-GENE × PER-STRATUM (top 25 by row-coverage) ===\n\n")
print(as.data.frame(head(wide_rna, 25)), row.names = FALSE)

cat("\n\n=== PER-GENE × PER-METH-SOURCE (top 25) ===\n\n")
print(as.data.frame(head(wide_src, 25)), row.names = FALSE)


# ---- Annotated long table joining gene-aggregate + per-pair details ----
detail <- merge(all_df[, c("gene","rna_stratum","meth_source","rna_fc","rna_fdr",
                            "meth_fc","meth_fdr")],
                by_gene[, c("gene","n_pairings","direction","is_cross_omics")],
                by = "gene")
detail <- detail[order(-detail$n_pairings, detail$gene, detail$rna_fdr), ]
fwrite(detail, file.path(OUT_DIR, "INVERSE_CONCORDANT_detail.tsv"), sep = "\t")
cat(sprintf("\nDetailed long-format rows: %d  ·  saved to INVERSE_CONCORDANT_detail.tsv\n",
            nrow(detail)))


# ---- Figure 1: bubble plot (top 30) ----
plot_df <- head(by_gene, 30) %>%
  mutate(label_rank = rank(-n_pairings) + rank(best_rna_fdr) * 0.001)
p1 <- ggplot(plot_df,
             aes(x = best_rna_fc, y = -log10(best_rna_fdr),
                 size = n_pairings, colour = direction)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(aes(label = gene), size = 3.5, fontface = "bold",
                  max.overlaps = 30, box.padding = 0.5) +
  scale_size(range = c(3, 12), name = "n_pairings\n(strata × meth source)") +
  scale_colour_manual(values = c(RNA_UP_meth_DOWN = "#D62828",
                                  RNA_DOWN_meth_UP = "#1F4E79"),
                     name = "Direction") +
  geom_vline(xintercept = 0, colour = "grey30") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  labs(title = "Top 30 inverse-concordant genes (RNA ↔ methylation)",
       subtitle = sprintf("Discovery scan across %d RNA strata × 2 meth sources",
                          length(tx_list)),
       x = "RNA log2FC (best across strata)",
       y = "-log10 RNA FDR (best across strata)") +
  theme_classic(base_size = 11) + theme(legend.position = "right")
ggsave(file.path(FIG_DIR, "09_inverse_concordant_top30.png"),
       p1, width = 11, height = 8, dpi = 200)
print(p1)


# ---- Figure 2: gene × stratum coverage heatmap ----
# Numeric matrix: rows = top N inverse genes, cols = strata; cell = sign(rna_fc) * -log10(p)
top_n <- 40
top_genes <- head(by_gene$gene, top_n)

mat_long <- as.data.table(all_df)[gene %in% top_genes,
                                   .(score = mean(sign(rna_fc) * -log10(pmax(rna_fdr, 1e-30)))),
                                   by = .(gene, rna_stratum)]
mat <- dcast(mat_long, gene ~ rna_stratum, value.var = "score", fill = 0)
mat <- mat[match(top_genes, mat$gene), ]
mat_m <- as.matrix(mat[, -1, with = FALSE]); rownames(mat_m) <- mat$gene

# Significance star matrix
sig_long <- as.data.table(all_df)[gene %in% top_genes,
                                   .(p = min(rna_fdr, na.rm = TRUE)),
                                   by = .(gene, rna_stratum)]
sig_w <- dcast(sig_long, gene ~ rna_stratum, value.var = "p", fill = NA_real_)
sig_w <- sig_w[match(top_genes, sig_w$gene), ]
sig_mat <- as.matrix(sig_w[, -1, with = FALSE])
star_mat <- ifelse(is.na(sig_mat), "",
            ifelse(sig_mat < 0.001, "***",
            ifelse(sig_mat < 0.01,  "**",
            ifelse(sig_mat < 0.05,  "*", ""))))

# Row annotation: cross-omics flag + n_pairings
row_ann <- data.frame(
  n_pairings = by_gene$n_pairings[match(top_genes, by_gene$gene)],
  cross_omics = ifelse(by_gene$is_cross_omics[match(top_genes, by_gene$gene)], "yes", "no"),
  direction = by_gene$direction[match(top_genes, by_gene$gene)],
  row.names = top_genes
)

library(pheatmap)
pheatmap(mat_m,
         display_numbers = star_mat, fontsize_number = 10,
         color = colorRampPalette(c("#1F4E79","#FFFFFF","#D62828"))(50),
         na_col = "grey95",
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_row = row_ann,
         annotation_colors = list(
           cross_omics = c(yes = "#D62828", no = "#CCCCCC"),
           direction   = c(RNA_UP_meth_DOWN = "#D62828",
                            RNA_DOWN_meth_UP = "#1F4E79")),
         main = sprintf("Top %d inverse-concordant genes — RNA stratum coverage\n(cell = sign(RNA_logFC) × -log10 FDR; stars = sig stratum)",
                        top_n),
         filename = file.path(FIG_DIR, "09_inverse_concordant_stratum_heatmap.png"),
         width = 11, height = 12)
cat("Wrote 09_inverse_concordant_stratum_heatmap.png\n")


# ---- Figure 3: per-source presence matrix (Stouffer vs mCSEA) ----
src_long <- as.data.table(all_df)[gene %in% top_genes,
                                   .(present = 1),
                                   by = .(gene, meth_source, rna_stratum)]
src_long$combo <- paste(src_long$meth_source, src_long$rna_stratum, sep = "\n")
mat3 <- dcast(src_long, gene ~ combo, value.var = "present", fill = 0)
mat3 <- mat3[match(top_genes, mat3$gene), ]
mat3_m <- as.matrix(mat3[, -1, with = FALSE]); rownames(mat3_m) <- mat3$gene

pheatmap(mat3_m,
         color = colorRampPalette(c("#F5F5F5","#1F4E79"))(2),
         legend_breaks = c(0, 1), legend_labels = c("ns", "sig FDR<0.05"),
         cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = ifelse(mat3_m == 1, "x", ""),
         fontsize_number = 11, color_number = "white",
         main = sprintf("Inverse-concordance presence per (meth_source × RNA_stratum) — top %d", top_n),
         filename = file.path(FIG_DIR, "09_inverse_concordant_presence_matrix.png"),
         width = 12, height = 12)
cat("Wrote 09_inverse_concordant_presence_matrix.png\n")

