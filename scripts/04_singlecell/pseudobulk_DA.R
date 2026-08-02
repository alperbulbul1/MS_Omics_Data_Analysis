#!/usr/bin/env Rscript
# DIFFERENTIAL ABUNDANCE (cell-type composition), analysed SEPARATELY from differential state.
# Question: does the PROPORTION of a cell type differ between MS and HC?  (propeller-style:
# variance-stabilising transform of proportions + limma with the paired design.)
# This is conceptually distinct from pseudobulk DE and is reported separately, as recommended.
# It also matters for interpretation: a composition shift can produce a bulk-tissue fold-change
# with no within-cell-type expression change (and vice versa).
suppressPackageStartupMessages({ library(limma) })
DIR <- "__MS_GEO_ROOT__/Poster_v2/figures/pseudobulk_proper"

for (tag in c("coarse","fine")) {
  C <- read.csv(file.path(DIR,sprintf("PBC_%s_coldata.csv",tag)), check.names=FALSE)
  # complete MS/HC pairs only
  fp <- names(which(tapply(C$condition, C$batch_pair, function(z) all(c("MS","HC") %in% z))))
  C  <- C[C$batch_pair %in% fp, ]
  # proportions per donor
  tot <- tapply(C$n_cells, C$donor, sum)
  C$prop <- C$n_cells / tot[C$donor]
  P <- tapply(C$prop, list(C$cell_type, C$donor), function(z) if(length(z)) z[1] else 0)
  P[is.na(P)] <- 0
  don <- colnames(P)
  md  <- C[match(don, C$donor), c("donor","condition","batch_pair","sex","age")]
  grp <- factor(md$condition, levels=c("HC","MS")); pr <- factor(as.character(md$batch_pair))
  design <- model.matrix(~ pr + grp); cc <- ncol(design)
  # propeller's variance-stabilising transform for proportions
  Y <- asin(sqrt(P))
  fit <- eBayes(lmFit(Y, design), robust=TRUE)
  tt  <- topTable(fit, coef=cc, number=Inf, sort.by="P")
  cat(sprintf("\n=========== DIFFERENTIAL ABUNDANCE (%s: %d cell types, %d donors, %d pairs) ===========\n",
              tag, nrow(P), ncol(P), nlevels(pr)))
  cat("  cell_type        meanProp_MS  meanProp_HC     logFC(asin)     p        FDR\n")
  for (g in rownames(tt)) {
    ms <- mean(P[g, md$condition=="MS"]); hc <- mean(P[g, md$condition=="HC"])
    r <- tt[g,]
    cat(sprintf("  %-16s %10.4f %12.4f %13.4f  %.2e  %6.4f%s\n",
                g, ms, hc, r$logFC, r$P.Value, r$adj.P.Val, ifelse(r$adj.P.Val<0.05," **","")))
  }
  write.csv(cbind(cell_type=rownames(tt), tt), file.path(DIR,sprintf("DA_%s.csv",tag)), row.names=FALSE)
}
cat("\ndone -> DA_coarse.csv / DA_fine.csv\n")
