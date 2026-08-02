#!/usr/bin/env Rscript
# AGGREGATION SENSITIVITY: does the pseudobulk aggregation unit change the conclusions?
#
#   S1  SUM of raw integer UMI counts   -> edgeR-QL / voom      (muscat standard; already run)
#   S2  MEAN of per-cell CP10K (= CPM/100), log2 -> limma-trend  (normalise each cell, then average;
#                                                                 every cell weighted equally)
#   S3  SUM of per-cell CP10K, log2            -> limma-trend    (normalise then sum)
#
# For 10x 3' data TPM is not identifiable (no usable effective length from 3'-biased reads), so
# CP10K/CPM is the correct normalised unit; CPM = CP10K * 100 and the log2 fold-changes are
# identical under either scaling.
# Same paired design, same granularities, same three FDR scopes as the count-based run.
suppressPackageStartupMessages({ library(limma) })
DIR <- "__MS_GEO_ROOT__/Poster_v2/figures/pseudobulk_proper"
T1    <- c("ITGB2","LXN","CD79B","IKZF1","SH3BP4"); SUG <- "HLA-E"
AUX   <- c("CASP6","CASP8","DGKQ","MX1","IFIT1","NUP210","RUNX3")
PANEL <- c(T1,SUG,AUX,"CTSZ","CHL1","THRB","ITGAL","IFI44L","RPAP2","SLAMF1","PCNP",
           "STAT3","TYK2","ICAM1","MOSPD3","FOXP3")
BIO <- read.csv(file.path(DIR,"gene_biotypes.csv"), row.names=1)

complete_pairs <- function(s) {
  fp <- names(which(tapply(s$condition, s$batch_pair, function(z) all(c("MS","HC") %in% z))))
  s[s$batch_pair %in% fp, ]
}

fit_ct <- function(Ml, s, ct, label) {
  grp <- factor(s$condition, levels=c("HC","MS"))
  if (sum(grp=="MS") < 3 || sum(grp=="HC") < 3) return(NULL)
  pr <- factor(as.character(s$batch_pair)); if (nlevels(pr) < 3) return(NULL)
  design <- model.matrix(~ pr + grp); if (nrow(design)-ncol(design) < 2) return(NULL)
  E <- Ml[, s$sample, drop=FALSE]
  keep <- rowMeans(E > 0) >= 0.5 | rownames(E) %in% PANEL
  # protein-coding filter (biotype from mygene.info), as in the count-based run
  pc <- rownames(E) %in% rownames(BIO)[BIO$biotype=="protein-coding"]
  keep <- keep & (pc | rownames(E) %in% PANEL)
  y <- E[keep, , drop=FALSE]; if (nrow(y) < 50) return(NULL)
  aw  <- arrayWeights(y, design=design)
  fit <- eBayes(lmFit(y, design, weights=aw), trend=TRUE, robust=TRUE)
  tt  <- topTable(fit, coef=ncol(design), number=Inf, sort.by="none")
  out <- data.frame(gene=rownames(tt), logFC=tt$logFC, PValue=tt$P.Value)
  out$FDR_local <- p.adjust(out$PValue,"BH")
  cbind(scheme=label, cell_type=ct, n_ms=sum(grp=="MS"), n_hc=sum(grp=="HC"),
        n_pairs=nlevels(pr), out, row.names=NULL)
}

run_level <- function(tag, unit, label, naive) {
  f <- file.path(DIR, sprintf("PBN_%s_%s.csv", tag, unit))
  if (!file.exists(f)) { cat("missing:",f,"\n"); return(NULL) }
  M <- as.matrix(read.csv(f, row.names=1, check.names=FALSE))
  C <- read.csv(file.path(DIR,sprintf("PBN_%s_coldata.csv",tag)), check.names=FALSE)
  C <- C[match(colnames(M), C$sample), ]
  if (naive) C <- C[C$group != "MS1_nat", ]
  Ml <- log2(M + 1)                       # log2 for limma-trend
  acc <- list()
  for (ct in sort(unique(C$cell_type)))
    acc[[ct]] <- fit_ct(Ml, complete_pairs(C[C$cell_type==ct,]), ct, label)
  do.call(rbind, acc[!sapply(acc,is.null)])
}

res <- rbind(
  run_level("coarse","meanCP10K","S2_meanCP10K_coarse8",  FALSE),
  run_level("coarse","meanCP10K","S2_meanCP10K_coarse8_naive", TRUE),
  run_level("fine",  "meanCP10K","S2_meanCP10K_fine25",   FALSE),
  run_level("donor", "meanCP10K","S2_meanCP10K_wholePBMC",FALSE),
  run_level("coarse","sumCP10K", "S3_sumCP10K_coarse8",   FALSE),
  run_level("fine",  "sumCP10K", "S3_sumCP10K_fine25",    FALSE),
  run_level("donor", "sumCP10K", "S3_sumCP10K_wholePBMC", FALSE))
res$FDR_global <- NA
for (k in unique(res$scheme)) { i <- res$scheme==k; res$FDR_global[i] <- p.adjust(res$PValue[i],"BH") }
write.csv(res, file.path(DIR,"pseudobulk_norm_compare.csv"), row.names=FALSE)

for (k in unique(res$scheme)) {
  r <- res[res$scheme==k,]; p <- r[r$gene %in% PANEL,]; p$FDR_panel <- p.adjust(p$PValue,"BH")
  cat(sprintf("\n======== %s ========\n", k))
  cat(sprintf("cell types=%d tests=%d | global FDR<0.05=%d | local FDR<0.05=%d | panel FDR<0.05=%d\n",
              length(unique(r$cell_type)), nrow(r), sum(r$FDR_global<0.05,na.rm=TRUE),
              sum(r$FDR_local<0.05,na.rm=TRUE), sum(p$FDR_panel<0.05,na.rm=TRUE)))
  q <- p[order(p$PValue),][1:min(8,nrow(p)),]
  for (i in seq_len(nrow(q))) with(q[i,],
    cat(sprintf("   %-8s %-9s %+8.4f  p=%.2e  local=%.3f global=%.3f panel=%.3f%s\n",
        gene, cell_type, logFC, PValue, FDR_local, FDR_global, FDR_panel,
        ifelse(FDR_panel<0.05|FDR_global<0.05,"  **",""))))
}
cat("\ndone -> pseudobulk_norm_compare.csv\n")
