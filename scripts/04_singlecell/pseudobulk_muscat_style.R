#!/usr/bin/env Rscript
# Literature-standard pseudobulk differential-state analysis (Crowell/muscat; Squair 2021),
# following the canonical decision scheme:
#   1. cells clustered/annotated in a common space  -> the study's own published annotation
#   2. SUM raw integer UMI counts per sample x cell type
#   3. separate count matrix per cell type
#   4. filter low-cell sample-clusters (min_cells=10) and low-expressed genes (filterByExpr)
#   5. sample-level model with edgeR quasi-likelihood F-test (+ limma-voom as second engine)
#   6. sample-level covariates in the design (donor pair; sex where not absorbed by pairing)
#   7. report logFC + LOCAL FDR (within cell type) and GLOBAL FDR (all gene x cell type)
# Differential abundance is analysed separately (see pseudobulk_DA.R) - it is a different question.
suppressPackageStartupMessages({ library(edgeR); library(limma) })
DIR <- "__MS_GEO_ROOT__/Poster_v2/figures/pseudobulk_proper"
T1    <- c("ITGB2","LXN","CD79B","IKZF1","SH3BP4"); SUG <- "HLA-E"
AUX   <- c("CASP6","CASP8","DGKQ","MX1","IFIT1","NUP210","RUNX3")
PANEL <- c(T1,SUG,AUX,"CTSZ","CHL1","THRB","ITGAL","IFI44L","RPAP2","SLAMF1","PCNP",
           "STAT3","TYK2","ICAM1","MOSPD3","FOXP3")

complete_pairs <- function(s) {
  fp <- names(which(tapply(s$condition, s$batch_pair, function(z) all(c("MS","HC") %in% z))))
  s[s$batch_pair %in% fp, ]
}

# ---- one cell type, one engine ----
fit_ct <- function(cnt, s, ct, engine, label) {
  grp <- factor(s$condition, levels=c("HC","MS"))
  if (sum(grp=="MS") < 3 || sum(grp=="HC") < 3) return(NULL)
  pr <- factor(as.character(s$batch_pair)); if (nlevels(pr) < 3) return(NULL)
  design <- model.matrix(~ pr + grp)                 # paired: pair blocks donor-level baseline
  if (nrow(design) - ncol(design) < 2) return(NULL)
  cc <- ncol(design)
  y <- DGEList(counts=cnt[, s$sample, drop=FALSE], group=grp)
  keep <- filterByExpr(y, design=design); keep[rownames(y) %in% PANEL] <- TRUE
  y <- y[keep, , keep.lib.sizes=FALSE]; if (nrow(y) < 50) return(NULL)
  y <- normLibSizes(y, method="TMM")
  if (engine == "edgeR_QL") {
    y <- estimateDisp(y, design, robust=TRUE)
    fit <- glmQLFit(y, design, robust=TRUE)
    tt  <- topTags(glmQLFTest(fit, coef=cc), n=Inf, sort.by="none")$table
    out <- data.frame(gene=rownames(tt), logFC=tt$logFC, PValue=tt$PValue)
  } else {
    v   <- voomLmFit(y, design, sample.weights=TRUE, plot=FALSE)
    fit <- eBayes(v, robust=TRUE)
    tt  <- topTable(fit, coef=cc, number=Inf, sort.by="none")
    out <- data.frame(gene=rownames(tt), logFC=tt$logFC, PValue=tt$P.Value)
  }
  out$FDR_local <- p.adjust(out$PValue, "BH")          # local: within this cell type
  cbind(analysis=label, engine=engine, cell_type=ct,
        n_ms=sum(grp=="MS"), n_hc=sum(grp=="HC"), n_pairs=nlevels(pr),
        n_cells=sum(s$n_cells), out, row.names=NULL)
}

run_level <- function(tag, label, naive_only) {
  M <- as.matrix(read.csv(file.path(DIR,sprintf("PBC_%s_matrix.csv",tag)), row.names=1, check.names=FALSE))
  C <- read.csv(file.path(DIR,sprintf("PBC_%s_coldata.csv",tag)), check.names=FALSE)
  C <- C[match(colnames(M), C$sample), ]
  stopifnot(all(abs(M-round(M)) < 1e-9))               # assert integer counts
  if (naive_only) C <- C[C$group != "MS1_nat", ]
  acc <- list()
  for (ct in sort(unique(C$cell_type))) {
    s <- complete_pairs(C[C$cell_type==ct, ])
    for (eng in c("edgeR_QL","limma_voom"))
      acc[[paste(ct,eng)]] <- fit_ct(M, s, ct, eng, label)
  }
  do.call(rbind, acc[!sapply(acc,is.null)])
}

cat("running literature-standard pseudobulk (summed raw counts)...\n")
res <- rbind(
  run_level("coarse","A_coarse8_allPairs",   FALSE),
  run_level("coarse","B_coarse8_naiveOnly",  TRUE),
  run_level("fine",  "C_fine25_allPairs",    FALSE),
  run_level("fine",  "D_fine25_naiveOnly",   TRUE),
  run_level("donor", "E_wholePBMC_allPairs", FALSE),
  run_level("donor", "F_wholePBMC_naiveOnly",TRUE))
# GLOBAL FDR: across all gene x cell-type tests within each analysis x engine
res$FDR_global <- NA
for (k in unique(paste(res$analysis,res$engine))) {
  i <- paste(res$analysis,res$engine)==k
  res$FDR_global[i] <- p.adjust(res$PValue[i], "BH")
}
write.csv(res, file.path(DIR,"pseudobulk_muscat_style.csv"), row.names=FALSE)

for (k in unique(paste(res$analysis,res$engine))) {
  r <- res[paste(res$analysis,res$engine)==k, ]
  p <- r[r$gene %in% PANEL, ]; p$FDR_panel <- p.adjust(p$PValue,"BH")
  cat(sprintf("\n============ %s ============\n", k))
  cat(sprintf("cell types=%d | tests=%d | GLOBAL FDR<0.05=%d | any-local FDR<0.05=%d\n",
              length(unique(r$cell_type)), nrow(r), sum(r$FDR_global<0.05,na.rm=TRUE),
              sum(r$FDR_local<0.05,na.rm=TRUE)))
  cat(sprintf("candidate panel: %d tests | panel-FDR<0.05=%d\n", nrow(p), sum(p$FDR_panel<0.05,na.rm=TRUE)))
  q <- p[p$PValue<0.05, ]; q <- q[order(q$PValue), ]
  if (nrow(q)) { cat("  candidates with nominal p<0.05 (logFC | p | local | global | panel):\n")
    for (i in seq_len(min(12,nrow(q)))) with(q[i,],
      cat(sprintf("    %-8s %-9s %+7.3f  %.2e  %6.3f  %6.3f  %6.3f%s\n", gene, cell_type, logFC,
          PValue, FDR_local, FDR_global, FDR_panel,
          ifelse(FDR_local<0.05 | FDR_global<0.05 | FDR_panel<0.05, "  **", "")))) }
  else cat("  (no candidate with nominal p<0.05)\n")
}
cat("\ndone -> pseudobulk_muscat_style.csv\n")
