#!/usr/bin/env Rscript
# Donor-level pseudobulk for the brain (Jaekel) and CSF/blood (Beltran) cohorts, processed the same
# way as the Kaufmann cohort so all three are comparable.
#   Jaekel : SUM of raw integer UMI counts per patient x cell type -> edgeR-QL and limma-voom
#            (4 MS / 5 HC patients; region blocks collapsed to the patient)
#   Beltran: only TPM-like normalised values deposited -> MEAN per twin x cell type, log2,
#            limma-trend (4 MS / 4 HC). NOT a count model - stated as such.
# Both cohorts are severely under-replicated; results are reported with that caveat.
suppressPackageStartupMessages({ library(edgeR); library(limma) })
DIR <- "__MS_GEO_ROOT__/Poster_v2/figures/pseudobulk_proper"
T1 <- c("ITGB2","LXN","CD79B","IKZF1","SH3BP4"); SUG <- "HLA-E"
AUX <- c("CASP6","CASP8","DGKQ","MX1","IFIT1","NUP210","RUNX3")
PANEL <- c(T1,SUG,AUX,"CTSZ","CHL1","THRB","ITGAL","IFI44L","RPAP2","SLAMF1","PCNP",
           "STAT3","TYK2","ICAM1","MOSPD3","FOXP3")
MIN_DON <- 3
acc <- list()

## ---- Jaekel: raw counts ----
M <- as.matrix(read.csv(file.path(DIR,"PBC_jakel_matrix.csv"), row.names=1, check.names=FALSE))
C <- read.csv(file.path(DIR,"PBC_jakel_coldata.csv"), check.names=FALSE)
C <- C[match(colnames(M), C$sample), ]
stopifnot(all(abs(M-round(M))<1e-9))
cat(sprintf("Jaekel: %d genes x %d pseudobulk | %d MS / %d HC patients | %d cell types\n",
  nrow(M),ncol(M),length(unique(C$donor[C$condition=="MS"])),
  length(unique(C$donor[C$condition=="HC"])),length(unique(C$cell_type))))
for (ct in sort(unique(C$cell_type))) {
  s <- C[C$cell_type==ct,]
  nms <- length(unique(s$donor[s$condition=="MS"])); nhc <- length(unique(s$donor[s$condition=="HC"]))
  if (nms < MIN_DON || nhc < MIN_DON) next
  grp <- factor(s$condition, levels=c("HC","MS")); design <- model.matrix(~grp)
  y <- DGEList(M[, s$sample, drop=FALSE], group=grp)
  keep <- filterByExpr(y, design=design)   # NO force-keep: unexpressed panel genes must not be tested
  y <- y[keep,,keep.lib.sizes=FALSE]; if (nrow(y)<50) next
  y <- normLibSizes(y,"TMM")
  y <- estimateDisp(y, design, robust=TRUE); fit <- glmQLFit(y, design, robust=TRUE)
  tt <- topTags(glmQLFTest(fit, coef=2), n=Inf, sort.by="none")$table
  acc[[paste("J_edgeR",ct)]] <- data.frame(cohort="Jaekel_brain", engine="edgeR_QL", cell_type=ct,
    gene=rownames(tt), logFC=tt$logFC, PValue=tt$PValue, FDR_local=p.adjust(tt$PValue,"BH"),
    n_ms=nms, n_hc=nhc, row.names=NULL)
  v <- voomLmFit(y, design, sample.weights=TRUE, plot=FALSE); v <- eBayes(v, robust=TRUE)
  t2 <- topTable(v, coef=2, number=Inf, sort.by="none")
  acc[[paste("J_voom",ct)]] <- data.frame(cohort="Jaekel_brain", engine="limma_voom", cell_type=ct,
    gene=rownames(t2), logFC=t2$logFC, PValue=t2$P.Value, FDR_local=p.adjust(t2$P.Value,"BH"),
    n_ms=nms, n_hc=nhc, row.names=NULL)
}

## ---- Beltran: normalised -> limma-trend ----
Mb <- as.matrix(read.csv(file.path(DIR,"PBN_beltran_meanNorm.csv"), row.names=1, check.names=FALSE))
Cb <- read.csv(file.path(DIR,"PBN_beltran_coldata.csv"), check.names=FALSE)
Cb <- Cb[match(colnames(Mb), Cb$sample), ]
cat(sprintf("Beltran: %d genes x %d pseudobulk | %d MS / %d HC donors | %d cell types\n",
  nrow(Mb),ncol(Mb),length(unique(Cb$donor[Cb$condition=="MS"])),
  length(unique(Cb$donor[Cb$condition=="HC"])),length(unique(Cb$cell_type))))
Lb <- log2(Mb + 1)
for (ct in sort(unique(Cb$cell_type))) {
  s <- Cb[Cb$cell_type==ct,]
  nms <- length(unique(s$donor[s$condition=="MS"])); nhc <- length(unique(s$donor[s$condition=="HC"]))
  if (nms < MIN_DON || nhc < MIN_DON) next
  grp <- factor(s$condition, levels=c("HC","MS")); design <- model.matrix(~grp)
  if (nrow(design)-ncol(design) < 2) next
  E <- Lb[, s$sample, drop=FALSE]
  keep <- rowMeans(E>0) >= 0.5
  y <- E[keep,,drop=FALSE]; if (nrow(y)<50) next
  fit <- eBayes(lmFit(y, design), trend=TRUE, robust=TRUE)
  tt <- topTable(fit, coef=2, number=Inf, sort.by="none")
  acc[[paste("B",ct)]] <- data.frame(cohort="Beltran_CSF", engine="limma_trend", cell_type=ct,
    gene=rownames(tt), logFC=tt$logFC, PValue=tt$P.Value, FDR_local=p.adjust(tt$P.Value,"BH"),
    n_ms=nms, n_hc=nhc, row.names=NULL)
}

res <- do.call(rbind, acc[!sapply(acc,is.null)]); rownames(res) <- NULL
res$FDR_global <- NA
for (k in unique(paste(res$cohort,res$engine))) {
  i <- paste(res$cohort,res$engine)==k; res$FDR_global[i] <- p.adjust(res$PValue[i],"BH")
}
write.csv(res, file.path(DIR,"pseudobulk_brain_csf.csv"), row.names=FALSE)

for (k in unique(paste(res$cohort,res$engine))) {
  r <- res[paste(res$cohort,res$engine)==k,]
  p <- r[r$gene %in% PANEL,]; p$FDR_panel <- p.adjust(p$PValue,"BH")
  cat(sprintf("\n=========== %s ===========\n", k))
  cat(sprintf("cell types=%d tests=%d | global FDR<0.05=%d | panel FDR<0.05=%d\n",
      length(unique(r$cell_type)), nrow(r), sum(r$FDR_global<0.05,na.rm=TRUE), sum(p$FDR_panel<0.05,na.rm=TRUE)))
  cat("  -- LXN and ITGB2, all cell types --\n")
  for (g in c("LXN","ITGB2")) { x <- p[p$gene==g,]; x <- x[order(x$PValue),]
    if (!nrow(x)) { cat(sprintf("    %-6s not quantified\n",g)); next }
    for (i in seq_len(min(4,nrow(x)))) with(x[i,],
      cat(sprintf("    %-6s %-14s logFC=%+8.3f p=%.3e local=%.3f panel=%.3f (%dMS/%dHC)%s\n",
          gene, cell_type, logFC, PValue, FDR_local, FDR_panel, n_ms, n_hc,
          ifelse(FDR_panel<0.05|FDR_local<0.05,"  **","")))) }
  q <- p[p$PValue<0.05,]; q <- q[order(q$PValue),]
  cat(sprintf("  -- all candidates with nominal p<0.05: %d --\n", nrow(q)))
  if (nrow(q)) for (i in seq_len(min(10,nrow(q)))) with(q[i,],
    cat(sprintf("    %-8s %-14s logFC=%+8.3f p=%.3e panel=%.3f%s\n", gene, cell_type, logFC, PValue,
        FDR_panel, ifelse(FDR_panel<0.05,"  **",""))))
}
cat("\ndone -> pseudobulk_brain_csf.csv\n")
