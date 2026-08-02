#!/usr/bin/env Rscript
# =============================================================================
# run_all_methylation_combat.R — ALL methylation datasets, batch-corrected
#   M-value sources:
#     • 6 IDAT datasets (minfi preprocessIllumina): IDAT6_Preprocessed/Combined_..._Strict_M.csv
#     • GSE106648 (beta-only, no IDAT deposited): series-matrix beta → M
#     • GSE40360  (beta-only brain WM, no IDAT):  series-matrix beta → M
#   Merge common 450K probes → ComBat (batch=dataset, mod=~condition) → limma → Tier-1
# =============================================================================
suppressPackageStartupMessages({ library(sva); library(limma)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19) })
msg <- function(...) cat(sprintf("[all-meth] %s\n", paste0(...)))
MD  <- "__MS_GEO_ROOT__/Methylation_Data"
INV1<- c("HLA-E","ITGB2","LXN","CD79B","IKZF1","SH3BP4")

rdM <- function(p){ x<-read.csv(p, check.names=FALSE, row.names=1); as.matrix(x) }
M_idat <- rdM(file.path(MD,"IDAT6_Preprocessed","Combined_Methylation_Strict_M.csv"))
M_648  <- rdM(file.path(MD,"GSE106648_betaonly_M.csv"))
M_403  <- rdM(file.path(MD,"GSE40360_betaonly_M.csv"))
msg("IDAT6 M: ",nrow(M_idat),"×",ncol(M_idat)," | GSE106648: ",nrow(M_648),"×",ncol(M_648),
    " | GSE40360: ",nrow(M_403),"×",ncol(M_403))

common <- Reduce(intersect, list(rownames(M_idat),rownames(M_648),rownames(M_403)))
msg("common 450K probes across all sources: ", length(common))
M <- cbind(M_idat[common,,drop=FALSE], M_648[common,,drop=FALSE], M_403[common,,drop=FALSE])
M <- M[is.finite(rowSums(M)), , drop=FALSE]
msg("merged M: ", nrow(M)," probes × ", ncol(M)," samples")

# metadata
om <- read.csv(file.path(MD,"Combined_Methylation_Metadata.csv"), stringsAsFactors=FALSE)
om$base <- sub("__.*$","",om$dataset)
om <- om[match(colnames(M), om$sample_id), ]
stopifnot(all(!is.na(om$condition)))
batch <- factor(om$base); cond <- factor(om$condition, levels=c("HC","MS"))
msg("datasets(batches): ", paste(table(batch),names(table(batch)),collapse="  "))
msg("condition: ", sum(cond=="MS")," MS / ", sum(cond=="HC")," HC across ", nlevels(batch)," datasets")

# ── ComBat (preserve condition) ─────────────────────────────────────────────
mod <- model.matrix(~ cond)
msg("ComBat...")
Mc <- ComBat(dat=M, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
write.csv(data.frame(Probe=rownames(Mc), Mc, check.names=FALSE),
          file.path(MD,"AllMeth_ComBat_M.csv"), row.names=FALSE)
write.csv(data.frame(sample_id=colnames(Mc), dataset=as.character(batch),
                     condition=as.character(cond)),
          file.path(MD,"AllMeth_ComBat_Metadata.csv"), row.names=FALSE)
msg("saved AllMeth_ComBat_M.csv (",nrow(Mc),"×",ncol(Mc),")")

# ── PCA before/after ────────────────────────────────────────────────────────
rv<-function(x) apply(x,1,var)
pca<-function(mat){ v<-rv(mat); top<-order(v,decreasing=TRUE)[1:min(15000,length(v))]
  pr<-prcomp(t(mat[top,])); list(x=pr$x[,1],y=pr$x[,2],pv=round(100*pr$sdev[1:2]^2/sum(pr$sdev^2),1)) }
png(file.path(MD,"AllMeth_batch_PCA_before_after.png"),width=1700,height=1700,res=150)
par(mfrow=c(2,2),mar=c(4,4,3,1))
bc<-rainbow(nlevels(batch)); names(bc)<-levels(batch); cc<-c(HC="#1976D2",MS="#D32F2F")
for(st in list(list(M,"BEFORE ComBat"),list(Mc,"AFTER ComBat"))){
  p<-pca(st[[1]])
  plot(p$x,p$y,col=bc[as.character(batch)],pch=19,cex=0.8,xlab=paste0("PC1 (",p$pv[1],"%)"),
       ylab=paste0("PC2 (",p$pv[2],"%)"),main=paste0(st[[2]]," — by DATASET"))
  legend("topright",names(bc),col=bc,pch=19,cex=0.6,bty="n")
  plot(p$x,p$y,col=cc[as.character(cond)],pch=19,cex=0.8,xlab=paste0("PC1 (",p$pv[1],"%)"),
       ylab=paste0("PC2 (",p$pv[2],"%)"),main=paste0(st[[2]]," — by CONDITION"))
  legend("topright",names(cc),col=cc,pch=19,cex=0.8,bty="n")
}
dev.off(); msg("saved PCA QC")

# ── limma MS-vs-HC on ComBat M + Tier-1 ─────────────────────────────────────
fit<-eBayes(lmFit(Mc, model.matrix(~cond)))
tt<-topTable(fit,coef=2,number=Inf,sort.by="none"); tt$Probe<-rownames(tt)
ann<-getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
tt$gene<-ann[tt$Probe,"UCSC_RefGene_Name"]
ex<-do.call(rbind,lapply(seq_len(nrow(tt)),function(i){
  gs<-unique(strsplit(tt$gene[i],";")[[1]]); gs<-gs[nzchar(gs)]; if(!length(gs))return(NULL)
  data.frame(gene=gs,logFC=tt$logFC[i],P=tt$P.Value[i],FDR=tt$adj.P.Val[i],Probe=tt$Probe[i],stringsAsFactors=FALSE)}))
best<-ex[order(ex$FDR),]; best<-best[!duplicated(best$gene),]
write.csv(best,file.path(MD,"AllMeth_ComBat_limma_DMP_byGene.csv"),row.names=FALSE)
msg("genome-wide sig probes FDR<0.05: ",sum(tt$adj.P.Val<0.05),
    " | genes(best-probe) FDR<0.05: ",sum(best$FDR<0.05))
cat("\n=== INV-Tier-1 (ALL datasets, IDAT+beta merged, ComBat) ===\n")
for(g in INV1){ r<-best[best$gene==g,]
  if(nrow(r)) cat(sprintf("  %-8s logFC=%+.3f FDR=%.2e %s (%s)\n",g,r$logFC[1],r$FDR[1],
                          ifelse(r$FDR[1]<0.05,"SIG","ns"),r$Probe[1]))
  else cat(sprintf("  %-8s absent\n",g)) }
msg("DONE")
