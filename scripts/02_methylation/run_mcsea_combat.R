#!/usr/bin/env Rscript
# mCSEA on the NEW ComBat-corrected, IDAT-reprocessed 8-dataset M-values
#   Input : Methylation_Data/AllMeth_ComBat_M.csv + AllMeth_ComBat_Metadata.csv
#   Method: mCSEA promoter + gene-body GSEA on MS-vs-HC ranked probes
#   Compare vs old combined-cohort mCSEA (06_mCSEA_promoter.tsv)
suppressPackageStartupMessages({ library(mCSEA); library(limma) })
msg <- function(...) cat(sprintf("[mcsea-combat] %s\n", paste0(...)))
MD <- "__MS_GEO_ROOT__/Methylation_Data"
INV1 <- c("HLA-E","ITGB2","LXN","CD79B","IKZF1","SH3BP4")
T2A  <- c("CASP6","CASP8","DGKQ","MX1","IFIT1","NUP210","RUNX3")

M <- as.matrix(read.csv(file.path(MD,"AllMeth_ComBat_M.csv"), check.names=FALSE, row.names=1))
meta <- read.csv(file.path(MD,"AllMeth_ComBat_Metadata.csv"), stringsAsFactors=FALSE)
meta <- meta[match(colnames(M), meta$sample_id), ]
pheno <- data.frame(condition=factor(meta$condition, levels=c("HC","MS")),
                    row.names=colnames(M))
msg(nrow(M)," probes × ", ncol(M)," samples | ", sum(pheno$condition=="MS")," MS / ",
    sum(pheno$condition=="HC")," HC")

msg("rankProbes (M-values, MS vs HC)...")
myRank <- rankProbes(M, pheno, refGroup="HC", typeInput="M", explanatory="condition")
msg("ranked ", length(myRank), " probes")

msg("mCSEA promoters + genes (450k, minCpGs=5)...")
res <- mCSEA(myRank, M, pheno, regionsTypes=c("promoters","genes"),
             platform="450k", minCpGs=5)
saveRDS(res, file.path(MD,"AllMeth_ComBat_mCSEA.rds"))
for(rt in c("promoters","genes")){
  d <- res[[rt]]
  if(is.null(d)) next
  d$gene <- rownames(d)
  write.csv(d, file.path(MD, paste0("AllMeth_ComBat_mCSEA_",rt,".csv")), row.names=FALSE)
  msg(rt,": ", nrow(d)," regions, ", sum(d$padj<0.05,na.rm=TRUE)," at padj<0.05")
}

prom <- res$promoters; gb <- res$genes
look <- function(g, d){ if(!is.null(d) && g %in% rownames(d)){
    r<-d[g,]; sprintf("NES=%+.2f padj=%.2e %s", r$NES, r$padj, ifelse(r$padj<0.05,"SIG","ns"))
  } else "absent" }
cat("\n=== mCSEA on NEW ComBat data — INV-Tier-1 ===\n")
cat(sprintf("%-8s %-34s %s\n","gene","PROMOTER","GENE-BODY"))
for(g in INV1) cat(sprintf("%-8s %-34s %s\n", g, look(g,prom), look(g,gb)))
cat("\n--- Tier-2 aux INV ---\n")
for(g in T2A) cat(sprintf("%-8s %-34s %s\n", g, look(g,prom), look(g,gb)))
msg("DONE")
