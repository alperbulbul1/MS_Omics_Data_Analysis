suppressPackageStartupMessages({library(data.table);library(limma)})
S <- "__MS_GEO_ROOT__/Methylation_Data"
E <- "__MS_GEO_ROOT__/Expression_Data"
CAND <- c("ITGB2","IKZF1","CD79B","LXN","SH3BP4","RUNX3","CASP6","CASP8","DGKQ",
          "MX1","IFIT1","NUP210","CTSZ","CHL1","ICAM1","HLA-E")
dt <- fread(file.path(E,"Corrected_Batch_Corrected_Expression.csv"))
g <- dt[[1]]; M <- as.matrix(dt[,-1]); rownames(M) <- g; rm(dt); gc()
md <- fread(file.path(S,"rna_sex_persample.csv"))
md <- md[sex %in% c("F","M")]
md <- md[sample_id %in% colnames(M)]
M <- M[, md$sample_id, drop=FALSE]
cat(sprintf("[rna-sex] %d genes x %d samples with sex\n", nrow(M), ncol(M)))
cond <- factor(md$condition, levels=c("HC","MS")); ds <- factor(md$dataset); sx <- factor(md$sex)
print(table(cond, sx)); print(chisq.test(table(cond, sx)))
keep <- rowSums(is.finite(M)) == ncol(M); M <- M[keep,]
cat(sprintf("[rna-sex] finite genes: %d\n", nrow(M)))
dA <- model.matrix(~ 0 + ds + cond); dB <- model.matrix(~ 0 + ds + cond + sx)
cA <- tail(colnames(dA),1); cB <- colnames(dB)[which(colnames(dB)==paste0("cond","MS"))]
fA <- eBayes(lmFit(M,dA)); fB <- eBayes(lmFit(M,dB))
tA <- topTable(fA,coef=cA,number=Inf,sort.by="none"); tB <- topTable(fB,coef=cB,number=Inf,sort.by="none")
cat(sprintf("[rna-sex] genes FDR<0.05  A(no sex): %d | B(+sex): %d | shared: %d\n",
    sum(tA$adj.P.Val<0.05), sum(tB$adj.P.Val<0.05), sum(tA$adj.P.Val<0.05 & tB$adj.P.Val<0.05)))
cat(sprintf("[rna-sex] logFC Pearson r = %.5f | p-value Spearman = %.5f\n",
    cor(tA$logFC,tB$logFC), cor(tA$P.Value,tB$P.Value,method="spearman")))
i <- rownames(tA) %in% CAND
o <- data.frame(gene=rownames(tA)[i], logFC_nosex=tA$logFC[i], FDR_nosex=tA$adj.P.Val[i],
                logFC_sex=tB$logFC[i], FDR_sex=tB$adj.P.Val[i])
o <- o[order(o$FDR_nosex),]
cat("\n=========== CANDIDATES ===========\n"); print(format(o,digits=3), row.names=FALSE)
cat(sprintf("\nsign changes: %d | sig without sex: %d | with sex: %d\n",
    sum(sign(o$logFC_nosex)!=sign(o$logFC_sex)), sum(o$FDR_nosex<0.05), sum(o$FDR_sex<0.05)))
write.csv(o, file.path(S,"rna_sex_sensitivity_candidates.csv"), row.names=FALSE)
