suppressPackageStartupMessages({library(data.table);library(IlluminaHumanMethylation450kanno.ilmn12.hg19)})
b <- fread("__MS_GEO_ROOT__/Methylation_Data/New_Datasets/GSE88824_beta_raw.csv")
p <- b[[1]]; B <- as.matrix(b[,-1]); rownames(B) <- p
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
xp <- rownames(ann)[ann$chr=="chrX"]; yp <- rownames(ann)[ann$chr=="chrY"]
xp <- intersect(xp, rownames(B)); yp <- intersect(yp, rownames(B))
cat(sprintf("chrX probes %d | chrY probes %d | samples %d\n", length(xp), length(yp), ncol(B)))
Ymed <- apply(B[yp,,drop=FALSE], 2, median, na.rm=TRUE)
Xint <- apply(B[xp,,drop=FALSE], 2, function(v) mean(v>0.3 & v<0.7, na.rm=TRUE))
d <- data.frame(sample=colnames(B), chrY_median=round(Ymed,3), chrX_intermediate=round(Xint,3))
d <- d[order(d$chrY_median),]
print(d, row.names=FALSE)
# females: low chrY signal, HIGH chrX intermediate fraction
d$sex <- ifelse(d$chrX_intermediate > 0.35, "F", "M")
cat(sprintf("\ncall: %d F / %d M\n", sum(d$sex=="F"), sum(d$sex=="M")))
write.csv(d, "__MS_GEO_ROOT__/Methylation_Data/gse88824_sex.csv", row.names=FALSE)
