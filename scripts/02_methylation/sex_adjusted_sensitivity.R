#!/usr/bin/env Rscript
# Reviewer 3: is the methylation layer's MS-vs-HC signal sensitive to sex?
#
# Replicates the published AllMeth analysis exactly (run_all_methylation_combat.R): limma on the
# saved ComBat M-value matrix, gene level = best probe per gene by FDR. The ONLY difference between
# the two models compared here is the sex term; both are fitted on the SAME 448 samples for which
# sex is deposited, so the comparison isolates sex and is not confounded by a change of sample set.
suppressPackageStartupMessages({ library(data.table); library(limma)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19) })
msg <- function(...) cat(sprintf("[sex-sens] %s\n", paste0(...)))

MD <- "__MS_GEO_ROOT__/Methylation_Data"
SC <- "__MS_GEO_ROOT__/Methylation_Data"
CAND <- c("ITGB2","IKZF1","CD79B","LXN","SH3BP4","RUNX3","CASP6","CASP8","DGKQ",
          "MX1","IFIT1","NUP210","CTSZ","CHL1","ICAM1","HLA-E")

msg("reading ComBat M matrix (2.2 GB CSV)...")
dt <- fread(file.path(MD, "AllMeth_ComBat_M.csv"), showProgress = FALSE)
probes <- dt[[1]]
M <- as.matrix(dt[, -1]); rownames(M) <- probes; rm(dt); gc()
msg("M: ", nrow(M), " probes x ", ncol(M), " samples")

meta <- fread(file.path(SC, "meth_sex.csv"))
meta <- meta[match(colnames(M), meta$sample_id), ]
stopifnot(all(meta$sample_id == colnames(M)))

keep <- meta$sex %in% c("F", "M")
msg("samples with deposited sex: ", sum(keep), " / ", length(keep))
Mk <- M[, keep, drop = FALSE]; mk <- meta[keep, ]
rm(M); gc()

cond <- factor(mk$condition, levels = c("HC", "MS"))
sex  <- factor(mk$sex, levels = c("F", "M"))
msg("condition: ", sum(cond == "MS"), " MS / ", sum(cond == "HC"), " HC")
print(table(cond, sex))
cat("\nchi-square test of sex x condition independence:\n")
print(chisq.test(table(cond, sex)))

ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# gene-level table exactly as the published pipeline builds it: best probe per gene by FDR
by_gene <- function(fit, coef) {
  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  tt$Probe <- rownames(tt)
  g <- ann[tt$Probe, "UCSC_RefGene_Name"]
  ex <- do.call(rbind, lapply(seq_len(nrow(tt)), function(i) {
    gs <- unique(strsplit(g[i], ";")[[1]]); gs <- gs[nzchar(gs)]
    if (!length(gs)) return(NULL)
    data.frame(gene = gs, logFC = tt$logFC[i], P = tt$P.Value[i],
               FDR = tt$adj.P.Val[i], Probe = tt$Probe[i], stringsAsFactors = FALSE)
  }))
  ex <- ex[order(ex$FDR), ]
  ex[!duplicated(ex$gene), ]
}

msg("model A: ~ condition        (same 448 samples)")
fitA <- eBayes(lmFit(Mk, model.matrix(~ cond)))
A <- by_gene(fitA, 2)

msg("model B: ~ condition + sex  (same 448 samples)")
fitB <- eBayes(lmFit(Mk, model.matrix(~ cond + sex)))
B <- by_gene(fitB, 2)

ttA <- topTable(fitA, coef = 2, number = Inf, sort.by = "none")
ttB <- topTable(fitB, coef = 2, number = Inf, sort.by = "none")
msg("genome-wide probes FDR<0.05 -- A: ", sum(ttA$adj.P.Val < 0.05),
    " | B: ", sum(ttB$adj.P.Val < 0.05))
msg("Pearson r of condition logFC across all ", nrow(ttA), " probes: ",
    sprintf("%.5f", cor(ttA$logFC, ttB$logFC)))
msg("Spearman r of condition p-values: ",
    sprintf("%.5f", cor(ttA$P.Value, ttB$P.Value, method = "spearman")))

cmp <- merge(A[, c("gene","logFC","FDR","Probe")],
             B[, c("gene","logFC","FDR","Probe")],
             by = "gene", suffixes = c("_nosex", "_sex"))
cat("\n================ CANDIDATE GENES ================\n")
out <- cmp[cmp$gene %in% CAND, ]
out <- out[order(out$FDR_nosex), ]
print(format(out, digits = 3), row.names = FALSE)

cat("\nsign changes among candidates: ",
    sum(sign(out$logFC_nosex) != sign(out$logFC_sex)), "\n")
cat("candidates significant (FDR<0.05) without sex: ", sum(out$FDR_nosex < 0.05),
    " | with sex: ", sum(out$FDR_sex < 0.05), "\n")

cat("\ngenome-wide gene-level concordance:\n")
cat("  genes FDR<0.05 without sex: ", sum(cmp$FDR_nosex < 0.05), "\n")
cat("  genes FDR<0.05 with sex:    ", sum(cmp$FDR_sex < 0.05), "\n")
cat("  overlap:                    ",
    sum(cmp$FDR_nosex < 0.05 & cmp$FDR_sex < 0.05), "\n")

write.csv(cmp, file.path(SC, "sex_sensitivity_bygene.csv"), row.names = FALSE)
write.csv(out, file.path(SC, "sex_sensitivity_candidates.csv"), row.names = FALSE)
msg("done")
