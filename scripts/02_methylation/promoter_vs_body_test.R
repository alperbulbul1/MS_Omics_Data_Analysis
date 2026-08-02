#!/usr/bin/env Rscript
# Reviewer point 2, direct compartment test.
#
# mCSEA is a region-level ENRICHMENT test and needs a minimum number of CpGs per region, so genes
# with few promoter probes (e.g. CD79B, 3 probes) cannot be tested at all. That is a limitation of
# the test, not evidence about the gene. Here every gene in the discovery pool is instead given an
# explicit promoter-only and gene-body-only statistic by classifying each probe with the Illumina
# annotation (UCSC_RefGene_Group) and combining probe p-values within each compartment:
#     promoter  = TSS1500, TSS200, 5'UTR, 1stExon
#     gene body = Body, 3'UTR
# Combination is Stouffer's method on the per-probe p-values, signed by the probe logFC, matching
# the gene-level aggregation already used in the pipeline; the mean probe logFC is reported per
# compartment. This assigns a compartment to all 82 genes regardless of probe count.
suppressPackageStartupMessages({
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19); library(minfi)
})
RES <- "__MS_GEO_ROOT__/Methylation/results"
OUT <- "__MS_GEO_ROOT__/Poster_v2/figures"

ann <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
ann <- ann[, c("Name","UCSC_RefGene_Name","UCSC_RefGene_Group")]
# explode the semicolon-separated gene/group pairs
sp_g <- strsplit(ann$UCSC_RefGene_Name, ";", fixed=TRUE)
sp_r <- strsplit(ann$UCSC_RefGene_Group, ";", fixed=TRUE)
n    <- lengths(sp_g)
keep <- n > 0
map  <- data.frame(Probe = rep(ann$Name[keep], n[keep]),
                   Gene  = unlist(sp_g[keep]),
                   Group = unlist(sp_r[keep]), stringsAsFactors = FALSE)
PROM <- c("TSS1500","TSS200","5'UTR","1stExon"); BODY <- c("Body","3'UTR")
map$compartment <- ifelse(map$Group %in% PROM, "promoter",
                   ifelse(map$Group %in% BODY, "body", NA))
map <- unique(map[!is.na(map$compartment), c("Probe","Gene","compartment")])
cat(sprintf("annotation: %d probe-gene-compartment records | %d genes\n",
            nrow(map), length(unique(map$Gene))))

pool <- read.delim(file.path(RES,"INVERSE_CONCORDANT_by_gene.tsv"))$gene
cat(sprintf("discovery pool: %d genes\n", length(pool)))

stouffer <- function(p, s) {                      # signed Stouffer over probes
  p <- pmin(pmax(p, 1e-300), 1-1e-16)
  z <- qnorm(p/2, lower.tail=FALSE) * sign(s)
  zc <- sum(z)/sqrt(length(z))
  list(z=zc, p=2*pnorm(abs(zc), lower.tail=FALSE))
}

strata <- list("T cells"="01_tcells_meth_DMP.tsv", "Whole blood DMF"="02_wb_dmf_meth_DMP.tsv",
               "Whole blood OCR"="03_wb_ocrelizumab_meth_DMP.tsv",
               "T cells remission"="04_tcells_remission_meth_DMP.tsv",
               "Combined"="05_combined_meth_DMP.tsv")
out <- list()
for (sn in names(strata)) {
  d <- read.delim(file.path(RES, strata[[sn]]))
  rownames(d) <- d$Probe
  mm <- map[map$Gene %in% pool & map$Probe %in% d$Probe, ]
  for (g in unique(mm$Gene)) {
    sub <- mm[mm$Gene==g, ]
    row <- list(stratum=sn, gene=g)
    for (cmp in c("promoter","body")) {
      pr <- sub$Probe[sub$compartment==cmp]
      if (!length(pr)) { row[[paste0(cmp,"_n")]] <- 0L; row[[paste0(cmp,"_logFC")]] <- NA
                         row[[paste0(cmp,"_p")]] <- NA; next }
      dd <- d[pr, , drop=FALSE]
      st <- stouffer(dd$P.Value, dd$logFC)
      row[[paste0(cmp,"_n")]]     <- length(pr)
      row[[paste0(cmp,"_logFC")]] <- mean(dd$logFC)
      row[[paste0(cmp,"_p")]]     <- st$p
    }
    out[[length(out)+1]] <- as.data.frame(row, stringsAsFactors=FALSE)
  }
}
res <- do.call(rbind, out)
# BH within stratum x compartment
for (sn in unique(res$stratum)) {
  i <- res$stratum==sn
  res$promoter_fdr[i] <- p.adjust(res$promoter_p[i], "BH")
  res$body_fdr[i]     <- p.adjust(res$body_p[i], "BH")
}
write.csv(res, file.path(OUT,"promoter_vs_body_by_gene.csv"), row.names=FALSE)
cat(sprintf("\nwrote promoter_vs_body_by_gene.csv (%d gene x stratum rows)\n", nrow(res)))

T1 <- c("ITGB2","LXN","CD79B","IKZF1","SH3BP4")
cat("\n=================== TIER-1: PROMOTER vs GENE BODY ===================\n")
for (g in T1) {
  x <- res[res$gene==g, ]
  if (!nrow(x)) { cat(sprintf("\n%s: not annotated\n", g)); next }
  cat(sprintf("\n### %s   (promoter probes n=%d, body probes n=%d)\n",
              g, max(x$promoter_n), max(x$body_n)))
  for (i in seq_len(nrow(x))) with(x[i,], {
    pl <- if (promoter_n>0) sprintf("logFC=%+7.4f p=%.3g FDR=%.3g", promoter_logFC, promoter_p, promoter_fdr) else "no promoter probes"
    bl <- if (body_n>0)     sprintf("logFC=%+7.4f p=%.3g FDR=%.3g", body_logFC, body_p, body_fdr)             else "no body probes"
    cat(sprintf("   %-18s promoter[%d]: %-42s | body[%d]: %s\n", stratum, promoter_n, pl, body_n, bl))
  })
}
cat("\n=================== POOL-WIDE SUMMARY (per stratum) ===================\n")
for (sn in unique(res$stratum)) {
  x <- res[res$stratum==sn, ]
  cat(sprintf("%-18s genes=%3d | promoter FDR<0.05: %3d | body FDR<0.05: %3d | promoter-only signal: %3d\n",
      sn, nrow(x), sum(x$promoter_fdr<0.05,na.rm=TRUE), sum(x$body_fdr<0.05,na.rm=TRUE),
      sum(x$promoter_fdr<0.05 & (is.na(x$body_fdr) | x$body_fdr>=0.05), na.rm=TRUE)))
}
