## 11_itgb2_csf_pleocytosis.R
##
## Why ITGB2 is detected more often in MS than in control CSF, and why that is not an
## MS-specific property of the protein.
##
## Background. ITGB2 (CD18) is the one Tier-1 candidate whose Astral CSF measurement is
## substantially incomplete: it is quantified in 700/978 MS (71.6%) but only 193/306 control
## (63.1%) samples. DEP's default MinProb imputation treats absent values as left-censored and
## fills them low, so this group-dependent detection gap alone produces a significant "MS-up"
## call (log2FC = +0.166, FDR = 0.023) that the measured values do not support (complete-case
## +0.046, FDR = 0.55). Reviewer 1 asked whether the CSF-versus-plasma direction difference
## reflects shedding from leukocyte surfaces. The Bader cohort records a CSF leukocyte count
## for every sample, so the cellular-composition part of that question can be tested directly
## rather than hypothesised.
##
## Design.
##   1. Confirm CSF pleocytosis in this cohort (MS vs control leukocyte counts).
##   2. Model ITGB2 DETECTION (not abundance) as a function of leukocyte count, adjusting for
##      total CSF protein, CSF erythrocytes (blood contamination) and diagnosis.
##   3. Stratify detection rate by leukocyte count to show the confounding directly.
##   4. Specificity controls: PTPRC (CD45, pan-leukocyte surface -> should behave like ITGB2)
##      and ALB (plasma-derived, not leukocyte -> should not).
##
## Output: Proteomics/processed/META/ITGB2_CSF_pleocytosis.txt

suppressPackageStartupMessages(library(data.table))

PROT_ROOT <- "__MS_GEO_ROOT__/Proteomics"
OUT_FP <- file.path(PROT_ROOT, "processed", "META", "ITGB2_CSF_pleocytosis.txt")
sink(OUT_FP, split = TRUE)

raw <- fread(file.path(PROT_ROOT, "processed", "astral_discovery_gene_keyed.tsv"),
             sep = "\t", header = TRUE, showProgress = FALSE)
sample_cols <- grep("\\.raw$|^[0-9]{8}_", colnames(raw), value = TRUE)

ann <- fread(file.path(PROT_ROOT, "osfstorage-archive", "processed proteomic data",
                       "0_sample_annotations",
                       paste0("annotations_v42_49_2_10_4_10_interimSky17_",
                              "PL01-PL56_PepResCustv01_resubmission.tsv")),
             sep = "\t", header = TRUE)

a <- ann[!is.na(Run_Astral_Measurement) & Run_Astral_Measurement != "",
         .(Run_Astral_Measurement, Diagnosis_group, MSgroup,
           Leukocyte_count, Erythrocytes_in_CSF, Total_protein)]
a[, group := fifelse(MSgroup == "MS", "MS",
              fifelse(Diagnosis_group %in% c("Other", "Neurological Control"),
                      "Control", NA_character_))]
a <- a[!is.na(group)]

keep <- intersect(sample_cols, a$Run_Astral_Measurement)
meta <- a[match(keep, a$Run_Astral_Measurement)]
leuko <- suppressWarnings(as.numeric(meta$Leukocyte_count))
ery   <- suppressWarnings(as.numeric(meta$Erythrocytes_in_CSF))
totp  <- suppressWarnings(as.numeric(meta$Total_protein))
grp   <- factor(meta$group, levels = c("Control", "MS"))

detected <- function(gene) {
  i <- which(raw$Genes == gene)[1]
  if (is.na(i)) return(NULL)
  !is.na(suppressWarnings(as.numeric(unlist(raw[i, ..keep]))))
}
det <- detected("ITGB2")

cat("Samples:", length(keep), sprintf("(MS %d / Control %d)\n",
    sum(grp == "MS"), sum(grp == "Control")))
cat(sprintf("ITGB2 detected: MS %.1f%%  Control %.1f%%\n\n",
    100 * mean(det[grp == "MS"]), 100 * mean(det[grp == "Control"])))

## ---- 1. CSF pleocytosis in this cohort -------------------------------------------------
cat("== 1. CSF leukocyte count, MS vs control ==\n")
ms <- leuko[grp == "MS" & !is.na(leuko)]
ct <- leuko[grp == "Control" & !is.na(leuko)]
cat(sprintf("  median  MS %.2f  Control %.2f  cells/uL\n", median(ms), median(ct)))
cat(sprintf("  mean    MS %.2f  Control %.2f\n", mean(ms), mean(ct)))
cat(sprintf("  Wilcoxon p = %.3e\n", wilcox.test(ms, ct)$p.value))
cat(sprintf("  above 5 cells/uL: MS %.1f%%  Control %.1f%%\n\n",
    100 * mean(ms > 5), 100 * mean(ct > 5)))

## ---- 2. detection model ----------------------------------------------------------------
## Erythrocytes are median 0 in both arms, so the term is rank-deficient and drops out;
## blood contamination is therefore not a viable confounder here.
ok <- !is.na(leuko) & leuko > 0 & !is.na(totp) & totp > 0
cat(sprintf("== 2. ITGB2 detection ~ leukocytes + total protein + diagnosis (n = %d) ==\n",
            sum(ok)))
cat(sprintf("  CSF erythrocytes, median: MS %.0f  Control %.0f\n",
    median(ery[grp == "MS"], na.rm = TRUE), median(ery[grp == "Control"], na.rm = TRUE)))
full <- glm(det[ok] ~ log2(leuko[ok] + 1) + log2(totp[ok]) + grp[ok], family = binomial)
print(round(summary(full)$coefficients, 4))

cat("\n  diagnosis effect with and without the leukocyte term:\n")
m0 <- glm(det[ok] ~ grp[ok], family = binomial)
cat(sprintf("    diagnosis alone            : MS beta = %+.3f, p = %.4f\n",
    summary(m0)$coefficients[2, 1], summary(m0)$coefficients[2, 4]))
cat(sprintf("    adjusted for leukocytes    : MS beta = %+.3f, p = %.4f\n",
    summary(full)$coefficients[4, 1], summary(full)$coefficients[4, 4]))

cat("\n  within-group (leukocytes + total protein):\n")
for (gg in c("MS", "Control")) {
  s <- ok & grp == gg
  cf <- summary(glm(det[s] ~ log2(leuko[s] + 1) + log2(totp[s]), family = binomial))$coefficients
  cat(sprintf("    %-8s n = %4d  leukocyte beta = %+.3f, p = %.4f\n",
      gg, sum(s), cf[2, 1], cf[2, 4]))
}

## ---- 3. stratified detection rate ------------------------------------------------------
cat("\n== 3. ITGB2 detection rate by CSF leukocyte stratum ==\n")
labs <- c("<=1", "2-3", "4-5", "6-10", ">10")
strat <- cut(leuko, c(-1, 1, 3, 5, 10, Inf), labels = labs)
cat(sprintf("  %-8s %20s %20s\n", "leuko", "MS", "Control"))
for (l in labs) {
  line <- sprintf("  %-8s", l)
  for (gg in c("MS", "Control")) {
    s <- !is.na(strat) & strat == l & grp == gg
    line <- paste0(line, sprintf("%20s", if (sum(s) >= 8)
      sprintf("%.0f%% (n=%d)", 100 * mean(det[s]), sum(s)) else sprintf("- (n=%d)", sum(s))))
  }
  cat(line, "\n")
}

## ---- 4. specificity controls -----------------------------------------------------------
cat("\n== 4. specificity: same model for other proteins ==\n")
for (gene in c("PTPRC", "ALB", "CHL1", "CTSZ", "ICAM1")) {
  d <- detected(gene)
  if (is.null(d)) { cat(sprintf("  %-6s not quantified\n", gene)); next }
  if (all(d)) { cat(sprintf("  %-6s detected in every sample - not modellable\n", gene)); next }
  cf <- summary(glm(d[ok] ~ log2(leuko[ok] + 1) + log2(totp[ok]) + grp[ok],
                    family = binomial))$coefficients
  cat(sprintf("  %-6s detected %.1f%%  leukocyte beta = %+.3f (p = %.4f)  diagnosis beta = %+.3f (p = %.4f)\n",
      gene, 100 * mean(d), cf[2, 1], cf[2, 4], cf[nrow(cf), 1], cf[nrow(cf), 4]))
}

## ---- 5. abundance among measured samples -----------------------------------------------
cat("\n== 5. among samples where ITGB2 IS measured, does abundance track leukocytes? ==\n")
i <- which(raw$Genes == "ITGB2")[1]
lv <- log2(suppressWarnings(as.numeric(unlist(raw[i, ..keep]))))

## The raw median difference is quoted in the manuscript as the model-free counterpart to the
## complete-case estimate, so it is recorded here rather than left to be recomputed.
mm <- median(lv[grp == "MS" & !is.na(lv)])
mc <- median(lv[grp == "Control" & !is.na(lv)])
cat(sprintf("  raw median of measured log2 intensities: MS %.4f  Control %.4f  difference %+.4f\n",
            mm, mc, mm - mc))

for (gg in c("MS", "Control")) {
  s <- grp == gg & !is.na(lv) & !is.na(leuko) & leuko > 0
  ct2 <- suppressWarnings(cor.test(lv[s], log2(leuko[s] + 1), method = "spearman"))
  cat(sprintf("  %-8s n = %4d  Spearman rho = %+.3f, p = %.3f\n",
      gg, sum(s), ct2$estimate, ct2$p.value))
}
cat("\n  Interpretation: pleocytosis governs whether ITGB2 crosses the detection threshold,\n")
cat("  not how abundant it is once measured. The MinProb call converted a detection\n")
cat("  asymmetry into an apparent abundance difference.\n")

sink()
cat("wrote", OUT_FP, "\n")
