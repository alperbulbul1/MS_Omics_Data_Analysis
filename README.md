# MS multi-omics re-analysis — analysis code

Code for *"Inverse-concordant DNA-methylation × transcription integration across four omic layers
prioritises ITGB2 and IKZF1 in multiple sclerosis"* (IJMS, under revision).

This repository holds the **analysis pipeline only**. Manuscript-authoring utilities (LaTeX/DOCX
patching, citation renumbering, figure relabelling) are deliberately excluded; 89 of the project's
529 R/Python files are published here.

---

## 1. Point it at your machine

Every script was written against absolute paths. They ship with two placeholders, which
`configure.sh` substitutes in place:

```bash
./configure.sh /path/to/MS_GEO_data /path/to/python
./configure.sh --check
```

`--check` fails if any placeholder remains **or** if an author-specific path survives, so a
misconfigured tree cannot be run by accident. Re-running with different values re-substitutes.

## 2. Environment

`env/r-packages.txt` and `env/requirements.txt` pin the versions that produced the reported
results, captured from the live session rather than transcribed from the paper.

**Read `env/OPTIONAL_AND_FALLBACKS.md` before assuming a package is required.** Eleven packages are
imported somewhere but were *not installed* when the results were produced, and the scripts took a
documented fallback path. The most consequential:

| package | what actually ran |
|---|---|
| `DEP` | not installed — the proteomic differential tests ran on the DEP-equivalent functions in `Proteomics/r_notebooks/helpers.R` (`filter_missval_R`, `vsn_with_fallback`, `moderated_t_safe`, `dep_completecase_de`) |
| `wateRmelon` | not installed — BMIQ normalisation of the beta-only methylation series fell back to quantile normalisation |
| `maxprobes` | not installed — cross-reactive probe removal used the built-in list |

## 3. Data

No data is redistributed. `docs/DATA.md` lists every accession with the layer it feeds. In short:
29 primary datasets — 14 bulk-transcriptomic GEO series, 8 Illumina 450K/EPIC methylation series
plus GSE173787 (WGBS), 3 single-cell series, PRIDE PXD064570 and PXD045058, MassIVE MSV000096790,
and the published UK Biobank-PPP Olink summary statistics.

## 4. Layout and run order

```
scripts/
  00_data/           acquisition, harmonisation, matrix assembly
  01_transcriptome/  per-stratum limma DE + pan-tissue combined model
  02_methylation/    minfi/ComBat preprocessing, DMP, mCSEA, gene-level aggregation
  03_proteomics/     CSF and brain differential abundance (complete-case)
  04_singlecell/     cohort processing (GSE118257 / GSE127969 / GSE144744) + donor-level pseudobulk
  05_integration/    STRING PPI + g:Profiler pathway analysis, data-source table, omics inventory
  06_figures/        one script per manuscript figure (figure1_workflow.py … figure7_string_network.py)
```

`docs/RUN_ORDER.md` gives the dependency order with the inputs and outputs of each step, and flags
the steps that are **not safe to re-run blindly** — several overwrite their own inputs.

## 5. Scope

This release contains **only** code on the path to the results the paper reports. Dataset-search
scripts from the project's exploratory phase, code for data that never entered the paper
(E-MTAB-69, CELLxGENE census queries), and QC-only plotting were removed alongside the
superseded analyses. Analyses
that the revision replaced are not included — notably the MinProb-imputed proteomic variants
(withdrawn because `ITGB2` is detected in 71.6 % of MS but only 63.1 % of control CSF samples,
so left-censored imputation manufactured a significant MS-up call the measured values do not
support), the dependence-aware methylation weighting experiments, and the legacy v1
acquisition track. The proteomic layer runner `00_run_all.R` invokes the complete-case
scripts (`*cc_*_completecase.R`), which are the reported analysis.

## 6. Known reproducibility limits

Stated plainly, because they are real and a reader will otherwise find them:

1. **Two upstream artefacts cannot be regenerated from the scripts.**
   `Expression_Data/Corrected_Metadata_ComBat.csv` was reduced from 552 to 472 rows after the
   matching expression matrix was written, and no script performs that reduction;
   `Global_Harmonized_Metadata.csv` was later overwritten down to 150 rows, so the 552-sample
   combined matrix cannot be rebuilt from it. Both files are inputs to the transcriptome layer.
   The analysed strata under `Stratified_Analyses/Expression/` are unaffected and self-consistent:
   they hold 462 samples across 13 series, and every reported RNA result derives from them.
2. **Two harmonisation steps are network-dependent.** `harmonize_microarray_v2.py` resolves
   probe → gene symbol maps live from GEO platform SOFT files and from mygene.info, so the gene
   universe can drift with the date of the run.
3. **Several scripts overwrite their own inputs** and are not idempotent — `RUN_ORDER.md` marks
   each one.
4. **One step streams a ~16 GB matrix** (`04_singlecell/`); budget memory accordingly.
5. **Four intermediate tables have no producing script anywhere**, so the steps that read them
   cannot be re-derived from this release. Each was traced and is listed rather than hidden:
   - `Poster_v2/figures/scrna_WILCOXON_v1.tsv` — the cell-level Wilcoxon scan feeding
     `06_figures/figure5_singlecell.py` and several Supplementary Table S2 values (it holds the
     quoted IKZF1 / T-cell FDR 9.4 × 10⁻⁷¹). Regenerated by an ad-hoc session after the Kaufmann
     UMAP repair; no script survived.
   - `Poster_v2/figures/scrna_PSEUDOBULK_COMPARISON.tsv` — feeds workbook sheet
     `08_CellLevel_vs_Donor`. `04_singlecell/run_pseudobulk_reanalysis.py` produces its 192-row
     parent table and is shipped, but the filter-and-rename step that reduces it to the 27 rows
     used was never saved.
   - `Poster_v2/figures/COMBINED_pantissue_proper_DEG.csv` — panel G of Figure 2. Note this is a
     different table from `07_pan_tissue_DE.tsv`, which `01_transcriptome/07_total_combined_de.R`
     does produce (10,744 vs 7,116 genes; 6,597 shared).
   - `Expression_Data/Corrected_Metadata_ComBat.csv` — see limitation 1 above.

## 7. Citation

If you use this code, please cite the paper. Data citations belong to the original depositors,
listed per accession in `docs/DATA.md`.
