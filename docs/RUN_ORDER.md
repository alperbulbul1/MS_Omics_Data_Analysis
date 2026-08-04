# Run order

Layer order: 00 -> 01/02 -> 03/04 -> 05 -> 06. Steps marked **[!]** are not safe to re-run blindly.

Only steps present in this release are listed. Filenames were cleaned for readability; the figure
generators are named after the manuscript figure they produce.

## Bulk RNA (transcriptome): raw GEO series -> two-track harmonisation -> ComBat -> per-stratum limma DE -> pan-tissue combined model -> Figure 2 / inverse-concordance inputs

**Download first**
- Microarray track, 6 GEO series (series_matrix.txt.gz + platform SOFT annotation): GSE21942 (GPL570), GSE43591 (GPL570), GSE38010 (GPL570), GSE103005 (GPL10558), GSE138064 (GPL17586), GSE190847 (GPL23126)
- RNA-seq track, 9 GEO series with deposited gene-level quantifications in supplementary / RAW.tar: GSE137143, GSE172009, GSE173789, GSE207680, GSE209596, GSE211358, GSE211739, GSE214334, GSE288904
- GSE66573 (10th RNA-seq series) - present in harmonized_v2 and added to the whole-blood stratum by volcano_NEWDATA_perstratum.R, but ABSENT from Corrected_Metadata_ComBat.csv and therefore absent from every Figure 2 panel; the manuscript Methods say 'nine RNA-seq series'
- ArrayExpress E-MTAB-69 (split into E-MTAB-69-CSF and E-MTAB-69-Blood) - MS vs OND lymphocyte panels; used only by volcano_NEWDATA_perstratum.R and dge_OND_modes.R, not by Figure 2
- GSM4071601_56514a-CD14.genes.results.txt.gz - a manually re-downloaded replacement for a corrupt member of the GSE137143 RAW.tar, hardcoded as a PATCH entry in harmonize_rnaseq_v3.py; must be documented as a separate acquisition step
- GEO family.soft.gz per series (Expression_Data/<GSE>/<GSE>_family.soft.gz) - required by run_stratified_omics.py to assign tissue/label strata; a missing SOFT file causes the series to be dropped from stratification with no error
- Used_Methylation_and_RNAseq_Dataset_Inventory.xlsx - local curation workbook read by rerun_verified_case_control.py; not a public artefact, must be shipped with the code
- Live web services at run time: NCBI GEO acc.cgi platform SOFT tables and mygene.info (harmonize_microarray_v2.py, build_matrices.py) - results are date-dependent
- Datasets appearing in legacy artefacts but not in the final cohort - GSE146383 (only in the overwritten Global_Harmonized_Metadata.csv), GSE235357/GSE247181/GSE255952/GSE130478 (targets in download_expression.py). Do not list these as manuscript inputs.

1. `01_transcriptome/harmonize_rnaseq_v3.py` (python) **[!]**  
   RNA-seq track harmonisation: reads 10 deposited RNA-seq supplementary/RAW.tar quantifications, maps ENSG/RefSeq/symbol -> HGNC, assigns MS/HC from series_matrix !Sample_characteristics_ch1 (falls back to !Sample_title regex).
   → <GSE>_symbol_matrix.csv, rnaseq_v2_metadata.csv
   > Contains a hardcoded PATCH dict re-pointing GSM4071601 of GSE137143 to an externally re-downloaded clean file, because the RAW.tar member is corrupt. That file is NOT obtainable from the GSE137143 RAW.tar and must be shipped as an acquisition instruction or the run silently loses

6. `01_transcriptome/harmonize_microarray_v2.py` (python) **[!]**  
   Microarray track harmonisation: parses 6 GEO array series_matrix expression tables plus E-MTAB-69, maps probe -> HGNC symbol (GPL570 via mygene, GPL10558/GPL17586/GPL23126 via live GEO platform SOFT download), assigns MS/HC from m
   → <GSE>_symbol_matrix.csv, microarray_v2_metadata.csv
   > NETWORK-DEPENDENT AND NON-DETERMINISTIC: probe->symbol maps are fetched live from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi (180 s timeout) and from mygene.info. Re-running at a later date can produce a different gene universe. The release must freeze the resolved probe->sym

11. `01_transcriptome/merge_harmonized_v2.py` (python) **[!]**  
   Merges the two harmonisation tracks on gene symbols present in >=50% of datasets into one symbol x sample matrix + combined metadata.
   → Global_Harmonized_v2_Expression.csv, Global_Harmonized_v2_Metadata.csv
   > On-disk output is STALE: Global_Harmonized_v2_Metadata.csv holds 727 samples / 16 datasets and contains NO E-MTAB-69 rows, while microarray_v2_metadata.csv (rewritten later, Jun 14 21:08) does contain E-MTAB-69-CSF and E-MTAB-69-Blood. Re-running changes the file. Nothing current

16. `01_transcriptome/build_matrices.py` (python) **[!]**  
   LEGACY v1 acquisition/parsing of individual GEO series into an expression matrix; superseded by the harmonize_*_v2/v3 track.
   → Recovered_Expression.csv, Recovered_Metadata.csv
   > Legacy. Retained only because the still-load-bearing Combined_Expression_Pre_ComBat.csv descends from it.

21. `01_transcriptome/build_global_matrices.py` (python) **[!]**  
   LEGACY v1 global matrix build from SOFT headers + platform tables.
   → Global_Harmonized_Expression.csv, Global_Harmonized_Metadata.csv
   > Its on-disk output was OVERWRITTEN (Jun 14 10:09) and now contains only 150 samples / 6 datasets (including GSE146383, which is in no other file). It therefore can no longer regenerate the 552-sample Combined_* matrix that the canonical chain depends on. Provenance is broken at t

26. `01_transcriptome/merge_and_correct.py` (python) **[!]**  
   LEGACY: merges Global_Harmonized_* with Recovered_* and writes the combined pre-ComBat matrix that the canonical stratified chain still consumes.
   → Combined_Expression_Pre_ComBat.csv, Combined_Metadata.csv
   > Cannot be re-run to reproduce the current Combined_* (552 samples / 15 datasets) because its Global_Harmonized_* input has been shrunk to 150/6. The 552-sample matrix on disk is effectively an orphan artefact.

31. `01_transcriptome/rerun_verified_case_control.py` (python) **[!]**  
   Filters Combined_* down to the datasets marked verified case/control in the inventory workbook.
   → Combined_Expression_Pre_ComBat.csv, Combined_Metadata.csv, Verified_Expression_Datasets.csv
   > DESTRUCTIVE IN-PLACE OVERWRITE: it reads EXPR_PRE / EXPR_META and writes back to the SAME paths. Running it twice, or running it on an already-filtered matrix, silently re-filters. It is not idempotent with respect to provenance and must be guarded in the release.

35b. `01_transcriptome/restore_zeroed_series.py` (python) **[!]**  
   MUST RUN BEFORE 36. Combined_Expression_Pre_ComBat.csv (16 March) carries four series as exactly zero across all 8,860 genes - GSE190847 (121 samples), GSE137143 (80), GSE172009 (8), GSE207680 (5) - because it predates the June re-harmonisation in harmonized_v2/. Restores them from the per-series harmonised matrices, mapping sample labels to GSM accessions via SOFT title / description / supplementary-filename. Asserts that no all-zero column survives and that the two zero-free strata (PBMC, whole blood) are byte-identical to the input.
   → Combined_Expression_Pre_ComBat_REPAIRED.csv
   > Skipping this step does not error: an all-zero series forms its own ComBat and limma batch, contributes a within-batch case-control contrast of exactly zero, and silently shrinks every effect while inflating residual df. The B-cell stratum loses 63.4% of its pooling weight. Step 36 now raises rather than continuing if any dataset is all-zero after normalisation.

36. `01_transcriptome/correct_and_normalize.py` (python) **[!]**  
   Per-dataset scale detection (raw counts vs log2), CPM+log2 for counts, within-dataset quantile normalisation, merge on genes present in >=75% of datasets, then neuroCombat batch correction with condition preserved; emits PCA valid
   → Corrected_Expression_Pre_ComBat.csv, Corrected_Metadata_ComBat.csv, Corrected_Batch_Corrected_Expression.csv
   > MUST NOT BE RE-RUN AS-IS. Its output Corrected_Metadata_ComBat.csv was hand-edited afterwards (Jun 14 22:04) to drop all 80 GSE137143 samples (552 -> 472 rows, 15 -> 14 datasets); the matching matrix Corrected_Expression_Pre_ComBat.csv (Jun 14 10:09) still carries all 553 columns

41. `01_transcriptome/run_expression_subgroup_limma.R` (R) **[!]**  
   Per-stratum limma worker invoked by run_stratified_omics.py. Fits ~ condition + dataset on the UNCORRECTED matrix (explicit in-code comment: 'use uncorrected mat for lmFit'), eBayes(robust=TRUE, trend=TRUE); separately writes a re
   → DGE_Results_MS_vs_HC.csv, Batch_Corrected_Expression.csv, Summary.csv
   > Takes paths as argv, no hardcoded paths - the cleanest script in the layer. Note the methodological split from Transcriptome/r_notebooks: this one tests on the uncorrected matrix with batch in the design; the r_notebooks re-test on the removeBatchEffect output. The two DE result 

46. `01_transcriptome/run_stratified_omics.py` (python) **[!]**  
   Defines the strata (cell/tissue case-control and label-context) by parsing GEO SOFT sample text, writes per-stratum metadata.csv, and shells out to run_expression_subgroup_limma.R for each; also drives the methylation strata.
   → , , 
   > main() re-runs EVERY stratum (expression and methylation). But on disk only cell_tissue_case_control_t_cells was re-run (Jun 14 22:04-22:05); all other strata are Jun 14 11:42-11:43. So the current state was produced by an ad-hoc partial invocation that does not exist as code. Al

51. `01_transcriptome/helpers.R` (R) **[!]**  
   Shared library for notebooks 01-08: PROJ_ROOT/OUT_DIR constants, CROSS_OMICS/RECURRING/PAPER_TOP gene panels, load_stratum(), run_limma_stratum() (eBayes trend+robust on the batch-corrected matrix), tx_volcano_gg().
   > Carries the OBSOLETE gene panels: CROSS_OMICS = LXN, SH3BP4, CHL1, CTSZ, RPAP2, PCNP, THRB - i.e. the old 7-gene cross-omics set. It does not contain ITGB2 or IKZF1, the current tier-1 pair. Every is_cross_omics column in Transcriptome/results/*.tsv is therefore wrong relative to

55. `01_transcriptome/00_run_all.R` (R) **[!]**  
   Driver that sources 01-08 in order.
   > BROKEN as documented. Its own header says 'Rscript 00_run_all.R', but line 1 of the body is HERE <- dirname(normalizePath(sys.frame(1)$ofile)), which under Rscript fails immediately with 'Error in sys.frame(1): not that many frames on the stack' (verified). It only works if itsel

59. `01_transcriptome/01_pbmc_de.R` (R) **[!]**  
   PBMC stratum limma DE re-run on the batch-corrected matrix.
   → 01_pbmc_DE.tsv, 01_pbmc_volcano.png, 01_pbmc_volcano.pdf
   > On-disk output verified consistent with the script (9 columns, lowercase 'gene', values match Unified_All_Assays_Long.tsv exactly for ITGB2/LXN/IKZF1).

64. `01_transcriptome/02_tcells_de.R` (R) **[!]**  
   T-cell stratum limma DE.
   → 02_tcells_DE.tsv, 02_tcells_volcano.png
   > THE ON-DISK 02_tcells_DE.tsv WAS NOT PRODUCED BY THIS SCRIPT. It is a hand-converted copy of Stratified_Analyses/Expression/cell_tissue_case_control_t_cells/DGE_Results_MS_vs_HC.csv: identical byte count (820,378), identical values, but CSV column order transposed to TSV, 7 colum

69. `01_transcriptome/03_bcells_de.R` (R)  
   B-cell stratum limma DE.
   → 03_bcells_DE.tsv, 03_bcells_volcano.png

73. `01_transcriptome/04_brainwm_de.R` (R)  
   Brain white-matter stratum limma DE (18 MS / 12 HC).
   → 04_brainwm_DE.tsv, 04_brainwm_volcano.png

77. `01_transcriptome/05_whole_blood_de.R` (R)  
   Whole-blood stratum limma DE.
   → 05_whole_blood_DE.tsv, 05_whole_blood_volcano.png

81. `01_transcriptome/06_pbmc_ifnb_de.R` (R) **[!]**  
   IFN-beta PBMC treatment-context stratum limma DE (45 MS-on-IFNb / 8 HC, single dataset so batch_method='none').
   → 06_pbmc_ifnb_DE.tsv, 06_pbmc_ifnb_volcano.png
   > The script labels this MS vs HC. Per the tiering decisions this stratum is a treatment-response contrast, not a disease contrast; the manuscript caption for Figure 2C now says so explicitly but the script's own text does not.

86. `01_transcriptome/07_total_combined_de.R` (R) **[!]**  
   Pan-tissue combined model: loads the 5 case-control strata, intersects genes, cbinds the matrices, and fits a SINGLE limma model ~ tissue + condition (model.matrix(~ tissue_f + condition_f), coefficient = last column), eBayes tren
   → 07_pan_tissue_DE.tsv, 07_pan_tissue_volcano.png, 07_pan_tissue_volcano.pdf
   > STALE AND CONTRADICTS THE MANUSCRIPT. Output is Jun 14 11:42 but its T-cell input was rebuilt Jun 14 22:04 - the output predates its own input. More seriously, it is NOT the source of Figure 2 panel G: 07_pan_tissue_DE.tsv gives IKZF1 logFC -0.029, FDR 0.573 (not significant), wh

91. `01_transcriptome/08_cross_stratum_master.R` (R) **[!]**  
   Cross-stratum summary: pulls the per-stratum TSVs, builds a gene x stratum logFC/FDR matrix, pheatmap with FDR stars, and a multi-stratum validated gene list.
   → CrossStratum_R_Summary.tsv, MultiStratum_Validated_R.tsv, CrossStratum_R_Heatmap.png
   > CONTAINS A SILENT FILENAME BUG: the assays list points the pan-tissue entry at OUT_DIR/'07_total_DE.tsv', but script 07 writes '07_pan_tissue_DE.tsv'. The file-exists guard turns this into a silent message('SKIP missing') - verified: CrossStratum_R_Summary.tsv contains only 6 str

96. `06_figures/figure2_rna_volcanoes.py` (python) **[!]**  
   MANUSCRIPT FIGURE 2 (image2.png, verified by md5 against Poster_v2/figures/per_celltype_volcanoes_INV.png). Renders panels A-F from the Stratified_Analyses per-stratum DGE CSVs and panel G from the pan-tissue DEG csv; holds the cu
   → per_celltype_volcanoes_INV.png, image2.png
   > This _figure_generators copy is CANONICAL and current (Aug 2). The sibling Poster_v2/make_volcano_v_INV.py is STALE: it still has the old 6-gene INV_TIER1 {HLA-E,ITGB2,LXN,CD79B,IKZF1,SH3BP4}, a 2x3 layout with no panel G, and the old suptitle. Do not ship the Poster_v2/ copy. No

**Notes:** BLOCKERS - things that would ship code that does not reproduce the paper.

1. ~~MISSING GENERATOR FOR FIGURE 2 PANEL G~~ — **RESOLVED**: panel G now reads `Transcriptome/results/07_pan_tissue_DE.tsv` (producer `01_transcriptome/07_total_combined_de.R`). Historical description follows. Poster_v2/figures/COMBINED_pantissue_proper_DEG.csv (Jun 15 01:03) is read verbatim by make_volcano_v_INV.py, but NO script in the project writes it. A full-tree grep for COMBINED_pantissue / COMBINED_allstrata / COMBINED_metaanalysis / properComBat matches only make_volcano_v_INV.py and its .pre* backups. Four sibling artefacts (COMBINED_allstrata_DEG.csv, COMBINED_metaanalysis_DEG.csv, COMBINED_pantissue_DEG.csv, COMBINED_properComBat_DEG.csv, all Jun 14-15) are likewise orph

## DNA methylation (array + WGBS): raw-IDAT/beta preprocessing, the two parallel ComBat M-value matrices, per-stratum limma DMP, mCSEA region enrichment, gene-level Stouffer aggregation (scripts 12-15), promoter-vs-body compartment test, and the Figure 3 generator.

**Download first**
- GSE130029 - Illumina 450K IDATs (23 samples: 12 MS / 11 HC); one of the six IDAT series
- GSE130030 - Illumina 450K IDATs (24 samples: 10 MS / 14 HC)
- GSE189255 - Illumina 450K IDATs, CD4 T cells (14 samples: 6 MS / 8 HC)
- GSE189256 - Illumina 450K IDATs, CD14 monocytes (14 samples: 6 MS / 8 HC)
- GSE219293 - Illumina 450K IDATs (47 samples: 29 MS / 18 HC)
- GSE88824 - Illumina 450K IDATs, whole blood + sorted cells (27 samples: 13 MS / 14 HC)
- GSE106648 - peripheral blood 450K, NO IDATs deposited; depositor methylated/unmethylated signal-intensity matrices (279 samples: 140 MS / 139 HC; 59% of the 475-sample merge)
- GSE40360 - brain white matter 450K, NO IDATs deposited; Illumina GenomeStudio export (47 samples: 28 MS / 19 HC)
- GSE151017 - BAL 450K (76 samples: 32 MS / 44 HC). PRESENT in Combined_Methylation_Batch_Corrected.csv (Track A, 548 samples) but deliberately EXCLUDED from AllMeth_ComBat_M.csv (Track B, 475 samples). The release must state this asymmetry.
- GSE173787 - sorted-immune WGBS (hg38), used for the orthogonal promoter-window validation

1. `00_data/select_methylation_case_control.py` (python) **[!]**  
   Filters the curated GEO inventory spreadsheet down to the case-control methylation series; emits the target list every downloader reads.
   → Methylation_Target_Datasets.csv
   > seconds; pure table filtering

6. `00_data/download_methylation.py` (python) **[!]**  
   Scrapes and downloads GEO supplementary files (IDATs, series matrices, signal-intensity matrices) for every accession in Methylation_Target_Datasets.csv.
   →  (raw IDAT + supplementary files)
   > network-bound, tens of GB; ThreadPoolExecutor; not deterministic if GEO changes hosting layout

11. `02_methylation/build_methylation_matrix.py` (python) **[!]**  
   Parses the downloaded per-series matrices, harmonises sample metadata (condition MS/HC, cell_type, dataset) and stacks them into the pre-batch probe x sample matrix.
   → Combined_Methylation_Pre_Batch.csv, Combined_Methylation_Metadata.csv
   > 446 MB output; MIN_FEATURE_OVERLAP=5000 hardcoded; condition/cell-type assignment is keyword-matching on free-text GEO characteristics, so GEO metadata drift silently changes group labels

16. `02_methylation/preprocess_methylation_arrays.R` (R) **[!]**  
   minfi IDAT pipeline: RGChannelSet -> detection-p and bead-count QC -> preprocessIllumina -> SNP / cross-reactive / sex-chromosome probe removal -> per-dataset beta + M -> merged Combined_Methylation_Strict_{Beta,M,Metadata}.csv in
   → Combined_Methylation_Strict_{Beta,M,M_BatchCorrected,Metadata}.csv, Dataset_QC_Summary.csv, Combined_Methylation_Strict_M.csv (second invocation, 6-series subset, 149 samples)
   > hours; multi-GB outputs. Run TWICE with different argv: once with default out_dir (9-series Strict_Array_Preprocessed, includes BAL) and once with out_dir=IDAT6_Preprocessed over a 6-series metadata subset. Neither invocation is recorded anywhere in the repo - the argv must be re

21. `02_methylation/normalize_beta_only.R` (R) **[!]**  
   Handles the series with no deposited IDATs: within-dataset BMIQ (Type I/II probe-bias correction), >20% missingness filter, SNP + sex-chromosome removal, beta -> M conversion.
   → Normalized_{Beta,M}.csv, BetaOnly_QC_Summary.csv
   > wateRmelon is absent, so BMIQ is skipped and the script silently falls back to plain quantile normalisation - the manuscript Methods explicitly claim BMIQ. This is a substantive, not cosmetic, mismatch.

26. `01_transcriptome/run_stratified_omics.py` (python) **[!]**  
   Stratum definition + dispatcher. Builds the cell_tissue_case_control_* and label_context_case_control_* sample subsets, writes each stratum's metadata.csv, and shells out to run_methylation_subgroup_limma.R per stratum.
   → metadata.csv, Methylation_Subgroup_{Definitions,Results_Summary}.csv, Methylation_Sample_Subgroups.csv
   > Prefers Strict_Array_Preprocessed and falls back to Combined_Methylation_Pre_Batch.csv - the fallback would produce a DIFFERENT stratified matrix without any error. Also drives the expression layer in the same run.

31. `02_methylation/run_methylation_subgroup_limma.R` (R) **[!]**  
   Per-stratum worker: top-variance probe prefilter, ComBat/removeBatchEffect for QC, limma DMP with batch covariate, DMRcate DMRs, mCSEA promoters + genes, missMethyl GO/KEGG. Writes each stratum's Batch_Corrected_M.csv.
   → Batch_Corrected_M.csv  [100,000 probes x n], DMP_Results_MS_vs_HC.csv, Promoter_Results_mCSEA.csv
   > CRITICAL: the current file has MAX_PROBES <- 450000L, but every shipped Batch_Corrected_M.csv has exactly 100,000 probes (script mtime Mar 17 10:19 POSTDATES the matrices, Mar 17 00:36). Re-running would produce 4.5x more probes and different BH-FDRs than the manuscript reports. 

36. `01_transcriptome/helpers.R` (R) **[!]**  
   Shared library for notebooks 01-15: path constants, load_meth_stratum(), load_combined_meth(), run_limma_meth(), annotate_probes_to_genes(), probe_to_gene_stouffer(), meth_volcano_gg(), the CROSS_OMICS / METH_TOP panels.
   → (library - no files)
   > THE VERIFIED-STALE FILE. Canonical Methylation/r_notebooks/helpers.R (Aug 2 02:12) adds d$SE and the se_logFC / ivw_logFC / min_logFC / max_logFC columns plus a 22-line provenance comment; MS_GEO_pipeline/scripts/02_methylation/Methylation__r_notebooks/helpers.R (May 17) lacks al

41. `02_methylation/01_tcells_meth_dmp.R` (R) **[!]**  
   T-cell case-control stratum: limma DMP on the stratum M-matrix with study as covariate, probe->gene annotation, Stouffer gene aggregation, probe + gene volcanoes.
   → 01_tcells_meth_DMP.tsv, 01_tcells_meth_gene.tsv, 01_tcells_meth_volcano_{probe,gene}.png
   > Contains a live BiocManager::install() call for the 450k annotation package - must be removed or guarded in a release.

46. `02_methylation/02_wb_dmf_meth_dmp.R` (R)  
   Whole-blood dimethyl-fumarate stratum, same pipeline as 01. Source of the manuscript's 2,650/16,667 DMF gene-level count and the THRB / CTSZ / CHL1 anchors.
   → 02_wb_dmf_meth_DMP.tsv, 02_wb_dmf_meth_gene.tsv, 02_wb_dmf_meth_volcano_{probe,gene}.png

50. `02_methylation/03_wb_ocrelizumab_meth_dmp.R` (R)  
   Whole-blood ocrelizumab stratum (99 gene-level DMPs, 12 mCSEA promoters in the manuscript).
   → 03_wb_ocrelizumab_meth_DMP.tsv, 03_wb_ocrelizumab_meth_gene.tsv, 03_wb_ocrelizumab_meth_volcano_{probe,gene}.png

54. `02_methylation/04_tcells_remission_meth_dmp.R` (R)  
   T-cell remission stratum. Not shown in Figure 3 (omitted by the figure generator) but its DMP/gene TSVs are consumed by scripts 12-15.
   → 04_tcells_remission_meth_DMP.tsv, 04_tcells_remission_meth_gene.tsv, 04_tcells_remission_meth_volcano_{probe,gene}.png

58. `02_methylation/05_combined_meth_dmp.R` (R) **[!]**  
   TRACK A pan-tissue combined stratum. limma on all MS/HC samples of Combined_Methylation_Batch_Corrected.csv with dataset as batch covariate, then Stouffer gene aggregation. This is the stratum whose gene-level FDRs the manuscript 
   → 05_combined_meth_DMP.tsv, 05_combined_meth_gene.tsv, 05_combined_meth_volcano_{probe,gene}.png
   > Loads a 332 MB matrix column-selectively via fread(select=); ~30 s read plus a large lmFit. Header comment says '~549 samples across 17 GSE series' - the actual matrix is 548 samples from 9 series; the comment is wrong and should not be shipped as documentation.

63. `02_methylation/06_mcsea_promoter_analysis.R` (R) **[!]**  
   mCSEA promoter and gene-body region enrichment on the same 9-series combined matrix. Source of every mCSEA promoter q-value in the manuscript (ITGB2 q=0.016, IKZF1 q=0.030, LXN q=0.0024, CASP8 q=0.019, SH3BP4 q=0.18).
   → 06_mCSEA_promoter.tsv, 06_mCSEA_gene_body.tsv
   > Reloads the full 332 MB matrix independently of script 05 (no caching). Contains a live BiocManager::install("mCSEA"). Gene-body call is wrapped in try() and fails silently on error.

68. `02_methylation/07_brainwm_rna_vs_meth.R` (R) **[!]**  
   Brain white-matter RNA x methylation inverse-concordance scan; supplies the Brain WM panel of Figure 3.
   → 07_BrainWM_RNA_vs_Meth_concordance.tsv, 07_BrainWM_RNA_vs_Meth_concordance.png
   > CROSS-LAYER DEPENDENCY: both inputs live under Proteomics/processed/ and are produced by the transcriptome/proteomics layer, not here. It has a stopifnot(file.exists(...)) guard, so it hard-fails if that layer has not been run. It does NOT use the brain_wm stratum directory, whic

73. `02_methylation/08_cross_stratum_meth_master.R` (R) **[!]**  
   Assembles the 01-05 gene-level TSVs into a cross-stratum summary table and pheatmap of the cross-omics + meth-top panels.
   → CrossStratum_Meth_R_Summary.tsv, CrossStratum_Meth_R_Heatmap.png
   > Silently skips any missing stratum TSV (message + next) rather than failing, so a partial run yields a quietly incomplete summary. This table is cited for the CTSZ / LXN DMF FDRs.

78. `02_methylation/09_inverse_concordance_scan.R` (R) **[!]**  
   The discovery step: scans all 7 transcriptome strata against the two methylation sources (05 gene-level Stouffer, 06 mCSEA promoter NES) for joint-FDR<0.05 opposite-sign genes. Produces the 82-gene inverse-concordant discovery poo
   → INVERSE_CONCORDANT_full_pairings.tsv, INVERSE_CONCORDANT_by_gene.tsv, INVERSE_CONCORDANT_by_gene_by_stratum.tsv
   > CROSS-LAYER: requires the full Transcriptome/results/ set. Uses mean_logFC = NES * 0.1 as a pseudo-effect-size for the mCSEA source - an arbitrary rescale that only affects the sign test, but it is undocumented in the manuscript. INVERSE_CONCORDANT_by_gene.tsv is the '82 genes' r

83. `02_methylation/10_inverse_proteomics_validation.R` (R) **[!]**  
   Looks the top-40 inverse-concordant genes up in every R proteomics result table and builds the gene x assay validation heatmap.
   → INV_proteomics_validation_long.tsv, INV_proteomics_validation_by_gene.tsv, 10_inverse_proteomics_heatmap.png
   > STALE INPUT PATHS. It reads CSF_Astral_R_DEP_results.tsv / CSF_timsTOF_R_DEP_results.tsv, i.e. the pre-revision imputed DEP outputs. The proteomics layer moved to complete-case *_CC_*.tsv this round. This script has NOT been updated and would re-introduce the rejected MinProb-imp

88. `02_methylation/promoter_vs_body_test.R` (R) **[!]**  
   Reviewer point 2: explicit promoter-only vs gene-body-only signed-Stouffer statistic for all pool genes in all 5 strata, classifying probes via UCSC_RefGene_Group. Bypasses mCSEA's minimum-region-size limitation.
   → promoter_vs_body_by_gene.csv
   > MUST be in the release: it is the sole source of the Section 2.4 facts that CD79B has only 3 promoter probes and that SH3BP4's gene-body FDR=0.016 beats its promoter FDR=0.13 in the DMF stratum. No copy exists under MS_GEO_pipeline/. Hardcoded RES and OUT paths.

93. `02_methylation/update_S2_promoter_flag.py` (python) **[SUPERSEDED — do not run]**  
   Adds the methylation-compartment classification (promoter-confirmed / promoter-only / composite-only), mCSEA promoter and gene-body padj columns to Supplementary Table S2 and runs the promoter-only sensitivity count.
   → Supplementary_Tables_IJMS_v2.xlsx
   > SUPERSEDED by step 104 (`rebuild_S2_94genes.py`). It patches columns onto the older 82-gene
   > sheet and produces the withdrawn '6 of 82 promoter-anchored / 76 composite-only' figures.
   > Running it regenerates the superseded v2 workbook; it is retained only for provenance.

98. `02_methylation/15_genelevel_weighting_corrected.R` (R) **[!]**  
   CANONICAL gene-level weighting analysis. Recomputes, per stratum, the published unweighted-Stouffer Z, the consistent 1/SE-weighted effect size, the IVW secondary estimator, the per-gene probe range, and the dependence-aware Z wit
   → 15_genelevel_weighting_corrected.tsv
   > Memory-heavy: builds a dense n x n hat matrix and a full residual matrix per stratum, including the 548-sample combined matrix; MAXPR=40 caps probes per gene in the pairwise-correlation step, so rho_sum is extrapolated for genes with >40 probes. Hardcoded setwd() and RES. NOT in 

103. `02_methylation/update_S2_weighting.py` (python) **[SUPERSEDED — do not run]**  
   Writes the six weighting columns promised in Methods and in the Reviewer-1 reply into Supplementary Table S2: CpG probes (n), unweighted mean logFC, 1/SE-weighted, inverse-variance, probe range, gene-level FDR - all from the 05_co
   → Supplementary_Tables_IJMS_v2.xlsx (MODIFIED IN PLACE)
   > SUPERSEDED by step 104 (`rebuild_S2_94genes.py`), which builds the whole sheet in one pass.
   > Retained for provenance only. It was also NOT IDEMPOTENT: it appended columns at
   > ws.max_column+1 and appended a sentence to the A1 caption on every invocation, so a second
   > run duplicated six columns and duplicated the caption text.

104. `02_methylation/rebuild_S2_94genes.py` (python) **[!]**  
   CANONICAL Supplementary Table S2 builder. Rebuilds the whole sheet in one pass for the current
   94-gene inverse-concordant pool: per-gene statistics, the methylation-compartment flag
   (promoter-confirmed / promoter-only / composite-only), the mCSEA promoter and gene-body padj
   columns, the promoter-only survival flag, and the CpG-count and weighting columns from the
   05_combined stratum.
   → Supplementary_Tables_IJMS_v3.xlsx
   > Run AFTER steps 06 (mCSEA) and 98 (`15_genelevel_weighting_corrected.R`), whose outputs it
   > reads, and instead of steps 93 and 103. Idempotent: it deletes and recreates the sheet rather
   > than appending, and asserts the pool is exactly 94 genes. Reports the promoter-only
   > sensitivity analysis to stdout: 8 of 94 survive, including both Tier-1 genes.

108. `02_methylation/run_all_methylation_combat.R` (R) **[!]**  
   TRACK B canonical builder. Merges the 6-series IDAT M-matrix with the two beta-only series on common 450K probes, runs sva::ComBat with condition preserved, and runs limma MS-vs-HC with best-probe-per-gene aggregation.
   → AllMeth_ComBat_M.csv  [VERIFIED 270,123 probes x 475 samples, 8 series, BAL excluded],  231 HC], AllMeth_ComBat_limma_DMP_byGene.csv
   > 2.2 GB output; reads three multi-GB CSVs with base read.csv (not fread) - very slow and memory-hungry. BLOCKING GAP: nothing in the repository produces GSE106648_betaonly_M.csv or GSE40360_betaonly_M.csv. normalize_beta_only.R writes Normalized_Beta_Only/GSE*/Normalized_M.csv ins

113. `02_methylation/run_mcsea_combat.R` (R) **[!]**  
   mCSEA promoter + gene-body enrichment on the 8-series ComBat matrix, as a cross-check of the Track A mCSEA result.
   → AllMeth_ComBat_mCSEA.rds, AllMeth_ComBat_mCSEA_promoters.csv, AllMeth_ComBat_mCSEA_genes.csv
   > Loads the full 2.2 GB matrix into a dense R matrix; expect very high peak RAM. Its INV1/T2A gene lists are the OLD 6-gene Tier-1 (HLA-E, ITGB2, LXN, CD79B, IKZF1, SH3BP4) - console labelling only, no effect on the written tables.

118. `06_figures/figure3_methylation.py` (python) **[!]**  
   FIGURE 3 GENERATOR. Five gene-level methylation volcanoes on a common logFC axis: panels A-D from the Track A per-stratum/combined gene TSVs plus the Brain WM concordance table, panel E from the Track B 475-sample ComBat result, p
   → methylation_v_INV.png, image3.png (md5-verified identical)
   > DIRECTION REVERSAL - the canonical copy is the one INSIDE MS_GEO_pipeline/scripts/_figure_generators/ (mtime Aug 2 15:20, INV_TIER1={ITGB2,IKZF1}, TIER2_AUX_INV includes CD79B/LXN/SH3BP4). The sibling Poster_v2/make_methylation_v_INV.py is STALE (Jun 7, 6-gene Tier-1 including HL

**Notes:** TWO PARALLEL MATRICES - definitive attribution (all verified by reading write calls and by counting rows/columns on disk, not inferred):

(A) Methylation_Data/Combined_Methylation_Batch_Corrected.csv - 34,196 probes x 548 samples, NINE series (includes GSE151017 BAL). SOLE WRITER: methylation_analysis_pipeline.R line 292 (mtime Mar 16 11:17 matches its Analysis_Results run). CONSUMERS: notebooks 05 and 06 (which produce the reported gene-level FDRs and mCSEA q-values), 14 and 15, and helpers.R::load_combined_meth. Feeds Figure 3 panels A-D.
  DANGER: run_methylation_limma.R (line 43) writes th

## Proteomics (Proteomics/r_notebooks/ complete-case differential abundance + Poster_v2/build_RDEP_CC_adapters.py figure adapters)

**Download first**
- Bader & Mann CSF Orbitrap Astral DIA-MS proteome - raw spectra PRIDE PXD064570; the analysis does NOT use PRIDE, it uses the OSF-deposited PROCESSED matrix from that study (Proteomics/osfstorage-archive/) reduced to Proteomics/processed/astral_discovery_gene_keyed.tsv (83 MB). No script in the repository performs that reduction.
- Bader & Mann CSF timsTOF DIA proteome - raw spectra PRIDE PXD045058; analysis input is Proteomics/processed/timsTOF_gene_mapped.tsv (213 MB), again derived from the OSF processed data with no script on record.
- Bader & Mann sample annotation table: 'osfstorage-archive/processed proteomic data/0_sample_annotations/annotations_v42_49_2_10_4_10_interimSky17_PL01-PL56_PepResCustv01_resubmission.tsv'. Required by 01cc, 02cc, 03 and 11. Script 11 additionally needs its clinical columns Leukocyte_count, Erythrocytes_in_CSF and Total_protein - these exist only in this OSF file, not in the PRIDE deposition, and the local copy sits in a 0700 directory. The release must give the exact OSF URL/DOI and note any access conditions on the clinical fields.
- Wang & Julien 2026 brain white-matter proteome (Nat Commun, DOI 10.1038/s41467-025-68118-0; raw spectra MassIVE MSV000096790, processed tables on Figshare). Local file Proteomics/Magliozzi2026_BrainProteomics_S1.xlsx, sheet '3-filtered-detection-proteins'. 04cc consumes the derived Proteomics/processed/Magliozzi_S1_sheet3.tsv.
- UK Biobank-PPP Olink plasma MS association statistics (Jacobs et al. 2024, Ann Clin Transl Neurol) - published supplementary, local copies Proteomics/blood_raw/ACN3-11-698-s002.xlsx / jacobs_supp.zip, flattened to Proteomics/blood_raw/Jacobs2024_UKB_primary_DE.tsv. Consumed directly by the Figure 4 and Figure 6 generators; no extraction script exists.
- GEO GSE32915 (Pegram, NK8+) and GSE78244 (Hellberg, CD4 T) - downloaded at runtime by GEOquery in the superseded scripts 05 and 06 only; not needed for any current manuscript figure.
- Hand-curated published-proteomics literature meta-cohort: Proteomics/MS_Serum_Proteomics_CaseControl_Expanded.xlsx, MS_BCell_OtherImmune_Proteomics.xlsx, MS_TCell_Proteomics_CaseControl.xlsx. These are the authors' own curation, not downloadable data - they must ship WITH the release or the '30-study meta-cohort' statement is unreproducible.

1. `01_transcriptome/helpers.R` (R) **[!]**  
   Shared library sourced by every r_notebooks script. Defines PROJ_ROOT/PROT_ROOT/OUT_DIR/FIG_DIR/CACHE_DIR, the candidate gene sets (CROSS_OMICS/RECURRING/PAPER_TOP/ECM_FAMILY), and the DEP-EQUIVALENT functions that actually run in
   > Not standalone; must be sourced with cwd = Proteomics/r_notebooks. Hardcodes PROJ_ROOT at line 21.

5b. `03_proteomics/dep_bh_equivalence_check.R` (R)  
   VERIFICATION, not part of the reported pipeline. Runs DEP itself on the CSF Astral matrix with DEP::impute omitted, and compares it with the reported complete-case limma result. Establishes that the two are the same test: identical fold changes (Pearson r = 1.000), identical p-value ranking (Spearman rho = 0.9999), and, once DEP's p-values are BH-adjusted, 955 significant proteins against the 941 reported, agreeing on 1,983 of 1,995 shared proteins (99.4%).
   → Proteomics/processed/META/CSF_Astral_DEP_BH_check.tsv
   > DEP::test_diff hard-codes fdrtool for multiple-testing adjustment and exposes no argument to change it; on its own default it calls only 35 proteins. Every layer of this study reports Benjamini-Hochberg, which is why limma is used directly. DEP was removed from Bioconductor at release 3.23; install with remotes::install_github("arnesmits/DEP"). Two pitfalls are handled in the script: make_se and normalize_vsn both assume raw intensities and will silently compress fold changes on this already-log matrix.

5. `03_proteomics/01cc_csf_astral_completecase.R` (R) **[!]**  
   CANONICAL Astral CSF (Bader & Mann) MS-vs-control differential abundance, complete-case (no imputation). Matches annotation to .raw run columns, builds MS/Control groups from MSgroup + Diagnosis_group, gene-dedups by max variance,
   → CSF_Astral_CC_results.tsv, CSF_Astral_CC_volcano.png, CSF_Astral_CC_volcano.pdf
   > 83 MB input matrix read with fread; a few minutes and several GB RAM. Must be run with cwd = Proteomics/r_notebooks (source("helpers.R") is relative). Verified output: 1,995 genes, 941 at FDR<0.05 - exactly the manuscript's Section 2.5 numbers. The volcano SUBTITLE still says 'vs

10. `03_proteomics/02cc_csf_timstof_completecase.R` (R) **[!]**  
   CANONICAL timsTOF CSF (same Bader & Mann study, second instrument) complete-case differential abundance. Identical method to 01cc; discovers the timsTOF run column in the annotation table by regex.
   → CSF_timsTOF_CC_results.tsv, CSF_timsTOF_CC_volcano.png, CSF_timsTOF_CC_volcano.pdf
   > 213 MB input matrix (fread nThread=4); the heaviest script in the layer, high peak RAM. Verified output: 1,069 genes, 759 at FDR<0.05 - matches the manuscript. Same wrong 'vsn + MinProb' figure subtitle. dep_completecase_de is copy-pasted from 01cc rather than shared via helpers.

15. `03_proteomics/04cc_magliozzi_brain_completecase.R` (R) **[!]**  
   CANONICAL brain white-matter proteome (cited in the manuscript as Wang & Julien 2026 Nat Commun, DOI 10.1038/s41467-025-68118-0, but named 'Magliozzi' throughout the code) complete-case differential abundance across 4 region contr
   → Magliozzi_CC_MS_CTX_vs_ODC_CTX.tsv, Magliozzi_CC_MS_NAWM_vs_ODC_WM.tsv, Magliozzi_CC_MS_WML_vs_ODC_WM.tsv
   > Fast (3,504 genes x 37 samples). IMPORTANT: unlike 01cc/02cc it does NOT apply an explicit '>=2 measured values per group' guard - with min group n=5 the 50% filter implies >=3 measured, so the Methods claim holds implicitly, but the code does not enforce it. Also unlike the supe

20. `03_proteomics/11_itgb2_csf_pleocytosis.R` (R) **[!]**  
   NEW this revision (Reviewer 1). Models ITGB2 DETECTION (not abundance) in Astral CSF as a function of CSF leukocyte count, adjusting for total CSF protein, CSF erythrocytes and diagnosis; stratified detection rates; specificity co
   → ITGB2_CSF_pleocytosis.txt
   > Self-contained: does NOT source helpers.R, hardcodes PROT_ROOT at line 28, so it can run from any cwd. Requires the annotation columns Leukocyte_count, Erythrocytes_in_CSF and Total_protein - clinical metadata that lives only in the OSF-deposited annotation TSV, not in PRIDE. Use

25. `03_proteomics/build_RDEP_CC_adapters.py` (python) **[!]**  
   CANONICAL figure-adapter builder. Rewrites the six complete-case result tables into the legacy column layout (Genes/log2FC_MSvsCtrl|log2FC/t/pval/FDR for CSF; Gene/log2FC/t/p/FDR for brain) under Proteomics/processed/RDEP_CC/, whi
   → Astral_RDEP.tsv, timsTOF_RDEP.tsv, Magliozzi2026_MS_CTX_vs_ODC_CTX.tsv
   > Seconds. Hardcodes both input and output absolute paths (lines 23-24). MUST run after 01cc/02cc/04cc and BEFORE the Figure 4/6 generators. On-disk timestamps confirm this order was respected (CC tables 02 Aug 00:02, RDEP_CC 00:11, figures 15:21).

30. `03_proteomics/03_csf_cross_platform_meta.R` (R) **[!]**  
   SUPERSEDED / STILL IMPUTED. Astral x timsTOF cross-platform integration by sva::ComBat(batch=platform, mod=~group) then limma. It was NOT converted to complete-case: lines 110-111 still call impute_minprob_R on both platforms. Its
   → CSF_combined_R_ComBat_DE.tsv, CSF_combined_R_PCA.png, CSF_combined_R_volcano.png
   > Loads both large matrices simultaneously - highest memory footprint of the layer. If the release ships this script it must be labelled as the MinProb-era cross-platform meta, or converted, otherwise it contradicts Methods 4.5 ('Missing intensities were not imputed').

35. `03_proteomics/05_t_lineage_meta.R` (R) **[!]**  
   SUPERSEDED / ORPHANED. T-lineage microarray meta of GSE32915 (Pegram, NK8+) + GSE78244 (Hellberg, CD4 T, unstimulated subset): GEOquery -> quantile normalise -> avereps -> ComBat(batch=study) -> limma. The current Figure 6 generat
   → T_lineage_R_combined_DE.tsv, T_lineage_R_pca.png, T_lineage_R_volcano.png
   > Network-dependent (GEOquery); results cached as _cache/<GSE>_eset.rds. Feeds script 09 only.

40. `03_proteomics/06_pegram_gse32915_de.R` (R) **[!]**  
   SUPERSEDED / ORPHANED. Standalone limma DE of GSE32915 (4 MS / 4 control, technical replicates). Feeds script 09 only; no current manuscript figure.
   → Pegram_R_DE_gene.tsv, Pegram_R_DE_volcano.png
   > Network-dependent; shares the GSE32915 cache with 05.

45. `03_proteomics/07_brainwm_rna_meth_rerun.R` (R) **[!]**  
   SUPERSEDED / CROSS-LAYER. R-side brain white-matter RNA DE plus an inverse-concordance plot against a methylation gene-level table produced by the methylation layer. Belongs to the proteomics notebook series only by filing acciden
   → BrainWM_R_RNA_DE.tsv, BrainWM_R_RNA_volcano.png, BrainWM_R_meth_volcano.png
   > Its methylation input predates the corrected gene-level weighting (Methylation script 15) - cross-check with the methylation layer owner before shipping. Feeds 09/10 only.

50. `03_proteomics/08_per_group_consistency.R` (R) **[!]**  
   SUPERSEDED-ish / CURATION LAYER. Splits the hand-curated published-proteomics meta-cohort into 29 tissue x cell-type groups and writes per-group consistency tables plus a summary. This is the code behind the manuscript's '30-study
   → consistency_R.tsv x29, _group_summary_R.tsv, all_groups_grid_R.png
   > Its input all_DEPs_with_group.csv is built by Proteomics/scripts/01_load_studies_by_group.py from three hand-curated xlsx workbooks; that CSV carries 64 distinct study_id values across 29 group directories, which does not obviously equal the '30 studies' the manuscript claims - t

55. `03_proteomics/09_cross_assay_lxn.R` (R) **[!]**  
   SUPERSEDED - READS MinProb TABLES. Cross-assay tracking grid over 10 result tables. Lines 26-38 point at CSF_Astral_R_DEP_results.tsv, CSF_timsTOF_R_DEP_results.tsv, CSF_combined_R_ComBat_DE.tsv and the four Magliozzi_R_DEP_* tabl
   → CrossAssay_R_summary.tsv, CrossAssay_R_grid.png
   > Either repoint at the *_CC_* tables or exclude from the release. Do not ship unchanged as 'the code that produced the paper'.

60. `03_proteomics/10_master_validation.R` (R) **[!]**  
   SUPERSEDED - READS MinProb TABLES. Builds the cross-omics master table and the 'triple validated' (CSF protein x brain RNA x brain methylation) intersection. Line 57 reads CSF_Astral_R_DEP_results.tsv (MinProb).
   → CrossOmics_Master_Table_R.tsv, Triple_Validated_R.tsv, CrossOmics_Master_Heatmap_R.png
   > Triple_Validated_R.tsv is pulled into the supplementary workbook MS_GEO_Master_Results.xlsx as sheet 'Prot_TripleValidated' by MS_GEO_pipeline/scripts/_manuscript_utils/rebuild_master_excel.py - so MinProb-era proteomics is currently embedded in that workbook.

65. `01_transcriptome/00_run_all.R` (R) **[!]**  
   SUPERSEDED DRIVER. Sources scripts 08, 04, 06, 05, 01, 02, 03, 07, 09, 10 in a dependency-aware order. It predates the revision: it does NOT know about 01cc/02cc/04cc/11 or build_RDEP_CC_adapters.py, and it runs the MinProb varian
   > Uses sys.frame(1)$ofile to setwd, so it only works under source(), not Rscript's usual invocation. Must be replaced by a new driver in the release: helpers -> 01cc -> 02cc -> 04cc -> 11 -> build_RDEP_CC_adapters.py.

**Notes:** SNAPSHOT STATUS FOR THIS LAYER. MS_GEO_pipeline/scripts/04_proteomics/Proteomics__r_notebooks/ is byte-identical to the canonical directory for all 14 files it contains (00_run_all.R, 01-10 .R, all 10 .ipynb, helpers.R, build_specs.py, _build_ipynb.py, README.md - verified by diff/cmp). So for proteomics the snapshot is not corrupted, it is INCOMPLETE: it is missing 01cc, 02cc, 04cc, 11 and build_RDEP_CC_adapters.py, i.e. every script that produces the numbers the revised manuscript reports. Shipping that directory would ship the rejected MinProb analysis in full and none of the corrected one.

## Single-cell / donor-level pseudobulk (scRNA-seq layer)

**Download first**
- GSE144744 (Kaufmann 2021, PBMC 10x Chromium 3'/CITE-seq) - GSE144744_RNA_normalised.tar.gz: RNA_normalised/matrix.mtx (16.1 GB, 529,680,471 non-zeros, 15,354 genes x 497,705 cells), genes.tsv, barcodes.tsv, plus cell_meta.csv.gz (donor, basictype, cluster_names, group, cohort, batch_pair, nCount_RNA, sex, age_sampling). Also GSE144744_saverx.tar.gz (466 MB, SaverX-imputed 425-gene panel) for the UMAP. Note the whole tree calls this cohort 'Ramesh2020' - it is Kaufmann.
- GSE118257 (Jaekel 2019, brain snRNA-seq) - expr.txt.gz (raw integer UMI counts) and anno.txt.gz
- GSE127969 (Beltran 2019, MS-discordant monozygotic twin CSF + PBMC) - deposited TPM-like normalised values only; raw counts are NOT recoverable, which is why Beltran is modelled with limma-trend rather than a count model
- GSE180759 (Absinta 2021, brain WM snRNA) - used by Methylation/r_notebooks/13
- GSE138266 (Schafflick 2020, PBMC/CSF) - the local h5ad is corrupt; script 13 silently falls back to pre-computed stats CSVs from an earlier run, which are an undocumented input
- GSE239626 (Vitamin D MS trial, blood) - loaded by script 11 but explicitly excluded from script 13 (no HC arm)
- mygene.info query results backing Poster_v2/figures/pseudobulk_proper/gene_biotypes.csv (protein-coding filter); no query script exists
- Methylation/results/INVERSE_CONCORDANT_by_gene.tsv - cross-layer input from the methylation layer, consumed by scripts 11, 12, 13 and cell_level_3cohorts.py

1. `04_singlecell/download_singlecell_datasets.py` (python) **[!]**  
   Downloads the CELLxGENE/GEO single-cell source datasets into SingleCell_CELLxGENE/data/. Acquisition step only - no statistics.
   → , , 
   > Network-bound; GSE144744 RNA_normalised alone is 16.1 GB (matrix.mtx). Total data/ tree is 18 GB.

6. `04_singlecell/process_GSE118257_jakel_brain.py` (python) **[!]**  
   Builds the Jaekel 2019 brain snRNA AnnData (QC, annotation carried from the depositor's anno.txt) used by every downstream single-cell step.
   → adata_jakel.h5ad
   > expr.txt.gz is a dense text matrix; peak RAM several GB.

11. `04_singlecell/process_GSE127969_beltran_csf.py` (python)  
   Builds the Beltran 2019 MS-discordant twin CSF+PBMC AnnData. Deposited values are TPM-like normalised; raw counts are not recoverable (this constrains all downstream Beltran modelling).
   → adata_beltran.h5ad

15. `04_singlecell/plot_GSE144744_kaufmann_celltypes.py` (python)  
   Re-renders the Kaufmann UMAP figure panels from the existing h5ad without re-running the UMAP embedding. Safe to run; cosmetic only.
   → fig*.png

19. `02_methylation/11_inverse_scRNA_validation.py` (python) **[!]**  
   First-pass inverse-concordant panel re-test across seven local AnnData datasets (Wilcoxon MS-vs-HC per gene x dataset x cell type, BH + Cohen's d). Superseded by script 13 for the per-study/per-tissue results but still the produce
   → INV_scRNA_validation_long.tsv, INV_scRNA_validation_by_gene.tsv, 11_inverse_scRNA_heatmap.png
   > Loads seven h5ads sequentially; several GB peak. Script header documents known defects fixed in script 13 (VitD wrongly treated as MS-vs-HC; Beltran CSF and PBMC lumped; brain studies pooled across tissues).

24. `02_methylation/13_perstudy_scRNA_validation.py` (python) **[!]**  
   CANONICAL per-study, tissue- and cell-type-aware scRNA validation of the CO7 + top-30 inverse-concordant panel. Fixes the three defects in script 11. Produces the INV_scRNA_* tables named in the task.
   → INV_scRNA_per_study_long.tsv, INV_scRNA_per_study_by_gene.tsv, INV_scRNA_per_gene_per_tissue.tsv
   > Schafflick PBMC h5ad is corrupt and the script falls back to pre-computed CSV stats from an earlier run - that fallback CSV is an undocumented external input.

29. `02_methylation/12_celltype_4layer_master.py` (python) **[!]**  
   Assembles the per-cell-type scRNA breakdown across all seven datasets and the 4-layer (RNA / methylation / proteomics / scRNA) integrated master table. Reads only saved TSVs, no h5ads.
   → MASTER_4layer_validation.tsv, 12_celltype_scRNA_heatmap_CO7.png, 12_celltype_scRNA_heatmap_INVERSE.png
   > Dated 2026-05-17: predates the tier-1 demotion (CD79B/LXN/SH3BP4 -> Tier-2) and the complete-case proteomics switch, so its Proteomics/processed/META reads may now resolve to superseded files.

34. `04_singlecell/build_pb_lognorm.py` (python) **[!]**  
   Streaming pass 1: donor x 8-lineage pseudobulk as the MEAN of the deposited log-normalised Kaufmann values (for limma-trend Track 1).
   → PB_lognorm_matrix.csv, PB_lognorm_coldata.csv
   > Streams the 16.1 GB matrix.mtx (529,680,471 non-zeros) in 25M-row chunks via pandas.read_csv. Long-running (tens of minutes to hours on a laptop); RAM dominated by the chunk plus a 15354 x n_groups float64 accumulator.

39. `04_singlecell/build_pb_counts.py` (python) **[!]**  
   Streaming pass 3: recovers exact raw integer UMI counts from the deposited log-normalised matrix (count = expm1(x) * nCount_RNA / 1e4) and sums them per donor x cluster / donor x lineage / donor (whole PBMC). This is the muscat-st
   → PBC_fine_matrix.csv, PBC_fine_coldata.csv, PBC_coarse_matrix.csv
   > Third full 16.1 GB stream; builds three accumulators simultaneously. Asserts integrality of the reconstruction and aborts if it fails.

44. `04_singlecell/build_pb_norm.py` (python) **[!]**  
   Streaming pass 4: mean and sum of per-cell CP10K per donor x cell type (the S2 and S3 aggregation units).
   → PBN_coarse_meanCP10K.csv, PBN_coarse_sumCP10K.csv, PBN_fine_meanCP10K.csv
   > Fourth full 16.1 GB stream. Writes ~770 MB of CSV (PBN_fine_meanCP10K.csv alone is 290 MB). PBN_donor_meanCP10K.csv carries the IKZF1 result the manuscript headlines.

49. `04_singlecell/build_pb_brain_csf.py` (python) **[!]**  
   Donor-level pseudobulk for the two non-Kaufmann cohorts: Jaekel = SUM of raw integer UMIs per patient x cell type (region blocks collapsed to patient); Beltran = MEAN of the deposited normalised values per twin x cell type (raw co
   → PBC_jakel_matrix.csv, PBC_jakel_coldata.csv, PBN_beltran_meanNorm.csv
   > Densifies the whole Beltran matrix (a.raw.X.toarray()); loads the Jaekel dense text matrix. Multi-GB peak but no streaming.

54. `04_singlecell/pseudobulk_muscat_style.R` (R) **[!]**  
   CANONICAL count-model differential state (S1): summed raw UMIs, paired ~batch_pair + condition, edgeR quasi-likelihood F-test and limma-voomLmFit, over six analyses (coarse-8 / fine-25 / whole-PBMC x all-pairs / naive-only), with 
   → pseudobulk_muscat_style.csv
   > Writes a 134 MB CSV. Reads the 50 MB PBC_fine_matrix.csv with read.csv (slow); estimateDisp/glmQLFit run per cell type x analysis.

59. `04_singlecell/pseudobulk_DA.R` (R) **[!]**  
   Differential abundance (cell-type composition) with a propeller-style arcsine-sqrt transform and the same paired limma design. Source of the manuscript's 'all BH-FDR > 0.19' composition claim.
   → DA_coarse.csv, DA_fine.csv
   > Fast (coldata only). Reads sex/age columns that only PBC_*_coldata.csv carries.

64. `04_singlecell/pseudobulk_norm_compare.R` (R) **[!]**  
   CANONICAL normalisation-based aggregation sensitivity (S2 mean-CP10K, S3 sum-CP10K) on log2 scale with limma-trend + arrayWeights, protein-coding restricted. Produces the whole-PBMC IKZF1 result the manuscript headlines and the S3
   → pseudobulk_norm_compare.csv
   > Reads the 290 MB and 275 MB PBN_fine CSVs with read.csv; peak RAM several GB. Writes a 78 MB CSV. Verified: S2_meanCP10K_wholePBMC has exactly 11,409 genes and IKZF1 logFC -0.0780 / p = 0.003582, matching the manuscript.

69. `04_singlecell/pseudobulk_brain_csf.R` (R) **[!]**  
   Donor-level pseudobulk DE for the Jaekel brain (edgeR-QL + voom on summed counts) and Beltran CSF/PBMC (limma-trend on mean normalised values) cohorts, min 3 donors per arm.
   → pseudobulk_brain_csf.csv
   > Fast. Its 17 MB output is NOT read by make_pseudobulk_excel.py, so the reviewer workbook does not contain the brain/CSF donor-level results.

74. `06_figures/build_pseudobulk_workbook.py` (python) **[!]**  
   Assembles the reviewer-facing donor-level pseudobulk workbook AND computes the panel-scope BH-FDR that the manuscript quotes (IKZF1 panel FDR 0.047). It is an analysis step, not just formatting - the 0.047 figure exists nowhere el
   → Pseudobulk_DonorLevel_Results.xlsx
   > Loads the 134 MB and 78 MB result CSVs into memory simultaneously. Its hardcoded README text still asserts the OLD five-gene tier-1 framing and calls ITGB2's tier-1 anchoring 'not single-cell', which is now inconsistent with the revised tier-1 = ITGB2 + IKZF1.

79. `06_figures/figure5_singlecell.py` (python) **[!]**  
   CANONICAL manuscript Figure 5 generator (composite single-cell panel A/B/C). Holds the CURRENT tier lists: INV_TIER1 = ITGB2 + IKZF1, SUGGESTIVE = HLA-E, TIER2_INV_AUX includes LXN, CD79B, SH3BP4. Edited in place 2026-08-02.
   → singlecell_combined_figure_v3_INV.png, image5.png (md5 verified identical)
   > Loads all three h5ads; several GB peak. This file exists ONLY under MS_GEO_pipeline/scripts/_figure_generators/ - there is no Poster_v2 sibling, so for this one file the suspect snapshot IS the source of truth.

**Notes:** TWO MANUSCRIPT-REPORTED OUTPUTS HAVE NO PRODUCER SCRIPT. (1) Poster_v2/figures/scrna_WILCOXON_v1.tsv is the source of every cell-level number in Section 2.6 and of Figure 5 Panel B (793 rows verified; IKZF1/t_cells FDR 9.389e-71 matches the quoted 9.4e-71). An exhaustive grep of the whole tree finds only readers - make_combined_singlecell_figure_v3_INV.py, _apply_v130_multipletest.py, _cleanup_batch1.py - and no writer. A sibling scrna_WILCOXON_v1.CORRUPT_BACKUP.tsv (2026-06-28 20:24) shows it was regenerated at 20:36 by an ad-hoc, unsaved session. This script must be rewritten to the Methods 

## Integration, network/pathway analysis and manuscript figures (MS_GEO_pipeline/scripts/05_integration_pathway/, Methylation/r_notebooks/09-10, Poster_v2/promoter_vs_body_test.R, MS_GEO_pipeline/scripts/_figure_generators/make_*.py)

**Download first**
- STRING v12 physical-interaction network, queried live at https://string-db.org/api (tsv/network, tsv/interaction_partners, tsv/enrichment; species 9606, network_type=physical) by ppi_analysis.py. No local cache of the API version - the release must ship PPI/ppi_string_physical_full.tsv, ppi_string_partners_top10.tsv, ppi_string_enrichment_*.tsv and ppi_hub_ranking.tsv as frozen artefacts and document the query date (archived run 2026-06-14).
- GeneMANIA physical-interaction network, queried live at https://genemania.org/json/network_data by ppi_analysis.py. Returned zero edges on the archived run (PPI/ppi_genemania_physical.tsv is a 1-byte file) - document this, or the reviewer will assume a failed re-run.
- g:Profiler gOSt API (organism hsapiens; GO:BP, KEGG, REAC, WP; g_SCS correction), queried live by make_string_network_v7_panel38.py. Manuscript pins this to "gOSt API, August 2025"; the API is not versioned locally. Frozen result appears to be PPI/ppi_gprofiler_17.csv - ship it and repoint the script's cache path away from /tmp/string17_gprofiler.csv.
- Enrichr gene-set libraries downloaded live by gseapy in run_ora_inv97.py: GO_Biological_Process_2023, GO_Molecular_Function_2023, GO_Cellular_Component_2023, KEGG_2021_Human, Reactome_2022, WikiPathway_2023_Human, MSigDB_Hallmark_2020, DisGeNET, GWAS_Catalog_2023, Jensen_DISEASES. Library snapshot dates are not recorded anywhere; the archived outputs in ORA_INV97/ are dated 2026-06-11.
- UK Biobank-PPP Olink plasma summary statistics (Jacobs et al. 2024, PMID 38282238) as Proteomics/blood_raw/Jacobs2024_UKB_primary_DE.tsv - published summary statistics, read directly by the Figure 4 and Figure 6 generators. Redistribution status must be checked; document acquisition from the paper's supplement.
- PRIDE PXD064570 (Bader & Mann CSF Orbitrap Astral) and PXD045058 (Bader & Mann CSF timsTOF), and MassIVE MSV000096790 / Figshare (Wang & Julien region-resolved brain WM proteome) - upstream of the RDEP_CC tables that Figures 4 and 6 consume.
- GEO series feeding this layer indirectly via Transcriptome/results and Methylation/results, and directly via the single-cell h5ad inputs to Figure 5: GSE118257 (Jaekel brain snRNA-seq), GSE144744 (Kaufmann PBMC 10x, stored under the misleading blood_Ramesh2020_UMAP/ directory name), GSE127969 (Beltran CSF/PBMC twins).
- Illumina 450K annotation package IlluminaHumanMethylation450kanno.ilmn12.hg19 (v0.6.1) - Bioconductor annotation data consumed by promoter_vs_body_test.R; must be version-pinned because probe-to-RefGene-group mapping determines the promoter/body split.

1. `02_methylation/09_inverse_concordance_scan.R` (R) **[!]**  
   Strict RNA x methylation inverse-concordance discovery scan: 7 RNA strata x 2 methylation sources (combined-cohort Stouffer gene-level from nb05, mCSEA promoter NES from nb06); requires FDR<0.05 in both layers with opposite sign. 
   → INVERSE_CONCORDANT_full_pairings.tsv, INVERSE_CONCORDANT_by_gene.tsv, INVERSE_CONCORDANT_by_gene_by_stratum.tsv
   > Fast (<1 min). MUST be run with cwd = Methylation/r_notebooks because it does source("helpers.R") with a relative path. helpers.R hardcodes PROJ_ROOT=<project root>.

6. `02_methylation/10_inverse_proteomics_validation.R` (R) **[!]**  
   Looks the top-40 inverse-concordant genes up in the R proteomics DE tables (CSF Astral, CSF timsTOF, CSF combined, T-lineage meta, Pegram NK8, 4 Magliozzi/Wang brain contrasts) and builds a gene x assay validation heatmap.
   → INV_proteomics_validation_long.tsv, INV_proteomics_validation_by_gene.tsv, 10_inverse_proteomics_heatmap.png
   > Fast. REVISION-CRITICAL STALENESS: it reads the *_R_DEP_results.tsv / Magliozzi_R_DEP_* tables, i.e. the MinProb-imputed DEP outputs that this revision REJECTED. It has NOT been repointed at Proteomics/processed/META/*_CC_*.tsv. Its outputs on disk are from 2026-06-14 11:43 and a

11. `02_methylation/promoter_vs_body_test.R` (R) **[!]**  
   Reviewer-1 point-2 compartment test. Classifies every 450K probe by UCSC_RefGene_Group into promoter (TSS1500/TSS200/5'UTR/1stExon) vs body (Body/3'UTR) and computes a signed-Stouffer promoter-only and body-only statistic per gene
   → promoter_vs_body_by_gene.csv
   > Fast. Two caveats: (a) it uses the 450K annotation ONLY, so EPIC-exclusive probes in the DMF/ocrelizumab (GSE219293, EPIC) stratum get no compartment assignment even though EPICanno is installed; (b) its console TIER-1 block is hardcoded to the OLD 5-gene tier-1 c("ITGB2","LXN","

16. `05_integration/ppi_analysis.py` (python) **[!]**  
   Builds the STRING physical + GeneMANIA PPI network for the candidate panel (CO7 + top-30 inverse-concordant from nb09 + multi-stratum TX + scRNA-validated + curated MS markers), ranks hubs by degree, runs Louvain communities. Its 
   → ppi_string_physical_full.tsv, ppi_string_partners_top10.tsv, ppi_string_enrichment_all.tsv
   > NON-DETERMINISTIC AND MUST NOT BE RE-RUN BLIND. (1) It queries the live STRING and GeneMANIA web services, so a re-run picks up whatever STRING release is current, not the v12 the manuscript cites. (2) ~90 sequential per-gene REST calls with time.sleep(0.15). (3) GeneMANIA return

21. `05_integration/export_used_omics_inventory.py` (python) **[!]**  
   Builds the per-dataset methylation + RNA-seq inventory workbook (GSE-level sample counts, tissue, platform) - the ancestor of Supplementary Table S1 (Data Sources).
   → Used_Methylation_and_RNAseq_Dataset_Inventory.xlsx
   > On-disk output is 2026-03-15, long before the revision re-verified sample counts against GEO/PRIDE. The manuscript's Supplementary Table S1 is now built by Poster_v2/build_data_sources_table.py, not by this script. Ship as provenance only; do not present it as the source of S1.

26. `06_figures/figure_constants.py` (python) **[!]**  
   Single source of truth for gene panels (INV_TIER1, TIER2_PROT, CO7, CANDIDATE_PANEL), cohort sizes and proteomic dataset metadata. Updated this round: INV_TIER1 = ["ITGB2","IKZF1"].
   → (importable constants only)
   > EXISTS ONLY under _figure_generators/ - there is no copy in Poster_v2/, so the Poster_v2 duplicates of make_extra_figures.py / make_figs4to7_v3.py / make_proteomics_v4.py / make_proteomics_v5_with_blood.py / make_intersection_heatmap_v4.py / make_workflow_v2.py cannot even import

31. `06_figures/figure1_workflow.py` (python) **[!]**  
   FIGURE 1 (figures/image1.png) - the four-track workflow schematic. Pure matplotlib drawing, no data input; all counts and gene names are literal text in the script.
   → image1.png (md5 d1df1ae8aa7d1a297d2c73b7e7af3624, verified identical)
   > CURRENT: Figure 1 reports Tier-1 (2: ITGB2, IKZF1), Tier-2 auxiliary (10, including CD79B), "both tier-1 genes recovered" for the ComBat matrix, and the corrected 79-donor single-cell total.

36. `06_figures/figure2_rna_volcanoes.py` (python) **[!]**  
   FIGURE 2 (figures/image2.png) - seven per-stratum RNA volcano panels (PBMC, T cells, IFN-b PBMC, B cells, whole blood, brain WM, pan-tissue) with tier-coloured candidate labels.
   → per_celltype_volcanoes_INV.png -> image2.png (md5 d1a0a7af5dcd44e98fca86c60b58ca00, verified)
   > Tier lists are CURRENT (INV_TIER1={ITGB2,IKZF1}; CD79B in TIER2_AUX_INV). Resolved: panel G now reads Transcriptome/results/07_pan_tissue_DE.tsv. (Historically it read Poster_v2/figures/COMBINED_pantissue_proper_DEG.csv (1.3 MB, 2026-06-15) and a project-wide grep finds NO script that writes it - the only files mentioning it are this generator 

41. `06_figures/figure3_methylation.py` (python) **[!]**  
   FIGURE 3 (figures/image3.png) - gene-level methylation volcanoes per stratum (A-D), combined IDAT/ComBat analysis (E) and the inverse-concordant discovery pool (F).
   → methylation_v_INV.png -> image3.png (md5 5379c4cf72a35d793eb5c58cbf4bda97, verified)
   > Tier code and documentation are CURRENT (INV_TIER1={ITGB2,IKZF1}; CD79B in the 10-gene TIER2_AUX_INV set). It reads the *_meth_gene.tsv tables that were regenerated on 2026-08-02 02:12 by the new helpers.R; I verified mean_logFC and adj.P

46. `06_figures/figure4_proteomics.py` (python) **[!]**  
   FIGURE 4 (figures/image4.png) - CSF Astral/timsTOF volcanoes, UK Biobank-PPP plasma, four region-resolved brain contrasts, and the panel-H directional-consistency heatmap across seven compartments.
   → proteomics_v_INV.png -> image4.png (md5 4539934e269dff222562724e1643aea7, verified)
   > CURRENT: repointed to RDEP_CC and INV_TIER1=["ITGB2","IKZF1"]. It still loads Unified_All_Assays_Long.tsv at line 71 into prot_u, but I traced every one of the seven panel-H columns to an explicit RDEP_CC branch - the prot_u fallback at line 185 is unreachable dead code. The Unif

51. `06_figures/figure5_singlecell.py` (python) **[!]**  
   FIGURE 5 (figures/image5.png) - single-cell validation: cohort UMAPs for Jaekel brain, Kaufmann PBMC and Beltran CSF/PBMC twins, plus the tie-corrected Wilcoxon log2FC heatmap with BH-FDR asterisks.
   → singlecell_combined_figure_v3_INV.png -> image5.png (md5 1b377e4ab585d9a84a52ad096fc7bf75, verified)
   > Loads ~170 MB of h5ad; slowest generator. Tier code CURRENT; docstring lines 10-11 still list the old 6-gene tier-1 - cosmetic. NAMING TRAP: cohort 2 is correctly labelled "Kaufmann 2021 / GSE144744" in the figure but its `short` key is "Ramesh_2020_PBMC" and its path is blood_Ra

56. `06_figures/figure6_intersection_heatmap.py` (python) **[!]**  
   FIGURE 6 (figures/image6.png) - the 17-gene x 22-readout cross-modal intersection matrix, blocked into inverse-concordant Tier-1, Tier-2 auxiliary and Tier-2 non-concordant proteomic anchors.
   → intersection_heatmap_v_INV.png -> image6.png (md5 5907b04af867285a756a0e3986708a70, verified)
   > CURRENT: INV_TIER1=["ITGB2","IKZF1"], CD79B moved into TIER2_AUX_INV, all seven proteomic columns routed through PROT_RAW (RDEP_CC). I verified the Unified_All_Assays_Long.tsv fallback branch is used ONLY for the 7 RNA and 4 methylation columns, never for proteomics - so the stal

61. `06_figures/figure7_string_network.py` (python) **[!]**  
   FIGURE 7 (figures/image7.png) - panel A: STRING v12 physical network of the 38-gene panel (17 named candidates + connected canonical MS markers, isolated named candidates drawn below the rule); panel B: g:Profiler g:SCS over-repre
   → string_network_v6.png -> image7.png (md5 aa23a62c199b9cd05b2ff272971577a9, verified)
   > THIS, NOT make_string_network_v6.py, IS THE FIGURE-7 GENERATOR - it deliberately writes to string_network_v6.png ("overwrite the embedded Fig 7", line 13). I re-executed its panel arithmetic against the archived PPI tables and reproduced the caption exactly: 58 hub genes -> 31 co

**Notes:** CANONICAL SOURCE. For every figure generator, _figure_generators/ is newer than the Poster_v2/ sibling and must be used. I diffed all 44 make_*.py pairs: 27 identical, 15 differ (in every differing case _figure_generators is newer - e.g. make_intersection_heatmap_v_INV.py 2026-08-02 vs 2026-06-03, make_proteomics_v_INV.py 2026-08-02 vs 2026-06-03), and 3 exist only in _figure_generators (make_combined_singlecell_figure_v3_INV.py, make_string_network_v7_panel38.py, _figure_constants.py). Since _figure_constants.py has no Poster_v2 copy, six Poster_v2 duplicates cannot even import. For the five 
