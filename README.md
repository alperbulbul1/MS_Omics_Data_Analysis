# MS multi-omics re-analysis — analysis code

Analysis code for the study prioritising *ITGB2* and *IKZF1* in multiple sclerosis by
inverse-concordant DNA-methylation × transcription integration across four omic layers.

Every script ships with the placeholders `__MS_GEO_ROOT__` (project root) and `__PYTHON_BIN__`
(interpreter); `./configure.sh /path/to/data /path/to/python` substitutes both in place, and
`./configure.sh --check` fails if any placeholder or author-specific path survives.

---

## `scripts/00_data/` — acquisition

| file | what it does |
|---|---|
| `download_bulk_rnaseq.py` | Downloads the bulk-transcriptomic GEO series (727 usable expression samples). |
| `download_methylation.py` | Downloads the methylation series listed in `Methylation_Target_Datasets.csv`, resuming partial files. |
| `select_methylation_case_control.py` | Screens the candidate methylation studies for MS-vs-control designs and writes the target list. |

## `scripts/01_transcriptome/` — bulk RNA differential expression

| file | what it does |
|---|---|
| `harmonize_rnaseq_v3.py` | Harmonises the RNA-seq series: per-dataset gene-ID → symbol mapping and metadata-based MS/HC assignment. |
| `harmonize_microarray_v2.py` | Harmonises the array series: platform probe → HGNC symbol via GEO platform SOFT annotation and mygene, `max` per symbol. |
| `merge_harmonized_v2.py` | Merges the RNA-seq and microarray tracks into one symbol × sample matrix. |
| `build_matrices.py` | Builds the per-series expression matrices with Ensembl → symbol resolution. |
| `build_global_matrices.py` | Builds the combined global matrix, parsing SOFT headers for MS/HC labels without loading the full tables. |
| `merge_and_correct.py` | Merges the recovered supplementary datasets into the global matrix and re-runs neuroComBat and limma. |
| `correct_and_normalize.py` | Per-dataset normalisation and ComBat batch correction with disease status protected; prefers the repaired matrix and raises on any all-zero dataset. |
| `rerun_verified_case_control.py` | Re-runs the case-control analyses against the verified dataset inventory. |
| `00_run_all.R` | Runs every transcriptome step in order. |
| `01_pbmc_de.R` … `05_whole_blood_de.R` | Per-stratum limma-trend differential expression: PBMC, T cells, B cells, brain white matter, whole blood. |
| `06_pbmc_ifnb_de.R` | IFN-β-versus-baseline PBMC contrast. This is a **treatment-response** contrast and is deliberately excluded from the inverse-concordance disease scan. |
| `07_total_combined_de.R` | Pan-tissue combined model (`~ tissue + condition`) over the pooled case-control strata. |
| `08_cross_stratum_master.R` | Cross-stratum summary table. |
| `run_expression_subgroup_limma.R` | Generic subgroup limma runner (`meta.csv matrix.csv out_dir [precorrected]`). |
| `run_stratified_omics.py` | Drives the per-stratum runs across both layers, preferring the IDAT-preprocessed methylation input and falling back to the combined matrix. |
| `infer_sex_from_expression.py` | Infers sample sex from XIST and Y-linked gene expression where the deposit does not state it. |
| `sex_adjusted_sensitivity_rna.R` | Sex-adjusted sensitivity re-run of the RNA layer for the candidate panel. |
| `helpers.R` | Shared paths and utilities for the transcriptome scripts. |

## `scripts/02_methylation/` — DNA methylation

| file | what it does |
|---|---|
| `preprocess_methylation_arrays.R` | minfi preprocessing of the IDAT-based array series. |
| `normalize_beta_only.R` | Normalisation path for the series deposited as beta/signal matrices rather than IDATs. |
| `build_methylation_matrix.py` | Assembles the per-series matrices onto common probes. |
| `run_all_methylation_combat.R` | ComBat batch correction across all methylation series with disease status protected. |
| `00_run_all.R` | Runs every methylation step in order. |
| `01_tcells_meth_dmp.R` … `05_combined_meth_dmp.R` | Per-stratum differentially methylated positions: T cells, whole-blood dimethyl fumarate, whole-blood ocrelizumab, T-cell remission, combined. |
| `06_mcsea_promoter_analysis.R`, `run_mcsea_combat.R` | mCSEA promoter and gene-body region enrichment on the ComBat-corrected M-values. |
| `07_brainwm_rna_vs_meth.R` | Brain white-matter RNA-versus-methylation concordance. |
| `08_cross_stratum_meth_master.R` | Cross-stratum methylation summary. |
| **`09_inverse_concordance_scan.R`** | **The core filter.** Scans for genes significant in both layers with opposite directions across matched tissue strata. The IFN-β PBMC contrast is excluded here by design. |
| `10_inverse_proteomics_validation.R` | Proteomic anchoring of the inverse-concordant pool. |
| `11_inverse_scRNA_validation.py` | Single-cell anchoring of the inverse-concordant pool. |
| `12_celltype_4layer_master.py` | Four-layer × cell-type master matrix. |
| `13_perstudy_scRNA_validation.py` | Per-study single-cell breakdown. |
| `15_genelevel_weighting_corrected.R` | Gene-level signed-Stouffer aggregation with the unweighted probe mean as the reported effect, plus the 1/SE-weighted and inverse-variance alternatives. |
| `promoter_vs_body_test.R` | Direct promoter-versus-gene-body compartment test, for genes where mCSEA has too few CpGs to form a testable region. |
| **`rebuild_S2_94genes.py`** | Rebuilds Supplementary Table S2 for the current 94-gene pool in one pass: per-gene statistics, methylation-compartment flag (promoter-confirmed / promoter-only / composite-only) and the promoter-only sensitivity analysis. **Supersedes `update_S2_promoter_flag.py` and `update_S2_weighting.py`, which patched an older 82-gene sheet.** |
| `run_methylation_subgroup_limma.R` | Generic subgroup limma runner for the methylation layer. |
| `infer_sex_GSE88824_chrXY.R` | Infers sex for GSE88824, the one array series without deposited sex, from chrX/chrY signal. |
| `sex_adjusted_sensitivity.R` | Sex-adjusted sensitivity re-run of the methylation layer. |
| `helpers.R` | Shared paths and utilities for the methylation scripts. |

## `scripts/03_proteomics/` — CSF and brain differential abundance

| file | what it does |
|---|---|
| `00_run_all.R` | Runs the proteomic steps in order; invokes the complete-case scripts, which are the reported analysis. |
| **`01cc_csf_astral_completecase.R`** | CSF Orbitrap Astral DIA-MS differential abundance. Tested with DEP **without imputation**, with DEP's `fdrtool` q-values replaced by Benjamini–Hochberg. |
| **`02cc_csf_timstof_completecase.R`** | The same for the timsTOF CSF platform. |
| `04cc_magliozzi_brain_completecase.R` | Brain white-matter proteome, complete-case. |
| `dep_bh_equivalence_check.R` | Shows that the reported result is equivalent to running DEP itself at the level of effect estimates, and that the only material difference is `fdrtool` versus BH. |
| `03_csf_cross_platform_meta.R` | Meta-analysis across the two CSF platforms. |
| `05_t_lineage_meta.R` | T-lineage proteomic meta-analysis. |
| `06_pegram_gse32915_de.R` | GSE32915 differential expression. |
| `07_brainwm_rna_meth_rerun.R` | Brain white-matter RNA/methylation re-run. |
| `08_per_group_consistency.R` | Per-group consistency checks. |
| `09_cross_assay_lxn.R` | Cross-assay *LXN* comparison. |
| `10_master_validation.R` | Master proteomic validation table. |
| `11_itgb2_csf_pleocytosis.R` | Tests why *ITGB2* is detected more often in MS than in control CSF, and shows this tracks leukocyte count rather than an MS-specific property of the protein. |
| `build_RDEP_CC_adapters.py` | Builds the figure-adapter tables from the complete-case proteomic output. |
| `helpers.R` | Shared paths and utilities for the proteomic scripts. |

## `scripts/04_singlecell/` — single-cell processing and donor-level pseudobulk

| file | what it does |
|---|---|
| `download_singlecell_datasets.py` | Downloads the single-cell h5ad files from CELLxGENE Discover. |
| `process_GSE118257_jakel_brain.py` | Jäkel et al. 2019 MS brain snRNA-seq, using the deposited final clustering. |
| `process_GSE127969_beltran_csf.py` | Beltrán et al. 2019 CSF and PBMC scRNA-seq from monozygotic twins discordant for MS. |
| `process_GSE144744_kaufmann_pbmc.py` | Kaufmann et al. PBMC and CSF scRNA-seq: 62 donors, 497,706 cells. |
| `process_GSE144744_kaufmann_genes.py` | Companion pass over the full count matrix for per-gene extraction. |
| `plot_GSE144744_kaufmann_celltypes.py` | Re-renders the cell-type UMAPs from the cached object. |
| `build_pb_counts.py` | Donor × cell-type pseudobulk as the **sum of raw integer UMI counts** (the muscat standard). |
| `build_pb_lognorm.py` | Pseudobulk from the deposited log-normalised matrix. |
| `build_pb_norm.py` | Normalisation-based pseudobulk, as a sensitivity comparison. |
| `build_pb_brain_csf.py` | Pseudobulk inputs for the brain and CSF/blood cohorts on the same footing. |
| `pseudobulk_muscat_style.R` | Donor-level differential-state analysis following the muscat/Squair decision rules. |
| `pseudobulk_brain_csf.R` | The same for the brain and CSF/blood cohorts. |
| `pseudobulk_norm_compare.R` | Aggregation-sensitivity comparison: does the pseudobulk unit change the conclusions? |
| `pseudobulk_DA.R` | Differential **abundance** (cell-type proportions), analysed separately from differential state. |
| `run_pseudobulk_reanalysis.py` | Donor-level MS-versus-HC re-analysis for the Tier-1 candidates, addressing cell-level p-value inflation. |

## `scripts/05_integration/` — network, pathway and inventory

| file | what it does |
|---|---|
| `ppi_analysis.py` | STRING physical-only protein–protein interaction query plus GeneMANIA, for the candidate panel. |
| `build_data_sources_table.py` | Generates Supplementary Table S1, the per-study data-source catalogue. |
| `export_used_omics_inventory.py` | Exports the dataset inventory actually used by the methylation and expression layers. |

## `scripts/06_figures/` — one script per manuscript figure

| file | what it does |
|---|---|
| `figure_constants.py` | Single source of truth for the figures: canonical gene panels, tier membership and dataset metadata. |
| `figure1_workflow.py` | Figure 1, the pipeline schematic. |
| `figure2_rna_volcanoes.py` | Figure 2, per-stratum RNA volcanoes plus the tissue-adjusted pan-tissue panel. |
| `figure3_methylation.py` | Figure 3, gene-level methylation panels. |
| `figure4_proteomics.py` | Figure 4, proteomic validation across the CSF, brain and plasma compartments. |
| `figure5_singlecell.py` | Figure 5, single-cell assessment across the three cohorts. |
| `figure6_intersection_heatmap.py` | Figure 6, the four-layer × tissue intersection matrix. |
| `figure7_string_network.py` | Figure 7, the STRING physical-interaction network and pathway enrichment. |
| `build_pseudobulk_workbook.py` | Assembles the donor-level pseudobulk results into one Excel workbook. |

---

## Supporting files

| file | what it does |
|---|---|
| `configure.sh` | Substitutes `__MS_GEO_ROOT__` and `__PYTHON_BIN__` across `scripts/`; `--check` verifies no placeholder or author path remains. |
| `docs/RUN_ORDER.md` | Dependency order with each step's inputs and outputs, and which steps are not safe to re-run blindly. |
| `docs/DATA.md` | Every accession with the layer it feeds. No data is redistributed. |
| `env/r-packages.txt`, `env/requirements.txt` | Package versions captured from the session that produced the reported results. |
