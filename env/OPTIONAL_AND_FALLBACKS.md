# Packages a script imports that are NOT present in the recorded environment.
# For each, what actually executed when the reported results were produced.

DEP (R/Bioconductor)
    -> NOT INSTALLED. Referenced in Proteomics/r_notebooks/01_,02_,01cc_,02cc_ and Proteomics/r_scripts/01_,03_. In the legacy 01_/02_ scripts the guard is `have_DEP <- requireNamespace("DEP", quietly=TRUE)` -> FALSE. In the CANONICAL complete-case scripts 01cc_csf_astral_completecase.R:50 and 02cc_csf_timstof_completecase.R:39 it is hard-coded `have_DEP <- FALSE  # COMPLETE-CASE: force the limma path`, so the `library(DEP)` branch is unreachable dead code. What actually produced every reported proteomics number: Proteomics/r_notebooks/helpers.R -> filter_missval_R (>=50% valid per condition, mirrors DEP::filter_proteins) + vsn_with_fallback (vsn::justvsn 3.80.0, median-centre if the matrix is already log-scale) + moderated_t_safe (limma 3.68.3 lmFit/eBayes), plus a >=2-observed-values-per-group estimability filter, NO imputation. DEP must NOT be listed as a runtime dependency of the release.

wateRmelon (R/Bioconductor)
    -> NOT INSTALLED. normalize_beta_only.R:83 `has_bmiq <- requireNamespace("wateRmelon", quietly=TRUE)` -> FALSE, printing 'wateRmelon not found - using quantile normalization fallback'. BMIQ Type-I/Type-II probe-bias correction therefore NEVER ran; the deposited normalised betas are quantile-normalised via preprocessCore 1.74.0. This is what the manuscript's methylation numbers rest on.

maxprobes (R, GitHub-only)
    -> NOT INSTALLED. preprocess_methylation_arrays.R:26 `requireNamespace("maxprobes")` -> FALSE, so load_xreactive_probes() falls through to the SNP-in-probe proxy: rownames of getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) with non-empty Probe_rs. The cross-reactive-probe exclusion actually applied is the 450K SNP-in-probe list, not the Chen/McCartney cross-reactive list.

openxlsx (R/CRAN)
    -> NOT INSTALLED. Only imported by Proteomics/r_notebooks/04_magliozzi_brain_dep.R:26 and Proteomics/r_scripts/03_brain_proteomics_DEP.R:34 - both SUPERSEDED by 04cc_magliozzi_brain_completecase.R, which reads the same Magliozzi table without openxlsx. The canonical brain outputs (processed/META/Magliozzi_CC_*.tsv, 2026-08-01) came from 04cc_. Ship 04cc_ and openxlsx is not needed; if the legacy 04_ is shipped for provenance it will error on load.

DMRcate (R/Bioconductor)
    -> NOT INSTALLED. `library(DMRcate)` at methylation_analysis_pipeline.R:37 and run_methylation_subgroup_limma.R:36 - an UNGUARDED library() call, so those two legacy monolithic scripts abort at startup on this machine. They are not the canonical methylation path: the reported DMP/mCSEA/gene-level results come from Methylation/r_notebooks/01-10 + 15_genelevel_weighting_corrected.R, none of which touch DMRcate. No DMR-level result is reported in the manuscript. Exclude both legacy scripts from the release (or ship them under a clearly-marked legacy/ with a note).

missMethyl (R/Bioconductor)
    -> NOT INSTALLED. Guarded: methylation_analysis_pipeline.R:616 and run_methylation_subgroup_limma.R:464 do `if (!requireNamespace("missMethyl")) stop("missMethyl not installed")`, i.e. the gometh/GO-of-DMP block hard-stops. No missMethyl-derived enrichment appears in the manuscript; pathway enrichment is g:Profiler (gprofiler-official 1.0.0 / gOSt API) instead.

EnhancedVolcano (R/Bioconductor)
    -> NOT INSTALLED. Unguarded library() in Proteomics/r_scripts/01_,02_,03_ only - the pre-r_notebooks generation of proteomics scripts, all superseded. Every published volcano is drawn by ggplot2 4.0.3 (helpers.R::dep_volcano_gg) or by matplotlib 3.10.8 in scripts/_figure_generators/. Not a release dependency.

python-louvain / `community` (Python)
    -> NOT INSTALLED. `import community as community_louvain` at MS_GEO_pipeline/scripts/05_integration_pathway/ppi_analysis.py:271 and _figure_generators/make_pathway_network_v2.py:167, v3.py:177 - each wrapped in try/except ImportError with fallback `networkx.community.greedy_modularity_communities`. The PPI communities and the Figure 7 network hulls as published are GREEDY-MODULARITY partitions, not Louvain. If the manuscript or a figure caption says 'Louvain', that is a mismatch; the release should either pin python-louvain or correct the wording to greedy modularity.

combat / pyComBat (Python)
    -> NOT INSTALLED and NO fallback. `from combat.pycombat import pycombat` at MS_GEO_pipeline/scripts/04_proteomics/Proteomics__scripts/10_meta_analysis_combine_batch_correct.py:362 raises ImportError, so that script cannot complete its ComBat step. Its output Proteomics/processed/META/T_lineage_combined_DE.tsv (2026-05-02) was SUPERSEDED on 2026-05-16 by T_lineage_R_combined_DE.tsv from Proteomics/r_notebooks/05_t_lineage_meta.R, which uses R sva::ComBat 3.60.0. Exclude script 10 from the release as superseded rather than pinning pyComBat.

ComplexHeatmap (R/Bioconductor)
    -> NOT INSTALLED, and NOT imported by any analysis script in the canonical set - the grep over all 134 .R files returns zero library()/requireNamespace() hits. It appears only in the environment audit, not in code. No action needed; do not put it in the manifest.

MSstats, msigdbr (R/Bioconductor)
    -> NOT INSTALLED. They appear only in the wish-list installer Proteomics/r_scripts/00_install_packages.R (which also claims 'Tested on R 4.3+/4.4+ with Bioconductor 3.18+/3.19+'). No analysis script imports them. That installer is aspirational, not a record of the environment - it must not be shipped as the release's dependency spec.
