# Data sources

No data is redistributed. Download each accession from its repository before running the
corresponding layer. Sample counts are the post-QC numbers used in the paper.

| layer | accession | cohort | tissue | platform | MS | HC |
|---|---|---|---|---|---|---|
| Bulk transcriptomics | GSE103005 | Yang et al. 2022 | Whole blood | HT-12 v4 (Array) | 8 | 12 |
| Bulk transcriptomics | GSE138064 | Hagan et al. 2020 | PBMC + WB | HTA-2.0 (Array) | 64 | 8 |
| Bulk transcriptomics | GSE172009 | Unpublished | CD4 T | NovaSeq (RNA-seq) | 4 | 4 |
| Bulk transcriptomics | GSE173789 | Aktas et al. 2022 | B cells | HiSeq X (RNA-seq) | 23 | 14 |
| Bulk transcriptomics | GSE190847 | Kular et al. 2023 | B cells | Clariom-D (Array) | 93 | 28 |
| Bulk transcriptomics | GSE207680 | Schirmer et al. 2023 | Brain / WM | NextSeq (RNA-seq) | 2 | 3 |
| Bulk transcriptomics | GSE209596 | Unpublished | T cells | HiSeq (RNA-seq) | 37 | 40 |
| Bulk transcriptomics | GSE211358 | Glanzman et al. 2022 | B cells | NovaSeq (RNA-seq) | 8 | 7 |
| Bulk transcriptomics | GSE214334 | Unpublished | Brain / WM | NovaSeq (RNA-seq) | 11 | 7 |
| Bulk transcriptomics | GSE21942 | Kemppinen 2011 | PBMC + WB | U133+2 (Array) | 14 | 15 |
| Bulk transcriptomics | GSE288904 | 2025 release | Whole blood | NovaSeq (RNA-seq) | 8 | 25 |
| Bulk transcriptomics | GSE38010 | Han et al. 2012 | Brain / WM | U133+2 (Array) | 5 | 2 |
| Bulk transcriptomics | GSE43591 | Ottoboni et al. 2013 | Whole blood | U133+2 (Array) | 10 | 10 |
| Bulk transcriptomics | GSE66573 | Spurlock et al. 2015 | Whole blood | RNA-seq (HiSeq 2500) | 6 | 8 |
| DNA methylation | GSE106648 | Kular et al. 2018 | Whole blood | 450K | 140 | 139 |
| DNA methylation | GSE130029 | Maltby et al. 2020 | T cells (CD4) | 450K | 12 | 11 |
| DNA methylation | GSE130030 | Maltby et al. 2020 | T cells (CD4/CD8) | 450K | 10 | 14 |
| DNA methylation | GSE189255 | Souren et al. 2022 | T cells (CD4) — Remission | 450K | 6 | 8 |
| DNA methylation | GSE189256 | Souren et al. 2022 | Monocytes (CD14) | 450K | 6 | 8 |
| DNA methylation | GSE219293 | Bos et al. 2023 | Whole blood — DMF / Ocrelizumab | EPIC | 29 | 18 |
| DNA methylation | GSE40360 | Huynh et al. 2014 | Brain / WM | 450K | 28 | 19 |
| DNA methylation | GSE88824 | Kular et al. 2018 (sub) | T cells / WB | 450K | 13 | 14 |
| DNA methylation | GSE173787 | Ma et al. 2021 | Sorted immune cells (CD19 B, CD4, CD8 T, CD14 mono) | WGBS (sequencing) | 72 | 57 |
| Single-cell RNA-seq | GSE118257 | Jäkel et al. 2019 | Brain / WM (snRNA-seq) | 10x Chromium | — | — |
| Single-cell RNA-seq | GSE144744 | Kaufmann et al. 2021 | PBMC | 10x Chromium 3′ (scRNA-seq + CITE-seq) | — | — |
| Single-cell RNA-seq | GSE127969 | Beltrán et al. 2019 | CSF + PBMC (twin pairs) | 10x Chromium (5′) | 5 | 5 |
| Proteomics (CSF Astral) | PXD064570 | Bader, Mann 2026 | CSF | Orbitrap Astral DIA | 978 | 306 |
| Proteomics (CSF timsTOF) | PXD045058 | Bader, Mann 2026 | CSF | Bruker timsTOF dia-PASEF | 1536 | 2363 |
| Proteomics (brain WM) | MSV000096790 (MassIVE) | Wang, Julien 2026 | Brain (cortex, NAWM, lesion, lesion-edge) | Bruker timsTOF dia-PASEF | 20-32 | 12-20 |
| Proteomics (UKB plasma) | UK Biobank-PPP (Olink Explore) | Jacobs et al. 2024 | Plasma | Olink Explore (antibody) | 407 | 39979 |

## Repositories

- **NCBI GEO** — all `GSE*` accessions: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSExxxxx
- **PRIDE** — `PXD064570` (Bader & Mann CSF Orbitrap Astral), `PXD045058` (same study, timsTOF CSF)
- **MassIVE** — `MSV000096790` (Wang & Julien region-resolved brain white-matter proteome; processed tables on Figshare)
- **UK Biobank-PPP** — published Olink Explore association statistics (Jacobs et al. 2024, PMID 38282238);
  individual-level data require approval via the UK Biobank Access Management System.

## Notes

- Two methylation series, **GSE106648** and **GSE40360**, deposit no IDATs; they enter from the
  depositor-processed matrices and therefore bypass the detection-p probe filter applied to the
  other six array series. See Methods 4.3.
- **Two series are downloaded and harmonised but enter no reported analysis**, and are therefore
  not among the 29 datasets above. They are fetched by `00_data/...download_bulk_rnaseq.py`
  because `harmonize_rnaseq_v3.py` processes them upstream of the point at which they are
  dropped, so the pipeline will not run without them:
  - `GSE137143` — 81 samples, 77 MS / 4 HC. Four controls cannot support a case-control
    contrast; under ComBat the series is a batch↔condition confound. Removed from the study
    during revision.
  - `GSE211739` — 10 samples, 4 MS / 6 HC, iPSC-derived oligodendrocytes. Present in the
    472-sample batch-corrected discovery matrix but in no analysis stratum.
