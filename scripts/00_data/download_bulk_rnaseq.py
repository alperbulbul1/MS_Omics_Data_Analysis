#!/usr/bin/env python3
"""Consolidated, reproducible download of the bulk-transcriptomic series.
727 usable expression samples across the 16 series listed below. Only 14 are declared in the
paper's Data Availability: GSE137143 and GSE211739 are processed upstream but enter no reported
analysis (see the notes on each) (the GEO catalog lists 1,076, but GSE211358's 32 are spread over 3 sub-experiment files →
15 in the MS-vs-HD matrix). 472 are used in the treatment-naive MS-vs-HC discovery; the remainder are
treatment-experienced arms (e.g. GSE138064 IFN-β time-course) analysed as separate per-stratum comparisons (Expression_Data/Combined_Metadata.csv).
Each entry is (n_usable_expression_samples, n_used_in_treatment_naive_discovery, tissue/first-author).

This supersedes the ad-hoc 6-GSE batch scripts: it enumerates the series
cohort explicitly so the RNA layer is reproducible from a single script.
For each GSE it pulls the depositor-provided series_matrix + processed supplementary tables
from the NCBI GEO FTP mirror (raw FASTQ/CEL RAW.tar tarballs are intentionally skipped).
"""
import os, re, urllib.parse, requests
from bs4 import BeautifulSoup

DEST = "__MS_GEO_ROOT__/Expression_Data"

# ── The definitive 15 bulk-RNA series (GSE : (n_deposited, n_used, tissue/author)) ──
RNA_SERIES = {
    "GSE190847": (121, 121, "B cells (Kular et al. 2023)"),
    # Downloaded and harmonised, but excluded from every reported analysis: 77 MS / 4 HC gives a
    # batch<->condition confound under ComBat, so it enters no stratum and is not in the paper.
    # Retained here because harmonize_rnaseq_v3.py processes it upstream of that exclusion.
    "GSE137143": (81, 0, "CD4/CD8 T + CD14 monocytes (Mexhitaj et al. 2019) [EXCLUDED from analyses]"),
    "GSE209596": (77, 77, "blood mTreg/mTeff (unpublished)"),
    "GSE138064": (227, 72, "PBMC, IFN-β cohort (Hagan et al. 2020)"),
    "GSE173789": (37, 37, "B cells (Aktas et al. 2022)"),
    "GSE288904": (37, 33, "neutrophils (2025 release)"),
    "GSE21942": (29, 29, "PBMC (Kemppinen et al. 2011)"),
    "GSE43591": (20, 20, "T cells (Ottoboni et al. 2013)"),
    "GSE103005": (20, 20, "whole blood (Yang et al. 2022)"),
    "GSE214334": (18, 18, "brain NAWM (unpublished)"),
    "GSE211358": (15, 15, "B cells (Glanzman et al. 2022)"),
    "GSE211739": (10, 10, "iPSC-derived oligodendrocytes (unpublished)"),
    # NOTE: GSE211739 is present in the batch-corrected discovery matrix (472 samples) but
    # appears in no analysis stratum and is not listed in the paper's Data Availability.
    # It is retained here only so that the 472-sample matrix can be rebuilt exactly.
    "GSE172009": (8, 8, "CD4+ T cells (unpublished)"),
    "GSE38010": (7, 7, "brain white matter (Han et al. 2012)"),
    "GSE207680": (6, 5, "cortex (Schirmer et al. 2023)"),
    "GSE66573": (14, 14, "whole blood, RRMS (Tetreault/Pawlowski 2015; FPKM)"),
}
# Candidate GSEs screened OUT during curation (NOT in the final 15) — kept for provenance:
EXCLUDED = ["GSE235357", "GSE247181", "GSE255952", "GSE130478"]

def download(url, out):
    if os.path.exists(out):
        print(f"  · exists, skip: {os.path.basename(out)}"); return True
    try:
        r = requests.get(url, stream=True, timeout=90); r.raise_for_status()
        with open(out, "wb") as f:
            for ch in r.iter_content(1024 * 1024): f.write(ch)
        print(f"  ✓ {os.path.basename(out)}"); return True
    except Exception as e:
        print(f"  ✗ {url}: {e}")
        if os.path.exists(out): os.remove(out)
        return False

def fetch(gse):
    nnn = gse[:-3] + "nnn"
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{nnn}/{gse}"
    # 1) series matrix
    download(f"{base}/matrix/{gse}_series_matrix.txt.gz",
             os.path.join(DEST, f"{gse}_series_matrix.txt.gz"))
    # 2) processed supplementary tables (skip RAW.tar)
    try:
        r = requests.get(f"{base}/suppl/", timeout=30)
        if r.status_code == 200:
            for a in BeautifulSoup(r.text, "html.parser").find_all("a"):
                h = a.get("href")
                if (h and gse in h and "RAW.tar" not in h
                        and h.endswith((".gz", ".txt", ".csv", ".tsv", ".xlsx"))):
                    download(urllib.parse.urljoin(f"{base}/suppl/", h), os.path.join(DEST, h))
    except Exception as e:
        print(f"  suppl list failed: {e}")

if __name__ == "__main__":
    os.makedirs(DEST, exist_ok=True)
    print(f"Downloading the definitive {len(RNA_SERIES)} bulk-RNA series "
          f"({sum(v[0] for v in RNA_SERIES.values())} samples deposited; "
          f"{sum(v[1] for v in RNA_SERIES.values())} used in analysis) → {DEST}\n")
    for gse, (dep, use, note) in RNA_SERIES.items():
        print(f"[{gse}]  GEO={dep}  used={use}  {note}")
        fetch(gse)
    print(f"\nDone. Excluded candidates (NOT in the 15): {EXCLUDED}")
