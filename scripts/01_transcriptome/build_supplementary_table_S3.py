#!/usr/bin/env python3
"""Assemble Supplementary Table S3 (per-series sex and age composition), Methods 4.2.

WHY THIS EXISTS. Methods cites Supplementary Table S3 for the sex and age composition of every
cohort, but no script assembled it: the delivered CSV was written by hand, so its per-series
percentages could not be checked against the per-sample calls they summarise.

WHAT IS DERIVED AND WHAT IS NOT. The bulk-RNA and methylation rows are computed here from the
per-sample tables that the sex pipeline produces:

    Methylation_Data/rna_sex_persample.csv   <- 01_transcriptome/build_rna_sex_persample.py
    Methylation_Data/meth_sex.csv            <- 02_methylation/build_meth_sex.py

so every percentage in those rows is recomputed rather than transcribed. The single-cell,
proteomic and UK Biobank rows are NOT derivable from anything in this repository: their sex and age
composition comes from the source publications and the depositors' phenotype tables. Those rows are
carried below as declared constants with the source named in the row itself, so a reader can see
which numbers were computed and which were quoted.

Age is likewise only available where the depositor recorded it; "not deposited" is written rather
than an imputed value.
"""
import os
import sys

import pandas as pd

ROOT = "__MS_GEO_ROOT__"
RNA_SEX = os.path.join(ROOT, "Methylation_Data", "rna_sex_persample.csv")
METH_SEX = os.path.join(ROOT, "Methylation_Data", "meth_sex.csv")
EXPR_META = os.path.join(ROOT, "Expression_Data", "Combined_Metadata.csv")
METH_META = os.path.join(ROOT, "Methylation_Data", "AllMeth_ComBat_Metadata.csv")
OUT = os.path.join(ROOT, "Poster_v2", "IJMS_submission_latex_v2",
                   "Supplementary_Table_S3_sex_age.csv")

# Rows that cannot be derived from this repository. Each carries its own provenance.
QUOTED = [
    dict(Layer="Single cell", Series="GSE118257 (Jaekel brain)", Source="publication"),
    dict(Layer="Single cell", Series="GSE127969 (Beltran twins)", Source="publication"),
    dict(Layer="Single cell", Series="GSE144744 (Kaufmann PBMC)", Source="publication"),
    dict(Layer="Proteomics", Series="Bader CSF timsTOF", Source="depositor phenotype table"),
    dict(Layer="Proteomics", Series="Bader CSF Astral", Source="depositor phenotype table"),
    dict(Layer="Proteomics", Series="Wang & Julien brain WM", Source="publication"),
    dict(Layer="Proteomics (external)", Series="UK Biobank-PPP", Source="publication"),
]


def summarise(sex_tbl, meta_tbl, id_col, ds_col, cond_col, layer, source_label):
    s = pd.read_csv(sex_tbl)
    m = pd.read_csv(meta_tbl)
    idc = next(c for c in s.columns if c.lower() in ("sample_id", "sample", "gsm"))
    sexc = next(c for c in s.columns if c.lower() == "sex")
    d = m.merge(s[[idc, sexc]], left_on=id_col, right_on=idc, how="left")

    rows = []
    for ds, g in d.groupby(ds_col):
        n = len(g)
        ms = int((g[cond_col].astype(str).str.upper() == "MS").sum())
        hc = n - ms
        known = g[sexc].isin(["F", "M"])
        f, mm = int((g[sexc] == "F").sum()), int((g[sexc] == "M").sum())
        if known.sum() == 0:
            sex = "not resolvable"
        else:
            sex = f"{f}F / {mm}M ({100 * f / (f + mm):.0f}% F)"
            if known.sum() < n:
                sex += f"; {n - int(known.sum())} unresolved"
        rows.append(dict(Layer=layer, Series=ds, n=n, MS=ms, HC=hc,
                         Sex_source=source_label, Sex=sex, Age="not deposited"))
    return rows


def main():
    rows = []
    rows += summarise(RNA_SEX, EXPR_META, "sample_id", "dataset", "condition",
                      "Bulk RNA", "deposited or inferred (see per-sample table)")
    rows += summarise(METH_SEX, METH_META, "sample_id", "dataset", "condition",
                      "Methylation", "deposited (GSE88824 inferred from chrX/chrY)")

    for q in QUOTED:
        rows.append(dict(Layer=q["Layer"], Series=q["Series"], n="", MS="", HC="",
                         Sex_source=q["Source"],
                         Sex="see source publication", Age="see source publication"))

    d = pd.DataFrame(rows)
    d.to_csv(OUT, index=False)
    n_der = sum(1 for r in rows if r["Sex"] != "see source publication")
    print(f"  {len(d)} rows: {n_der} computed from the per-sample tables, "
          f"{len(QUOTED)} quoted from sources")
    print(d.head(8).to_string(index=False))
    print(f"\n  wrote {OUT}")
    print("  NOTE: the quoted rows carry placeholders. Fill their Sex/Age from the sources before "
          "submission, or keep the delivered hand-curated values for those rows.")


if __name__ == "__main__":
    sys.exit(main())
