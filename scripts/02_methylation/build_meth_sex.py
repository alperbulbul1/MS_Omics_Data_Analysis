#!/usr/bin/env python3
"""Build Methylation_Data/meth_sex.csv, the per-sample sex table sex_adjusted_sensitivity.R reads.

WHY THIS EXISTS. sex_adjusted_sensitivity.R (Reviewer 3) fits the AllMeth MS-vs-HC model with and
without a sex term on the samples for which sex is known, and reads that sex from meth_sex.csv. No
script produced that file: it was assembled by hand, so the released sensitivity analysis could not
be re-run. This rebuilds it from the deposited records.

WHERE SEX COMES FROM. For each of the eight array series in the AllMeth matrix, the GEO series
matrix carries a `!Sample_characteristics_ch1` line of the form "Sex: F" / "gender: male". Those
deposited labels are used verbatim. GSE88824 deposits no sex, and its calls come instead from
infer_sex_GSE88824_chrXY.R, which classifies on chrX/chrY methylation signal and writes
gse88824_sex.csv; if that file is absent those samples are simply left unsexed rather than guessed,
and the sensitivity analysis drops them, which is the behaviour its own comment describes.

The output has one row per sample with columns sample_id, dataset, sex (F/M/NA). Samples whose
deposited value is neither male nor female are written as NA rather than dropped, so the count of
sexed samples is visible in the file itself.
"""
import glob
import gzip
import os
import re
import sys

import pandas as pd

ROOT = "__MS_GEO_ROOT__"
MD = os.path.join(ROOT, "Methylation_Data")
META = os.path.join(MD, "AllMeth_ComBat_Metadata.csv")
GSE88824_SEX = os.path.join(MD, "gse88824_sex.csv")
OUT = os.path.join(MD, "meth_sex.csv")

SEX_KEY = re.compile(r"^\s*(sex|gender)\s*:\s*(.+?)\s*$", re.I)


def norm_sex(v):
    v = str(v).strip().lower()
    if v.startswith("f") or v == "female":
        return "F"
    if v.startswith("m") or v == "male":
        return "M"
    return None


def sex_from_series_matrix(path):
    """GSM -> F/M from the deposited !Sample_characteristics_ch1 lines."""
    gsms, rows = [], []
    with gzip.open(path, "rt", errors="ignore") as fh:
        for ln in fh:
            if ln.startswith("!Sample_geo_accession"):
                gsms = [x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
            elif ln.startswith("!Sample_characteristics_ch1"):
                rows.append([x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]])
            elif ln.startswith("!series_matrix_table_begin"):
                break
    out = {}
    for row in rows:
        for gsm, cell in zip(gsms, row):
            m = SEX_KEY.match(cell)
            if m:
                s = norm_sex(m.group(2))
                if s:
                    out[gsm] = s
    return out


def main():
    meta = pd.read_csv(META)
    calls = {}
    for gse in sorted(meta.dataset.unique()):
        hits = glob.glob(f"{MD}/**/{gse}_series_matrix.txt.gz", recursive=True)
        if not hits:
            print(f"  {gse}: no series matrix, skipped")
            continue
        s = sex_from_series_matrix(hits[0])
        calls.update(s)
        print(f"  {gse}: {len(s)} deposited sex calls")

    if os.path.exists(GSE88824_SEX):
        g = pd.read_csv(GSE88824_SEX)
        col = next((c for c in g.columns if c.lower() in ("sex", "predicted_sex")), None)
        idc = next((c for c in g.columns if "sample" in c.lower() or c.lower().startswith("gsm")), g.columns[0])
        n = 0
        for _, r in g.iterrows():
            s = norm_sex(r[col]) if col else None
            if s:
                calls.setdefault(str(r[idc]), s)
                n += 1
        print(f"  GSE88824: {n} inferred calls from infer_sex_GSE88824_chrXY.R")
    else:
        print(f"  GSE88824: {os.path.basename(GSE88824_SEX)} absent -> those samples left unsexed "
              "(run infer_sex_GSE88824_chrXY.R first to include them)")

    meta["sex"] = [calls.get(str(s)) for s in meta.sample_id]
    meta[["sample_id", "dataset", "sex"]].to_csv(OUT, index=False)
    n_sexed = meta.sex.isin(["F", "M"]).sum()
    print(f"\n  wrote {OUT}")
    print(f"  {len(meta)} samples, {n_sexed} with sex "
          f"({(meta.sex == 'F').sum()}F / {(meta.sex == 'M').sum()}M)")
    per = meta.groupby("dataset").sex.apply(lambda s: s.isin(["F", "M"]).sum())
    for k, v in per.items():
        print(f"    {k}: {v}/{(meta.dataset == k).sum()}")


if __name__ == "__main__":
    sys.exit(main())
