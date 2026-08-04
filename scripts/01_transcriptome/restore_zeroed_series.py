#!/usr/bin/env python3
"""restore_zeroed_series.py
=========================
Repair step that must run BEFORE correct_and_normalize.py.

THE DEFECT THIS FIXES
---------------------
Expression_Data/Combined_Expression_Pre_ComBat.csv (16 March) carries four series as exactly zero
across every gene:

    GSE190847   121 samples   B-cell stratum      63.4% of that stratum's pooling weight
    GSE137143    80 samples   T cells (52) + monocytes (28)
    GSE172009     8 samples   T-cell stratum
    GSE207680     5 samples   Brain-WM stratum

correct_and_normalize.py reads that March file and propagates the zeros into
Corrected_Expression_Pre_ComBat.csv, which is the matrix entering batch correction and every
per-stratum limma design. The June re-harmonisation (harmonized_v2/*) has all four intact; the
pipeline was simply never repointed at it.

WHY IT MATTERS STATISTICALLY
----------------------------
With the per-stratum design ~ condition + batch, an all-zero series forms a batch whose within-batch
MS-vs-HC contrast is exactly zero. It cannot flip the sign of an effect, but it dilutes the pooled
estimate by a fixed, gene-independent factor and inflates the residual degrees of freedom. Effects
come out too small and p-values too significant. In the T-cell stratum the dilution factor was 0.917
and was identical to three decimals for every gene checked.

WHAT THIS SCRIPT DOES
---------------------
Replaces those four series' columns in the March matrix with their real values from the per-series
harmonised matrices, then writes the repaired matrix. Sample identifiers differ between the two
sources, so each series is mapped explicitly rather than positionally:

    GSE190847, GSE137143   harmonised columns are already GSM accessions -> identity
    GSE172009              columns HC1..RRMS4  -> SOFT !Sample_title prefix (HC1_CD4 -> HC1)
    GSE207680              columns R2100154... -> series-matrix !Sample_description line that holds
                           the run identifiers (there are two such lines; the generic one is skipped)

VALIDATION
----------
Two checks, both run at the end and both fatal on failure:
  * no all-zero column survives;
  * the two strata that never contained zeros - PBMC and whole blood - are byte-identical to the
    input, which is what distinguishes a repair from a reshape.

Downstream, correct_and_normalize.py carries an assertion that refuses to proceed if any dataset is
all-zero after normalisation, so this defect cannot recur silently.
"""
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd

ED = "__MS_GEO_ROOT__/Expression_Data"
SRC = f"{ED}/Combined_Expression_Pre_ComBat.csv"
OUT = f"{ED}/Combined_Expression_Pre_ComBat_REPAIRED.csv"
ZEROED = ["GSE137143", "GSE190847", "GSE172009", "GSE207680"]
CLEAN_CHECK = ["GSE21942", "GSE138064", "GSE103005", "GSE288904", "GSE43591"]


def soft_field(ds, field):
    """Map GSM -> first value of `field` from the SOFT family file."""
    for p in (f"{ED}/{ds}/{ds}_family.soft.gz", f"{ED}/{ds}_family.soft.gz"):
        if not os.path.exists(p):
            continue
        cur, out = None, {}
        with gzip.open(p, "rt", errors="ignore") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("^SAMPLE"):
                    cur = line.split("=", 1)[1].strip()
                elif cur and line.startswith(field):
                    out.setdefault(cur, line.split("=", 1)[1].strip())
        return out
    return {}


def run_id_map(ds):
    """GSE207680: the run identifiers live in the !Sample_description line that matches R\\d+.

    The series matrix carries two !Sample_description lines; the first is generic prose
    ("rRNA-depleted total RNA") and only the second holds R2100154 etc.
    """
    for p in (f"{ED}/{ds}/{ds}_series_matrix.txt.gz", f"{ED}/{ds}_series_matrix.txt.gz"):
        if not os.path.exists(p):
            continue
        gsm = desc = None
        with gzip.open(p, "rt", errors="ignore") as fh:
            for line in fh:
                parts = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
                if parts[0] == "!Sample_geo_accession":
                    gsm = parts[1:]
                elif parts[0] == "!Sample_description" and any(
                    re.fullmatch(r"R\d+", x) for x in parts[1:]
                ):
                    desc = parts[1:]
        if gsm and desc:
            return dict(zip(desc, gsm))
    return {}


def build_maps():
    maps = {}
    for ds in ("GSE137143", "GSE190847"):
        cols = list(pd.read_csv(f"{ED}/harmonized_v2/{ds}_symbol_matrix.csv", nrows=0).columns)[1:]
        maps[ds] = {c: c for c in cols}
    maps["GSE172009"] = {v.split("_")[0]: g for g, v in soft_field("GSE172009", "!Sample_title").items()}
    maps["GSE207680"] = run_id_map("GSE207680")
    return maps


def main():
    expr = pd.read_csv(SRC, index_col=0)
    meta = pd.read_csv(f"{ED}/Combined_Metadata.csv")
    meta["base"] = meta["dataset"].astype(str).str.split("__").str[0]
    before = {d: expr[[c for c in expr.columns
                       if c in set(meta.loc[meta.base == d, "sample_id"].astype(str))]].copy()
              for d in CLEAN_CHECK}

    maps = build_maps()
    replaced = 0
    for ds in ZEROED:
        h = pd.read_csv(f"{ED}/harmonized_v2/{ds}_symbol_matrix.csv", index_col=0)
        h = h[~h.index.duplicated(keep="first")]
        h = h.rename(columns={c: maps[ds].get(c.replace(".counts.htseq", ""), maps[ds].get(c, c))
                              for c in h.columns})
        want = set(meta.loc[meta.base == ds, "sample_id"].astype(str))
        cols = [c for c in expr.columns if c in want and c in h.columns]
        missing = sorted(want - set(cols))
        assert not missing, f"{ds}: eslesmeyen {len(missing)} ornek, ilk: {missing[:3]}"
        expr[cols] = h.reindex(index=expr.index)[cols].values
        replaced += len(cols)
        print(f"  {ds:<11} {len(cols):>4} sutun geri yuklendi "
              f"(medyan ornekler-arasi SD {np.nanmedian(h[cols].std(axis=1)):.3f})")

    zero = int((expr.abs().sum(axis=0) == 0).sum())
    assert zero == 0, f"hala {zero} tamamen-sifir sutun var"
    for d, old in before.items():
        assert old.equals(expr[old.columns]), f"{d} degismis olmamaliydi"
    print(f"\n  toplam {replaced} sutun | tamamen-sifir sutun: 0 | "
          f"dokunulmayan {len(CLEAN_CHECK)} seri degismedi")

    expr.to_csv(OUT)
    print(f"  yazildi -> {OUT}")
    print("  correct_and_normalize.py bu dosyayi girdi olarak kullanmalidir.")


if __name__ == "__main__":
    sys.exit(main())
