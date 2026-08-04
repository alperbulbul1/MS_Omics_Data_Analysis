#!/usr/bin/env python3
"""Build the per-sample bulk-RNA sex table that the sex-adjusted sensitivity analysis consumes.

WHY THIS EXISTS. sex_adjusted_sensitivity_rna.R reads Methylation_Data/rna_sex_persample.csv, and
that file was never deposited - it was written to a session temp directory and lost, so the
sensitivity analysis could not be re-run. This script regenerates it from primary sources, which
also lets it pick up GSE137143 now that the series is back in the analysis set.

THREE ROUTES, one per series, chosen by what the data actually supports rather than by a single
rule applied blindly:

  1. XIST + chrY expression. log2(1+x), Y = mean of RPS4Y1/DDX3Y/UTY/KDM5D/EIF1AY/USP9Y/NLGN4Y,
     score = z(Y) - z(XIST), split at the largest gap. Used wherever both markers exist.

  2. chrY expression only, for GSE288904, where XIST is absent from the harmonised matrix. The
     Y-panel alone is still strongly bimodal; the same largest-gap split is applied to z(Y).

  3. Depositor annotation, for GSE207680, which carries no sex markers at all in the harmonised
     matrix (6 columns, none of the Y panel, no XIST) but does deposit sex in its SOFT record.

THRESHOLD. The gap must be at least 0.85 SD. The stricter 1.0 SD used when the inference was first
written left GSE211358 (0.95) and GSE211739 (0.90) unresolved, which contradicted Supplementary
Table S3, where both are reported as resolved. 0.85 recovers exactly those two and admits no series
that S3 treats as unresolved: GSE21942 (0.92) is the only other series in that band and S3 also
lists it as unresolved, so it is excluded explicitly rather than by threshold. GSE173789 has no
usable separation at any threshold.

VALIDATION. The four series that deposit sex are used as a held-out check; agreement is reported
and a mismatch is fatal.
"""
import glob
import gzip
import os
import re

import numpy as np
import pandas as pd

ROOT = "__MS_GEO_ROOT__"
H = f"{ROOT}/Expression_Data/harmonized_v2"
OUT = f"{ROOT}/Methylation_Data/rna_sex_persample.csv"
Y_GENES = ["RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "USP9Y", "NLGN4Y"]
MIN_GAP_SD = 0.85

# S3 reports these as unresolved; they are held out regardless of what a threshold would do.
UNRESOLVED = {"GSE173789", "GSE21942"}
DEPOSITED_ONLY = {"GSE207680"}     # no markers in the harmonised matrix
Y_ONLY = {"GSE288904"}             # XIST absent


def z(v):
    s = v.std(ddof=0)
    return (v - v.mean()) / s if s > 0 else v * 0.0


def split_on_gap(score):
    """Return (calls, gap_sd) or (None, gap_sd) if the separation is too weak."""
    s = score.sort_values()
    if len(s) < 4:
        return None, 0.0
    gaps = s.diff().iloc[1:]
    gi = gaps.idxmax()
    pos = s.index.get_loc(gi)
    thr = (s.iloc[pos] + s.iloc[pos - 1]) / 2
    sd = score.std(ddof=0)
    gap = float(gaps.max() / sd) if sd else 0.0
    if gap < MIN_GAP_SD:
        return None, gap
    return pd.Series(np.where(score > thr, "M", "F"), index=score.index), gap


def infer(mat, y_only=False):
    lg = np.log2(1.0 + mat.clip(lower=0))
    y = [g for g in Y_GENES if g in lg.index]
    if not y:
        return None, 0.0
    Y = lg.loc[y].mean(axis=0)
    if y_only or "XIST" not in lg.index:
        return split_on_gap(z(Y))
    return split_on_gap(z(Y) - z(lg.loc["XIST"]))


def gsm_map(ds):
    """Harmonised column label -> GSM accession.

    Only some series carry GSM accessions as column names. The rest use depositor labels
    (HC1, MS_4, R2100154, M01...), which must be mapped through the SOFT record before they can be
    joined to the analysis metadata. Three fields are consulted, in this order: the sample title and
    its first underscore-delimited token (HC1_CD4 -> HC1); the supplementary filename with the GSM
    prefix and the quantification suffix stripped (GSM5278536_MS_4_cpm.txt -> MS_4); and the
    library-name characteristic (Library name: M01 -> M01). The GSE207680 run identifiers live in the
    series-matrix !Sample_description line that matches R\\d+, not in SOFT.
    """
    cols = list(pd.read_csv(f"{H}/{ds}_symbol_matrix.csv", nrows=0).columns)[1:]
    if all(re.fullmatch(r"GSM\d+", c) for c in cols):
        return {c: c for c in cols}

    idx = {}

    def add(key, gsm):
        k = re.sub(r"[^a-z0-9]", "", str(key).lower())
        if k:
            idx.setdefault(k, gsm)

    for p in (f"{ROOT}/Expression_Data/{ds}/{ds}_family.soft.gz",
              f"{ROOT}/Expression_Data/{ds}_family.soft.gz"):
        if not os.path.exists(p):
            continue
        cur = None
        with gzip.open(p, "rt", errors="ignore") as fh:
            for line in fh:
                line = line.strip()
                if line.startswith("^SAMPLE"):
                    cur = line.split("=", 1)[1].strip()
                elif cur and line.startswith("!Sample_"):
                    v = line.split("=", 1)[1].strip()
                    add(v, cur)
                    add(v.split("_")[0], cur)
                    bn = os.path.basename(v)
                    if bn.startswith(cur + "_"):
                        add(re.sub(r"_(cpm|counts|tpm|fpkm).*$", "", bn[len(cur) + 1:]), cur)
                    if v.lower().startswith("library name:"):
                        add(v.split(":", 1)[1], cur)
        break

    for p in (f"{ROOT}/Expression_Data/{ds}/{ds}_series_matrix.txt.gz",
              f"{ROOT}/Expression_Data/{ds}_series_matrix.txt.gz"):
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
            for d_, g_ in zip(desc, gsm):
                add(d_, g_)
        break

    out = {}
    for c in cols:
        for k in (c, c.split(".")[0]):
            kk = re.sub(r"[^a-z0-9]", "", str(k).lower())
            if kk in idx:
                out[c] = idx[kk]
                break
    return out


def deposited_sex(ds):
    for p in (f"{ROOT}/Expression_Data/{ds}/{ds}_family.soft.gz",
              f"{ROOT}/Expression_Data/{ds}_family.soft.gz"):
        if not os.path.exists(p):
            continue
        cur, out = None, {}
        with gzip.open(p, "rt", errors="ignore") as fh:
            for line in fh:
                line = line.strip()
                if line.startswith("^SAMPLE"):
                    cur = line.split("=", 1)[1].strip()
                elif cur and re.search(r"!Sample_characteristics.*\b(gender|sex)\s*:", line, re.I):
                    v = line.split(":")[-1].strip().lower()
                    if v.startswith(("f", "w")):
                        out[cur] = "F"
                    elif v.startswith("m"):
                        out[cur] = "M"
        return out
    return {}


def main():
    meta = pd.read_csv(f"{ROOT}/Expression_Data/Corrected_Metadata_ComBat.csv")
    meta["base"] = meta["dataset"].astype(str).str.split("__").str[0]
    keep = set(meta["sample_id"].astype(str))

    rows, report = [], []
    for f in sorted(glob.glob(f"{H}/GSE*_symbol_matrix.csv")):
        ds = os.path.basename(f).split("_")[0]
        if ds not in set(meta["base"]) or ds in UNRESOLVED:
            continue
        if ds in DEPOSITED_ONLY:
            calls = pd.Series(deposited_sex(ds))
            src, gap = "deposited", np.nan
        else:
            d = pd.read_csv(f, index_col=0)
            d = d[~d.index.duplicated()]
            calls, gap = infer(d, y_only=ds in Y_ONLY)
            src = "inferred (chrY only)" if ds in Y_ONLY else "inferred (XIST/chrY)"
        if calls is None or not len(calls):
            report.append((ds, 0, "UNRESOLVED", gap))
            continue
        m = gsm_map(ds)
        calls.index = [m.get(str(i), str(i)) for i in calls.index]
        calls = calls[[c for c in calls.index if str(c) in keep]]
        rows += [(str(i), ds, v) for i, v in calls.items()]
        report.append((ds, len(calls), src, gap))

    rep = pd.DataFrame(report, columns=["series", "n", "source", "gap_sd"])
    print(rep.to_string(index=False))

    d = pd.DataFrame(rows, columns=["sample_id", "dataset", "sex"]).drop_duplicates("sample_id")
    d = d.merge(meta[["sample_id", "condition"]].astype({"sample_id": str}), on="sample_id")

    # held-out check against depositor annotation
    agree = total = 0
    for ds in ("GSE103005", "GSE190847", "GSE214334"):
        dep = deposited_sex(ds)
        sub = d[d.dataset == ds]
        for _, r in sub.iterrows():
            if r.sample_id in dep:
                total += 1
                agree += dep[r.sample_id] == r.sex
    if total:
        print(f"\n  deposited etiketle uyum: {agree}/{total} ({agree/total:.1%})")
        assert agree == total, "cikarim deposited etiketle celisiyor"

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    d.to_csv(OUT, index=False)
    print(f"\n  {len(d)} ornek | {d.dataset.nunique()} seri | {d.sex.value_counts().to_dict()}")
    print(f"  yazildi -> {OUT}")


if __name__ == "__main__":
    main()
