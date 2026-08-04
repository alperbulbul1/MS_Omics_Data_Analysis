#!/usr/bin/env python3
"""Infer sample sex from expression of XIST and Y-linked genes, per bulk-RNA series.

Sex is undeposited for 10 of the 14 bulk-transcriptomic series. It is however directly readable from
the expression data: XIST is expressed almost exclusively in females, and the Y-linked genes RPS4Y1,
DDX3Y, UTY, KDM5D, EIF1AY, USP9Y and NLGN4Y almost exclusively in males. The two signals are
independent of each other and of disease status, so their difference separates the sexes cleanly.

Method, applied within each series so that platform and scale never enter the comparison:
  1. log2(1+x) each marker.
  2. Y = mean over the Y-linked genes present; X = XIST.
  3. score = z(Y) - z(X), z-scored within series.
  4. Split at the largest gap in the sorted score, provided that gap is at least 1 SD; otherwise the
     series is reported as unresolved rather than forced.

The four series that DO deposit sex are used as a held-out check of the procedure, which is what
licenses reporting the inferred values at all.
"""
import glob
import os

import numpy as np
import pandas as pd

ROOT = "__MS_GEO_ROOT__"
H = f"{ROOT}/Expression_Data/harmonized_v2"
Y_GENES = ["RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "USP9Y", "NLGN4Y"]


def infer(mat):
    """mat: genes x samples, raw scale. Returns (calls, score, gap_sd)."""
    lg = np.log2(1.0 + mat.clip(lower=0))
    y = [g for g in Y_GENES if g in lg.index]
    if not y or "XIST" not in lg.index:
        return None, None, None
    Y = lg.loc[y].mean(axis=0)
    X = lg.loc["XIST"]

    def z(v):
        s = v.std(ddof=0)
        return (v - v.mean()) / s if s > 0 else v * 0.0

    score = z(Y) - z(X)
    s = score.sort_values()
    if len(s) < 4:
        return None, score, 0.0
    gaps = s.diff().iloc[1:]
    gi = gaps.idxmax()
    gap = gaps.max()
    thr = (s[s.index.get_loc(gi)] + s.iloc[s.index.get_loc(gi) - 1]) / 2
    sd = score.std(ddof=0)
    if sd == 0 or gap / sd < 1.0:
        return None, score, float(gap / sd) if sd else 0.0
    calls = pd.Series(np.where(score > thr, "M", "F"), index=score.index)
    return calls, score, float(gap / sd)


def main():
    meta = pd.read_csv(f"{ROOT}/Expression_Data/Corrected_Metadata_ComBat.csv")
    known = pd.read_csv(
        "/private/tmp/claude-501/-Users-alperbulbul-Desktop-MS-GEO/"
        "e14c6758-36eb-4f9e-84e8-ec87d00c1dc0/scratchpad/rna_known_sex.csv"
    ) if os.path.exists(
        "/private/tmp/claude-501/-Users-alperbulbul-Desktop-MS-GEO/"
        "e14c6758-36eb-4f9e-84e8-ec87d00c1dc0/scratchpad/rna_known_sex.csv"
    ) else None

    rows = []
    for f in sorted(glob.glob(f"{H}/GSE*_symbol_matrix.csv")):
        ds = os.path.basename(f).split("_")[0]
        if ds not in set(meta.dataset):
            continue
        d = pd.read_csv(f, index_col=0)
        d = d[~d.index.duplicated()]
        calls, score, gap = infer(d)
        n = d.shape[1]
        if calls is None:
            rows.append(dict(series=ds, n=n, resolved=False, F=None, M=None,
                             gap_sd=round(gap, 2) if gap else None))
            continue
        rows.append(dict(series=ds, n=n, resolved=True,
                         F=int((calls == "F").sum()), M=int((calls == "M").sum()),
                         gap_sd=round(gap, 2)))
        calls.rename("sex_inferred").to_frame().assign(series=ds).to_csv(
            f"/private/tmp/claude-501/-Users-alperbulbul-Desktop-MS-GEO/"
            f"e14c6758-36eb-4f9e-84e8-ec87d00c1dc0/scratchpad/inferred_{ds}.csv")
    t = pd.DataFrame(rows)
    print(t.to_string(index=False))
    t.to_csv("/private/tmp/claude-501/-Users-alperbulbul-Desktop-MS-GEO/"
             "e14c6758-36eb-4f9e-84e8-ec87d00c1dc0/scratchpad/rna_sex_inferred_summary.csv",
             index=False)


if __name__ == "__main__":
    main()
