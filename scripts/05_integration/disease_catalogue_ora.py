#!/usr/bin/env python3
"""Hypergeometric disease-catalogue over-representation for the inverse-concordant pool.

Methods 4.8: "Disease-catalogue overlap of the full 94-gene discovery pool used a hypergeometric
test against the genes testable in both discovery layers." No script in the release performed it;
the earlier ORA_INV97/run_ora_inv97.py was written for the 97-gene pool and read its inputs from
/tmp files that no script produced.

This version derives BOTH inputs from the canonical result tables, so the test is reproducible:

  foreground  Methylation/results/INVERSE_CONCORDANT_by_gene.tsv          (the pool)
  background  union of the per-stratum RNA gene universes, intersected with the methylation
              gene universe from 15_genelevel_weighting_corrected.tsv    ("testable in both layers")

The test is gseapy's offline hypergeometric against Enrichr disease libraries, with the background
passed explicitly - not Enrichr's default whole-genome background, which would inflate every term.

NOTE ON THE BACKGROUND SIZE. The manuscript reports 8,111 genes. This script recomputes the
intersection from the current tables and prints what it gets; if the two differ, the printed value
is the reproducible one and the discrepancy needs resolving before the number is quoted again.
Nothing here hardcodes a background size.
"""
import glob
import os
import sys

import pandas as pd

ROOT = "__MS_GEO_ROOT__"
POOL = os.path.join(ROOT, "Methylation", "results", "INVERSE_CONCORDANT_by_gene.tsv")
METH = os.path.join(ROOT, "Methylation", "results", "15_genelevel_weighting_corrected.tsv")
RNA_GLOB = os.path.join(ROOT, "Transcriptome", "results", "0*_DE.tsv")
OUT = os.path.join(ROOT, "ORA_INV97")

LIBS = {"DisGeNET": "DisGeNET",
        "GWAS_Catalog": "GWAS_Catalog_2023",
        "Jensen_DISEASES": "Jensen_DISEASES"}


def gene_col(d):
    return (d["gene"] if "gene" in d.columns else d.iloc[:, 0]).astype(str)


def main():
    pool = sorted(set(gene_col(pd.read_csv(POOL, sep="\t"))))
    meth = set(gene_col(pd.read_csv(METH, sep="\t")))

    rna = set()
    files = sorted(glob.glob(RNA_GLOB))
    if not files:
        raise FileNotFoundError(f"no per-stratum DE tables matched {RNA_GLOB}")
    for f in files:
        rna |= set(gene_col(pd.read_csv(f, sep="\t")))

    bg = sorted(rna & meth)
    print(f"  foreground (inverse-concordant pool) : {len(pool)} genes")
    print(f"  RNA universe over {len(files)} strata      : {len(rna)} genes")
    print(f"  methylation universe                 : {len(meth)} genes")
    print(f"  background (testable in both layers) : {len(bg)} genes")
    missing = [g for g in pool if g not in set(bg)]
    if missing:
        print(f"  WARNING: {len(missing)} pool genes are absent from the background: {missing[:10]}")

    try:
        import gseapy as gp
    except ImportError:
        print("\n  gseapy is not installed; install it to run the enrichment "
              "(pip install gseapy). Gene lists above are still valid.")
        return 0

    os.makedirs(OUT, exist_ok=True)
    for name, lib in LIBS.items():
        try:
            enr = gp.enrich(gene_list=pool, gene_sets=lib, background=bg, outdir=None)
            r = enr.results.sort_values("Adjusted P-value").copy()
            # gseapy does not report the term size in this version, and the manuscript quotes
            # overlap/term-size pairs. Recover it from the library, counted against the same
            # background the test used, so the denominator matches the test.
            sets = gp.parser.download_library(lib, "Human")
            bgset = set(bg)
            r["n_hits"] = [len(str(g).split(";")) for g in r["Genes"]]
            r["term_size_in_background"] = [len(set(sets.get(term, [])) & bgset) for term in r["Term"]]
            r["Overlap"] = [f"{a}/{b}" for a, b in zip(r.n_hits, r.term_size_in_background)]
            r.to_csv(os.path.join(OUT, f"ORA_{name}.csv"), index=False)
            top = r.head(3)
            print(f"\n  {name}: {len(r)} terms, {(r['Adjusted P-value'] < 0.05).sum()} at q < 0.05")
            ov = "Overlap" if "Overlap" in r.columns else (
                "Genes" if "Genes" in r.columns else None)
            for _, x in top.iterrows():
                extra = f" {x[ov]}" if ov else ""
                print(f"    {str(x['Term'])[:52]:<54}{extra:<12} q={x['Adjusted P-value']:.3g}")
        except Exception as e:                       # noqa: BLE001 - report, do not abort the rest
            print(f"  {name}: failed ({type(e).__name__}: {e})")
    print(f"\n  wrote per-library tables to {OUT}")


if __name__ == "__main__":
    sys.exit(main())
