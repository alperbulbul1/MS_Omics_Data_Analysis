#!/usr/bin/env python3
"""figure_constants.py — single source of truth for the MS multi-omics manuscript figures.

Canonical gene panels, proteomic-dataset metadata, and cohort sizes. Every ``make_*.py``
figure generator should import from here instead of hardcoding values, so a figure can
never silently disagree with the analysis outputs or with another figure.

Canonical references: ``README.md`` · ``figure1_workflow.py`` · ``figure4_proteomics.py``.

Import pattern (works regardless of the current working directory)::

    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from figure_constants import INV_TIER1, CANDIDATE_PANEL, PROTEOMICS, prot_label, prot_n
"""

# ── Gene panels ──────────────────────────────────────────────────────────────
# THE five INV-concordant Tier-1 candidates: RNA x methylation inverse-concordant,
# BH-FDR<0.05 in *both* layers, opposite direction of effect. This is the only set
# that should ever be labelled "Tier-1" (README + every *_INV figure).
# HLA-E REMOVED: its bulk-RNA down-regulation did not survive proper cross-dataset
# quantile-normalised ComBat (the prior significance was a normalization artifact).
INV_TIER1 = ["ITGB2","IKZF1"]   # SH3BP4 reclassified to Tier-2 auxiliary

# Secondary proteomic-anchor tier shown alongside INV_TIER1 on the proteomic panels.
TIER2_PROT = ["CTSZ", "CHL1", "ICAM1", "FOXP3", "ITGAL"]

# The seven cross-omics candidates ("CO7") used by the pathway-network figures.
# DISTINCT from INV_TIER1 — 7 cross-omics genes != 6 Tier-1. Do not conflate the two.
CO7 = ["LXN", "SH3BP4", "CHL1", "CTSZ", "RPAP2", "PCNP", "THRB"]

# Broad candidate-highlight panel, formerly (mis)named "TIER1" in v4/v5/figs4to7.
# These are NOT the 6 Tier-1 genes — they are the wider set of candidates annotated on
# the volcano/heatmap panels, so legends must say "candidates", never "Tier-1".
# The order is canonical: the 17-/16-gene variants used elsewhere are strict prefixes
# of this list (see CANDIDATE_PANEL[:N]).
CANDIDATE_PANEL = [
    "ITGB2", "CTSZ", "CHL1", "LXN", "THRB", "ITGAL", "CD79B", "IFI44L",
    "IKZF1", "SH3BP4", "RPAP2", "SLAMF1", "PCNP", "STAT3", "TYK2", "ICAM1", "CASP6",
    "MX1", "IFIT1", "MOSPD3", "FOXP3", "NUP210", "DAXX", "CASP8", "DGKQ",
]

# ── Cohort sizes (canonical — figure1_workflow.py / README) ────────────
COHORTS = {
    "rna":         {"n_series": 14, "n_samples": 472},   # ComBat discovery matrix (291 MS / 181 HC); the 7 strata that produce the reported results hold 462 samples from 13 series
    "methylation": {"n_series": 9, "n_combat": 475},   # 8 arrays + GSE173787 WGBS; AllMeth ComBat = 475 samples
    "scrna":       {"n_cohorts": 3, "n_donors": 81},   # Jäkel 9 + Kaufmann 62 (GSE144744) + Beltrán 10
}

# ── Proteomic datasets (canonical labels — figure4_proteomics.py) ──────────
# CRITICAL CORRECTION: the "timsTOF" DE table is the SECOND PLATFORM of the CSF
# Bader & Mann 2024 cohort (same samples as Astral, different instrument). It is NOT a
# brain dataset and NOT Wang & Julien. Earlier v4/v5/figs4to7 scripts mislabelled it
# "Brain timsTOF — Wang & Julien 2025, n=52"; that is wrong on tissue, citation and n.
PROTEOMICS = {
    "csf_astral":  {"compartment": "CSF",   "platform": "Astral DIA-MS",
                    "citation": "Bader & Mann 2024", "n_ms": 978,  "n_hc": 306},
    "csf_timstof": {"compartment": "CSF",   "platform": "timsTOF DIA",
                    "citation": "Bader & Mann 2024", "n_ms": 1536, "n_hc": 2363},
    "blood_ukb":   {"compartment": "Blood", "platform": "UK Biobank Olink",
                    "citation": "Jacobs 2024", "n_ms": 407, "n_hc": 39979},
    # Region-resolved brain white matter (CTX / NAWM / WML), n=8 per group.
    # UNRESOLVED PROVENANCE: the data files are named ``Magliozzi2026_*``, the _INV
    # figure labels this "Wang & Julien 2025", and README says "Wang & Julien 2026".
    # The author/year is inconsistent across the repo — the authors must confirm which
    # is correct. Left at the _INV value here pending that decision.
    "brain":       {"compartment": "Brain WM", "platform": "region-resolved DIA",
                    "citation": "Wang & Julien 2025", "n_per_group": 8},
}


def prot_label(key):
    """'Compartment - Platform (Citation)' label for a proteomic dataset key."""
    d = PROTEOMICS[key]
    return f"{d['compartment']} · {d['platform']} ({d['citation']})"


def prot_n(key):
    """'n=<MS> MS / <HC> HC' (or 'n=<k>/group') string for a proteomic dataset key."""
    d = PROTEOMICS[key]
    if "n_per_group" in d:
        return f"n={d['n_per_group']}/group"
    return f"n={d['n_ms']:,} MS / {d['n_hc']:,} HC"


# ── Hard-coded display statistics still pending data-driven regeneration ──────
# These values are overlaid as text on poster panels in make_extra_figures.py. They are
# centralised here so there is ONE place to audit them, but they were transcribed from a
# manuscript draft and are NOT yet read from the analysis outputs.
#   TODO: regenerate each from its source table and drop this dict.
FIG_STATS = {
    # Astral vs timsTOF CSF cross-platform concordance (both Bader & Mann CSF cohort).
    # SOURCE: Proteomics/.../Cell2026_concordance_Astral_vs_timsTOF.* — VERIFY.
    "csf_astral_timstof_concordance": {"pearson_r": 0.449, "n": 1762},
    # Per-tissue scRNA MS-vs-HC effect sizes. SOURCE: SingleCell per-tissue DE table — VERIFY.
    "hla_e_csf_cd8":  {"cohens_d": -0.87, "fdr": 2.9e-17},
    # CTSZ multi-layer FDRs. SOURCE: respective per-layer DE tables — VERIFY.
    "ctsz_pbmc_rna":  {"fdr": 8.6e-5},
    "ctsz_meth_wb_dmf": {"fdr": 7.9e-15},
    "ctsz_csf_astral": {"fdr": 4.5e-10},
    "hla_e_meth_combined": {"fdr": 8.1e-8},
}
