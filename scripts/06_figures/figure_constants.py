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
# Canonical evidence hierarchy used by every manuscript figure and caption.
# Tier-1 requires inverse RNA × methylation concordance plus an independent
# donor-level single-cell or proteomic anchor.  The criterion has no exceptions.
INV_TIER1 = ["ITGB2", "IKZF1"]

# HLA-E is retained for biological follow-up because of its strong cell-level signal,
# but its bulk-RNA arm did not survive cross-dataset normalisation; it is not Tier-1.
SUGGESTIVE = ["HLA-E"]

# Genes passing the inverse RNA × methylation screen without a qualifying orthogonal
# anchor.  CD79B, LXN and SH3BP4 therefore belong here, not in Tier-1.
TIER2_AUX_INV = [
    "CD79B", "LXN", "SH3BP4", "CASP6", "CASP8", "DGKQ", "MX1", "IFIT1",
    "NUP210", "RUNX3",
]

# Strong proteomic candidates without a qualifying inverse RNA × methylation pairing.
# FOXP3 was removed from this group on revision: it is not quantified in ANY of the seven
# proteomic compartments (both CSF instruments, all four brain-region contrasts, UK Biobank-PPP),
# so it could not be a "strong proteomic candidate", which is what defines this group. It is
# retained in the STRING display as a canonical MS immune context gene, which is the role it
# actually plays (the IKZF1-RUNX3/FOXP3-STAT1-STAT3 axis).
TIER2_PROT = ["CTSZ", "CHL1", "ICAM1", "ITGAL"]

# The strict 17-gene panel used for pathway enrichment.  HLA-E is shown separately as
# suggestive wherever the underlying assay contains it and is never counted in a tier.
TIERED_PANEL = INV_TIER1 + TIER2_AUX_INV + TIER2_PROT
DISPLAY_PANEL = INV_TIER1 + SUGGESTIVE + TIER2_AUX_INV + TIER2_PROT

# Backwards-compatible alias for scripts that use the wider displayed candidate set.
CANDIDATE_PANEL = DISPLAY_PANEL

# ── Cohort sizes (canonical — figure1_workflow.py / README) ────────────
COHORTS = {
    "rna":         {"n_series": 14, "n_samples": 472},   # ComBat discovery matrix (291 MS / 181 HC); the 7 strata that produce the reported results hold 462 samples from 13 series
    "methylation": {"n_series": 9, "n_combat": 475},   # 8 arrays + GSE173787 WGBS; AllMeth ComBat = 475 samples
    "scrna":       {"n_cohorts": 3, "n_donors": 79},   # Jäkel 9 + Kaufmann 62 (GSE144744) + Beltrán 8 unique donors
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
    # PROVENANCE RESOLVED: this is Wang & Julien 2026 (ref9 in the manuscript; raw spectra
    # MassIVE MSV000096790). The on-disk tables and the 04cc script are still named
    # ``Magliozzi*`` for historical reasons; that is a filename, not an attribution.
    "brain":       {"compartment": "Brain WM", "platform": "region-resolved DIA",
                    "citation": "Wang & Julien 2026", "n_per_group": 8},
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
