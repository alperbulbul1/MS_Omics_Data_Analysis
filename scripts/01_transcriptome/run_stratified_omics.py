import gzip
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


ROOT = Path("__MS_GEO_ROOT__")
OUT_DIR = ROOT / "Stratified_Analyses"

EXPR_META_PATH = ROOT / "Expression_Data" / "Corrected_Metadata_ComBat.csv"
EXPR_MATRIX_PATH = ROOT / "Expression_Data" / "Corrected_Expression_Pre_ComBat.csv"

# Methylation: prefer IDAT-preprocessed, fall back to a combined matrix
METH_DATA_DIR  = ROOT / "Methylation_Data"
METH_STRICT_DIR = METH_DATA_DIR / "Strict_Array_Preprocessed"
METH_BETAONLY_DIR = METH_DATA_DIR / "Normalized_Beta_Only"

# Combined metadata: strict (IDAT) preferred, fall back to PythonPipeline metadata
METH_META_PATH = (
    METH_STRICT_DIR / "Combined_Methylation_Strict_Metadata.csv"
    if (METH_STRICT_DIR / "Combined_Methylation_Strict_Metadata.csv").exists()
    else METH_DATA_DIR / "Combined_Methylation_Metadata.csv"
)
# Combined M-value matrix: strict preferred, then pre-batch from Python pipeline
METH_MATRIX_PATH = (
    METH_STRICT_DIR / "Combined_Methylation_Strict_M.csv"
    if (METH_STRICT_DIR / "Combined_Methylation_Strict_M.csv").exists()
    else METH_DATA_DIR / "Combined_Methylation_Pre_Batch.csv"
)
# Beta matrix for plotCpg / mCSEA
METH_BETA_PATH = (
    METH_STRICT_DIR / "Combined_Methylation_Strict_Beta.csv"
    if (METH_STRICT_DIR / "Combined_Methylation_Strict_Beta.csv").exists()
    else Path("")
)

# The two R workers ship alongside this file in the release, not under the data root, so resolve
# them from this script's own location. Reading them from ROOT pointed every subprocess call at a
# path that does not exist once configure.sh has substituted the data directory.
HERE = Path(__file__).resolve().parent
EXPR_R = HERE / "run_expression_subgroup_limma.R"
METH_R = HERE.parent / "02_methylation" / "run_methylation_subgroup_limma.R"
for _w in (EXPR_R, METH_R):
    if not _w.exists():
        raise FileNotFoundError(f"R worker not found: {_w}")

MIN_CASES = 4
MIN_CONTROLS = 4

LABEL_ORDER = [
    "Responder",
    "Non-responder",
    "IFNb",
    "DMF",
    "Fingolimod",
    "Glatiramer",
    "Ocrelizumab",
    "Untreated",
    "Relapse",
    "Remission",
    "Smoker",
    "Non-smoker",
]


@dataclass
class AnalysisDefinition:
    omics: str
    analysis_type: str
    subgroup: str
    cell_group: str
    label: str
    sample_ids: list[str]
    datasets: list[str]
    ms_samples: int
    hc_samples: int

    @property
    def slug(self) -> str:
        parts = [self.analysis_type, self.cell_group]
        if self.label:
            parts.append(self.label)
        text = "__".join(parts)
        return re.sub(r"[^A-Za-z0-9]+", "_", text.strip().lower()).strip("_")


def normalize_text(text: str) -> str:
    return f" {str(text or '').lower().replace('_', ' ')} "


def parse_soft_samples(soft_path: Path) -> dict[str, str]:
    samples: dict[str, list[str]] = {}
    current = None
    prefixes = (
        "!Sample_title = ",
        "!Sample_characteristics_ch1 = ",
        "!Sample_source_name_ch1 = ",
        "!Sample_description = ",
    )
    with gzip.open(soft_path, "rt", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith("^SAMPLE = "):
                current = line.split("=", 1)[1].strip()
                samples[current] = []
                continue
            if current is not None:
                if line.startswith("^"):
                    current = None
                    continue
                for prefix in prefixes:
                    if line.startswith(prefix):
                        samples[current].append(line.split("=", 1)[1].strip())
                        break
    return {sample_id: " ".join(parts) for sample_id, parts in samples.items()}


def primary_cell_group(text: str) -> str:
    text = normalize_text(text)
    if any(token in text for token in [" white matter ", " frontal lobe ", " brain ", " cortex ", " corpus callosum "]):
        return "Brain / WM"
    if " bronchoalveolar " in text or " bal " in text:
        return "BAL"
    if any(token in text for token in [" cd14 ", " monocyte "]):
        return "Monocytes"
    if any(token in text for token in [" cd19", " b cell ", " b cells "]):
        return "B cells"
    if any(token in text for token in [" cd4", " cd8", " t cell", " t cells ", " tcell "]):
        return "T cells"
    if " pbmc " in text or " peripheral blood mononuclear cells " in text:
        return "PBMC"
    if any(token in text for token in [" whole blood ", " wholeblood ", " peripheral blood ", " tissue: blood "]):
        return "Whole blood"
    return "Unspecified"


def extract_labels(text: str) -> list[str]:
    text = normalize_text(text)
    patterns = {
        "Responder": [" responder ", " complete responder ", " partial responder "],
        "Non-responder": [" non responder ", " nonresponder ", " non-responder "],
        "IFNb": [" interferon ", " ifn-b ", " ifnb ", " ifn beta ", " ifn-beta ", " rebif ", " avonex "],
        "DMF": [" dimethyl fumarate ", " dmf ", " tecfidera "],
        "Fingolimod": [" fingolimod ", " gilenya "],
        "Glatiramer": [" glatiramer ", " copaxone "],
        "Ocrelizumab": [" ocrevus ", " ocrelizumab "],
        "Untreated": [" untreated ", " treatment: untreated ", " wash out ", " no treatment "],
        "Relapse": [" relapse "],
        "Remission": [" remission "],
        "Smoker": [" smoking status: ever smoker ", " smoker ", " smoking_status: s"],
        "Non-smoker": [" smoking status: never smoker ", " nonsmoker ", " smoking_status: ns"],
    }
    labels = [label for label, tokens in patterns.items() if any(token in text for token in tokens)]
    return labels or ["Not reported"]


def read_expression_sample_metadata() -> pd.DataFrame:
    meta = pd.read_csv(EXPR_META_PATH)
    cache: dict[str, dict[str, str]] = {}
    rows = []
    for base_gse in sorted(meta["dataset"].astype(str).str.split("__").str[0].unique()):
        soft_path = ROOT / "Expression_Data" / base_gse / f"{base_gse}_family.soft.gz"
        if not soft_path.exists():
            soft_path = ROOT / "Expression_Data" / f"{base_gse}_family.soft.gz"
        if not soft_path.exists():
            continue
        cache[base_gse] = parse_soft_samples(soft_path)
        subset = meta[meta["dataset"].astype(str).str.startswith(base_gse)].copy()
        for _, row in subset.iterrows():
            raw_text = cache[base_gse].get(str(row["sample_id"]), "")
            rows.append(
                {
                    "sample_id": row["sample_id"],
                    "dataset": row["dataset"],
                    "condition": row["condition"],
                    "raw_text": raw_text,
                    "cell_group": primary_cell_group(raw_text),
                    "labels": "; ".join(extract_labels(raw_text)),
                }
            )
    df = pd.DataFrame(rows)
    return df


def read_methylation_sample_metadata() -> pd.DataFrame:
    if not METH_META_PATH.exists():
        raise FileNotFoundError(f"Methylation metadata not found: {METH_META_PATH}")
    meta = pd.read_csv(METH_META_PATH)
    meta = meta.copy()
    if "raw_text" not in meta.columns:
        meta["raw_text"] = ""

    # Cell group: prefer cell_type column (always correct from merge), fall back to raw_text parsing
    def _assign_cell_group(row):
        ct = str(row.get("cell_type", "")).lower().strip()
        # Direct mapping from cell_type column
        ct_map = {
            "whole_blood": "Whole blood",
            "pbmc": "PBMC",
            "cd4_tcell": "T cells", "cd4_t_cells": "T cells",
            "cd8_t_cells": "T cells",
            "cd14_monocytes": "Monocytes", "monocytes": "Monocytes",
            "b_cells": "B cells",
            "nk_cells": "NK cells",
            "neutrophils": "Neutrophils",
            "brain_nawm": "Brain / WM",
            "bronchoalveolar_lavage": "BAL",
        }
        if ct in ct_map:
            return ct_map[ct]
        # Fall back to parsing raw_text
        return primary_cell_group(str(row.get("raw_text", "")))

    meta["cell_group"] = meta.apply(_assign_cell_group, axis=1)
    meta["labels"] = meta["raw_text"].fillna("").map(lambda text: "; ".join(extract_labels(text)))
    # Keep base_dataset for batch covariate
    if "base_dataset" not in meta.columns and "dataset" in meta.columns:
        meta["base_dataset"] = meta["dataset"].str.split("__").str[0]
    keep_cols = [c for c in ["sample_id", "dataset", "base_dataset", "condition", "raw_text", "cell_group", "labels"] if c in meta.columns]
    return meta[keep_cols]


def label_present(series: pd.Series, label: str) -> pd.Series:
    pattern = rf"(?:^|; ){re.escape(label)}(?:$|;)"
    return series.fillna("").str.contains(pattern, regex=True)


def build_definitions(df: pd.DataFrame, omics: str) -> list[AnalysisDefinition]:
    definitions: list[AnalysisDefinition] = []

    for cell_group in sorted(group for group in df["cell_group"].dropna().unique() if group != "Unspecified"):
        subset = df[df["cell_group"] == cell_group].copy()
        ms_samples = int((subset["condition"] == "MS").sum())
        hc_samples = int((subset["condition"] == "HC").sum())
        if ms_samples >= MIN_CASES and hc_samples >= MIN_CONTROLS:
            definitions.append(
                AnalysisDefinition(
                    omics=omics,
                    analysis_type="cell_tissue_case_control",
                    subgroup=cell_group,
                    cell_group=cell_group,
                    label="",
                    sample_ids=subset["sample_id"].astype(str).tolist(),
                    datasets=sorted(subset["dataset"].astype(str).unique()),
                    ms_samples=ms_samples,
                    hc_samples=hc_samples,
                )
            )

    for cell_group in sorted(group for group in df["cell_group"].dropna().unique() if group != "Unspecified"):
        cell_subset = df[df["cell_group"] == cell_group].copy()
        for label in LABEL_ORDER:
            ms_subset = cell_subset[(cell_subset["condition"] == "MS") & label_present(cell_subset["labels"], label)].copy()
            if len(ms_subset) < MIN_CASES:
                continue
            datasets = sorted(ms_subset["dataset"].astype(str).unique())
            hc_pool = cell_subset[(cell_subset["condition"] == "HC") & (cell_subset["dataset"].astype(str).isin(datasets))].copy()
            hc_labeled = hc_pool[label_present(hc_pool["labels"], label)].copy()
            hc_subset = hc_labeled if len(hc_labeled) >= MIN_CONTROLS else hc_pool
            if len(hc_subset) < MIN_CONTROLS:
                continue
            combined = pd.concat([ms_subset, hc_subset], ignore_index=True).drop_duplicates(subset="sample_id", keep="first")
            definitions.append(
                AnalysisDefinition(
                    omics=omics,
                    analysis_type="label_context_case_control",
                    subgroup=f"{label} in {cell_group}",
                    cell_group=cell_group,
                    label=label,
                    sample_ids=combined["sample_id"].astype(str).tolist(),
                    datasets=datasets,
                    ms_samples=int((combined["condition"] == "MS").sum()),
                    hc_samples=int((combined["condition"] == "HC").sum()),
                )
            )

    return definitions


def run_analysis(defn: AnalysisDefinition, sample_meta: pd.DataFrame, matrix_path: Path) -> dict:
    omics_dir = OUT_DIR / defn.omics
    out_dir = omics_dir / defn.slug
    out_dir.mkdir(parents=True, exist_ok=True)

    subset_meta = sample_meta[sample_meta["sample_id"].astype(str).isin(defn.sample_ids)].copy()
    subset_meta.to_csv(out_dir / "metadata.csv", index=False)

    if defn.omics == "Expression":
        cmd = ["Rscript", str(EXPR_R), str(out_dir / "metadata.csv"), str(matrix_path), str(out_dir)]
    else:
        beta_arg = str(METH_BETA_PATH) if METH_BETA_PATH and METH_BETA_PATH.exists() else ""
        cmd = ["Rscript", str(METH_R), str(out_dir / "metadata.csv"), str(matrix_path), str(out_dir), beta_arg, "450K"]
    run = subprocess.run(cmd, capture_output=True, text=True)
    (out_dir / "stdout.log").write_text(run.stdout or "")
    (out_dir / "stderr.log").write_text(run.stderr or "")

    record = {
        "omics": defn.omics,
        "analysis_type": defn.analysis_type,
        "subgroup": defn.subgroup,
        "cell_group": defn.cell_group,
        "label": defn.label,
        "datasets": ", ".join(defn.datasets),
        "ms_samples": defn.ms_samples,
        "hc_samples": defn.hc_samples,
        "status": "ok" if run.returncode == 0 else "failed",
        "out_dir": str(out_dir),
    }
    if run.returncode != 0:
        record["error"] = (run.stderr or run.stdout or "").strip()[:4000]
        return record

    summary_path = out_dir / "Summary.csv"
    if summary_path.exists():
        summary_df = pd.read_csv(summary_path)
        for _, row in summary_df.iterrows():
            record[str(row["metric"])] = row["value"]

    if defn.omics == "Expression":
        res_path = out_dir / "DGE_Results_MS_vs_HC.csv"
        if res_path.exists():
            res = pd.read_csv(res_path)
            if not res.empty and "adj.P.Val" in res.columns:
                res["adj.P.Val"] = pd.to_numeric(res["adj.P.Val"], errors="coerce")
                top = res.dropna(subset=["adj.P.Val"]).nsmallest(5, "adj.P.Val")[["Gene", "adj.P.Val"]]
                record["top_hits"] = "; ".join(f"{g} ({p:.2e})" for g, p in top.itertuples(index=False))
    else:
        res_path = out_dir / "DMP_Results_MS_vs_HC.csv"
        if res_path.exists():
            res = pd.read_csv(res_path)
            if not res.empty and "adj.P.Val" in res.columns:
                res["adj.P.Val"] = pd.to_numeric(res["adj.P.Val"], errors="coerce")
                probe_col = "Probe" if "Probe" in res.columns else res.columns[0]
                top = res.dropna(subset=["adj.P.Val"]).nsmallest(5, "adj.P.Val")[[probe_col, "adj.P.Val"]]
                record["top_probe_hits"] = "; ".join(f"{g} ({p:.2e})" for g, p in top.itertuples(index=False))
        # New filename: Promoter_Results_mCSEA.csv (updated pipeline)
        for pname in ["Promoter_Results_mCSEA.csv", "Promoter_Results_MS_vs_HC.csv"]:
            promoter_path = out_dir / pname
            if promoter_path.exists():
                promoter = pd.read_csv(promoter_path)
                padj_col = next((c for c in ["padj", "adj.P.Val", "P.adjust"] if c in promoter.columns), None)
                gene_col  = next((c for c in ["Gene", "gene", "pathway"] if c in promoter.columns), None)
                if padj_col and gene_col and not promoter.empty:
                    promoter[padj_col] = pd.to_numeric(promoter[padj_col], errors="coerce")
                    top_p = promoter.dropna(subset=[padj_col]).nsmallest(5, padj_col)[[gene_col, padj_col]]
                    record["top_promoter_hits"] = "; ".join(f"{g} ({p:.2e})" for g, p in top_p.itertuples(index=False))
                break
        # Gene body results
        gene_path = out_dir / "Gene_Results_mCSEA.csv"
        if gene_path.exists():
            gene_res = pd.read_csv(gene_path)
            padj_col2 = next((c for c in ["padj", "adj.P.Val"] if c in gene_res.columns), None)
            gene_col2  = next((c for c in ["Gene", "gene"] if c in gene_res.columns), None)
            if padj_col2 and gene_col2 and not gene_res.empty:
                gene_res[padj_col2] = pd.to_numeric(gene_res[padj_col2], errors="coerce")
                top_g = gene_res.dropna(subset=[padj_col2]).nsmallest(5, padj_col2)[[gene_col2, padj_col2]]
                record["top_gene_body_hits"] = "; ".join(f"{g} ({p:.2e})" for g, p in top_g.itertuples(index=False))

    return record


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    expr_samples = read_expression_sample_metadata()
    meth_samples = read_methylation_sample_metadata()
    expr_samples.to_csv(OUT_DIR / "Expression_Sample_Subgroups.csv", index=False)
    meth_samples.to_csv(OUT_DIR / "Methylation_Sample_Subgroups.csv", index=False)

    expr_defs = build_definitions(expr_samples, "Expression")
    meth_defs = build_definitions(meth_samples, "Methylation")

    pd.DataFrame([d.__dict__ | {"slug": d.slug, "datasets": ", ".join(d.datasets), "sample_ids": "; ".join(d.sample_ids)} for d in expr_defs]).to_csv(
        OUT_DIR / "Expression_Subgroup_Definitions.csv", index=False
    )
    pd.DataFrame([d.__dict__ | {"slug": d.slug, "datasets": ", ".join(d.datasets), "sample_ids": "; ".join(d.sample_ids)} for d in meth_defs]).to_csv(
        OUT_DIR / "Methylation_Subgroup_Definitions.csv", index=False
    )

    results = []
    for defn in expr_defs:
        results.append(run_analysis(defn, expr_samples, EXPR_MATRIX_PATH))
    for defn in meth_defs:
        results.append(run_analysis(defn, meth_samples, METH_MATRIX_PATH))

    results_df = pd.DataFrame(results)
    expr_results = results_df[results_df["omics"] == "Expression"].copy()
    meth_results = results_df[results_df["omics"] == "Methylation"].copy()
    expr_results.to_csv(OUT_DIR / "Expression_Subgroup_Results_Summary.csv", index=False)
    meth_results.to_csv(OUT_DIR / "Methylation_Subgroup_Results_Summary.csv", index=False)

    with pd.ExcelWriter(OUT_DIR / "Stratified_Analysis_Summary.xlsx", engine="openpyxl") as writer:
      expr_samples.to_excel(writer, sheet_name="Expr_Samples", index=False)
      meth_samples.to_excel(writer, sheet_name="Meth_Samples", index=False)
      pd.DataFrame([d.__dict__ | {"slug": d.slug, "datasets": ", ".join(d.datasets), "sample_ids": "; ".join(d.sample_ids)} for d in expr_defs]).to_excel(
          writer, sheet_name="Expr_Subgroups", index=False
      )
      pd.DataFrame([d.__dict__ | {"slug": d.slug, "datasets": ", ".join(d.datasets), "sample_ids": "; ".join(d.sample_ids)} for d in meth_defs]).to_excel(
          writer, sheet_name="Meth_Subgroups", index=False
      )
      expr_results.to_excel(writer, sheet_name="Expr_Results", index=False)
      meth_results.to_excel(writer, sheet_name="Meth_Results", index=False)

    print(f"Finished {len(expr_defs)} expression and {len(meth_defs)} methylation subgroup analyses")


if __name__ == "__main__":
    main()
