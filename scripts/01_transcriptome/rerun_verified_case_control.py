import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import correct_and_normalize

ROOT = Path("__MS_GEO_ROOT__")
WORKBOOK = ROOT / "Used_Methylation_and_RNAseq_Dataset_Inventory.xlsx"

EXPR_DIR = ROOT / "Expression_Data"
EXPR_META = EXPR_DIR / "Combined_Metadata.csv"
EXPR_PRE = EXPR_DIR / "Combined_Expression_Pre_ComBat.csv"

METH_DIR = ROOT / "Methylation_Data"
METH_META = METH_DIR / "Combined_Methylation_Metadata.csv"
METH_PRE = METH_DIR / "Combined_Methylation_Pre_Batch.csv"


def load_verified_datasets(sheet_name: str) -> pd.DataFrame:
    df = pd.read_excel(WORKBOOK, sheet_name=sheet_name, header=1)
    status_col = "Case / Control Usage Status"
    keep = (
        df["Used MS / Case"].fillna(0).astype(int).ge(2)
        & df["Used HC / Control"].fillna(0).astype(int).ge(2)
        & ~df[status_col].fillna("").str.contains("nonstandard", case=False, regex=False)
    )
    verified = df.loc[keep].copy()
    verified["GSE ID"] = verified["GSE ID"].astype(str).str.strip()
    return verified


def filter_expression() -> None:
    verified = load_verified_datasets("Transcriptome_Used")
    keep_datasets = set(verified["GSE ID"])

    meta = pd.read_csv(EXPR_META)
    expr = pd.read_csv(EXPR_PRE, index_col=0)

    meta = meta[meta["condition"].isin(["MS", "HC"])].copy()
    meta["base_gse"] = meta["dataset"].astype(str).str.split("__", n=1).str[0]
    meta = meta[meta["base_gse"].isin(keep_datasets)].copy()
    meta = meta[meta["sample_id"].isin(expr.columns)].copy()
    meta = meta.drop_duplicates(subset="sample_id", keep="first")

    sample_order = [sample for sample in meta["sample_id"] if sample in expr.columns]
    meta = meta.set_index("sample_id").loc[sample_order].reset_index()
    expr = expr[sample_order].copy()
    expr = expr.loc[expr.var(axis=1) > 0]

    expr.to_csv(EXPR_PRE)
    meta.drop(columns=["base_gse"]).to_csv(EXPR_META, index=False)
    verified.to_csv(ROOT / "Verified_Expression_Datasets.csv", index=False)

    print(
        f"Expression rerun input ready: {len(meta)} samples, "
        f"{meta['dataset'].nunique()} datasets, {expr.shape[0]} genes"
    )


def filter_methylation() -> None:
    verified = load_verified_datasets("Methylation_Used")
    keep_datasets = set(verified["GSE ID"])

    meta = pd.read_csv(METH_META)
    mat = pd.read_csv(METH_PRE, index_col=0)

    meta = meta[meta["condition"].isin(["MS", "HC"])].copy()
    meta["base_gse"] = meta["dataset"].astype(str).str.split("__", n=1).str[0]
    meta = meta[meta["base_gse"].isin(keep_datasets)].copy()
    meta = meta[meta["sample_id"].isin(mat.columns)].copy()
    meta = meta.drop_duplicates(subset="sample_id", keep="first")

    sample_order = [sample for sample in meta["sample_id"] if sample in mat.columns]
    meta = meta.set_index("sample_id").loc[sample_order].reset_index()
    mat = mat[sample_order].copy()
    mat = mat.loc[mat.var(axis=1) > 0]

    mat.to_csv(METH_PRE)
    meta.drop(columns=["base_gse"]).to_csv(METH_META, index=False)
    verified.to_csv(ROOT / "Verified_Methylation_Datasets.csv", index=False)

    print(
        f"Methylation rerun input ready: {len(meta)} samples, "
        f"{meta['dataset'].nunique()} datasets, {mat.shape[0]} probes"
    )


def main() -> None:
    filter_expression()
    print("\n" + "=" * 60)
    print("Running correct_and_normalize pipeline on filtered expression data...")
    print("=" * 60)
    correct_and_normalize.run()
    filter_methylation()


if __name__ == "__main__":
    main()
