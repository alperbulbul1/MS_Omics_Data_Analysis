import argparse
import gzip
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


TARGETS_PATH = "__MS_GEO_ROOT__/Methylation_Target_Datasets.csv"
BASE_DIR = "__MS_GEO_ROOT__/Methylation_Data"
MIN_FEATURE_OVERLAP = 5000


def sanitize_label(value: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9]+", "_", value.strip().lower()).strip("_")
    return clean or "unspecified"


def classify_condition(text: str) -> str:
    text = f" {text.lower()} "
    hc_terms = [
        " healthy ",
        " healthy donor ",
        " healthy subject ",
        " healthy control ",
        " control ",
        " hc ",
        " diagnosis: h ",
        " diagnosis: healthy ",
        " non-inflammatory ",
        " non ms ",
    ]
    ms_terms = [
        " multiple sclerosis ",
        " ms patient ",
        " ms patients ",
        " disease status: multiple sclerosis ",
        " rrms ",
        " spms ",
        " relapsing-remitting ",
        " relapsing remitting ",
        " secondary progressive ",
        " primary progressive ",
        " case ",
    ]

    is_hc = any(term in text for term in hc_terms)
    is_ms = any(term in text for term in ms_terms)

    if is_hc and not is_ms:
        return "HC"
    if is_ms and not is_hc:
        return "MS"
    return "Unknown"


def detect_cell_type(text: str) -> str:
    text = text.lower()
    mapping = [
        ("bronchoalveolar_lavage", ["bronchoalveolar lavage", "bal cells", "alveolar macrophages"]),
        ("vitd3_toldc", ["vitd3-toldc", "vitd3 toldc"]),
        ("toldc", ["toldc", "tolerogenic dendritic"]),
        ("mature_dendritic", [" mature dendritic ", " mdc "]),
        ("cd19_bcell", ["cd19", "b cell", "b-cell"]),
        ("cd14_monocyte", ["cd14", "monocyte"]),
        ("cd8_tcell", ["cd8+", " cd8 t", "cd8 t cell"]),
        ("cd4_tcell", ["cd4+", " cd4 t", "cd4 t cell"]),
        ("whole_blood", ["whole blood", "peripheral blood", "blood"]),
        ("pbmc", ["pbmc", "mononuclear"]),
    ]
    padded = f" {text} "
    for label, terms in mapping:
        if any(term in padded for term in terms):
            return label
    return "unspecified"


def parse_soft_metadata(soft_path: str) -> dict[str, dict[str, str]]:
    metadata: dict[str, dict[str, str]] = {}
    current_sample = None
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
                current_sample = line.split("=", 1)[1].strip()
                metadata[current_sample] = {"parts": []}
                continue

            if current_sample is None:
                continue

            if line.startswith("^"):
                current_sample = None
                continue

            for prefix in prefixes:
                if line.startswith(prefix):
                    metadata[current_sample]["parts"].append(line.split("=", 1)[1].strip())
                    break

    for sample_id, sample_meta in metadata.items():
        text = " ".join(sample_meta["parts"])
        metadata[sample_id] = {
            "sample_id": sample_id,
            "condition": classify_condition(text),
            "cell_type": detect_cell_type(text),
            "raw_text": text,
        }
    return metadata


def read_series_matrix(matrix_path: str) -> pd.DataFrame:
    df = pd.read_csv(
        matrix_path,
        sep="\t",
        comment="!",
        compression="gzip",
        low_memory=False,
    )

    if df.empty:
        raise ValueError("Series matrix is empty")

    df.columns = [str(col).strip().strip('"') for col in df.columns]
    first_col = df.columns[0]
    df[first_col] = df[first_col].astype(str).str.strip().str.strip('"')
    df = df.set_index(first_col)
    df = df.loc[~df.index.duplicated(keep="first")]
    df = df[[col for col in df.columns if not col.startswith("Unnamed:")]]
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(axis=0, how="all").dropna(axis=1, how="all")
    return df


def read_soft_sample_tables(soft_path: str, valid_samples: set[str]) -> pd.DataFrame:
    tables: list[pd.DataFrame] = []
    current_sample = None
    in_table = False
    header = None
    rows = []

    def flush_table(sample_id: str, table_header: list[str] | None, table_rows: list[list[str]]) -> None:
        if sample_id not in valid_samples or not table_header or not table_rows:
            return

        table = pd.DataFrame(table_rows, columns=table_header)
        table.columns = [str(col).strip() for col in table.columns]
        id_candidates = [col for col in table.columns if col.lower() in {"id_ref", "id", "probe", "cpg", "ilmnid"}]
        if not id_candidates:
            return
        id_col = id_candidates[0]

        value_candidates = [
            col
            for col in table.columns
            if col.lower() in {"value", "beta_value", "beta", "m_value", "methylation_level", "signal"}
        ]
        if not value_candidates:
            value_candidates = [col for col in table.columns if col != id_col]
        if not value_candidates:
            return

        value_col = value_candidates[-1]
        sub = table[[id_col, value_col]].copy()
        sub.columns = ["Probe", sample_id]
        sub[sample_id] = pd.to_numeric(sub[sample_id], errors="coerce")
        sub["Probe"] = sub["Probe"].astype(str).str.strip().str.strip('"')
        sub = sub.dropna(subset=[sample_id]).drop_duplicates(subset=["Probe"])
        if not sub.empty:
            tables.append(sub.set_index("Probe"))

    with gzip.open(soft_path, "rt", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")

            if line.startswith("^SAMPLE = "):
                flush_table(current_sample, header, rows)
                current_sample = line.split("=", 1)[1].strip()
                in_table = False
                header = None
                rows = []
                continue

            lowered = line.lower().strip()
            if lowered == "!sample_table_begin":
                in_table = True
                header = None
                rows = []
                continue

            if lowered == "!sample_table_end":
                flush_table(current_sample, header, rows)
                in_table = False
                header = None
                rows = []
                continue

            if not in_table:
                continue

            parts = line.split("\t")
            if header is None:
                header = parts
            else:
                if len(parts) == len(header):
                    rows.append(parts)

    flush_table(current_sample, header, rows)

    if not tables:
        return pd.DataFrame()

    combined = pd.concat(tables, axis=1, join="inner")
    combined = combined.loc[~combined.index.duplicated(keep="first")]
    combined = combined.dropna(axis=0, how="all").dropna(axis=1, how="all")
    return combined


def read_sample_supplementary_tables(gse_dir: str, metadata: pd.DataFrame) -> pd.DataFrame:
    sample_to_condition = metadata.set_index("sample_id")["condition"].to_dict()
    sample_files = sorted(
        [
            path
            for path in Path(gse_dir).glob("*.gz")
            if path.name.endswith((".txt.gz", ".tsv.gz", ".csv.gz"))
            and "_family.soft.gz" not in path.name
            and "_series_matrix.txt.gz" not in path.name
        ]
    )
    if not sample_files:
        return pd.DataFrame()

    tables: list[pd.DataFrame] = []
    for path in sample_files:
        sample_match = re.search(r"(GSM\d+)", path.name)
        if not sample_match:
            continue
        sample_id = sample_match.group(1)
        if sample_to_condition.get(sample_id) not in {"MS", "HC"}:
            continue

        try:
            df = pd.read_csv(path, sep="\t", compression="gzip", dtype=str, low_memory=False)
        except Exception:
            try:
                df = pd.read_csv(path, sep=",", compression="gzip", dtype=str, low_memory=False)
            except Exception:
                continue

        if df.empty or df.shape[1] < 2:
            continue

        numeric = df.apply(pd.to_numeric, errors="coerce")
        numeric_cols = [col for col in numeric.columns if numeric[col].notna().sum() > 0]
        if not numeric_cols:
            continue
        value_col = numeric_cols[-1]

        non_numeric_cols = [col for col in df.columns if col != value_col]
        if len(non_numeric_cols) >= 3:
            feature = (
                df[non_numeric_cols[:3]]
                .fillna("")
                .astype(str)
                .agg(":".join, axis=1)
                .str.replace(r":+$", "", regex=True)
            )
        else:
            feature = df[non_numeric_cols[0]].fillna("").astype(str)

        sub = pd.DataFrame({"Probe": feature, sample_id: pd.to_numeric(df[value_col], errors="coerce")})
        sub = sub.dropna(subset=[sample_id])
        sub["Probe"] = sub["Probe"].astype(str).str.strip()
        sub = sub[sub["Probe"] != ""].drop_duplicates(subset=["Probe"])
        if not sub.empty:
            tables.append(sub.set_index("Probe"))

    if not tables:
        return pd.DataFrame()

    combined = pd.concat(tables, axis=1, join="inner")
    combined = combined.loc[~combined.index.duplicated(keep="first")]
    combined = combined.dropna(axis=0, how="all").dropna(axis=1, how="all")
    return combined


def beta_to_m_values(df: pd.DataFrame, eps: float = 1e-6) -> pd.DataFrame:
    clipped = df.clip(lower=eps, upper=1 - eps)
    return np.log2(clipped / (1 - clipped))


def prepare_subdataset(
    gse_id: str,
    matrix: pd.DataFrame,
    meta_df: pd.DataFrame,
    cell_type: str,
) -> tuple[str, pd.DataFrame, pd.DataFrame] | None:
    subset_meta = meta_df[meta_df["cell_type"] == cell_type].copy()
    subset_meta = subset_meta[subset_meta["condition"].isin(["MS", "HC"])]

    counts = subset_meta["condition"].value_counts()
    if counts.get("MS", 0) < 2 or counts.get("HC", 0) < 2:
        return None

    valid_samples = [sample for sample in subset_meta["sample_id"] if sample in matrix.columns]
    if len(valid_samples) < 4:
        return None

    subset_meta = subset_meta[subset_meta["sample_id"].isin(valid_samples)].copy()
    subset_meta["dataset"] = f"{gse_id}__{sanitize_label(cell_type)}"

    subset_matrix = matrix[valid_samples].copy()
    subset_matrix = subset_matrix.dropna(axis=0, how="all")
    subset_matrix = subset_matrix.loc[subset_matrix.notna().sum(axis=1) >= max(4, int(0.8 * len(valid_samples)))]
    if subset_matrix.empty:
        return None

    return subset_meta["dataset"].iloc[0], subset_matrix, subset_meta


def load_dataset(gse_id: str) -> list[tuple[str, pd.DataFrame, pd.DataFrame]]:
    gse_dir = os.path.join(BASE_DIR, gse_id)
    soft_path = os.path.join(gse_dir, f"{gse_id}_family.soft.gz")
    matrix_path = os.path.join(gse_dir, f"{gse_id}_series_matrix.txt.gz")

    if not os.path.exists(soft_path):
        print(f"Skipping {gse_id}: missing soft file")
        return []

    try:
        metadata = parse_soft_metadata(soft_path)
    except Exception as exc:
        print(f"Skipping {gse_id}: failed to parse soft metadata ({exc})")
        return []
    meta_df = pd.DataFrame(metadata.values())
    if meta_df.empty:
        print(f"Skipping {gse_id}: no sample metadata")
        return []

    valid_sample_ids = set(meta_df["sample_id"])
    matrix = pd.DataFrame()
    matrix_source = "none"
    if os.path.exists(matrix_path):
        try:
            matrix = read_series_matrix(matrix_path)
            matrix_source = "series_matrix"
        except Exception as exc:
            print(f"{gse_id}: series matrix parse failed ({exc}), trying soft sample tables")

    if matrix.empty:
        try:
            matrix = read_soft_sample_tables(soft_path, valid_sample_ids)
            if not matrix.empty:
                matrix_source = "soft_sample_table"
        except Exception as exc:
            print(f"Skipping {gse_id}: soft sample table parse failed ({exc})")
            return []

    if matrix.empty:
        try:
            matrix = read_sample_supplementary_tables(gse_dir, meta_df)
            if not matrix.empty:
                matrix_source = "sample_supplementary"
        except Exception as exc:
            print(f"Skipping {gse_id}: sample supplementary parse failed ({exc})")
            return []

    if matrix.empty:
        print(f"Skipping {gse_id}: matrix empty")
        return []

    matrix.columns = [str(col).strip().strip('"') for col in matrix.columns]
    valid_meta = meta_df[meta_df["sample_id"].isin(matrix.columns)].copy()
    valid_meta = valid_meta[valid_meta["condition"].isin(["MS", "HC"])]
    if valid_meta.empty:
        print(f"Skipping {gse_id}: no case/control samples matched to matrix")
        return []

    if matrix_source == "sample_supplementary":
        expected_case_control = meta_df[meta_df["condition"].isin(["MS", "HC"])]["sample_id"].nunique()
        observed = valid_meta["sample_id"].nunique()
        if expected_case_control and observed < max(4, int(np.ceil(expected_case_control * 0.8))):
            print(
                f"Skipping {gse_id}: sample supplementary download incomplete "
                f"({observed}/{expected_case_control} case-control sample files present)"
            )
            return []

    matrix = matrix[[col for col in matrix.columns if col in set(valid_meta["sample_id"])]]
    matrix = matrix.loc[~matrix.index.duplicated(keep="first")]

    sample_min = np.nanmin(matrix.values)
    sample_max = np.nanmax(matrix.values)
    if sample_min >= -1e-6 and sample_max <= 1.000001:
        matrix = beta_to_m_values(matrix)

    cell_types = sorted(valid_meta["cell_type"].fillna("unspecified").unique())
    outputs: list[tuple[str, pd.DataFrame, pd.DataFrame]] = []

    if len(cell_types) > 1:
        for cell_type in cell_types:
            prepared = prepare_subdataset(gse_id, matrix, valid_meta, cell_type)
            if prepared is not None:
                outputs.append(prepared)
    else:
        valid_meta = valid_meta.copy()
        valid_meta["dataset"] = gse_id
        counts = valid_meta["condition"].value_counts()
        if counts.get("MS", 0) >= 2 and counts.get("HC", 0) >= 2:
            valid_samples = [sample for sample in valid_meta["sample_id"] if sample in matrix.columns]
            sub_matrix = matrix[valid_samples].copy()
            sub_matrix = sub_matrix.loc[sub_matrix.notna().sum(axis=1) >= max(4, int(0.8 * len(valid_samples)))]
            outputs.append((gse_id, sub_matrix, valid_meta[valid_meta["sample_id"].isin(valid_samples)]))

    if outputs:
        print(f"Loaded {gse_id}: {len(outputs)} usable cohort(s)")
    else:
        print(f"Skipping {gse_id}: no usable cohort after sample filtering")
    return outputs


def select_common_features(
    matrices: dict[str, pd.DataFrame],
    min_fraction: float,
) -> pd.Index:
    feature_counts = pd.Series(dtype=float)
    for matrix in matrices.values():
        feature_counts = feature_counts.add(pd.Series(1, index=matrix.index), fill_value=0)

    threshold = max(1, int(np.ceil(len(matrices) * min_fraction)))
    common = feature_counts[feature_counts >= threshold].index
    return pd.Index(common)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("gse_ids", nargs="*")
    parser.add_argument("--min-datasets-fraction", type=float, default=0.75)
    parser.add_argument("--max-probes", type=int, default=450000)
    args = parser.parse_args()

    targets = pd.read_csv(TARGETS_PATH)
    gse_ids = sorted(targets["gse_id"].dropna().astype(str).unique())
    if args.gse_ids:
        allowed = set(args.gse_ids)
        gse_ids = [gse_id for gse_id in gse_ids if gse_id in allowed]

    matrices: dict[str, pd.DataFrame] = {}
    metadata_frames: list[pd.DataFrame] = []

    for gse_id in gse_ids:
        for dataset_id, matrix, meta_df in load_dataset(gse_id):
            if matrix.empty:
                continue
            matrices[dataset_id] = matrix
            metadata_frames.append(meta_df)

    if not matrices:
        raise SystemExit("No methylation datasets were successfully loaded")

    common_features = select_common_features(matrices, args.min_datasets_fraction)
    if common_features.empty:
        raise SystemExit("No common methylation probes found across datasets")

    compatible_matrices = {
        dataset_id: matrix
        for dataset_id, matrix in matrices.items()
        if matrix.index.intersection(common_features).size >= min(MIN_FEATURE_OVERLAP, len(common_features))
    }
    dropped = sorted(set(matrices) - set(compatible_matrices))
    if dropped:
        print(f"Dropping feature-incompatible datasets: {', '.join(dropped)}")
        matrices = compatible_matrices
        if not matrices:
            raise SystemExit("All datasets were feature-incompatible after overlap filtering")
        common_features = select_common_features(matrices, args.min_datasets_fraction)
        if common_features.empty:
            raise SystemExit("No common methylation probes remained after overlap filtering")

    harmonized = []
    for dataset_id, matrix in matrices.items():
        sub = matrix.reindex(common_features)
        harmonized.append(sub)

    final_matrix = pd.concat(harmonized, axis=1)

    # Global sample-level overlap check: keep only probes with real (non-NA)
    # data in >= 80% of ALL samples across all datasets combined.
    n_total_samples = final_matrix.shape[1]
    real_coverage = final_matrix.notna().sum(axis=1)
    min_samples = max(4, int(0.80 * n_total_samples))
    probes_before = len(final_matrix)
    final_matrix = final_matrix.loc[real_coverage >= min_samples]
    print(
        f"Global sample-level overlap filter: {probes_before} → {len(final_matrix)} probes "
        f"(kept probes with real data in ≥80% of {n_total_samples} samples)"
    )

    # Impute remaining sporadic NAs with per-probe row mean (vectorised)
    row_means = final_matrix.mean(axis=1)
    final_matrix = final_matrix.T.fillna(row_means).T

    if args.max_probes and len(final_matrix.index) > args.max_probes:
        variances = final_matrix.var(axis=1).sort_values(ascending=False)
        keep = variances.head(args.max_probes).index
        final_matrix = final_matrix.loc[keep]
        print(f"Retained top {len(keep)} variable probes for tractable downstream modeling")

    final_meta = pd.concat(metadata_frames, ignore_index=True)
    final_meta = final_meta[final_meta["sample_id"].isin(final_matrix.columns)].copy()

    final_matrix = final_matrix.astype(np.float32)
    final_matrix.to_csv(os.path.join(BASE_DIR, "Combined_Methylation_Pre_Batch.csv"))
    final_meta.to_csv(os.path.join(BASE_DIR, "Combined_Methylation_Metadata.csv"), index=False)

    print(
        "Saved harmonized methylation matrix:",
        final_matrix.shape,
        "samples:",
        len(final_meta),
        "datasets:",
        final_meta["dataset"].nunique(),
    )


if __name__ == "__main__":
    main()
