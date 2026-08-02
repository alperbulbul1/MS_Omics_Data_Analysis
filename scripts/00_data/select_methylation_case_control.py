import pandas as pd


INPUT_PATH = "__MS_GEO_ROOT__/Filtered_Naive_Methylation_MS_Datasets_with_Tissues.xlsx"
OUTPUT_PATH = "__MS_GEO_ROOT__/Methylation_Target_Datasets.csv"
ADDITIONAL_PATH = "__MS_GEO_ROOT__/Additional_Methylation_Candidates.csv"


def contains_any(text: str, terms: list[str]) -> bool:
    return any(term in text for term in terms)


def main() -> None:
    df = pd.read_excel(INPUT_PATH)

    organism = df["Organism"].fillna("").astype(str).str.lower()
    assay = df["Assay Type"].fillna("").astype(str).str.lower()
    text = (
        df["Accession Title"].fillna("").astype(str)
        + " "
        + df["Sample Information"].fillna("").astype(str)
        + " "
        + df.get("Tissue / Source Information", pd.Series("", index=df.index)).fillna("").astype(str)
    ).str.lower()

    immune_terms = [
        "blood",
        "pbmc",
        "mononuclear",
        "b cell",
        "t cell",
        "monocyte",
        "cd4",
        "cd8",
        "cd14",
        "cd19",
        "pbl",
        "lymphocyte",
        "dendritic",
    ]
    assay_terms = [
        "humanmethylation450",
        "methylation450",
        "methylationepic",
        "infinium methylationepic",
    ]
    control_terms = ["healthy", "control", "hc", "healthy donor", "healthy subject"]
    disease_terms = ["multiple sclerosis", " ms ", "ms patients", "ms patient", "rrms", "spms"]

    selected = []
    for _, row in df.iterrows():
        gse_id = str(row.get("GEO IDs", "")).strip()
        if not gse_id or gse_id == "nan":
            continue

        row_text = text.loc[row.name]
        row_assay = assay.loc[row.name]
        row_organism = organism.loc[row.name]

        is_human = "human" in row_organism or "sapiens" in row_organism
        is_array_methylation = contains_any(row_assay, assay_terms)
        is_immune = contains_any(row_text, immune_terms)
        has_case_control_signal = contains_any(row_text, control_terms) and contains_any(row_text, disease_terms)

        if not (is_human and is_array_methylation and is_immune and has_case_control_signal):
            continue

        selected.append(
            {
                "gse_id": gse_id,
                "assay_type": row.get("Assay Type", ""),
                "title": row.get("Accession Title", ""),
                "sample_size": row.get("Sample Size", ""),
            }
        )

    out_df = pd.DataFrame(selected).drop_duplicates(subset=["gse_id"]).sort_values("gse_id")

    if pd.io.common.file_exists(ADDITIONAL_PATH):
        extra = pd.read_csv(ADDITIONAL_PATH)
        extra = extra.rename(columns={"technology": "assay_type", "title": "title"})
        extra = extra.assign(sample_size=pd.NA)
        extra = extra[["gse_id", "assay_type", "title", "sample_size"]]
        out_df = (
            pd.concat([out_df, extra], ignore_index=True)
            .drop_duplicates(subset=["gse_id"], keep="first")
            .sort_values("gse_id")
            .reset_index(drop=True)
        )

    out_df.to_csv(OUTPUT_PATH, index=False)

    print(f"Saved {len(out_df)} methylation case/control targets to {OUTPUT_PATH}")
    if not out_df.empty:
        print(out_df.to_string(index=False))


if __name__ == "__main__":
    main()
