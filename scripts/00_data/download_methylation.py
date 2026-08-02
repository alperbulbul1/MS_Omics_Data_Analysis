import argparse
import gzip
import os
import urllib.parse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import requests
from bs4 import BeautifulSoup


TARGETS_PATH = "__MS_GEO_ROOT__/Methylation_Target_Datasets.csv"
DEST_DIR = "__MS_GEO_ROOT__/Methylation_Data"
USED_META_PATH = "__MS_GEO_ROOT__/Methylation_Data/Combined_Methylation_Metadata.csv"


def is_complete_file(path: str) -> bool:
    if not os.path.exists(path):
        return False
    if os.path.getsize(path) == 0:
        return False
    if path.endswith(".gz"):
        try:
            with gzip.open(path, "rb") as handle:
                while handle.read(1024 * 1024):
                    pass
            return True
        except Exception:
            return False
    return True


def extract_sample_supplementary_urls(soft_path: str) -> dict[str, list[str]]:
    urls: dict[str, list[str]] = {}
    if not os.path.exists(soft_path):
        return urls

    current_sample = None
    try:
        with gzip.open(soft_path, "rt", encoding="utf-8", errors="ignore") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if line.startswith("^SAMPLE = "):
                    current_sample = line.split("=", 1)[1].strip()
                    urls.setdefault(current_sample, [])
                    continue
                if line.startswith("!Sample_supplementary_file"):
                    url = line.split("=", 1)[1].strip()
                    if current_sample and (url.startswith("ftp://") or url.startswith("https://")):
                        urls.setdefault(current_sample, []).append(url)
    except Exception as exc:
        print(f"Could not parse sample supplementary URLs from {soft_path}: {exc}")
    return {sample_id: sorted(set(sample_urls)) for sample_id, sample_urls in urls.items()}


def load_used_sample_ids() -> dict[str, set[str]]:
    if not os.path.exists(USED_META_PATH):
        return {}
    meta = pd.read_csv(USED_META_PATH)
    meta["dataset"] = meta["dataset"].astype(str)
    meta["base_dataset"] = meta["dataset"].str.split("__").str[0]
    used: dict[str, set[str]] = {}
    for base_dataset, subset in meta.groupby("base_dataset"):
        used[base_dataset] = set(subset["sample_id"].astype(str))
    return used


def download_file(url: str, output_path: str) -> bool:
    print(f"Downloading {url} -> {output_path}")
    try:
        response = requests.get(url, stream=True, timeout=90)
        response.raise_for_status()
        with open(output_path, "wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
        return True
    except Exception as exc:
        print(f"Failed: {url} ({exc})")
        if os.path.exists(output_path):
            os.remove(output_path)
        return False


def download_gse(gse_id: str, include_idat: bool = False, used_sample_ids: set[str] | None = None) -> None:
    nnn_dir = f"{gse_id[:-3]}nnn"
    gse_dir = os.path.join(DEST_DIR, gse_id)
    os.makedirs(gse_dir, exist_ok=True)

    soft_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{nnn_dir}/{gse_id}/soft/{gse_id}_family.soft.gz"
    soft_path = os.path.join(gse_dir, f"{gse_id}_family.soft.gz")
    if not is_complete_file(soft_path):
        download_file(soft_url, soft_path)

    matrix_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{nnn_dir}/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz"
    matrix_path = os.path.join(gse_dir, f"{gse_id}_series_matrix.txt.gz")
    if not is_complete_file(matrix_path):
        download_file(matrix_url, matrix_path)

    suppl_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{nnn_dir}/{gse_id}/suppl/"
    try:
        response = requests.get(suppl_url, timeout=45)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")
        for link in soup.find_all("a"):
            href = link.get("href")
            if not href or gse_id not in href:
                continue
            if "RAW.tar" in href:
                continue
            if not href.endswith((".txt", ".txt.gz", ".csv", ".csv.gz", ".tsv", ".tsv.gz", ".xlsx", ".xls")):
                continue

            file_url = urllib.parse.urljoin(suppl_url, href)
            file_path = os.path.join(gse_dir, href)
            if not is_complete_file(file_path):
                download_file(file_url, file_path)
    except Exception as exc:
        print(f"Supplementary listing failed for {gse_id}: {exc}")

    supplementary_by_sample = extract_sample_supplementary_urls(soft_path)
    for sample_id, sample_urls in supplementary_by_sample.items():
        if used_sample_ids is not None and sample_id not in used_sample_ids:
            continue
        for file_url in sample_urls:
            filename = os.path.basename(urllib.parse.urlparse(file_url).path)
            allowed_suffixes = (".txt", ".txt.gz", ".csv", ".csv.gz", ".tsv", ".tsv.gz")
            if include_idat:
                allowed_suffixes = allowed_suffixes + (".idat", ".idat.gz")
            if not filename.endswith(allowed_suffixes):
                continue
            file_path = os.path.join(gse_dir, filename)
            if not is_complete_file(file_path):
                download_file(file_url.replace("ftp://", "https://"), file_path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("gse_ids", nargs="*")
    parser.add_argument("--max-workers", type=int, default=4)
    parser.add_argument("--include-idat", action="store_true")
    parser.add_argument("--only-used-samples", action="store_true")
    args = parser.parse_args()

    os.makedirs(DEST_DIR, exist_ok=True)
    targets = pd.read_csv(TARGETS_PATH)
    gse_ids = sorted(targets["gse_id"].dropna().astype(str).unique())
    if args.gse_ids:
        allowed = set(args.gse_ids)
        gse_ids = [gse for gse in gse_ids if gse in allowed]
    print(f"Downloading {len(gse_ids)} methylation datasets into {DEST_DIR}")

    used_samples_by_dataset = load_used_sample_ids() if args.only_used_samples else {}

    with ThreadPoolExecutor(max_workers=args.max_workers) as pool:
        futures = []
        for gse_id in gse_ids:
            futures.append(
                pool.submit(
                    download_gse,
                    gse_id,
                    args.include_idat,
                    used_samples_by_dataset.get(gse_id),
                )
            )
        for future in futures:
            future.result()


if __name__ == "__main__":
    main()
