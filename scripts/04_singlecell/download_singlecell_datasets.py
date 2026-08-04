"""
download_singlecell_datasets.py
--------------------------
Download the MS single-cell h5ad files from CELLxGENE Discover into data/.

Two modes:
  (A) If metadata/ms_datasets_cellxgene.csv exists (produced by 01_query_...py),
      iterate its rows and download each h5ad via the Discover datasets endpoint.
  (B) Fallback: use the hand-curated list in CURATED_DATASETS below (see
      README.md for provenance). Useful if the Census API is unreachable.

Downloads land in data/<collection_short>/<dataset_id>.h5ad

Requires: requests, tqdm (pip install requests tqdm)
"""

from __future__ import annotations
import os
import sys
import csv
import json
import time
import logging
from pathlib import Path

import requests
from tqdm import tqdm

# In the authors' tree this file sat in SingleCell_CELLxGENE/scripts/, so parent.parent was the
# single-cell data root. In the release it resolves to <repo>/scripts/, which would stream the
# ~18 GB single-cell tree into the git checkout and put it where no downstream script looks.
# Use the same placeholder the sibling pseudobulk scripts use.
ROOT = Path("__MS_GEO_ROOT__") / "SingleCell_CELLxGENE"
DATA = ROOT / "data"
META = ROOT / "metadata"
LOGS = ROOT / "logs"
for p in (DATA, META, LOGS):
    p.mkdir(exist_ok=True)

logging.basicConfig(
    filename=LOGS / "02_download.log",
    filemode="w",
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)

API = "https://api.cellxgene.cziscience.com/curation/v1"

# Hand-curated list of MS single-cell studies in CELLxGENE Discover.
# Collection IDs are the stable UUIDs CELLxGENE assigns to each publication.
# Each record has:
#   short      -> short name used as subfolder
#   collection -> CELLxGENE collection UUID
#   compartment-> blood, brain, CSF
#   pmid / doi / first_author / year / title   (from the peer-reviewed paper)
# Verify / refresh via: https://cellxgene.cziscience.com/collections/<collection>
CURATED_DATASETS = [
    {
        "short": "Schafflick2020_CSF_Blood",
        "collection": "2b02dff7-e427-4cdc-96fe-c36c6e41bdc0",  # Schafflick 2020
        "compartment": "blood + CSF",
        "tissue": "PBMC, CSF",
        "pmid": "31980637",
        "doi": "10.1038/s41467-019-14118-w",
        "first_author": "Schafflick D",
        "year": 2020,
        "journal": "Nature Communications",
        "title": "Integrated single cell analysis of blood and cerebrospinal fluid leukocytes in multiple sclerosis",
    },
    {
        "short": "Ostkamp2023_CSF_compendium",
        "collection": "aa5fcb2a-d8fe-4f2c-be26-7b587f4a4cf2",  # Ostkamp compendium
        "compartment": "blood + CSF",
        "tissue": "PBMC, CSF",
        "pmid": "36536441",
        "doi": "10.1186/s12974-022-02667-9",
        "first_author": "Ostkamp P",
        "year": 2022,
        "journal": "J Neuroinflammation",
        "title": "Integrated single-cell transcriptomics of CSF cells in treatment-naive MS",
    },
    {
        "short": "Schirmer2019_MS_cortex",
        "collection": "180bff9c-c8a5-4539-b13b-ddbc00d643e6",
        "compartment": "brain",
        "tissue": "cortical grey + white matter",
        "pmid": "31316211",
        "doi": "10.1038/s41586-019-1404-z",
        "first_author": "Schirmer L",
        "year": 2019,
        "journal": "Nature",
        "title": "Neuronal vulnerability and multilineage diversity in multiple sclerosis",
    },
    {
        "short": "Absinta2021_MS_lesions",
        "collection": "a72afd53-ab92-4511-88da-252fb0e26b9a",
        "compartment": "brain",
        "tissue": "subcortical white-matter lesions",
        "pmid": "34497421",
        "doi": "10.1038/s41586-021-03892-7",
        "first_author": "Absinta M",
        "year": 2021,
        "journal": "Nature",
        "title": "A lymphocyte-microglia-astrocyte axis in chronic active MS",
    },
    {
        "short": "Jakel2019_OL_heterogeneity",
        "collection": "2a79d190-a41e-4408-88c8-ac5c4d03c0fc",
        "compartment": "brain",
        "tissue": "white matter (OL lineage)",
        "pmid": "30747918",
        "doi": "10.1038/s41586-019-0903-2",
        "first_author": "Jakel S",
        "year": 2019,
        "journal": "Nature",
        "title": "Altered human oligodendrocyte heterogeneity in MS",
    },
    {
        "short": "Macnair2024_WM_glial",
        "collection": "",  # filled in by API lookup
        "compartment": "brain",
        "tissue": "white matter lesions",
        "pmid": "39657672",
        "doi": "10.1016/j.neuron.2024.11.016",
        "first_author": "Macnair W",
        "year": 2024,
        "journal": "Neuron",
        "title": "snRNA-seq stratifies MS patients into distinct white matter glial responses",
    },
    {
        "short": "Lerma-Martin2024_spatial_MS",
        "collection": "",
        "compartment": "brain",
        "tissue": "subcortical MS lesions (snRNA-seq + Visium)",
        "pmid": "39424983",
        "doi": "10.1038/s41593-024-01796-z",
        "first_author": "Lerma-Martin C",
        "year": 2024,
        "journal": "Nature Neuroscience",
        "title": "Cell type mapping reveals tissue niches and interactions in subcortical MS lesions",
    },
]


def resolve_collection(coll_uuid: str) -> dict:
    """GET /collections/{id} and return the full record (datasets + links)."""
    r = requests.get(f"{API}/collections/{coll_uuid}", timeout=30)
    r.raise_for_status()
    return r.json()


def download_one(url: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        log.info("Skip (exists): %s", dest)
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        with open(dest, "wb") as fh, tqdm(
            total=total, unit="B", unit_scale=True, desc=dest.name
        ) as bar:
            for chunk in r.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
                bar.update(len(chunk))


def main():
    # Prefer the full discovered list if present
    discovered = META / "ms_datasets_cellxgene.csv"
    if discovered.exists():
        log.info("Using discovered dataset list: %s", discovered)
        with open(discovered) as fh:
            rows = list(csv.DictReader(fh))
        for row in rows:
            dsid = row["dataset_id"]
            coll = row.get("collection_name", "unknown").replace("/", "_")[:60]
            # Discover serves h5ad at this stable path:
            url = f"https://datasets.cellxgene.cziscience.com/{dsid}.h5ad"
            dest = DATA / coll / f"{dsid}.h5ad"
            try:
                download_one(url, dest)
            except Exception as e:
                log.error("Failed %s: %s", dsid, e)
        return

    # Fallback: curated list -> resolve each collection, then download each asset
    log.info("Discovered list missing; using CURATED_DATASETS")
    manifest = []
    for rec in CURATED_DATASETS:
        if not rec["collection"]:
            log.warning("No collection UUID for %s - please fill manually", rec["short"])
            continue
        try:
            coll = resolve_collection(rec["collection"])
        except Exception as e:
            log.error("Resolve failed %s: %s", rec["short"], e)
            continue
        for ds in coll.get("datasets", []):
            dsid = ds["dataset_id"]
            # Pick the h5ad asset
            h5_url = None
            for asset in ds.get("assets", []):
                if asset.get("filetype") == "H5AD":
                    h5_url = asset.get("url")
                    break
            if not h5_url:
                h5_url = f"https://datasets.cellxgene.cziscience.com/{dsid}.h5ad"
            dest = DATA / rec["short"] / f"{dsid}.h5ad"
            manifest.append({
                "short": rec["short"], "dataset_id": dsid,
                "title": ds.get("title"), "cell_count": ds.get("cell_count"),
                "h5ad_url": h5_url, "collection_id": rec["collection"],
                "pmid": rec["pmid"], "doi": rec["doi"],
            })
            try:
                download_one(h5_url, dest)
            except Exception as e:
                log.error("Download failed %s: %s", dsid, e)
    # Write manifest
    if manifest:
        out = META / "downloaded_manifest.csv"
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(manifest[0].keys()))
            w.writeheader()
            w.writerows(manifest)
        print(f"Wrote {out}")


if __name__ == "__main__":
    main()
