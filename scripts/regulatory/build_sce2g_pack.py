"""Build the bundled scE2G regulatory-links pack: download each public scE2G_multiome cell-type file, compact it
(compact_sce2g.py), delete the raw file, and write regulatory_links/sce2g/manifest.json.

Source: scE2G, Sheth, Qiu et al. (Engreitz lab), https://data.broadinstitute.org/scE2Gpredictions/
    PBMC_9_Cell_Types (9 cell types), BMMC_22_Cell_Types (22), Pancreatic_Islets (9); the "merged" files are skipped.

    python scripts/regulatory/build_sce2g_pack.py <hg38ToHg19.over.chain.gz> [--work DIR]
Re-running skips cell types already compacted.
"""
import gzip
import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import compact_sce2g  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "regulatory_links" / "sce2g"
BASE = "https://data.broadinstitute.org/scE2Gpredictions/"

# (folder, file cell-type token, label, tissue group, lineage)
CELL_TYPES = [
    ("PBMC_9_Cell_Types", "17_B", "B (PBMC)", "PBMC", "B"),
    ("PBMC_9_Cell_Types", "10_cDC", "cDC (PBMC)", "PBMC", "Myeloid"),
    ("PBMC_9_Cell_Types", "11_CD14.Mono.1", "CD14 Mono1 (PBMC)", "PBMC", "Myeloid"),
    ("PBMC_9_Cell_Types", "12_CD14.Mono.2", "CD14 Mono2 (PBMC)", "PBMC", "Myeloid"),
    ("PBMC_9_Cell_Types", "13_CD16.Mono", "CD16 Mono (PBMC)", "PBMC", "Myeloid"),
    ("PBMC_9_Cell_Types", "20_CD4.N1", "CD4 Naive (PBMC)", "PBMC", "T/NK"),
    ("PBMC_9_Cell_Types", "22_CD4.M", "CD4 Memory (PBMC)", "PBMC", "T/NK"),
    ("PBMC_9_Cell_Types", "24_CD8.CM", "CD8 CM (PBMC)", "PBMC", "T/NK"),
    ("PBMC_9_Cell_Types", "25_NK", "NK (PBMC)", "PBMC", "T/NK"),
    ("BMMC_22_Cell_Types", "B1_B", "B1 B (BM)", "Bone marrow", "B"),
    ("BMMC_22_Cell_Types", "Naive_CD20_pos_B", "Naive B (BM)", "Bone marrow", "B"),
    ("BMMC_22_Cell_Types", "Transitional_B", "Transitional B (BM)", "Bone marrow", "B"),
    ("BMMC_22_Cell_Types", "Plasma_cell", "Plasma cell (BM)", "Bone marrow", "B"),
    ("BMMC_22_Cell_Types", "CD4_pos_T_activated", "CD4 T activated (BM)", "Bone marrow", "T/NK"),
    ("BMMC_22_Cell_Types", "CD4_pos_T_naive", "CD4 T naive (BM)", "Bone marrow", "T/NK"),
    ("BMMC_22_Cell_Types", "CD8_pos_T", "CD8 T (BM)", "Bone marrow", "T/NK"),
    ("BMMC_22_Cell_Types", "CD8_pos_T_naive", "CD8 T naive (BM)", "Bone marrow", "T/NK"),
    ("BMMC_22_Cell_Types", "NK", "NK (BM)", "Bone marrow", "T/NK"),
    ("BMMC_22_Cell_Types", "ILC", "ILC (BM)", "Bone marrow", "T/NK"),
    ("BMMC_22_Cell_Types", "CD14_pos_Mono", "CD14 Mono (BM)", "Bone marrow", "Myeloid"),
    ("BMMC_22_Cell_Types", "CD16_pos_Mono", "CD16 Mono (BM)", "Bone marrow", "Myeloid"),
    ("BMMC_22_Cell_Types", "cDC2", "cDC2 (BM)", "Bone marrow", "Myeloid"),
    ("BMMC_22_Cell_Types", "pDC", "pDC (BM)", "Bone marrow", "Myeloid"),
    ("BMMC_22_Cell_Types", "HSC", "HSC (BM)", "Bone marrow", "Progenitor"),
    ("BMMC_22_Cell_Types", "Lymph_prog", "Lymph prog (BM)", "Bone marrow", "Progenitor"),
    ("BMMC_22_Cell_Types", "G_M_prog", "GM prog (BM)", "Bone marrow", "Progenitor"),
    ("BMMC_22_Cell_Types", "ID2_hi_myeloid_prog", "Myeloid prog (BM)", "Bone marrow", "Progenitor"),
    ("BMMC_22_Cell_Types", "MK_E_prog", "MK E prog (BM)", "Bone marrow", "Progenitor"),
    ("BMMC_22_Cell_Types", "Proerythroblast", "Proerythroblast (BM)", "Bone marrow", "Erythroid"),
    ("BMMC_22_Cell_Types", "Erythroblast", "Erythroblast (BM)", "Bone marrow", "Erythroid"),
    ("BMMC_22_Cell_Types", "Normoblast", "Normoblast (BM)", "Bone marrow", "Erythroid"),
    ("Pancreatic_Islets", "beta", "Beta (islets)", "Pancreatic islets", "Endocrine"),
    ("Pancreatic_Islets", "alpha", "Alpha (islets)", "Pancreatic islets", "Endocrine"),
    ("Pancreatic_Islets", "delta", "Delta (islets)", "Pancreatic islets", "Endocrine"),
    ("Pancreatic_Islets", "PP", "PP (islets)", "Pancreatic islets", "Endocrine"),
    ("Pancreatic_Islets", "acinar", "Acinar (islets)", "Pancreatic islets", "Exocrine"),
    ("Pancreatic_Islets", "duct", "Duct (islets)", "Pancreatic islets", "Exocrine"),
    ("Pancreatic_Islets", "endothelial", "Endothelial (islets)", "Pancreatic islets", "Other"),
    ("Pancreatic_Islets", "immune", "Immune (islets)", "Pancreatic islets", "Other"),
    ("Pancreatic_Islets", "mesenchymal", "Mesenchymal (islets)", "Pancreatic islets", "Other"),
]


def file_id(label):
    return "".join(c if c.isalnum() else "_" for c in label).strip("_").replace("__", "_")


def count_rows(path):
    with gzip.open(path, "rt") as f:
        return sum(1 for ln in f if not ln.startswith("#")) - 1


def main():
    chain = sys.argv[1]
    work = Path(sys.argv[sys.argv.index("--work") + 1]) if "--work" in sys.argv else ROOT.parent / "LYNXgwas-regulatory" / "raw"
    work.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    entries = []
    for folder, token, label, group, lineage in CELL_TYPES:
        fid = file_id(label)
        out = OUT / f"{fid}.tsv.gz"
        url = f"{BASE}{folder}/scE2G/pairs.E2G.scE2G_multiome.res.{token}.tsv.gz"
        if not out.is_file():
            raw = work / f"{fid}.raw.tsv.gz"
            t0 = time.time()
            for attempt in range(1, 4):
                try:
                    urllib.request.urlretrieve(url, raw)
                    break
                except Exception as e:
                    print(f"  download failed ({e}); retry {attempt}", flush=True)
                    time.sleep(20 * attempt)
            else:
                sys.exit(f"could not download {url}")
            print(f"{label}: downloaded {raw.stat().st_size / 1e6:.0f} MB in {time.time() - t0:.0f} s", flush=True)
            compact_sce2g.main(str(raw), label, chain, str(out))
            raw.unlink()
        md5 = hashlib.md5(out.read_bytes()).hexdigest()
        entries.append({"id": fid, "label": label, "group": group, "lineage": lineage, "file": out.name,
                        "pairs": count_rows(out), "bytes": out.stat().st_size, "md5": md5, "source": url})
        print(f"  {label}: {entries[-1]['pairs']:,} pairs, {entries[-1]['bytes'] / 1e6:.2f} MB", flush=True)
    manifest = {
        "resource": "scE2G",
        "title": "scE2G enhancer-gene links (scE2G_multiome)",
        "citation": "Sheth MU, Qiu W-L, et al. Mapping enhancer-gene regulatory interactions from single-cell data. "
                    "bioRxiv 2024.11.23.624931 (Engreitz lab). Predictions: https://data.broadinstitute.org/scE2Gpredictions/",
        "threshold": "E2G.Score.qnorm >= 0.177 (scE2G_multiome recall threshold)",
        "builds": "GRCh37 (lifted with the UCSC hg38ToHg19 chain) and the original GRCh38",
        "columns": compact_sce2g.KEEP_OUT,
        "cell_types": entries,
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=1))
    print(f"ALL DONE: {len(entries)} cell types, {sum(e['bytes'] for e in entries) / 1e6:.1f} MB", flush=True)


if __name__ == "__main__":
    main()
