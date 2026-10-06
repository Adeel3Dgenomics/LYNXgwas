"""Build the LYNXgwas example-data pack (results only) from a processed workspace.

The pack is one zip that the app's "Add example data" installer downloads (from Zenodo) and merges
into a user's workspace. It holds what is needed to browse every result, and nothing that the app can
rebuild or that we may not redistribute:

  included   project settings and metadata, the viewer data (data/), heritability, regulatory
             enrichment, the sample-size check, and for every locus the LATEST run of each tool
             (the one the app shows) without its bulky inputs; the cross-dataset runs
             (output/multi_locus/<job>) and their joint fine-mapping runs (output/mvsusie/<run>);
             the Roadmap histone peaks (regulatory_data/)
  left out   the GWAS summary statistics themselves (not ours to redistribute), LD matrices,
             reference-panel subsets and other intermediates, superseded runs

Absolute paths of the machine that built the pack are replaced by placeholders (@@LYNX_WORKSPACE@@ and
friends, see PLACEHOLDERS); pack.json lists every file that contains one, so the installer only
rewrites those. Every project is marked "results_only" in its project.json.

  python scripts/export_example_data.py --workspace <ws> --jobs <job-id> [<job-id> ...] \
      --map "<abs path>=@@LYNX_WORKSPACE@@" ... --out <dir> --version 1 [--dry-run]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import time
import zipfile
from collections import Counter
from pathlib import Path

PACK_FORMAT = 1
PLACEHOLDERS = {
    "@@LYNX_WORKSPACE@@": "the workspace the pack is installed into",
    "@@LYNX_APP@@": "the app folder (the workspace for a pip install)",
    "@@LYNX_RESOURCES@@": "<workspace>/resources (gene annotation)",
    "@@LYNX_REF_PANEL@@": "<workspace>/refpanels (reference panel, not included)",
    "@@LYNX_GWAS@@": "<workspace>/example-data/gwas (original summary statistics, not included)",
}
TEXT_EXT = {".json", ".properties", ".tsv", ".txt", ".log", ".r", ".yaml", ".yml", ".js", ".csv",
            ".md", ".key", ".snps", ".cojo", ".master", ".html"}
# Bulky inputs/intermediates inside a run folder: rebuilt by the tool whenever it is run again
RUN_EXCLUDE = re.compile(r"^(finemap\.ld|finemap\.z|ref_bigsnp\..*|susie_input\.tsv|cojo_input\.ma|ma_audit\.tsv|"
                         r"magma_snploc\.txt|magma_input\.pval|magma_step\d+\.genes\.annot\.txt|.*\.cma\.cojo|"
                         r".*\.(bed|bim|fam|bk|rds|ld))$", re.I)
MVSUSIE_KEEP = {"meta.json", "result.json", "manifest.json", "run.log", "extract.txt", "ukb_snps.tsv"}
PROJECT_TOP_SKIP = {".fingerprint_cache", "plink_subsets", "ld_results", "tmp", "plots", "loci_analysis"}
PROJECT_TOP_SKIP_RE = re.compile(r"(\.bak|\.hg18|_done)$")
# data/locus_<n>.js duplicate data/locus_<n>.json for opening the viewer straight from disk (file://);
# the app serves the .json through its API, so the copies are left out
DATA_EXCLUDE = re.compile(r"^locus_\d+\.js$")


def latest_runs(runs_dir: Path) -> list[Path]:
    """The newest run of each tool, chosen exactly as LocalServer.latestRuns does."""
    best: dict[str, tuple[int, Path]] = {}
    for d in runs_dir.iterdir() if runs_dir.is_dir() else []:
        prov, res = d / "provenance.json", d / "result.tsv"
        if not (d.is_dir() and prov.is_file() and res.is_file()):
            continue
        try:
            p = json.loads(prov.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        tool = p.get("tool", "")
        ts = p.get("timestamp") if isinstance(p.get("timestamp"), (int, float)) else int(d.stat().st_mtime * 1000)
        if tool not in best or ts > best[tool][0]:
            best[tool] = (ts, d)
    return [d for _, d in best.values()]


def select_files(ws: Path, jobs: list[str]) -> list[tuple[Path, str]]:
    """(source file, path inside the pack) for everything that goes into the pack."""
    out: list[tuple[Path, str]] = []

    def add_tree(src: Path, arc: str, keep=lambda f: True):
        for root, _, files in os.walk(src):
            for f in files:
                if keep(f):
                    p = Path(root) / f
                    out.append((p, arc + "/" + p.relative_to(src).as_posix()))

    projects = sorted(p for p in (ws / "projects").iterdir() if (p / "config.properties").is_file())
    for P in projects:
        base = f"projects/{P.name}"
        for e in sorted(P.iterdir()):
            if e.name in PROJECT_TOP_SKIP or PROJECT_TOP_SKIP_RE.search(e.name):
                continue
            if e.is_file():
                out.append((e, f"{base}/{e.name}"))
            elif e.is_dir():
                add_tree(e, f"{base}/{e.name}", keep=(lambda f: not DATA_EXCLUDE.match(f)) if e.name == "data" else (lambda f: True))
        for L in sorted((P / "loci_analysis").iterdir()):
            la = f"{base}/loci_analysis/{L.name}"
            if (L / "build.log").is_file():
                out.append((L / "build.log", f"{la}/build.log"))
            for name in ("consistency_summary.json", "consistency_report.tsv"):
                if (L / "ld" / name).is_file():
                    out.append((L / "ld" / name, f"{la}/ld/{name}"))
            if (L / "failures").is_dir():
                add_tree(L / "failures", f"{la}/failures")
            for run in latest_runs(L / "runs"):
                add_tree(run, f"{la}/runs/{run.name}", keep=lambda f: not RUN_EXCLUDE.match(f))

    for job in jobs:
        jd = ws / "output" / "multi_locus" / job
        if not (jd / "job.json").is_file():
            sys.exit(f"cross-dataset run {job} not found in {jd}")
        add_tree(jd, f"output/multi_locus/{job}")
    for md in sorted((ws / "output" / "mvsusie").iterdir()):
        meta = md / "meta.json"
        if meta.is_file() and json.loads(meta.read_text(encoding="utf-8")).get("job") in jobs:
            add_tree(md, f"output/mvsusie/{md.name}", keep=lambda f: f in MVSUSIE_KEEP)
    if (ws / "regulatory_data").is_dir():
        add_tree(ws / "regulatory_data", "regulatory_data")
    return out


def path_variants(prefix: str) -> list[str]:
    """The spellings an absolute path takes in config, JSON and logs: C:/x, C:\\x, C:\\\\x, /c/x."""
    fwd = prefix.replace("\\", "/").rstrip("/")
    drive, rest = (fwd[0], fwd[2:]) if re.match(r"^[A-Za-z]:/", fwd) else ("", fwd)
    v = {fwd, fwd.replace("/", "\\"), fwd.replace("/", "\\\\")}
    if drive:
        for d in (drive.upper(), drive.lower()):
            v |= {f"{d}:{rest}", f"{d}:{rest}".replace("/", "\\"), f"{d}:{rest}".replace("/", "\\\\"), f"/{d.lower()}{rest}"}
    return sorted(v, key=len, reverse=True)


def rewrite(text: str, rules: list[tuple[str, str]]) -> str:
    for old, new in rules:
        if old in text:
            text = text.replace(old, new)
    return text


def project_overrides(arc: str, text: str) -> str:
    """config.properties and project.json of a results-only project."""
    if arc.endswith("/config.properties"):
        # point the GFF3 at the copy the installer fetches from GENCODE, wherever it was on this machine
        text = re.sub(r"(?m)^gff3\.file=.*$", "gff3.file=@@LYNX_RESOURCES@@/gencode.v37lift37.annotation.gff3", text)
    elif arc.count("/") == 2 and arc.endswith("/project.json"):
        d = json.loads(text)
        d["results_only"] = True
        text = json.dumps(d, indent=2) + "\n"
    return text


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workspace", required=True)
    ap.add_argument("--jobs", nargs="+", required=True, help="cross-dataset runs to include")
    ap.add_argument("--map", action="append", default=[], help='"<absolute path>=<placeholder>", most specific first')
    ap.add_argument("--blank", action="append", default=[], help="machine-only paths to clear (e.g. a local R library)")
    ap.add_argument("--leak", action="append", default=[], help="strings that must not remain (e.g. the user name)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--version", required=True)
    ap.add_argument("--title", default="LYNXgwas article data: 30 GWAS datasets, 6 diseases (results only)")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    ws = Path(a.workspace).resolve()
    rules: list[tuple[str, str]] = []
    for m in a.map:
        path, ph = m.rsplit("=", 1)
        if ph not in PLACEHOLDERS:
            sys.exit(f"unknown placeholder {ph}")
        rules += [(v, ph) for v in path_variants(path)]
    for path in a.blank:
        rules += [(v, "") for v in path_variants(path)]
    rules.sort(key=lambda r: len(r[0]), reverse=True)

    t0 = time.time()
    files = select_files(ws, a.jobs)
    total = sum(p.stat().st_size for p, _ in files)
    by_top = Counter()
    for p, arc in files:
        parts = arc.split("/")
        key = "runs" if "runs" in parts else ("data" if len(parts) > 2 and parts[2] == "data" else parts[0] if parts[0] != "projects" else "project other")
        by_top[key] += p.stat().st_size
    print(f"{len(files):,} files, {total / 1e9:.2f} GB uncompressed ({time.time() - t0:.0f}s)")
    for k, v in by_top.most_common():
        print(f"  {k:15s} {v / 1e9:6.2f} GB")
    if a.dry_run:
        return

    out_dir = Path(a.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    name = f"lynxgwas-example-data-v{a.version}"
    zpath = out_dir / f"{name}.zip"
    tmp = zpath.with_suffix(".zip.part")
    projects = sorted({arc.split("/")[1] for _, arc in files if arc.startswith("projects/")})
    rewritten, leaks = [], Counter()
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED, compresslevel=9, allowZip64=True) as z:
        for i, (p, arc) in enumerate(files):
            if Path(arc).suffix.lower() in TEXT_EXT or arc.endswith("/config.properties"):
                raw = p.read_bytes()
                try:
                    text = raw.decode("utf-8")
                except UnicodeDecodeError:
                    z.write(p, arc)
                    continue
                new = project_overrides(arc, rewrite(text, rules))
                for s in a.leak:
                    if s.lower() in new.lower():
                        leaks[arc.split("/")[-1] if "runs" not in arc else "runs/*/" + arc.split("/")[-1]] += 1
                if "@@LYNX_" in new:
                    rewritten.append(arc)
                z.writestr(zipfile.ZipInfo.from_file(p, arc), new.encode("utf-8"), zipfile.ZIP_DEFLATED, 9)
            else:
                z.write(p, arc)
            if i % 5000 == 0:
                print(f"  {i:,}/{len(files):,}", flush=True)
        pack = {
            "format": PACK_FORMAT, "name": name, "version": str(a.version), "title": a.title,
            "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "projects": projects, "cross_dataset_runs": a.jobs,
            "files": len(files), "uncompressed_bytes": total,
            "placeholders": PLACEHOLDERS, "rewrite_files": rewritten,
            "gff3": {"file": "gencode.v37lift37.annotation.gff3",
                     "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh37_mapping/gencode.v37lift37.annotation.gff3.gz"},
        }
        z.writestr("pack.json", json.dumps(pack, indent=1))
        here = Path(__file__).resolve().parent / "example_data"
        for extra in ("README.md", "LICENSE.txt"):
            if (here / extra).is_file():
                z.write(here / extra, extra)
    if leaks:
        print("WARNING: leak strings remain in:", dict(leaks.most_common(20)))
    tmp.replace(zpath)
    md5, sha = hashlib.md5(), hashlib.sha256()
    with open(zpath, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            md5.update(chunk); sha.update(chunk)
    info = {"file": zpath.name, "size": zpath.stat().st_size, "md5": md5.hexdigest(), "sha256": sha.hexdigest(),
            "uncompressed_bytes": total, "files": len(files), "projects": len(projects), "version": str(a.version)}
    (out_dir / f"{name}.checksums.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
