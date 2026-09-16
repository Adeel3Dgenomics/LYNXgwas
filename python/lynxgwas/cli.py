"""Launcher for LYNXgwas: resolves Java/PLINK/GFF3, then runs the bundled Java app."""
from __future__ import annotations

import argparse
import gzip
import importlib.resources
import json
import platform
import shutil
import subprocess
import sys
import time
import urllib.request
import webbrowser
from pathlib import Path

SERVER_URL = "http://localhost:8765/"
PLINK_DOWNLOAD_PAGE = "https://www.cog-genomics.org/plink/1.9/"
GENCODE_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/"
    "gencode.v44.basic.annotation.gff3.gz"
)


def bundled_path(*parts: str) -> Path:
    root = importlib.resources.files("lynxgwas") / "_bundled"
    for part in parts:
        root = root / part
    return Path(str(root))


def config_path(workspace: Path) -> Path:
    return workspace / ".lynxgwas-config.json"


def load_config(workspace: Path) -> dict:
    p = config_path(workspace)
    if p.exists():
        try:
            return json.loads(p.read_text())
        except (json.JSONDecodeError, OSError):
            return {}
    return {}


def save_config(workspace: Path, cfg: dict) -> None:
    config_path(workspace).write_text(json.dumps(cfg, indent=2))


def ensure_workspace_assets(workspace: Path) -> None:
    workspace.mkdir(parents=True, exist_ok=True)
    web = bundled_path("web")
    for name in ("index.html", "viewer.html"):
        dest = workspace / name
        if not dest.exists():
            shutil.copy2(web / name, dest)
    dest_assets = workspace / "assets"
    if not dest_assets.exists():
        shutil.copytree(web / "assets", dest_assets)
    dest_images = workspace / "docs" / "images"
    if not dest_images.exists():
        dest_images.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(web / "docs" / "images", dest_images)
    dest_tools = workspace / "tools"
    if not dest_tools.exists():
        shutil.copytree(bundled_path("tools"), dest_tools)


def require_java() -> None:
    if shutil.which("java") is None:
        print(
            "lynxgwas needs a Java 11+ runtime on PATH.\n"
            "Install one from https://adoptium.net/ and re-run.",
            file=sys.stderr,
        )
        sys.exit(1)


def stage_plink(workspace: Path, resolved: Path) -> None:
    bin_dir = workspace / "bin"
    bin_dir.mkdir(exist_ok=True)
    dest_name = "plink.exe" if platform.system() == "Windows" else "plink"
    dest = bin_dir / dest_name
    if dest.resolve() != resolved.resolve():
        shutil.copy2(resolved, dest)
        if platform.system() != "Windows":
            dest.chmod(0o755)


def resolve_plink(workspace: Path, cli_arg: str | None) -> Path | None:
    cfg = load_config(workspace)
    if cli_arg:
        p = Path(cli_arg).expanduser().resolve()
        cfg["plink"] = str(p)
        save_config(workspace, cfg)
        stage_plink(workspace, p)
        return p
    if cfg.get("plink") and Path(cfg["plink"]).exists():
        p = Path(cfg["plink"])
        stage_plink(workspace, p)
        return p
    found = shutil.which("plink") or shutil.which("plink.exe")
    if found:
        p = Path(found)
        cfg["plink"] = str(p)
        save_config(workspace, cfg)
        stage_plink(workspace, p)
        return p
    print(f"\nPLINK 1.9 (required for LD + fine-mapping) was not found.")
    try:
        answer = input(
            f"Enter a path to your PLINK 1.9 executable, or press Enter to open the "
            f"download page ({PLINK_DOWNLOAD_PAGE}): "
        ).strip()
    except EOFError:
        print(
            "No terminal input available. Pass --plink <path-to-plink> (LD/fine-mapping "
            "features will be unavailable without it).",
            file=sys.stderr,
        )
        return None
    if answer:
        p = Path(answer).expanduser().resolve()
        if p.exists():
            cfg["plink"] = str(p)
            save_config(workspace, cfg)
            stage_plink(workspace, p)
            return p
        print(f"'{p}' does not exist — continuing without PLINK.")
        return None
    webbrowser.open(PLINK_DOWNLOAD_PAGE)
    print("Re-run lynxgwas with --plink <path> once you've downloaded it.")
    return None


def resolve_gff3(workspace: Path, cli_arg: str | None) -> Path | None:
    cfg = load_config(workspace)
    if cli_arg:
        p = Path(cli_arg).expanduser().resolve()
        cfg["gff3"] = str(p)
        save_config(workspace, cfg)
        return p
    if cfg.get("gff3") and Path(cfg["gff3"]).exists():
        return Path(cfg["gff3"])
    resources_dir = workspace / "resources"
    default_dest = resources_dir / "gencode.v44.basic.annotation.gff3"
    if default_dest.exists():
        cfg["gff3"] = str(default_dest)
        save_config(workspace, cfg)
        return default_dest
    print("\nA GENCODE GFF3 gene-annotation file (required for the gene track) was not found.")
    try:
        answer = input(
            "Enter a path to an existing GFF3 file, or press Enter to download "
            f"GENCODE v44 basic annotation now ({GENCODE_URL}): "
        ).strip()
    except EOFError:
        print(
            "No terminal input available. Pass --gff3 <path> once you have one.",
            file=sys.stderr,
        )
        return None
    if answer:
        p = Path(answer).expanduser().resolve()
        if p.exists():
            cfg["gff3"] = str(p)
            save_config(workspace, cfg)
            return p
        print(f"'{p}' does not exist — continuing without a GFF3 annotation.")
        return None
    resources_dir.mkdir(exist_ok=True)
    gz_path = resources_dir / "gencode.v44.basic.annotation.gff3.gz"
    try:
        print("Downloading GENCODE annotation (~45MB)...")
        urllib.request.urlretrieve(GENCODE_URL, gz_path)
        with gzip.open(gz_path, "rb") as src, open(default_dest, "wb") as dst:
            shutil.copyfileobj(src, dst)
        gz_path.unlink()
        cfg["gff3"] = str(default_dest)
        save_config(workspace, cfg)
        print(f"Saved to {default_dest}")
        return default_dest
    except Exception as exc:
        print(f"Download failed ({exc}). Re-run with --gff3 <path> once you have one.")
        return None


def ensure_at_least_one_project(workspace: Path, gff3: Path | None) -> None:
    """The Java app refuses to start its server with zero projects/*/config.properties
    (Main.discoverProjects/main). Seed a harmless placeholder so first-run users reach
    the home page and can use the "+ New project" wizard; they can delete it from there."""
    projects_dir = workspace / "projects"
    if projects_dir.is_dir():
        for child in projects_dir.iterdir():
            if (child / "config.properties").exists():
                return
    placeholder = projects_dir / "getting-started"
    placeholder.mkdir(parents=True, exist_ok=True)
    placeholder_gwas = placeholder / "placeholder_gwas.tsv"
    placeholder_gwas.write_text("chrom\tpos\tp\tea\tnea\n")
    lines = [
        f"gwas.file={placeholder_gwas.as_posix()}",
        f"gff3.file={gff3.as_posix() if gff3 else ''}",
        "sample.n=1",
    ]
    (placeholder / "config.properties").write_text("\n".join(lines) + "\n")


def build_classpath() -> str:
    classes = bundled_path("classes")
    lib = bundled_path("lib")
    jars = sorted(str(j) for j in lib.glob("*.jar"))
    sep = ";" if platform.system() == "Windows" else ":"
    return sep.join([str(classes), *jars])


def main() -> None:
    parser = argparse.ArgumentParser(prog="lynxgwas", description="Run LYNXgwas locally.")
    parser.add_argument("--workspace", default=str(Path.home() / "lynxgwas-workspace"))
    parser.add_argument("--project", dest="project", default=None)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--plink", default=None)
    parser.add_argument("--gff3", default=None)
    parser.add_argument("--no-browser", action="store_true")
    args, extra = parser.parse_known_args()

    require_java()
    workspace = Path(args.workspace).expanduser().resolve()
    ensure_workspace_assets(workspace)
    resolve_plink(workspace, args.plink)
    gff3 = resolve_gff3(workspace, args.gff3)
    ensure_at_least_one_project(workspace, gff3)

    java_args = ["java", "-Xmx8g", "-cp", build_classpath(), "Main"]
    if args.project:
        java_args.append(f"--project={args.project}")
    if args.all:
        java_args.append("--all")
    java_args.extend(extra)

    print(f"Starting LYNXgwas in {workspace} ...")
    proc = subprocess.Popen(java_args, cwd=workspace)
    if not args.no_browser:
        time.sleep(2)
        webbrowser.open(SERVER_URL)
    proc.wait()


if __name__ == "__main__":
    main()
