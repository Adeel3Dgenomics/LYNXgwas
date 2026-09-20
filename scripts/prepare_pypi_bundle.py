"""Regenerate python/lynxgwas/_bundled/ from the current Java build output.

Run this after `build.bat` and before building the PyPI wheel/sdist. It never
edits src/, bin/, or lib/ — it only copies from them.
"""
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BUNDLED = ROOT / "python" / "lynxgwas" / "_bundled"


def reset(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True)


def copy_classes() -> None:
    dest = BUNDLED / "classes"
    reset(dest)
    src = ROOT / "bin"
    for item in src.rglob("*.class"):
        # bin/ also holds this project's standalone tests/*Test.class (compiled there by
        # run_tests.bat, since there's no separate build-tool output directory) — those have no
        # business in the shipped runtime bundle, so skip anything named "*Test.class" or nested
        # under a class named "*Test$..." (its inner classes).
        if item.stem.split("$")[0].endswith("Test"):
            continue
        rel = item.relative_to(src)
        target = dest / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(item, target)


def copy_lib() -> None:
    dest = BUNDLED / "lib"
    reset(dest)
    for jar in (ROOT / "lib").glob("*.jar"):
        shutil.copy2(jar, dest / jar.name)


def copy_web() -> None:
    dest = BUNDLED / "web"
    reset(dest)
    shutil.copy2(ROOT / "index.html", dest / "index.html")
    shutil.copy2(ROOT / "viewer.html", dest / "viewer.html")
    shutil.copy2(ROOT / "gene_constellation.html", dest / "gene_constellation.html")
    # README-only images (screenshots) don't belong in the runtime app bundle
    shutil.copytree(ROOT / "assets", dest / "assets", ignore=shutil.ignore_patterns("screenshots"))
    # index.html/viewer.html reference the logo/favicon at this exact relative path
    shutil.copytree(ROOT / "docs" / "images", dest / "docs" / "images")


def copy_tools() -> None:
    dest = BUNDLED / "tools"
    reset(dest)
    for yaml_file in (ROOT / "tools").glob("*.yaml"):
        shutil.copy2(yaml_file, dest / yaml_file.name)


def main() -> None:
    if not (ROOT / "bin").exists():
        raise SystemExit("bin/ not found — run build.bat first")
    copy_classes()
    copy_lib()
    copy_web()
    copy_tools()
    n_classes = len(list((BUNDLED / "classes").rglob("*.class")))
    n_jars = len(list((BUNDLED / "lib").glob("*.jar")))
    print(f"Bundled {n_classes} class files, {n_jars} jar(s), web assets, and tool descriptors "
          f"into {BUNDLED}")


if __name__ == "__main__":
    main()
