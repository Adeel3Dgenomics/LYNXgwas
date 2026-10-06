# Example data

A fresh LYNXgwas install has no projects. The home page offers the article's example data — the 30
GWAS datasets of the validation corpus (6 diseases × 5, GRCh37) with every analysis already run — as
one download from Zenodo. One click on **Download** does the rest; **Not now** hides the banner, and
the toolbar's **Example data** button brings it back.

## What the user gets

The 30 projects, the two cross-dataset runs over all of them (the Locus Serpent Plot and the Gene
Constellation) with their joint fine-mapping (mvSuSiE) runs, and the Roadmap histone peaks. Everything
can be browsed: the viewer, every tool's latest result at every locus, the summary pages, heritability,
enrichment, the cross-dataset views.

The projects are **results only**: the GWAS summary statistics are not included (they belong to the
consortia that published them) and neither are the LD matrices and other intermediates the app rebuilds.
Such a project carries `"results_only": true` in its `project.json`:

- it always counts as up to date (its inputs are missing by design, not changed), and the home page shows
  it as *Example data*, without the Reprocess and QC buttons;
- the server refuses, with an explanation, everything that needs the summary statistics: building the
  analysis inputs, running a tool, heritability runs, the sample-size correction, QC and locus edits
  (`RESULTS_ONLY_BLOCKED` in `LocalServer`);
- setting the project's GWAS file to the original summary statistics (Edit) turns it back into a normal
  project, which can then be reprocessed.

## How the install works

`ExampleDataInstaller` (desktop app; endpoints under `/api/example-data/`):

1. **Find the pack**: the latest version of the Zenodo record in `config/example_data.json`
   (`zenodo_record`, or `DEFAULT_ZENODO_RECORD` in the code), via
   `https://zenodo.org/api/records/<id>/versions/latest`. A `url` with `size` and `md5`/`sha256` can be
   configured instead. **Install from a file** skips this and uses a zip already on disk.
2. **Check disk space** for the download, the unpacked pack and the gene annotation.
3. **Download** into `example-data/<file>.part`, continuing a partial download when the server allows
   it, retrying on connection errors. Only HTTPS redirects are followed.
4. **Verify** the size and the MD5 (Zenodo's) and/or SHA-256; a damaged file is deleted.
5. **Unpack** into `example-data/staging/`, refusing any entry outside it or outside the pack's known
   folders (`projects/`, `output/`, `regulatory_data/`), and stopping if it unpacks to more than it declares.
6. **Fill in paths**: the pack stores machine paths as placeholders (`@@LYNX_WORKSPACE@@`, …);
   `pack.json` lists the files that contain one, and only those are rewritten.
7. **Gene annotation**: the projects use GENCODE v37lift37 (GRCh37). If `resources/` lacks it, it is
   downloaded from GENCODE (69 MB, 1.7 GB unpacked).
8. **Merge**: projects, cross-dataset runs and joint fine-mapping runs are moved into the workspace. Any
   with a name already present is skipped and reported, never replaced. Regulatory peak files are only
   added where missing. The pip launcher's untouched `getting-started` placeholder is removed.
9. The cross-dataset runs are loaded without a restart; the downloaded zip is deleted.

Cancel keeps the download for the next try. The state (installed version, banner dismissed) is in
`example-data/state.json`.

## Building a new pack

`scripts/export_example_data.py` builds the pack from a processed workspace: the latest run of each tool
at each locus without its bulky inputs, the viewer data (without the `.js` duplicates of the `.json`
files), the chosen cross-dataset runs and their mvSuSiE runs, and the regulatory peaks. Machine paths
become placeholders (`--map`), paths that only make sense on the building machine are cleared
(`--blank`), and `--leak` warns about any file that still names the machine.

Version 1 was built with:

```
T=C:/Users/<user>/AppData/Local/Temp
python scripts/export_example_data.py \
  --workspace $T/claude/scz-full-demo/workspace \
  --jobs d5924980-27fb-497a-8443-1ebc69d076e1 ef832895-6854-454d-9265-bcd5f7acf9f5 \
  --map "$T/claude/scz-full-demo/workspace=@@LYNX_WORKSPACE@@" \
  --map "$T/claude/scz-full-demo/ref1000g=@@LYNX_REF_PANEL@@" \
  --map "$T/claude/scz-full-demo=@@LYNX_GWAS@@" \
  --map "$T/claude/disease-demos=@@LYNX_GWAS@@" \
  --map "C:/Users/<user>/Documents/mygithub/LYNXgwas/resources=@@LYNX_RESOURCES@@" \
  --map "C:/Users/<user>/Documents/mygithub/LYNXgwas=@@LYNX_APP@@" \
  --map "C:/Users/<USER~1>/AppData/Local/Temp/claude/scz-full-demo/workspace=@@LYNX_WORKSPACE@@" \
  --blank "C:/Users/<user>/R/mvsusie-lib" \
  --leak <user> --leak "<USER~1>" \
  --out ../LYNXgwas-example-data --version 1
```

It gave `lynxgwas-example-data-v1.zip`: 801 MB, 3.25 GB unpacked, 111,923 files, 30 projects.
Publish a new pack as a new version of the same Zenodo record; the app always takes the latest.

To test without Zenodo, serve the zip locally and point a workspace's `config/example_data.json` at it
(`url`, `file`, `size`, `md5`); a second instance can run next to the first with `-Dlynxgwas.port=8790`.
