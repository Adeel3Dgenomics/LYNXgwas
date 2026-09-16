<p align="center">
  <img src="https://raw.githubusercontent.com/AlsammanAlsamman/LYNXgwas/main/docs/images/icon.png" alt="LYNXgwas icon" width="120">
</p>

<h1 align="center">LYNXgwas</h1>
<p align="center"><b>L</b>ocus ana<b>Y</b>sis and ge<b>N</b>omic e<b>X</b>plorer</p>

<p align="center">
<a href="https://pypi.org/project/lynxgwas/"><img src="https://img.shields.io/pypi/v/lynxgwas.svg" alt="PyPI"></a>
<a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"></a>
</p>

A local GWAS visualization and fine-mapping platform. Point it at your GWAS summary statistics
and it identifies loci, recovers rsIDs, computes LD, runs a full fine-mapping suite, and gives you
an interactive Manhattan/LD/gene-track viewer per locus — all in your browser, entirely on your
own machine. No data ever leaves your computer, and no cloud account or upload step is involved.

## Why LYNXgwas

Going from a raw GWAS summary-statistics file to a set of well-characterized, fine-mapped loci
usually means stitching together several separate tools by hand — a clumping step, a liftover or
rsID lookup, a manual LD calculation, one-off scripts for each fine-mapping method, and finally
some plotting code to actually look at the result. LYNXgwas wraps that whole workflow into one
local application with a UI: define a project once, and loci identification, LD, annotation, and
fine-mapping all run through the same pipeline and land in the same interactive viewer — including
comparing the same locus across multiple GWAS datasets, and cross-ancestry fine-mapping across
projects at once.

## Screenshots

**Home page** — manage multiple GWAS projects, see status and locus/SNP counts at a glance:

![LYNXgwas home page](https://raw.githubusercontent.com/AlsammanAlsamman/LYNXgwas/main/assets/screenshots/homepage.png)

**Locus viewer** — Manhattan plot, gene track, and LD triangle for a single locus:

![LYNXgwas locus viewer](https://raw.githubusercontent.com/AlsammanAlsamman/LYNXgwas/main/assets/screenshots/viewer.png)

## Install

```bash
pip install lynxgwas
```

Requires a Java 11+ runtime on your machine ([Adoptium](https://adoptium.net/) is a good source
if you don't have one).

## Quick start

```bash
lynxgwas
```

This starts the local server and opens `http://localhost:8765/` in your browser. On first run it
will ask for (or offer to download) two things the pipeline needs:

- **PLINK 1.9** — for LD computation and reference-panel subsetting (`--plink <path>` to skip the prompt)
- **A GENCODE GFF3 annotation** — for the gene track (`--gff3 <path>` to skip the prompt)

Both are remembered afterward. From the home page, use **+ New project** to point LYNXgwas at
your own GWAS summary statistics and walk through the setup wizard.

## What it does

**Loci & annotation**
- Automatic loci identification from GWAS summary stats via PLINK clumping (or supply your own `loci.txt`)
- rsID recovery from a local dbSNP VCF, with optional NCBI/gnomAD API completion for anything left unmatched
- Gene-track annotation from a GENCODE GFF3, with nearest-gene lookup outside the plotted window
- Custom annotation tracks from your own TSV files (point/bar/flag styles) via an in-viewer wizard
- Genome-wide QC triage: genomic inflation factor (λ<sub>GC</sub>), with an optional LDSC-intercept comparison to separate genuine polygenicity from confounding

**LD & fine-mapping**
- LD computation and pairwise LD triangles per locus (PLINK reference-panel subsetting under the hood)
- **SuSiE** — Bayesian sum-of-single-effects fine-mapping with credible sets
- **FINEMAP**-style Wakefield ABF fine-mapping for single-causal-variant loci
- **GCTA-COJO** — conditional & joint SNP selection with an automatic reliability/artifact check
- **coloc** — colocalization against a second trait's summary stats (PP.H0–H4)
- **GWAMA** — meta-analysis pass-through
- **SuSiEx** — cross-ancestry joint fine-mapping across multiple LYNXgwas projects at once
- **Locus Matrix** — compare the same locus region across many GWAS datasets side by side

**Viewer & output**
- Interactive Manhattan plot, gene track, and LD triangle, all zooming together
- Live locus editing: resize, split, create, and reorder loci without re-running the whole pipeline
- Per-locus and whole-genome PDF export; full Excel export of all annotated SNPs
- Multi-project management from a single home page, with per-project staleness tracking

## Requirements

| Requirement | Why | How it's handled |
|---|---|---|
| Java 11+ | Runs the backend | Install separately (bundled Java classes only) |
| PLINK 1.9 | LD, clumping, reference panel subsetting | Prompted for a path, or auto-downloaded |
| GENCODE GFF3 | Gene track annotation | Prompted for a path, or auto-downloaded |
| A PLINK-format reference panel (`.bed`/`.bim`/`.fam`) | LD computation | You supply this in the project wizard |

## Documentation

Full architecture, Java package reference, API endpoints, and configuration keys:
[`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md)

## License

[MIT](LICENSE)
