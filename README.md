<p align="center">
  <img src="docs/images/icon.png" alt="LYNXgwas" width="200">
</p>

<h1 align="center">LYNXgwas</h1>

<p align="center">
  <strong>Locus Analysis and Genomic Explorer</strong><br>
  Interactive GWAS locus visualization with LD computation and gene annotation
</p>

<p align="center">
  <img src="https://img.shields.io/badge/language-Java-orange" alt="Java">
  <img src="https://img.shields.io/badge/platform-Windows%20%7C%20Linux%20%7C%20macOS-blue" alt="Platform">
  <img src="https://img.shields.io/badge/license-MIT-green" alt="License">
</p>

---

## Overview

LYNXgwas is a local-first GWAS locus visualization tool that takes summary statistics, locus definitions, and a GFF3 gene annotation, then produces an interactive browser-based explorer with:

- **Manhattan plots** per locus with gene tracks
- **LD heatmap triangles** computed via PLINK against a reference panel (e.g. 1000 Genomes)
- **Gabriel-style LD block detection**
- **Gene annotation** overlays from GENCODE GFF3
- **Genome-wide skyline** overview across all loci
- **Live progress tracking** via a built-in HTTP server

## Quick Start

### Prerequisites

- **Java 11+** (JDK for building, JRE for running)
- **PLINK 1.9** (optional, for LD computation)
- A reference panel in PLINK binary format (`.bed/.bim/.fam`)

### Build & Run

```bash
# Compile
build.bat

# Configure your paths in config.properties, then:
run.bat
```

The viewer opens automatically at `http://localhost:8080/index.html`.

## Configuration

Edit `config.properties` to set input paths and parameters:

| Key | Description |
|-----|-------------|
| `gwas.file` | GWAS summary statistics (tab-delimited) |
| `loci.file` | Locus definitions (chr, start, end) |
| `gff3.file` | GENCODE GFF3 annotation file |
| `ref.panel.path` | PLINK bfile prefix for LD reference panel |
| `ref.panel.population` | Population label (e.g. `EUR`, `EAS`) |
| `ld.enabled` | Enable/disable LD computation |
| `output.dir` | Output directory |

## Project Structure

```
LYNXgwas/
├── src/                  # Java source files
│   ├── Main.java         # Pipeline orchestrator
│   ├── GwasParser.java   # GWAS summary stats streamer
│   ├── LociParser.java   # Locus definitions parser
│   ├── GffParser.java    # GFF3 gene annotation parser
│   ├── LdCalculator.java # Pairwise LD computation
│   ├── JsonExporter.java # JSON output for the viewer
│   ├── LocalServer.java  # Built-in HTTP server
│   └── ...
├── index.html            # Interactive viewer (single-page app)
├── config.properties     # Runtime configuration
├── build.bat             # Compile script
├── run.bat               # Run script
└── docs/images/          # Project assets
```

## Output

After a run, `output/` contains:

- `index.html` — Interactive viewer with all locus panels
- `data/` — Per-locus JSON files consumed by the viewer
- `plots/` — Generated plot assets

## License

MIT
