# LYNXgwas — Architecture & Reference

Full technical reference for LYNXgwas: Java backend packages, frontend pages, pipeline stages, API
endpoints, per-project layout, and configuration keys. For installation and a quick start, see the
main [README](../README.md).

---

## Table of Contents

1. [Architecture Overview](#architecture-overview)
2. [Directory Structure](#directory-structure)
3. [Java Backend — Source Packages](#java-backend--source-packages)
4. [Frontend — Pages and UI Forms](#frontend--pages-and-ui-forms)
5. [Pipeline Stages and Data Flow](#pipeline-stages-and-data-flow)
6. [Analysis Tools (YAML Descriptors)](#analysis-tools-yaml-descriptors)
7. [API Endpoints Reference](#api-endpoints-reference)
8. [Per-Project Directory Layout](#per-project-directory-layout)
9. [Configuration Reference](#configuration-reference)
10. [Build and Run](#build-and-run)
11. [Key Design Decisions](#key-design-decisions)

---

## Architecture Overview

```
Browser (index.html / viewer.html)
        │  HTTP (port 8765)
        ▼
LocalServer.java  ──────────────────────────────────────────────────────
   │  orchestrates                                                       │
   ├── Main.runPipeline()          4-phase pipeline per project          │
   │     Phase 1: GwasParser, LociParser, GffParser, GenomeSkyline      │
   │     Phase 2: PlinkSubsetter  (PLINK 1.9 subprocess)                │
   │     Phase 3: LdCalculator   (PLINK 1.9 subprocess)                 │
   │     Phase 4: JsonExporter   (per-locus JSON + manifest)            │
   │                                                                      │
   ├── RsidPipeline               rsID recovery (HTSJDK + optional API)  │
   ├── LociIdentifier             GWAS → PLINK clump → loci.txt          │
   ├── LocusUpdater               live resize / split / create loci      │
   └── PluginEngine               run R/binary analysis tools per locus  │
                                                                          │
projects/{id}/data/               ◄──── static JSON served directly ─────┘
```

---

## Directory Structure

```
LYNXgwas/
│
├── src/                          # Java source (all compiled to bin/)
│   ├── Main.java
│   ├── Config.java
│   ├── LocalServer.java
│   ├── GwasParser.java
│   ├── LociParser.java
│   ├── GffParser.java
│   ├── LdCalculator.java
│   ├── PlinkSubsetter.java
│   ├── GenomeSkyline.java
│   ├── SnpAnnotator.java
│   ├── JsonExporter.java
│   ├── LocusUpdater.java
│   ├── LociMutationService.java
│   ├── ProjectMetadata.java
│   ├── ProgressTracker.java
│   ├── Locus.java
│   ├── LocusOutput.java
│   ├── Snp.java
│   ├── Gene.java / Transcript.java / Exon.java
│   │
│   ├── rsid/                     # rsID recovery module
│   │   ├── RsidPipeline.java
│   │   ├── RsidRecovery.java
│   │   ├── RsidMatcher.java
│   │   ├── RsidDetector.java
│   │   ├── RsidApiCompleter.java
│   │   ├── NcbiDbSnpProvider.java
│   │   ├── GnomadProvider.java
│   │   ├── RsidApiCache.java
│   │   ├── RateLimiter.java
│   │   ├── GlobalConfig.java
│   │   ├── RsidProgress.java
│   │   ├── RsidApiProvider.java
│   │   ├── DbSnpRecord.java
│   │   └── MatchResult.java
│   │
│   ├── loci/                     # Loci identification module
│   │   ├── LociIdentifier.java
│   │   └── LociProgress.java
│   │
│   ├── analysis/                 # Fine-mapping / analysis tools module
│   │   ├── PluginEngine.java
│   │   ├── BaseStepPipeline.java
│   │   ├── ToolDescriptor.java
│   │   ├── InputContractWriter.java
│   │   ├── OutputContractValidator.java
│   │   ├── LocusGwasExtractor.java
│   │   ├── SnpMatcher.java
│   │   ├── AlleleHarmonizer.java
│   │   ├── LdMatrixComputer.java
│   │   ├── LdGwasDiagnostic.java
│   │   ├── StepManifest.java
│   │   ├── ContentHasher.java
│   │   ├── StableSnpId.java
│   │   ├── AnalysisColumnProvider.java
│   │   ├── SusieAdapter.java
│   │   ├── FinemapAdapter.java
│   │   ├── CojoAdapter.java
│   │   ├── ColocAdapter.java
│   │   ├── GwamaAdapter.java
│   │   └── SusiexAdapter.java     # cross-project — not a PluginEngine tool
│   │
│   └── export/                   # Excel export module
│       ├── ExcelExporter.java
│       ├── XlsxWriter.java
│       ├── ExportRegistry.java
│       ├── ColumnSpec.java
│       ├── SnpContext.java
│       ├── LocusContext.java
│       ├── SnpColumnProvider.java
│       └── LocusColumnProvider.java
│
├── index.html                    # Home page + project wizard (single file)
├── viewer.html                   # Per-locus interactive viewer (single file)
│
├── tools/                        # Analysis tool descriptors (YAML)
│   ├── susie_finemapping.yaml
│   ├── finemap.yaml
│   ├── cojo_conditional.yaml
│   ├── coloc.yaml
│   └── gwama_meta.yaml
│
├── projects/                     # One subdirectory per project (runtime)
│   └── {project_id}/
│       ├── config.properties
│       ├── project.json
│       ├── annotations.yaml
│       ├── loci.txt
│       ├── loci_detail.tsv
│       ├── data/                 # Generated JSON files
│       ├── plink_subsets/        # PLINK .bed/.bim/.fam per locus
│       └── ld_results/           # Temp LD files (deleted after reading)
│
├── config/
│   └── global.json               # Shared ref panels + SNP databases
│
├── resources/
│   └── gencode.v37.annotation.gff3
│
├── input/                        # GWAS summary statistics files
├── lib/                          # HTSJDK + dependencies (JAR files)
├── bin/                          # Compiled Java classes (output of build.bat)
│
├── build.bat                     # Compile all Java sources → bin/
├── run.bat                       # Start server on port 8765
├── checks.json                   # Validation rules / quality checks
└── assets/
    └── lynxgwas-icon.svg
```

---

## Java Backend — Source Packages

### Core (`src/*.java`)

| File | Role |
|---|---|
| `Main.java` | Entry point. Parses CLI args (`--all`, `--project=`, `--server`). Runs `runPipeline()` for one project or `runMultiProject()` for all. Also contains `PipelineResult` data class. |
| `Config.java` | Loads `config.properties`, resolves all column names, LD settings, dataset metadata. Auto-enables LD if `ref.panel.path` is set. `applyJson()` used by wizard POST. `writeProperties()` is the single canonical config writer. |
| `LocalServer.java` | Built-in HTTP server (Java `com.sun.net.httpserver`, port 8765). Registers all API endpoints. Per-project state held in `ConcurrentHashMap<String, ProjectState>`. Background threads per project for pipeline, rsID, loci, analysis. |
| `GwasParser.java` | Streams GWAS TSV once. Builds per-chromosome interval lookup (sorted by `paddedStart`) so SNPs are correctly assigned to loci even if the GWAS file is not sorted by position. |
| `LociParser.java` | Reads `loci.txt` (columns: `meta_chr`, `meta_start`, `meta_end`). Returns `List<Locus>`. |
| `GffParser.java` | Parses GENCODE GFF3. Builds an in-memory gene/transcript/exon tree. `overlapping()` returns genes within a padded window. `nearestGeneNames()` finds closest genes outside the window. |
| `LdCalculator.java` | For each locus: (A) selects a 100-SNP window around the top SNP from the PLINK subset BIM, (B) runs `plink --r2 --ld-snp` to get r² of all SNPs vs. the top SNP → `r2ByPos` map, (C) runs `plink --r2 square` on the window → pairwise matrix filtered to GWAS SNP positions → `LdTriangle`. All temp files are deleted after reading. Runs in parallel via `ExecutorService`. |
| `PlinkSubsetter.java` | Runs `plink --bfile --chr --from-bp --to-bp --make-bed` to cut a region from the reference panel per locus. Also runs `plink --clump` for loci identification. `findPlink()` searches PATH and common install locations. |
| `GenomeSkyline.java` | Builds a genome-wide binned Manhattan signal array. Written to `data/skyline.json`. Powers the mini-overview strip in the viewer. |
| `SnpAnnotator.java` | Queries NCBI dbSNP E-utilities to fetch rsIDs for lead SNPs that lack them. Used when no rsID column is mapped and no top SNP file is provided. |
| `GwasQc.java` | Genome-wide inflation/confounding triage: streams the full GWAS file once, converts each p-value to a chi-square statistic (via `StatsUtil.qnorm`), and reports the genomic inflation factor λ<sub>GC</sub> (flagged above 1.1). Accepts an optional user-supplied LDSC intercept to derive the attenuation ratio `(intercept-1)/(mean_chi2-1)`, which separates genuine polygenicity from confounding. Not LD-pruned — a quick triage signal, not a substitute for a real LDSC run. Cached to `data/qc.json` by GWAS-file content hash. |
| `StatsUtil.java` | Small self-contained statistics helpers (no external math dependency): `qnorm()` (inverse normal CDF, Acklam's approximation) and `zToP()`/`erfc()` (accurate two-sided p-value from a z-score, avoiding the 1-CDF cancellation that floors extreme p-values). `deriveBetaSeFromOrP()` back-derives beta/SE from OR+p-value for GWAS files that report only the former. |
| `JsonExporter.java` | Serializes `LocusOutput` to JSON (hand-built, no library). Writes `locus_N.json`, `locus_N.js` (JSONP wrapper), and `manifest.json`. `exportManifest()` writes the project-level locus index. |
| `LocusUpdater.java` | Handles live locus operations from the viewer: `update()` (resize boundaries), `create()` (new locus), `split()` (divide into sub-loci), `validateSplit()` (cross-region LD check). Each operation re-streams GWAS, re-runs PLINK, and re-exports JSON for the affected locus. |
| `LociMutationService.java` | Higher-level service wrapping `LocusUpdater`; handles the merge + re-export + manifest rebuild after mutations. |
| `ProjectMetadata.java` | Stores `project.json` (loci count, SNP count, fingerprints, rsID status). `checkStaleness()` compares fingerprints of GWAS, loci, GFF3, config, annotations, and pipeline version. |
| `ProgressTracker.java` | Thread-safe phased progress (phase name, step count, advance). Polled by `/api/project/{id}/progress`. |
| `Locus.java` | `index`, `chr`, `start`, `end`, `paddedStart`, `paddedEnd`, `snps`. `chrInt()` / `chrToInt()` for numeric sorting. |
| `LocusOutput.java` | Full locus output model: `gwasSnps`, `genes`, `ldTriangle`, `topSnp`, `refPanel`, `locusContext` (prev/next navigation), etc. `LdTriangle` is a nested class with `snps[]`, `matrix[][]`, `topSnpRank`. |
| `Snp.java` | `id`, `chr`, `pos`, `pvalue`, `ea`, `nea`, `beta`, `se`, `oddsRatio`, `r2`, `maf`, `infoScore`, `sampleN`. |

### `src/rsid/` — rsID Recovery

| File | Role |
|---|---|
| `RsidPipeline.java` | Orchestrates full rsID recovery: local dbSNP query → optional API completion → patch locus JSONs in-place. |
| `RsidRecovery.java` | Uses HTSJDK to query tabix-indexed dbSNP VCF files per-locus region. Returns `DbSnpRecord` list. |
| `RsidMatcher.java` | Matches GWAS SNPs to dbSNP records by chr:pos + allele verification (forward/reverse strand). |
| `RsidDetector.java` | Auto-detects rsID column in GWAS file header by scanning for `rs\d+` patterns in first few rows. |
| `RsidApiCompleter.java` | For SNPs not matched locally: queries NCBI Variation Services and/or gnomAD GraphQL API. Manages rate limiting and disk cache. |
| `NcbiDbSnpProvider.java` | NCBI E-utilities / Variation Services HTTP client. |
| `GnomadProvider.java` | gnomAD GraphQL HTTP client. |
| `RsidApiCache.java` | Disk-backed cache for API responses (avoid re-querying on repeated runs). |
| `RateLimiter.java` | Token-bucket rate limiter shared across API providers. |
| `GlobalConfig.java` | Loads/saves `config/global.json`. Manages reference panels and SNP database registrations shared across all projects. |

### `src/loci/` — Loci Identification

| File | Role |
|---|---|
| `LociIdentifier.java` | Full loci pipeline: (1) stream GWAS to find candidate SNPs by p-value threshold, (2) match candidates to reference panel by chr:pos, (3) write a PLINK-compatible summary stats file, (4) run `plink --clump` (r² 0.6 default), (5) parse clumped intervals, (6) merge intervals within `mergeDistanceBp`. Writes `loci.txt` + `loci_detail.tsv`. |
| `LociProgress.java` | Progress data class for loci identification (step names, counts, done flag). |

### `src/analysis/` — Analysis Tools Engine

| File | Role |
|---|---|
| `PluginEngine.java` | Discovers tool YAML descriptors from `tools/`, runs them per locus. Manages job lifecycle, cancellation, and result storage. |
| `ToolDescriptor.java` | Parses a tool YAML into a descriptor: `tool`, `label`, `language`, `command`, `requires`, `params`, `output_mapping`. |
| `BaseStepPipeline.java` | Abstract pipeline for each tool run: (1) extract GWAS for locus, (2) harmonize alleles, (3) match to ref panel, (4) compute LD matrix, (5) run LD-GWAS diagnostic, (6) write input contract, (7) execute tool subprocess, (8) validate output contract, (9) store results. |
| `LocusGwasExtractor.java` | Extracts SNPs for a single locus into `base/locus_gwas.tsv`. Also back-derives missing beta/SE from OR+p-value (`StatsUtil.deriveBetaSeFromOrP`) and recomputes a literal/floored p=0 from beta/se via the chi-square survival function, so GWAS files that report only OR, or whose source tool's own p-value underflowed, still produce usable SNPs. |
| `AlleleHarmonizer.java` | Aligns GWAS alleles to reference panel alleles (handles strand flips). Palindromic SNPs (A/T, C/G) can't be strand-resolved from allele letters alone, so they're resolved by comparing the ref panel's A1 frequency against the GWAS's own effect-allele frequency, and dropped (not guessed) when either MAF is too close to 0.5 to be informative or the frequencies don't clearly support one orientation. |
| `SnpMatcher.java` | Matches GWAS SNPs to reference panel variants for the analysis input. |
| `LdMatrixComputer.java` | Computes the full LD r/r² matrix for analysis tools (dense square matrix, written to `ld_r.matrix`). |
| `LdGwasDiagnostic.java` | DENTIST-style LD–GWAS consistency check. Flags SNPs whose z-score is inconsistent with LD neighbours. Writes `consistency_report.tsv`. |
| `InputContractWriter.java` | Writes all input files (`harmonized_gwas.tsv`, `matched_ref.txt`, `ld_r.matrix`, `ld_snp_order.txt`) expected by tool scripts. |
| `OutputContractValidator.java` | Validates tool output files against the `output_mapping` spec. |
| `SusieAdapter.java` | Adapts SuSiE output columns (`chr`, `pos`, `susie_pip`, `susie_cs`, `susie_cs_coverage`) into the annotation system. Uses `susieR::susie_rss` on a bigsnpr-computed, Ledoit-Wolf-shrunk LD matrix. Requires sample size N to be set. |
| `FinemapAdapter.java` | Single-causal-variant Wakefield ABF fine-mapping (chr, pos, `finemap_pip`, `finemap_log10bf`) — a lightweight approximation, not the real multi-causal FINEMAP shotgun-stochastic-search binary. Falls back to ref-panel allele frequency when the GWAS row's own MAF is missing or out of range. Requires sample size N to be set. |
| `CojoAdapter.java` | Prepares GCTA `--cojo-slct` input with allele/freq orientation verified against the ref panel, runs conditional & joint selection plus a `--cojo-cond` pass for conditional p-values on every SNP, and auto-retries at a looser collinearity threshold when the joint estimates look like artifacts (effect inflation, sign flips). Adapts output columns (`cojo_pJ`, `cojo_bJ`, `cojo_bJ_se`, `cojo_pC`, `cojo_bC`, `cojo_bC_se`, `cojo_selected`, plus `chr`/`pos`) into the annotation system, and writes a reliability verdict + artifact flags to `cojo_diag.json`. |
| `ColocAdapter.java` | Joins this locus's own harmonized GWAS (trait 1) against an externally supplied second trait's summary stats file (trait 2, auto-detects chr/pos/ea/nea/beta/se/maf column names), harmonizes alleles, and generates an R script calling `coloc::coloc.abf`. Writes a raw `SNP.PP.H4`/`SNP.PP.H3`/`SNP.PP.H0` output that the generic `PluginEngine.mapOutput()` maps per `coloc.yaml`, plus a `coloc_summary.json` with the standard PP.H4 interpretation scale (&ge;0.75 strong / &ge;0.50 moderate / 0.10&ndash;0.50 weak / &lt;0.10 none). Optional `restrict_to_finemapped` param narrows the tested SNP set to a prior SuSiE/FINEMAP/COJO run's high-confidence SNPs at this locus instead of the whole region. |
| `GwamaAdapter.java` | Writes a single-cohort GWAMA input file (`MARKERNAME`/`EA`/`NEA`/`BETA`/`SE`/`N`/`EAF`/`STRAND`) from `harmonized_gwas.tsv`. Only meta-analyzes one cohort today (a genomic-control pass-through) — real multi-cohort combination needs a cross-project runner like `SusiexAdapter`. |
| `SusiexAdapter.java` | Cross-ancestry fine-mapping via the external [SuSiEx](https://github.com/getian107/SuSiEx) tool. Unlike the tools above this is **not** a per-locus, single-project `tools/*.yaml` plugin — it takes harmonized GWAS + the matched-ref PLINK subset from *multiple projects* at once (driven by `LocalServer#susiexRun`, not `PluginEngine`) and jointly fine-maps a shared locus across ancestries. Output-schema-tolerant parser (SuSiEx's exact output columns vary by version): looks for a PIP-like and an id-like column by name rather than a fixed position. |
| `StepManifest.java` | Records which pipeline steps completed, their timestamps and content hashes (for caching/skip logic). |
| `ContentHasher.java` | SHA-256 hashing of input files to detect changes between runs. |
| `StableSnpId.java` | Generates stable SNP identifiers (chr:pos:a1:a2 canonical form) for cross-step matching. |
| `AnalysisColumnProvider.java` | Registers analysis tool output columns in the export system. |

### `src/export/` — Excel Export

| File | Role |
|---|---|
| `ExcelExporter.java` | Orchestrates full Excel export: iterates loci, collects SNP+locus rows, writes XLSX. |
| `XlsxWriter.java` | Low-level XLSX writer (uses Apache POI-style byte construction without POI dependency). |
| `ExportRegistry.java` | Central registry of all exportable columns (GWAS stats + annotation columns). |
| `ColumnSpec.java` | Column definition: `id`, `label`, `type`, `scope` (per_snp / per_locus). |
| `SnpContext.java` | All per-SNP data bundled for export (SNP + locus metadata + annotations). |
| `LocusContext.java` | All per-locus data bundled for export. |
| `SnpColumnProvider.java` | Provides standard GWAS SNP column values for export. |
| `LocusColumnProvider.java` | Provides locus-level column values for export. |

---

## Frontend — Pages and UI Forms

### `index.html` — Home Page

Single-file frontend (~1,800 lines). No external dependencies except inline D3.js usage for the skyline.

#### Project Table

The main table shows all projects with columns:
- Project name (click to open viewer)
- Loci count, SNP count
- Annotation sources (badge count)
- Status: `up_to_date` | `needs_reprocessing` | `processing` | `never_processed`
- rsID status: `present` | `recovered` | `not_recovered` | `none`
- Loci source: `identified` | `manual` | `none`

Actions per row: **Open**, **Reprocess**, **Edit**, **Delete**, **Get Loci**, **Add rsIDs**, **Analysis**

---

#### Project Wizard (New / Edit Project)

Triggered by **+ New project** or **Edit**. 5 steps:

**Step 1 — Basics**
- Project name (text input; auto-generates project ID as slug)
- Description (optional text input)
- Project ID (text input, locked when editing)

**Step 2 — Files**
- GWAS summary statistics path (text + Browse button → `/pick-file` dialog)
- Loci definition file path (text + Browse; optional — can be generated later)
- GFF3 gene annotation path (text + Browse)
- Reference panel PLINK prefix (text + Browse; `.bed/.bim/.fam` extension stripped automatically)

**Step 3 — Column mapping**
- Auto-detected from GWAS file header via `/api/peek-file-header`
- Required dropdowns (5): Chromosome, Position, P-value, Effect allele (A1), Other allele (A2)
- Optional dropdowns (7): rsID, Variant ID, Beta, Odds Ratio, Standard Error, Sample size N, MAF, Info score

**Step 4 — Parameters**

*Locus parameters:*
- Locus padding (number, default 200,000 bp)
- Max SNPs per locus (number, default 5,000)

*LD computation:*
- Enable LD checkbox (requires ref panel)
- LD triangle boundary — SNPs each side (number, default 100)
- Parallel LD jobs (number, default 4)

*Reference panel:*
- Population label (text, e.g. EUR, EAS)

*Dataset metadata (required for analysis tools):*
- Sample size N (number, required for COJO/SuSiE/FINEMAP)
- Trait type (dropdown: Not specified / Quantitative / Binary)
- N cases (number)
- N controls (number)
- Effect type (dropdown: Auto-detect / Beta / Odds Ratio / Log Odds Ratio)
- Genome build (dropdown: GRCh37/hg19 / GRCh38/hg38)
- Ancestry / Population (text, e.g. EUR, EAS, Hispanic)

**Step 5 — Preview**
- Summary table of all settings before submit
- Submit button: **Process** (new) or **Save & reprocess** (edit)

---

#### Get Loci Modal

Triggered by **Get Loci** button. Parameters:
- Seed p-value threshold (number, default 5e-8)
- Lead p-value threshold (number, default 5e-8)
- Merge distance (kb, default 250)
- Min SNPs per locus (number, default 5)

Progress polling via `/api/loci-progress?project_id=`.

---

#### Add rsIDs Modal

Triggered by **Add rsIDs** button. Parameters:
- SNP database selector (populated from `global.json` SNP databases)
- Enable API completion checkbox (NCBI + gnomAD)

Progress polling via `/api/rsid-progress?project_id=`.

---

#### Resources & Settings Panel

Accessed from the top-right settings icon.

*Reference panels tab:*
- Table of registered panels (path, population, build)
- Add panel form: path (Browse), population label, genome build
- Remove button per panel

*SNP databases tab:*
- Table of registered dbSNP VCF folders
- Add database form: label, folder path (Browse), genome build, file pattern
- Remove button per database

*API settings tab:*
- NCBI API toggle + API key field + rate limit
- gnomAD toggle
- Restrict to lead SNPs toggle

Saved to `config/global.json` via `POST /api/global-config`.

---

#### Analysis Panel

Per-project analysis (triggered from **Analysis** button).

*Tool list view:*
- Lists all discovered tools from `tools/*.yaml`
- Each tool shows: label, description, status, last-run results
- **Run** button opens tool parameter form
- Base-status table: one row per locus with dots for extract/match/harmonize/LD completion, the LD&ndash;GWAS
  diagnostic verdict (`pass` / `warn` / `high_warn`, the last flagging likely reference-panel ancestry
  mismatch), and a **Log** link that opens that locus's `build.log` (narrates each base pipeline step)

*Tool run form:*
- Dynamically generated from tool YAML `params` list
- Field types: `string` (text input), `int`/`float` (number input), `boolean` (checkbox), `select` (dropdown)
- **Run** triggers `POST /api/project/{id}/analysis/run-tool`

*Locus-level analysis view (opened from viewer):*
- Same tool list but scoped to a single locus
- Shows per-locus run log and results
- Results columns (PIP, CS, p-values) are added as annotation tracks

---

#### SuSiEx Modal (cross-ancestry fine-mapping)

Triggered by the **SuSiEx** toolbar button (top-level, not per-project — it spans projects).
- Pick 2+ projects (checkboxes), then a locus dropdown per selected project (populated from that
  project's manifest) — the locus that corresponds to the shared region being fine-mapped.
- Params: max causal signals (L), p-value threshold, MAF threshold, PLINK path, and the local path
  to `SuSiEx.py` (external tool, not bundled).
- Run polls `/api/susiex-progress` then renders `/api/susiex-result` as a SNP × PIP × credible-set
  table, sorted by PIP.

---

### `viewer.html` — Interactive Locus Viewer

Single-file viewer (~4,500 lines). Uses D3.js for all rendering.

#### Main Layout (per locus panel)

From top to bottom:
1. **Panel header**: locus name, chromosome:start-end, top SNP rsID, ref panel population, prev/next navigation
2. **Annotation tracks**: configurable tracks above the Manhattan plot (point, bar, flag styles)
3. **Manhattan plot**: -log10(p) vs position; SNPs colored by LD r² (red→blue gradient); hover tooltip
4. **Context sketch**: mini overview of chromosome with locus position indicator
5. **Gene track**: GENCODE genes (protein-coding = blue, non-coding = green); exon structures; hover gene name
6. **LD triangle**: pairwise r² heatmap (200×200 max); Gabriel-style block outlines; top SNP marker

#### Viewer Controls

- **Loci sidebar**: searchable locus list (by gene name, rsID, chr:pos); click to jump
- **Search bar**: find locus by gene or rsID
- **Track toggles**: Context / Genes / LD / Annotations
- **Zoom**: scroll-wheel on x-axis; all tracks zoom together
- **Locus resize**: drag left/right handles on context sketch → calls `/api/project/{id}/locus/{n}/update`
- **Split locus**: button in panel header → split dialog with LD cross-region validation
- **New locus**: draw region on chromosome sketch → calls `/api/project/{id}/locus/create`
- **Reorder Loci**: renumbers all loci to match genomic (chr:pos) order → `POST /api/project/{id}/reorder-loci`.
  Manually created loci keep their original append-order index otherwise, so this is a manual, reversible action
  rather than something the viewer does automatically — navigation order always follows the server's genomic
  sort regardless.

#### PDF Export Modal (`#pdf-config-modal`)

- Export mode: **per-locus panels** or **whole-genome overview**
- Per-locus mode: resolution (1×–4×), page layout (portrait/landscape), loci-per-page, track selection, locus selection (all/current/range)
- Whole-genome mode: chromosomes-per-page (2/3/4), and a squeeze option that compresses the
  p&lt;5×10⁻¹⁰ region to 60% of plot height (with a visible axis break) so extreme peaks don't
  flatten the rest of the signal
- **Export** button → generates client-side PDF via canvas

#### Annotation Wizard (in viewer)

5-step inline wizard for adding annotation tracks:
1. Select TSV file (Browse)
2. Map join columns (chr, pos, allele columns)
3. Choose columns to display
4. Configure track style (point / bar / flag) and Manhattan channel (color / size / shape)
5. Save → writes `annotations.yaml` via `POST /api/project/{id}/annotation-config`

---

## Pipeline Stages and Data Flow

```
Input files                   Pipeline phase          Output
───────────────────────────────────────────────────────────────────────────
loci.txt (meta_chr/start/end) ─┐
GWAS summary stats (.tsv)     ─┤─ Phase 1: Parse ─► Locus list + SNP lists
GFF3 annotation (.gff3)        │                     Gene tree
                               │                     Skyline data
                               │                     Top SNP per locus
                               └─────────────────────────────────────────

Top SNP + Locus regions        ─┬─ Phase 2: PLINK ─► plink_subsets/locus_N.bed
Reference panel (.bed/.bim/.fam)─┘  subset            plink_subsets/locus_N.bim
                                                       plink_subsets/locus_N.fam

plink_subsets/locus_N.*        ─┬─ Phase 3: LD ───► r2ByPos map per locus
Top SNP BIM IDs                ─┘  computation        LdTriangle (pairwise matrix)
                                                       [temp .ld files deleted]

All of the above               ─── Phase 4: Export ► data/locus_N.json
                                                      data/locus_N.js (JSONP)
                                                      data/manifest.json
                                                      data/skyline.json
```

### GwasParser — How SNPs Are Assigned to Loci

The parser streams the GWAS file once regardless of sort order:
1. Builds a per-chromosome map of loci sorted by `paddedStart`
2. For each GWAS row, binary-searches the sorted list to find the upper bound where `paddedStart > pos`
3. Scans back through loci with `paddedStart <= pos` and collects those where `paddedEnd >= pos`

This correctly assigns SNPs even when the GWAS file is unsorted within chromosomes.

---

## Analysis Tools (YAML Descriptors)

Tool descriptors live in `tools/*.yaml`. New tools are auto-discovered at startup.

### Schema

```yaml
tool: <id>              # unique tool identifier
version: "1.0"
label: "Display Name"
description: "..."
language: R | python | binary
command: "source('{run_dir}/script.R')"   # {run_dir} = locus working dir
requires:               # input files the BaseStepPipeline must provide
  - harmonized_gwas     # harmonized_gwas.tsv
  - matched_ref         # matched_ref.txt + ld_r.matrix + ld_snp_order.txt
  - ld_r                # ld_r.matrix only
params:
  - name: param_name
    type: int | float | string | boolean | select
    default: "value"
    label: "UI label"
    description: "hint text"
    options: [a, b, c]  # for select type only
output_mapping:
  - raw: column_name_in_output    # column in tool's output TSV
    as:  annotation_column_id    # name used in annotation system
    scope: per_snp | per_credible_set
    type: double | int | string
```

### Available Tools

| File | Tool | Language | Outputs |
|---|---|---|---|
| `susie_finemapping.yaml` | SuSiE Fine-Mapping | R | `susie_pip`, `susie_cs`, `susie_cs_coverage` |
| `finemap.yaml` | FINEMAP (ABF) | R | `finemap_pip`, `finemap_log10bf`, `finemap_cs` |
| `cojo_conditional.yaml` | COJO Conditional & Joint (GCTA `--cojo-slct`) | R | `cojo_pJ`, `cojo_bJ`, `cojo_bJ_se`, `cojo_pC`, `cojo_bC`, `cojo_bC_se`, `cojo_selected` |
| `coloc.yaml` | Colocalisation (coloc) | R | `coloc_pp_h4`, `coloc_pp_h3`, `coloc_pp_h0` |
| `gwama_meta.yaml` | GWAMA Meta-Analysis (single-cohort pass-through today) | binary | `gwama_beta`, `gwama_se`, `gwama_p`, `gwama_direction`, `gwama_i2` |

Not a `tools/*.yaml` plugin — see below: **SuSiEx** (cross-ancestry fine-mapping across multiple projects at once), driven by its own `/api/susiex-*` endpoints and the **SuSiEx** toolbar button in `index.html`, since it needs more than one project's artifacts simultaneously and doesn't fit the one-project-one-locus `PluginEngine` model. Requires an external [SuSiEx](https://github.com/getian107/SuSiEx) install (path supplied in the run form) — LYNXgwas prepares each selected project/locus's summary stats in SuSiEx's expected column layout and points `--ld_file` straight at that project's already-computed `matched_ref` PLINK subset.

---

## API Endpoints Reference

All endpoints served by `LocalServer.java` on port 8765.

### Project Management

| Endpoint | Method | Description |
|---|---|---|
| `/api/projects` | GET | List all projects with status, counts, staleness |
| `/api/projects/process` | POST | `{"id":"x"}` or `{"all_stale":true}` — trigger pipeline |
| `/api/project/{id}/progress` | GET | Pipeline progress (phase, step, done flag, completed loci indices) |
| `/api/project/{id}/config` | GET | Read project config as JSON |
| `/api/project/{id}/config` | POST | Write project config (wizard submit) |
| `/api/project/{id}/manifest` | GET | Locus index (names, coords, top SNPs) |
| `/api/project/{id}/locus/{n}` | GET | Full locus JSON (SNPs, genes, LD, annotations) |
| `/api/project/{id}/qc` | GET | Genomic inflation factor (λ<sub>GC</sub>), cached by GWAS-file content hash |
| `/api/project/{id}/qc` | POST | `{"ldsc_intercept":1.02}` — same as GET plus the attenuation ratio against a user-supplied LDSC intercept |
| `/api/project-delete` | POST | `{"id":"x"}` — delete project directory |

### Locus Mutations (from viewer)

| Endpoint | Method | Description |
|---|---|---|
| `/api/project/{id}/locus/{n}/update` | POST | Resize locus boundaries; re-streams GWAS + re-runs PLINK |
| `/api/project/{id}/locus/create` | POST | Create new locus from a chr:start-end region |
| `/api/project/{id}/locus/{n}/split` | POST | Split locus into named sub-regions |
| `/api/project/{id}/locus/{n}/validate-split` | POST | Check cross-region LD before confirming split |
| `/api/project/{id}/reorder-loci` | POST | Renumber all loci to genomic (chr:pos) order |

### rsID Recovery

| Endpoint | Method | Description |
|---|---|---|
| `/api/rsid-recover` | POST | Start rsID recovery (`project_id`, `db_id`, `use_api`) |
| `/api/rsid-progress` | GET | `?project_id=x` — recovery progress |
| `/api/rsid-detect` | GET | `?project_id=x` — auto-detect rsID column in GWAS |

### Loci Identification

| Endpoint | Method | Description |
|---|---|---|
| `/api/loci-identify` | POST | Start loci identification (`project_id`, thresholds, distances) |
| `/api/loci-progress` | GET | `?project_id=x` — identification progress |

### Analysis Tools

| Endpoint | Method | Description |
|---|---|---|
| `/api/project/{id}/analysis/base-status` | GET | Status of base pipeline steps (harmonize, LD, diagnostic) |
| `/api/project/{id}/analysis/run-base` | POST | Run base analysis pipeline for a locus |
| `/api/project/{id}/analysis/tools` | GET | List available tool descriptors |
| `/api/project/{id}/analysis/run-tool` | POST | Run a tool (`tool_id`, `locus_id`, params) |
| `/api/project/{id}/analysis/tool-status` | GET | Status of a running tool job |
| `/api/project/{id}/analysis/tool-result` | GET | Get tool output columns |
| `/api/project/{id}/analysis/cancel-tool` | POST | Cancel a running tool job |
| `/api/project/{id}/analysis/locus-log/{locusId}` | GET | Read a locus's `build.log` from the base pipeline |

### SuSiEx (Cross-Ancestry Fine-Mapping)

| Endpoint | Method | Description |
|---|---|---|
| `/api/susiex-run` | POST | `{"members":[{"project_id","locus_index"}, ...], "params":{...}}` — at least 2 members. Auto-builds each member's base-pipeline artifacts if missing, then runs the external SuSiEx tool. |
| `/api/susiex-progress` | GET | `?job=<id>` — `{"status":"running"\|"done"\|"error", "error":...}` |
| `/api/susiex-result` | GET | `?job=<id>` — `{"ok","n_credible_sets","rows":[{"snp_id","chr","pos","pip","cs_id"}, ...]}` |

### Annotations and Export

| Endpoint | Method | Description |
|---|---|---|
| `/api/project/{id}/annotation-config` | GET | Read `annotations.yaml` |
| `/api/project/{id}/annotation-config` | POST | Write `annotations.yaml` |
| `/api/annotation-file` | GET | `?path=...` — read a TSV annotation file as JSON |
| `/api/export-excel` | POST | Start Excel export; returns job ID |
| `/api/export-progress` | GET | Export progress |

### Utilities

| Endpoint | Method | Description |
|---|---|---|
| `/api/global-config` | GET | Read `config/global.json` |
| `/api/global-config` | POST | Write `config/global.json` |
| `/api/peek-file-header` | GET | `?path=...` — read column names from TSV header |
| `/api/detect-rsid-column` | GET | `?path=...` — scan file for rsID column |
| `/pick-file` | GET | Opens native file picker dialog; returns `{"path":"..."}` |

---

## Per-Project Directory Layout

```
projects/{project_id}/
│
├── config.properties         # All pipeline configuration (see below)
├── project.json              # Metadata: loci count, SNP count, fingerprints, rsID status
├── annotations.yaml          # Annotation track configuration (written by wizard)
├── loci.txt                  # Locus definitions: meta_chr, meta_start, meta_end
├── loci_detail.tsv           # Lead SNPs, p-values, clump sizes per locus
│
├── data/                     # Served directly as static files
│   ├── manifest.json         # [{locusIndex, chr, start, end, topSnp, locusName}, ...]
│   ├── skyline.json          # Genome-wide binned signal
│   ├── locus_1.json          # Full locus data for locus 1
│   ├── locus_1.js            # JSONP wrapper: window.LOCUS_DATA[1] = {...}
│   ├── locus_2.json
│   └── ...
│
├── plink_subsets/            # PLINK binary subsets (one per locus)
│   ├── locus_1.bed
│   ├── locus_1.bim
│   ├── locus_1.fam
│   └── ...
│
└── ld_results/               # Temporary LD computation directory
    │                         # (files deleted after reading; directory persists)
    └── consistency_report.tsv  # DENTIST-style LD diagnostic (if analysis run)
```

### `locus_N.json` Schema

```json
{
  "id": "uuid",
  "locus_index": 1,
  "locus_name": "Locus 1",
  "chr": "1",
  "start": 992819,
  "end": 1191870,
  "padded_start": 792819,
  "padded_end": 1391870,
  "ref_panel": "EAS",
  "top_snp": { "id": "rs123", "chr": "1", "pos": 1050000, "pvalue": 5e-12, ... },
  "gwas_snps": [ { "id": "rs...", "pos": ..., "pvalue": ..., "r2": 0.95, ... }, ... ],
  "genes": [ { "gene_name": "GENE1", "strand": "+", "start": ..., "transcripts": [...] } ],
  "nearest_genes": ["GENE1", "GENE2"],
  "ld_triangle": {
    "snp_count": 201,
    "top_snp_rank": 100,
    "snps": [ { "bim_id": "rs...", "pos": ..., "is_top_snp": false } ],
    "matrix": [[1.0, 0.8, ...], ...]
  },
  "locus_context": {
    "prev_locus": { "index": 0, "chr": "1", "start": ..., "distance_bp": 500000 },
    "next_locus": { "index": 2, ... }
  }
}
```

---

## Configuration Reference

### Per-project: `config.properties`

| Key | Default | Required | Description |
|---|---|---|---|
| `gwas.file` | — | Yes | Path to GWAS summary statistics (tab-delimited TSV) |
| `loci.file` | `input/loci.txt` | No | Path to loci definitions (`meta_chr`, `meta_start`, `meta_end`). Generated by "Get Loci" if absent. |
| `gff3.file` | `resources/gencode.v37.annotation.gff3` | Yes | Path to GENCODE GFF3 annotation |
| `ref.panel.path` | — | No | PLINK bfile prefix (without `.bed`/`.bim`/`.fam`) |
| `ref.panel.population` | `EAS` | No | Population label shown in UI |
| `col.chr` | `chrom` | Yes | GWAS column name for chromosome |
| `col.pos` | `pos` | Yes | GWAS column name for position (integer) |
| `col.pvalue` | `p` | Yes | GWAS column name for p-value |
| `col.ea` | `ea` | Yes | GWAS column name for effect allele |
| `col.nea` | `nea` | Yes | GWAS column name for non-effect allele |
| `col.rsid` | — | No | GWAS column name for rsID (if present in file) |
| `col.varid` | `varid` | No | GWAS column name for variant ID (chr:pos fallback) |
| `col.beta` | — | No | Beta effect size column |
| `col.or` | — | No | Odds ratio column |
| `col.se` | — | No | Standard error column |
| `col.n` | — | No | Sample size per SNP column |
| `col.maf` | — | No | Minor allele frequency column |
| `col.info` | — | No | Imputation info score column |
| `locus.padding` | `200000` | No | Bp added each side of locus for display |
| `ld.enabled` | `false` | No | Enable LD computation (auto-true if `ref.panel.path` set) |
| `ld.triangle.boundary` | `100` | No | SNPs each side of top SNP in LD triangle |
| `ld.r2.threshold` | `0.0` | No | Minimum r² to report |
| `ld.parallel.jobs` | `4` | No | Parallel PLINK LD jobs |
| `threads` | `4` | No | Pipeline threads |
| `max.snps.per.locus` | `5000` | No | Downsample loci exceeding this |
| `split.ld.threshold` | `0.2` | No | r² threshold for cross-region LD warning when splitting |
| `split.min.distance.bp` | `250000` | No | Minimum bp gap between split sub-loci |
| `sample.n` | `0` | **Yes** | Total sample size. Required by COJO, SuSiE, and FINEMAP; the wizard blocks Step 3 without it. |
| `n.cases` / `n.controls` | `0` | **Yes** | Case/control counts, enforced by the wizard alongside `sample.n` (at least one of the pair must be set) |
| `trait.type` | — | No | `quantitative` or `binary` |
| `effect.type` | — | No | `beta`, `OR`, or `logOR` |
| `genome.build` | `GRCh37` | No | `GRCh37` or `GRCh38` |
| `ancestry` | — | No | Population label for analysis tools (e.g. `EAS`, `EUR`) |

### Global: `config/global.json`

```json
{
  "reference_panels": [
    { "id": "...", "path": "...", "population": "EAS", "build": "GRCh37" }
  ],
  "snp_databases": [
    { "id": "dbsnp_hg19", "label": "dbSNP (hg19)", "build": "hg19",
      "folder": "/path/to/vcf/", "file_pattern": "{chr}.vcf.gz" }
  ],
  "ncbi_api": {
    "enabled": false, "api_key": "", "rate_limit_per_sec": 3,
    "restrict_to_lead_snps": true
  }
}
```

---

## Build and Run

### Prerequisites

- Java 11+ (JDK)
- PLINK 1.9 (on PATH or discoverable)
- Reference panel in PLINK binary format (`.bed/.bim/.fam`)
- GENCODE GFF3 annotation file

### Build

```bat
build.bat
```

Compiles all Java sources to `bin/`. Also copies `index.html`, `viewer.html`, and assets to `output/`.

### Run

```bat
run.bat                       # Start server, open http://localhost:8765/
run.bat --all                 # Process all stale projects, then start server
run.bat --project=myproject   # Process one project, then start server
```

The server stays running until Ctrl+C. The browser UI is served at `http://localhost:8765/`.

---

## Key Design Decisions

| Decision | Reason |
|---|---|
| **Single-file frontend** (`index.html`, `viewer.html`) | No build toolchain required; easy to deploy and diff |
| **GWAS parser uses interval lookup, not sorted streaming** | GWAS files are often not sorted by position within chromosomes; the interval lookup assigns SNPs correctly regardless of file order |
| **PLINK temp files deleted after reading** | `ld_results/` contains only transient files; the permanent LD data lives in `data/locus_N.json` |
| **Concurrent per-project threads** | Each project's pipeline, loci identification, rsID recovery, and analysis jobs run in independent background threads; different projects can run simultaneously |
| **Fingerprint-based staleness** | SHA-256 of input files + annotation config + pipeline version; avoids unnecessary reprocessing |
| **Tool descriptors as YAML** | New analysis tools can be added without recompiling Java; the engine auto-discovers `tools/*.yaml` at startup |
| **JSONP wrapper (`locus_N.js`)** | Enables static file serving without CORS issues when opening from a local filesystem |
| **rsID patching in-place** | rsID recovery patches existing locus JSONs without a full pipeline rerun; fingerprint updated to avoid re-patching |
