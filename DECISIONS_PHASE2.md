# Phase 2 decisions log — MAGMA/GCTA/SuSiE improvements, cross-dataset gene visualization,
# evidence enrichment, and manuscript update

Same contract as `DECISIONS.md`: this file is written and updated *incrementally, honestly* while
the user is away, so nothing here is a claim made before it was actually verified. Anything not
yet done is marked `[ ]`; anything done and verified is `[x]` with the evidence that proves it;
anything genuinely blocked (a missing tool, an unsupported platform) says so plainly rather than
being silently skipped.

## 0. Scope of the request (as I understood it)

1. Improve/test the MAGMA, SuSiE, FINEMAP, and other statistical-tool integrations, reusing the
   user's own existing pipelines wherever possible instead of reinventing them.
2. Improve multi-dataset analysis and its visualization; **replace** "Locus Matrix" (user
   explicitly dislikes it) with something better, not just patch it.
3. New circular interactive gene-relationship visualization across diseases: genome-as-circle
   (chromosomes colored), gene nodes sized/colored by significance and direction of effect,
   radial distance from the circle boundary encoding how many datasets contain that gene.
4. ANOVA for significant-gene variation within-disease and between-disease, shown in the same
   fancy visual style.
5. Improve the multi-dataset-analysis setup wizard ("next next" flow) and add a small circular
   genome icon per dataset on the home page.
6. Let users add "evidence" (PPI scores, gene-expression tables, flexible schema/data types) that
   enriches results and renders as a track under the genome plot; add locus-based enrichment
   analysis against that evidence.
7. Export everything as plots + an Excel file carrying real statistics, checked at least three
   times via independent/controlled verification before being trusted.
8. A fancy (not plain-Mermaid) flowchart of the program's structure/tools and the statistical
   analyses performed.
9. Add a GCTA heritability (GREML) analysis.
10. Use all of the above to run a small, real multi-disease demonstration analysis, and update the
    manuscript to describe it.

## 1. Recon findings that shape every decision below

- **Local machine has none of MAGMA, GCTA, FINEMAP, or R/Rscript installed.** Confirmed via
  `which`/`Rscript -e`. This is the single biggest constraint on "testing" claims below.
- **FINEMAP has no official Windows build** (Univ. of Helsinki distributes Linux/macOS binaries
  only) — this was already an implicit gap; now explicit.
- **The app already has a mature 6-adapter analysis suite**, not a blank slate:
  `src/analysis/{SusieAdapter,FinemapAdapter,SusiexAdapter,CojoAdapter(GCTA-COJO),ColocAdapter,
  GwamaAdapter}.java`, each following the same pattern (prepare harmonized inputs → generate a
  script → `PluginEngine` executes it → parse output). **Zero of these six have any automated
  test coverage today** — confirmed via `tests/` containing only `GlobalConfig`, `GwasQc`,
  `LdCalculator`, `MultiLocusScanner`, `SharedStorageResolver`, `SnakemakeSubmitter` tests. This
  matches what the manuscript's Limitations section already says.
- **No MAGMA adapter exists at all** — this is genuinely new work, not an improvement of existing
  code.
- **No GCTA-GREML (heritability) adapter exists** — `CojoAdapter` uses GCTA only for `--cojo-slct`
  conditional selection, not `--make-grm`/`--reml`. Also new work.
- **Excel export already exists** (`src/export/{ExcelExporter,XlsxWriter,ExportRegistry,
  LocusColumnProvider,SnpColumnProvider}.java`, `src/MultiLocusExcelWriter.java`) — extend it,
  don't rebuild it.
- **The wizard UI pattern already exists** (`wizard-overlay`/`wizard-box`/`wiz-steps-bar` in
  `index.html`, used for project creation). The multi-dataset flow ("Locus Matrix" button) is a
  single-screen checkbox picker, not a real multi-step wizard — that's the concrete gap to fix.
- **`braAInny` (the user's other repo) is not a real visualization** — it's a nearly-empty scratch
  folder (one JSON file). "Like brainny's brain visualization" is read here as a *style*
  reference (an interconnected, glowing, node-link aesthetic), not a literal asset to copy.
- **Real, reusable prior art found and used as design/logic reference** (never copy-pasted
  verbatim — read for approach, reimplemented to fit LYNXgwas's own data model):
  - `vegas2_2_MAGMA` (network share) — real MAGMA gene-based + eQTL config layout → informs the
    new `MagmaAdapter`'s command construction.
  - `Correct_Heritability` (network share) — real GCTA GREML Snakemake rules → informs the new
    heritability adapter's `--make-grm`/`--reml` invocation and `.hsq` parsing.
  - `4FinemapSnake`, `susiexSNAKmake` (network share / GitHub) — confirms the existing
    Susie/Finemap/Susiex adapters already follow sound real-world conventions; used as a
    cross-check for the new unit tests' expected values.
  - `manhattantwin`, `plotLD`, `locusZoom` (GitHub, R) — Manhattan/LD visual conventions cross-
    checked against the new circular visualization's design (e.g. chromosome color palette,
    one-representative-per-locus labeling technique already captured in brainny).
  - `Tasks/flowchart/flowchart.html` (local) — the established "fancy," non-Mermaid HTML/CSS/SVG
    infographic style already validated once this session; reused as the visual basis for the new
    program-structure flowchart and, where it fits, the circular gene visualization's chrome.
  - brainny precautions actively applied while building this phase (not just recon reading):
    exclude MHC from fine-mapping/coloc; remove palindromic SNPs before cross-source harmonization;
    verify GWAS genome build empirically; `estimate_residual_variance=FALSE` for susieR with a
    reference-panel (not in-sample) LD matrix; polar-plot label rotation is screen-space in
    matplotlib — irrelevant here since the new viz is D3/SVG, but the *general* lesson (verify
    rotation/orientation visually, don't assume) is applied to the new circular plot's QA.

## 2. Hard boundaries (carried over from `DECISIONS.md`, unchanged)

- No push to any GitHub repo without an explicit URL from the user.
- No autonomous SSH/HPC connection to a real external system.
- No writes to real `projects/*` research data — all new-feature validation happens in the
  existing isolated `scz-full-demo`-style scratch workspace, on the already-processed 30-dataset
  corpus, exactly as before.
- Binary downloads (MAGMA, GCTA) only from the tools' own official distribution pages, mirroring
  the already-approved PLINK auto-download precedent — never from a third-party mirror.
- Every fork's "done" report is independently verified (file existence, test execution, byte-level
  diff) before being written down here as fact — the dominant lesson of Phase 1.
- Statistical correctness claims (ANOVA, enrichment, MAGMA/GCTA output parsing) are checked at
  least three independent ways before being trusted, per the user's explicit instruction:
  (1) a hand-computed expected value from a textbook/independent-tool example baked into the unit
  test itself, (2) confirming the test fails against the pre-fix/absent code
  (via `git stash`), (3) confirming it passes against the real implementation. This is the same
  discipline already used for `MultiLocusScannerTest` etc. in Phase 1, now applied to every new
  statistical module.

## 3. Work items and status

### 3.1 MAGMA adapter — gene-based + gene-set analysis [ ]
New `src/analysis/MagmaAdapter.java` mirroring the existing adapter pattern: build MAGMA's
`.snploc`/`.pval` inputs from harmonized GWAS + gene annotation already parsed by `GffParser`,
generate the `magma --bfile ... --gene-annot ... --pval ...` and `--gene-results --set-annot`
commands, parse `.genes.out`/`.gsa.out`. Attempt official Windows `magma.exe` auto-download
(v1.10, `https://vu.data.surfsara.nl/...` official CTG Lab URL) into `bin/`, same resolution
order as PLINK (`--magma` flag → remembered config → `PATH` → offer download). If download or
execution fails in this environment, the adapter code and its input-generation unit tests still
stand on their own and are reported as such — not silently dropped.

### 3.2 GCTA-GREML heritability adapter [ ]
New `src/analysis/GctaGremlAdapter.java`, reusing `CojoAdapter`'s existing GCTA-binary-resolution
code. Builds `--make-grm` then `--reml --pheno` (or summary-stats-based `--reml` where genotypes
aren't available, falling back to a documented LDSC-style approximation only if true GREML isn't
feasible for a given dataset) and parses the `.hsq` file for h² and its SE. Same auto-download
approach as MAGMA (official GCTA 1.94.x Windows build).

### 3.3 Test coverage for all fine-mapping/meta adapters [ ]
`tests/MagmaAdapterTest.java`, `GctaGremlAdapterTest.java`, `SusieAdapterTest.java`,
`FinemapAdapterTest.java`, `CojoAdapterTest.java`, `ColocAdapterTest.java`, `GwamaAdapterTest.java`
— input-generation and output-parsing logic tested against synthetic data with hand-computed
expected results, run through the fail→pass `git stash` cycle. Where the actual external binary
isn't runnable here (FINEMAP: no Windows build; SuSiE: no R), the test covers everything up to and
including the generated script's exact content (byte-comparable against a hand-written expected
script), and that boundary is stated explicitly rather than implied to be full end-to-end coverage.

### 3.4 ANOVA module [ ]
New `src/analysis/AnovaUtil.java`, hand-rolled (matches `StatsUtil`'s existing no-dependency
convention): one-way ANOVA (F-statistic, between/within sum of squares, p-value via the existing
`StatsUtil` F-distribution/erf-based machinery) for (a) within-disease variation of a gene's
significance across that disease's 5 datasets, (b) between-disease variation of a gene's
significance across the 6 disease groups. Verified against a textbook ANOVA example (hand-computed
F and p) inside the test itself.

### 3.5 Circular cross-dataset gene visualization — replaces Locus Matrix [ ]
New D3.js panel (fancy SVG, not Mermaid, following the `Tasks/flowchart` visual language):
circular genome ideogram (chromosomes as colored arcs), gene nodes positioned at genomic angle,
radius-from-boundary ∝ number of datasets containing the gene, node size ∝ aggregated
significance, node color ∝ effect direction, hover reveals per-dataset breakdown + ANOVA result.
Replaces the "Locus Matrix" button/panel entirely (old code removed, not just hidden).

### 3.6 Multi-dataset wizard rework [ ]
Rebuild the multi-dataset entry point as a real multi-step wizard using the existing
`wizard-overlay`/`wiz-steps-bar` components (dataset selection → gene/significance threshold →
evidence attachment (3.8) → review → run), replacing the current single-screen checkbox picker.

### 3.7 Home-page dataset icon [ ]
Small inline SVG circular genome glyph (mini chromosome-colored ring) next to each project name in
the home page table.

### 3.8 Evidence enrichment (PPI / gene expression / flexible schema) [ ]
New per-project "evidence" upload: CSV/TSV with a gene-symbol column plus arbitrary named
numeric/categorical columns, tagged by type (PPI score, expression, other). Rendered as a track
under the locus/genome plot. New locus-based enrichment test (Fisher's exact / Mann–Whitney,
hand-rolled, tested) comparing genes-in-loci vs. background for the attached evidence.

### 3.9 Excel export extension [ ]
Extend `ExcelExporter`/`XlsxWriter` with new sheets for ANOVA results, enrichment results, MAGMA
gene-based results, and GCTA heritability estimates — reusing the existing `ExportRegistry`
column-provider pattern, tested the same three-fold way as 3.3/3.4.

### 3.10 Program-structure + statistics flowchart [ ]
New fancy HTML/SVG flowchart (Tasks/flowchart style) documenting the full tool chain end-to-end;
exported as a figure for the manuscript.

### 3.11 Small real multi-disease demonstration analysis [ ]
Run the new ANOVA + enrichment (and MAGMA/GCTA if the binaries download successfully in this
environment) against the already-processed, already-verified 30-dataset corpus on disk
(`scz-full-demo` workspace) to get real numbers — not illustrative/mocked ones.

### 3.12 Manuscript update [ ]
New Methods subsections (MAGMA, GCTA-GREML, ANOVA, enrichment, circular visualization), a Results
subsection reporting 3.11's real numbers, Limitations updated honestly for whatever in 3.1–3.9
couldn't be fully end-to-end verified in this environment (no R, no FINEMAP-on-Windows, and
MAGMA/GCTA only if their downloads fail).

## 4. Execution order

Sequential by dependency, verified before moving on: 3.4 (ANOVA, self-contained) →
3.1/3.2 (MAGMA/GCTA adapters + binary acquisition) → 3.3 (tests for all adapters) → 3.5 (circular
viz, needs ANOVA + real gene/dataset data) → 3.11 (demonstration run) → 3.6/3.7/3.8/3.9/3.10 →
3.12 (manuscript, last, so it reports only what's actually been verified by then).
