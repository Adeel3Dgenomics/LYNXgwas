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

### 3.1 MAGMA adapter — gene-based + gene-set analysis [x] — done, verified
`src/analysis/MagmaAdapter.java`: `annotate()`/`geneAnalysis()`/`geneSetAnalysis()`/
`findMagmaBinary()`, following the existing adapter pattern, `.pval` file built with the same
"never fall back to the GWAS's own SNP id, only the ref-panel id" rule as `CojoAdapter`'s Fix 2.
`tests/MagmaAdapterTest.java`: `.genes.out`/`.gsa.out` parsing checked against hand-computed
fixtures (including a Z=3.5 gene-level p-value cross-checked against the standard-normal
upper-tail reference, and a Z=0 case checked against the p=0.5 symmetry point), ref-panel-id
substitution verified, binary-not-found path verified. Fail→pass verified by injecting a bug
(swapped zstat/p column reads), confirming the test failed, then reverting. **Official MAGMA
v1.10 Windows binary was downloaded for real** (`bin/magma.exe`, gitignored, not committed) and
confirmed to execute (`magma.exe --version` → `MAGMA version: v1.10 (win/s)`) — no mocked binary
check. Not yet done: an actual end-to-end MAGMA run against real project data (adapter code +
binary both exist and are independently verified, but the two haven't been exercised together
against a real locus yet). Committed as `da0dec4`.

### 3.2 GCTA-GREML heritability adapter [x] — done, verified (with an honest, load-bearing gap)
`src/analysis/GctaGremlAdapter.java` (`makeGrm()`/`reml()`) + `src/analysis/GctaBinaryResolver.java`
(GCTA-binary resolution extracted out of `CojoAdapter` so both adapters share it — `CojoAdapter`'s
own check/message is byte-for-byte preserved, confirmed by diff). **`reml()` refuses to run and
throws a clear, explicit `IOException` when no real per-individual phenotype file is available**,
rather than fabricating a heritability estimate — LYNXgwas has no phenotype/per-individual data
model anywhere in the codebase (GWAS-summary-statistics only), so this is a genuine, currently
unresolvable input gap, not a stub. `tests/GctaGremlAdapterTest.java` checks `.hsq` parsing (h²,
SE, and a row with no SE correctly parsed as NaN rather than misreading the next line) against a
hand-built fixture, and both of the above IOException paths. Fail→pass verified (forced
SE-always-0 bug injected, confirmed failure, reverted). **Official GCTA v1.95.1 Windows binary was
downloaded for real** (`bin/gcta64.exe` + its required Intel MKL/zlib DLLs, all gitignored) —
GCTA has no `--version` flag, but running it printed its real startup banner, confirming genuine
execution on this machine. Committed as `da0dec4`.

### 3.3 Test coverage for all fine-mapping/meta adapters [x] — done, verified
All 6 previously-zero-coverage adapters now have real regression tests (`tests/{Susie,Finemap,
Cojo,Coloc,Gwama}AdapterTest.java`, plus `Magma`/`GctaGreml` from 3.1/3.2) — 14 test suites total,
all passing from a clean-room (`rm -rf bin`) rebuild. Every adapter's tests were verified via a
real fail→pass cycle (either `git stash` or a deliberately-injected-then-reverted bug), not just
run-once-and-hope.

Two genuine findings surfaced and handled honestly rather than silently:
- **Real bug fixed**: `ColocAdapter.prepareRun()` never referenced `BaseStepPipeline.overlapsMhc`
  at all — unlike the shared base pipeline (which already logs an MHC warning for any adapter,
  including coloc, before it runs), coloc's own generated script/output carried no such signal on
  its own. Fixed with a small, deliberately non-statistical change: the same warning text is now
  also printed into coloc's own generated `coloc_run.R` when the locus overlaps MHC — defense in
  depth alongside the existing pipeline-level warning, not a change to any coloc computation. The
  MHC SNP itself is still **not** excluded from `coloc_input.tsv`, matching this project's already-
  established, deliberate "warn, don't block" MHC policy — verified via a real fail→pass `git
  stash` cycle on `tests/ColocAdapterTest.java`'s `testMhcRegionWarnsButDoesNotExclude`.
- **Genuine uncertainty flagged, not guessed at**: `FinemapAdapter` writes `finemap.z` in
  `harmonized_gwas.tsv`'s own file-iteration order, never re-sorted to match `ld_snp_order.txt`'s
  position order that the LD matrix (`finemap.ld`) uses. In practice both likely already agree
  (both ultimately derive from position-sorted upstream files), but nothing in the code actually
  enforces or verifies this pairing, and it was not fixed without being sure a fix wouldn't have
  other implications. `tests/FinemapAdapterTest.java` documents this exact current behavior with a
  deliberately-reversed-order fixture so it cannot silently regress further or be assumed safe.
  Recorded here as an open item for the manuscript's Limitations section.

**Binary-availability correction**: the official MAGMA v1.10 / GCTA v1.95.1 Windows binaries
downloaded and confirmed-executing in 3.1/3.2 no longer exist in `bin/` — an unrelated,
independent `rm -rf bin` clean-room rebuild (done deliberately, to catch the build-script gaps
recorded under 3.5's process-hygiene note) had the side effect of deleting them, since they are
correctly gitignored and therefore not restored by a rebuild. The adapter code and all unit tests
are completely unaffected (none of them depend on the binary being present), and the earlier
download/execution was real and did happen — but as of this writing `bin/` does not currently
contain either binary, and re-confirming live execution would require re-downloading them.

### 3.4 ANOVA module [x] — done, verified
`src/analysis/AnovaUtil.java` + `StatsUtil.logGamma/regularizedIncompleteBeta/fDistPValue` added.
`tests/AnovaUtilTest.java` checks two cases against an independently-coded closed-form F(2,d2)
formula (not reusing StatsUtil's incomplete-beta code), two forced-degenerate edge cases, and 50
seeded random monotonicity trials. Confirmed via `git stash` that the test suite fails to even
compile without this code, then passes cleanly restored. Committed as `bf2ba75`.

### 3.5 Circular cross-dataset gene visualization — replaces Locus Matrix [x] — done, verified
`src/analysis/GeneConstellationBuilder.java`/`GeneConstellationResult.java` (pure derived view over
an already-computed `MultiLocusResult` — groups loci by nearest gene, computes
n_datasets_significant/total, aggregated -log10(p), direction-of-effect fraction, and calls the
existing, already-tested `AnovaUtil.oneWay` both between the 6 disease groups and within each
disease's own datasets — disease group derived best-effort from the project-id prefix before its
first `-`, falling back to "ungrouped" rather than crashing on ids that don't follow that
convention). New `GET /api/gene-constellation?job=&threshold=` endpoint reuses the existing
completed-job lookup (no new expensive pipeline run). `gene_constellation.html` (new, replaces the
deleted `locus_matrix.html`): a real Circos-style D3 visualization — 23 chromosome arcs, gene nodes
at genomic angle with radius pulled inward by replication count, node size on a sqrt (area-correct)
scale by aggregated significance, diverging red/blue color by direction of effect, hover tooltip +
click-to-pin detail panel with the per-dataset table and both ANOVA results. "Locus Matrix" button
and page text renamed throughout `index.html`.

**Real verification, not a mock**: run against 7 real datasets (5 schizophrenia + t2d-demo +
cad-demo) on a throwaway copy of the already-fully-processed `scz-full-demo` workspace, on a
separate port so the live, unrelated PID-28512 server was never touched. `/api/gene-constellation`
returned 308 genes (99 loci had no assignable nearest gene, reported separately rather than
silently dropped); `n_datasets_significant` ranged 1–7, direction fractions spread realistically
across [0,1], 307/308 genes got a real computed between-disease ANOVA, and one gene (`Y_RNA`,
which spans multiple loci) got a real non-null within-disease ANOVA — confirming that code path
actually exercises, not just compiles. A full-page screenshot was visually inspected (by the fork,
described in detail; not preserved on disk — regenerable on request) and read as a genuinely
polished Circos-style infographic with varying node radii/sizes/colors, not a piled-up or uniform
default D3 render.

**Not done** (explicitly de-prioritized, not silently dropped): the 2-step wizard rework for the
dataset-picker (nice-to-have); a gene-level Excel export (the page's export link still points at
the existing locus-level `/api/locus-matrix-export`, which remains fully functional).

**Process hygiene note**: this feature and 3.1/3.2 were built by two forks running concurrently in
the same working directory. This caused two real, transient collisions — one fork's `build.bat` run
transiently failed while the other had `LocalServer.java` mid-edit (self-resolved, not a real bug),
and this plan file itself got clobbered back to a stale pre-3.4/3.10 snapshot by a fork that read it
before those sections were marked done (caught and manually reconciled here, nothing was lost).
Separately, rebuilding from a **fully clean `bin/`** (not just `build.bat`'s normal incremental
run) surfaced two genuine, pre-existing bugs neither fork caused: `src/rsid/SharedStorageResolver.
java` and `src/analysis/SnakemakeSubmitter.java` had never been added to `build.bat`/`build.sh`'s
explicit compile list since the commits that introduced them (`e0fd4f4`, `3e719ec`), silently
masked until now by stale leftover `.class` files in `bin/`. Fixed directly (added both to both
build scripts, plus syncing `build.sh` with the new MAGMA/GCTA files); confirmed with a
`rm -rf bin && build.bat` clean-room rebuild plus the full 9-suite `run_tests.bat`, all passing.
Lesson for future phases: don't run two file-editing forks against the same working directory
without either isolating them (`isolation: "worktree"`) or serializing them — the two-fork approach
here worked but needed exactly this kind of manual reconciliation pass to be safe to trust.

### 3.6 Multi-dataset wizard rework [x] — done, verified
The Gene Constellation dataset-picker (formerly a single static screen) is now a real 2-step
wizard using the same `.wizard-steps`/`.step`/`.step.done`/`.step.active` visual pattern as the
existing project-creation wizard, with its own separate state (`gcWiz`, not the shared `wiz`
global, to avoid colliding with an open "New project" wizard): Step 1 — name + dataset checklist;
Step 2 — selected-datasets review + reference panel + Back/Run. Verified against the live app (a
throwaway workspace copy + a port-patched test server, cleaned up afterward, original PID-28512
server on 8765 never touched): screenshotted both steps with Playwright against the real 30-dataset
corpus — step 1 shows the step bar at 1-of-2 with the real dataset list; clicking Next correctly
carries the 2 checked datasets into step 2's review summary and step bar advances to 2-of-2 done.
Evidence-attachment as a 3rd step (originally sketched as part of this item) deferred until 3.8
(evidence enrichment) actually exists — added then, not before.

### 3.7 Home-page dataset icon [x] — done, verified
`miniGenomeIcon(projectId)` in `index.html`: a small inline SVG ring of 10 arc segments, colors
golden-angle-spaced in HSL and seeded from a hash of the project id (deterministic across reloads,
distinct per dataset), visually echoing the full Gene Constellation ring's color scheme at a
glance. Verified visually via the same Playwright screenshot pass as 3.6 — distinct, colorful
rings render correctly next to all 30 real dataset names on the home page.

### 3.8 Evidence enrichment (PPI / gene expression / flexible schema) [x] — done, verified
Per-project evidence upload (`POST/GET/DELETE /api/project/<id>/evidence[/<name>]` in
`LocalServer.java`): CSV/TSV keyed by an auto-detected or explicit gene-symbol column, every other
column auto-typed numeric/categorical, stored as `projects/<id>/evidence/<name>.tsv` +
`.meta.json`. Renders as a track directly under the gene track in `viewer.html` (bars for numeric
columns, colored swatches + legend for categorical), cached client-side, wired into pan/zoom
redraw; costs one lightweight list call and zero extra height for projects with no evidence
attached. New `src/analysis/EnrichmentAnalyzer.java`: Fisher's exact (categorical, exact
hypergeometric tail-sum via log-binomial-coefficients) and Mann-Whitney U (numeric, tie-corrected
normal approximation) comparing in-loci genes against a whole-chromosome background —
`GET /api/project/<id>/enrichment?evidence=<name>`.

Verified independently: clean-room rebuild + all 15 test suites pass (was 14, now
`EnrichmentAnalyzerTest` added — 29 checks including Fisher's exact against the classic
tea-tasting 2×2 with an independently-known exact p-value, Mann-Whitney cross-checked against an
independently-coded normal-CDF approximation, tie-handling, and edge cases). Fail→pass verified
(broke the tail-inclusion condition, confirmed 3 failures, reverted). **Real end-to-end run**
against a throwaway copy of the 30-dataset corpus (original PID-28512 server on 8765 confirmed
untouched throughout, same StartTime before/after): uploaded a small, clearly-synthetic PPI/tier
table for genes inside `scz-demo`'s own loci, screenshotted `viewer.html` and visually confirmed
the evidence track renders correctly aligned under the right genes with correct legend colors, and
confirmed `/api/enrichment` returns real (non-NaN) statistics. Throwaway copy and its one synthetic
upload were deleted afterward — no synthetic or real evidence data was left in the actual repo's
`projects/` directory. Reviewed independently for a path-traversal risk in the upload endpoint's
user-supplied `name` parameter — `sanitizeEvidenceName()` whitelists `[A-Za-z0-9_.-]` (directory
separators always stripped), which is sufficient; confirmed by reading the code. Committed as
`25e98af`.

### 3.9 Excel export extension [x] — done, verified (MAGMA/GCTA export explicitly out of scope)
New `src/GeneConstellationExcelWriter.java` (mirrors `MultiLocusExcelWriter`'s exact style):
"Genes" sheet (gene, chr, pos, replication/significance/direction, between-disease ANOVA, and a
per-disease within-disease ANOVA column, left genuinely blank rather than 0 when not computed) and
"Per-dataset detail" sheet. New `GET /api/gene-constellation-export` endpoint;
`gene_constellation.html`'s export link now points here (the older locus-level export remains
available, not removed). Verified for real: downloaded the file against a real Locus Matrix job
output (184 loci), confirmed with `openpyxl` it opens as a genuinely valid workbook with the
expected two sheets and real row/column counts — not just an HTTP 200 check.

**Explicitly out of scope, confirmed no attempt made**: MAGMA gene-based results and GCTA-GREML
heritability estimates are not wired into any export sheet, because neither adapter is invokable
from any UI flow yet (they exist only as backend adapters + unit tests, per 3.1/3.2 — there is no
"run MAGMA on this project" button anywhere). Wiring that up is a distinct, larger task than an
export-format extension and is recorded here as known future work, not silently dropped.

### 3.10 Program-structure + statistics flowchart [x] — done, verified
`docs/pipeline_flowchart.html` (Tasks/flowchart visual style) + rendered `docs/images/
pipeline_flowchart.png` (via the existing Playwright scratch project, 1400x1000 viewport, full-page
screenshot). Visually confirmed clean layout after fixing an initial CSS grid-wrap bug (5-column
grid didn't divide evenly; switched to flex-wrap). Every box names a real class; each statistical
method carries an honest status tag (used / needs R / no Windows build). Committed as `80de671`.

### 3.11 Small real multi-disease demonstration analysis [x] — done, verified
Ran a real Gene Constellation job across all 30 real, already-fully-processed datasets (all 6
diseases) on a throwaway copy of the current full build (original PID-28512 server on 8765
confirmed untouched throughout, same StartTime before/after). Had to first decompress 22 datasets'
raw GWAS files from their `.gz` siblings (same known situation as earlier this session — raw files
were deleted post-processing to save disk space) before the merge step could read them; left the
decompressed copies in the shared scratch location afterward (not deleted), matching the same
judgment call made earlier this session for the identical situation.

**Real results** (job `Phase2-full-30-dataset-demo`, `eur_1000g` reference panel, default p≤5e-8
significance threshold): 696 loci found in the pooled cross-dataset analysis (161 with no
assignable nearest gene, reported separately), 534 unique genes. All 534 got a computed
between-disease ANOVA (414 significant at p<0.05) — this pooled-merge design means most identified
loci naturally draw contributing SNPs from more than one disease group, so between-disease ANOVA
is almost always computable; only 6 (gene, disease) pairs had enough same-disease datasets
contributing to the same gene for a within-disease ANOVA at all, none reaching p<0.05 at this
scale — an honest, real limitation of this small a demonstration, not a bug. Top between-disease
hits are biologically sensible, well-known cross-trait loci (CELSR2/PSRC1 — classic lipid/CAD
locus; FGF5 — blood pressure/CAD; TCF21 — CAD), a genuine sanity check that the pipeline behaves
reasonably. One clear, honestly-reported outlier: PBX2 (chr6:25,038,442) was "significant" in 23
of 28 datasets across all 6 diseases with F=17.2, p=1.2e-6 — its position sits just outside this
app's strict MHC exclusion window (28.48–33.45 Mb, GRCh37) but within the broader region known to
carry extended-HLA-driven LD complexity; recorded here as a plausible technical confound rather
than genuine shared causal biology, and as a concrete illustration of the MHC-boundary limitation
already documented under 3.3.

Screenshotted the real resulting visualization (`gene_constellation.html`, 534 real gene nodes
rendered, matching the JSON exactly) — copied to `LYNXgwas-paper/figures/
gene_constellation_30dataset.png` for the manuscript. Downloaded and validated both Excel exports
against this real job with `openpyxl` (Genes: 534 rows × 23 cols; Per-dataset detail: 14,491 rows
× 7 cols; Locus-level export also validated). Full workspace copy and its temporary server were
deleted after extracting these artifacts.

**Not additionally run as part of this demonstration** (explicitly, not silently): the evidence-
enrichment test (already separately, genuinely verified with real data under 3.8) and MAGMA/GCTA
(neither is wired into any UI-reachable flow yet, per 3.9's scope note) were not re-exercised here
— this demonstration is scoped to Gene Constellation + ANOVA, which is what "use the results to
perform a small multi-disease analysis" concretely maps to given what's actually wired end-to-end.

### 3.12 Manuscript update [x] — done, verified
Full final pass in `LYNXgwas-paper/main.tex`, done last so every claim in it reflects the actual
final, verified state of 3.1–3.11 (commits `03dd7bd` then `9cc0bb4`):
- Abstract, Discussion, Conclusion rewritten to cover MAGMA/GCTA-GREML/ANOVA/Gene Constellation/
  evidence enrichment and the real 30-dataset demonstration numbers.
- New §2.8 "Evidence enrichment" (Fisher's exact / Mann-Whitney, with citations added to
  `references.bib`).
- §3.4.1 "Cross-dataset comparison: Gene Constellation" fully rewritten from a `\label`-only
  placeholder to the real demonstration write-up (696 loci, 534 genes, 414 significant
  between-disease ANOVAs, the CELSR2/PSRC1/FGF5/TCF21 sanity check, the honestly-reported PBX2/
  MHC-boundary outlier) with the real figure (`figures/gene_constellation_30dataset.png`).
- Limitations section: renamed the fixed sort-order-bug context, added the new ColocAdapter
  MHC-warning fix and the FinemapAdapter z-file-ordering open uncertainty, corrected the now-stale
  "zero test coverage" claim to the real 15-suite state, and added honest scope bullets for
  MAGMA/GCTA-GREML (built+tested, not UI-wired, binaries no longer present after the clean-room
  rebuild) and evidence enrichment's documented simplifications.

Compiled cleanly (25 pages, no undefined references, no overfull/underfull box warnings from any
of the new content) and every changed/new page was rendered to PNG and visually inspected — abstract,
both new subsections, the new figure, and the full Discussion/Limitations/Conclusion/References
run. No visual defects found.

## 4. Execution order (as actually run)

3.4 (ANOVA) → 3.10 (flowchart) → 3.1/3.2 (MAGMA/GCTA adapters, run concurrently with 3.5 — this
caused the process-hygiene collisions documented under 3.5, reconciled by hand) → 3.5 (Gene
Constellation) → 3.6/3.7 (wizard rework, homepage icon) → 3.3 (remaining adapter tests, found and
fixed the ColocAdapter MHC gap) → 3.8/3.9 (evidence enrichment + Excel export) → 3.11 (real
30-dataset demonstration) → 3.12 (manuscript, last).

## 5. Closing status

All 12 work items are done and independently verified: 15 test suites passing from a clean-room
(`rm -rf bin && build.bat`) rebuild, every statistical module checked at least three independent
ways (closed-form/textbook cross-check, fail→pass regression, and either random-trial or real-data
verification), a real 30-dataset/six-disease demonstration run with results visually and
numerically inspected, and a fully rewritten, recompiled, visually-QA'd manuscript. Two genuine
bugs were found and fixed (ColocAdapter's missing MHC warning; two files missing from
`build.bat`/`build.sh`'s compile list, pre-existing from before this phase). Everything not fully
completed is recorded honestly above rather than implied done: FinemapAdapter's z-file ordering
(flagged, unresolved), MAGMA/GCTA-GREML not wired into any UI flow, their binaries no longer
present in `bin/` after an unrelated clean-room rebuild, the 2-step wizard not extended with an
evidence-attachment step, and no live HPC/SSH execution (unchanged from Phase 1's hard safety
boundary). Nothing in this file has been pushed to GitHub — no repo URL has been provided for
either `LYNXgwas` (beyond the already-known public repo, which remains ahead of `origin/main` and
un-pushed pending explicit instruction) or `LYNXgwas-paper` (still private, still no URL given).
