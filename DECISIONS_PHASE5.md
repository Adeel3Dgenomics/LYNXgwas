# Phase 5 decisions log — explicit disease field, region-level cross-dataset ANOVA,
# effect-size ANOVA, viewer toggle, white-canvas re-theme

Same contract as prior phases: written and updated incrementally and honestly. `[ ]` not done,
`[x]` done and verified, blockers stated plainly.

## 0. Scope of the request

1. Disease should be an explicit field the user provides per dataset, not inferred from the project
   id's naming convention.
2. A new cross-dataset analysis, structurally similar to Gene Constellation, comparing SNPs/genes/
   regions across and within diseases for variation in significance *and* effect size (OR/beta), by
   aligning genomic regions across datasets rather than only by shared nearest-gene label.
3. A "fancy circular" visualization for this, either standalone or as a mode the user switches to from
   within the existing Gene Constellation viewer (the user's own stated preference).
4. Re-theme Gene Constellation's canvas from dark to white/light, adjusting every color tuned
   specifically for the old dark background.

## 1. What already exists (checked before writing new code, not assumed)

- `Config.diseaseName` already exists as a per-project field (`disease.name` in config.properties) —
  it was simply never wired into the disease-grouping logic, which instead infers a group from the
  project id's prefix (`GeneConstellationBuilder.diseaseGroupOf`, "substring before first '-'").
- Cross-dataset **region-level pooling already exists**: `MultiLocusScanner`/`MultiLocusMerger` already
  pool loci across every dataset by genomic position into `MultiLocusResult.LocusRow` (chr/start/end),
  each carrying a `Map<datasetId, DatasetLocusStat>` with that dataset's own best SNP, p-value, and
  beta/OR for that pooled region. `GeneConstellationBuilder` takes this already-region-pooled structure
  and *collapses* it further to one node per nearest-gene. So "aligning regions across datasets" for
  significance is not new work — it's already the intermediate representation Gene Constellation itself
  is built from. What's genuinely new: (a) exposing that region-level granularity as its own view
  instead of only ever collapsing it to genes, and (b) an ANOVA on effect size, which does not exist at
  any granularity today (only -log10(p) is ANOVA'd currently).
- True **SNP-level** cross-dataset alignment (the same exact variant tested as significant in multiple
  datasets, not just each dataset's own single best SNP per region) does not exist and would require
  re-scanning raw GWAS rows per region across every dataset — comparable cost to the regulatory-
  enrichment SNP scan. Scoped out of this phase (section 6), not silently dropped.
- Gene Constellation's canvas already uses CSS custom properties (`--canvas-bg-1/2`, edge stroke, a
  glow filter) that are hard-coded dark regardless of the page's own light/dark theme — confirmed by
  reading the actual CSS, not guessed.

## 2. Design decisions resolved here

- **Explicit disease field, backward compatible**: `diseaseGroupOf`/`resolveEid` gain an overload that
  takes the project's `Config.diseaseName` and prefers it (trimmed, lower-cased for grouping) when
  non-empty; the existing id-prefix heuristic remains the fallback for projects that never set it
  (including the existing 30-dataset corpus), so nothing already built breaks.
- **Region-level node = one `LocusRow`, not one gene**: a new `RegionConstellationBuilder` mirrors
  `GeneConstellationBuilder`'s structure but does not collapse rows sharing a nearest gene into one
  node — each pooled region is its own node, labeled by its nearest gene when one exists and by
  `chrN:start-end` otherwise (so the 161-of-30-dataset "no assignable nearest gene" loci, currently
  reported separately and excluded from Gene Constellation entirely, get a place in this view instead
  of being dropped).
- **Effect-size ANOVA uses log(effect), not raw OR**: raw odds ratios are not symmetric around 1 (an OR
  of 0.5 and 2.0 are equal-magnitude opposite effects), so ANOVA needs `ln(OR)` (or `beta` directly,
  already on an additive scale) as its input, mirroring how `-log10(p)` is used for significance.
- **Honest, stated limitation, not hidden**: this phase does *not* attempt cross-study allele/effect-
  direction harmonization (a real, much bigger undertaking — this app's existing `AlleleHarmonizer`
  only harmonizes a single project's own SNPs against its own reference panel, not one dataset's
  effect-allele convention against another's). Effect-size ANOVA is computed on each dataset's own
  reported `ln(OR)`/`beta` as-is. This is scientifically fine for asking "how much does the *magnitude*
  of effect vary across/within disease groups" but not for claiming a consistent risk/protective
  direction across datasets — documented in the manuscript-facing limitations the same way the
  existing "background gene set" and "single-test-per-column" simplifications already are.
- **One viewer, a mode toggle** (the user's own preference): `gene_constellation.html` gains a
  Gene/Region toggle rather than being a separate page — same D3 circular layout code, fed either
  gene-nodes or region-nodes, since both shapes are structurally identical for rendering purposes
  (id/chr/pos/significance-count/direction-fraction/per-dataset breakdown/ANOVA summaries).
- **Re-theme scope**: only the canvas background and the specific colors proven (by reading the CSS)
  to have been tuned for it — `--canvas-bg-1/2`, the co-significance edge stroke, the node's dark-only
  glow filter, and the chromosome-ring stroke fallback. `d3.interpolateRdBu` (the risk/protective node
  fill scale) already works on any background and is left unchanged.

## 3. Work items and status

### 3.1 Explicit disease field [x] — done, verified
`MultiLocusResult.DatasetInfo` gains a `diseaseName` field, populated from `cfg.diseaseName` in
`MultiLocusScanner.scan()`. `GeneConstellationBuilder.diseaseGroupOf` gains an explicit-name-preferred
overload used for cross-dataset grouping (ANOVA/edges) — a free-text name is a perfectly good grouping
key there. **Deliberately not applied to `RegulatoryPeakIndex.resolveEid`**: its `DISEASE_TO_EID` table
is a small, fixed, short-code vocabulary (`scz`/`alz`/`t2d`/`cad`/`ra`/`ibd`) tied to which reference
epigenome to use, not a free grouping key — routing a user-typed name like "Schizophrenia" through it
would silently break the existing id-prefix-based regulatory-track lookup for no benefit, since a free
name has no natural mapping to that fixed vocabulary. Left unchanged, deliberately, not overlooked.
New Project wizard UI already had a "Disease / Trait name" field (`disease.name`) — no UI change needed.

### 3.2 `RegionConstellationBuilder` — region-level + effect-size ANOVA [x] — done, verified
New class, deliberately not a rewrite of `GeneConstellationBuilder` (different node identity), directly
reusing `AnovaUtil.oneWay` (already tested) and `GeneConstellationBuilder.effectValue` for the ln(effect)
scale. Also extended `GeneConstellationBuilder` itself with the same effect-size ANOVA at the gene
level (`betweenDiseaseAnovaEffect`/`withinDiseaseAnovaEffect` on `GeneEntry`), since "genes, regions...
for variation, OR" asked for effect-size comparison at both granularities, not only the new one.

**Real finding, discovered while testing, not decided in advance**: within-disease ANOVA is only
computed at the *gene* level, never at the *region* level. A gene can span several pooled regions, so
one dataset can contribute several observations to a within-disease-by-dataset group; a region is by
definition a single `LocusRow`, so every dataset contributes at most one observation to it — "group by
dataset within a disease" therefore always has exactly one observation per group (zero residual
degrees of freedom) at region granularity, not merely when data is sparse. The first version of this
code computed it anyway and a test caught the always-null result; removed the always-degenerate
computation entirely rather than leaving dead code whose only possible output is `null`, and documented
why directly in `RegionConstellationBuilder`'s class javadoc.

### 3.3 Result model + endpoint [x] — done, verified
`GeneConstellationResult` gained `RegionEntry` and a `regions` list, populated inside the *same*
`GeneConstellationBuilder.build()` call (calls `RegionConstellationBuilder.build(mlr, threshold)`
internally) rather than a second endpoint/mode parameter — one API response now carries both views, so
the frontend toggle (3.4) needs no second network round-trip. `LocalServer.geneConstellation`'s existing
job-running/caching infrastructure is unchanged and fully reused.

### 3.4 Viewer: Gene/Region toggle [x] — done, verified
`gene_constellation.html` gained a Genes/Regions header toggle. `normalizeForMode()` adapts either
`data.genes` or `data.regions` into one common shape so the entire existing circular-layout renderer,
tooltip, and detail panel are reused completely unchanged for both modes — not a second renderer.
Region mode always shows zero co-significance edges (a gene-level construct, out of scope for regions
per 3.2's own finding) rather than a fabricated edge list. Detail panel extended with the new
effect-size ANOVA sections (both modes) and a plain-language note when within-disease ANOVA is
structurally unavailable (region mode only).

### 3.5 White-canvas re-theme [x] — done, verified
`--canvas-bg-1/2` → light (`#FAFBFC`/`#FFFFFF`, was `#11141c`/`#1b2033`). Every color actually tuned
for the old dark backdrop (checked by reading the CSS/JS, not assumed) was re-tuned alongside it: the
co-significance edge stroke (`#E8ECF5` → `#8C9BC4`, otherwise near-invisible on white), the node's glow
filter (was a blur-behind-the-shape halo that only reads as a "glow" against a dark backdrop — replaced
with a proper drop shadow, the light-background equivalent), and the node stroke fallback (`#fff` →
`#1B2430`, otherwise invisible against the new white canvas for the opposite reason it was invisible
against the old dark one). `d3.interpolateRdBu` (the risk/protective node fill scale) needed no change,
confirmed by reading its usage, not assumed.

### 3.6 Testing [x] — done, verified
`tests/RegionConstellationBuilderTest.java` (19 cases): region grouping including a no-nearest-gene
row, explicit-disease-name grouping (deliberately using non-hyphenated dataset ids so id-prefix
inference cannot silently make the test pass for the wrong reason), between-disease ANOVA on both
scales cross-checked against a direct `AnovaUtil.oneWay` call on the same hand-derived groups, and the
within-disease-always-empty-at-region-level finding from 3.2. Two real fail-before-fix/pass-after-fix
cycles run during this phase: the explicit-disease-name wiring (reverted to id-prefix-only, the
between-disease-ANOVA assertion failed as expected) and the always-null within-disease computation
(the test's *own* wrong initial expectation caught the structural reality, not a code bug — see 3.2).
Clean-room `rm -rf bin && build.bat` + full `run_tests.bat` (including both new suites) pass together.
Visual verification used a Playwright fixture test (network-intercepted, no real backend needed) driving
the actual production `gene_constellation.html` through both modes, a node-click, and the detail panel —
zero console/page errors, screenshots confirm the toggle, title, node count, center label, and both
ANOVA sections all update correctly between modes. Also ran one real (not fixture) locus-matrix job on
the article server against `scz-demo`+`t2d-demo`: real `genes` (238) and `edges` (36) came back non-empty
as expected, but `regions` came back empty — because that server process predates this phase's compiled
code (same restart-required situation as every other backend-code change this session; not a bug in the
new code, confirmed by reading the response and recognizing the stale-process pattern rather than
assuming success). **Not yet re-verified against a real job on a restarted server** — pending the same
user go-ahead to restart as every other backend deployment this session.

## 4. Hard boundaries (carried over, unchanged)

No push without explicit instruction; every fork's report independently verified before being trusted
(n/a this phase — done directly); the real production/article servers are never restarted without
asking first.

## 5. Closing status

All four scope items (section 0) are done and independently verified at the code level: disease is now
an explicit, wired-through field for cross-dataset grouping; a real region-level cross-dataset view
exists alongside the gene-level one, with effect-size (ln OR/beta) ANOVA added at both granularities on
top of the pre-existing significance ANOVA; one viewer with a mode toggle, per the user's own stated
preference, rather than a second page; and the canvas is white with every dark-tuned color re-tuned
alongside it, not just the background swapped. 19 new tests (`RegionConstellationBuilderTest`) plus the
existing suite all pass together after a clean-room rebuild. A real structural finding — within-disease
ANOVA is impossible at region granularity by construction, not just sometimes null — was caught by a
test's own wrong assumption during development and fixed in the implementation, not papered over.

**Deployment status, honestly reported, not implied complete**: like every prior phase's backend
changes, this is compiled Java code the already-running production and article servers have not loaded
yet (confirmed directly: a real job on the still-running article server returns real genes/edges but
empty regions). The frontend toggle/re-theme are static-file changes and are already live on both
servers. Restarting the backend on each is the user's call, not automatically done here.
