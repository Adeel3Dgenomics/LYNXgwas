# Phase 7 decisions log — Locus Serpent Plot + category analysis

Same contract as the earlier `DECISIONS*.md` files: written as the work happened, `[x]` = done and
verified, limitations stated plainly.

## 0. Scope of the request

1. Integrate the Locus Serpent Plot (prototype: `prototypes/loci_snake/loci_snake.html`, simulated data)
   into LYNXgwas and add it to the home page.
2. Runs work like Gene Constellation: the user selects datasets, picks a reference panel, and runs it
   on real data.
3. Add a category-based analysis: which dataset category (disease, ancestry, ...) best explains the
   clustering of datasets by consistently significant loci, and which loci drive it.

## 1. Design

- [x] **No second pipeline.** A Serpent run is a normal Locus Matrix job (merge -> identify -> scan)
  tagged `kind: "serpent"` (`LocusMatrixJobMeta.kind`, `POST /api/locus-matrix-run {"kind":"serpent"}`).
  Both views read the same `MultiLocusResult`, so every finished run links to the other view.
- [x] **`GET /api/serpent-plot?job=`** returns the derived view from `SerpentPlotBuilder` (pure,
  read-only), cached per job because the category analysis runs permutations.
- [x] **Dataset metadata** added to `MultiLocusResult.DatasetInfo`: `ancestry` (Config.ancestry, which
  falls back to the ref-panel population), `sampleN`, and free-form `categories` from new
  `category.<name>=<value>` config lines (read and written by `Config`).
- [x] **Tiers**: 0 `p >= 5e-3` or missing, 1 `[5e-5, 5e-3)`, 2 `[5e-8, 5e-5)`, 3 `< 5e-8`, the
  plot's grey/green/yellow/red steps.
- [x] **Offset**: dataset best-SNP position minus the position of the locus's overall best SNP.
- [x] **Distance**: cosine distance between datasets' tier vectors (all loci).
- [x] **Clustering**: UPGMA (Lance-Williams), K = best mean silhouette in [2, n/2].
- [x] **Sparse datasets excluded from clustering/tests**: datasets with fewer than 3 genome-wide
  significant loci. First real run without this rule: silhouette chose K=2 (one X-chromosome-only
  dataset vs everyone), ARI vs disease 0.00. A simulation on the real tier matrix (scipy/sklearn,
  independent of the Java code) showed that excluding the 5 sparse datasets gives K=8, ARI 0.70.
  They are still drawn, in their own "Too few sig. loci" band, and named on the page.
- [x] **Category tests**, per category (Disease, Ancestry, Sample-size tertile as a power control,
  any custom `category.*`):
  PERMANOVA (R^2, pseudo-F, 999 label permutations) on the same distance matrix; adjusted Rand
  index of the category vs the K clusters; per level, cohesion (1 - mean within distance),
  separation, and **consistently and specifically significant loci**: p < 5e-8 in >= 50% (and >= 2)
  of the level's datasets, more often than elsewhere, two-sided Fisher p < 0.05 (uncorrected).
  The Fisher filter was added after the first real run listed PBX2 (the extended-MHC artifact,
  significant almost everywhere) as RA-specific with Fisher p = 1.0.
- [x] A category whose only replicated level is one group (e.g. 25 EUR vs 1 multi-ancestry) is still
  computed but flagged as uninformative, not as "no effect".

## 2. Verification

- [x] `tests/SerpentPlotBuilderTest.java` (added to `run_tests.bat`/`.sh`): tier boundaries;
  PERMANOVA == centroid sums of squares for Euclidean 1-D points; permutation p vs the exact 2/70;
  ARI identical/relabeled/hand example (0.8/3.3); UPGMA heights by hand; silhouette by hand;
  end-to-end synthetic corpus (disease R^2 = 1, ancestry R^2 = 0, ARI = 1, consistent loci exact,
  Fisher p equals `EnrichmentAnalyzer`); outlier robustness; everywhere-significant locus excluded.
  The outlier test failed against the first cluster-count rule (fixed K = sqrt(n) non-singletons),
  confirming it is not vacuous. Full suite passes after a clean `build.bat`.
- [x] End-to-end on the real 31-dataset / 7-disease corpus (GRCh37, EUR 1000G panel), started
  through the same API the home-page button uses: 704 loci; serpent JSON in ~0.3 s.

## 3. Real-data result (31 datasets, 26 clustered)

- Disease explains **76%** of the shared-loci variation (PERMANOVA p = 0.001); clusters vs disease
  ARI **0.70** (K = 8): all 5 SCZ datasets in one cluster, all 5 CAD in one, RA 3 of 4, IBD's three
  large studies together; one mixed cluster (AD + the two smaller IBD studies + sarcoidosis).
- Sample size explains 15% (p = 0.021): part of the structure is statistical power. The page says so.
- Ancestry: uninformative in this corpus (25 EUR vs 1 multi-ancestry).
- Consistent loci per disease are the expected ones: CAD PCSK9, CELSR2, MIA3, PLPP3; IBD NOD2,
  IL23R, IKZF1; RA PTPN22, ANKRD55; T2D TCF7L2, HHEX, KCNJ11, WFS1, PPARG, CDKAL1; AD CD2AP, MS4A,
  CR1; SCZ AKT3, VRK2, TCF4, CACNA1C, NT5C2.

## 4. Known limitations

- ~~Runs are in memory and disappear on restart~~ fixed, see section 5.
- Whole-genome view of hundreds of loci is dense; large runs open filtered to loci significant in
  >= 2 diseases, with chromosome zoom for detail.
- Fisher p-values for consistent loci are not multiplicity-corrected; the clusters and R^2 are
  exploratory, not a genetic-correlation estimate; datasets of one disease can share samples.
- The pip bundle (`python/lynxgwas/_bundled`) is not refreshed yet; `cli.py` already lists
  `serpent_plot.html`.

## 5. Persisted runs (Gene Constellation and Serpent)

- [x] Every Locus Matrix run is saved in its own folder, `output/multi_locus/<jobId>/`, next to the
  `merged.tsv` the pipeline already wrote there: `job.json` (name, kind, project ids, dataset names,
  ref panel, created, status running/done/error, error message) is written at start and updated at
  the end; `result.json` (`MultiLocusResult.toJson()`) is written when the run finishes. Both are
  written to a temp file and renamed, so a crash cannot leave a truncated file.
- [x] At startup `LocalServer` reloads every folder with a `job.json` (`MultiLocusResult.fromJson`):
  finished runs come back fully viewable; a run still `running` when LYNXgwas stopped comes back as
  "Interrupted" with a delete-and-rerun message; folders from before persistence (only `merged.tsv`)
  are skipped because their result was never saved.
- [x] Delete removes the whole run folder.
- [x] `MultiLocusResult.toJson()` now escapes control characters (it previously escaped only `\` and
  `"`, so a tab or newline in a dataset's extra column would have produced invalid JSON) and includes
  dataset categories.
- [x] `tests/MultiLocusResultPersistenceTest.java`: toJson -> fromJson -> toJson byte-identical
  (NaN, categories, ancestry, sample size, awkward characters), and the Serpent view built from the
  reloaded result is identical to the original.
