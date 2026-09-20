# Phase 3 decisions log — regulatory/functional genomics integration (Option C) + real
# cross-disease enrichment analysis

Same contract as `DECISIONS.md`/`DECISIONS_PHASE2.md`: written and updated incrementally and
honestly. `[ ]` not done, `[x]` done and verified, blockers stated plainly.

## 0. Scope of the request

1. Implement **Option C** (multi-mark stacked regulatory tracks) from the four locus-level view
   concepts presented on the design canvas — real feature, not a mockup.
2. Apply it to the real 30-dataset/6-disease corpus using **real downloaded regulatory data**, not
   synthetic peaks.
3. Run a real enrichment analysis and report it in the manuscript as evidence of what LYNXgwas's
   new capability adds.
4. Save storage space as much as possible while doing this.
5. Generate a plot of the single most important result and add it to the manuscript.

## 1. Real data source (verified before writing any code)

**NIH Roadmap Epigenomics**, consolidated narrowPeak histone ChIP-seq peak calls, hosted at
`https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/` — confirmed live via
a real directory listing and `HEAD` requests (not assumed from memory). Chosen over ENCODE's raw
experiment accessions for this task because Roadmap already publishes one *consolidated* peak call
per reference epigenome per mark — exactly the shape this feature needs — with a small, uniform,
well-documented naming scheme (`<EID>-<mark>.narrowPeak.gz`).

**Disease → reference epigenome (EID) mapping**, each independently cross-checked against at least
one real source before use (not guessed):
| Disease | EID | Epigenome | Marks used |
|---|---|---|---|
| Schizophrenia | E073 | Brain, Dorsolateral Prefrontal Cortex | H3K27ac, H3K4me1, H3K4me3 |
| Alzheimer's disease | E073 | Brain, Dorsolateral Prefrontal Cortex (shared with SCZ — both neuropsychiatric/neurodegenerative, same reference tissue) | H3K27ac, H3K4me1, H3K4me3 |
| Type 2 diabetes | E098 | Pancreas | H3K27ac, H3K4me1, H3K4me3 |
| Coronary artery disease | E065 | Aorta | H3K27ac, H3K4me1, H3K4me3 |
| Rheumatoid arthritis | E034 | Primary T cells, peripheral blood | H3K27ac, H3K4me1, H3K4me3 |
| Crohn's / IBD | E106 | Sigmoid colon | H3K27ac, H3K4me1, H3K4me3 |

All 15 files (5 tissues × 3 marks) confirmed to actually exist via a live directory listing before
committing to this table. Real gzipped sizes checked via `HEAD`: 1.3–3.2 MB each, ~25 MB total for
all 15 — confirmed cheap before downloading anything.

**Honest limitation, stated up front, not discovered later**: a reference epigenome is a real,
publicly generated ChIP-seq dataset from real human tissue, but it is not the GWAS cohort's own
tissue, and one epigenome stands in for an entire disease category (e.g., "Aorta" for all of
coronary artery disease). This is a standard, defensible proxy in the field (the same logic behind
tools like GREGOR, LDSC-SEG, and FUMA's tissue-enrichment modules), but it is a proxy, not ground
truth, and is written up as such in the manuscript rather than overclaimed.

## 2. Storage-saving design (addressing the explicit request)

- **narrowPeak, not bigWig**: peak calls (~1–3 MB/file) instead of continuous signal tracks
  (hundreds of MB–GB each) — a 100–1000x reduction for the same biological question ("where are the
  active regulatory elements"), since this feature only needs peak locations, not signal shape.
- **Stored once per tissue, not once per project**: 30 datasets map to only 5 unique tissues (SCZ
  and Alzheimer's share E073). Files live in one shared `regulatory_data/<EID>/` location the app
  resolves by disease, not copied into each of the 30 `projects/<id>/evidence/` folders — a ~6x
  reduction over the naive per-project-copy approach on top of the format saving above.
- Downloaded files are **kept gzipped on disk**; parsed once into an in-memory interval index at
  query time (matching this app's existing `GwasCatalogLocalIndex`/`GlobalSearchIndex` convention of
  parse-once-cache-in-memory rather than pre-expanding to a larger on-disk form).
- Intermediate per-dataset analysis artifacts from the cross-disease enrichment run (section 5) are
  deleted after the final summary numbers/plot are extracted and verified, keeping only: the 15
  source peak files (~25 MB), the final summary JSON/CSV, and the generated figure.

## 3. Work items and status

### 3.0 Real data download [x] — done, verified
All 15 files (5 tissues × 3 marks) downloaded for real from the live Roadmap Epigenomics server into
`regulatory_data/<EID>/` (gitignored, matching `gwascatalog_data/`'s existing convention). Verified
non-corrupt: `zcat`'d and confirmed real 10-column narrowPeak rows with plausible chromosome/position
values (e.g. E073-H3K27ac: 206,746 real peaks). Total on-disk size 31 MB for all 15 files — confirms
the storage-saving design target from section 2 (narrowPeak instead of bigWig, one shared copy per
tissue instead of per project) up front, before any feature code exists to consume them.
Citation for the manuscript: Roadmap Epigenomics Consortium et al. (2015), *Nature* 518:317–330.

### 3.1 Backend: regulatory peak ingestion + interval index [x] — done, verified
`src/analysis/RegulatoryPeakIndex.java`: lazy per-(EID,mark) in-memory interval index over the
narrowPeak files, same normalizeChr+binary-search convention as `GwasParser`/`MultiLocusScanner`;
holds the disease→EID and EID→tissue-name tables; handles missing/corrupt files gracefully (logs,
returns empty rather than throwing). `LocalServer.java` exposes it via
`GET /api/project/{id}/locus/{n}/regulatory`. Covered by `tests/RegulatoryPeakIndexTest.java`
(overlap-boundary, missing-file, unknown-EID/mark, corrupt-gzip cases) — confirmed fail-before-fix
(broke the overlap boundary check `p.end > start` vs `>=`) / pass-after-fix.

### 3.2 Frontend: Option C stacked multi-mark tracks in viewer.html [x] — done, verified
`viewer.html`: one labeled row per mark (H3K27ac / H3K4me1 / H3K4me3), stacked directly under the
gene track above the evidence track, matching the design canvas's Option C mockup; renders as
`{marks:{}}` (i.e. nothing) for any project whose disease has no mapped tissue, fully
backward-compatible with every pre-Phase-3 project. Fetched once per locus load via the new
`/locus/{n}/regulatory` endpoint with a 4s timeout guard.

### 3.3 Real locus-based regulatory enrichment test [x] — done, verified
`src/analysis/RegulatoryEnrichmentAnalyzer.java`: real interval-overlap Fisher's exact test —
reuses `EnrichmentAnalyzer.fisherExactTwoSided(a,b,c,d)` directly (confirmed via `grep`, not
reimplemented), adding only the interval-overlap classification step and an odds ratio with
Haldane-Anscombe correction for zero cells. Exposed via
`GET /api/project/{id}/regulatory-enrichment?threshold=5e-8`. Covered by
`tests/RegulatoryEnrichmentAnalyzerTest.java` — confirmed fail-before-fix (broke the odds-ratio
formula `a*c/(b*d)` vs. the correct `a*d/(b*c)`) / pass-after-fix.

### 3.4 Real 30-dataset run [x] — done, verified
Ran the real `/api/project/{id}/regulatory-enrichment` endpoint against all 30 real, already-fully-
processed datasets on a throwaway copy of the `scz-full-demo` workspace (original PID-32928
production server confirmed untouched throughout, same start time before/after). 4 of 30 datasets
(`alz-src-belloy2024`, `t2d-src-saxena2007`, `t2d-src-wood2016`, `ra-src-shigesi2025`) have
`n_foreground_total` < 20 genome-wide-significant SNPs (matching what Phase 2 already documented
about these specific datasets being small/underpowered) and were excluded from pooled statistics
before computing anything, not after seeing an inconvenient number — reported individually, not
silently dropped. Raw per-dataset counts (foreground-in-peak/total, background-in-peak/total) summed
across each disease's remaining datasets, then ONE Fisher's exact + odds ratio computed per disease
per mark from the pooled counts (statistically sounder than averaging 30 separate p-values/ORs).

**Cross-checked two independent ways**: (1) the Java implementation's own output for `scz-demo`
H3K27ac (OR=2.129777, p=2.422244e-147) recomputed from the identical raw counts using Python's
`scipy.stats.fisher_exact` — exact agreement to 6 decimal places; (2) the earlier delegated build
already cross-checked one project's foreground-in-peak count against an independent from-scratch
`awk` script (17/255 exact match) before this 30-dataset run even started.

**Real result — a consistent, strong, real enrichment signal across all 6 diseases**: pooling 26
included datasets, GWAS-significant SNPs are 1.6–4.6x more likely to overlap an active regulatory
mark in the disease-relevant tissue than background SNPs, in every one of 6 diseases × 3 marks (18/18
comparisons), every one below $p=1\times10^{-4}$ (several underflow below double precision, i.e.
$p<10^{-300}$ exactly, not literally zero). Overall pooled (all 6 diseases, 180,942 foreground SNPs
vs. 269,015,376 background SNPs): H3K27ac OR=2.63, H3K4me1 OR=2.85, H3K4me3 OR=3.28, all
$p\to 0$. Per-disease range: SCZ/T2D/CAD/ALZ cluster around OR 1.6–2.8; RA and IBD (immune/gut
tissues) show the strongest enrichment, OR up to 4.6 (RA, H3K4me3).

### 3.5 Most-important-result plot [x] — done, verified
`figures/regulatory_enrichment_6disease.png`: grouped bar chart, one group per disease, one bar per
mark, real pooled odds ratios with significance stars, generated directly from the verified 3.4
numbers (not illustrative). Visually confirmed: every bar clears the OR=1 reference line, all
significant, tells the "consistent across independent diseases" story the real numbers actually
support — chosen after seeing the numbers, not decided in advance.

### 3.6 Manuscript update [x] — done, verified
`LYNXgwas-paper/main.tex`: new Methods §2.9 "Regulatory and functional genomics integration"
(data source, storage-saving rationale, `tab:tissuemap` disease→EID→tissue table, Option C viewer
description, enrichment-test methodology explicitly noting reuse of `EnrichmentAnalyzer`'s existing
Fisher's-exact function); new Results §3.4.2 "Regulatory-element enrichment across all six diseases"
with the real pooled per-disease/overall numbers from 3.4 and the `regulatory_enrichment_6disease.png`
figure from 3.5; a new Limitations bullet covering the reference-epigenome-as-proxy caveat and the
pre-stated small-N exclusion threshold; Abstract, Conclusion, and Acknowledgments updated to reference
the new capability and headline result. `references.bib` gained the Roadmap Epigenomics citation
(verified via WebSearch: Nature 518(7539):317–330, 2015, DOI 10.1038/nature14248). Full
`pdflatex`+`bibtex`+`pdflatex`×2 compile cycle run clean (no undefined references/citations); visual
QA performed by rendering the title/author page, the Methods §2.9 page, both new Results pages
(Gene Constellation figure + regulatory-enrichment figure), and the Limitations/Conclusion pages to
PNG and inspecting them directly — all cross-references, the new table, and both figures render
correctly. Committed as `f848216` in `LYNXgwas-paper` (not pushed — no push instruction given for
this phase).

## 4. Hard boundaries (carried over, unchanged)

No push without explicit instruction beyond what's already been authorized; no writes to real
`projects/*` beyond the evidence/peak-index additions this feature itself requires; every fork's
report independently verified before being trusted; the real production server (currently PID —
recheck at execution time) is never restarted without asking first, exactly as established in
Phase 2's search-feature rollout.

## 5. Closing status

All five items in the original scope of request (section 0) are done and independently verified:
Option C is a real, working feature (not a mockup); the real 30-dataset/6-disease corpus was analyzed
end-to-end using real downloaded Roadmap Epigenomics data; the enrichment analysis and its headline
result are written up in the manuscript as concrete evidence of what this capability adds; storage was
minimized at every layer (narrowPeak not bigWig, one shared copy per tissue not per project, gzipped
on disk, intermediate per-dataset artifacts deleted after the final summary/figure were extracted); and
the most-important-result plot is generated and embedded in the manuscript. A clean-room rebuild
(`rm -rf bin && build.bat`) and the full existing test suite (including both new Phase-3 suites) were
re-run and confirmed passing immediately before closing this phase. LYNXgwas and LYNXgwas-paper commits
for this phase are local only — neither repo was pushed, matching every hard boundary above.
