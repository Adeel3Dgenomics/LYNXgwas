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

### 3.1 Backend: regulatory peak ingestion + interval index [ ]
New ingestion for BED/narrowPeak-format evidence (chrom, start, end, name, score, strand,
signalValue, pValue, qValue, peak) — a genuinely different shape from the existing gene-symbol-keyed
evidence tables (3.8 in Phase 2), so this is new code, not a reuse of `EnrichmentAnalyzer`'s existing
per-gene path. Peaks resolved per-project by the project's disease (derived the same
best-effort way `GeneConstellationBuilder.diseaseGroupOf` already does) against the shared
`regulatory_data/<EID>/` location, not re-uploaded per project.

### 3.2 Frontend: Option C stacked multi-mark tracks in viewer.html [ ]
One labeled row per mark (H3K27ac / H3K4me1 / H3K4me3), stacked under the gene track, matching the
design canvas's Option C mockup, restricted to whatever peaks overlap the current locus's window.

### 3.3 Real locus-based regulatory enrichment test [ ]
For each dataset: of its genome-wide-significant SNPs (foreground) vs. a background SNP set from the
same GWAS file, what fraction overlaps a peak of each mark, and is that enrichment significant
(Fisher's exact, reusing `EnrichmentAnalyzer`'s existing exact hypergeometric implementation with a
new interval-overlap classifier in place of its existing per-gene-value lookup)?

### 3.4 Real 30-dataset run [ ]
Run 3.3 against all 30 already-processed datasets using their disease-mapped tissue's real peaks.
Aggregate per-disease and overall. This is the analysis the manuscript will report.

### 3.5 Most-important-result plot [ ]
Chosen only after seeing the real 3.4 numbers (not decided in advance) — most likely a cross-disease
bar chart of regulatory-overlap enrichment (odds ratio or fold-enrichment with confidence interval,
one bar per disease, per mark or combined) if the real result supports that story; documented
honestly if the real numbers suggest a different, more accurate framing instead.

### 3.6 Manuscript update [ ]
New Methods subsection (regulatory data source, tissue mapping table, enrichment test), new Results
subsection with the real 3.4 numbers, the 3.5 figure, and a Limitations addition for the
reference-epigenome-as-proxy caveat from section 1.

## 4. Hard boundaries (carried over, unchanged)

No push without explicit instruction beyond what's already been authorized; no writes to real
`projects/*` beyond the evidence/peak-index additions this feature itself requires; every fork's
report independently verified before being trusted; the real production server (currently PID —
recheck at execution time) is never restarted without asking first, exactly as established in
Phase 2's search-feature rollout.
