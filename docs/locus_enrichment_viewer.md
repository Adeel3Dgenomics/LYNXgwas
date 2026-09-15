# Locus Enrichment Viewer — Design Discussion

Living design doc for evolving the current **Locus Matrix** (`locus_matrix.html`) into a
**Locus Enrichment Viewer**: a cross-dataset, multi-evidence-type workspace for known-loci-based
comparison (GWAS layers today; gene enrichment, pathway, and fine-mapping/PIP layers planned).

Status: **design discussion — no code written yet.** This file is updated as the discussion
progresses; see the changelog at the bottom for what changed and when.

---

## 1. Why change it

The current Locus Matrix (`locus_matrix.html` + `MultiLocusMerger`/`MultiLocusScanner`/
`MultiLocusResult` + `/api/locus-matrix-*` endpoints in `LocalServer.java`) is a
dataset × locus **presence** grid: rows are LYNXgwas projects, columns are loci, and a cell is
colored by chromosome if the dataset has any SNP in that locus, or gray/white-striped if not.
The actually useful numbers (`best_p`, `beta`/`or`, `se`, `n_p5e5`, `n_p5e8`, nearest gene) only
show up in a hover tooltip, one cell at a time.

Two limitations drove this redesign:

1. **It's meta-analysis under the hood, not locus comparison.** `MultiLocusMerger` pools the
   best p-value per position across every selected dataset's raw GWAS file and *re-clumps* that
   pooled signal into new loci via `LociIdentifier`. That's fine for finding a shared significant
   region, but it's not "compare these datasets at their own already-identified loci."
2. **It only knows one evidence type.** SNP-level GWAS stats. There's no way to bring in gene-level
   enrichment, pathway enrichment, or fine-mapping (PIP) evidence alongside the GWAS comparison.

## 2. Scope for this pass

Explicitly **in scope** now:
- Cross-dataset GWAS layer comparison, aligned on each dataset's own pre-identified loci
  (not re-clumped).
- A Gene Enrichment tab (GO Biological Process / Molecular Function / Cellular Component,
  plus a Protein–Protein Interaction sub-tab) — UI skeleton and data model; the actual
  compute-vs-import decision is still open (see §5).

Explicitly **deferred** ("not now" per discussion):
- SNP-PIP / fine-mapping layers (whether pulled from an existing project's SuSiE/FINEMAP/COJO
  run, or uploaded externally).
- Pathway (`snp_pathway` / `gene_pathway`) layers.
- Cell-level rendering detail for "multiple genes in one locus" (best-gene-only vs. a ranked
  list) — deferred until the Gene Enrichment tab is actually being built.

Decided and locked in:
- **GWAS layer ingestion is existing-projects-only.** No raw-file quick-add / one-off column
  mapping for this view — a dataset must already be a configured LYNXgwas project.

## 3. Locus alignment model

Replaces `MultiLocusMerger`'s pool-and-reclump behavior for this view.

- User selects one or more **existing projects that already have identified loci**
  (`loci.txt`/`loci_detail.tsv` present). Projects without loci are excluded from the picker —
  loci are never auto-computed for this view; "if the loci isn't defined, it isn't detected."
- **Merge algorithm**: union-merge each selected project's own locus intervals wherever they
  overlap (or fall within a merge-distance threshold — see open question in §5). A merged
  locus's displayed span is the min-start/max-end envelope of its contributing intervals. This
  reuses the same interval-merge logic `LociIdentifier` already applies within a single project,
  just applied across datasets.
- Each merged locus tracks, per contributing dataset: whether it's *supported* there, plus that
  dataset's own best-SNP/p/beta/se/etc. (same per-cell stats the matrix already computes today).
- **Aggregate summary value**: one blended stat per merged locus across whichever datasets
  support it (e.g. mean of −log10(p), or an N-weighted mean beta). This is explicitly a
  *summary*, not a formal meta-analysis p-value — `tools/gwama_meta.yaml` already exists for
  real fixed/random-effects meta-analysis if that's ever wanted instead of this.

## 4. Layout

Three regions, deliberately reusing patterns `viewer.html` already has working rather than
inventing new UI conventions.

**Left — loci sidebar.** Same pattern as the existing per-project viewer sidebar: searchable
list of merged loci, each row showing a support count (e.g. "4/6 datasets") and the best signal
across supporting datasets.

**Top — genome-wide overview strip.** One band per merged locus, ordered by genomic position,
color-coded by chromosome, showing which datasets support each locus. Key change from today's
matrix: **not** proportional-to-bp width. Today's proportional width is exactly why small loci
visually disappear next to megabase-scale ones. Switching to equal-width columns trades literal
physical proportionality for legibility, since the strip's job is "which loci, supported by
whom," not physical distance. This strip has its own tab bar; the first tab is **Across-GWAS**
(functionally today's matrix, redrawn with the new width model).

**Bottom — locus detail workspace**, driven by whichever locus is selected in the sidebar or the
overview strip. Tabbed:
- **GWAS** tab — per-dataset regional/Manhattan comparison for just the selected locus (the
  per-locus counterpart to the overview strip's Across-GWAS tab).
- **Gene Enrichment** tab — sub-tabs for the three standard GO aspects (Biological Process,
  Molecular Function, Cellular Component) plus a Protein–Protein Interaction sub-tab
  (STRING-based).
- *(reserved, not built now)* Pathway tab, SNP-PIP/fine-mapping tab. The tab container should be
  generic enough that adding these later doesn't require restructuring.

## 5. Open questions

Unresolved as of the last update — needed before implementation can start on the affected part:

1. **Gene Enrichment data source.** Does LYNXgwas *compute* GO/PPI enrichment itself (calling an
   enrichment method/API against the genes in each locus), or *import* enrichment results the
   user already ran elsewhere (STRING/g:Profiler/etc. export)? This decides whether we're
   building a computation pipeline or an ingestion+display feature — changes the backend shape
   substantially.
2. **Merge distance for cross-dataset locus alignment.** Reuse each project's own
   `merge.distance.bp` config value, or introduce a separate matrix-level setting the user sets
   per run?
3. **Aggregate summary weighting.** N-weighted mean, or unweighted mean, for the blended
   per-locus summary stat described in §3?

Deferred questions (revisit when their features come into scope):
- SNP-PIP source(s): pull from existing project analysis runs, external upload, or both.
- Multi-gene-per-locus cell rendering: best gene only vs. a small ranked list.

## 6. Changelog

- **2026-08-20** — Initial version. Captures the discussion that moved this from "improve the
  existing Locus Matrix visualization" to "generalize it into a multi-layer Locus Enrichment
  Viewer." Locked: existing-projects-only GWAS ingestion. Deferred: PIP layers, pathway layers,
  multi-gene cell rendering. Three open questions recorded in §5.
