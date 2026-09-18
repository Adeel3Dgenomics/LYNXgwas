# Autonomous execution log — limitations fixes + HPC/shared-storage features

Started 2026-09-17, while the user is away from the computer. This file records the plan and every
non-obvious decision made while executing it, per the user's explicit request. Updated as work
progresses — check timestamps/checkmarks below for current status, not just the plan text.

## Scope, as requested

1. Fix each limitation listed in the manuscript's Limitations section, one by one, where a real
   autonomous fix is safe and verifiable.
2. Add: shared-storage options for reference data (private GitHub repo, shared lab folder, cloud-synced
   folder) so users don't have to re-download PLINK/GFF3/reference panels every time.
3. Add: optional HPC/cluster execution — user provides SSH access, the app generates and submits a
   Snakemake workflow remotely, polls job status on a configurable interval (5–20 min), using
   user-specified module names for cluster-side tools. Must be **opt-in** (local execution stays the
   default) and safe.
4. Add: a GUI setup/configuration wizard for all of the above, re-runnable any time from the home page.
5. Reprocess all 6 diseases × 5 datasets (30 total) with **every** genome-wide-significant locus fully
   processed (not just the ~14-locus illustrative subset used for the manuscript figures).

## Hard safety boundary set for this session (decision, not negotiable)

**No autonomous SSH connection to any real, external, production system** (including the user's own
HPC cluster referenced elsewhere in this environment) will be attempted while unsupervised. Reasons:
- The user's own standing instructions elsewhere in this environment describe a cluster incident where
  a single bare command on a login node once cascaded into hung sessions — i.e. there is a specific,
  known-bad failure mode for unsupervised SSH/cluster interaction with this exact user.
- Autonomous, unattended use of SSH credentials against a real external system is exactly the class of
  action the general operating instructions require a human present for.

Consequence: the HPC feature is built completely (config, GUI, Snakemake generation, SSH invocation via
the system `ssh`/`scp` binaries using key-based auth only, job status polling) and tested as far as
possible **without** a live external target (config round-trip, generated-file correctness, and a
loopback test against `localhost` if an SSH server is available on this machine). Connecting it to the
user's actual cluster is left as a manual step for when they're back and can supervise the first real
run.

## Per-limitation decisions

| # | Limitation (from manuscript) | Decision | Status |
|---|---|---|---|
| 1 | Parallel LD memory scaling (OOM crash) | **Done.** Found the actual root cause while implementing this: `--ld-window-kb` for the index-SNP r2 scan (Step B) scaled directly with the locus's *raw genomic span*, uncapped — an 8.8Mb MHC-sized locus drove PLINK into a multi-thousand-kb window against tens of thousands of markers. Added `MAX_LD_WINDOW_KB` (2000) cap — bounds the actual memory/runtime driver directly, not just a band-aid on thread count. Also added a lighter defense-in-depth: the worker pool drops to at most 2 concurrent jobs (from whatever `ld.parallel.jobs` is configured) whenever any locus in the batch spans over 4Mb. Verified with `tests/LdCalculatorTest.java`. | ✅ done, tested |
| 2 | Fine-mapping method fidelity (ABF vs real FINEMAP) | **Deferred.** Integrating the real FINEMAP binary is a substantial new external-tool integration (binary distribution, licensing check, new adapter), not a safe same-day autonomous change. Left as future work. | deferred |
| 3 | Meta-analysis scope (GWAMA single-cohort only) | **Deferred.** True multi-cohort meta-analysis is a statistical-correctness-critical feature; implementing it without domain review risks silently wrong results, worse than the current honest limitation. Left as future work. | deferred |
| 4 | QC signal not LD-pruned | **Done, but scoped down from the original plan.** True `PLINK --indep-pairwise` pruning needs genotypes/a reference panel, which would break `GwasQc`'s explicit no-external-tool design goal (many projects run QC without ever configuring a reference panel). Implemented **distance-based** pruning instead (keep at most 1 SNP per 250kb per chromosome) — a genotype-free proxy that still substantially reduces dense-LD-block inflation of the median, honestly documented as approximate rather than true LD-based pruning. Also fixed a real cache-staleness bug found while implementing this: the QC result cache was keyed only by GWAS-file hash, so a pre-fix cached result would have been served forever after this change (same file, same hash) — added a cache-format marker so stale caches are recomputed. Verified with `tests/GwasQcTest.java`. | ✅ done, tested |
| 5 | External dependency management (no bundled ref panel/GFF3) | **Done** — see "Shared storage + HPC feature build" below for the full implementation (`SharedStorageResolver`, GUI section, test endpoint). This doesn't eliminate the need for a reference panel/GFF3 somewhere, but removes the "re-download per machine/project" pain the limitation was about. | ✅ done, tested |
| 6 | No automated test suite | **Started, not complete.** This project has no build tool (plain `javac`, no Maven/Gradle), so a full JUnit setup would be new infrastructure; added a dependency-free `tests/` directory (standalone classes with `main()` methods, PASS/FAIL output, non-zero exit on failure) plus `run_tests.bat`/`run_tests.sh`, and wrote real regression tests for every fix made today (3 test files, 8 assertions total, all passing). This is a genuine start, not a comprehensive suite — most existing modules (rsID recovery, fine-mapping adapters, export) still have zero coverage. | 🟡 started |
| 7 | Single-user, single-machine design | **Deferred.** Adding auth/multi-tenancy is a product-direction decision (which auth model, is shared deployment even wanted), not something to decide unilaterally. Left as future work. | deferred |
| 8 | Windows-only build tooling | **Done.** Added `build.sh`/`run.sh`/`run_tests.sh` mirroring the `.bat` files file-for-file. Verified: shell syntax checked clean (`sh -n`); the actual `javac`/`java` invocations could only be exercised via git-bash on this Windows machine, which calls the *Windows* JDK (semicolon classpath separator) even though the script correctly uses `:` per the POSIX/Linux JDK convention — so the compile step fails *in this exact test environment* for a reason that has nothing to do with the script being wrong, and it remains genuinely unverified on a real Linux/Mac JDK. Flagged honestly rather than claimed as fully tested. | ✅ (partially verified — see caveat) |
| 9 | Locus Matrix coverage inconsistency (found during validation, not yet a doc'd limitation) | **Root cause found and fixed** — actually two compounding bugs: (a) `MultiLocusScanner.scanOne` assumed the raw GWAS file is sorted by chr:pos (a sequential sweep that never backtracks `nextLocusIdx`), unlike `GwasParser` which explicitly handles unsorted files via interval lookup — rewritten to reuse that same interval-lookup approach; (b) chromosome label mismatches across datasets in the same comparison (e.g. one file using "1", another "chr1") silently failed to match — added `normalizeChr()` used consistently on both sides of the lookup. Both fixes verified with a real regression test (`tests/MultiLocusScannerTest.java`) that fails against the pre-fix code and passes after. | ✅ done, tested |

## New features

### Shared storage for reference data
Options implemented: (a) a plain shared filesystem path (works for any network/lab-shared drive or a
locally-synced cloud folder — Google Drive Desktop, Dropbox, etc. all present as an ordinary local path
once synced, so this one config option covers all of them without needing a real OAuth integration);
(b) a git repository URL, cloned/pulled on demand (covers private GitHub repos with an SSH deploy key or
HTTPS token the user already has configured in their own git credential store — this app shells out to
the user's own `git`, it does not handle credentials itself). Decision: **not** building a bespoke Google
Drive API/OAuth integration — out of scope for a same-day autonomous change and the synced-folder path
covers the real use case with far less risk.

### HPC/Snakemake execution
Opt-in, configured via a new "Cluster (HPC)" settings tab: SSH host, username, path to an existing
private key (never a password field — key-based auth only), remote working directory, a free-form
module-name map (tool → `module load <name>`), and a poll-interval range (5–20 minutes, the app should
pick within this range based on the submitted job's expected size rather than a fixed value). The app
generates a Snakefile wrapping the existing per-locus pipeline (PLINK subset → LD → export), copies it
and the minimal required inputs to the remote working directory via `scp`, submits with
`ssh <host> "cd <remotedir> && snakemake --cores <n> --jobs <m>"` in the background, and polls job
status on the configured interval via `squeue`/`qstat` (whichever the target scheduler uses — both
attempted, first one that returns a sane result wins). All remote commands are invoked via
`ProcessBuilder` with argument arrays (never a concatenated shell string), so no path/hostname the user
types can be interpreted as extra shell syntax.

### GUI setup wizard
Extends the existing home-page "Resources & Settings" panel (already present in `index.html`) with two
new tabs: "Shared storage" and "Cluster (HPC)", both re-openable and re-editable at any time — this is
not a one-shot first-run wizard, it's just more settings tabs, which is simpler, more consistent with
the existing UI, and means "reconfigure later" comes for free instead of needing separate wizard-restart
logic.

## Status: shared storage + HPC feature build — done

A first attempt at delegating this to a background agent reported "in progress" without having
actually written any code — caught by the coordinator independently checking `git log`/working-tree
state rather than trusting the report, and corrected by implementing it directly instead, in five
concrete, compiled-and-tested steps (each committed separately so the history shows real progress,
not just a final summary):

1. `src/rsid/GlobalConfig.java` — added `SharedStorage` and `HpcConfig` nested classes + fields,
   extended `toJson()`/`parse()` to round-trip them (commit `f844afe`). Both default to fully inert
   (`mode="none"`, `enabled=false`) — zero behavior change for existing installs. Verified with a
   round-trip unit test (`tests/GlobalConfigTest.java`).
2. `src/rsid/SharedStorageResolver.java` — resolves a reference-data path against shared storage
   when it doesn't exist directly; `ensureGitClone()`/`test()` shell out to the user's own `git` via
   `ProcessBuilder` argument arrays, never a shell string, never touching credentials. New
   `POST /api/shared-storage/test` endpoint in `LocalServer.java` (commit `e0fd4f4`). Verified with
   a unit test covering every branch that needs no network access.
3. Two new `<h4>` sections ("Shared Storage", "Cluster (HPC)") added to the Resources & Settings
   panel in `index.html`, following the exact existing GET-mutate-POST pattern (commit `c1fc971`).
   **Verified with a real headless-browser test (Playwright)** against an isolated server instance
   on a temporary port (8765 was held by the still-running 30-dataset reprocessing job the whole
   time — confirmed undisturbed before and after; the `PORT` constant was temporarily changed,
   tested, then reverted, with a clean `git diff` confirming the revert before committing). That
   test caught a real bug before it shipped: `saveSharedStorage()` called `loadResources()` without
   awaiting it, so `testSharedStorage()` could write its result into a DOM node that had already
   been silently replaced by the (still in-flight) re-render — fixed by awaiting `loadResources()`
   in both `saveSharedStorage()` and `saveHpc()`.
4. `src/analysis/SnakemakeSubmitter.java` — generates a Snakefile wrapping the existing PLINK
   subset/LD commands with a `module load` prefix per configured tool; builds `scp`/`ssh` argument
   arrays (never shell strings); validates every remote-facing config value against shell
   metacharacters *twice* — once implicitly via `ProcessBuilder`'s argument-array form (blocks
   local shell injection) and once explicitly via `validateNoShellMetachars()` (blocks the value
   from changing what the *remote* shell does, which the local argument-array form can't prevent
   on its own); picks a poll interval that scales with job size within the configured
   `[min,max]` range. **Per the hard safety constraint, `pollStatus()`/an actual submit path were
   never invoked against any real host** — only the pure logic (Snakefile content, interval
   selection, argument-array construction, metachar validation, wrapper-script content) is
   unit-tested (`tests/SnakemakeSubmitterTest.java`, 7 cases, all passing).

**Explicitly not done, by design**: no "Submit job" button anywhere in the UI — HPC config can be
saved and validated, but nothing in this pass actually dispatches work to a cluster. Wiring a real
submit action into the per-project pipeline UI, and testing `pollStatus()`/the scp+ssh submit path
against a real (or realistic loopback) target, are left as explicit future work for when the user
is present to supervise the first real connection.

All 7 test files pass via `run_tests.bat`/`run_tests.sh` as of this update.

## Full 30-dataset, all-loci reprocessing

Launched as a background task (see task log below) reusing the existing scz-full-demo infrastructure
(1000G EUR reference panel, gzip-compressed raw files decompressed on demand). Unlike the manuscript
validation (which processed ~14 illustrative loci per dataset), this run processes **every**
genome-wide-significant locus found by clumping, for **all 30 datasets** (not just the 6 primaries) —
i.e., the 24 datasets previously used only for the lightweight Locus Matrix scan now also get full
PLINK-subset + LD + annotation treatment. Given locus counts already observed (178/130/70/70/70/70 for
primaries alone), this is expected to take many hours; `ld.parallel.jobs` kept conservative (2, not 1
and not the crash-prone 4) with a hard per-locus size cap (skip/flag anything above ~4 Mb) to avoid
repeating the earlier OOM crash while not being fully serial for a job this large.

**Update — the first full run's "success" report was wrong, caught by direct verification, not
by trusting the job's own output.** `full_reprocess.py` finished all 30 datasets and wrote
`full_reprocess_summary.json` with `"status": "ok"` for every one. Spot-checking actual file
counts in `<project>/data/` (not the JSON summary) showed the 6 primary (`*-demo`) datasets were
genuinely correct (real, proportional file counts), but **all 24 secondary (`*-src-*`) datasets had
zero output files** despite being marked done. Root cause, confirmed via
`curl http://localhost:8765/api/project/<id>/progress` directly against the still-running server:
every one of the 24 secondary configs has a leftover `col.varid=varid` default that doesn't match
that dataset's real header (each has a different actual column set), which the pipeline's
`Config.validate()` treats as a hard error for a *configured-but-absent* optional column — so
every secondary dataset's processing failed immediately. Separately, `full_reprocess.py`'s
polling loop only checked `done`, not whether the reported `phase` was actually an error string,
so it recorded every one of these failures as success. Both are now understood and a fix +
re-run for the 24 affected datasets is in progress (see task log). Lesson applied going forward
in this file: a job's own "done"/"ok" self-report is not verification — checking the actual
artifact (file count, in this case) is.

**Second update — the fix-and-rerun batch also crashed partway, root cause different and smaller.**
After the `col.varid` fix, re-verified live against the running server that `cad-src-aragam2022`
genuinely processes correctly (watched real progress, not just a status flag). The corrected batch
then silently stopped advancing for 30+ minutes with java's memory completely static (not GC noise —
byte-identical across checks), which looked like it could be a repeat of the original OOM failure
mode. Direct investigation (querying `/api/project/<id>/progress` live, checking a stdout capture
file for a Python traceback) found the real cause: the Java pipeline had actually finished
`cad-src-nam2026` correctly (`done:true`, real files confirmed on disk), but the orchestrator's
cleanup step (`os.remove()` on the decompressed raw file, to re-gzip it) hit
`PermissionError: WinError 32` because Windows still held a lock on the file — unlike POSIX, an
open handle blocks deletion outright. That exception was unhandled and killed the whole batch,
stopping progress on the 8 datasets still queued behind it, even though the actual GWAS pipeline
work was fine. Fixed with a retry-with-backoff-then-warn-and-continue wrapper around that one
`os.remove()` call (never let a housekeeping/cleanup step abort real, already-succeeded work), added
a missing done-marker for `cad-src-aragam2022` (it was verified correct manually outside the normal
per-dataset flow, so never got one, and would otherwise have been wastefully reprocessed from
scratch), and restarted the orchestrator — it correctly skipped everything already marked done and
resumed exactly where it left off.

**Final result — genuinely complete, independently verified.** All 30 datasets now have real,
non-zero, loci-proportional output confirmed by direct `ls <project>/data/` file counts (not API
status, not the orchestrator's own summary JSON alone) across every one of the 6 disease groups:

| Dataset | Files | Dataset | Files | Dataset | Files |
|---|---|---|---|---|---|
| scz-demo | 358 | cad-demo | 154 | ra-demo | 134 |
| scz-src-clozuk2018 | 284 | cad-src-aragam2022 | 422 | ra-src-eyre2012 | 30 |
| scz-src-pgc1-2011 | 10 | cad-src-hartiala2021 | 166 | ra-src-glanville2021 | 24 |
| scz-src-pgc2-2014 | 192 | cad-src-nam2026 | 148 | ra-src-shigesi2025 | 8 |
| scz-src-sweden2013 | 24 | cad-src-nikpay2015 | 88 | ra-src-stahl2010 | 14 |
| t2d-demo | 274 | alz-demo | 34 | ibd-demo | 186 |
| t2d-src-diagramv3 | 20 | alz-src-belloy2024 | 4 | ibd-src-franke2010 | 118 |
| t2d-src-mahajan2018 | 120 | alz-src-eadb2026 | 64 | ibd-src-liu2015 | 280 |
| t2d-src-saxena2007 | 4 | alz-src-morenograu2019 | 6 | ibd-src-verma2024 | 10 |
| t2d-src-wood2016 | 4 | alz-src-nicolas2025 | 30 | ibd-src-zorina2023 | 18 |

**Total: 3,228 output files across all 30 datasets.** Datasets with very few/zero loci
(`t2d-src-saxena2007`, `t2d-src-wood2016`, `alz-src-belloy2024`) are genuine results, not failures —
these are small/underpowered/narrowly-scoped studies (an early 2007 T2D GWAS, an insulin-secretion-
focused sub-study, and an X-chromosome-only AD scan respectively), correctly finding few or no
genome-wide-significant autosomal loci. Java server remained stable throughout (memory fluctuated
5–8.6GB across the run but never crashed after the two fixes above), no stray processes left
running unexpectedly.

## What this file will NOT do

- Will not push anything to GitHub without a repo URL (still not provided).
- Will not attempt any live HPC/cluster connection.
- Will not touch the real research project data under LYNXgwas's own `projects/` folder.
- Will not update the manuscript's claims about a fix until that fix has actually been rebuilt and
  verified working, not just written.
