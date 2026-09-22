# Phase 6 decisions log — GitHub project save/extract with safe token handling

Same contract as prior phases: written and updated incrementally and honestly. `[ ]` not done,
`[x]` done and verified, blockers stated plainly.

## 0. Scope of the request

A project (not just shared reference data, which already has git support) should be pushable to and
pullable from a GitHub repository, including a private one, with a GitHub icon in the UI and an
explicit, carefully designed way to supply a token: an "ask every time, nothing stored" mode and a
persistent "stored the extremely safe way" mode, with the user choosing between them.

## 1. What already exists (checked first, not assumed)

`SharedStorageResolver`/`GlobalConfig.SharedStorage` already git-clones/pulls a *reference-data*
repository (PLINK panel, GFF3), shelling out to the user's own `git` binary via an argument-array
`ProcessBuilder` and never touching a credential directly — that pattern is reused, not reinvented,
for the new project-level sync. This machine's git is already configured with
`credential.helper=manager` (Git Credential Manager), which on Windows stores credentials in Windows
Credential Manager, DPAPI-encrypted and tied to the Windows login — a real, audited, OS-level secure
store already sitting there unused by the reference-data feature (which relies on ambient credentials
already being configured, but never explicitly sets one up).

## 2. Security design (resolved here, the part that actually matters)

- **LYNXgwas's own code and config files never contain a token, in either mode.** For the persistent
  mode, the token is handed to `git credential approve` over **stdin** (never a CLI argument, never
  written to a LYNXgwas config file) — this delegates real storage entirely to the already-configured
  OS credential store. LYNXgwas can ask "is a credential already resolvable for github.com?" via
  `git credential fill`, but that check discards the returned password immediately and only ever
  reports a boolean back to the UI — never the value.
- **For the one-time (unsaved) mode**, the token must reach `git` without ever touching disk as a
  persistent artifact and without appearing as a plain command-line argument (visible to other
  processes on the same machine via a process listing). The mechanism: a small, fixed, *secret-free*
  helper script (safe to leave on disk permanently, contains no token, ships as a template) is pointed
  to via `GIT_ASKPASS`; the actual token is passed only through that one child process's environment
  block (`ProcessBuilder.environment()`), which is not written to any file and disappears when the
  process exits. This is materially safer than the two common alternatives it deliberately avoids:
  embedding the token in the remote URL (persists in `.git/config`, readable by `git remote -v`
  indefinitely after) or passing it as a bare CLI argument (visible in a process listing for the
  process's lifetime).
- **Honest limitation, stated up front**: an environment variable on a child process is not a perfect
  secret on a multi-user or compromised machine (a sufficiently privileged local process/administrator
  can still read another process's environment block) — this is a real, known limitation shared by
  essentially every CLI credential-injection pattern, not something specific to this implementation,
  and is documented rather than glossed over. The stored-credential mode's actual security rests on
  Windows Credential Manager, not on anything LYNXgwas does itself.
- **Scope-limited to github.com over HTTPS**: the repo URL is validated as `https://github.com/...`
  before anything runs. This is a real safety boundary (predictable credential host, no accidental
  operation against an arbitrary/malicious git remote), not just a convenience.
- **No destructive default behavior**: pull/"extract" refuses to clone into a non-empty directory that
  isn't already that same repo's own clone (mirroring the precaution this session has already applied
  to other clone operations) rather than overwriting anything silently.

## 3. Work items and status

### 3.1 `GithubTokenStore` — credential save/forget/status via git's own store [x]
`git credential approve`/`reject`/`fill`, all via stdin, none of it ever logged or returned to a
caller as plaintext beyond the single in-process call that discards it immediately. Verified by
`tests/GithubTokenStoreTest.java` (8/8 assertions) against a real, disposable git "store" credential
helper — fail-before-fix confirmed on the "store file contains expected username" assertion.

### 3.2 `GithubProjectSync` — push/pull a whole project directory [x]
Reuses the `SharedStorageResolver` ProcessBuilder/timeout conventions. Push: init-if-needed, add,
commit (skip cleanly if nothing changed), push. Pull ("extract"): clone if empty/missing, pull if
already the same repo's clone, refuse otherwise. Verified by `tests/GithubProjectSyncTest.java`
(16/16 assertions) against a real local `git init --bare` repo, plus a real end-to-end run against a
throwaway private GitHub repo (see 3.6).

### 3.3 One-time token injection via `GIT_ASKPASS` + child-process environment [x]
Static, secret-free askpass helper script; token only ever lives in one child process's environment
block for the duration of one git invocation. Verified against the real throwaway GitHub repo: pushed
successfully using a token extracted via `git credential fill` and passed as a one-time token; direct
inspection of `.git/config` afterward showed only the plain HTTPS URL, no embedded credential; no
leftover `lynxgwas_askpass_*.cmd` temp file remained (cleanup in the `finally` block).

### 3.4 New `LocalServer` endpoints [x]
Implemented as: `GET /api/github/token-status` (boolean only), `POST /api/github/token` (save) /
`DELETE /api/github/token` (forget), `POST /api/project/{id}/github/push` (project-scoped push),
`POST /api/github/import` (top-level; takes a target `project_id` and clones/pulls into
`projects/<project_id>`). Names differ slightly from the original sketch (`import` instead of
`pull-into`, and it's top-level rather than project-scoped since the project directory may not exist
yet at import time) — noted here rather than silently diverging from the plan.

### 3.5 UI: GitHub icon + settings [x]
An inline-SVG GitHub octocat icon on each project row's "Save" action, and an "Import from GitHub"
button (with the same icon) in the home-page toolbar; a "GitHub (save / extract projects)" settings
section with the two credential modes described in section 2 (radio choice, "Ask me each time"
selected by default — nothing assumed/saved silently). Verified live via Playwright against the
running production server: settings section renders with the correct default selection, the import
modal opens and its fields render, and a project row's GitHub "Save" button opens the push modal with
its fields — screenshots captured for all three
(`github_settings_panel.png`, `github_import_modal.png`, `github_push_modal.png`). A stray leftover
lock-emoji reference in the settings copy (from an earlier draft, before the icon was finalized as the
GitHub octocat) was caught during this check and corrected to describe the actual button.

### 3.6 Testing [x]
Unit tests against a local `git init --bare` repository (no real network dependency, matching this
project's existing convention of not depending on a real external service in the automated suite) for
push/pull correctness, plus a real end-to-end check against one real, throwaway private GitHub repo
(`AlsammanAlsamman/lynxgwas-test-project-sync`) covering: push via ambient credential, clone-extract
with real content verification, a second push+pull round-trip, and the one-time-token path. A
dedicated check confirmed the token never appears in `.git/config` in either the pushed project or the
extracted clone. **Known follow-up**: the throwaway test repo could not be deleted programmatically
(`DELETE /repos/.../lynxgwas-test-project-sync` → HTTP 403, the ambient credential lacks the
`delete_repo` scope) — it still exists on the user's GitHub account (private) and needs manual
deletion.

## 4. Hard boundaries (carried over, unchanged)

No push without explicit instruction beyond what this feature itself is for; the real production
server is never restarted without asking first; nothing here changes existing shared-storage behavior.
