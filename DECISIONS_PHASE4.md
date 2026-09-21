# Phase 4 decisions log — local/cloud LLM Agent panel

Same contract as `DECISIONS.md`/`DECISIONS_PHASE2.md`/`DECISIONS_PHASE3.md`: written and updated
incrementally and honestly. `[ ]` not done, `[x]` done and verified, blockers stated plainly.

## 0. Scope of the request

A robot-icon button on the home page opens an interactive agent panel. The agent can talk to the
user, read the project/tool configuration, create/edit projects, and run the pipeline/analysis tools
on the user's behalf — making the tool "AI-friendly": an LLM should be able to drive LYNXgwas through
its own YAML/JSON configuration surface instead of the human clicking through the wizard by hand. The
user can point the agent at either (a) a locally-run open-source LLM on their own machine, or (b) a
cloud API of their choice, by supplying a key.

## 1. Architecture decisions (resolved here, not left implicit)

- **Wire protocol: OpenAI-compatible Chat Completions with tool-calling** (`POST {base_url}/chat/
  completions`, `tools`/`tool_calls` schema). One client implementation this way transparently covers
  Ollama (its own `/v1` OpenAI-compatible endpoint), LM Studio (`/v1`), llama.cpp's `server` binary
  (`/v1`), and cloud providers using the same shape (OpenAI itself, Groq, Together, Fireworks, etc.) —
  both the "local open-source model" and "bring your own API key" cases the request asks for, from a
  single code path. Anthropic's native Messages API uses a different request/response shape; out of
  scope for this MVP, left as a documented extension point behind the same `LlmClient` interface
  rather than hard-coded around.
- **Two configured modes**, both ultimately just populate `{base_url, api_key, model}`:
  - *Local*: default `base_url=http://localhost:11434/v1` (Ollama's own default port), `api_key`
    blank. `GET /api/agent/detect-local` probes `localhost:11434/api/tags` (Ollama) and
    `localhost:1234/v1/models` (LM Studio) so the UI can offer a live dropdown of models the user
    already has pulled, instead of asking them to type a model name blind.
  - *Cloud API*: user supplies `base_url` (default `https://api.openai.com/v1`), `api_key`, and
    `model` by hand.
  Persisted to `<workspace>/.agent_config.json` (gitignored); the API key is never logged or echoed
  back in any response body other than the config-read endpoint itself.
- **Capability = a fixed tool registry, each tool a thin loop-back HTTP call to this same server's
  existing REST API** — `list_projects`, `get_project_config`, `create_or_update_project`,
  `run_pipeline`, `get_pipeline_progress`, `list_analysis_tools`, `run_analysis_tool`,
  `get_project_manifest`. The agent can therefore never do anything a human user couldn't already do
  by clicking through the existing UI — no new capability or attack surface is introduced beyond what
  already exists; it's automation of the existing surface, not a new one. Deliberately **no
  delete-project tool** in this MVP — destructive actions stay a manual, explicit UI action outside
  the agent's reach.
- **Transparency**: every tool call the agent makes is returned to the frontend as a visible "step"
  (name, arguments, one-line result summary) in call order — the user sees exactly what the agent did,
  never a hidden action. The tool-call loop is capped at 6 rounds per user message to bound runaway
  looping.
- **UI**: a robot-icon button in the home page header opens a bottom slide-up panel — chat/step log,
  a text input, and a settings gear for the local/cloud model configuration described above. One
  blocking request per user turn (`POST /api/agent/chat` returns the full set of steps plus the final
  reply together) rather than token-level streaming — a stated MVP simplification, not an oversight.

## 2. Work items and status

### 2.1 `MiniJson` — minimal hand-rolled JSON parser [x] — done, verified
`src/analysis/MiniJson.java`: recursive-descent parser (objects/arrays/strings incl. `\u` escapes/
numbers/booleans/null) plus a symmetric `encode()` serializer reused by `LlmClient`,
`AgentToolRegistry`, and `LocalServer.agentChat` instead of duplicating a third hand-rolled encoder.
`tests/MiniJsonTest.java`: 24 cases including a real-shaped OpenAI tool-call response (nested
objects/arrays, a JSON-encoded string as a value, a Windows-style backslash path). Fail-before-fix
verified (broke the `\\` escape case, both re-parse assertions failed; reverted, passed again).

### 2.2 `LlmClient` — OpenAI-compatible chat-completions HTTP client [x] — done, verified
`src/analysis/LlmClient.java`: `java.net.http.HttpClient`-based; takes `{base_url, api_key, model}`
plus a message list and tool schema, returns the assistant's text and/or parsed tool calls. Exercised
for real (not just compiled) in the section 2.7 end-to-end run below.

### 2.3 `AgentToolRegistry` — tool schemas + dispatch [x] — done, verified
`src/analysis/AgentToolRegistry.java`: 8 tool schemas, dispatched as real HTTP calls back to
`http://localhost:<LocalServer.PORT>/...` — reusing the exact same code path a human click would hit,
not reimplementing project-creation/pipeline-running logic. Deliberately handles the one real
correctness trap in this design: `POST /api/project/{id}/config` applies its body onto a *fresh*
`Config` (missing fields silently revert to hard-coded defaults), so `create_or_update_project` reads
the current config first and merges the model's fields onto it before writing back — see the class's
own javadoc. `tests/AgentToolRegistryTest.java` (28 cases) against a small fake HTTP backend on an
ephemeral port (never collides with a real server): every tool, plus the merge behavior specifically —
fail-before-fix verified (removed the merge step, the two "untouched field preserved" assertions
failed as expected; restored, passed again).

### 2.4 `AgentOrchestrator` — the tool-calling loop [x] — done, verified
`src/analysis/AgentOrchestrator.java`: message-history + tool-calling loop, capped at 6 rounds, records
every tool call as a `Step`. Exercised for real in the section 2.7 end-to-end run.

### 2.5 New `LocalServer` endpoints [x] — done, verified
`POST /api/agent/chat`, `GET`/`POST /api/agent/config`, `GET /api/agent/detect-local` (probes
`localhost:11434/api/tags` for Ollama and `localhost:1234/v1/models` for LM Studio with a ~1s timeout).
Compiles clean with the rest of `LocalServer.java`; exercised for real below.

### 2.6 Frontend: robot-icon agent panel in `index.html` [x] — done, verified
Bottom slide-up panel: settings (local/cloud radio, detected-local-model dropdown or cloud
base_url/api_key/model fields), a chat/step log (user bubbles, assistant bubbles, a distinct line per
tool call with a ✓/✗ marker and its one-line summary), and an input box. Confirmed rendering and
working via the Playwright pass below, including live in the existing project table once a project is
created.

### 2.7 Testing [x] — done, verified
- `MiniJsonTest` (2.1) and `AgentToolRegistryTest` (2.3): both wired into `build.bat`/`build.sh` and
  `run_tests.bat`/`run_tests.sh`; a full clean-room `rm -rf bin && build.bat` + `run_tests.bat` run
  passed all suites, old and new, together.
- **Real end-to-end run** (not just unit tests): built a tiny scripted stand-in HTTP server
  (`fake_llm.py`) implementing the OpenAI chat-completions response shape with a hand-scripted,
  deterministic tool-call sequence (`create_or_update_project` → `run_pipeline` → final text), and a
  second, fully isolated copy of the real LYNXgwas server (patched `PORT=8777`, scratch workspace, a
  minimal placeholder project so `Main` would start) — the real production server (PID 32928, port
  8765) and the article-demo server (PID on port 8766) were both confirmed still healthy before and
  after and were never touched. Pointed the agent config at the fake LLM and sent one real chat message
  through `/api/agent/chat`; verified:
  - the exact expected two tool calls executed in order, both `ok:true`;
  - a real `projects/agent-e2e-demo/config.properties` appeared on disk with exactly the 5 fields the
    (scripted) model supplied and `Config`'s own defaults for everything else;
  - `run_pipeline` genuinely started a real background job (`GET .../progress` showed a real, honest
    error once it hit the test's intentionally-nonexistent `loci.txt`, i.e. real execution, not a stub).
- **Playwright UI pass** against that same live test server: robot icon visible and clickable, panel
  opens, settings panel toggles, a sent message renders as a user bubble, both tool-call steps render
  with their summaries, the final assistant reply renders, and the newly-created project is visible
  live in the existing home-page project table behind the panel — screenshot inspected directly.
- All test-only processes (fake LLM, test server) stopped and the scratch workspace deleted afterward.
- **Honest limitation, stated up front, not glossed over**: no real local open-source LLM runtime
  (e.g. an actual installed Ollama with a downloaded model) is available in this sandboxed dev
  environment, and this phase does not download one. The scripted stand-in server exercises the
  identical wire protocol a real Ollama/LM Studio/OpenAI endpoint would use, so the *mechanism*
  (tool-calling loop, real side-effects, UI wiring, the config merge correctness) is genuinely verified
  end-to-end; real model behavior/prompt quality/tool-selection judgment against an actual model is not
  verified in this environment.
- **Important deployment caveat, discovered and reported rather than silently left**: `index.html`/
  `viewer.html` are static files the live production server re-reads from disk on every request, so
  the new robot-icon UI is already visible on the real running production server (port 8765) the
  moment these files were edited, with no restart needed for that part. The three new `/api/agent/*`
  *endpoints*, however, live in `LocalServer.class`, which the JVM only loads once at process start —
  the already-running production server process was started before this phase's code existed, so it
  is still running the old compiled class and returns 404 for `/api/agent/config` right now. The UI is
  live but non-functional on the real server until it is restarted; per this project's own established
  rule this session, the live production server is never restarted without asking the user first.

## 3. Hard boundaries (carried over, unchanged)

No push without explicit instruction; no writes to real `projects/*` beyond what the agent's own
tools legitimately do through the existing API surface; every tool call logged, never hidden; the real
production server is never restarted without asking first.

## 4. Follow-up: making local-model setup easier (same session)

Explicit user request: even installing/finding a local model was too much friction — asked for a
public, stable link and an easier way to get a model running.

### 4.1 `OllamaPuller` — one-click model download via Ollama's own `/api/pull` [x] — done, verified
`src/analysis/OllamaPuller.java`: drives Ollama's real, documented, streaming `POST /api/pull`
endpoint in a background thread, tracked via a simple `PullStatus` (status text, completed/total
bytes, done, error) polled by the UI — this cannot and does not install Ollama itself (a real system
installer, intentionally left to the user via the official link, never silently automated).
`tests/OllamaPullerTest.java` (8 cases) against a fake streaming NDJSON server: incremental progress,
an Ollama-style in-stream `{"error":...}` case (bad model name — reported as HTTP 200, not a failed
status code), and a connection-refused case (nothing listening) — all reach `done=true` without
hanging. Fail-before-fix verified (broke `completed`-byte tracking; the "final completed == total"
assertion failed as expected; reverted, passed again).

### 4.2 New `LocalServer` endpoints [x] — done, verified
`POST /api/agent/pull-model` (body `{"model": "..."}`, proxies to `http://localhost:11434/api/pull`),
`GET /api/agent/pull-progress`. Compiles clean; exercised for real below.

### 4.3 Settings-panel UI [x] — done, verified
When Ollama isn't detected: a direct link to `https://ollama.com/download` (the real, official,
stable installer page — Windows/Mac/Linux) with a one-line explanation. When it is detected: a short,
curated list (Llama 3.2 3B, Qwen 2.5 7B — both documented by Ollama as tool-calling-capable, not the
full model catalog, so a first-time user isn't asked to choose among hundreds) each with a "Download"
button that starts the pull, polls progress live (MB downloaded / total, percent), and on success
re-detects and auto-selects the new model — no terminal required anywhere in this flow.

### 4.4 Real-world verification — an unplanned but valuable discovery [x]
While testing, `netstat`/`tasklist` revealed a **real, already-installed Ollama** (`ollama.exe`,
version 0.34.2, PID confirmed) genuinely running on this machine with a real model already pulled
(`gemma4:26b`, ~18.6 GB, tool-calling-capable) — something a much earlier check in this same session
had confirmed was *not* present, so it must have been installed at some point in between. This
upgraded the verification from "protocol-correct against a scripted stand-in" (section 2.7) to a
**genuine live test against a real local model**:
- `GET /api/agent/detect-local` against the real Ollama correctly returned `available:true` and the
  real model name.
- Configured the real production agent (`POST /api/agent/config`) to use it, then sent a real,
  deliberately read-only question ("how many projects do I have?") through `/api/agent/chat` on the
  live production server. The real model correctly called the real `list_projects` tool and received
  real data back (49 real projects) — genuine end-to-end tool-calling against a real model, not a
  mock. **Reported honestly, not glossed over**: the model then also called `list_analysis_tools`
  unprompted and its final answer addressed the wrong question (it described analysis tools instead
  of stating the project count) — a real model-judgment limitation for this specific small model on
  this specific prompt, not a bug in the tool-calling mechanism itself, which behaved correctly
  throughout. This is exactly the kind of finding worth surfacing rather than hiding behind a vague
  "it works."
- A real bogus-model-name pull was (accidentally, harmlessly) exercised against the real Ollama via
  the new endpoint and correctly surfaced Ollama's own real error ("pull model manifest: file does not
  exist") through to `GET /api/agent/pull-progress` — confirming the error path works against the real
  API, not just the fake one.
- The production server was restarted twice this session (each time only after asking first and
  getting an explicit go-ahead) to pick up first the core agent backend, then this follow-up; both
  restarts verified via clean-room rebuild + full test suite (`OllamaPullerTest` included) passing
  beforehand, and a live endpoint check afterward. It is currently live and configured to use the
  user's real local model.

## 5. Closing status

All work items across both the original scope (section 0) and the same-session easier-setup follow-up
(section 4) are done and independently verified: a real, working tool-calling agent, driven by either
a locally-run open-source model or a user-supplied cloud API key through one shared OpenAI-compatible
client, with a fixed, transparent, loop-back-only tool surface (list/read/create/edit projects, run
the pipeline, run an analysis tool, check progress), a one-click local-model download via Ollama's own
API, and a working chat/step panel in the existing home page. A clean-room rebuild and the full test
suite (34 pre-existing tests plus 24 `MiniJsonTest`, 27 `AgentToolRegistryTest`, and 8
`OllamaPullerTest` new cases) pass together. Verification went further than originally planned: what
started as a scripted-stand-in-only end-to-end run (section 2.7) was upgraded mid-phase by a genuine
discovery — a real Ollama with a real model already running on this machine (section 4.4) — into an
actual real-model test, including one real, honestly-reported model-judgment miss (see 4.4), not just
protocol-level plumbing.

**Activated on the real production server, twice, each time only after asking first**: the server was
restarted once for the core agent backend and once more for the easier-setup follow-up, since both
live in compiled Java code the running JVM does not hot-reload (unlike the static HTML/JS, which took
effect immediately). Both restarts were preceded by a clean-room rebuild + full passing test suite and
followed by a live endpoint check. The production server is currently live, has the complete feature
set, and is configured to use the user's own already-installed local model.

