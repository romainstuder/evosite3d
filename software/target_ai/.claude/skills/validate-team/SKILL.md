---
name: validate-team
description: Multi-agent target validation. Fans out specialist subagents (clinical, druggability, biology, landscape) in parallel next to the deterministic scoring script, then reconciles them. Use when the user asks for a team/multi-agent/second-opinion validation, or a deeper validation than /validate across several targets and diseases.
---

# Multi-Agent Target Validation

```
/validate-team BRAF,KRAS — cutaneous melanoma,non-small cell lung carcinoma
```

Same input format as `/validate`: `TARGETS — DISEASES`, separated by `—`, `--` or `vs`.

## Architecture

```
            ┌──────────────┐
  input ──▶ │ ot-resolver  │  (haiku)  IDs + ambiguity check
            └──────┬───────┘
                   │ confirmed IDs
   ┌───────────────┼───────────────┬───────────────┬────────────────┐
   ▼               ▼               ▼               ▼                ▼
 script        clinical-      druggability-    biology-       landscape-
 (baseline)    evidence       evidence         evidence       scout (web)
   │               │               │               │                │
   └─ CSV/HTML ────┴──── dossiers in results/agents/{run_id}/ ──────┘
                   │
            ┌──────▼───────────────┐
            │ evidence-synthesizer │  → consensus.md
            └──────────────────────┘
```

Subagents can't spawn subagents, so **you (the main session) are the orchestrator**.
The script stays the reproducible baseline. The agents give an independent,
evidence-cited second opinion.

## Instructions

Given `$ARGUMENTS`:

1. **Parse** the target and disease lists. Set
   `run_id = {YYYYMMDD-HHMM}_{targets joined by _}` (lowercase).

2. **Resolve.** Spawn `ot-resolver` with both lists. If any item comes back
   `ambiguous`, `broad` or `not_found`, show the options and ask the user before going
   on. Don't guess.

3. **Fan out, all in a single message so they run concurrently:**
   - Bash: `python3 open_targets_client.py validate "ENSG…,ENSG…" "EFO_…,EFO_…"`
     (pass the **resolved IDs** so the script and the agents score the same entities)
   - Spawn `clinical-evidence`, `druggability-evidence`, `biology-evidence` and
     `landscape-scout`. Give each the same brief:
     ```
     run_id: <run_id>
     pairs:
       - BRAF (ENSG00000157764) × cutaneous melanoma (EFO_0000389)
       - ...
     write your dossier to: results/agents/<run_id>/<agent-name>.json
     ```

4. **Find the baseline outputs.** The script prints `Open results/validation_{slug}.html`.
   The CSV is `results/scores_{slug}.csv`. The slug comes from resolved names in
   composite order, so read it from the output rather than predicting it.

5. **Synthesize.** Once all four agents have reported (a failed agent is fine; note it),
   spawn `evidence-synthesizer` with the `run_id`, the CSV path, and the dossier
   directory.

6. **Present:**
   - The synthesizer's ranked table and headline findings, as is
   - Discrepancies of 2 points or more, if any. These matter most to the user.
   - Point the user to `results/validation_{slug}.html` (interactive matrix) and
     `results/agents/{run_id}/consensus.md` (multi-agent report)

## Rules

- Every rule in `CLAUDE.md` still applies: no curl, no `python3 -c`, single-line
  script calls.
- Don't do evidence collection in the main session. That's what the agents are for,
  and doing it here defeats the context isolation.
- For a single target, `/score-target` is cheaper. Suggest it when the input is
  one target with one disease and the user hasn't asked for the team.
