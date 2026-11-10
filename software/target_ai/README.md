# 🧬 Target Validation Tool for Claude Code

An AI-powered drug target validation assistant that scores and ranks therapeutic targets against diseases using real-time data from Open Targets Platform and supplementary biomedical databases.

## What It Does

Given one or more **gene targets** and **diseases**, this tool:

1. **Resolves** gene symbols → Ensembl IDs and disease names → EFO/MONDO IDs
2. **Fetches** evidence from Open Targets Platform (genetic associations, clinical trials, tractability, pathways, safety)
3. **Scores** each target–disease pair on three dimensions (0–5 each):
   - **Clinical Evidence** — genetic and somatic evidence, plus trials of this target's drugs _in this disease_
   - **Druggability** — small-molecule/antibody tractability, approved and clinical-stage drugs
   - **Pathway/Biology** — pathway membership, tissue expression, animal models
4. **Ranks** by weighted composite score with compact narrative
5. **Generates** an interactive HTML scoring matrix plus Markdown report, CSV, and raw JSON

An optional **multi-agent mode** (`/validate-team`) runs specialist subagents next to the
scoring script to give an independent, evidence-cited second opinion. See
[Multi-Agent Mode](#multi-agent-mode).

## Quick Start

### Prerequisites

- [Claude Code](https://docs.anthropic.com/en/docs/claude-code/overview) installed (`npm install -g @anthropic-ai/claude-code`)
- An Anthropic API key or Claude Pro/Max subscription
- Python 3.10+ (for the helper script; Claude Code can also query APIs directly)
- No additional API keys needed — Open Targets Platform is free and open

### Setup

```bash
# Clone or download this project
git clone https://github.com/romainstuder/target-ai target-validation-tool
cd target-validation-tool

# Start Claude Code
claude
```

That's it. Claude Code reads the `CLAUDE.md` file on startup and becomes your target validation assistant.

### Usage — Commands

Once inside Claude Code, use these commands:

#### `/validate` — Multi-target × multi-disease scoring

Runs the deterministic scoring script (`validate.md`).

```
/validate BRAF,KRAS — cutaneous melanoma,non-small cell lung carcinoma
/validate PCSK9,HMGCR,ANGPTL3 — coronary artery disease,familial hypercholesterolemia
/validate CDK4,CDK6 — breast carcinoma,ovarian carcinoma
/validate JAK1,JAK2,JAK3,TYK2 — rheumatoid arthritis
```

Prefer specific disease terms ("cutaneous melanoma", not "cancer"): Open Targets scores
each ontology level separately, so umbrella terms miss evidence filed under their
children.

#### `/score-target` — Deep dive on one target

SWOT analysis, sub-scores, and supplementary web research for a single target (`SKILL.md`).

```
/score-target PCSK9 — coronary artery disease, familial hypercholesterolemia
```

#### `/validate-team` — Multi-agent validation

The script plus parallel specialist subagents, reconciled into a consensus report. See
[Multi-Agent Mode](#multi-agent-mode).

```
/validate-team BRAF,KRAS — cutaneous melanoma,non-small cell lung carcinoma
```

### Usage — Natural Language

You can also just talk to Claude naturally:

```
> What are the best drug targets for type 2 diabetes? Score the top 10.

> Is TREM2 a good target for Alzheimer's? Compare it against BACE1 and BIN1.

> I'm working on a grant proposal for targeting GPR84 in inflammatory bowel disease.
  Can you validate this target and identify the key risks?

> Score these oncology targets for pancreatic cancer: KRAS, TP53, CDKN2A, SMAD4
```

### Usage — Data access

Two ways to reach Open Targets data, no API key needed either way:

**MCP server** (`open-targets-platform`, registered in `.mcp.json`) — for search, ID
resolution, and any one-off lookup. Claude calls its tools directly during a
conversation: `search_entities`, `get_open_targets_graphql_schema`,
`get_type_dependencies`, `query_open_targets_graphql`, `batch_query_open_targets_graphql`.
No script needed for these.

**Python script** (`open_targets_client.py`) — for the two commands that must produce
deterministic, reproducible, file-based output:

```bash
# Full validation — scores every target × disease pair and writes results/
python3 open_targets_client.py validate "BRAF,KRAS" "cutaneous melanoma,non-small cell lung carcinoma"

# Pathway visualization — interactive Reactome map
python3 open_targets_client.py pathways "PCSK9,HMGCR,ANGPTL3"
```

Targets and diseases can be given as names or as IDs (`ENSG…`, `EFO_…`, `MONDO_…`).
Names resolve to the top Open Targets search hit, so pass IDs when a name is ambiguous.

The script stops with an `OpenTargetsError` if any Open Targets query fails (HTTP or
GraphQL error). It never scores a failed query as "no evidence".

## Multi-Agent Mode

`/validate-team` runs the scoring script and a team of
[Claude Code subagents](https://docs.anthropic.com/en/docs/claude-code/sub-agents) in
parallel. Your main Claude Code session acts as the orchestrator, because subagents
can't start other subagents.

```
            ┌──────────────┐
  input ──▶ │ ot-resolver  │  IDs + ambiguity check (asks you if a term is broad)
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

| Agent                   | Role                                                                                | Tools               |
| ----------------------- | ----------------------------------------------------------------------------------- | ------------------- |
| `ot-resolver`           | Resolves names to IDs; flags ambiguous or umbrella disease terms (runs on Haiku)    | Open Targets search |
| `clinical-evidence`     | Scores Clinical: genetics, trials, approvals for this indication                    | Open Targets MCP    |
| `druggability-evidence` | Scores Druggability: tractability, modalities, chemical matter, safety risks        | Open Targets MCP    |
| `biology-evidence`      | Scores Pathway/Biology: disease-relevant pathways, expression, mouse models         | Open Targets MCP    |
| `landscape-scout`       | Recent trials, competitors, patents, publications (does not score)                  | Web search          |
| `evidence-synthesizer`  | Compares agent and script scores, flags gaps of 2 or more points, writes the report | Reads files only    |

- **The script's scores stay the official ones.** The agents never edit the CSV or HTML.
  They give a second opinion in which every claim cites an Open Targets field or a URL.
- **Each agent covers one dimension for all pairs**, using batched queries, so the
  number of agents doesn't grow with the size of the input.
- **Each agent writes a JSON dossier** following `dossier_schema.md`. The synthesizer
  merges them into `results/agents/{run_id}/consensus.md`, with an agent composite
  score next to the script's composite for each pair.
- **A run costs several times more than `/validate`**: six agent runs instead of one
  script run. Use `/score-target` for a single target and disease.

Agent definitions live in `.claude/agents/`, and the orchestration steps in
`.claude/skills/validate-team/SKILL.md`. Start `claude` from this directory so they're
picked up.

## Scoring Methodology

### Dimensions & Weights

| Dimension         | Weight | Data Sources                                                           |
| ----------------- | ------ | ---------------------------------------------------------------------- |
| Clinical Evidence | 40%    | GWAS Catalog, ClinVar, ClinGen, UniProt, clinical trials (ChEMBL)      |
| Druggability      | 35%    | Open Targets tractability, ChEMBL drugs, structural data               |
| Pathway/Biology   | 25%    | Reactome, literature mining (EuropePMC), expression (Expression Atlas) |

### Composite Score

```
Composite = (Clinical × 0.40) + (Druggability × 0.35) + (Pathway × 0.25)
```

### How Drugs Are Counted

- **Druggability** counts any approved or clinical-stage drug acting on the target, in
  any indication.
- **Clinical** counts only this target's drugs that have been trialled **in this
  disease**: an approval for NSCLC doesn't raise a melanoma score. Matching uses the
  exact disease term, so an approval filed under a broader parent term (e.g.
  "melanoma" for "cutaneous melanoma") isn't counted.
- Drug stages come from Open Targets' `maxClinicalStage` (`APPROVAL`, `PHASE_4`,
  `PHASE_3`, `PHASE_2_3`, …). Withdrawn drugs count as Phase 3.

### Confidence Levels

| Level      | Criteria                                                |
| ---------- | ------------------------------------------------------- |
| **High**   | Composite ≥ 3.5 and all three dimensions scored above 0 |
| **Medium** | Composite 2.0–3.49, or ≥ 3.5 with a dimension scored 0  |
| **Low**    | Composite < 2.0                                         |

## Output Files

After each validation run, results are saved to the `results/` directory. File names encode both targets and diseases:

| File                                      | Contents                                                |
| ----------------------------------------- | ------------------------------------------------------- |
| `validation_{targets}_vs_{diseases}.html` | **★ Interactive HTML scoring matrix** — open in browser |
| `validation_{targets}_vs_{diseases}.md`   | Markdown report with tables and narrative               |
| `scores_{targets}_vs_{diseases}.csv`      | Machine-readable scores for downstream analysis         |
| `raw_data_{targets}_vs_{diseases}.json`   | Raw API responses for reproducibility and auditing      |

Example: `validation_braf_kras_vs_cutaneous_melanoma_non-small_cell_lung_carcinoma.html`

Names use the resolved target symbols and disease names **in the order you gave them**,
so re-running the same validation overwrites the same files. A file can still get a new
name if Open Targets renames a disease (e.g. "coronary artery disease" became
"coronary artery disorder").

`/validate-team` also writes to `results/agents/{run_id}/`: one JSON dossier per
evidence agent, plus `consensus.md`.

## Project Structure

```
target_ai/
├── CLAUDE.md                          # Core instructions (read by Claude Code on startup)
├── README.md                          # This file
├── validate.md                        # /validate skill
├── SKILL.md                           # /score-target skill (single-target deep dive)
├── dossier_schema.md                  # JSON contract for multi-agent evidence dossiers
├── open_targets_client.py             # Python API client, scoring engine + HTML generator
├── report_template.html               # HTML template for interactive scoring matrix
├── .claude/
│   ├── agents/                        # Subagents for /validate-team
│   │   ├── ot-resolver.md
│   │   ├── clinical-evidence.md
│   │   ├── druggability-evidence.md
│   │   ├── biology-evidence.md
│   │   ├── landscape-scout.md
│   │   └── evidence-synthesizer.md
│   └── skills/validate-team/SKILL.md  # /validate-team orchestration steps
└── results/                           # Generated outputs (committed as examples)
    ├── validation_*_vs_*.html         #   ★ Interactive HTML (primary deliverable)
    ├── validation_*_vs_*.md
    ├── scores_*_vs_*.csv
    ├── raw_data_*_vs_*.json
    └── agents/{run_id}/               #   /validate-team dossiers + consensus.md
```

> **Note:** If a slash command shows "Unknown skill", describe the request in plain
> language ("validate PCSK9 against coronary artery disease"). `CLAUDE.md` is always
> loaded and routes it to the right command.

## Customization

### Adjusting Weights

Edit `composite_score()` in `open_targets_client.py`, and the matching weights table in
`CLAUDE.md` so the agents use the same weights:

```
# Default: Clinical-heavy (industry standard)
Composite = (Clinical × 0.40) + (Druggability × 0.35) + (Pathway × 0.25)

# Academic / early discovery: Biology-heavy
Composite = (Clinical × 0.25) + (Druggability × 0.25) + (Pathway × 0.50)

# Repurposing focus: Druggability-heavy
Composite = (Clinical × 0.30) + (Druggability × 0.50) + (Pathway × 0.20)
```

### Adding New Scoring Criteria

The rubric tables in `CLAUDE.md` are what the agents and the deep-dive skill score
against. The script implements them in `score_clinical()`, `score_druggability()` and
`score_pathway()`, so change both together.

### Adding New Data Sources

Describe additional APIs or databases in `CLAUDE.md` under "Data Sources" and Claude Code will query them alongside Open Targets.

## Limitations

- **Scores are heuristic, not predictive.** They summarize available evidence but don't guarantee clinical success.
- **Data freshness.** Open Targets releases quarterly. Very recent trial results may require supplementary web search.
- **Bias toward well-studied targets.** Under-researched targets will score lower due to data scarcity, not necessarily due to poor biology.
- **Somatic mutation counts as strong genetics.** A somatic-mutation score above 0.5 gives Clinical 4 even with no trials (e.g. KRAS × cutaneous melanoma), where a reviewer might score 2.
- **Well-validated targets saturate at 5.0.** Approved drugs max out Clinical and Druggability, so the ranking barely separates established targets. Use `/validate-team` or `/score-target` to tell them apart.
- **Not a substitute for expert review.** Use this tool to prioritize and generate hypotheses, not as a final decision-maker.

## License

MIT — use freely for academic and commercial research.

## Acknowledgments

- [Open Targets Platform](https://platform.opentargets.org/) for the comprehensive, open-access target–disease association data
- [Anthropic Claude Code](https://docs.anthropic.com/en/docs/claude-code/overview) for the agentic coding framework
