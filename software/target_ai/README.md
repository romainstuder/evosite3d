# 🧬 Target Validation Tool for Claude Code

An AI-powered drug target validation assistant that scores and ranks therapeutic targets against diseases using real-time data from Open Targets Platform and supplementary biomedical databases.

## What It Does

Given one or more **gene targets** and **diseases**, this tool:

1. **Resolves** gene symbols → Ensembl IDs and disease names → EFO IDs
2. **Fetches** evidence from Open Targets Platform (genetic associations, clinical trials, tractability, pathways, safety)
3. **Scores** each target–disease pair on three dimensions (0–5 each):
   - **Clinical Evidence** — GWAS hits, rare variants, clinical trial history
   - **Druggability** — small-molecule/antibody tractability, existing chemical matter
   - **Pathway/Biology** — pathway membership, tissue expression, animal models
4. **Ranks** by weighted composite score with compact narrative
5. **Generates** an interactive HTML scoring matrix plus Markdown report, CSV, and raw JSON

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

#### `validate` — Multi-target × multi-disease scoring

```
validate BRAF,KRAS — melanoma,lung cancer
validate PCSK9,HMGCR,ANGPTL3 — cardiovascular disease
validate CDK4,CDK6 — breast cancer, ovarian cancer
validate EGFR,HER2,HER3,MET — non-small cell lung cancer
compare-targets JAK1,JAK2,JAK3,TYK2 — rheumatoid arthritis
```

#### `find-targets` — Discover top targets for a disease

```
find-targets Alzheimer's disease
find-targets Crohn's disease 20
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

### Usage — Python Script (standalone)

The included Python script can also be run independently:

```bash
# Search for a target or disease
python3 open_targets_client.py search-target BRAF
python3 open_targets_client.py search-disease "breast cancer"

# Find top targets for a disease
python3 open_targets_client.py find-targets "Alzheimer disease"
python3 open_targets_client.py find-targets "Alzheimer disease" 20

# Full validation
python3 open_targets_client.py validate "BRAF,KRAS" "cutaneous melanoma,non-small cell lung carcinoma"

# Individual lookups
python3 open_targets_client.py target-info ENSG00000157764
python3 open_targets_client.py known-drugs ENSG00000157764
python3 open_targets_client.py association ENSG00000157764 EFO_0000389

# Pathway visualization
python3 open_targets_client.py pathways "PCSK9,HMGCR,ANGPTL3"
```

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

### Confidence Levels

| Level      | Criteria                                                    |
| ---------- | ----------------------------------------------------------- |
| **High**   | Composite ≥ 3.5 with data from ≥ 3 independent source types |
| **Medium** | Composite 2.0–3.49 or limited independent sources           |
| **Low**    | Composite < 2.0 or based primarily on literature mining     |

## Output Files

After each validation run, results are saved to the `results/` directory. File names encode both targets and diseases:

| File                                      | Contents                                                |
| ----------------------------------------- | ------------------------------------------------------- |
| `validation_{targets}_vs_{diseases}.html` | **★ Interactive HTML scoring matrix** — open in browser |
| `validation_{targets}_vs_{diseases}.md`   | Markdown report with tables and narrative               |
| `scores_{targets}_vs_{diseases}.csv`      | Machine-readable scores for downstream analysis         |
| `raw_data_{targets}_vs_{diseases}.json`   | Raw API responses for reproducibility and auditing      |

Example: `validation_braf_kras_vs_cutaneous_melanoma_non-small_cell_lung_carcinoma.html`

## Project Structure

```
target-validation-tool/
├── CLAUDE.md                          # Core instructions (read by Claude Code on startup)
├── README.md                          # This file
├── open_targets_client.py             # Python API client, scoring engine + HTML generator
├── report_template.html               # HTML template for interactive scoring matrix
└── results/                           # Generated outputs (gitignored)
    ├── validation_*_vs_*.html         #   ★ Interactive HTML (primary deliverable)
    ├── validation_*_vs_*.md
    ├── scores_*_vs_*.csv
    └── raw_data_*_vs_*.json
```

> **Note:** All commands (`validate`, `find-targets`, `score-target`, `compare-targets`)
> are handled via `CLAUDE.md` instructions — just type naturally and Claude will run the
> right Python command.

## Customization

### Adjusting Weights

Edit the composite formula in `CLAUDE.md` under "Step 4 — Compute Composite Score":

```
# Default: Clinical-heavy (industry standard)
Composite = (Clinical × 0.40) + (Druggability × 0.35) + (Pathway × 0.25)

# Academic / early discovery: Biology-heavy
Composite = (Clinical × 0.25) + (Druggability × 0.25) + (Pathway × 0.50)

# Repurposing focus: Druggability-heavy
Composite = (Clinical × 0.30) + (Druggability × 0.50) + (Pathway × 0.20)
```

### Adding New Scoring Criteria

Add new sub-scores to the rubric tables in `CLAUDE.md`. Claude Code will automatically incorporate them.

### Adding New Data Sources

Describe additional APIs or databases in `CLAUDE.md` under "Data Sources" and Claude Code will query them alongside Open Targets.

## Limitations

- **Scores are heuristic, not predictive.** They summarize available evidence but don't guarantee clinical success.
- **Data freshness.** Open Targets releases quarterly. Very recent trial results may require supplementary web search.
- **Bias toward well-studied targets.** Under-researched targets will score lower due to data scarcity, not necessarily due to poor biology.
- **Not a substitute for expert review.** Use this tool to prioritize and generate hypotheses, not as a final decision-maker.

## License

MIT — use freely for academic and commercial research.

## Acknowledgments

- [Open Targets Platform](https://platform.opentargets.org/) for the comprehensive, open-access target–disease association data
- [Anthropic Claude Code](https://docs.anthropic.com/en/docs/claude-code/overview) for the agentic coding framework
