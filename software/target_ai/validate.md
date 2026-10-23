---
name: validate
description: Run a full target validation analysis scoring multiple gene targets against diseases. Use when the user wants to validate, score, or rank drug targets against diseases. Triggers on phrases like "validate targets", "score targets", "rank targets", "target validation".
---

# Target Validation Command

Run a full target validation analysis.

## Usage

```
/validate BRAF,KRAS — melanoma,lung cancer
```

The argument format is: `TARGET1,TARGET2,... — DISEASE1,DISEASE2,...`
Use `—` (em-dash), `--`, or `vs` as the separator between targets and diseases.
Targets can be gene symbols (BRAF) or Ensembl IDs (ENSG00000157764).
Diseases can be common names (breast cancer) or EFO IDs (EFO_0000305).

## Instructions

Given the user's input in `$ARGUMENTS`:

1. **Parse inputs.** Split on the separator to extract target list and disease list. Trim whitespace. Handle both gene symbols and Ensembl IDs.

2. **(Optional) Sanity-check identifiers.** If a target or disease name looks ambiguous, use the `open-targets-platform` MCP tool `search_entities` to confirm the canonical Ensembl/EFO ID with the user before running the script — `validate` (step 6) will resolve them anyway, but this catches ambiguity early.

3. **Run the script.** Skip manual data fetching — `open_targets_client.py validate` (step 6) does resolution, fetching, and scoring internally in one deterministic pass. Only reach for the MCP tools (`query_open_targets_graphql`, `batch_query_open_targets_graphql`) if the user wants a field the script doesn't cover.

4. **Score.** The script applies the 0–5 rubrics from CLAUDE.md for Clinical, Druggability, and Pathway dimensions and computes the weighted Composite score — don't recompute this by hand.

5. **Rank and present.** The script already sorts by Composite score descending and prints the results table plus narrative. Print that output directly; don't recompute or reformat it.

6. **Generate outputs.** Run the Python script to produce all files:

   ```bash
   python open_targets_client.py validate "TARGET1,TARGET2" "DISEASE"
   ```

   This creates four files in `results/`:
   - `validation_report.md` — Markdown report
   - `scores.csv` — CSV for downstream analysis
   - `raw_data.json` — Raw API data
   - `validation_{disease}.html` — **Interactive HTML scoring matrix**

   The HTML file is the primary deliverable. Tell the user to open it in their browser.

Print the table and narrative directly in the conversation as well.
