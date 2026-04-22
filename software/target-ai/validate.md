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

2. **Resolve identifiers.** For each target and disease, query the Open Targets Platform search API to get canonical Ensembl Gene IDs and EFO Disease IDs. Confirm with the user if any resolution is ambiguous (multiple hits with similar scores).

3. **Fetch data.** For every (target, disease) pair, run the GraphQL queries defined in CLAUDE.md:

   - Association scores (overall + per datatype + per datasource)
   - Target tractability and safety liabilities
   - Known drugs for the target–disease pair
   - Pathway annotations, GO terms, tissue expression

4. **Score.** Apply the 0–5 rubrics from CLAUDE.md for Clinical, Druggability, and Pathway dimensions. Compute the weighted Composite score.

5. **Rank and present.** Sort by Composite score descending. Output the full results table in Markdown. Follow with a narrative analysis for each pair.

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
