---
name: score-target
description: Deep-dive validation of a single gene target across one or more diseases. Use when the user wants detailed scoring, a deep-dive, SWOT analysis, or comprehensive profile of a single drug target. Triggers on "score target", "deep dive", "profile target", "analyze target".
---

# Deep-Dive Target Scoring

Perform an in-depth validation of a single target across one or more diseases.

## Usage

```
/score-target PCSK9 — cardiovascular disease, familial hypercholesterolemia
```

## Instructions

Given `$ARGUMENTS`:

1. **Parse.** Extract one target and one or more diseases from the input.

2. **Resolve.** Map to Ensembl and EFO IDs via the Open Targets search API.

3. **Deep data collection.** Go beyond the standard queries — also fetch:

   - Genetic constraint (pLI, LOEUF from gnomAD via Open Targets)
   - Mouse phenotypes (IMPC data from Open Targets)
   - Tissue expression (baseline expression from Open Targets)
   - Pharmacogenomics data if available
   - Safety liabilities and adverse events

4. **Supplementary search.** Use web search to find:

   - Latest clinical trial updates on ClinicalTrials.gov
   - Recent publications (last 2 years) about this target
   - Patent landscape (any composition-of-matter patents?)
   - Competitive landscape (who else is developing drugs against this target?)

5. **Score.** Apply the standard rubrics, but also provide sub-scores:

   - Clinical: Genetic evidence sub-score + Trial evidence sub-score
   - Druggability: Structural sub-score + Chemical matter sub-score + Modality sub-score
   - Pathway: Pathway centrality sub-score + Expression sub-score + Animal model sub-score

6. **Output.** Present a detailed profile card for the target including:

   - Summary table with all scores
   - Strengths / Weaknesses / Opportunities / Threats (SWOT) analysis
   - Comparison across diseases if multiple provided
   - Recommended next experiments

7. **Generate HTML report.** Run the Python script to produce the interactive HTML matrix:

   ```bash
   python open_targets_client.py validate "TARGET" "DISEASE1,DISEASE2"
   ```

   This creates `results/validation_{disease}.html` — the interactive scoring matrix.
   Tell the user to open it in their browser.

8. **Save** narrative to `results/{TARGET}_deep_dive.md` alongside the HTML.
