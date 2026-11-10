# Evidence Dossier Schema

Shared contract between the `/validate-team` orchestrator and its subagents. Every
evidence agent writes exactly one JSON file to
`results/agents/{run_id}/{agent_name}.json` with this shape:

```json
{
  "agent": "clinical-evidence",
  "dimension": "clinical",
  "generated": "2026-11-03T14:05:00Z",
  "pairs": [
    {
      "target": "BRAF",
      "target_id": "ENSG00000157764",
      "disease": "cutaneous melanoma",
      "disease_id": "EFO_0000389",
      "proposed_score": 5,
      "rationale": "One or two sentences tying the evidence to the rubric row chosen.",
      "evidence": [
        {
          "claim": "Vemurafenib approved for BRAF V600E melanoma",
          "source": "OT:drugAndClinicalCandidates",
          "strength": "strong"
        }
      ],
      "gaps": ["No rare-variant burden data returned"]
    }
  ],
  "notes": "Optional free text: caveats, queries that failed, fields that were empty."
}
```

## Rules

- `dimension` is one of `clinical`, `druggability`, `pathway`, `landscape`.
- `proposed_score` is an integer 0–5 taken from the rubric in `CLAUDE.md`, or `null`
  for the `landscape` agent (which does not score).
- `source` is either an Open Targets field path prefixed `OT:` (e.g.
  `OT:target.tractability`) or a full URL for web evidence. No source → no claim.
- `strength` is `strong`, `moderate`, or `weak`.
- Include one entry in `pairs` for every target–disease pair in the brief, even when
  nothing was found. In that case use score 0 and explain it in `gaps`.
- If unsure between two rubric rows, pick the lower one and say why in `gaps`.
- Never fabricate evidence. Leave out anything you cannot cite.
