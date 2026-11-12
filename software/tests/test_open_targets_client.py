"""Tests for software/target_ai/open_targets_client.py (pure scoring logic)."""

from ..target_ai.open_targets_client import (
    _make_slug,
    composite_score,
    confidence_level,
    format_markdown_table,
    format_narrative,
    score_clinical,
    score_druggability,
    score_pathway,
)

# ── composite_score ────────────────────────────────────────────


def test_composite_score_weights() -> None:
    # 5*0.40 + 5*0.35 + 5*0.25 = 5
    assert composite_score(5, 5, 5) == 5.0


def test_composite_score_zero() -> None:
    assert composite_score(0, 0, 0) == 0.0


def test_composite_score_partial() -> None:
    # 4*0.40 + 3*0.35 + 2*0.25 = 1.6 + 1.05 + 0.5 = 3.15
    assert composite_score(4, 3, 2) == 3.15


# ── confidence_level ───────────────────────────────────────────


def test_confidence_level_high_requires_score_and_sources() -> None:
    assert confidence_level(4.0, 3) == "High"


def test_confidence_level_high_score_alone_not_enough() -> None:
    # Three sources required for High.
    assert confidence_level(4.0, 2) == "Medium"


def test_confidence_level_medium_threshold() -> None:
    assert confidence_level(2.0, 1) == "Medium"


def test_confidence_level_low() -> None:
    assert confidence_level(1.0, 5) == "Low"


# ── score_clinical ─────────────────────────────────────────────


def test_score_clinical_no_data_returns_zero() -> None:
    score, reason = score_clinical({}, {}, "EFO_0000305")
    assert score == 0
    assert "No clinical or genetic" in reason


def test_score_clinical_phase_iv_drug_returns_5() -> None:
    drugs = {
        "data": {
            "target": {
                "drugAndClinicalCandidates": {
                    "count": 1,
                    "rows": [
                        {
                            "maxClinicalStage": "Phase IV",
                            "drug": {"name": "TestDrug"},
                        }
                    ],
                }
            }
        }
    }
    score, _ = score_clinical({}, drugs, "EFO_0000305")
    assert score == 5


def test_score_clinical_phase_iii_returns_4() -> None:
    drugs = {
        "data": {
            "target": {
                "drugAndClinicalCandidates": {
                    "count": 1,
                    "rows": [{"maxClinicalStage": "Phase III", "drug": {"name": "X"}}],
                }
            }
        }
    }
    score, _ = score_clinical({}, drugs, "EFO_0000305")
    assert score == 4


def test_score_clinical_strong_genetic_signal_returns_4() -> None:
    assoc = {
        "data": {
            "disease": {
                "associatedTargets": {
                    "rows": [
                        {
                            "datatypeScores": [
                                {"id": "genetic_association", "score": 0.7},
                            ]
                        }
                    ]
                }
            }
        }
    }
    score, reason = score_clinical(assoc, {}, "EFO_0000305")
    assert score == 4
    assert "Genetic" in reason


# ── score_druggability ─────────────────────────────────────────


def test_score_druggability_no_data_returns_zero() -> None:
    score, _ = score_druggability({}, {})
    assert score == 0


def test_score_druggability_approved_drug_returns_5() -> None:
    target_info = {
        "data": {"target": {"tractability": [{"modality": "SM", "value": True, "label": "good"}]}}
    }
    drugs = {
        "data": {
            "target": {
                "drugAndClinicalCandidates": {
                    "count": 1,
                    "rows": [{"maxClinicalStage": "Approved"}],
                }
            }
        }
    }
    score, _ = score_druggability(target_info, drugs)
    assert score == 5


def test_score_druggability_tractable_only_returns_3() -> None:
    target_info = {
        "data": {
            "target": {
                "tractability": [{"modality": "SM", "value": True, "label": "Clinical Precedence"}]
            }
        }
    }
    score, _ = score_druggability(target_info, {})
    assert score == 3


# ── score_pathway ──────────────────────────────────────────────


def test_score_pathway_no_data_returns_zero() -> None:
    score, _ = score_pathway({}, {})
    assert score == 0


def test_score_pathway_many_pathways_returns_5() -> None:
    target_info = {
        "data": {
            "target": {
                "pathways": [{"pathway": f"pw{i}"} for i in range(6)],
            }
        }
    }
    assoc = {
        "data": {
            "disease": {
                "associatedTargets": {
                    "rows": [{"datatypeScores": [{"id": "rna_expression", "score": 0.5}]}]
                }
            }
        }
    }
    score, _ = score_pathway(target_info, assoc)
    assert score == 5


# ── _make_slug ─────────────────────────────────────────────────


def test_make_slug_simple() -> None:
    results = [{"target": "BRAF", "disease": "Melanoma"}]
    assert _make_slug(results) == "braf_vs_melanoma"


def test_make_slug_strips_apostrophes_and_spaces() -> None:
    results = [{"target": "Foo's Gene", "disease": "Some Disease"}]
    assert _make_slug(results) == "foos_gene_vs_some_disease"


def test_make_slug_deduplicates() -> None:
    results = [
        {"target": "BRAF", "disease": "Melanoma"},
        {"target": "BRAF", "disease": "Melanoma"},  # duplicate
        {"target": "KRAS", "disease": "Melanoma"},
    ]
    assert _make_slug(results) == "braf_kras_vs_melanoma"


# ── Output formatters ──────────────────────────────────────────


def test_format_markdown_table_ranks_results() -> None:
    results = [
        {
            "target": "BRAF",
            "disease": "melanoma",
            "clinical": 5,
            "druggability": 4,
            "pathway": 3,
            "composite": 4.25,
            "confidence": "High",
        }
    ]
    md = format_markdown_table(results)
    assert "| 1 | **BRAF** | melanoma | 5 | 4 | 3 | **4.25** | High |" in md
    assert "Composite" in md


def test_format_narrative_includes_reasons() -> None:
    results = [
        {
            "target": "BRAF",
            "disease": "melanoma",
            "clinical": 5,
            "druggability": 4,
            "pathway": 3,
            "safety": 5,
            "composite": 4.25,
            "confidence": "High",
            "clinical_reason": "Approved drug",
            "druggability_reason": "SM tractable",
            "pathway_reason": "5 pathways",
            "safety_reason": "no liabilities",
        }
    ]
    out = format_narrative(results)
    assert "BRAF×melanoma" in out
    assert "Approved drug" in out
    assert "no liabilities" in out
