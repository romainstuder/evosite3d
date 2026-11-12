"""Tests for scripts/extract_seq_prob.py."""

from scripts.extract_seq_prob import parse_probabilities


def test_parse_probabilities_picks_correct_node_block(tmp_path) -> None:
    # The parser turns on at "Prob distribution at node N, by site"
    # and turns off at "Prob distribution at node N+1, by site".
    # Each captured row is split into tokens; columns 3..23 hold the per-AA probabilities.
    # The format used elsewhere in the codebase encodes each as "<aa>(0.NN)".
    probs = " ".join(f"{aa}(0.05)" for aa in "ACDEFGHIKLMNPQRSTVWY")
    rst = tmp_path / "rst.txt"
    rst.write_text(
        "Prob distribution at node 5, by site\n"
        f"  100 extra extra {probs}\n"
        "Prob distribution at node 6, by site\n"
        f"  101 extra extra {probs}\n"  # belongs to node 6, must be skipped
    )

    result = parse_probabilities(str(rst), node=5)

    assert list(result.keys()) == ["100"]
    assert len(result["100"]) == 20


def test_parse_probabilities_no_match_returns_empty(tmp_path) -> None:
    rst = tmp_path / "rst.txt"
    rst.write_text("nothing of interest here\n")

    assert parse_probabilities(str(rst), node=1) == {}
