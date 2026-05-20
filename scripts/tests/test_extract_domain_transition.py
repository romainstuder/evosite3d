"""Tests for scripts/extract_domain_transition.py."""

from scripts.extract_domain_transition import extract_transitions


def test_emits_only_target_branch_and_position(tmp_path, capsys) -> None:
    # The script's hard-coded targets:
    #   target_branches  = ["208..209", "209..210", "209..265", "208..230"]
    #   target_position  = ["198", "204", "222", "234", "251", "270", "272", "307", "318"]
    text = (
        "Summary of changes along branches\n"
        "Branch 1: 208..209\n"
        "198 some data\n"  # branch 208..209, pos 198 -> matches
        "999 ignored\n"  # not a target position
        "Branch 2: 999..999\n"  # not a target branch
        "198 should be ignored\n"
        "List of extant and reconstructed sequences\n"  # section ends
        "198 also ignored\n"
    )
    inp = tmp_path / "in.txt"
    inp.write_text(text)

    extract_transitions(str(inp))

    out = capsys.readouterr().out.splitlines()
    assert out == ["208..209 198 198 some data"]


def test_no_section_yields_nothing(tmp_path, capsys) -> None:
    inp = tmp_path / "in.txt"
    inp.write_text("Branch 1: 208..209\n198 data\n")

    extract_transitions(str(inp))
    assert capsys.readouterr().out == ""
