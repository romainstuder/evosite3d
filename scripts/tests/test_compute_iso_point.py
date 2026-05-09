"""Tests for scripts/compute_iso_point.py."""

from scripts.compute_iso_point import compute_properties


def test_compute_properties_emits_pi_and_mw(tmp_path, capsys) -> None:
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">prot1\nMKWVTFISLLLLFSSAYS\n")

    compute_properties(str(fasta))

    out = capsys.readouterr().out.strip().split("\t")
    assert out[0] == "prot1"
    pI = float(out[1])
    MW = float(out[2])
    # Just sanity-check the output ranges to keep the test BioPython-version stable.
    assert 0.0 < pI < 14.0
    assert MW > 1000.0  # 18 residues -> well over 1 kDa


def test_compute_properties_strips_gaps(tmp_path, capsys) -> None:
    fasta_gapped = tmp_path / "gapped.fasta"
    fasta_plain = tmp_path / "plain.fasta"
    fasta_gapped.write_text(">prot1\nMKW-VTFI--SLLLLFSSAYS\n")
    fasta_plain.write_text(">prot1\nMKWVTFISLLLLFSSAYS\n")

    compute_properties(str(fasta_gapped))
    gapped = capsys.readouterr().out

    compute_properties(str(fasta_plain))
    plain = capsys.readouterr().out

    assert gapped == plain
