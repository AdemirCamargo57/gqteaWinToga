"""Excitation files: Gaussian CIS/TDA/TD output and Multiwfn's plain-text format."""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import excitationData as ed
import gaussianWavefunction as gw
from nto_fixtures import FIXTURES, fixture

SRC = Path(__file__).resolve().parent.parent

RESTRICTED_TEXT = """Excited State   1   1   5.7945
       5 ->  6         0.70642

Excited State   2   1   7.8943
       5 ->  7         0.63860
       5 ->  8         0.30006

Excited State   3   1   7.8943
       5 ->  7        -0.30006
       5 ->  8         0.63860
       4 <-  8         0.01000
"""

UNRESTRICTED_TEXT = """Excited State   1      0     2.07774
       600B -> 601B        -0.676085
       598A -> 602A         0.454805
       600A -> 603A         0.416875

 Excited State   2      0     2.07792
       599B -> 601B         0.561757
       598A -> 601A         0.496762
       600B -> 603B         0.468381
"""


def gaussian_block(energies):
    lines = [" Excitation energies and oscillator strengths:", " "]
    for n, e in enumerate(energies, 1):
        lines.append(f" Excited State   {n}:      Singlet-A      {e:.4f} eV  200.00 nm  f=0.0100  <S**2>=0.000")
        lines.append("       2 ->  3         0.70000")
        lines.append("")
    return "\n".join(lines)


@pytest.fixture(scope="module")
def uracil():
    return ed.read_excitation_file(fixture("uracil")[1])


class TestGaussianOutput:
    def test_uracil_states(self, uracil):
        assert uracil.source == "gaussian" and not uracil.unrestricted and uracil.blocks == 1
        assert [s.number for s in uracil.states] == [1, 2, 3]
        s3 = uracil.state(3)
        assert (s3.label, s3.multiplicity) == ('Singlet-A"', 1)
        assert (s3.energy_ev, s3.wavelength_nm, s3.oscillator, s3.s_squared) == (6.0180, 206.02, 0.0, 0.0)
        assert s3.title() == '3: Singlet-A" 6.0180 eV 206.02 nm f=0.0000'

    def test_uracil_configuration_counts_match_multiwfn(self, uracil):
        # Multiwfn: "MO pairs: 954 / 1824 / 846"
        assert [len(s.configurations) for s in uracil.states] == [954, 1824, 846]

    def test_uracil_normalisation_terms_match_multiwfn(self, uracil):
        x2, y2 = uracil.state(3).norm_terms()
        assert x2 == pytest.approx(0.500203, abs=5e-7) and y2 == pytest.approx(0.000196, abs=5e-7)

    def test_three_digit_index_without_a_space(self, uracil):
        c = next(c for c in uracil.state(1).configurations if c.virtual == 103 and c.occupied == 10)
        assert c.coefficient == 0.00017 and c.spin == ed.SPIN_ALPHA and not c.deexcitation
        assert c.line > 0

    def test_counts_and_geometry(self, uracil):
        assert (uracil.n_basis, uracil.n_independent, uracil.n_alpha, uracil.n_beta) == (128, 128, 29, 29)
        assert uracil.atomic_numbers.tolist() == [6, 6, 6, 6, 1, 1, 1, 1, 7, 7, 8, 8]
        assert uracil.coordinates_angstrom.shape == (12, 3)

    def test_singlets_triplets_and_deexcitations(self):
        data = ed.read_excitation_file(fixture("h2co_td5050")[1])
        assert [s.multiplicity for s in data.states] == [3, 1, 3, 3, 3, 1, 1, 1]
        assert any(c.deexcitation for s in data.states for c in s.configurations)

    def test_cis_has_no_deexcitations(self):
        data = ed.read_excitation_file(fixture("h2co_cis")[1])
        assert len(data.states) == 3
        assert not any(c.deexcitation for s in data.states for c in s.configurations)

    def test_open_shell(self):
        data = ed.read_excitation_file(fixture("ch3")[1])
        assert data.unrestricted and (data.n_alpha, data.n_beta) == (5, 4)
        s1 = data.state(1)
        assert (s1.label, s1.multiplicity) == ('2.004-E"', 0)
        assert {c.spin for c in s1.configurations} == {ed.SPIN_ALPHA, ed.SPIN_BETA}

    def test_unknown_symmetry_label(self):
        data = ed.read_excitation_file(fixture("benzene")[1])
        assert data.state(1).label == "Singlet-?Sym" and data.state(1).multiplicity == 1

    def test_last_block_is_used(self, tmp_path):
        path = tmp_path / "opt.log"
        path.write_text(" Entering Gaussian System\n" + gaussian_block([5.0, 6.0]) + "\n stuff\n"
                        + gaussian_block([4.5, 5.5, 6.5]) + "\n Normal termination\n")
        data = ed.read_excitation_file(path)
        assert data.blocks == 2 and [s.energy_ev for s in data.states] == [4.5, 5.5, 6.5]

    def test_state_lookup(self, uracil):
        with pytest.raises(ValueError, match="state 7.*1-3"):
            uracil.state(7)


class TestPlainText:
    def test_restricted(self, tmp_path):
        path = tmp_path / "exc.txt"
        path.write_text(RESTRICTED_TEXT)
        data = ed.read_excitation_file(path)
        assert data.source == "text" and not data.unrestricted and len(data.states) == 3
        s3 = data.state(3)
        assert (s3.multiplicity, s3.energy_ev, len(s3.configurations)) == (1, 7.8943, 3)
        assert s3.configurations[2].deexcitation and s3.configurations[2].coefficient == 0.01
        assert data.n_basis is None and data.coordinates_angstrom is None
        assert s3.title() == "3: Singlet 7.8943 eV"

    def test_unrestricted(self, tmp_path):
        path = tmp_path / "exc.txt"
        path.write_text(UNRESTRICTED_TEXT)
        data = ed.read_excitation_file(path)
        assert data.unrestricted and data.state(1).multiplicity == 0
        c = data.state(1).configurations[0]
        assert (c.occupied, c.virtual, c.spin, c.coefficient) == (600, 601, ed.SPIN_BETA, -0.676085)


class TestErrors:
    @pytest.mark.parametrize("text, message", [
        ("hello\nworld\n", "no excited states"),
        ("Excited State   1   1   5.0\n   5A -> 6B   0.7\n", "spin-flip"),
        ("Excited State   1   1   5.0\n   5 -> 6   0.7\n   5 -> 6   0.1\n", "appears twice"),
        ("Excited State   1   1   5.0\n\nExcited State   2   1   6.0\n   5 -> 6   0.7\n", "state 1 lists no"),
        ("Excited State   1   1   5.0\n   5A -> 6A   0.7\n   4 -> 6   0.1\n", "mixes"),
    ])
    def test_bad_text(self, tmp_path, text, message):
        path = tmp_path / "bad.txt"
        path.write_text(text)
        with pytest.raises(ValueError, match=message):
            ed.read_excitation_file(path)

    def test_gaussian_output_without_excitations(self, tmp_path):
        path = tmp_path / "sp.log"
        path.write_text(" Entering Gaussian System, Link 0=g09\n SCF Done:  E(RHF) = -1.0\n Normal termination\n")
        with pytest.raises(ValueError, match="CIS/TDA/TD"):
            ed.read_excitation_file(path)

    def test_missing_file(self, tmp_path):
        with pytest.raises(ValueError, match="Cannot read"):
            ed.read_excitation_file(tmp_path / "none.log")


class TestAgainstWavefunction:
    @pytest.mark.parametrize("name", list(FIXTURES))
    def test_every_fixture_pair_matches(self, name):
        fchk, log = fixture(name)
        assert ed.check_against_wavefunction(ed.read_excitation_file(log), gw.read_fchk(fchk)) == []

    def test_other_molecule(self):
        data = ed.read_excitation_file(fixture("uracil")[1])
        with pytest.raises(ValueError, match="atoms"):
            ed.check_against_wavefunction(data, gw.read_fchk(fixture("h2co_cis")[0]))

    def test_shifted_geometry(self):
        fchk, log = fixture("h2co_cis")
        wfn = gw.read_fchk(fchk)
        wfn.coordinates = wfn.coordinates + np.array([1e-3 / 0.529177210544, 0.0, 0.0])
        with pytest.raises(ValueError, match="same job"):
            ed.check_against_wavefunction(ed.read_excitation_file(log), wfn)

    def test_restricted_output_with_unrestricted_fchk(self):
        data = ed.read_excitation_file(fixture("h2co_cis")[1])
        with pytest.raises(ValueError, match="unrestricted"):
            ed.check_against_wavefunction(data, gw.read_fchk(fixture("ch3")[0]))

    def test_no_mo_coefficients(self):
        fchk, log = fixture("h2co_cis")
        wfn = gw.read_fchk(fchk)
        wfn.alpha_mo = None
        with pytest.raises(ValueError, match="MO coefficients"):
            ed.check_against_wavefunction(ed.read_excitation_file(log), wfn)

    def test_restricted_open_shell_is_refused(self):
        fchk, log = fixture("h2co_cis")
        wfn = gw.read_fchk(fchk)
        wfn.method = "ROB3LYP"
        with pytest.raises(ValueError, match="restricted open-shell"):
            ed.check_against_wavefunction(ed.read_excitation_file(log), wfn)

    def test_index_beyond_the_orbitals(self, tmp_path):
        wfn = gw.read_fchk(fixture("h2co_cis")[0])
        path = tmp_path / "exc.txt"
        path.write_text(f"Excited State   1   1   5.0\n    8 -> {wfn.n_mo + 1}   0.70\n")
        with pytest.raises(ValueError, match="state 1.*line 2"):
            ed.check_against_wavefunction(ed.read_excitation_file(path), wfn)

    def test_several_blocks_give_a_note(self, tmp_path):
        fchk, log = fixture("h2co_cis")
        text = Path(log).read_text()
        start = text.index(" Excitation energies and oscillator strengths:")
        path = tmp_path / "twice.log"
        path.write_text(text[:start] + text[start:] + text[start:])
        notes = ed.check_against_wavefunction(ed.read_excitation_file(path), gw.read_fchk(fchk))
        assert len(notes) == 1 and "last" in notes[0]


def test_sibling_excitation_file(tmp_path):
    fchk = tmp_path / "mol.fchk"
    fchk.write_text("")
    assert ed.sibling_excitation_file(fchk) is None
    (tmp_path / "mol.out").write_text("")
    assert ed.sibling_excitation_file(fchk) == tmp_path / "mol.out"
    (tmp_path / "mol.log").write_text("")
    assert ed.sibling_excitation_file(fchk) == tmp_path / "mol.log"


def test_does_not_import_toga():
    code = "import sys, excitationData; print('toga' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], cwd=SRC, capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "False"
