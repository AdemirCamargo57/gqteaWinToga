"""QTAIM analysis: agreement with Multiwfn 2026.7.15, report and outputs (qtaimAnalysis)."""
import subprocess
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

import qtaimAnalysis as qa
import qtaimTopology as qt
from qtaim_fixtures import REFERENCE_FIXTURES, fchk, load_reference

SRC = Path(__file__).resolve().parent.parent
_RESULTS = {}


def results_for(name):
    if name not in _RESULTS:
        _RESULTS[name] = qa.QTAIMAnalyser().run(str(fchk(name)))
    return _RESULTS[name]


def prop_vector(cp):
    p = cp.properties
    return np.array([p.rho, p.laplacian, p.g, p.k, p.v, p.h, *p.eigenvalues, p.ellipticity, p.eta])


def match(results, ref_cps):
    """Our CP for each Multiwfn CP (nearest position)."""
    ours = np.array([c.position for c in results.critical_points])
    out = {}
    for r in ref_cps:
        d = np.linalg.norm(ours - r.position, axis=1)
        out[r.index] = (results.critical_points[int(np.argmin(d))], float(d.min()))
    return out


@pytest.mark.parametrize("name", REFERENCE_FIXTURES)
class TestAgreementWithMultiwfn:
    def test_same_critical_points(self, name):
        ref, _ = load_reference(name)
        res = results_for(name)
        assert Counter(r.kind for r in ref) == Counter(c.kind for c in res.critical_points)
        for r in ref:
            cp, distance = match(res, ref)[r.index]
            assert cp.kind == r.kind
            assert distance < 1e-7, f"CP {r.index} {r.kind} is {distance:.2e} bohr from Multiwfn's"

    def test_same_properties(self, name):
        ref, _ = load_reference(name)
        res = results_for(name)
        for r in ref:
            cp, _ = match(res, ref)[r.index]
            ours = prop_vector(cp)
            # Multiwfn prints 10 significant digits; ellipticity and eta with 6 decimals
            assert ours[:9] == pytest.approx(r.props[:9], rel=1e-7, abs=1e-10), f"CP {r.index} {r.kind}"
            assert ours[9:] == pytest.approx(r.props[9:], abs=2e-6), f"CP {r.index} {r.kind}"

    def test_same_bond_atoms(self, name):
        ref, _ = load_reference(name)
        res = results_for(name)
        for r in ref:
            if r.kind != qt.BCP:
                continue
            cp, _ = match(res, ref)[r.index]
            assert cp.assignment[0] == "pair" and frozenset(cp.assignment[1:]) == r.atoms

    def test_same_path_lengths(self, name):
        ref, ref_paths = load_reference(name)
        res = results_for(name)
        mapping = match(res, ref)
        ours = sorted((p.start_cp, p.end_cp, round(p.length, 5)) for p in res.bond_paths + res.ring_paths)
        theirs = sorted((mapping[a][0].index, mapping[b][0].index, round(length, 5)) for a, b, length in ref_paths)
        assert len(ours) == len(theirs)
        for (a1, b1, l1), (a2, b2, l2) in zip(ours, theirs):
            assert (a1, b1) == (a2, b2) and l1 == pytest.approx(l2, abs=5e-4)

    def test_poincare_hopf(self, name):
        assert results_for(name).poincare_hopf() == (1, True)


class TestReport:
    def test_header_and_sections(self):
        text = qa.QTAIMAnalyser().report_text(results_for("water_631gd_6d.fchk"))
        lines = text.splitlines()
        assert lines[0] == "# gQTEA QTAIM Topology Analysis"
        keys = [line[2:25].strip() for line in lines if line.startswith("# ") and line[2:25].strip()]
        assert keys[2:] == ["file", "method / basis", "density", "electrons", "atoms", "search_modes",
                            "cp_search", "path_generation", "engine", "units"]
        assert "# density                Total SCF Density" in lines
        assert "# search_modes           nuclei, atom pairs, triangles, tetrahedra" in lines
        for heading in ("CRITICAL POINT SUMMARY", "CRITICAL POINTS", "DENSITY PROPERTIES AT CRITICAL POINTS (a.u.)",
                        "BOND PATHS (each (3,-1) to its two (3,-3); lengths along the path; ratio = total / ncp_distance)"):
            assert heading in lines
        assert "  Poincare-Hopf  3 - 2 + 0 - 0 = 1   satisfied" in lines
        assert "RING PATHS" not in text   # no cage point in water

    def test_one_row_per_cp_and_assignments(self):
        text = qa.QTAIMAnalyser().report_text(results_for("water_631gd_6d.fchk"))
        block = text.split("CRITICAL POINTS\n")[1].split("\n\n")[0].splitlines()[1:]
        assert [row.split()[-1] for row in block] == ["O1", "H2", "H3", "O1-H2", "O1-H3"]
        assert "3.5679964" in text   # rho at the O-H BCP, Multiwfn 0.3567996402

    def test_bond_path_rows(self):
        text = qa.QTAIMAnalyser().report_text(results_for("water_631gd_6d.fchk"))
        rows = text.split("ratio = total / ncp_distance)\n")[1].split("\n\n")[0].splitlines()[1:]
        assert len(rows) == 2
        total, distance, ratio = (float(x) for x in rows[0].split()[-3:])
        assert distance == pytest.approx(0.923339, abs=1e-5) and 1.0 <= ratio < 1.001   # NCP-NCP distance
        assert total == pytest.approx((0.29034 + 1.45472) * qt.BOHR_RADIUS_ANGSTROM, abs=1e-4)

    def test_ring_paths_section_for_a_cage(self):
        text = qa.QTAIMAnalyser().report_text(results_for("bicyclopentane_631gd.fchk"))
        rows = text.split("RING PATHS ((3,+1) to (3,+3))\n")[1].split("\n\n")[0].splitlines()[1:]
        assert len(rows) == 3
        assert "(3,+1) 3   (3,+3) 1" in text

    def test_notes_and_written_lines(self, tmp_path):
        res = results_for("ch3_ub3lyp_631gd.fchk")
        written = qa.QTAIMAnalyser().write_outputs(res, str(tmp_path))
        text = Path(written[0]).read_text(encoding="utf-8")
        assert Path(written[0]).name == "ch3_ub3lyp_631gd_QTAIM.txt"
        assert "NOTE: Open-shell wavefunction" in text
        assert f"Written: {written[0]}" in text


class TestRunValidation:
    def test_ecp_note(self, monkeypatch):
        import gaussianWavefunction as gw
        real = gw.read_fchk

        def with_ecp(path):
            w = real(path)
            w.nuclear_charges = w.nuclear_charges.copy()
            w.nuclear_charges[0] = 6.0
            return w
        monkeypatch.setattr(gw, "read_fchk", with_ecp)
        res = qa.QTAIMAnalyser().run(str(fchk("water_631gd_6d.fchk")), modes=(2, 3))
        assert any("ECP atoms (O1)" in n for n in res.notes)

    def test_missing_cps_give_a_poincare_hopf_note(self):
        res = qa.QTAIMAnalyser().run(str(fchk("benzene_631gd.fchk")), modes=(2,))
        assert res.poincare_hopf() == (12, False)
        assert any("Poincare-Hopf sum is 12" in n for n in res.notes)

    @pytest.mark.parametrize("modes, message", [((), "at least one"), ((2, 6), "Unknown CP search mode")])
    def test_invalid_modes(self, modes, message):
        with pytest.raises(ValueError, match=message):
            qa.QTAIMAnalyser().run(str(fchk("water_631gd_6d.fchk")), modes=modes)

    def test_invalid_parameters(self):
        with pytest.raises(ValueError, match="step"):
            qa.QTAIMAnalyser().run(str(fchk("water_631gd_6d.fchk")), paths=qt.PathParameters(step=0))

    def test_unreadable_file(self, tmp_path):
        with pytest.raises(ValueError, match="Cannot read"):
            qa.QTAIMAnalyser().run(str(tmp_path / "missing.fchk"))

    def test_progress_is_reported_and_may_raise(self):
        seen = []

        def progress(fraction, message):
            seen.append(fraction)
            raise RuntimeError("must not abort")
        res = qa.QTAIMAnalyser().run(str(fchk("water_631gd_6d.fchk")), progress=progress)
        assert len(res.critical_points) == 5
        assert seen[0] == 0.0 and seen[-1] == 1.0 and seen == sorted(seen)

    def test_write_error_is_a_value_error(self, tmp_path):
        blocker = tmp_path / "file"
        blocker.write_text("x")
        with pytest.raises(ValueError, match="Could not write"):
            qa.QTAIMAnalyser().write_outputs(results_for("water_631gd_6d.fchk"), str(blocker / "sub"))

    def test_resolve_output_dir(self, tmp_path):
        f = tmp_path / "a.fchk"
        assert qa.resolve_output_dir(str(f), "  ") == str(tmp_path.resolve())
        assert qa.resolve_output_dir(str(f), "C:/out") == "C:/out"

    def test_numpy_and_numba_give_the_same_topology(self):
        pytest.importorskip("numba")
        a = qa.QTAIMAnalyser().run(str(fchk("benzene_631gd.fchk")), use_numba=False)
        b = qa.QTAIMAnalyser().run(str(fchk("benzene_631gd.fchk")), use_numba=True)
        assert [c.kind for c in a.critical_points] == [c.kind for c in b.critical_points]
        assert np.allclose([c.position for c in a.critical_points], [c.position for c in b.critical_points],
                           atol=1e-9)
        assert a.engine == "numpy" and b.engine == "numba (compiled)"


def test_modules_do_not_import_toga():
    code = "import sys, qtaimTopology, qtaimAnalysis; print('toga' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], cwd=SRC, capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "False"
