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


class TestPlaneMap:
    @pytest.fixture(scope="class")
    def benzene_map(self):
        res = results_for("benzene_631gd.fchk")
        w = res.wavefunction
        plane = qt.plane_from_atoms(w.coordinates, *qt.default_plane_atoms(w.coordinates), symbols=w.atom_symbols())
        return res, plane, qa.QTAIMAnalyser().plane_map(res, plane, qa.PlaneMapOptions(points=60))

    def test_every_cp_of_a_planar_molecule_is_shown(self, benzene_map):
        res, _, m = benzene_map
        assert Counter(p[0] for p in m.points) == Counter({qt.NCP: 12, qt.BCP: 12, qt.RCP: 1})
        assert {p[3] for p in m.points if p[0] == qt.NCP} == {f"C{i}" for i in range(1, 7)} | {f"H{i}" for i in range(7, 13)}

    def test_lines(self, benzene_map):
        _, _, m = benzene_map
        kinds = Counter(k for k, _ in m.lines)
        assert kinds["gradient"] == 12 * 24 and kinds["bond_path"] == 24 and kinds["interbasin"] == 24
        for kind, pts in m.lines:
            assert kind in qa.MAP_LINE_KINDS and pts.shape[1] == 2 and len(pts) >= 2
            assert pts[:, 0].min() >= m.extent[0] - 0.1 and pts[:, 0].max() <= m.extent[1] + 0.1

    def test_gradient_lines_start_next_to_a_nucleus(self, benzene_map):
        _, _, m = benzene_map
        nuclei = np.array([(p[1], p[2]) for p in m.points if p[0] == qt.NCP])
        for kind, pts in m.lines:
            if kind == "gradient":
                d = np.linalg.norm(nuclei - pts[0], axis=1).min()
                assert d == pytest.approx(qa.GRADIENT_LINE_START_RADIUS * qt.BOHR_RADIUS_ANGSTROM, abs=1e-3)

    def test_grid_and_contours(self, benzene_map):
        _, _, m = benzene_map
        assert m.grid.shape == (60, 60) and np.all(m.grid > 0)
        assert m.contour_levels and set(m.contour_levels) <= set(qa.BADER_CONTOURS)
        assert max(m.contour_levels) <= m.grid.max()

    def test_extent_covers_atoms_plus_margin(self, benzene_map):
        res, plane, m = benzene_map
        st = plane.to_plane(res.wavefunction.coordinates) * qt.BOHR_RADIUS_ANGSTROM
        assert m.extent[0] == pytest.approx(st[:, 0].min() - 2.0) and m.extent[3] == pytest.approx(st[:, 1].max() + 2.0)

    def test_manifest_is_json_and_complete(self, benzene_map):
        import json
        entry = json.loads(json.dumps(qa.plane_manifest(benzene_map[2])))
        assert entry["type"] == "plane_map"
        assert set(entry) == {"type", "title", "xlabel", "ylabel", "extent", "grid", "contours", "lines", "points"}
        assert len(entry["grid"]) == 60 and len(entry["grid"][0]) == 60
        assert entry["xlabel"] == "u (Angstrom)" and "plane of C1, C2, C3" in entry["title"]

    def test_switches_and_far_plane(self):
        res = results_for("water_631gd_6d.fchk")
        far = qa.QTAIMAnalyser().plane_map(res, qt.axis_plane("XY", 5.0),
                                           qa.PlaneMapOptions(points=20, contours=False, gradient_lines=False,
                                                              interbasin=False))
        assert far.grid is None and far.lines == [] and far.points == []
        assert qa.plane_manifest(far)["grid"] is None

    def test_water_in_its_own_plane(self):
        res = results_for("water_631gd_6d.fchk")
        # water lies in the YZ plane of its fchk (x = 0)
        m = qa.QTAIMAnalyser().plane_map(res, qt.axis_plane("YZ"), qa.PlaneMapOptions(points=30))
        assert Counter(p[0] for p in m.points) == Counter({qt.NCP: 3, qt.BCP: 2})
        assert qa.plane_manifest(m)["xlabel"] == "y (Angstrom)"

    @pytest.mark.parametrize("kwargs, message", [({"points": 5}, "20 to 1000"), ({"margin_angstrom": -1}, "margin"),
                                                 ({"lines_per_atom": 2}, "lines per atom"),
                                                 ({"show_distance_bohr": 0}, "distance")])
    def test_option_validation(self, kwargs, message):
        with pytest.raises(ValueError, match=message):
            qa.PlaneMapOptions(**kwargs).validate()


class TestPlanes:
    def test_plane_from_atoms(self):
        coords = np.array([[0, 0, 0], [2, 0, 0], [0, 3, 0.0]])
        p = qt.plane_from_atoms(coords, 0, 1, 2, symbols=["O", "H", "H"])
        assert np.allclose(p.u, [1, 0, 0]) and np.allclose(p.v, [0, 1, 0]) and np.allclose(p.normal, [0, 0, 1])
        assert np.allclose(p.origin, coords.mean(axis=0)) and p.label == "plane of O1, H2, H3"

    def test_collinear_atoms_are_rejected(self):
        coords = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [0, 1, 0.0]])
        with pytest.raises(ValueError, match="collinear"):
            qt.plane_from_atoms(coords, 0, 1, 2)
        assert qt.default_plane_atoms(coords) == (0, 1, 3)
        with pytest.raises(ValueError, match="three different atoms"):
            qt.plane_from_atoms(coords, 0, 0, 2)

    def test_linear_molecule_has_no_default_plane(self):
        with pytest.raises(ValueError, match="non-collinear"):
            qt.default_plane_atoms(np.array([[0, 0, 0], [0, 0, 1], [0, 0, 2.0]]))

    def test_axis_planes(self):
        p = qt.axis_plane("xz", 1.0)
        assert np.allclose(p.u, [1, 0, 0]) and np.allclose(p.v, [0, 0, 1]) and np.allclose(p.origin, [0, 1, 0])
        assert p.axis_names == ("x", "z") and "y = 0.529177" in p.label
        with pytest.raises(ValueError, match="Unknown plane"):
            qt.axis_plane("AB")
