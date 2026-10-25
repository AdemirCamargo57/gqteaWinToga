"""molecularOrbitals: the MO catalogue, orbital selections, grids, cubes and the energy files.

Orbital values are checked against Gaussian's cubegen (tests/data/orbitals/
*.cube.gz, regenerated only by tests/data/orbitals/make_reference.py) and
against the SCF density, Sum occ |phi|^2 = rho, for restricted, unrestricted
and restricted-open files.
"""
import dataclasses
import gzip
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import cubeFile as cf
import gaussianWavefunction as gw
import molecularOrbitals as mo
import molecularSurface as ms

SRC = Path(__file__).resolve().parent.parent
DATA = Path(__file__).resolve().parent / "data"
WATER = DATA / "fchk" / "water_631gd_5d.fchk"
CH3_U = DATA / "fchk" / "ch3_ub3lyp_631gd.fchk"
CH3_RO = DATA / "orbitals" / "ch3_rob3lyp_631gd.fchk"
COARSE = mo.preset_spacing("Coarse")


@pytest.fixture(scope="module")
def water():
    return mo.load_orbitals(WATER)


@pytest.fixture(scope="module")
def ch3_u():
    return mo.load_orbitals(CH3_U)


@pytest.fixture(scope="module")
def ch3_ro():
    return mo.load_orbitals(CH3_RO)


def test_the_module_does_not_import_toga():
    code = "import sys, molecularOrbitals; assert 'toga' not in sys.modules, 'toga imported'"
    subprocess.run([sys.executable, "-c", code], cwd=SRC, check=True)


class TestCatalogue:
    def test_restricted(self, water):
        assert water.kind == mo.RESTRICTED and water.spins == (mo.ALPHA,) and not water.unrestricted
        assert water.n_mo == 18 and water.homo == {mo.ALPHA: 5} and water.lumo(mo.ALPHA) == 6
        occupancies = [o.occupancy for o in water.orbitals[mo.ALPHA]]
        assert occupancies == [2.0] * 5 + [0.0] * 13
        labels = [o.label for o in water.orbitals[mo.ALPHA][:8]]
        assert labels == ["HOMO-4", "HOMO-3", "HOMO-2", "HOMO-1", "HOMO", "LUMO", "LUMO+1", "LUMO+2"]
        homo = water.orbital(mo.ALPHA, 5)
        assert homo.energy == pytest.approx(-0.49756346) and homo.energy_ev == pytest.approx(-0.49756346 * 27.211386245988)
        assert water.gap(mo.ALPHA) == pytest.approx(0.21457438 + 0.49756346)
        assert water.spin_text(mo.ALPHA) == ""

    def test_unrestricted(self, ch3_u):
        assert ch3_u.kind == mo.UNRESTRICTED and ch3_u.spins == (mo.ALPHA, mo.BETA)
        assert ch3_u.homo == {mo.ALPHA: 5, mo.BETA: 4}
        assert [o.occupancy for o in ch3_u.orbitals[mo.ALPHA][:6]] == [1, 1, 1, 1, 1, 0]
        assert [o.occupancy for o in ch3_u.orbitals[mo.BETA][:6]] == [1, 1, 1, 1, 0, 0]
        assert ch3_u.orbital(mo.BETA, 5).label == "LUMO" and ch3_u.orbital(mo.ALPHA, 5).label == "HOMO"
        assert ch3_u.spin_text(mo.BETA) == " beta"
        energies = gw.read_fchk(CH3_U).beta_orbital_energies
        assert ch3_u.orbital(mo.BETA, 5).energy == pytest.approx(energies[4])

    def test_restricted_open(self, ch3_ro):
        assert ch3_ro.kind == mo.RESTRICTED_OPEN and ch3_ro.spins == (mo.ALPHA,)
        assert [o.occupancy for o in ch3_ro.orbitals[mo.ALPHA][:7]] == [2, 2, 2, 2, 1, 0, 0]
        assert ch3_ro.homo[mo.ALPHA] == 5 and ch3_ro.orbital(mo.ALPHA, 5).label == "HOMO"

    def test_summary_names_the_frontier_orbitals(self, water, ch3_u):
        assert water.summary().startswith("RHF/6-31G(d), 10 electrons, 18 MOs (restricted).")
        assert "HOMO 5: -13.539 eV" in water.summary() and "gap 19." in water.summary()
        assert "Alpha HOMO 5" in ch3_u.summary() and "Beta HOMO 4" in ch3_u.summary()

    def test_a_file_without_orbitals_is_refused(self):
        wfn = dataclasses.replace(gw.read_fchk(WATER), alpha_mo=None)
        with pytest.raises(ValueError, match="no molecular orbitals"):
            mo.orbital_catalogue(wfn)
        wfn = dataclasses.replace(gw.read_fchk(WATER), alpha_orbital_energies=None)
        with pytest.raises(ValueError, match="Alpha Orbital Energies"):
            mo.orbital_catalogue(wfn)

    def test_inconsistent_energies_and_electrons_are_refused(self):
        wfn = gw.read_fchk(WATER)
        with pytest.raises(ValueError, match="17 alpha orbital energies for 18 MOs"):
            mo.orbital_catalogue(dataclasses.replace(wfn, alpha_orbital_energies=wfn.alpha_orbital_energies[:-1]))
        with pytest.raises(ValueError, match="do not fit"):
            mo.orbital_catalogue(dataclasses.replace(wfn, n_alpha=19, n_beta=19))
        u = gw.read_fchk(CH3_U)
        with pytest.raises(ValueError, match="Beta Orbital Energies"):
            mo.orbital_catalogue(dataclasses.replace(u, beta_orbital_energies=None))

    def test_unknown_orbital_number(self, water):
        with pytest.raises(ValueError, match="MO 19 does not exist"):
            water.orbital(mo.ALPHA, 19)

    def test_labels(self):
        assert [mo.orbital_label(n, 3) for n in (1, 3, 4, 6)] == ["HOMO-2", "HOMO", "LUMO", "LUMO+2"]


class TestOrbitalSpec:
    @pytest.mark.parametrize("text, expected", [
        ("5", [5]),
        ("HOMO", [5]),
        ("lumo", [6]),
        ("HOMO-1, LUMO+1", [4, 7]),
        ("HOMO-2:LUMO+2", [3, 4, 5, 6, 7, 8]),
        ("LUMO+1:HOMO", [5, 6, 7]),
        ("2-4 9", [2, 3, 4, 9]),
        ("1-3, 3-4", [1, 2, 3, 4]),
        ("10:LUMO", [6, 7, 8, 9, 10]),
    ])
    def test_valid(self, water, text, expected):
        assert mo.parse_orbital_spec(text, water) == expected

    def test_spin_matters_for_homo_and_lumo(self, ch3_u):
        assert mo.parse_orbital_spec("HOMO", ch3_u, mo.ALPHA) == [5]
        assert mo.parse_orbital_spec("HOMO", ch3_u, mo.BETA) == [4]

    @pytest.mark.parametrize("text, message", [
        ("", "Enter the orbitals"),
        ("HOMO-9", "MO -4 does not exist"),
        ("19", "MO 19 does not exist"),
        ("SOMO", "not an orbital"),
        ("5-2", "not a valid MO range"),
        ("0-3", "not a valid MO range"),
        ("HOMO:", "two orbitals joined"),
        ("1:2:3", "two orbitals joined"),
    ])
    def test_invalid(self, water, text, message):
        with pytest.raises(ValueError, match=message):
            mo.parse_orbital_spec(text, water)


def _cubegen(name):
    path = DATA / "orbitals" / f"{name}.cube.gz"
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        cube = Path(tmp) / "ref.cube"
        cube.write_bytes(gzip.decompress(path.read_bytes()))
        return cf.read_cube_volumes(cube)[0]


class TestAgainstCubegen:
    @pytest.mark.parametrize("name, path, spin, number", [
        ("water_5d_mo5", WATER, mo.ALPHA, 5),
        ("water_5d_mo6", WATER, mo.ALPHA, 6),
        ("ch3_u_alpha5", CH3_U, mo.ALPHA, 5),
        ("ch3_u_beta5", CH3_U, mo.BETA, 5),
        ("ch3_ro_mo5", CH3_RO, mo.ALPHA, 5),
    ])
    def test_orbital_values_and_phases_match_cubegen(self, monkeypatch, name, path, spin, number):
        monkeypatch.setattr(mo, "SLAB_POINTS", 5000)  # several slabs
        ref = _cubegen(name)
        assert ref.orbital == number and np.allclose(ref.axes, ref.axes[0, 0] * np.eye(3))
        catalogue = mo.load_orbitals(path)
        values = mo.evaluate_orbitals(catalogue, spin, [number], ref.origin, float(ref.axes[0, 0]),
                                      ref.values.shape)[0]
        points = ms.grid_points(ref.origin, float(ref.axes[0, 0]), ref.values.shape)
        near = np.linalg.norm(points[:, None] - catalogue.wavefunction.coordinates[None], axis=2).min(axis=1) < 0.3
        diff = np.abs(values.ravel() - ref.values.ravel())
        # cubegen writes 6 significant digits and approximates near the nuclei.
        assert diff[~near].max() < 3e-5 and diff.max() < 1e-4
        strong = np.abs(ref.values) > 1e-3
        assert np.all(np.sign(values[strong]) == np.sign(ref.values[strong]))  # GaussView's phases


def _occupation_density(catalogue):
    """P = Sum_spin Sum_i n_i C_i C_i^T from the catalogue's occupancies (fchk basis)."""
    wfn = catalogue.wavefunction
    total = 0.0
    for spin in catalogue.spins:
        c = gw.mo_coefficients(wfn, spin)
        n = np.array([o.occupancy for o in catalogue.orbitals[spin]])
        total = total + (c * n) @ c.T
    return total


class TestGrids:
    @pytest.mark.parametrize("path", [WATER, CH3_U, CH3_RO])
    def test_occupancies_count_every_electron(self, path):
        import espIntegrals as ei
        catalogue = mo.load_orbitals(path)
        wfn = catalogue.wavefunction
        p_cart = wfn.transform.T @ _occupation_density(catalogue) @ wfn.transform
        assert np.trace(p_cart @ ei.overlap_matrix(wfn)) == pytest.approx(wfn.n_electrons, abs=1e-6)

    @pytest.mark.parametrize("path", [WATER, CH3_U])
    def test_occupancies_reproduce_the_fchk_density(self, path):
        catalogue = mo.load_orbitals(path)
        assert np.allclose(_occupation_density(catalogue), catalogue.wavefunction.density, atol=1e-7)

    def test_the_g09_rohf_fchk_density_is_the_alpha_density_only(self, ch3_ro):
        # Found while writing these tests: G09W's formchk writes the alpha
        # density as "Total SCF Density" for a ROHF job (it integrates to 5 of
        # the 9 electrons). The Orbitals tab never reads that section; the
        # occupancies above count all 9 electrons.
        wfn = ch3_ro.wavefunction
        c = wfn.alpha_mo[:, :5]
        assert np.allclose(c @ c.T, wfn.density, atol=1e-8)

    @pytest.mark.parametrize("path", [WATER, CH3_U, CH3_RO])
    def test_occupied_orbitals_on_the_grid_add_up_to_the_density(self, path):
        catalogue = mo.load_orbitals(path)
        total = None
        for spin in catalogue.spins:
            occupied = [o for o in catalogue.orbitals[spin] if o.occupancy > 0]
            volumes = mo.orbital_volumes(catalogue, spin, [o.number for o in occupied], 0.4, 0.05)
            part = sum(o.occupancy * v.values ** 2 for o, v in zip(occupied, volumes))
            total = part if total is None else total + part
            origin, shape = volumes[0].origin, volumes[0].values.shape
        points = ms.grid_points(origin, 0.4, shape)
        wfn = catalogue.wavefunction
        wfn = dataclasses.replace(wfn, cart_density=wfn.transform.T @ _occupation_density(catalogue) @ wfn.transform,
                                  packed=None)
        rho = gw.electron_density(wfn, points).reshape(shape)
        assert np.max(np.abs(total - rho)) < 1e-8 * max(1.0, rho.max())

    def test_the_shared_grid_is_large_enough_for_every_orbital(self, water):
        numbers = [5, 6, 8]
        volumes = mo.orbital_volumes(water, mo.ALPHA, numbers, COARSE, 0.02)
        assert len({v.values.shape for v in volumes}) == 1
        for volume in volumes:
            v = volume.values
            faces = [v[0], v[-1], v[:, 0], v[:, -1], v[:, :, 0], v[:, :, -1]]
            assert max(np.abs(f).max() for f in faces) < 0.002
            assert np.allclose(volume.axes, COARSE * np.eye(3))
        assert volumes[0].orbital == 5 and "MO 5 (HOMO), occupancy 2" in volumes[0].title

    def test_one_orbital_grid_is_not_larger_than_a_shared_one(self, water):
        alone = mo.orbital_volumes(water, mo.ALPHA, [5], COARSE, 0.02)[0]
        shared = mo.orbital_volumes(water, mo.ALPHA, [5, 18], COARSE, 0.02)[0]
        assert np.prod(alone.values.shape) <= np.prod(shared.values.shape)

    def test_bad_requests(self, water):
        with pytest.raises(ValueError, match="grid spacing"):
            mo.orbital_volumes(water, mo.ALPHA, [5], 0.0, 0.02)
        with pytest.raises(ValueError, match="isovalue"):
            mo.orbital_volumes(water, mo.ALPHA, [5], COARSE, -1)
        with pytest.raises(ValueError, match="at least one"):
            mo.orbital_volumes(water, mo.ALPHA, [], COARSE, 0.02)
        with pytest.raises(ValueError, match="million points"):
            mo.orbital_volumes(water, mo.ALPHA, [5], 0.01, 0.02)
        with pytest.raises(ValueError, match="MO 30"):
            mo.orbital_volumes(water, mo.ALPHA, [30], COARSE, 0.02)

    def test_density_volume_squares_the_amplitude(self, water):
        volume = mo.orbital_volumes(water, mo.ALPHA, [5], COARSE, 0.02)[0]
        density = mo.density_volume(water, volume, mo.ALPHA, 5)
        assert np.array_equal(density.values, volume.values ** 2) and density.orbital is None
        assert density.title.startswith("|MO|^2 5 (HOMO)") and np.array_equal(density.origin, volume.origin)

    def test_progress_errors_never_abort(self, water):
        def bad(fraction):
            raise RuntimeError("boom")
        assert mo.orbital_volumes(water, mo.ALPHA, [5], COARSE, 0.02, progress=bad)

    def test_presets(self):
        assert mo.preset_spacing("Medium") == pytest.approx(1 / 6) and mo.DEFAULT_GRID == "Medium"
        with pytest.raises(ValueError, match="grid quality"):
            mo.preset_spacing("Ultra")


class TestFiles:
    def test_cube_round_trip(self, water, tmp_path):
        paths = mo.write_cubes(water, str(tmp_path), mo.ALPHA, [5, 6], COARSE, 0.02)
        assert [Path(p).name for p in paths] == ["water_631gd_5d_MO5.cube", "water_631gd_5d_MO6.cube"]
        direct = mo.orbital_volumes(water, mo.ALPHA, [5, 6], COARSE, 0.02)
        for path, volume in zip(paths, direct):
            back = cf.read_cube_volumes(path)[0]
            assert back.values.shape == volume.values.shape
            assert np.allclose(back.values, volume.values, rtol=1e-5, atol=1e-10)
            assert np.allclose(back.origin, volume.origin)
            assert "MO" in back.title and back.title.isascii()

    def test_batches_give_the_same_cubes(self, water, tmp_path, monkeypatch):
        whole = mo.write_cubes(water, str(tmp_path / "a"), mo.ALPHA, [4, 5, 6], COARSE, 0.02)
        monkeypatch.setattr(mo, "VALUES_PER_BATCH", 1)  # one MO per batch
        batched = mo.write_cubes(water, str(tmp_path / "b"), mo.ALPHA, [4, 5, 6], COARSE, 0.02)
        assert [Path(p).read_bytes() for p in whole] == [Path(p).read_bytes() for p in batched]

    def test_unrestricted_density_cube_names(self, ch3_u, tmp_path):
        path = mo.write_cubes(ch3_u, str(tmp_path), mo.BETA, [5], COARSE, 0.02, density=True)[0]
        assert Path(path).name == "ch3_ub3lyp_631gd_MO5_beta_density.cube"
        assert cf.read_cube_volumes(path)[0].values.min() >= 0.0

    def test_too_many_cubes_are_refused(self, water, tmp_path, monkeypatch):
        monkeypatch.setattr(mo, "MAX_EXPORT_ORBITALS", 2)
        with pytest.raises(ValueError, match="at most 2"):
            mo.write_cubes(water, str(tmp_path), mo.ALPHA, [1, 2, 3], COARSE, 0.02)
        assert list(tmp_path.iterdir()) == []

    def test_write_errors_are_value_errors(self, water, tmp_path):
        blocker = tmp_path / "file"
        blocker.write_text("x")
        with pytest.raises(ValueError, match="Cannot write"):
            mo.write_cubes(water, str(blocker / "sub"), mo.ALPHA, [5], COARSE, 0.02)
        with pytest.raises(ValueError, match="Cannot write"):
            mo.write_energy_table(water, str(blocker / "sub"))

    def test_energy_table(self, water, ch3_u, tmp_path):
        path = mo.write_energy_table(water, str(tmp_path))
        assert Path(path).name == "water_631gd_5d_MO_energies.txt"
        text = Path(path).read_text(encoding="utf-8")
        assert text.startswith("# gQTEA Molecular Orbital Energies\n")
        assert "# wavefunction  restricted" in text and "# homo  5  -0.497563 Ha" in text
        assert "# mo spin occupancy label energy_hartree energy_ev" in text
        data = np.loadtxt(path, usecols=(0, 2, 4, 5), comments="#")
        assert data.shape == (18, 4) and data[4, 1] == 2 and data[5, 1] == 0
        assert np.allclose(data[:, 3], data[:, 2] * mo.HARTREE_TO_EV, atol=1e-4)
        u_text = mo.energy_table_text(ch3_u)
        assert "# homo_alpha  5" in u_text and "# homo_beta  4" in u_text and "# gap_beta" in u_text
        assert len([l for l in u_text.splitlines() if not l.startswith("#")]) == 2 * ch3_u.n_mo


class TestEnergyLevels:
    def test_entry(self, water, ch3_u, ch3_ro):
        entry = mo.energy_levels_entry(water)
        assert entry["type"] == "energy_levels" and len(entry["columns"]) == 1
        column = entry["columns"][0]
        assert column["label"] == "" and column["homo"] == 5 and column["lumo"] == 6
        assert len(column["levels"]) == 18
        lo, hi = entry["ylim"]
        energies = [l["energy_ev"] for l in column["levels"]]
        # HOMO-10 would be the O 1s core (MO 1, -559 eV): the window stops 20 eV below the HOMO instead.
        assert energies[0] < lo <= energies[4] - mo.LEVEL_SPAN_EV
        assert hi >= min(energies[15], energies[5] + mo.LEVEL_SPAN_EV)   # LUMO+10 = MO 16, or 20 eV up
        u = mo.energy_levels_entry(ch3_u)
        assert [c["label"] for c in u["columns"]] == ["Alpha", "Beta"] and u["columns"][1]["homo"] == 4
        assert u["ylim"][0] > -100.0   # the C 1s level (-277 eV) is outside the opening view
        assert mo.energy_levels_entry(ch3_ro)["columns"][0]["levels"][4]["occupancy"] == 1.0

    def test_the_plotter_accepts_and_describes_it(self, water, ch3_u, tmp_path):
        from plotter import PlotterBase
        path = mo.write_levels(mo.energy_levels_entry(ch3_u), tmp_path / mo.levels_filename("ch3"))
        figures = PlotterBase.load_json_plot_file(path)
        assert "42 orbital levels in 2 column(s)" in PlotterBase.describe_json_figures(figures)
        assert json.loads(Path(path).read_text())[0]["source"]["wavefunction"] == "unrestricted"

    @pytest.mark.parametrize("figure, message", [
        ({"type": "energy_levels"}, "'columns' must be a non-empty list"),
        ({"type": "energy_levels", "columns": [{"levels": []}]}, "column 1: 'levels'"),
        ({"type": "energy_levels", "columns": [{"levels": [{"number": 1}]}]}, "level 1 needs"),
        ({"type": "energy_levels", "columns": [{"levels": [{"number": 1, "energy_ev": "x", "occupancy": 2}]}]},
         "only numbers"),
        ({"type": "energy_levels", "columns": [{"levels": [{"number": 1, "energy_ev": -1, "occupancy": 2}],
                                                "homo": 1.5}]}, "'homo' must be an MO number"),
    ])
    def test_the_plotter_rejects_bad_entries(self, figure, message):
        from plotter import PlotterBase
        with pytest.raises(ValueError, match=message):
            PlotterBase.validate_plot_manifest([figure])

    def test_it_renders(self, ch3_ro, ch3_u):
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure
        from displayPlots import ENERGY_LEVEL_COLOURS, _degenerate_groups, draw_energy_levels
        for catalogue in (ch3_ro, ch3_u):
            figure = Figure()
            FigureCanvasAgg(figure)
            axes = figure.add_subplot()
            entry = mo.energy_levels_entry(catalogue)
            draw_energy_levels(axes, entry)
            figure.canvas.draw()
            assert tuple(axes.get_ylim()) == pytest.approx(entry["ylim"])
            texts = [t.get_text() for t in axes.texts]
            assert any(t.startswith("HOMO ") for t in texts) and any(t.startswith("gap ") for t in texts)
        # ROHF: MOs 3 and 4 are degenerate and share one slot; MO 5 is singly occupied (orange).
        groups = _degenerate_groups(mo.energy_levels_entry(ch3_ro)["columns"][0]["levels"])
        assert [l["number"] for l in groups[2]] == [3, 4]
        figure = Figure()
        FigureCanvasAgg(figure)
        axes = figure.add_subplot()
        draw_energy_levels(axes, mo.energy_levels_entry(ch3_ro))
        from matplotlib.colors import to_hex
        colours = {to_hex(c.get_colors()[0]) for c in axes.collections if c.get_colors().size}
        assert to_hex(ENERGY_LEVEL_COLOURS["singly"]) in colours
