import asyncio

import numpy as np
import pytest

import periodicPreOptimizer as ppo
import solventPacking
from sharedWallBox import SharedWallBox, SharedWallBoxUI

WATER = [["O", 0.0, 0.0, 0.1173], ["H", 0.0, 0.7572, -0.4692], ["H", 0.0, -0.7572, -0.4692]]
METHANOL = [["C", -0.0465, 0.6630, 0.0], ["O", -0.0465, -0.7570, 0.0], ["H", -1.0861, 0.9729, 0.0],
            ["H", 0.4366, 1.0816, 0.8904], ["H", 0.4366, 1.0816, -0.8904], ["H", 0.8624, -1.0857, 0.0]]
CO2 = [["C", 0.0, 0.0, 0.0], ["O", 0.0, 0.0, 1.16], ["O", 0.0, 0.0, -1.16]]

needs_numba = pytest.mark.skipif(not solventPacking.numba_available(), reason="numba not installed")


def make_builder(ab=16.0, c1=12.0, c2=14.0, periodic=True, engine="numpy", attempts=2500,
                 solutes=(False, False), comps=(1.0, 0.5), seeds=(11, 12), scales=(0.8, 0.75),
                 padding=0.3, clearance=0.4, rotate=(True, True), density=(0.95, 0.9)):
    b = SharedWallBox()
    b.shared_z = c1
    b.lattice_vectors = [ab, ab, c1 + c2]
    b.use_periodic_minimum_image = periodic
    b.packing_engine = engine
    for i, box in enumerate(b.boxes):
        box.update(spacing=padding, target_density=density[i], max_attempts=attempts, vdw_scale=scales[i],
                   random_seed=seeds[i], min_distance=clearance, composition_A=comps[i],
                   composition_B=1.0 - comps[i], use_solute=solutes[i], rotate=rotate[i],
                   coords_A=[a[:] for a in WATER], coords_B=[a[:] for a in METHANOL],
                   solute_coords=[a[:] for a in CO2])
    return b


def run(engine, **kwargs):
    b = make_builder(engine=engine, **kwargs)
    assert b.build_box()
    return b


def coords(data):
    return np.array([a[1:] for a in data], dtype=float)


CASES = [
    dict(),
    dict(solutes=(True, True)),
    dict(periodic=False, solutes=(False, True)),
    dict(comps=(0.0, 1.0), rotate=(False, True), seeds=(5, 9)),
    dict(ab=28.0, c1=22.0, c2=24.0, attempts=2000),          # real cell lists
    dict(density=(0.0, 0.95)),                                # an empty first region
]


@pytest.mark.parametrize("case", CASES)
def test_numpy_engine_reproduces_the_legacy_structure_exactly(case):
    legacy = run("legacy", **case)
    fast = run("numpy", **case)
    for lb, fb in zip(legacy.boxes, fast.boxes):
        assert (fb["count_solvent_A"], fb["count_solvent_B"]) == (lb["count_solvent_A"], lb["count_solvent_B"])
        assert fb["solvent_attempts_used"] == lb["solvent_attempts_used"]
    assert [a[0] for a in fast.combined_data] == [a[0] for a in legacy.combined_data]
    assert np.array_equal(coords(fast.combined_data), coords(legacy.combined_data))
    assert len(legacy.combined_data) > 30


@needs_numba
@pytest.mark.parametrize("case", CASES)
def test_numba_engine_reproduces_the_legacy_structure_exactly(case):
    legacy = run("legacy", **case)
    fast = run("numba", **case)
    assert np.array_equal(coords(fast.combined_data), coords(legacy.combined_data))


def test_regions_and_molecules_are_recorded_in_cell_coordinates():
    b = run("numpy", solutes=(True, True))
    assert len(b.molecules) == len(b.molecule_regions)
    total = sum(len(m[1]) for m in b.molecules)
    assert total == len(b.combined_data)
    system = b.packed_system()
    centres = ppo.molecule_centres(system.coordinates, system.masses, system.molecule_slices)
    regions = np.array(b.molecule_regions)
    assert np.all(centres[regions == 1, 2] <= 12.0) and np.all(centres[regions == 2, 2] >= 12.0)
    # The written structure keeps the original z shift (lowest atom at z = 0).
    assert min(a[3] for a in b.combined_data) == pytest.approx(0.0)


@pytest.mark.parametrize("periodic", [True, False])
@pytest.mark.parametrize("mode", [ppo.MODE_RIGID, ppo.MODE_FLEXIBLE])
def test_preoptimization_keeps_each_region_on_its_side(periodic, mode):
    b = run("numpy", solutes=(True, False), periodic=periodic)
    b.preoptimize = True
    b.preopt_mode = mode
    b.run_preoptimization()
    result = b.preopt_result
    assert result is not None, b.preopt_error
    assert result.report_initial["clashes"] > 0
    assert result.report_final["clashes"] == 0
    assert result.connectivity_preserved is True
    system = b.packed_system()
    centres = ppo.molecule_centres(result.coordinates, system.masses, system.molecule_slices)
    regions = np.array(b.molecule_regions)
    tol = 1e-3
    assert np.all(centres[regions == 1, 2] >= -tol) and np.all(centres[regions == 1, 2] <= 12.0 + tol)
    assert np.all(centres[regions == 2, 2] >= 12.0 - tol) and np.all(centres[regions == 2, 2] <= 26.0 + tol)
    assert not any("outside their region" in note for note in result.notes)
    if not periodic:
        xyz = result.coordinates
        assert xyz[:, :2].min() > 0.3 - 0.05 and xyz[:, :2].max() < 16.0 - 0.3 + 0.05


def test_wall_bounds_have_no_wall_at_the_interface():
    b = run("numpy", periodic=False)
    system = b.packed_system()
    low, high = b.preopt_wall_bounds(system)
    regions = np.repeat(b.molecule_regions, [stop - start for start, stop in system.molecule_slices])
    assert np.all(np.isinf(high[regions == 1, 2])) and np.all(low[regions == 1, 2] == 0.3)
    assert np.all(np.isinf(low[regions == 2, 2])) and np.all(high[regions == 2, 2] == 26.0 - 0.3)


def test_hooks():
    b = run("numpy", solutes=(False, True), scales=(0.7, 0.85), clearance=0.6)
    assert b.packing_vdw_scale() == pytest.approx(0.85)
    assert b.effective_preopt_contact_scale() == pytest.approx(0.90)
    assert b.preopt_solute_clearance() == pytest.approx(0.6)
    b.boxes[1]["solute_charge"] = -1
    b.boxes[0]["charge_A"] = 1
    b.boxes[1]["charge_B"] = 2
    expected = -1 + b.boxes[0]["count_solvent_A"] + 2 * b.boxes[1]["count_solvent_B"]
    assert b.net_charge() == expected
    b.multiplicity = 1
    assert b.compute_system_properties().atom_count == len(b.combined_data)


def test_outputs(tmp_path):
    b = run("numpy", solutes=(True, False))
    b.output_dir = str(tmp_path)
    b.multiplicity = 1
    b.preoptimize = True
    b.run_preoptimization()
    b.compute_system_properties()
    summary = b.write_outputs()
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "shared_wall_box.txt", "shared_wall_box.xyz", "shared_wall_box_preopt.xyz"]
    lines = (tmp_path / "shared_wall_box.xyz").read_text(encoding="utf-8").splitlines()
    assert lines[1] == "Shared-wall double solvent box generated by density-driven random insertion"
    assert "Box 1 summary" in summary and "Total number of electrons" in summary
    assert "--- Pre-optimization ---" in summary and "cell coordinates" in summary


def test_cancel_stops_without_combined_output():
    b = make_builder(attempts=50000, density=(1.3, 1.3))

    def progress(payload):
        if payload["accepted"] >= 2:
            b.cancel_requested = True

    assert b.build_box(progress) is False
    assert b.combined_data == []


# ----------------------------------------------------------------------
# UI plumbing, headless
# ----------------------------------------------------------------------
class Field:
    def __init__(self, value=""):
        self.value = value
        self.enabled = True
        self.text = ""


def test_read_charges_per_box():
    ui = SharedWallBoxUI.__new__(SharedWallBoxUI)
    SharedWallBox.__init__(ui)
    values = {"textInput_box1_solute_charge": "1", "textInput_box2_solute_charge": "",
              "textInput_box1_charge_A": "-1", "textInput_box2_charge_A": "0",
              "textInput_box1_charge_B": "", "textInput_box2_charge_B": "2"}
    for name, value in values.items():
        setattr(ui, name, Field(value))
    ui.read_charges()
    assert [(b["solute_charge"], b["charge_A"], b["charge_B"]) for b in ui.boxes] == [(1, -1, 0), (0, 0, 2)]
    ui.textInput_box2_charge_B.value = "two"
    with pytest.raises(ValueError, match="Box 2 solvent B charge"):
        ui.read_charges()


@pytest.mark.parametrize("value,fragment", [("", "not inferred"), ("0", "positive")])
def test_missing_multiplicity_points_to_the_field(value, fragment, monkeypatch):
    import toga

    class Recorded:
        def __init__(self, title, message):
            self.title, self.message = title, message

    class Dialogs:
        def __init__(self):
            self.shown = []

        async def dialog(self, d):
            self.shown.append(d)

    class FocusField(Field):
        focused = 0

        def focus(self):
            self.focused += 1

    class Tabs:
        current_tab = 2

    monkeypatch.setattr(toga, "ErrorDialog", Recorded)
    ui = SharedWallBoxUI.__new__(SharedWallBoxUI)
    SharedWallBox.__init__(ui)
    ui.main_window = Dialogs()
    ui.tabs = Tabs()
    ui.textInput_multiplicity = FocusField(value)
    assert asyncio.run(ui.read_multiplicity()) is False
    message = ui.main_window.shown[-1].message
    assert fragment in message and "Shared geometry section at the top" in message
    assert ui.tabs.current_tab == ui.SETUP_TAB
    assert ui.textInput_multiplicity.focused == 1
    ui.textInput_multiplicity.value = "3"
    assert asyncio.run(ui.read_multiplicity()) is True and ui.multiplicity == 3


def test_progress_bar_splits_the_two_regions():
    ui = SharedWallBoxUI.__new__(SharedWallBoxUI)
    SharedWallBox.__init__(ui)
    ui.lattice_vectors = [10.0, 10.0, 20.0]
    ui.progress_bar, ui.status_label = Field(0), Field()

    class View:
        def submit(self, snapshot, force=False):
            self.last = snapshot

    ui.construction_view = View()
    ui.preoptimize = False
    ui._apply_progress({"stage": "insertion", "label": "Packing Box 2", "attempts": 10, "max_attempts": 100,
                        "accepted": 5, "target": 10, "base": 50.0, "span": 50.0, "molecules": []})
    assert ui.progress_bar.value == pytest.approx(75.0)
    assert ui.status_label.text.startswith("Packing Box 2")
