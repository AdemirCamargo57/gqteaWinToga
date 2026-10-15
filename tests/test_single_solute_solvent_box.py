import asyncio

import numpy as np
import pytest

import periodicPreOptimizer as ppo
import solventPacking
from single_solute_solvent_box import SingleSoluteSolventBox, SingleSoluteSolventBoxUI

WATER = [["O", 0.0, 0.0, 0.1173], ["H", 0.0, 0.7572, -0.4692], ["H", 0.0, -0.7572, -0.4692]]
CO2 = [["C", 0.0, 0.0, 0.0], ["O", 0.0, 0.0, 1.16], ["O", 0.0, 0.0, -1.16]]

needs_numba = pytest.mark.skipif(not solventPacking.numba_available(), reason="numba not installed")


def make_builder(length=12.0, density=0.997, attempts=3000, seed=11, solute=None, rotate=True,
                 periodic=True, padding=0.0, clearance=0.0, engine="numpy", vdw_scale=0.8):
    b = SingleSoluteSolventBox()
    b.lattice_vectors = list(length) if isinstance(length, (list, tuple)) else [length] * 3
    b.target_density = density
    b.max_attempts = attempts
    b.vdw_scale = vdw_scale
    b.random_seed = seed
    b.spacing = padding
    b.min_distance = clearance
    b.rotate_solvent = rotate
    b.use_periodic_minimum_image = periodic
    b.packing_engine = engine
    b.solvent_coords = [atom[:] for atom in WATER]
    if solute is not None:
        b.include_solute = True
        b.solute_coords = [atom[:] for atom in solute]
    b.center_of_mass()
    return b


def run(engine, **kwargs):
    b = make_builder(engine=engine, **kwargs)
    b.build_box()
    return b


CASES = [
    dict(),
    dict(seed=3, rotate=False),
    dict(seed=5, periodic=False, padding=0.4),
    dict(seed=8, solute=CO2, clearance=0.7),
    dict(seed=13, length=(11.0, 13.0, 15.0), solute=CO2),
    # Large enough for a real cell list (> 3 cells per axis), periodic and not.
    dict(seed=21, length=30.0, attempts=2500),
    dict(seed=22, length=30.0, attempts=2500, periodic=False, solute=CO2),
]


@pytest.mark.parametrize("case", CASES)
def test_numpy_engine_reproduces_the_legacy_box_exactly(case):
    legacy = run("legacy", **case)
    fast = run("numpy", **case)
    assert legacy.solvent_insertion_count > 5
    assert fast.solvent_insertion_count == legacy.solvent_insertion_count
    assert fast.solvent_attempts_used == legacy.solvent_attempts_used
    assert [a[0] for a in fast.data] == [a[0] for a in legacy.data]
    # Bit-identical coordinates, not merely close.
    assert np.array_equal(np.array([a[1:] for a in fast.data], dtype=float),
                          np.array([a[1:] for a in legacy.data], dtype=float))


@needs_numba
@pytest.mark.parametrize("case", CASES)
def test_numba_engine_reproduces_the_legacy_box_exactly(case):
    legacy = run("legacy", **case)
    fast = run("numba", **case)
    assert fast.solvent_insertion_count == legacy.solvent_insertion_count
    assert np.array_equal(np.array([a[1:] for a in fast.data], dtype=float),
                          np.array([a[1:] for a in legacy.data], dtype=float))


def test_molecules_match_data_order():
    b = run("numpy", solute=CO2)
    flat = np.concatenate([m[2] for m in b.molecules])
    assert np.array_equal(flat, np.array([a[1:] for a in b.data], dtype=float))
    assert b.molecules[0][0] == "solute"
    assert all(m[0] == "solvent" for m in b.molecules[1:])
    assert len(b.molecules) == 1 + b.solvent_insertion_count


def test_build_reports_progress_and_can_be_cancelled():
    b = make_builder(attempts=20000, density=1.2)
    payloads = []

    def progress(payload):
        payloads.append(payload)
        if payload["accepted"] >= 3:
            b.cancel_requested = True

    b.build_box(progress)
    assert payloads[0]["stage"] == "insertion"
    assert b.solvent_attempts_used < 20000


def test_auto_engine_resolves():
    b = run("auto")
    assert b.engine_used in ("numpy", "numba")


# ----------------------------------------------------------------------
# Pre-optimization inside the builder
# ----------------------------------------------------------------------
@pytest.mark.parametrize("mode", [ppo.MODE_RIGID, ppo.MODE_FLEXIBLE])
def test_preoptimization_keeps_molecules_whole_and_inside(mode):
    # Packing at 0.6 leaves contacts the 0.9 relaxation must remove.
    b = run("numpy", length=14.0, solute=CO2, attempts=4000, vdw_scale=0.6)
    b.preopt_contact_scale = 0.9
    b.preoptimize = True
    b.preopt_mode = mode
    b.run_preoptimization()
    result = b.preopt_result
    assert result is not None, b.preopt_error
    assert result.report_initial["clashes"] > 0
    assert result.report_final["clashes"] < result.report_initial["clashes"]
    assert result.connectivity_preserved is True
    system = b.packed_system()
    final = np.array([a[1:] for a in b.preopt_data])
    centres = ppo.molecule_centres(final, system.masses, system.molecule_slices)
    assert np.all(centres >= 0.0) and np.all(centres < 14.0)
    if mode == ppo.MODE_RIGID:
        assert result.max_internal_deviation < 1e-10
    assert [a[0] for a in b.preopt_data] == [a[0] for a in b.data]


def test_default_contact_scale_gives_the_preoptimization_work_to_do():
    # Packing guarantees d >= vdw_scale*(ri+rj); relaxing to the same scale would be a no-op.
    b = run("numpy", length=14.0, attempts=4000)
    assert b.preopt_contact_scale is None
    assert b.effective_preopt_contact_scale() == pytest.approx(0.90)
    b.run_preoptimization()
    assert b.preopt_result.report_initial["clashes"] > 0
    assert b.preopt_result.report_final["clashes"] == 0
    b.vdw_scale = 0.95
    assert b.effective_preopt_contact_scale() == pytest.approx(0.95)


def test_preoptimization_problem_is_recorded_not_raised():
    b = run("numpy", length=12.0)
    b.lattice_vectors = [4.0, 4.0, 4.0]  # contact range no longer fits half the box
    b.run_preoptimization()
    assert b.preopt_result is None
    assert "half the" in b.preopt_error


# ----------------------------------------------------------------------
# Properties and output
# ----------------------------------------------------------------------
def test_net_charge_and_electrons():
    b = run("numpy", solute=CO2)
    b.solute_charge = -1
    b.solvent_charge = 0
    b.multiplicity = 2
    props = b.compute_system_properties()
    n = b.solvent_insertion_count
    assert props.atom_count == 3 + 3 * n
    assert props.electron_count == 22 + 10 * n + 1
    assert props.multiplicity_valid


def test_solvent_charge_scales_with_inserted_molecules():
    b = run("numpy")
    b.solvent_charge = 1
    b.multiplicity = 1
    assert b.net_charge() == b.solvent_insertion_count
    assert b.compute_system_properties().electron_count == 9 * b.solvent_insertion_count


def test_outputs_keep_original_files_and_add_preopt(tmp_path):
    b = run("numpy", solute=CO2)
    b.output_dir = str(tmp_path)
    b.multiplicity = 1
    b.preoptimize = True
    b.run_preoptimization()
    b.compute_system_properties()
    summary = b.write_outputs()

    names = sorted(p.name for p in tmp_path.iterdir())
    assert names == sorted([
        "single_solute_solvent_box.xyz", "single_solute_solvent_box.txt",
        "single_solute_solvent_box_preopt.xyz", "single_solute_solvent.txt",
        "single_solute_box_cmass.xyz",
    ])
    lines = (tmp_path / "single_solute_solvent_box.xyz").read_text().splitlines()
    assert int(lines[0]) == len(b.data)
    assert lines[1] == "Solvent box generated by random-insertion solvent box builder"
    assert lines[2] == f"{'C':<3s}{b.data[0][1]:>14.7f}{b.data[0][2]:>14.7f}{b.data[0][3]:>14.7f}"
    pre = (tmp_path / "single_solute_solvent_box_preopt.xyz").read_text().splitlines()
    assert int(pre[0]) == len(b.data) and "cell a b c = 12.000000" in pre[1]
    assert "Total number of electrons" in summary and "--- Pre-optimization ---" in summary
    assert (tmp_path / "single_solute_solvent_box.txt").read_text(encoding="utf-8") == summary
    assert "Total number of atoms" in (tmp_path / "single_solute_solvent.txt").read_text()


def test_summary_reports_skipped_preoptimization():
    b = run("numpy")
    b.preoptimize = True
    b.preopt_error = "box too small"
    assert "Pre-optimization not performed: box too small" in b.build_summary_text()


# ----------------------------------------------------------------------
# UI plumbing, headless
# ----------------------------------------------------------------------
class Field:
    def __init__(self, value=""):
        self.value = value
        self.enabled = True
        self.text = ""


class Window:
    def __init__(self):
        self.dialogs = []

    async def dialog(self, dialog):
        self.dialogs.append(dialog)


def stub_ui(**values):
    ui = SingleSoluteSolventBoxUI.__new__(SingleSoluteSolventBoxUI)
    SingleSoluteSolventBox.__init__(ui)
    ui.main_window = Window()
    defaults = {
        "textInput_lattice_vectors": "15 15 15", "textInput_spacing": "", "textInput_target_density": "0.997",
        "textInput_max_attempts": "", "textInput_vdw_scale": "", "textInput_seed": "4",
        "textInput_min_distance": "", "textInput_solute_charge": "", "textInput_solvent_charge": "",
        "textInput_multiplicity": "1", "textInput_preopt_iterations": "",
        "textInput_preopt_contact_scale": "",
    }
    defaults.update(values)
    for name, value in defaults.items():
        setattr(ui, name, Field(value))
    ui.switch = Field(True)
    ui.switch_minimum_image = Field(True)
    ui.switch_include_solute = Field(False)
    ui.switch_preoptimize = Field(True)
    ui.selection_preopt_mode = Field("Flexible molecules")
    ui.selection_engine = Field("NumPy cell list")
    ui.multi_line_text = Field("")
    return ui


def test_read_params_supplies_every_engine_option():
    ui = stub_ui(textInput_solvent_charge="-1", textInput_multiplicity="3", textInput_preopt_iterations="250",
                 textInput_preopt_contact_scale="0.95")
    assert asyncio.run(ui.read_params(None)) is True
    assert ui.packing_engine == "numpy"
    assert ui.preoptimize is True and ui.preopt_mode == ppo.MODE_FLEXIBLE
    assert ui.preopt_max_iterations == 250
    assert ui.preopt_contact_scale == pytest.approx(0.95)
    assert ui.solvent_charge == -1 and ui.solute_charge == 0
    assert ui.multiplicity == 3
    assert ui.rotate_solvent is True and ui.include_solute is False
    # Every keyword the pre-optimizer and the packer read is an attribute read_params set.
    for name in ["lattice_vectors", "vdw_scale", "spacing", "min_distance", "use_periodic_minimum_image",
                 "random_seed", "max_attempts", "target_density"]:
        assert getattr(ui, name) is not None


@pytest.mark.parametrize("field,value,message", [
    ("textInput_multiplicity", "", "not inferred"),
    ("textInput_multiplicity", "0", "positive"),
    ("textInput_solvent_charge", "0.5", "Solvent molecule charge"),
    ("textInput_preopt_iterations", "0", "iterations"),
    ("textInput_lattice_vectors", "10 -2 10", "positive"),
    ("textInput_seed", "abc", "random seed"),
    ("textInput_preopt_contact_scale", "-1", "contact scale"),
    ("textInput_preopt_contact_scale", "x", "contact scale"),
])
def test_read_params_rejects_bad_inputs(field, value, message, monkeypatch):
    import toga

    class Recorded:
        def __init__(self, title, message):
            self.title, self.message = title, message

    monkeypatch.setattr(toga, "ErrorDialog", Recorded)
    monkeypatch.setattr(toga, "InfoDialog", Recorded)
    ui = stub_ui(**{field: value})
    assert asyncio.run(ui.read_params(None)) is False
    assert ui.main_window.dialogs
    assert message in ui.main_window.dialogs[-1].message


class StubView:
    def __init__(self):
        self.snapshots = []

    def submit(self, snapshot, force=False):
        self.snapshots.append(snapshot)


def test_apply_progress_updates_bar_status_and_canvas():
    ui = stub_ui()
    asyncio.run(ui.read_params(None))
    ui.progress_bar = Field(0)
    ui.status_label = Field()
    ui.construction_view = StubView()
    ui._apply_progress({"stage": "insertion", "attempts": 50, "accepted": 10, "target": 20, "molecules": []})
    assert ui.progress_bar.value == pytest.approx(0.7 * 50.0)
    assert "accepted 10/20" in ui.status_label.text
    ui._apply_progress({"stage": "preopt", "iteration": 250, "max_iterations": 500, "energy": 1.5,
                        "molecules": []})
    assert ui.progress_bar.value == pytest.approx(85.0)
    assert [s.stage for s in ui.construction_view.snapshots] == ["Inserting solvent", "Pre-optimizing"]
