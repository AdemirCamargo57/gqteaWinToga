import asyncio

import numpy as np
import pytest

import periodicPreOptimizer as ppo
import solventPacking
from mixtureSolventBox import MixtureSolvationBox, MixtureSolventBoxUI

WATER = [["O", 0.0, 0.0, 0.1173], ["H", 0.0, 0.7572, -0.4692], ["H", 0.0, -0.7572, -0.4692]]
METHANOL = [["C", -0.0465, 0.6630, 0.0], ["O", -0.0465, -0.7570, 0.0], ["H", -1.0861, 0.9729, 0.0],
            ["H", 0.4366, 1.0816, 0.8904], ["H", 0.4366, 1.0816, -0.8904], ["H", 0.8624, -1.0857, 0.0]]
CO2 = [["C", 0.0, 0.0, 0.0], ["O", 0.0, 0.0, 1.16], ["O", 0.0, 0.0, -1.16]]

needs_numba = pytest.mark.skipif(not solventPacking.numba_available(), reason="numba not installed")


def make_builder(length=14.0, density=0.9, attempts=3000, seed=11, comp_a=0.6, solute=None, rotate=True,
                 periodic=True, padding=0.0, clearance=0.0, engine="numpy", vdw_scale=0.8):
    b = MixtureSolvationBox()
    b.lattice_vectors = list(length) if isinstance(length, (list, tuple)) else [length] * 3
    b.target_density = density
    b.max_attempts = attempts
    b.vdw_scale = vdw_scale
    b.random_seed = seed
    b.spacing = padding
    b.min_distance = clearance
    b.composition_A, b.composition_B = comp_a, 1.0 - comp_a
    b.rotate_solvent = rotate
    b.use_periodic_minimum_image = periodic
    b.packing_engine = engine
    b.coords_A = [a[:] for a in WATER]
    b.coords_B = [a[:] for a in METHANOL]
    if solute is not None:
        b.include_solute = True
        b.solute_coords = [a[:] for a in solute]
    b.prepare_molecules()
    b.estimate_target_molecule_count()
    return b


def run(engine, **kwargs):
    b = make_builder(engine=engine, **kwargs)
    b.build_box()
    return b


CASES = [
    dict(),
    dict(seed=3, comp_a=1.0),
    dict(seed=4, comp_a=0.0),
    dict(seed=5, rotate=False),
    dict(seed=6, periodic=False, padding=0.5),
    dict(seed=8, solute=CO2, clearance=0.7),
    dict(seed=21, length=30.0, attempts=2500),
    dict(seed=22, length=(28.0, 31.0, 33.0), attempts=2500, periodic=False, solute=CO2),
]


def coords_of(b):
    return np.array([a[1:] for a in b.data], dtype=float)


@pytest.mark.parametrize("case", CASES)
def test_numpy_engine_reproduces_the_legacy_box_exactly(case):
    legacy = run("legacy", **case)
    fast = run("numpy", **case)
    assert legacy.count_solvent_A + legacy.count_solvent_B > 5
    assert (fast.count_solvent_A, fast.count_solvent_B) == (legacy.count_solvent_A, legacy.count_solvent_B)
    assert fast.solvent_attempts_used == legacy.solvent_attempts_used
    assert [a[0] for a in fast.data] == [a[0] for a in legacy.data]
    assert np.array_equal(coords_of(fast), coords_of(legacy))


@needs_numba
@pytest.mark.parametrize("case", CASES)
def test_numba_engine_reproduces_the_legacy_box_exactly(case):
    legacy = run("legacy", **case)
    fast = run("numba", **case)
    assert np.array_equal(coords_of(fast), coords_of(legacy))


def test_unseeded_runs_still_complete():
    b = run("numpy", seed=None)
    assert b.count_solvent_A + b.count_solvent_B > 0


def test_molecules_carry_their_species_and_match_data_order():
    b = run("numpy", solute=CO2)
    kinds = [m[0] for m in b.molecules]
    assert kinds[0] == "solute"
    assert kinds.count("solvent_A") == b.count_solvent_A
    assert kinds.count("solvent_B") == b.count_solvent_B
    flat = np.concatenate([m[2] for m in b.molecules])
    assert np.array_equal(flat, coords_of(b))


def test_box_too_small_for_the_molecules_stops_cleanly():
    # Mixture bounding radii exclude vdW: water (~0.96 A) cannot be centred in a 1.5 A box.
    b = run("numpy", length=1.5, density=30.0, attempts=100)
    assert b.num_target_solvent_molecules > 0
    assert b.count_solvent_A + b.count_solvent_B == 0
    assert b.solvent_attempts_used == 1


@pytest.mark.parametrize("mode", [ppo.MODE_RIGID, ppo.MODE_FLEXIBLE])
def test_preoptimization_of_a_mixture(mode):
    b = run("numpy", length=16.0, solute=CO2, attempts=4000)
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
    assert np.all(centres >= 0.0) and np.all(centres < 16.0)
    if mode == ppo.MODE_RIGID:
        assert result.max_internal_deviation < 1e-10


def test_net_charge_counts_each_species():
    b = run("numpy", solute=CO2)
    b.solute_charge, b.solvent_A_charge, b.solvent_B_charge = 2, -1, 1
    assert b.net_charge() == 2 - b.count_solvent_A + b.count_solvent_B
    b.multiplicity = 1
    props = b.compute_system_properties()
    electrons = 22 + 10 * b.count_solvent_A + 18 * b.count_solvent_B - b.net_charge()
    assert props.electron_count == electrons


def test_outputs(tmp_path):
    b = run("numpy", solute=CO2)
    b.output_dir = str(tmp_path)
    b.multiplicity = 1
    b.preoptimize = True
    b.run_preoptimization()
    b.compute_system_properties()
    summary = b.write_outputs()
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "mixture_solvent_box.txt", "mixture_solvent_box.xyz", "mixture_solvent_box_preopt.xyz"]
    lines = (tmp_path / "mixture_solvent_box.xyz").read_text(encoding="utf-8").splitlines()
    assert int(lines[0]) == len(b.data)
    assert lines[1] == "Mixture solvent box generated by density-driven random insertion"
    assert "Inserted solvent A molecules" in summary and "Total number of electrons" in summary
    assert "--- Pre-optimization ---" in summary


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
    ui = MixtureSolventBoxUI.__new__(MixtureSolventBoxUI)
    MixtureSolvationBox.__init__(ui)
    ui.main_window = Window()
    defaults = {
        "textInput_lattice_vectors": "20 20 20", "textInput_spacing": "", "textInput_target_density": "0.9",
        "textInput_max_attempts": "", "textInput_vdw_scale": "", "textInput_seed": "",
        "textInput_composition_A": "70", "textInput_composition_B": "30", "textInput_min_distance": "",
        "textInput_solute_charge": "", "textInput_solvent_A_charge": "", "textInput_solvent_B_charge": "",
        "textInput_multiplicity": "1", "textInput_preopt_iterations": "", "textInput_preopt_contact_scale": "",
    }
    defaults.update(values)
    for name, value in defaults.items():
        setattr(ui, name, Field(value))
    ui.switch_rotate = Field(False)
    ui.switch_minimum_image = Field(True)
    ui.switch_include_solute = Field(True)
    ui.switch_preoptimize = Field(True)
    ui.selection_preopt_mode = Field("Rigid molecules")
    ui.selection_engine = Field("Legacy (reference, slow)")
    ui.multi_line_text = Field("")
    return ui


def test_read_params_supplies_every_option():
    ui = stub_ui(textInput_solvent_A_charge="1", textInput_solvent_B_charge="-1", textInput_solute_charge="2",
                 textInput_multiplicity="2")
    assert asyncio.run(ui.read_params(None)) is True
    assert (ui.solute_charge, ui.solvent_A_charge, ui.solvent_B_charge) == (2, 1, -1)
    assert ui.multiplicity == 2
    assert ui.rotate_solvent is False and ui.include_solute is True
    assert ui.packing_engine == "legacy" and ui.preoptimize is True
    assert ui.composition_A == pytest.approx(0.7)


@pytest.mark.parametrize("field,value,message", [
    ("textInput_multiplicity", "", "not inferred"),
    ("textInput_solvent_B_charge", "x", "Solvent B molecule charge"),
    ("textInput_lattice_vectors", "20 0 20", "positive"),
    ("textInput_seed", "1.5", "random seed"),
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
    assert message in ui.main_window.dialogs[-1].message
