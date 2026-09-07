"""Tests for the molecularViewer atom-picking geometry and selection rules.

These exercise the pure projection / nearest-atom math used by the
click-to-identify feature, plus the ordered pick list that feeds the Measure
tab. They construct known OpenGL-style matrices and feed them directly, so no
GL context (and no GUI) is required — matching the headless test convention in
tests/conftest.py.
"""
import numpy as np
import pytest

from molecularViewer import MolecularViewer, MolecularViewerUI


def _ortho_matrix(left, right, bottom, top, near, far):
    """glOrtho matrix in the column-major layout glGetDoublev returns.

    Stored so that ``obj_vec @ M`` reproduces the OpenGL transform.
    """
    return np.array(
        [
            [2.0 / (right - left), 0.0, 0.0, 0.0],
            [0.0, 2.0 / (top - bottom), 0.0, 0.0],
            [0.0, 0.0, -2.0 / (far - near), 0.0],
            [
                -(right + left) / (right - left),
                -(top + bottom) / (top - bottom),
                -(far + near) / (far - near),
                1.0,
            ],
        ],
        dtype=float,
    )


def _viewer_with_view():
    viewer = MolecularViewer()
    viewer.frames = []
    viewer.molecule_data = [
        ("C", (5.0, 0.0, 0.0)),
        ("O", (-5.0, 0.0, 0.0)),
        ("N", (0.0, 0.0, 0.0)),
    ]
    identity = np.eye(4)
    proj = _ortho_matrix(-10, 10, -10, 10, -10, 10)
    viewer._pick_view = (identity, proj, (0, 0, 200, 200))
    return viewer


def test_project_atom_maps_to_expected_pixels():
    viewer = _viewer_with_view()
    identity = np.eye(4)
    proj = _ortho_matrix(-10, 10, -10, 10, -10, 10)

    sx, sy, sz = viewer._project_atom((5.0, 0.0, 0.0), identity, proj, (0, 0, 200, 200))

    assert sx == pytest.approx(150.0)
    assert sy == pytest.approx(100.0)
    assert 0.0 <= sz <= 1.0


def test_pick_selects_atom_under_cursor():
    viewer = _viewer_with_view()

    # Cursor uses GLFW top-left origin; window height is 200.
    # C -> window (150, 100) -> cursor (150, 100)
    assert viewer._pick_atom_at(150, 100) == 0
    # O -> window (50, 100) -> cursor (50, 100)
    assert viewer._pick_atom_at(50, 100) == 1
    # N -> window (100, 100) -> cursor (100, 100)
    assert viewer._pick_atom_at(100, 100) == 2


def test_pick_empty_space_returns_none():
    viewer = _viewer_with_view()

    assert viewer._pick_atom_at(10, 10) is None


def test_pick_returns_none_without_captured_view():
    viewer = MolecularViewer()
    viewer._pick_view = None

    assert viewer._pick_atom_at(100, 100) is None


def test_pick_prefers_nearer_atom_on_overlap():
    viewer = MolecularViewer()
    viewer.frames = []
    # Both atoms project to the same pixel but differ in depth.
    viewer.molecule_data = [
        ("C", (0.0, 0.0, -5.0)),
        ("O", (0.0, 0.0, 5.0)),
    ]
    identity = np.eye(4)
    proj = _ortho_matrix(-10, 10, -10, 10, -10, 10)
    viewer._pick_view = (identity, proj, (0, 0, 200, 200))

    # O (z=5) maps to the smaller window-z (nearer), so it wins.
    assert viewer._pick_atom_at(100, 100) == 1


# ----------------------------------------------------------------------
# Ordered canvas selection feeding the Measure tab
# ----------------------------------------------------------------------
def _viewer_with_atoms(count=6):
    viewer = MolecularViewer()
    viewer.frames = []
    viewer.molecule_data = [("C", (float(i), 0.0, 0.0)) for i in range(count)]
    return viewer


def test_measurement_atom_count_per_measure_type():
    viewer = MolecularViewer()

    assert viewer.measurement_atom_count("Bond length") == 2
    assert viewer.measurement_atom_count("Bond angle") == 3
    assert viewer.measurement_atom_count("Dihedral angle") == 4
    assert viewer.measurement_atom_count("Atom coordinates") == 1
    assert viewer.measurement_atom_count("Nonsense") is None


def test_picks_append_in_click_order():
    viewer = _viewer_with_atoms()
    viewer.set_measurement_pick_capacity(3)

    viewer._toggle_picked_atom(4)
    viewer._toggle_picked_atom(0)
    viewer._toggle_picked_atom(2)

    # Click order is the measurement order: the 2nd click is the angle vertex.
    assert viewer.get_picked_atoms() == [4, 0, 2]


def test_clicking_a_selected_atom_removes_it():
    viewer = _viewer_with_atoms()
    viewer.set_measurement_pick_capacity(3)

    viewer._toggle_picked_atom(1)
    viewer._toggle_picked_atom(2)
    viewer._toggle_picked_atom(1)

    assert viewer.get_picked_atoms() == [2]


def test_picks_roll_the_oldest_out_when_full():
    viewer = _viewer_with_atoms()
    viewer.set_measurement_pick_capacity(2)

    viewer._toggle_picked_atom(0)
    viewer._toggle_picked_atom(1)
    viewer._toggle_picked_atom(2)

    assert viewer.get_picked_atoms() == [1, 2]

    viewer._toggle_picked_atom(3)

    assert viewer.get_picked_atoms() == [2, 3]


def test_shrinking_capacity_keeps_the_most_recent_picks():
    viewer = _viewer_with_atoms()
    viewer.set_measurement_pick_capacity(4)
    for index in (0, 1, 2, 3):
        viewer._toggle_picked_atom(index)

    viewer.set_measurement_pick_capacity(2)

    assert viewer.get_picked_atoms() == [2, 3]


def test_growing_capacity_keeps_existing_picks():
    viewer = _viewer_with_atoms()
    viewer.set_measurement_pick_capacity(2)
    viewer._toggle_picked_atom(0)
    viewer._toggle_picked_atom(1)

    viewer.set_measurement_pick_capacity(4)

    assert viewer.get_picked_atoms() == [0, 1]


def test_picked_atoms_text_is_one_based_and_ordered():
    viewer = _viewer_with_atoms()
    viewer.set_measurement_pick_capacity(3)
    viewer._toggle_picked_atom(4)
    viewer._toggle_picked_atom(0)

    assert viewer.picked_atoms_text() == "5,1"


def test_picked_atoms_text_is_empty_without_picks():
    viewer = _viewer_with_atoms()

    assert viewer.picked_atoms_text() == ""


def test_clear_picked_atoms_reports_whether_anything_was_cleared():
    viewer = _viewer_with_atoms()
    viewer.set_measurement_pick_capacity(2)

    assert viewer.clear_picked_atoms() is False

    viewer._toggle_picked_atom(0)

    assert viewer.clear_picked_atoms() is True
    assert viewer.get_picked_atoms() == []


def test_pick_changes_notify_the_ui_hook():
    viewer = _viewer_with_atoms()
    viewer.set_measurement_pick_capacity(2)
    calls = []
    viewer._notify_picked_atoms_changed = lambda: calls.append(viewer.get_picked_atoms())

    viewer._toggle_picked_atom(1)
    viewer._toggle_picked_atom(1)
    viewer.clear_picked_atoms()

    assert calls == [[1], [], []]


# ----------------------------------------------------------------------
# Measure-tab plumbing (stub widgets; no Toga window is created)
# ----------------------------------------------------------------------
class _StubInput:
    def __init__(self):
        self.value = ""


class _StubLoop:
    """Just enough event loop for the pick -> field -> measure hand-off."""

    def __init__(self):
        self.measurements = 0

    def is_closed(self):
        return False

    def call_soon_threadsafe(self, callback, *args):
        callback(*args)

    def create_task(self, coro):
        coro.close()  # never awaited here; we only assert it was scheduled
        self.measurements += 1


def _headless_ui(measure_type="Bond length", atoms=6):
    ui = MolecularViewerUI.__new__(MolecularViewerUI)
    MolecularViewer.__init__(ui)
    ui.loop = _StubLoop()
    ui.loading_label = None
    ui.measure_indices_input = _StubInput()
    ui.measurement_label = _StubInput()
    ui.frames = []
    ui.molecule_data = [("C", (float(i), 0.0, 0.0)) for i in range(atoms)]
    ui.set_measurement_pick_capacity(ui.measurement_atom_count(measure_type))
    return ui


def test_clicking_atoms_fills_the_measure_field():
    ui = _headless_ui("Bond angle")

    ui._toggle_picked_atom(2)
    assert ui.measure_indices_input.value == "3"
    ui._toggle_picked_atom(0)

    assert ui.measure_indices_input.value == "3,1"


def test_measurement_runs_once_the_selection_is_complete():
    ui = _headless_ui("Bond length")

    ui._toggle_picked_atom(0)
    assert ui.loop.measurements == 0

    ui._toggle_picked_atom(1)

    assert ui.measure_indices_input.value == "1,2"
    assert ui.loop.measurements == 1


def test_rolling_pick_remeasures_with_the_new_pair():
    ui = _headless_ui("Bond length")
    ui._toggle_picked_atom(0)
    ui._toggle_picked_atom(1)

    ui._toggle_picked_atom(4)

    assert ui.measure_indices_input.value == "2,5"
    assert ui.loop.measurements == 2


def test_clearing_the_selection_empties_the_field():
    ui = _headless_ui("Bond length")
    ui._toggle_picked_atom(0)
    ui._toggle_picked_atom(1)

    ui.clear_measure_selection(None)

    assert ui.measure_indices_input.value == ""
    assert ui.get_picked_atoms() == []


# ----------------------------------------------------------------------
# Canvas overlay text: the bare number, no label / prefix / unit
# ----------------------------------------------------------------------
def test_bond_length_overlay_is_the_bare_number():
    viewer = MolecularViewer()
    viewer.frames = []
    viewer.molecule_data = [("C", (0.0, 0.0, 0.0)), ("O", (1.5, 0.0, 0.0))]

    value, overlay = viewer.compute_measurement_value("Bond length", [0, 1])

    assert value == pytest.approx(1.5)
    assert overlay == "1.5000"


def test_bond_angle_overlay_is_the_bare_number():
    viewer = MolecularViewer()
    viewer.frames = []
    viewer.molecule_data = [
        ("H", (1.0, 0.0, 0.0)),
        ("O", (0.0, 0.0, 0.0)),
        ("H", (0.0, 1.0, 0.0)),
    ]

    value, overlay = viewer.compute_measurement_value("Bond angle", [0, 1, 2])

    assert value == pytest.approx(90.0)
    assert overlay == "90.000"


def test_dihedral_overlay_is_the_bare_number():
    viewer = MolecularViewer()
    viewer.frames = []
    # Planar arrangement, so the dihedral is exactly 0 or 180 degrees.
    viewer.molecule_data = [
        ("C", (1.0, 1.0, 0.0)),
        ("C", (0.0, 0.0, 0.0)),
        ("C", (1.0, 0.0, 0.0)),
        ("C", (2.0, 1.0, 0.0)),
    ]

    value, overlay = viewer.compute_measurement_value("Dihedral angle", [0, 1, 2, 3])

    assert overlay == f"{value:.3f}"
    assert not any(ch.isalpha() for ch in overlay)


def test_result_field_keeps_its_full_label():
    viewer = MolecularViewer()
    viewer.frames = []
    viewer.molecule_data = [("C", (0.0, 0.0, 0.0)), ("O", (1.5, 0.0, 0.0))]

    text = viewer.format_measurement_label("Bond length", [0, 1], 1.5)

    assert text == "Measurement: d(1,2) = 1.5000 Å"


def test_deselecting_back_to_empty_clears_the_field():
    ui = _headless_ui("Bond length")
    ui._toggle_picked_atom(3)
    assert ui.measure_indices_input.value == "4"

    ui._toggle_picked_atom(3)

    assert ui.measure_indices_input.value == ""
