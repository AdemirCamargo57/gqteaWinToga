"""Headless tests for the ESP viewer: colour maths and surface state (no GL window)."""
import numpy as np
import pytest

import espViewer as ev
from molecularViewer import MolecularViewer


class TestColours:
    def test_hue_sweep_end_points(self):
        rgb = ev.hue_colours(np.array([0.0, 0.25, 0.5, 0.75, 1.0]))
        assert np.allclose(rgb, [[1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 1, 1], [0, 0, 1]])

    def test_esp_zero_is_green_and_clipped(self):
        rgb = ev.esp_colours(np.array([-0.05, 0.0, 0.05, 0.2]), 0.05)
        assert np.allclose(rgb[0], [1, 0, 0]) and np.allclose(rgb[1], [0, 1, 0])
        assert np.allclose(rgb[2], [0, 0, 1]) and np.allclose(rgb[3], [0, 0, 1])

    def test_symmetric_range_is_max_abs(self):
        assert ev.symmetric_range(np.array([-0.02097, 0.01166])) == pytest.approx(0.02097)
        assert ev.symmetric_range(np.zeros(3)) == 1.0  # never a zero-width bar

    @pytest.mark.parametrize("value, text", [(-0.02072, "-2.072e-2"), (0.02072, "2.072e-2"),
                                             (-13.157, "-1.316e1"), (0.0, "0.000e0")])
    def test_gaussview_label(self, value, text):
        assert ev.gaussview_label(value) == text

    def test_colour_bar_text_names_both_quantities(self):
        left, right, caption = ev.colour_bar_text(0.02097, ev.UNIT_KCAL, 0.0004)
        assert left == "-1.316e1" and right == "1.316e1"
        assert "ESP (kcal/mol)" in caption and "rho = 0.000400 a.u." in caption

    def test_unknown_unit(self):
        with pytest.raises(ValueError):
            ev.unit_factor("eV")


class TestViewerState:
    @pytest.fixture
    def viewer(self):
        v = ev.ESPSurfaceViewer()
        v.set_molecule(["C", "C"], np.array([[0, 0, -0.67], [0, 0, 0.67]]))
        return v

    def test_base_viewer_hooks_are_no_ops(self):
        base = MolecularViewer()
        assert base.render_scene_extras() is None and base.draw_screen_overlay() is None
        assert base.window_title == "3D Molecular Viewer"

    def test_set_surface_builds_draw_data(self, viewer):
        verts = np.array([[3.0, 0, 0], [0, 3.0, 0], [0, 0, 3.0]])
        viewer.set_surface(verts, np.array([[0, 1, 2]]), np.eye(3), np.array([-0.02, 0.0, 0.01]), 0.0004)
        s = viewer.get_surface()
        assert s.vertices.dtype == np.float32 and s.indices.dtype == np.uint32
        assert s.half_range == pytest.approx(0.02)
        assert np.allclose(s.rgba[:, 3], 1.0)
        assert viewer.get_scene_radius() >= 3.0  # the camera frames the surface, not just the atoms
        assert viewer.colour_bar_text()[0] == "-2.000e-2"

    def test_opacity_updates_alpha(self, viewer):
        viewer.set_surface(np.zeros((3, 3)) + [[1, 0, 0]], np.array([[0, 1, 2]]), np.eye(3), np.zeros(3), 0.001)
        viewer.set_display_options(opacity=0.5)
        assert np.allclose(viewer.get_surface().rgba[:, 3], 0.5)

    def test_clear_surface(self, viewer):
        viewer.set_surface(np.eye(3), np.array([[0, 1, 2]]), np.eye(3), np.zeros(3), 0.001)
        viewer.clear_surface()
        assert viewer.get_surface() is None and viewer.colour_bar_text() is None


class TestOverlayText:
    def test_text_contrasts_with_the_background(self):
        assert ev.contrasting_text_colour((0.0, 0.0, 0.0)) == (1.0, 1.0, 1.0)
        assert ev.contrasting_text_colour((1.0, 1.0, 1.0)) == (0.0, 0.0, 0.0)
        assert ev.contrasting_text_colour((0.85, 0.85, 0.85)) == (0.0, 0.0, 0.0)
        assert ev.contrasting_text_colour((0.25, 0.25, 0.25)) == (1.0, 1.0, 1.0)


def test_custom_caption_replaces_the_esp_caption():
    v = ev.ESPSurfaceViewer()
    tri = np.array([[0, 1, 2]])
    verts = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0]])
    v.set_surface(verts, tri, np.tile([0, 0, 1.0], (3, 1)), [0.01, -0.01, 0.0], 0.01, caption="f- (a.u.)")
    assert v.colour_bar_text()[2] == "f- (a.u.)"
    v.set_surface(verts, tri, np.tile([0, 0, 1.0], (3, 1)), [0.01, -0.01, 0.0], 0.01)
    assert v.colour_bar_text()[2].startswith("ESP (a.u.)")


class TestSelectionFilter:
    """Show only selected: the surface and markers follow the selection (nearest-atom owners)."""

    # Two H2 molecules 10 A apart; one triangle (3 vertices) round each atom.
    XYZ = np.array([[0, 0, 0], [0.74, 0, 0], [10, 0, 0], [10.74, 0, 0]])

    @pytest.fixture
    def viewer(self):
        v = ev.ESPSurfaceViewer()
        v.set_molecule(["H"] * 4, self.XYZ)
        offsets = np.array([[0, 0.5, 0], [0, 0, 0.5], [0, -0.5, 0]])
        verts = np.vstack([atom + offsets for atom in self.XYZ])
        tris = np.arange(12).reshape(4, 3)
        values = np.array([0.01] * 3 + [0.02] * 3 + [-0.03] * 3 + [-0.04] * 3)
        markers = [ev.SurfaceMarker(tuple(self.XYZ[0] + [0, 0.5, 0]), "max", "global", "1"),
                   ev.SurfaceMarker(tuple(self.XYZ[3] + [0, 0.5, 0]), "min", "global", "2")]
        v.set_surface(verts, tris, np.tile([0, 0, 1.0], (12, 1)), values, 0.01, markers,
                      vertex_atoms=np.repeat(np.arange(4), 3), marker_atoms=[0, 3])
        return v

    def test_without_show_only_the_whole_surface_is_drawn(self, viewer):
        viewer.set_selected_atoms([0, 1])
        s = viewer.get_surface()
        assert len(s.vertices) == 12 and len(s.indices) == 12 and len(s.markers) == 2
        assert not viewer.filtering_by_selection()

    def test_show_only_keeps_the_selected_atoms_part(self, viewer):
        viewer.set_selected_atoms([0, 1])
        viewer.set_show_selected_only(True)
        s = viewer.get_surface()
        assert len(s.vertices) == 6 and s.indices.tolist() == [0, 1, 2, 3, 4, 5]
        assert np.all(s.vertices[:, 0] < 1.0)  # only round the first molecule
        assert [m.label for m in s.markers] == ["1"]  # extrema-file row numbers kept
        assert s.half_range == pytest.approx(0.04)  # whole-system colour range
        assert viewer.filtering_by_selection() and not viewer.surface_hidden_by_selection()
        # Framed round the selection (as set_projection asks), not the whole system.
        assert viewer.get_scene_radius(viewer.camera_focus_frame()) < 3.0

    def test_selection_changes_refilter(self, viewer):
        viewer.set_show_selected_only(True)
        viewer.set_selected_atoms([3])
        assert np.all(viewer.get_surface().vertices[:, 0] > 10.5)
        viewer.toggle_selected_atom(3)  # empty selection: everything again
        assert len(viewer.get_surface().vertices) == 12
        viewer.set_selected_atoms([1])
        viewer.set_show_selected_only(False)
        assert len(viewer.get_surface().vertices) == 12

    def test_a_selection_with_no_surface_hides_it_and_says_so(self, viewer):
        viewer.set_surface(np.eye(3), np.array([[0, 1, 2]]), np.eye(3), np.zeros(3), 0.01,
                           vertex_atoms=[2, 2, 3])
        viewer.set_selected_atoms([0])
        viewer.set_show_selected_only(True)
        assert viewer.get_surface() is None and viewer.surface_hidden_by_selection()
        viewer.clear_surface()
        assert not viewer.surface_hidden_by_selection()

    def test_a_surface_without_owners_is_never_filtered(self):
        v = ev.ESPSurfaceViewer()
        v.set_molecule(["H"] * 4, self.XYZ)
        v.set_selected_atoms([0])
        v.set_show_selected_only(True)
        v.set_surface(np.eye(3), np.array([[0, 1, 2]]), np.eye(3), np.zeros(3), 0.01)
        assert len(v.get_surface().vertices) == 3

    def test_owner_arrays_must_match(self, viewer):
        with pytest.raises(ValueError, match="vertex_atoms"):
            viewer.set_surface(np.eye(3), np.array([[0, 1, 2]]), np.eye(3), np.zeros(3), 0.01,
                               vertex_atoms=[0, 1])
        with pytest.raises(ValueError, match="marker_atoms"):
            viewer.set_surface(np.eye(3), np.array([[0, 1, 2]]), np.eye(3), np.zeros(3), 0.01,
                               markers=[ev.SurfaceMarker((0, 0, 0), "max", "global", "1")],
                               vertex_atoms=[0, 1, 2], marker_atoms=[])

    def test_same_structure_keeps_the_selection(self, viewer):
        viewer.set_selected_atoms([2, 3])
        viewer.set_molecule(["H"] * 4, self.XYZ + 1e-6)
        assert viewer.get_selected_atoms() == [2, 3]
        viewer.set_molecule(["H"] * 3, self.XYZ[:3])  # another molecule: a real load
        assert viewer.get_selected_atoms() == []

    def test_selection_listener_is_called_on_the_loop(self, viewer):
        calls = []

        class Loop:
            def is_closed(self):
                return False

            def call_soon_threadsafe(self, callback, *args):
                calls.append(callback)

        viewer.loop = Loop()
        viewer.selection_listener = lambda: None
        viewer.set_selected_atoms([0])
        assert calls == [viewer.selection_listener]


def test_atom_view_sets_style_and_numbers():
    v = ev.ESPSurfaceViewer()
    v.set_atom_view("CPK style", True)
    assert v.atom_display_style == "CPK style" and v.show_atom_numbers is True
    v.set_atom_view("vdW style", False)
    assert v.atom_display_style == "vdW style" and v.show_atom_numbers is False
    with pytest.raises(ValueError, match="atom style"):
        v.set_atom_view("Cartoon", False)
