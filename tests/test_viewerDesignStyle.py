"""The 3D viewer's "Design style" (Molecular Design Canvas look) and its label/style sync."""
from types import SimpleNamespace

import numpy as np
import pytest

import molecularViewer as mv
from molecularDesign import ATOM_CANVAS_COLORS, ATOM_LABEL_COLOR, DEFAULT_ATOM_COLOR
from molecularViewer import MolecularViewer, MolecularViewerUI

WATER = [("O", (0.0, 0.0, 0.0)), ("H", (0.96, 0.0, 0.0)), ("H", (-0.24, 0.93, 0.0))]
BONDS = [(0, 1), (0, 2)]


def _viewer(style=mv.DESIGN_STYLE):
    v = MolecularViewer()
    v.atom_display_style = style
    return v


@pytest.fixture
def gl(monkeypatch):
    """Record the GL calls render_frame makes (no context needed)."""
    calls = []
    names = ("glPushMatrix", "glPopMatrix", "glTranslatef", "glMaterialfv", "glDisable", "glEnable",
             "glLineWidth", "glColor3f", "glBegin", "glEnd", "glVertex3f", "glCullFace")
    for name in names:
        monkeypatch.setattr(mv, name, lambda *a, _n=name: calls.append((_n, a)))
    identity = np.eye(4)
    monkeypatch.setattr(mv, "glGetDoublev", lambda what: identity)
    monkeypatch.setattr(mv, "glGetIntegerv", lambda what: np.array([0, 0, 800, 600]))
    return calls


def _record(viewer):
    spheres, labels = [], []
    viewer.draw_sphere = lambda r: spheres.append(round(r, 6))
    viewer._draw_text_window = lambda pos, text, color, font, centred, view: labels.append((text, tuple(color), centred))
    return spheres, labels


class TestHelpers:
    def test_atom_styles(self):
        assert mv.ATOM_STYLES == ("Line style", "CPK style", "vdW style", "Design style")

    def test_hex_to_rgb(self):
        assert mv.hex_to_rgb("#ff8000") == pytest.approx((1.0, 128 / 255, 0.0))
        assert mv.hex_to_rgb("33dddd") == pytest.approx((0.2, 221 / 255, 221 / 255))
        with pytest.raises(ValueError, match="#rrggbb"):
            mv.hex_to_rgb("red")

    def test_rgb_to_hex_round_trip(self):
        assert mv.rgb_to_hex(mv.hex_to_rgb("#3a3a3a")) == "#3a3a3a"
        assert mv.rgb_to_hex((0.0, 0.0, 0.0)) == "#000000"

    def test_design_colours_are_the_canvas_colours(self):
        assert mv.design_atom_colour("C") == pytest.approx(mv.hex_to_rgb(ATOM_CANVAS_COLORS["C"]))
        assert mv.design_atom_colour("Xe") == pytest.approx(mv.hex_to_rgb(DEFAULT_ATOM_COLOR))

    def test_design_palette_follows_the_background(self):
        light = mv.design_palette((1.0, 1.0, 1.0))
        dark = mv.design_palette((0.0, 0.0, 0.0))
        assert sum(light["bond"]) < 1.0 < sum(dark["bond"])
        assert sum(light["number"]) < sum(dark["number"])
        assert light["symbol"] == dark["symbol"] == pytest.approx(mv.hex_to_rgb(ATOM_LABEL_COLOR))

    def test_screen_axes_identity_and_rotated(self):
        right, up, toward = mv.screen_axes(np.eye(4))
        assert np.allclose(right, [1, 0, 0]) and np.allclose(up, [0, 1, 0]) and np.allclose(toward, [0, 0, 1])
        # 90 degrees about z in OpenGL's column-major layout: object x appears as screen y
        rot = np.eye(4)
        rot[:2, :2] = [[0.0, 1.0], [-1.0, 0.0]]   # glGetDoublev layout (transpose of the math matrix)
        right, up, _ = mv.screen_axes(rot)
        assert np.allclose(right, [0, -1, 0]) and np.allclose(up, [1, 0, 0])

    def test_design_radius_is_the_same_for_every_element(self):
        v = _viewer()
        assert v.get_atom_radius("H") == v.get_atom_radius("O") == mv.DESIGN_ATOM_RADIUS
        assert not v.show_atoms_as_spheres() and not v.show_bonds_as_lines() and not v.show_bonds_as_cylinders()


class TestDesignRendering:
    def test_discs_outlines_and_single_bond_lines(self, gl):
        v = _viewer()
        v.background_color = (1.0, 1.0, 1.0)
        spheres, labels = _record(v)
        v.render_frame(WATER, BONDS)
        r = mv.DESIGN_ATOM_RADIUS
        assert spheres == [round(r * mv.DESIGN_OUTLINE_SCALE, 6), round(r, 6)] * 3
        begins = [a for n, a in gl if n == "glBegin"]
        assert begins == [(mv.GL_LINES,)]                           # all bonds in one batch
        assert sum(1 for n, _ in gl if n == "glVertex3f") == 2 * len(BONDS)
        colours = [a for n, a in gl if n == "glColor3f"]
        assert colours[0] == pytest.approx(mv.design_palette((1.0, 1.0, 1.0))["bond"])
        assert pytest.approx(mv.design_atom_colour("O")) in colours
        assert labels == []                                        # both label switches off

    def test_symbols_centred_and_numbers_beside(self, gl):
        v = _viewer()
        v.background_color = (0.0, 0.0, 0.0)
        v.show_atom_symbols = v.show_atom_numbers = True
        _, labels = _record(v)
        v.render_frame(WATER, BONDS)
        palette = mv.design_palette((0.0, 0.0, 0.0))
        symbols = [(t, c) for t, c, centred in labels if centred]
        numbers = [(t, c) for t, c, centred in labels if not centred]
        assert [t for t, _ in symbols] == ["O", "H", "H"] and [t for t, _ in numbers] == ["1", "2", "3"]
        assert all(c == pytest.approx(palette["symbol"]) for _, c in symbols)
        assert all(c == pytest.approx(palette["number"]) for _, c in numbers)

    def test_either_label_alone(self, gl):
        v = _viewer()
        v.show_atom_numbers = True
        _, labels = _record(v)
        v.render_frame(WATER, BONDS)
        assert [t for t, _, _ in labels] == ["1", "2", "3"]

    def test_label_positions_sit_in_front_of_and_beside_the_disc(self):
        v = _viewer()
        right, up, toward = np.eye(3)
        symbol, number = v.design_label_positions((1.0, 2.0, 3.0), 0.3, right, up, toward)
        assert symbol[2] > 3.0 + 0.3 and symbol[:2] == pytest.approx((1.0, 2.0))
        assert number[0] > 1.3 and number[1] > 2.3

    def test_other_styles_unchanged(self, gl):
        v = _viewer("Line style")
        v.show_atom_symbols = v.show_atom_numbers = True
        texts = []
        v._draw_text_3d = lambda pos, text, color=(1, 1, 0): texts.append(text)
        v.draw_bond_line = lambda *a: None
        v.render_frame(WATER, BONDS)
        assert texts == ["O-1", "H-2", "H-3"]


class TestViewSettingsSync:
    def _ui(self):
        ui = MolecularViewerUI.__new__(MolecularViewerUI)
        MolecularViewer.__init__(ui)
        ui.atom_display_style_selection = SimpleNamespace(value="Line style")
        ui.atom_numbers_switch = SimpleNamespace(value=False)
        ui.atom_symbols_switch = SimpleNamespace(value=False)
        ui.qtaim_controller = SimpleNamespace(synced=0)
        ui.qtaim_controller.sync_canvas_widgets = lambda: setattr(ui.qtaim_controller, "synced",
                                                                  ui.qtaim_controller.synced + 1)
        ui.set_status_message = lambda m: None
        return ui

    def test_apply_view_settings_updates_state_and_both_tabs(self):
        ui = self._ui()
        ui.apply_view_settings(style=mv.DESIGN_STYLE, symbols=True, numbers=True)
        assert ui.atom_display_style == mv.DESIGN_STYLE and ui.show_atom_symbols and ui.show_atom_numbers
        assert ui.atom_display_style_selection.value == mv.DESIGN_STYLE
        assert ui.atom_symbols_switch.value is True and ui.atom_numbers_switch.value is True
        assert ui.qtaim_controller.synced == 1

    def test_atom_scale_updates_the_viewer_and_the_display_field(self):
        ui = self._ui()
        ui.atom_scale_input = SimpleNamespace(value="1.0")
        ui.apply_view_settings(style=mv.DESIGN_STYLE)
        assert ui.atom_scale_input.value == "1.0"   # untouched unless the scale changes
        ui.apply_view_settings(atom_scale=1.44)
        assert ui.atom_scale_factor == 1.44 and ui.atom_scale_input.value == "1.44"
        assert ui.atom_display_style == mv.DESIGN_STYLE

    @pytest.mark.parametrize("bad", [0.0, -1.0, float("nan"), float("inf")])
    def test_bad_atom_scale_is_rejected(self, bad):
        ui = self._ui()
        with pytest.raises(ValueError, match="atom size"):
            ui.apply_view_settings(atom_scale=bad)
        assert ui.atom_scale_factor == 1.0

    @pytest.mark.parametrize("style", ["CPK style", mv.DESIGN_STYLE])
    def test_atom_scale_resizes_the_drawn_atoms(self, gl, style):
        v = _viewer(style)
        spheres, _ = _record(v)
        v.draw_bond_line = v.create_cylinder = lambda *a, **k: None
        v.render_frame(WATER, BONDS)
        normal = list(spheres)
        spheres.clear()
        v.atom_scale_factor = 1.5
        v.render_frame(WATER, BONDS)
        assert normal and spheres == pytest.approx([r * 1.5 for r in normal], abs=2e-6)

    def test_unknown_style_is_rejected(self):
        with pytest.raises(ValueError, match="atom style"):
            self._ui().apply_view_settings(style="Comic style")

    def test_display_tab_handlers_route_through_the_sync(self):
        ui = self._ui()
        ui.atom_symbols_switch.value = True
        ui.toggle_atom_symbols(None)
        ui.atom_display_style_selection.value = mv.DESIGN_STYLE
        ui.set_atom_display_style(None)
        assert ui.show_atom_symbols and ui.atom_display_style == mv.DESIGN_STYLE
        assert ui.qtaim_controller.synced == 2

    def test_display_tab_offers_the_design_style(self):
        import inspect
        assert "items=list(ATOM_STYLES)" in inspect.getsource(MolecularViewerUI._build_display_tab)
