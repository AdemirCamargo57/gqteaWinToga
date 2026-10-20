"""QTAIM tab of the 3D Molecular Viewer: scene layer, controller (headless) and integration."""
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import qtaimAnalysis as qa
import qtaimLayer as ql
import qtaimTopology as qt
from qtaim_fixtures import fchk

SRC = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def water_results():
    return qa.QTAIMAnalyser().run(str(fchk("water_631gd_6d.fchk")))


@pytest.fixture(scope="module")
def cage_results():
    return qa.QTAIMAnalyser().run(str(fchk("bicyclopentane_631gd.fchk")))


class RecordingViewer:
    def __init__(self):
        self.spheres, self.texts = [], []

    def draw_sphere(self, radius):
        self.spheres.append(radius)

    def _draw_text_3d(self, position, text, color=(1, 1, 0)):
        self.texts.append((text, color))


@pytest.fixture
def no_gl(monkeypatch):
    calls = []
    for name in ("glPushMatrix", "glPopMatrix", "glTranslatef", "glMaterialfv", "glDisable", "glEnable",
                 "glLineWidth", "glColor3f", "glBegin", "glEnd", "glVertex3f"):
        monkeypatch.setattr(ql, name, lambda *a, _n=name: calls.append(_n))
    return calls


class TestSceneLayer:
    def test_scene_from_results_is_in_angstrom(self, water_results):
        scene = ql.scene_from_results(water_results)
        assert len(scene.cps) == 5 and len(scene.bond_paths) == 4 and scene.ring_paths == ()
        cp = water_results.critical_points[3]
        assert np.allclose(scene.cps[3].position, cp.position * qt.BOHR_RADIUS_ANGSTROM)
        assert scene.bond_paths[0].dtype == np.float32

    def test_visible_cps_follow_the_type_switches(self, cage_results):
        layer = ql.QTAIMSceneLayer()
        layer.set_scene(ql.scene_from_results(cage_results))
        assert len(layer.visible_cps()) == 31
        layer.set_display_options(show_types={qt.NCP: False})
        assert {c.kind for c in layer.visible_cps()} == {qt.BCP, qt.RCP, qt.CCP}
        layer.set_display_options(show_types={qt.BCP: False, qt.RCP: False, qt.CCP: False})
        assert layer.visible_cps() == ()
        with pytest.raises(ValueError, match="Unknown critical point type"):
            layer.set_display_options(show_types={"(3,0)": True})

    def test_radius(self):
        layer = ql.QTAIMSceneLayer()
        assert layer.radius() == 0.0
        layer.set_scene(ql.QTAIMScene((ql.SceneCP(1, qt.NCP, (1.0, 0.0, 0.0)),),
                                      (np.array([[0, 0, 0], [0, 3.0, 0]], dtype=np.float32),), (), (0.0, 0.0, 0.0)))
        assert layer.radius() == pytest.approx(3.0)

    def test_draw_scene_spheres_highlight_labels_and_paths(self, water_results, no_gl):
        layer = ql.QTAIMSceneLayer()
        layer.set_scene(ql.scene_from_results(water_results))
        layer.set_display_options(highlight=4)
        viewer = RecordingViewer()
        layer.draw_scene(viewer)
        assert viewer.spheres.count(ql.HIGHLIGHT_RADIUS) == 1 and viewer.spheres.count(ql.CP_RADIUS) == 4
        assert [t for t, _ in viewer.texts] == ["1", "2", "3", "4", "5"]
        assert no_gl.count("glBegin") == 4 and no_gl.count("glEnd") == 4
        last_end = len(no_gl) - 1 - no_gl[::-1].index("glEnd")
        assert "glEnable" in no_gl[last_end:]   # lighting restored after the paths

    def test_draw_scene_honours_the_switches(self, cage_results, no_gl):
        layer = ql.QTAIMSceneLayer()
        layer.set_scene(ql.scene_from_results(cage_results))
        layer.set_display_options(bond_paths=False, labels=False, show_types={qt.NCP: False})
        viewer = RecordingViewer()
        layer.draw_scene(viewer)
        assert no_gl.count("glBegin") == 3      # the three ring paths only
        assert viewer.texts == [] and len(viewer.spheres) == 18

    def test_empty_layer_draws_nothing(self, no_gl):
        viewer = RecordingViewer()
        ql.QTAIMSceneLayer().draw_scene(viewer)
        assert viewer.spheres == [] and no_gl == []

    def test_layer_does_not_import_toga_or_the_viewer(self):
        code = "import sys, qtaimLayer; print('toga' in sys.modules, 'molecularViewer' in sys.modules)"
        out = subprocess.run([sys.executable, "-c", code], cwd=SRC, capture_output=True, text=True, check=True)
        assert out.stdout.strip() == "False False"
