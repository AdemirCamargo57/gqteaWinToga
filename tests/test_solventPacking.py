import math
import pathlib
import subprocess
import sys

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

import solventPacking as sp
from single_solute_solvent_box import SingleSoluteSolventBox

SRC = pathlib.Path(__file__).resolve().parents[1]
needs_numba = pytest.mark.skipif(not sp.numba_available(), reason="numba not installed")


def legacy_clash(builder, cand_xyz, cand_el, other_xyz, other_el, extra=0.0):
    a = [[e, *map(float, p)] for e, p in zip(cand_el, cand_xyz)]
    b = [[e, *map(float, p)] for e, p in zip(other_el, other_xyz)]
    return builder.molecules_clash(a, b, extra=extra)


def legacy_builder(lengths, periodic, scale):
    b = SingleSoluteSolventBox()
    b.lattice_vectors = list(lengths)
    b.use_periodic_minimum_image = periodic
    b.vdw_scale = scale
    return b


@pytest.mark.parametrize("periodic", [True, False])
@pytest.mark.parametrize("extra", [0.0, 0.6])
def test_numpy_clash_decisions_match_legacy(periodic, extra):
    rng = np.random.default_rng(42)
    lengths = np.array([9.0, 10.0, 11.0])
    builder = legacy_builder(lengths, periodic, 0.8)
    elements = ["O", "H", "C", "N"]
    agree = clashes = 0
    for _ in range(600):
        n, m = rng.integers(1, 5), rng.integers(1, 12)
        ce = list(rng.choice(elements, n))
        oe = list(rng.choice(elements, m))
        cx = rng.uniform(0, lengths, size=(n, 3))
        ox = rng.uniform(0, lengths, size=(m, 3))
        cr = np.array([builder.get_vdw_radius(e) for e in ce])
        orr = np.array([builder.get_vdw_radius(e) for e in oe])
        expected = legacy_clash(builder, cx, ce, ox, oe, extra)
        got = sp.atoms_clash_numpy(cx, cr, ox, orr, 0.8, lengths, periodic, extra)
        agree += expected == got
        clashes += expected
    assert agree == 600
    assert 50 < clashes < 550  # the sample exercises both outcomes


@needs_numba
def test_numba_clash_decisions_match_numpy():
    rng = np.random.default_rng(7)
    lengths = np.array([8.0, 8.0, 8.0])
    for periodic in (True, False):
        for _ in range(400):
            cx = rng.uniform(0, 8, size=(3, 3))
            ox = rng.uniform(0, 8, size=(rng.integers(0, 20), 3))
            cr = rng.uniform(1.0, 1.8, size=3)
            orr = rng.uniform(1.0, 1.8, size=len(ox))
            extra = float(rng.choice([0.0, 0.5]))
            assert (sp.atoms_clash(cx, cr, ox, orr, 0.8, lengths, periodic, extra, sp.ENGINE_NUMBA)
                    == sp.atoms_clash_numpy(cx, cr, ox, orr, 0.8, lengths, periodic, extra))


@pytest.mark.parametrize("periodic", [True, False])
def test_grid_returns_exactly_what_the_legacy_centre_filter_keeps(periodic):
    rng = np.random.default_rng(3)
    lengths = np.array([20.0, 13.0, 31.0])
    radius = 2.0
    grid = sp.MoleculeGrid(lengths, periodic, min_cell=2 * radius + sp.CENTRE_FILTER_MARGIN)
    centres = rng.uniform(0, lengths, size=(300, 3))
    for c in centres:
        grid.add(c, radius, np.zeros((1, 3)), np.ones(1))
    builder = legacy_builder(lengths, periodic, 0.8)
    for probe in rng.uniform(0, lengths, size=(200, 3)):
        expected = []
        for index, c in enumerate(centres):
            dx, dy, dz = builder.minimum_image_displacement(*(probe - c))
            if not math.sqrt(dx * dx + dy * dy + dz * dz) > radius + radius + 2.5:
                expected.append(index)
        assert grid.nearby(probe, radius) == expected


@needs_numba
@pytest.mark.parametrize("periodic", [True, False])
def test_numba_grid_kernel_matches_numpy_grid(periodic):
    rng = np.random.default_rng(11)
    lengths = np.array([28.0, 26.0, 30.0])
    grid = sp.MoleculeGrid(lengths, periodic, min_cell=6.5)
    assert not grid.covers_all
    template = np.array([[0.0, 0.0, 0.1173], [0.0, 0.7572, -0.4692], [0.0, -0.7572, -0.4692]])
    radii = np.array([1.52, 1.2, 1.2])
    for c in rng.uniform(0, lengths, size=(400, 3)):
        grid.add(c, 2.0, template + c, radii)
    outcomes = set()
    for c in rng.uniform(0, lengths, size=(500, 3)):
        expected = grid.clashes(template + c, radii, c, 2.0, 0.8, 0.0, sp.ENGINE_NUMPY)
        got = grid.clashes(template + c, radii, c, 2.0, 0.8, 0.0, sp.ENGINE_NUMBA)
        assert expected == got
        outcomes.add(expected)
    assert outcomes == {True, False}


def test_small_boxes_with_fewer_than_three_cells_do_not_duplicate():
    grid = sp.MoleculeGrid([10.0, 10.0, 10.0], True, min_cell=6.0)
    assert grid.shape == (1, 1, 1)
    grid.add(np.array([1.0, 1.0, 1.0]), 1.0, np.zeros((1, 3)), np.ones(1))
    assert grid.nearby(np.array([9.0, 9.0, 9.0]), 1.0) == [0]


def test_rotation_stream_matches_sequential_rotation_apply():
    template = np.array([[0.0, 0.0, 0.1173], [0.0, 0.7572, -0.4692], [0.0, -0.7572, -0.4692],
                         [1.3, -0.2, 0.9]])
    np.random.seed(99)
    expected = [Rotation.random().apply(template) for _ in range(1300)]
    np.random.seed(99)
    stream = sp.RotationStream(block=512)
    got = [template @ stream.next().T for _ in range(1300)]
    assert all(np.array_equal(e, g) for e, g in zip(expected, got))


def test_resolve_engine():
    assert sp.resolve_engine("numpy") == "numpy"
    assert sp.resolve_engine("legacy") == "legacy"
    assert sp.resolve_engine("auto") in ("numpy", "numba")
    with pytest.raises(ValueError):
        sp.resolve_engine("gpu")


def test_numba_is_only_imported_by_its_kernel_module():
    hits = [p.name for p in SRC.glob("*.py")
            if p.name != "espNumba.py" and "import numba" in p.read_text(encoding="utf-8", errors="ignore")
            or (p.name != "espNumba.py" and "from numba" in p.read_text(encoding="utf-8", errors="ignore"))]
    assert hits == ["solventPackingNumba.py"]


def test_module_does_not_import_toga():
    code = "import sys, solventPacking; print('toga' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, cwd=str(SRC))
    assert out.stdout.strip() == "False", out.stderr
