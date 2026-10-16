"""Native Hirshfeld charges against Multiwfn 2026.7.15 (built-in free-atom densities)."""
import subprocess
import sys
from collections import defaultdict
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

import gaussianWavefunction as gw
import hirshfeldCharges as hc
import proatomDensities as pad
from help import AtomicData
from fukui_fixtures import CDFT, MINUS, N, PLUS

DATA = Path(__file__).resolve().parent / "data" / "hirshfeld"
SRC = Path(hc.__file__).resolve().parent


def multiwfn_reference():
    ref = defaultdict(list)
    for line in (DATA / "multiwfn_hirshfeld.txt").read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            name, q = line.split()
            ref[name].append(float(q))
    return dict(ref)


def phenol_multiwfn_charges():
    """q(N), q(N+1), q(N-1) columns of Multiwfn's CDFT.txt (4 decimals)."""
    lines = CDFT.read_text().splitlines()
    start = next(i for i, l in enumerate(lines) if l.split()[:4] == ["Atom", "q(N)", "q(N+1)", "q(N-1)"])
    rows = []
    for line in lines[start + 1:]:
        parts = line.split()
        if len(parts) < 7:
            break
        rows.append([float(v) for v in parts[-7:-4]])
    return np.array(rows)


class TestProatoms:
    def test_h_to_ar_are_tabulated(self):
        assert hc.supported_elements() == ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
                                           "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar"]

    @pytest.mark.parametrize("symbol", list(pad.DENSITIES))
    def test_each_table_holds_its_electrons(self, symbol):
        r = np.asarray(pad.RADII_BOHR)
        rho = np.asarray(pad.DENSITIES[symbol])
        z = AtomicData.atomic_numbers[symbol]
        assert np.all(rho >= 0) and np.all(rho[r < 20.0] > 0)  # tight atoms underflow far out
        assert np.trapezoid(4 * np.pi * r ** 2 * rho, r) == pytest.approx(z, rel=3e-4)

    def test_interpolation_reproduces_table_and_vanishes_far_away(self):
        r = np.asarray(pad.RADII_BOHR)
        assert np.allclose(hc.proatom_density("C", r[::50]), np.asarray(pad.DENSITIES["C"])[::50], rtol=1e-12)
        assert hc.proatom_density("C", np.array([100.0]))[0] == 0.0


class TestGrid:
    def test_grid_integrates_the_electron_count(self):
        w = gw.read_fchk(N)
        grid = hc.molecular_grid(w.coordinates, w.atom_symbols())
        assert len(grid.points) == 13 * 75 * 434
        assert grid.weights @ gw.electron_density(w, grid.points) == pytest.approx(50.0, abs=2e-4)

    @pytest.mark.skipif(not gw.numba_available(), reason="numba not installed (optional)")
    def test_becke_weights_numba_matches_numpy(self):
        w = gw.read_fchk(DATA / "ch3scl.fchk")
        fast = hc.molecular_grid(w.coordinates, w.atom_symbols(), use_numba=True)
        slow = hc.molecular_grid(w.coordinates, w.atom_symbols(), use_numba=False)
        # Becke's step 0.5(1 - mu) with mu -> 1 amplifies rounding order: compare loosely,
        # then strictly on what the weights are for, the integrated populations.
        scale = np.abs(slow.weights).max()
        assert np.max(np.abs(fast.weights - slow.weights)) <= 1e-12 * scale
        p_fast = hc.hirshfeld_populations(w, fast, use_numba=True)
        p_slow = hc.hirshfeld_populations(w, slow, use_numba=False)
        assert np.max(np.abs(p_fast - p_slow)) < 1e-9


class TestAgainstMultiwfn:
    @pytest.mark.parametrize("name", sorted(multiwfn_reference()))
    def test_validation_molecules(self, name):
        q = hc.hirshfeld_charges(gw.read_fchk(DATA / f"{name}.fchk"))
        assert np.max(np.abs(q - multiwfn_reference()[name])) < 3e-4

    def test_phenol_three_states(self):
        ref = phenol_multiwfn_charges()
        for k, path in enumerate((N, PLUS, MINUS)):
            q = hc.hirshfeld_charges(gw.read_fchk(path))
            assert np.max(np.abs(q - ref[:, k])) < 2e-4  # Multiwfn prints 4 decimals here

    def test_charges_sum_to_the_molecular_charge(self):
        w = gw.read_fchk(PLUS)
        assert hc.hirshfeld_charges(w).sum() == pytest.approx(w.charge, abs=1e-10)


class TestRejections:
    def test_element_without_free_atom_density(self):
        w = gw.read_fchk(DATA / "lih.fchk")
        heavy = replace(w, atomic_numbers=np.array([26, 1]), nuclear_charges=np.array([26.0, 1.0]))
        with pytest.raises(ValueError, match="No free-atom density for Fe"):
            hc.hirshfeld_charges(heavy)

    def test_effective_core_potential(self):
        w = gw.read_fchk(DATA / "nah.fchk")
        with pytest.raises(ValueError, match="effective core potentials"):
            hc.hirshfeld_charges(replace(w, nuclear_charges=np.array([9.0, 1.0])))


@pytest.mark.parametrize("module", ["hirshfeldCharges", "proatomDensities", "densityNumba"])
def test_module_does_not_import_toga(module):
    code = f"import sys, {module}; print('toga' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, cwd=str(SRC), check=True)
    assert out.stdout.strip() == "False"
