"""Tests for the Mean Residence Time tool (``meanResidenceTime.py``).

Run from ``venv/src/``::

    python -m pytest tests/ -q

The Toga GUI is never instantiated -- ``MeanResidenceTimeCalculator`` is pure
numpy/IO, so every test drives it directly.

Expected values are analytic wherever possible.  For an occupancy train made of
runs of length ``n`` frames the two continuous time constants are known exactly:

* mean residence time      = <T>              = n*dt
* IMM survival integral    = <T^2>/(2<T>)     = (n+1)/2 * dt   (discrete form)

so the tests pin the equations rather than freezing whatever the code happens to
return.
"""
import asyncio
import glob
import math
import os

import numpy as np
import pytest

from meanResidenceTime import MeanResidenceTimeCalculator


@pytest.fixture
def calc():
    return MeanResidenceTimeCalculator()


# --------------------------------------------------------------------------- #
# Helpers                                                                       #
# --------------------------------------------------------------------------- #
def occupancy_from_runs(run_lengths, gap=3, n_objects=1):
    """Build a 1-object occupancy row: run, gap, run, gap, ... with leading gap."""
    row = [0] * gap
    for n in run_lengths:
        row.extend([1] * n)
        row.extend([0] * gap)
    return np.array([row] * n_objects, dtype=int)


def brute_force_origin_survival(occupancy, max_lag):
    """Reference implementation: scan every time origin explicitly.

    S(lag) = (number of (object, origin) pairs occupied continuously through lag)
             / (number of occupied origins)
    """
    occ = occupancy.astype(bool)
    n_obj, n_frames = occ.shape
    total_origins = occ.sum()
    out = np.zeros(max_lag, dtype=float)
    for lag in range(max_lag):
        survivors = 0
        for o in range(n_obj):
            for t0 in range(n_frames):
                if not occ[o, t0]:
                    continue
                if t0 + lag < n_frames and occ[o, t0:t0 + lag + 1].all():
                    survivors += 1
        out[lag] = survivors / total_origins
    return out


# --------------------------------------------------------------------------- #
# Origin-averaged (IMM / Garcia-Stiller) survival                               #
# --------------------------------------------------------------------------- #
class TestOriginAveragedSurvival:
    def test_matches_brute_force_origin_scan(self, calc):
        """The fast run-length form must equal an explicit scan over time origins."""
        rng = np.random.default_rng(20240905)
        occ = (rng.random((4, 120)) < 0.35).astype(int)
        max_lag = 25

        _, fast = calc.compute_origin_averaged_survival(occ, dt=1.0, max_lag=max_lag)
        slow = brute_force_origin_survival(occ, max_lag)

        assert np.allclose(fast, slow), "run-length survival disagrees with origin scan"

    def test_single_run_is_linear_ramp(self, calc):
        """One run of n frames gives S(l) = (n-l)/n exactly."""
        n = 10
        occ = np.array([[0, 0] + [1] * n + [0, 0]], dtype=int)
        _, S = calc.compute_origin_averaged_survival(occ, dt=1.0, max_lag=n)
        expected = np.array([(n - l) / n for l in range(n)])
        assert np.allclose(S, expected)

    def test_starts_at_one(self, calc):
        occ = occupancy_from_runs([4, 7, 2])
        _, S = calc.compute_origin_averaged_survival(occ, dt=0.5, max_lag=10)
        assert S[0] == pytest.approx(1.0)

    def test_integral_is_mean_residual_time(self, calc):
        """For uniform runs of n frames the integral is (n+1)/2 * dt."""
        n, dt = 8, 0.25
        occ = occupancy_from_runs([n, n, n, n])
        _, S = calc.compute_origin_averaged_survival(occ, dt=dt, max_lag=n)
        assert calc.discrete_integral(S, dt) == pytest.approx((n + 1) / 2 * dt)

    def test_empty_occupancy_is_safe(self, calc):
        occ = np.zeros((2, 20), dtype=int)
        t, S = calc.compute_origin_averaged_survival(occ, dt=1.0, max_lag=5)
        assert len(t) == len(S)
        assert np.all(S == 0.0)


# --------------------------------------------------------------------------- #
# The two continuous time constants are genuinely different                     #
# --------------------------------------------------------------------------- #
class TestResidenceVersusSurvivalIntegral:
    def test_uniform_runs_give_the_two_known_values(self, calc):
        n, dt = 6, 0.5
        occ = occupancy_from_runs([n] * 5)

        durations = calc.extract_event_durations(occ, dt, censor_boundary=True)
        residence = float(np.mean(durations))

        _, S = calc.compute_origin_averaged_survival(occ, dt=dt, max_lag=n)
        imm = calc.discrete_integral(S, dt)

        assert residence == pytest.approx(n * dt)          # <T>
        assert imm == pytest.approx((n + 1) / 2 * dt)      # <T^2>/(2<T>)

    def test_broad_distribution_separates_them(self, calc):
        """With a broad duration distribution the two must not coincide.

        This is the regression that motivated the change: the old code reported
        the mean event duration twice under two different names.
        """
        runs = [2] * 8 + [60]
        occ = occupancy_from_runs(runs, gap=4)
        dt = 1.0

        durations = calc.extract_event_durations(occ, dt, censor_boundary=True)
        residence = float(np.mean(durations))
        _, S = calc.compute_origin_averaged_survival(occ, dt=dt, max_lag=80)
        imm = calc.discrete_integral(S, dt)

        # Analytic values for a set of runs {n_i}:
        #   <T>          = sum(n_i) / count
        #   IMM integral = sum(n_i(n_i+1)/2) / sum(n_i)
        total = sum(runs)
        assert residence == pytest.approx(total / len(runs) * dt)
        assert imm == pytest.approx(sum(n * (n + 1) / 2 for n in runs) / total * dt)

        # ...and those are far apart, which is the whole point.
        assert imm > 2.5 * residence, (
            f"broad distribution should separate the estimators "
            f"(residence={residence}, imm={imm})"
        )


# --------------------------------------------------------------------------- #
# Boundary censoring                                                            #
# --------------------------------------------------------------------------- #
class TestBoundaryCensoring:
    OCC = np.array([[1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1]], dtype=int)

    def test_uncensored_keeps_truncated_runs(self, calc):
        durations = calc.extract_event_durations(self.OCC, dt=1.0, censor_boundary=False)
        assert sorted(durations.tolist()) == [2.0, 5.0, 5.0]

    def test_censored_drops_runs_touching_either_edge(self, calc):
        durations = calc.extract_event_durations(self.OCC, dt=1.0, censor_boundary=True)
        assert durations.tolist() == [2.0], "only the interior event is complete"

    def test_counts_censored_events(self, calc):
        assert calc.count_boundary_events(self.OCC) == 2

    def test_interior_only_is_unaffected(self, calc):
        occ = np.array([[0, 1, 1, 0, 1, 0]], dtype=int)
        assert calc.count_boundary_events(occ) == 0
        assert np.allclose(
            calc.extract_event_durations(occ, 1.0, censor_boundary=True),
            calc.extract_event_durations(occ, 1.0, censor_boundary=False),
        )

    def test_fully_occupied_row_yields_no_complete_events(self, calc):
        occ = np.ones((1, 10), dtype=int)
        assert calc.extract_event_durations(occ, 1.0, censor_boundary=True).size == 0
        assert calc.count_boundary_events(occ) == 1

    def test_censoring_removes_downward_bias(self, calc):
        """Truncated edge runs drag the mean down; censoring must fix that."""
        occ = np.array([[1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1]], dtype=int)
        biased = calc.extract_event_durations(occ, 1.0, censor_boundary=False).mean()
        clean = calc.extract_event_durations(occ, 1.0, censor_boundary=True).mean()
        assert clean > biased


# --------------------------------------------------------------------------- #
# Periodic boundary conditions                                                  #
# --------------------------------------------------------------------------- #
class TestPeriodicBoundaries:
    def test_pair_across_boundary_is_close_with_pbc(self, calc):
        box = (10.0, 10.0, 10.0)
        ref = np.array([[0.5, 5.0, 5.0]])
        obs = np.array([[[9.5, 5.0, 5.0]]])       # 1.0 A away through the wall

        _, occ_no = calc.build_occupancy(ref, obs, cutoff=2.0, cell_lengths=None)
        _, occ_pbc = calc.build_occupancy(ref, obs, cutoff=2.0, cell_lengths=box)

        assert occ_no[0, 0] == 0, "without PBC the wrapped pair looks 9 A apart"
        assert occ_pbc[0, 0] == 1, "with PBC the pair is 1 A apart"

    def test_box_shift_invariance(self, calc):
        """Translating every coordinate by a whole box must not change occupancy."""
        rng = np.random.default_rng(7)
        box = np.array([12.0, 11.0, 13.0])
        ref = rng.random((30, 3)) * box
        obs = rng.random((5, 30, 3)) * box

        _, base = calc.build_occupancy(ref, obs, cutoff=3.0, cell_lengths=tuple(box))
        _, shifted = calc.build_occupancy(
            ref + box, obs + box, cutoff=3.0, cell_lengths=tuple(box)
        )
        assert np.array_equal(base, shifted)

    def test_distance_never_exceeds_half_box(self, calc):
        rng = np.random.default_rng(3)
        box = np.array([8.0, 8.0, 8.0])
        ref = rng.random((40, 3)) * box
        obs = rng.random((3, 40, 3)) * box
        d, _ = calc.build_occupancy(ref, obs, cutoff=1.0, cell_lengths=tuple(box))
        assert d.max() <= math.sqrt(3) * 4.0 + 1e-9

    def test_cutoff_above_half_the_shortest_edge_is_rejected(self, calc):
        ref = np.zeros((3, 3))
        obs = np.zeros((1, 3, 3))
        with pytest.raises(ValueError, match="half"):
            calc.build_occupancy(ref, obs, cutoff=5.1, cell_lengths=(10.0, 12.0, 14.0))

    def test_cutoff_exactly_half_is_accepted(self, calc):
        ref = np.zeros((3, 3))
        obs = np.zeros((1, 3, 3))
        calc.build_occupancy(ref, obs, cutoff=5.0, cell_lengths=(10.0, 12.0, 14.0))

    def test_non_positive_edge_is_rejected(self, calc):
        ref = np.zeros((3, 3))
        obs = np.zeros((1, 3, 3))
        with pytest.raises(ValueError):
            calc.build_occupancy(ref, obs, cutoff=1.0, cell_lengths=(10.0, 0.0, 14.0))

    def test_parse_cell_lengths(self, calc):
        assert calc.parse_cell_lengths("10 11 12") == (10.0, 11.0, 12.0)
        assert calc.parse_cell_lengths("10,11,12") == (10.0, 11.0, 12.0)
        assert calc.parse_cell_lengths("") is None
        assert calc.parse_cell_lengths("   ") is None
        with pytest.raises(ValueError):
            calc.parse_cell_lengths("10 11")


# --------------------------------------------------------------------------- #
# Bounded intermittent integral                                                 #
# --------------------------------------------------------------------------- #
class TestIntermittentBounding:
    def test_default_max_lag_is_a_tenth_of_the_trajectory(self, calc):
        occ = (np.random.default_rng(1).random((2, 500)) < 0.3).astype(int)
        t, C, R, _ = calc.compute_intermittent_functions(occ, dt=1.0)
        assert len(t) == len(C) == len(R) == 50

    def test_explicit_max_lag_is_honoured(self, calc):
        occ = (np.random.default_rng(1).random((2, 500)) < 0.3).astype(int)
        t, C, _, _ = calc.compute_intermittent_functions(occ, dt=1.0, max_lag=17)
        assert len(t) == 17

    def test_max_lag_is_clamped_to_trajectory_length(self, calc):
        occ = np.ones((1, 8), dtype=int)
        t, _, _, _ = calc.compute_intermittent_functions(occ, dt=1.0, max_lag=999)
        assert len(t) <= 8

    def test_c_starts_at_one(self, calc):
        occ = (np.random.default_rng(5).random((3, 200)) < 0.4).astype(int)
        _, C, _, _ = calc.compute_intermittent_functions(occ, dt=1.0)
        assert C[0] == pytest.approx(1.0)

    def test_zero_crossing_truncation_stops_at_first_sign_change(self, calc):
        y = np.array([1.0, 0.6, 0.3, 0.1, -0.2, 0.4, -0.5])
        value, n_used, rule = calc.integrate_correlation(y, dt=1.0, rule="zero_crossing")
        assert n_used == 4
        assert value == pytest.approx(2.0)
        assert rule == "zero_crossing"

    def test_zero_crossing_falls_back_to_full_axis(self, calc):
        y = np.array([1.0, 0.5, 0.25])
        value, n_used, _ = calc.integrate_correlation(y, dt=2.0, rule="zero_crossing")
        assert n_used == 3
        assert value == pytest.approx(3.5)

    def test_full_rule_uses_everything(self, calc):
        y = np.array([1.0, -1.0, 1.0, -1.0])
        value, n_used, rule = calc.integrate_correlation(y, dt=1.0, rule="full")
        assert n_used == 4
        assert value == pytest.approx(0.0)
        assert rule == "full"

    def test_exponential_fit_recovers_known_tau(self, calc):
        tau, dt = 7.5, 0.5
        t = np.arange(60) * dt
        y = np.exp(-t / tau)
        value, _, rule = calc.integrate_correlation(y, dt=dt, rule="exponential")
        assert value == pytest.approx(tau, rel=1e-3)
        assert rule == "exponential"

    def test_unknown_rule_is_rejected(self, calc):
        with pytest.raises(ValueError):
            calc.integrate_correlation(np.array([1.0]), dt=1.0, rule="nonsense")


# --------------------------------------------------------------------------- #
# Uncertainty                                                                   #
# --------------------------------------------------------------------------- #
class TestUncertainty:
    def test_constant_sample_has_zero_uncertainty(self, calc):
        assert calc.block_average_uncertainty(np.full(50, 3.0)) == pytest.approx(0.0)

    def test_uncertainty_is_positive_for_scattered_data(self, calc):
        rng = np.random.default_rng(2)
        assert calc.block_average_uncertainty(rng.normal(5.0, 2.0, 400)) > 0.0

    def test_too_few_samples_gives_nan(self, calc):
        assert math.isnan(calc.block_average_uncertainty(np.array([1.0])))

    def test_uncertainty_shrinks_with_more_data(self, calc):
        rng = np.random.default_rng(11)
        small = calc.block_average_uncertainty(rng.normal(0.0, 1.0, 200))
        large = calc.block_average_uncertainty(rng.normal(0.0, 1.0, 20000))
        assert large < small


# --------------------------------------------------------------------------- #
# t* tolerance scan                                                             #
# --------------------------------------------------------------------------- #
class TestToleranceScan:
    def test_scan_returns_one_row_per_tolerance(self, calc):
        occ = occupancy_from_runs([3, 3, 3], gap=2)
        rows = calc.scan_tolerance(occ, dt=1.0, tolerances=[0, 1, 2, 5])
        assert [r["tolerance_frames"] for r in rows] == [0, 1, 2, 5]

    def test_residence_time_grows_with_tolerance(self, calc):
        """Bridging short gaps can only merge events, never split them."""
        occ = occupancy_from_runs([4, 4, 4, 4], gap=2)
        rows = calc.scan_tolerance(occ, dt=1.0, tolerances=[0, 1, 2, 3])
        values = [r["mean_residence_time"] for r in rows]
        assert all(b >= a for a, b in zip(values, values[1:])), values
        assert values[-1] > values[0], "a large t* must merge the runs"

    def test_tolerance_does_not_bridge_the_leading_gap(self, calc):
        occ = np.array([[0, 0, 1, 1, 0, 0]], dtype=int)
        bridged = calc.apply_tolerance_to_occupancy(occ, tolerance_frames=5)
        assert bridged.tolist() == occ.tolist()


# --------------------------------------------------------------------------- #
# End-to-end                                                                    #
# --------------------------------------------------------------------------- #
class TestEndToEnd:
    @staticmethod
    def write_xyz(path, frames):
        """frames: list of list of (symbol, (x, y, z))."""
        with open(path, "w", encoding="utf-8") as f:
            for i, frame in enumerate(frames):
                f.write(f"{len(frame)}\n")
                f.write(f"frame {i}\n")
                for sym, (x, y, z) in frame:
                    f.write(f"{sym} {x:.6f} {y:.6f} {z:.6f}\n")

    def make_trajectory(self, tmp_path, distances):
        """Atom 1 fixed at the origin; atom 2 placed at the given distance."""
        frames = [
            [("O", (0.0, 0.0, 0.0)), ("H", (d, 0.0, 0.0))] for d in distances
        ]
        path = tmp_path / "traj.xyz"
        self.write_xyz(path, frames)
        return str(path)

    def test_known_occupancy_pattern(self, tmp_path, calc):
        # in / in / out / in / in / in  with cutoff 2.0
        path = self.make_trajectory(tmp_path, [1.0, 1.0, 5.0, 1.0, 1.0, 1.0])
        res = calc.run(
            xyz_file=path, dt=1.0, dt_unit="fs", cutoff=2.0, tolerance_frames=0,
            reference_mode="single_atom", reference_definition="1",
            observed_mode="single_atom", observed_definition="2",
        )
        assert res.n_frames == 6
        assert res.occupancy_mean == pytest.approx(5 / 6)
        # both runs touch a trajectory edge, so both are censored
        assert res.n_censored_events == 2
        assert res.n_continuous_events == 0

    def test_interior_event_survives_censoring(self, tmp_path, calc):
        path = self.make_trajectory(tmp_path, [5.0, 1.0, 1.0, 1.0, 5.0])
        res = calc.run(
            xyz_file=path, dt=2.0, dt_unit="fs", cutoff=2.0, tolerance_frames=0,
            reference_mode="single_atom", reference_definition="1",
            observed_mode="single_atom", observed_definition="2",
        )
        assert res.n_continuous_events == 1
        assert res.n_censored_events == 0
        assert res.mean_residence_time == pytest.approx(6.0)      # 3 frames * 2.0

    def test_geometric_mode_does_not_need_masses(self, tmp_path, calc):
        """An element missing from the mass table must not break a mass-free mode."""
        frames = [
            [("Xx", (0.0, 0.0, 0.0)), ("Xx", (2.0, 0.0, 0.0)), ("H", (0.5, 0.0, 0.0))]
            for _ in range(4)
        ]
        path = tmp_path / "u.xyz"
        self.write_xyz(path, frames)
        res = calc.run(
            xyz_file=str(path), dt=1.0, dt_unit="fs", cutoff=3.0, tolerance_frames=0,
            reference_mode="geometric_center", reference_definition="1,2",
            observed_mode="single_atom", observed_definition="3",
        )
        assert res.n_frames == 4

    def test_center_of_mass_mode_still_reports_missing_mass(self, tmp_path, calc):
        """A genuinely unknown symbol must still be reported, clearly."""
        frames = [[("Xx", (0.0, 0.0, 0.0)), ("H", (1.0, 0.0, 0.0))] for _ in range(3)]
        path = tmp_path / "u2.xyz"
        self.write_xyz(path, frames)
        with pytest.raises(ValueError, match="mass"):
            calc.run(
                xyz_file=str(path), dt=1.0, dt_unit="fs", cutoff=3.0, tolerance_frames=0,
                reference_mode="center_of_mass", reference_definition="1,2",
                observed_mode="single_atom", observed_definition="2",
            )

    def test_pbc_changes_the_result(self, tmp_path, calc):
        """Same trajectory, with and without a box, must disagree."""
        path = self.make_trajectory(tmp_path, [9.5, 9.5, 9.5, 9.5])
        common = dict(
            xyz_file=path, dt=1.0, dt_unit="fs", cutoff=2.0, tolerance_frames=0,
            reference_mode="single_atom", reference_definition="1",
            observed_mode="single_atom", observed_definition="2",
        )
        without = calc.run(**common)
        with_box = calc.run(**common, cell_lengths=(10.0, 10.0, 10.0))
        assert without.occupancy_mean == pytest.approx(0.0)
        assert with_box.occupancy_mean == pytest.approx(1.0)

    def test_summary_reports_the_new_quantities(self, tmp_path, calc):
        path = self.make_trajectory(tmp_path, [5.0, 1.0, 1.0, 1.0, 5.0])
        res = calc.run(
            xyz_file=path, dt=1.0, dt_unit="fs", cutoff=2.0, tolerance_frames=0,
            reference_mode="single_atom", reference_definition="1",
            observed_mode="single_atom", observed_definition="2",
        )
        text = calc.format_results_text(res)
        for expected in (
            "Mean residence time",
            "IMM survival integral",
            "Censored",
            "Periodic boundary",
        ):
            assert expected in text, f"summary is missing '{expected}'"

    def test_vectorised_group_centres_match_the_loop(self, tmp_path, calc):
        rng = np.random.default_rng(4)
        n_frames, n_atoms = 6, 8
        frames = [
            [("C", tuple(rng.random(3) * 10)) for _ in range(n_atoms)]
            for _ in range(n_frames)
        ]
        path = tmp_path / "g.xyz"
        self.write_xyz(path, frames)
        symbols, coords = calc.read_xyz(str(path))
        masses = calc.get_masses(symbols)

        groups = [[0, 1, 2], [3, 4], [5, 6, 7]]
        got, _ = calc.build_observed_positions(
            coords, masses, "group_list_geometric", "1,2,3; 4,5; 6,7,8"
        )
        expected = np.stack([coords[:, g, :].mean(axis=1) for g in groups])
        assert np.allclose(got, expected)


# --------------------------------------------------------------------------- #
# Vectorised centre construction (must match the per-frame loop exactly)        #
# --------------------------------------------------------------------------- #
class TestCentreConstruction:
    @staticmethod
    def reference_loop_geometric(coords, indices):
        out = np.zeros((coords.shape[0], 3))
        for i in range(coords.shape[0]):
            out[i] = coords[i][indices].mean(axis=0)
        return out

    @staticmethod
    def reference_loop_com(coords, indices, masses):
        out = np.zeros((coords.shape[0], 3))
        m = masses[indices]
        for i in range(coords.shape[0]):
            out[i] = (coords[i][indices] * m[:, None]).sum(axis=0) / m.sum()
        return out

    @pytest.fixture
    def sample(self):
        rng = np.random.default_rng(99)
        coords = rng.random((25, 9, 3)) * 15.0
        masses = np.array([12.011, 1.008, 15.999, 1.008, 12.011,
                           15.999, 1.008, 32.06, 1.008])
        return coords, masses

    def test_geometric_reference_matches_loop(self, calc, sample):
        coords, masses = sample
        got = calc.build_reference_positions(coords, masses, "geometric_center", "1,3,5")
        assert np.allclose(got, self.reference_loop_geometric(coords, [0, 2, 4]))

    def test_com_reference_matches_loop(self, calc, sample):
        coords, masses = sample
        got = calc.build_reference_positions(coords, masses, "center_of_mass", "1,3,5")
        assert np.allclose(got, self.reference_loop_com(coords, [0, 2, 4], masses))

    def test_com_is_mass_weighted_not_geometric(self, calc, sample):
        """A heavy/light pair must not give the midpoint."""
        coords, masses = sample
        geo = calc.build_reference_positions(coords, masses, "geometric_center", "2,8")
        com = calc.build_reference_positions(coords, masses, "center_of_mass", "2,8")
        assert not np.allclose(geo, com)

    def test_group_com_matches_loop(self, calc, sample):
        coords, masses = sample
        got, labels = calc.build_observed_positions(
            coords, masses, "group_list_com", "1,2,3; 4,5; 6,7,8,9"
        )
        for k, group in enumerate([[0, 1, 2], [3, 4], [5, 6, 7, 8]]):
            assert np.allclose(got[k], self.reference_loop_com(coords, group, masses))
        assert len(labels) == 3

    def test_single_atom_reference_is_the_atom_itself(self, calc, sample):
        coords, masses = sample
        got = calc.build_reference_positions(coords, masses, "single_atom", "4")
        assert np.allclose(got, coords[:, 3, :])

    def test_atom_list_individual_shape_and_order(self, calc, sample):
        coords, masses = sample
        got, labels = calc.build_observed_positions(
            coords, masses, "atom_list_individual", "9,1,5"
        )
        assert got.shape == (3, coords.shape[0], 3)
        assert np.allclose(got[0], coords[:, 8, :])
        assert labels == ["Atom 9", "Atom 1", "Atom 5"]

    def test_single_atom_group_equals_the_atom(self, calc, sample):
        coords, masses = sample
        got, _ = calc.build_observed_positions(coords, masses, "group_list_geometric", "7")
        assert np.allclose(got[0], coords[:, 6, :])


# --------------------------------------------------------------------------- #
# UI plumbing that can be checked without a display                             #
# --------------------------------------------------------------------------- #
class FakeWidget:
    def __init__(self):
        self.text = ""
        self.value = 0


class FakeWindow:
    """Records which Toga dialog type was raised.

    Toga 0.5.3 dialogs keep their title/message inside the backend impl and
    expose neither, so the type name is what can be asserted portably.
    """

    def __init__(self):
        self.dialogs = []

    async def dialog(self, dlg):
        self.dialogs.append(type(dlg).__name__)
        return None


class TestUIPlumbing:
    @staticmethod
    def bare_ui():
        """A MeanResidenceTimeUI with stub widgets and no Toga window."""
        from meanResidenceTime import MeanResidenceTimeUI
        ui = MeanResidenceTimeUI.__new__(MeanResidenceTimeUI)
        ui.status_label = FakeWidget()
        ui.progress_bar = FakeWidget()
        ui.summary_output = FakeWidget()
        ui.window = FakeWindow()
        ui.results = None
        return ui

    def test_accepts_the_button_toga_passes(self):
        """gqteaWinToga registers the class itself as on_press, so __init__ is
        handed the Button. It must take *args like every other tool."""
        import inspect
        from meanResidenceTime import MeanResidenceTimeUI
        params = list(inspect.signature(MeanResidenceTimeUI.__init__).parameters.values())
        assert any(p.kind is inspect.Parameter.VAR_POSITIONAL for p in params), (
            "MeanResidenceTimeUI.__init__ must accept *args"
        )

    def test_clear_output_discards_stale_results(self):
        ui = self.bare_ui()
        ui.results = object()
        ui.clear_output(None)
        assert ui.results is None, (
            "Export after Clear must not silently write the previous run"
        )

    def test_clear_output_resets_the_progress_bar(self):
        ui = self.bare_ui()
        ui.progress_bar.value = 100
        ui.clear_output(None)
        assert ui.progress_bar.value == 0

    def test_show_error_uses_the_dialog_api(self):
        ui = self.bare_ui()
        asyncio.get_event_loop_policy().new_event_loop().run_until_complete(
            ui._show_error("Title", "Message")
        )
        assert ui.window.dialogs == ["ErrorDialog"]

    def test_show_info_uses_the_dialog_api(self):
        ui = self.bare_ui()
        asyncio.get_event_loop_policy().new_event_loop().run_until_complete(
            ui._show_info("Title", "Message")
        )
        assert ui.window.dialogs == ["InfoDialog"]

    def test_export_without_results_is_reported_not_crashed(self):
        ui = self.bare_ui()
        ui.calculator = MeanResidenceTimeCalculator()
        asyncio.get_event_loop_policy().new_event_loop().run_until_complete(
            ui.export_results(None)
        )
        assert ui.window.dialogs, "the user must be told there is nothing to export"

    def test_broken_helper_is_gone(self):
        """open_mean_residence_time_window() called __init__ with an 'app'
        keyword it never accepted, so it could only ever raise TypeError."""
        import meanResidenceTime
        assert not hasattr(meanResidenceTime, "open_mean_residence_time_window")

    def test_no_deprecated_toga_dialog_calls(self):
        """The rest of the codebase uses window.dialog(toga.XDialog(...))."""
        import inspect
        import meanResidenceTime
        source = inspect.getsource(meanResidenceTime)
        for deprecated in ("error_dialog(", "info_dialog(", "question_dialog("):
            assert deprecated not in source, f"deprecated Toga API: {deprecated}"


# --------------------------------------------------------------------------- #
# Shared reference data                                                         #
# --------------------------------------------------------------------------- #
class TestMassSource:
    def test_masses_come_from_the_shared_atomic_data_hub(self, calc):
        """help.py AtomicData is the project's reference-data hub; the tool must
        not carry a second, shorter copy of the periodic table."""
        from help import AtomicData
        assert calc.ATOMIC_MASSES is AtomicData.atomic_masses

    def test_covers_elements_the_old_local_table_missed(self, calc):
        for element in ("U", "Th", "Np", "Pu"):
            assert element in calc.ATOMIC_MASSES, element

    def test_common_organic_masses_are_unchanged(self, calc):
        # The elements that actually matter for these trajectories must not
        # have shifted when the table was swapped.
        for element, mass in [("H", 1.00784), ("C", 12.011),
                              ("N", 14.007), ("O", 15.999)]:
            assert calc.ATOMIC_MASSES[element] == pytest.approx(mass, abs=1e-6)


# --------------------------------------------------------------------------- #
# Mode labels and mode-driven help                                              #
# --------------------------------------------------------------------------- #
class TestModePresentation:
    @staticmethod
    def ui_class():
        from meanResidenceTime import MeanResidenceTimeUI
        return MeanResidenceTimeUI

    def test_labels_map_onto_the_calculator_modes(self):
        ui = self.ui_class()
        from meanResidenceTime import MeanResidenceTimeCalculator as Calc
        assert set(ui.REFERENCE_MODE_LABELS.values()) == {
            "single_atom", "geometric_center", "center_of_mass"
        }
        assert set(ui.OBSERVED_MODE_LABELS.values()) == {
            "single_atom", "atom_list_individual",
            "group_list_geometric", "group_list_com"
        }
        # every mass-dependent mode must be reachable from the UI
        reachable = set(ui.REFERENCE_MODE_LABELS.values()) | set(ui.OBSERVED_MODE_LABELS.values())
        assert Calc.MASS_DEPENDENT_MODES <= reachable

    def test_labels_are_human_readable_not_identifiers(self):
        ui = self.ui_class()
        for label in list(ui.REFERENCE_MODE_LABELS) + list(ui.OBSERVED_MODE_LABELS):
            assert "_" not in label, f"{label!r} is still a raw identifier"

    @pytest.mark.parametrize("mode", [
        "single_atom", "geometric_center", "center_of_mass",
        "atom_list_individual", "group_list_geometric", "group_list_com",
    ])
    def test_every_mode_has_a_hint_and_an_example(self, mode):
        ui = self.ui_class()
        hint = ui.mode_hint(mode)
        assert hint and isinstance(hint, str)
        assert ui.mode_example(mode)

    def test_group_modes_document_the_semicolon(self):
        ui = self.ui_class()
        for mode in ("group_list_geometric", "group_list_com"):
            assert ";" in ui.mode_example(mode), mode

    def test_single_atom_example_is_one_index(self):
        ui = self.ui_class()
        assert "," not in ui.mode_example("single_atom")
        assert ";" not in ui.mode_example("single_atom")

    def test_unknown_mode_falls_back_without_raising(self):
        ui = self.ui_class()
        assert isinstance(ui.mode_hint("no_such_mode"), str)

    def test_long_definition_is_abbreviated(self, calc=None):
        """Selecting every solvent H makes a definition thousands of chars long."""
        from meanResidenceTime import MeanResidenceTimeCalculator as C
        long_def = ",".join(str(i) for i in range(1, 300))
        out = C._abbreviate(long_def)
        assert len(out) < 200
        assert "299 entries in total" in out

    def test_short_definition_is_left_alone(self):
        from meanResidenceTime import MeanResidenceTimeCalculator as C
        assert C._abbreviate("1,2,3") == "1,2,3"

    def test_many_object_labels_are_summarised(self):
        from meanResidenceTime import MeanResidenceTimeCalculator as C
        labels = [f"Atom {i}" for i in range(1, 281)]
        out = C._abbreviate_labels(labels)
        assert len(out) == 8
        assert out[0] == "Atom 1"
        assert out[-1] == "Atom 280"
        assert "273 more" in out[5]

    def test_few_object_labels_are_listed_in_full(self):
        from meanResidenceTime import MeanResidenceTimeCalculator as C
        labels = [f"Atom {i}" for i in range(1, 5)]
        assert C._abbreviate_labels(labels) == labels

    def test_help_entry_exists_and_is_not_the_legacy_one(self):
        """The legacy tool keeps help_mrt; this tool needs its own text."""
        from help import HelpGqteaWin
        assert hasattr(HelpGqteaWin, "help_mrt_advanced")
        assert HelpGqteaWin.help_mrt_advanced != HelpGqteaWin.help_mrt
        text = HelpGqteaWin.help_mrt_advanced
        for topic in ("Cell lengths", "t*", "IMM survival integral",
                      "censored", "Mean residence time"):
            assert topic in text, f"help text does not cover {topic!r}"


# --------------------------------------------------------------------------- #
# Output directory                                                              #
# --------------------------------------------------------------------------- #
class TestOutputDirectory:
    def test_defaults_to_the_trajectory_folder(self, calc, tmp_path):
        xyz = tmp_path / "sub" / "traj.xyz"
        assert calc.resolve_output_dir(str(xyz), "") == str(xyz.parent)
        assert calc.resolve_output_dir(str(xyz), None) == str(xyz.parent)

    def test_explicit_choice_wins(self, calc, tmp_path):
        other = tmp_path / "elsewhere"
        other.mkdir()
        assert calc.resolve_output_dir(str(tmp_path / "t.xyz"), str(other)) == str(other)

    def test_missing_directory_is_rejected(self, calc, tmp_path):
        with pytest.raises(ValueError, match="directory"):
            calc.resolve_output_dir(str(tmp_path / "t.xyz"), str(tmp_path / "nope"))

    def test_export_writes_into_the_chosen_directory(self, calc, tmp_path):
        frames = [[("O", (0.0, 0.0, 0.0)), ("H", (d, 0.0, 0.0))]
                  for d in (5.0, 1.0, 1.0, 1.0, 5.0)]
        xyz = tmp_path / "traj.xyz"
        with open(xyz, "w", encoding="utf-8") as f:
            for i, fr in enumerate(frames):
                f.write(f"{len(fr)}\nframe {i}\n")
                for s, (x, y, z) in fr:
                    f.write(f"{s} {x} {y} {z}\n")
        out = tmp_path / "results"
        out.mkdir()

        res = calc.run(
            xyz_file=str(xyz), dt=1.0, dt_unit="fs", cutoff=2.0, tolerance_frames=0,
            reference_mode="single_atom", reference_definition="1",
            observed_mode="single_atom", observed_definition="2",
        )
        files = calc.export_curves(res, output_dir=str(out))
        for path in files.values():
            assert str(out) in path
            assert os.path.exists(path)


# --------------------------------------------------------------------------- #
# Plots                                                                         #
# --------------------------------------------------------------------------- #
class TestPlots:
    @staticmethod
    def bare_ui(tmp_path):
        from meanResidenceTime import MeanResidenceTimeUI
        ui = MeanResidenceTimeUI.__new__(MeanResidenceTimeUI)
        ui.output_dir = str(tmp_path)
        ui.saved_plot_data = []
        ui.saved_plot_files = []
        return ui

    @pytest.fixture
    def results(self, calc, tmp_path):
        occ = occupancy_from_runs([4, 9, 3, 7], gap=3)
        durations = calc.extract_event_durations(occ, 1.0, censor_boundary=True)
        t_s, S = calc.compute_origin_averaged_survival(occ, 1.0, max_lag=10)
        t_c, C, R, _ = calc.compute_intermittent_functions(occ, 1.0, max_lag=10)
        from meanResidenceTime import MRTResults
        return MRTResults(
            xyz_file=str(tmp_path / "t.xyz"), n_frames=occ.shape[1], n_atoms=2,
            dt=1.0, dt_unit="fs", cutoff=2.0, tolerance_frames=0, cell_lengths=None,
            reference_mode="single_atom", reference_definition="1",
            observed_mode="single_atom", observed_definition="2",
            n_observed_objects=1, occupancy_mean=0.4,
            n_continuous_events=len(durations), n_censored_events=0,
            mean_residence_time=float(durations.mean()),
            mean_residence_uncertainty=0.5,
            imm_survival_time=calc.discrete_integral(S, 1.0),
            event_durations=durations, survival_time_axis=t_s, origin_survival=S,
            intermittent_mrt=1.0, intermittent_lags_used=5,
            integration_rule="zero_crossing", max_lag=10,
            time_axis=t_c, intermittent_correlation=C, intermittent_relaxation=R,
            object_labels=["Atom 2"], tolerance_scan=[],
        )

    def test_records_the_three_expected_figures(self, tmp_path, results):
        ui = self.bare_ui(tmp_path)
        ui.save_result_plots(results)
        titles = [entry.get("title", "") for entry in ui.saved_plot_data]
        assert len(titles) == 3, titles
        joined = " ".join(titles).lower()
        assert "survival" in joined
        assert "intermittent" in joined
        assert "duration" in joined

    def test_writes_no_png_files(self, tmp_path, results):
        """The project recipe: figures go through the interactive viewer only."""
        ui = self.bare_ui(tmp_path)
        ui.save_result_plots(results)
        assert glob.glob(os.path.join(str(tmp_path), "*.png")) == []

    def test_intermittent_figure_overlays_both_curves(self, tmp_path, results):
        ui = self.bare_ui(tmp_path)
        ui.save_result_plots(results)
        overlay = [e for e in ui.saved_plot_data if "series" in e]
        assert overlay, "C(t) and R(t) should be overlaid as a multi-series figure"
        assert len(overlay[0]["series"]) == 2

    def test_no_events_does_not_crash(self, tmp_path, results):
        ui = self.bare_ui(tmp_path)
        results.event_durations = np.array([])
        ui.save_result_plots(results)
        assert ui.saved_plot_data, "survival/correlation figures should still be recorded"


# --------------------------------------------------------------------------- #
# Real progress reporting                                                       #
# --------------------------------------------------------------------------- #
class TestProgressReporting:
    @staticmethod
    def tiny_trajectory(tmp_path):
        path = tmp_path / "p.xyz"
        with open(path, "w", encoding="utf-8") as f:
            for i, d in enumerate([5.0, 1.0, 1.0, 1.0, 5.0, 1.0]):
                f.write(f"2\nframe {i}\nO 0 0 0\nH {d} 0 0\n")
        return str(path)

    @staticmethod
    def params(path):
        return dict(
            xyz_file=path, dt=1.0, dt_unit="fs", cutoff=2.0, tolerance_frames=0,
            reference_mode="single_atom", reference_definition="1",
            observed_mode="single_atom", observed_definition="2",
        )

    def test_progress_is_reported_and_monotonic(self, calc, tmp_path):
        seen = []
        calc.run(**self.params(self.tiny_trajectory(tmp_path)),
                 progress=lambda frac, msg: seen.append((frac, msg)))
        assert seen, "no progress was reported"
        fractions = [f for f, _ in seen]
        assert all(0.0 <= f <= 1.0 for f in fractions), fractions
        assert fractions == sorted(fractions), f"progress went backwards: {fractions}"
        assert fractions[-1] == pytest.approx(1.0), "progress must finish at 1.0"

    def test_every_step_carries_a_message(self, calc, tmp_path):
        seen = []
        calc.run(**self.params(self.tiny_trajectory(tmp_path)),
                 progress=lambda frac, msg: seen.append((frac, msg)))
        assert all(isinstance(m, str) and m for _, m in seen)

    def test_reading_is_reported_before_the_analysis(self, calc, tmp_path):
        seen = []
        calc.run(**self.params(self.tiny_trajectory(tmp_path)),
                 progress=lambda frac, msg: seen.append(msg.lower()))
        assert "read" in seen[0], seen[0]

    def test_a_failing_callback_does_not_break_the_run(self, calc, tmp_path):
        """Progress reporting is cosmetic; it must never abort an analysis."""
        def boom(frac, msg):
            raise RuntimeError("callback exploded")
        res = calc.run(**self.params(self.tiny_trajectory(tmp_path)), progress=boom)
        assert res.n_frames == 6

    def test_progress_is_optional(self, calc, tmp_path):
        assert calc.run(**self.params(self.tiny_trajectory(tmp_path))).n_frames == 6


# --------------------------------------------------------------------------- #
# Advanced options must be reachable from the GUI                               #
# --------------------------------------------------------------------------- #
class TestAdvancedOptionsExposed:
    def test_every_run_parameter_is_supplied_by_the_ui(self):
        """Anything run() accepts should be settable in the window, otherwise
        the engine has capabilities the user cannot reach."""
        import ast
        import inspect
        from meanResidenceTime import MeanResidenceTimeCalculator

        accepted = set(inspect.signature(MeanResidenceTimeCalculator.run).parameters)
        accepted -= {"self", "progress"}          # progress is wired separately

        source = inspect.getsource(
            __import__("meanResidenceTime").MeanResidenceTimeUI._collect_inputs
        )
        tree = ast.parse(source.lstrip())
        ret = [n for n in ast.walk(tree) if isinstance(n, ast.Return)][-1]
        supplied = {k.value for k in ret.value.keys if isinstance(k, ast.Constant)}

        missing = accepted - supplied
        assert not missing, f"not reachable from the GUI: {sorted(missing)}"

    def test_integration_rules_offered_match_the_engine(self):
        from meanResidenceTime import MeanResidenceTimeUI as UI
        assert set(UI.INTEGRATION_RULES) == {"zero_crossing", "exponential", "full"}

    def test_max_lag_blank_means_automatic(self, calc):
        from meanResidenceTime import MeanResidenceTimeUI as UI
        assert UI.parse_optional_int("") is None
        assert UI.parse_optional_int("   ") is None
        assert UI.parse_optional_int("250") == 250

    @pytest.mark.parametrize("bad", ["0", "-5", "abc", "1.5"])
    def test_max_lag_rejects_nonsense(self, bad):
        from meanResidenceTime import MeanResidenceTimeUI as UI
        with pytest.raises(ValueError):
            UI.parse_optional_int(bad)


# --------------------------------------------------------------------------- #
# Dead code must stay gone                                                      #
# --------------------------------------------------------------------------- #
class TestNoDeadCode:
    @pytest.mark.parametrize("name", [
        "geometric_center",                        # superseded by _geometric_center_series
        "center_of_mass",                          # superseded by _center_of_mass_series
        "compute_continuous_survival_from_events",  # the removed duplicate metric
    ])
    def test_superseded_calculator_methods_are_removed(self, calc, name):
        assert not hasattr(calc, name), f"{name} is dead code"

    def test_show_window_is_removed(self):
        """Its only caller was open_mean_residence_time_window, already deleted."""
        from meanResidenceTime import MeanResidenceTimeUI
        assert not hasattr(MeanResidenceTimeUI, "show_window")

    def test_no_oscillating_fake_progress_loop(self):
        """The bar used to call start() (indeterminate) while also assigning
        .value in a sawtooth loop -- two contradictory modes, no real signal."""
        import inspect
        import meanResidenceTime
        source = inspect.getsource(meanResidenceTime)
        assert "progress_direction" not in source
        assert "self.progress_bar.start()" not in source


# --------------------------------------------------------------------------- #
# Parsers (unchanged behaviour, pinned so the refactor cannot drift)            #
# --------------------------------------------------------------------------- #
class TestParsers:
    @pytest.mark.parametrize("text", ["1,2,3", "1 2 3", "1, 2 3"])
    def test_atom_list_accepts_mixed_separators(self, calc, text):
        assert calc.parse_atom_list(text) == [0, 1, 2]

    # --- range syntax ----------------------------------------------------
    def test_range_expands_to_the_same_list_as_the_explicit_form(self, calc):
        """'3-7,10,15-17' must mean exactly '3,4,5,6,7,10,15,16,17'."""
        assert (calc.parse_atom_list("3-7,10,15-17")
                == calc.parse_atom_list("3,4,5,6,7,10,15,16,17"))

    def test_range_is_inclusive_and_one_based(self, calc):
        assert calc.parse_atom_list("3-7") == [2, 3, 4, 5, 6]

    def test_single_value_range_is_one_atom(self, calc):
        assert calc.parse_atom_list("5-5") == [4]

    @pytest.mark.parametrize("text", [
        "3-7,10,15-17",
        "3-7 10 15-17",
        "3-7, 10, 15-17",
        " 3-7 ,10,  15-17 ",
    ])
    def test_ranges_accept_the_same_separators_as_plain_indices(self, calc, text):
        assert calc.parse_atom_list(text) == [2, 3, 4, 5, 6, 9, 14, 15, 16]

    def test_ranges_and_singles_can_be_mixed_in_any_order(self, calc):
        assert calc.parse_atom_list("10,3-5,1") == [9, 2, 3, 4, 0]

    def test_input_order_is_preserved(self, calc):
        """MRT labels objects in the order given, so ranges must not sort."""
        assert calc.parse_atom_list("9,1,5") == [8, 0, 4]
        assert calc.parse_atom_list("7-8,1-2") == [6, 7, 0, 1]

    def test_ranges_work_inside_group_definitions(self, calc):
        assert calc.parse_group_list("1-3; 5-7") == [[0, 1, 2], [4, 5, 6]]

    def test_ranges_and_singles_mix_inside_a_group(self, calc):
        assert calc.parse_group_list("1-3,8; 5") == [[0, 1, 2, 7], [4]]

    @pytest.mark.parametrize("bad", [
        "7-3",        # descending
        "0-5",        # zero lower bound
        "3-0",        # zero upper bound
        "-5",         # missing lower bound
        "3-",         # missing upper bound
        "3-5-7",      # too many bounds
        "a-b",        # not numbers
        "3-x",        # half not a number
        "1.5-3",      # not integers
        "--",         # nonsense
    ])
    def test_malformed_ranges_are_rejected(self, calc, bad):
        with pytest.raises(ValueError):
            calc.parse_atom_list(bad)

    def test_a_bad_token_rejects_the_whole_input(self, calc):
        """A malformed range must not be silently skipped."""
        with pytest.raises(ValueError):
            calc.parse_atom_list("1,2,7-3,9")

    def test_range_error_names_the_offending_token(self, calc):
        with pytest.raises(ValueError, match="7-3"):
            calc.parse_atom_list("1,7-3")

    def test_out_of_range_expansion_is_still_validated(self, calc):
        """Ranges are checked against the trajectory like any other index."""
        with pytest.raises(ValueError, match="out of range"):
            calc._validate_indices(calc.parse_atom_list("1-50"), n_atoms=10)

    def test_large_range_expands_fully(self, calc):
        assert calc.parse_atom_list("1-280") == list(range(280))

    def test_atom_list_is_one_based(self, calc):
        assert calc.parse_atom_list("5") == [4]

    @pytest.mark.parametrize("bad", ["", "   ", "0", "-3", "x"])
    def test_atom_list_rejects_invalid(self, calc, bad):
        with pytest.raises(ValueError):
            calc.parse_atom_list(bad)

    def test_group_list_splits_on_semicolons(self, calc):
        assert calc.parse_group_list("1,2,3; 4,5,6; 10 11 12") == [
            [0, 1, 2], [3, 4, 5], [9, 10, 11]
        ]

    def test_group_list_ignores_trailing_semicolon(self, calc):
        assert calc.parse_group_list("1,2; 3,4;") == [[0, 1], [2, 3]]

    def test_out_of_range_index_is_reported_one_based(self, calc):
        with pytest.raises(ValueError, match="99"):
            calc._validate_indices([98], n_atoms=10, label="test")
