"""UV-Vis spectrum of the excited states (NTO tab's Plot Spectrum): math, manifest, drawing, replot."""
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import excitationData as ed
import uvVisSpectrum as uv
from nto_fixtures import fixture
from plotter import PlotterBase

SRC = Path(__file__).resolve().parent.parent


def state(number, energy_ev, f, label="Singlet-A"):
    return ed.ExcitedState(number, label, 1, energy_ev, round(uv.NM_EV / energy_ev, 2), f, 0.0, ())


@pytest.fixture(scope="module")
def uracil_entry():
    data = ed.read_excitation_file(fixture("uracil")[1])
    return uv.spectrum_entry(data.states, source={"fchk": "uracil.fch", "excitation_file": "uracil.out"})


class TestBandShape:
    def test_gaussview_half_width_is_gaussian_sigma_0_4_ev(self):
        assert uv.sigma_from_hwhm(uv.DEFAULT_HWHM_EV) == pytest.approx(0.4, abs=1e-3)

    def test_peak_height_and_half_width(self):
        e0, f = 6.0, 0.1
        sigma = uv.sigma_from_hwhm(0.333)
        peak = uv.epsilon(np.array([uv.NM_EV / e0]), [e0], [f], sigma)[0]
        assert peak == pytest.approx(uv.EPSILON_PREFACTOR * f / (sigma * uv.EV_TO_WAVENUMBER), rel=1e-12)
        for e in (e0 - 0.333, e0 + 0.333):  # the half width is in energy, on either side
            assert uv.epsilon(np.array([uv.NM_EV / e]), [e0], [f], sigma)[0] == pytest.approx(peak / 2, rel=1e-9)

    def test_bands_add(self):
        lam = np.linspace(150, 300, 50)
        a = uv.epsilon(lam, [5.0], [0.2], 0.4)
        b = uv.epsilon(lam, [6.5], [0.05], 0.4)
        assert np.allclose(uv.epsilon(lam, [5.0, 6.5], [0.2, 0.05], 0.4), a + b)

    def test_range_covers_every_band(self):
        lo, hi = uv.wavelength_range([4.0, 7.0], 0.4)
        assert lo == pytest.approx(uv.NM_EV / (7.0 + 4 * 0.4)) and hi == pytest.approx(uv.NM_EV / (4.0 - 4 * 0.4))
        lo, hi = uv.wavelength_range([0.5], 0.4)  # never down to zero energy
        assert 0 < lo < uv.NM_EV / 0.5 < hi < math.inf

    def test_a_dark_state_only_needs_room_for_its_stick(self):
        lo, hi = uv.wavelength_range([3.35, 9.0], 0.4, strengths=[0.0, 0.2])
        assert hi == pytest.approx(uv.NM_EV / (3.35 - 0.4))     # 1 sigma beyond the dark triplet
        assert lo == pytest.approx(uv.NM_EV / (9.0 + 4 * 0.4))  # the bright band in full
        lo, hi = uv.wavelength_range([5.0], 0.4, strengths=[0.0])  # all dark: still a range
        assert lo < uv.NM_EV / 5.0 < hi


class TestEntry:
    def test_every_state_with_its_identifier_wavelength_and_strength(self, uracil_entry):
        e = uracil_entry
        assert e["type"] == "uv_vis" and e["title"] == "UV-Vis Spectrum"
        assert [s["number"] for s in e["states"]] == [1, 2, 3]
        assert e["states"][2] == {"number": 3, "label": 'Singlet-A"', "energy_ev": 6.018,
                                  "wavelength_nm": 206.02, "oscillator": 0.0}
        assert e["states"][1]["oscillator"] == 0.1317
        assert e["xlabel"] == "Wavelength (nm)" and e["ylabel_right"] == "Oscillator Strength"

    def test_settings_and_curve(self, uracil_entry):
        b = uracil_entry["broadening"]
        assert b == {"shape": "gaussian", "hwhm_ev": 0.333, "sigma_ev": pytest.approx(uv.sigma_from_hwhm(0.333)),
                     "prefactor": uv.EPSILON_PREFACTOR, "x": "wavenumber"}
        c = uracil_entry["curve"]
        assert len(c["wavelength_nm"]) == len(c["epsilon"]) == uv.CURVE_POINTS
        lo, hi = uracil_entry["xlim"]
        assert lo <= min(c["wavelength_nm"]) and max(c["wavelength_nm"]) <= hi
        assert all(lo < s["wavelength_nm"] < hi for s in uracil_entry["states"])
        # the curve is the band sum at its own wavelengths
        states = uracil_entry["states"]
        expected = uv.epsilon(np.array(c["wavelength_nm"]), [s["energy_ev"] for s in states],
                              [s["oscillator"] for s in states], b["sigma_ev"])
        assert np.allclose(c["epsilon"], expected)
        assert uracil_entry["source"] == {"fchk": "uracil.fch", "excitation_file": "uracil.out"}

    def test_triplets_are_kept_with_zero_strength(self):
        data = ed.read_excitation_file(fixture("h2co_td5050")[1])
        entry = uv.spectrum_entry(data.states)
        assert len(entry["states"]) == 8
        assert [s["label"].split("-")[0] for s in entry["states"]][:2] == ["Triplet", "Singlet"]
        assert entry["states"][0]["oscillator"] == 0.0

    def test_a_custom_half_width(self):
        entry = uv.spectrum_entry([state(1, 5.0, 0.1)], hwhm_ev=0.2)
        assert entry["broadening"]["hwhm_ev"] == 0.2
        assert max(entry["curve"]["epsilon"]) == pytest.approx(
            uv.EPSILON_PREFACTOR * 0.1 / (uv.sigma_from_hwhm(0.2) * uv.EV_TO_WAVENUMBER), rel=1e-4)

    @pytest.mark.parametrize("hwhm", [0.0, -0.1, float("nan")])
    def test_bad_half_width(self, hwhm):
        with pytest.raises(ValueError, match="half-width"):
            uv.spectrum_entry([state(1, 5.0, 0.1)], hwhm_ev=hwhm)

    def test_states_without_oscillator_strengths_are_refused(self):
        bare = ed.ExcitedState(1, "Singlet", 1, 5.0, None, None, None, ())
        with pytest.raises(ValueError, match="oscillator strengths"):
            uv.spectrum_entry([bare])
        with pytest.raises(ValueError, match="no excited states"):
            uv.spectrum_entry([])

    def test_write_and_reload_with_plotter(self, uracil_entry, tmp_path):
        path = uv.write_spectrum(uracil_entry, tmp_path / uv.spectrum_filename("uracil"))
        assert Path(path).name == "uracil_UVVis.json"
        assert json.loads(Path(path).read_text(encoding="utf-8")) == [json.loads(json.dumps(uracil_entry))]
        assert PlotterBase.load_json_plot_file(path)[0]["type"] == "uv_vis"
        with pytest.raises(ValueError, match="Could not write"):
            uv.write_spectrum(uracil_entry, tmp_path)  # a folder


class TestDrawing:
    def test_gaussview_layout(self, uracil_entry):
        import matplotlib.pyplot as plt
        from displayPlots import draw_uv_vis

        figure = plt.figure()
        axes = figure.gca()
        twin = draw_uv_vis(axes, uracil_entry)
        assert len(figure.axes) == 2 and twin is figure.axes[1]
        assert axes.get_title() == "UV-Vis Spectrum" and axes.get_xlabel() == "Wavelength (nm)"
        assert twin.get_ylabel() == "Oscillator Strength"
        assert len(axes.lines) == 1                              # the epsilon curve
        sticks = twin.collections[0].get_segments()
        assert len(sticks) == 3 and sticks[1][1][1] == pytest.approx(0.1317)
        assert [t.get_text() for t in twin.texts] == ['1: Singlet-A"', "2: Singlet-A'", '3: Singlet-A"']
        assert axes.get_xlim() == tuple(uracil_entry["xlim"])
        assert axes.get_ylim()[0] == 0 and twin.get_ylim()[0] == 0
        assert axes.get_facecolor()[:3] == pytest.approx((1.0, 1.0, 0.94), abs=0.01)
        plt.close(figure)

    def test_all_dark_states(self):
        import matplotlib.pyplot as plt
        from displayPlots import draw_uv_vis

        entry = uv.spectrum_entry([state(1, 4.0, 0.0, "Triplet-A")])
        figure = plt.figure()
        twin = draw_uv_vis(figure.gca(), entry)
        assert twin.get_ylim()[1] > 0 and figure.gca().get_ylim()[1] > 0  # no singular axes
        plt.close(figure)

    def test_viewer_renders_it(self, uracil_entry):
        import matplotlib.pyplot as plt
        import plotViewer

        plt.close("all")
        plotViewer.build_figures([json.loads(json.dumps(uracil_entry))])
        assert plt.figure(1).axes[0].get_title() == "UV-Vis Spectrum" and len(plt.figure(1).axes) == 2
        plt.close("all")


class TestPlotterManifests:
    def test_accepted_and_described(self, uracil_entry):
        figures = PlotterBase.validate_plot_manifest([json.loads(json.dumps(uracil_entry))])
        assert "UV-Vis Spectrum - 3 excited states" in PlotterBase.describe_json_figures(figures)

    @pytest.mark.parametrize("change, message", [
        (lambda f: f.pop("states"), "'states'"),
        (lambda f: f.update(states=[]), "'states'"),
        (lambda f: f["states"][1].update(oscillator="big"), "Figure 1, state 2"),
        (lambda f: f["states"][0].pop("wavelength_nm"), "Figure 1, state 1"),
        (lambda f: f["curve"]["epsilon"].pop(), "same length"),
        (lambda f: f.pop("curve"), "'curve'"),
    ])
    def test_rejected_with_a_message_naming_the_figure(self, uracil_entry, change, message):
        fig = json.loads(json.dumps(uracil_entry))
        change(fig)
        with pytest.raises(ValueError) as info:
            PlotterBase.validate_plot_manifest([fig])
        assert message in str(info.value)


def test_module_does_not_import_toga_or_matplotlib():
    code = "import sys, uvVisSpectrum; print('toga' in sys.modules, 'matplotlib' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], cwd=SRC, capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "False False"
