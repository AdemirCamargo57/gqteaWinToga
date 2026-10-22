"""UV-Vis absorption spectrum of a set of excited states (the NTO tab's Plot Spectrum).

Pure numpy, no Toga, no matplotlib: it builds a ``{"type": "uv_vis"}`` plot
manifest entry that ``displayPlots.draw_uv_vis`` draws in the interactive plot
viewer, and that Plotter's JSON plot type reopens later.

The broadened curve is GaussView's / Gaussian's ("UV/Visible Plots in
Gaussian"): one Gaussian band per state in wavenumbers,

    eps(nu) = sum_i 1.3062974e8 * f_i / sigma * exp(-((nu - nu_i) / sigma)^2)

with nu, nu_i and sigma in cm^-1 and eps in L mol^-1 cm^-1. GaussView's
"peak half-width at half height" is HWHM = sigma * sqrt(ln 2); its default
0.333 eV is Gaussian's sigma = 0.4 eV.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import numpy as np

EPSILON_PREFACTOR = 1.3062974e8   # L mol^-1 cm^-2 per unit oscillator strength (Gaussian's constant)
NM_EV = 1239.841984               # h c in eV nm (CODATA 2018): lambda (nm) = NM_EV / E (eV)
EV_TO_WAVENUMBER = 1e7 / NM_EV    # cm^-1 per eV, 8065.544; derived so E, lambda and nu agree exactly
DEFAULT_HWHM_EV = 0.333           # GaussView's default peak half-width at half height
CURVE_POINTS = 2000
RANGE_SIGMAS = 4.0                # the plotted range covers every band +/- this many sigma
DARK_MARGIN_SIGMAS = 1.0          # ... and every dark (f = 0) state's stick +/- this many
SPECTRUM_SUFFIX = "_UVVis.json"
TITLE = "UV-Vis Spectrum"


def sigma_from_hwhm(hwhm_ev: float) -> float:
    """Gaussian width parameter (eV) of a band whose half width at half height is ``hwhm_ev``."""
    return hwhm_ev / math.sqrt(math.log(2.0))


def epsilon(wavelength_nm: np.ndarray, energies_ev: Sequence[float], strengths: Sequence[float],
            sigma_ev: float) -> np.ndarray:
    """Molar absorptivity (L mol^-1 cm^-1) at each wavelength: the sum of every state's band."""
    wavenumber = 1e7 / np.asarray(wavelength_nm, dtype=float)                    # cm^-1
    centres = np.asarray(energies_ev, dtype=float) * EV_TO_WAVENUMBER
    sigma = sigma_ev * EV_TO_WAVENUMBER
    f = np.asarray(strengths, dtype=float)
    return (EPSILON_PREFACTOR * f / sigma
            * np.exp(-((wavenumber[:, None] - centres[None, :]) / sigma) ** 2)).sum(axis=1)


def wavelength_range(energies_ev: Sequence[float], sigma_ev: float,
                     strengths: Optional[Sequence[float]] = None) -> Tuple[float, float]:
    """(shortest, longest) wavelength (nm) to plot (never down to 0 eV).

    Every band with f > 0 is covered +/- RANGE_SIGMAS sigma; a dark state
    (f = 0, e.g. a triplet) has no band, so it only needs room for its stick
    (+/- DARK_MARGIN_SIGMAS sigma) and cannot stretch the axis over empty space.
    Without ``strengths`` every state counts as bright.
    """
    energies = list(energies_ev)
    if strengths is None:
        bright = energies
    else:
        bright = [e for e, f in zip(energies, strengths) if f > 0]
    margins = [(e, RANGE_SIGMAS) for e in bright] + [(e, DARK_MARGIN_SIGMAS) for e in energies]
    high = max(e + k * sigma_ev for e, k in margins)
    low = min(e - k * sigma_ev for e, k in margins)
    if low <= 0.5 * min(energies):
        low = 0.5 * min(energies)
    return NM_EV / high, NM_EV / low


def spectrum_entry(states, hwhm_ev: float = DEFAULT_HWHM_EV, title: str = TITLE,
                   source: Optional[Dict[str, str]] = None) -> dict:
    """The plot manifest entry for ``states`` (excitationData.ExcitedState): data, settings and curve."""
    if not (isinstance(hwhm_ev, (int, float)) and math.isfinite(hwhm_ev) and hwhm_ev > 0):
        raise ValueError(f"The peak half-width must be a positive number of eV, not {hwhm_ev}.")
    states = list(states)
    if not states:
        raise ValueError("There are no excited states to plot.")
    missing = [s.number for s in states if s.oscillator is None]
    if missing:
        raise ValueError(f"The excitation file gives no oscillator strengths (states "
                         f"{', '.join(map(str, missing[:5]))}{'...' if len(missing) > 5 else ''}); "
                         "plot the spectrum from the Gaussian output instead of a plain-text file.")
    sigma = sigma_from_hwhm(float(hwhm_ev))
    energies = [s.energy_ev for s in states]
    strengths = [s.oscillator for s in states]
    lo, hi = wavelength_range(energies, sigma, strengths)
    grid = np.linspace(lo, hi, CURVE_POINTS)
    return {
        "type": "uv_vis",
        "title": title,
        "xlabel": "Wavelength (nm)",
        "ylabel": "ε (L mol⁻¹ cm⁻¹)",
        "ylabel_right": "Oscillator Strength",
        "states": [{"number": s.number, "label": s.label, "energy_ev": s.energy_ev,
                    "wavelength_nm": s.wavelength_nm if s.wavelength_nm is not None else NM_EV / s.energy_ev,
                    "oscillator": s.oscillator} for s in states],
        "broadening": {"shape": "gaussian", "hwhm_ev": float(hwhm_ev), "sigma_ev": sigma,
                       "prefactor": EPSILON_PREFACTOR, "x": "wavenumber"},
        "curve": {"wavelength_nm": grid.tolist(), "epsilon": epsilon(grid, energies, strengths, sigma).tolist()},
        "xlim": [float(lo), float(hi)],
        "source": dict(source or {}),
    }


def spectrum_filename(stem: str) -> str:
    return stem + SPECTRUM_SUFFIX


def write_spectrum(entry: dict, path) -> str:
    """Write ``[entry]`` as a plot manifest (Plotter's JSON plot type reopens it); OSError -> ValueError."""
    path = Path(path)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps([entry], indent=1), encoding="utf-8")
    except OSError as exc:
        raise ValueError(f"Could not write {path}: {exc}") from exc
    return str(path)
