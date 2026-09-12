"""Fit-and-subtract for a residual solvent line that overlaps an end-group window.

Some of the most informative end-group windows sit right next to the residual
solvent line, and the line is enormous.  Its tail leaks into the window; naive
integration of the window can come out *negative*, because the leaking tail is
partly dispersive (a phasing artefact) and dispersion is negative on one side
of the line.

The fix is to model the interloper and subtract it before integrating:

* fit a phase-mixed Lorentzian, ``cos(phi) * absorption + sin(phi) *
  dispersion``, plus a local linear baseline,
* fit it on windows that *exclude* the end-group window itself (the solvent
  core on one side, a short anchor on the far side of the window), so the
  signal being measured never pulls the model,
* subtract the fitted line shape from the spectrum,
* integrate the end-group window on the corrected trace.

Fit across the window instead and the model absorbs the signal and reports that
there is nothing there.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

from .processing import Spectrum

__all__ = [
    "SolventLineFit",
    "mixed_lorentzian",
    "fit_solvent_line",
    "subtract_solvent_lines",
]


@dataclass
class SolventLineFit:
    """Fitted parameters for one interfering line."""

    center_ppm: float
    hwhm_ppm: float
    amplitude: float
    phase_rad: float
    baseline_offset: float
    baseline_slope: float
    applies_to: str = ""
    success: bool = True
    message: str = ""

    @property
    def params(self) -> np.ndarray:
        return np.array(
            [
                self.center_ppm,
                self.hwhm_ppm,
                self.amplitude,
                self.phase_rad,
                self.baseline_offset,
                self.baseline_slope,
            ],
            dtype=float,
        )

    def as_row(self) -> dict:
        return {
            "solvent_fit_center_ppm": self.center_ppm,
            "solvent_fit_hwhm_ppm": self.hwhm_ppm,
            "solvent_fit_phase_rad": self.phase_rad,
            "solvent_fit_ok": self.success,
        }


def mixed_lorentzian(params, x: np.ndarray, include_baseline: bool = True) -> np.ndarray:
    """``[center, hwhm, amplitude, phase, offset, slope]`` evaluated on ``x``.

    ``phase = 0`` is pure absorption, ``phase = pi/2`` pure dispersion.
    """
    center, hwhm, amplitude, phase, offset, slope = params
    delta = np.asarray(x, dtype=float) - center
    denom = hwhm * hwhm + delta * delta
    absorption = hwhm / denom
    dispersion = delta / denom
    line = amplitude * (np.cos(phase) * absorption + np.sin(phase) * dispersion)
    if include_baseline:
        line = line + offset + slope * np.asarray(x, dtype=float)
    return line


def fit_solvent_line(
    spectrum: Spectrum,
    fit_windows,
    center_bounds: tuple[float, float],
    max_hwhm_ppm: float = 0.2,
    applies_to: str = "",
    max_nfev: int = 20000,
) -> SolventLineFit:
    """Least-squares fit of one phase-mixed Lorentzian.

    Parameters
    ----------
    fit_windows
        Sequence of ``(low_ppm, high_ppm)`` pairs to fit on.  Give at least the
        solvent core plus one anchor on the far side of the window you are
        protecting, and do **not** include that window.
    center_bounds
        Hard ``(low, high)`` bounds on the line position, so the fit cannot
        wander onto a different peak.
    """
    mask = np.zeros_like(spectrum.ppm, dtype=bool)
    for low, high in fit_windows:
        mask |= spectrum.mask(low, high)
    x = spectrum.ppm[mask]
    y = spectrum.intensity[mask]
    if x.size < 8:
        return SolventLineFit(
            center_ppm=float(np.mean(center_bounds)),
            hwhm_ppm=float("nan"),
            amplitude=0.0,
            phase_rad=0.0,
            baseline_offset=0.0,
            baseline_slope=0.0,
            applies_to=applies_to,
            success=False,
            message=f"only {x.size} points inside the fit windows; nothing subtracted",
        )

    i0 = int(np.argmax(y))
    start = [
        float(np.clip(x[i0], *center_bounds)),
        0.01,
        float(abs(y[i0])) * 0.01,
        0.0,
        0.0,
        0.0,
    ]
    lower = [center_bounds[0], 1e-4, 0.0, -1.2, -np.inf, -np.inf]
    upper = [center_bounds[1], max_hwhm_ppm, np.inf, 1.2, np.inf, np.inf]
    start = [float(np.clip(v, lo, hi)) for v, lo, hi in zip(start, lower, upper)]

    result = least_squares(
        lambda p: mixed_lorentzian(p, x) - y,
        start,
        bounds=(lower, upper),
        max_nfev=max_nfev,
    )
    center, hwhm, amplitude, phase, offset, slope = (float(v) for v in result.x)
    return SolventLineFit(
        center_ppm=center,
        hwhm_ppm=hwhm,
        amplitude=amplitude,
        phase_rad=phase,
        baseline_offset=offset,
        baseline_slope=slope,
        applies_to=applies_to,
        success=bool(result.success),
        message=str(result.message),
    )


def subtract_solvent_lines(
    spectrum: Spectrum,
    config,
    include_baseline: bool = False,
) -> tuple[Spectrum, list[SolventLineFit]]:
    """Fit and subtract every line declared in the ``solvent_lines`` block.

    ``include_baseline=False`` (the default) subtracts only the line shape.
    The local linear term is part of the *fit*, where it soaks up whatever
    offset remains in the fit windows, but it is not extrapolated across the
    whole spectrum, which would tilt regions far from the fit.

    Returns the corrected spectrum and the fits, in declaration order.  A
    config with no ``solvent_lines`` block returns the spectrum unchanged.
    """
    lines = config.solvent_lines or []
    if not lines:
        return spectrum, []

    corrected = np.array(spectrum.intensity, dtype=float)
    fits: list[SolventLineFit] = []
    for line in lines:
        fit = fit_solvent_line(
            spectrum,
            [tuple(w) for w in line["fit_windows"]],
            tuple(line["center_bounds"]),
            max_hwhm_ppm=line.get("max_hwhm_ppm", 0.2),
            applies_to=line.get("applies_to", ""),
        )
        fits.append(fit)
        if fit.success and np.isfinite(fit.hwhm_ppm):
            corrected -= mixed_lorentzian(
                fit.params, spectrum.ppm, include_baseline=include_baseline
            )
    return spectrum.with_intensity(corrected), fits
