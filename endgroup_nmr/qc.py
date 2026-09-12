"""Shim-quality gate: the ``valley_ratio``.

The main-chain CH2 line is three or four orders of magnitude taller than every
signal worth measuring.  A Lorentzian has heavy tails and a poorly shimmed line
heavier ones, so the foot of the main-chain peak spills into the neighbouring
CH3 window and is integrated as chain ends.  The spectrum looks fine; the
number is inflated.

``valley_ratio`` measures that spill:

    valley_ratio = (intensity at the lowest point between the CH3 maximum and
                    the main-chain maximum) / (main-chain maximum)

On a well-shimmed spectrum the peaks resolve down to the noise and the ratio is
small; as the line broadens, the trough lifts off the baseline and the ratio
grows.  Rejecting spectra above a threshold throws away the ones whose CH3
integral is contaminated by main-chain tail, a failure no later arithmetic can
undo.

The threshold is empirical and belongs in the regions file, since it depends on
nucleus, field, temperature, solvent and how much CH3 you expect.  Measure a
handful of spectra you trust and a handful you know are bad, then put the
number where it separates them.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .processing import Spectrum

__all__ = ["ValleyRatio", "valley_ratio", "check_valley_ratio"]


@dataclass
class ValleyRatio:
    """Result of one shim-quality measurement.

    ``passed`` is ``None`` when no threshold was supplied, and ``None`` also
    when the ratio could not be measured, since an unmeasurable ratio is not a
    pass.
    """

    ratio: float
    peak_ppm: float
    peak_height: float
    valley_ppm: float
    valley_height: float
    window_peak_ppm: float
    window_peak_height: float
    threshold: float | None = None
    passed: bool | None = None
    note: str = ""

    def as_row(self) -> dict:
        """Flat dict for a results table."""
        return {
            "valley_ratio": self.ratio,
            "valley_ratio_threshold": self.threshold,
            "valley_ratio_passed": self.passed,
            "main_peak_ppm": self.peak_ppm,
            "valley_ppm": self.valley_ppm,
            "qc_note": self.note,
        }


def _argmax_in(spectrum: Spectrum, low: float, high: float) -> tuple[float, float]:
    mask = spectrum.mask(low, high)
    if not mask.any():
        raise ValueError(f"no spectrum points inside [{low}, {high}] ppm")
    ppm = spectrum.ppm[mask]
    y = spectrum.intensity[mask]
    i = int(np.argmax(y))
    return float(ppm[i]), float(y[i])


def valley_ratio(
    spectrum: Spectrum,
    main_chain_region: tuple[float, float],
    endgroup_region: tuple[float, float],
    threshold: float | None = None,
) -> ValleyRatio:
    """Measure the trough between the main-chain peak and the CH3 peak.

    Parameters
    ----------
    spectrum
        Phased, baseline-corrected, calibrated spectrum.
    main_chain_region, endgroup_region
        ``(low_ppm, high_ppm)`` for the tall reference line and for the window
        it might be leaking into.  With the shipped regions file these are the
        main-chain CH2 and the CH3 windows.
    threshold
        Optional maximum acceptable ratio.  When given, :attr:`ValleyRatio.passed`
        is set.

    Notes
    -----
    The trough is searched *between the two maxima*, not between the two window
    edges, so the measurement does not move when you widen a window by a few
    hundredths of a ppm.
    """
    peak_ppm, peak_height = _argmax_in(spectrum, *main_chain_region)
    win_ppm, win_height = _argmax_in(spectrum, *endgroup_region)

    low, high = sorted((peak_ppm, win_ppm))
    between = spectrum.mask(low, high)
    n_between = int(between.sum())

    note = ""
    if peak_height <= 0:
        ratio = float("nan")
        valley_ppm = float("nan")
        valley_height = float("nan")
        note = "main-chain peak height is not positive; check phase and baseline"
    elif n_between < 3:
        ratio = float("nan")
        valley_ppm = float("nan")
        valley_height = float("nan")
        note = (
            "the two maxima are less than three points apart; the windows "
            "probably overlap or the spectrum is badly under-sampled"
        )
    else:
        ppm_between = spectrum.ppm[between]
        y_between = spectrum.intensity[between]
        i = int(np.argmin(y_between))
        valley_ppm = float(ppm_between[i])
        valley_height = float(y_between[i])
        ratio = valley_height / peak_height

    passed: bool | None = None
    if threshold is not None and np.isfinite(ratio):
        passed = bool(ratio <= threshold)

    return ValleyRatio(
        ratio=ratio,
        peak_ppm=peak_ppm,
        peak_height=peak_height,
        valley_ppm=valley_ppm,
        valley_height=valley_height,
        window_peak_ppm=win_ppm,
        window_peak_height=win_height,
        threshold=threshold,
        passed=passed,
        note=note,
    )


def check_valley_ratio(spectrum: Spectrum, config) -> ValleyRatio | None:
    """Run :func:`valley_ratio` using the ``qc`` block of a regions config.

    Returns ``None`` when the config declares no ``valley_ratio`` block, so a
    minimal regions file stays usable.
    """
    block = (config.qc or {}).get("valley_ratio")
    if not block:
        return None
    main = config.region(block["main_chain_region"])
    ends = config.region(block["endgroup_region"])
    return valley_ratio(
        spectrum,
        (main.low, main.high),
        (ends.low, ends.high),
        threshold=block.get("max_ratio"),
    )
