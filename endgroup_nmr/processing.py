"""From a raw acquisition to a phased, baseline-corrected, ppm-calibrated trace.

The 1D workup is the standard nmrglue recipe, in this order:

1. **remove the digital filter** (Bruker only -- ``grpdly`` group delay, which
   otherwise wraps the first points of the FID and ruins the baseline),
2. **zero-fill** to the next power of two, then double it,
3. **FFT**,
4. **reverse the axis** -- see :func:`fourier_transform` for why this step is
   not optional on Bruker data,
5. **automatic phase correction** (ACME entropy minimisation),
6. **discard the imaginary channel**,
7. **polynomial baseline correction**,
8. **ppm calibration** against a line of known shift (usually the residual
   solvent line).

Every step is a separate function so a caller who has already phased their
data in the vendor software can skip straight to :func:`calibrate_ppm`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import nmrglue as ng

from .readers import Acquisition

__all__ = [
    "Spectrum",
    "fourier_transform",
    "autophase",
    "correct_baseline",
    "ppm_axis",
    "calibrate_ppm",
    "process",
]


@dataclass
class Spectrum:
    """A real 1D spectrum on a ppm axis.

    ``ppm`` is stored exactly as the vendor conventions produce it, i.e.
    **descending** (left edge = high ppm).  Nothing in this package assumes a
    direction: every consumer selects points with a ``low <= ppm <= high``
    mask, and :func:`endgroup_nmr.quantify.integrate_region` sorts before it
    integrates.
    """

    ppm: np.ndarray
    intensity: np.ndarray
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.ppm = np.asarray(self.ppm, dtype=float)
        self.intensity = np.asarray(self.intensity, dtype=float)
        if self.ppm.shape != self.intensity.shape:
            raise ValueError(
                f"ppm and intensity must have the same shape, got "
                f"{self.ppm.shape} and {self.intensity.shape}"
            )

    def mask(self, low: float, high: float) -> np.ndarray:
        """Boolean mask for the closed ppm interval ``[low, high]``."""
        if high < low:
            low, high = high, low
        return (self.ppm >= low) & (self.ppm <= high)

    def slice(self, low: float, high: float) -> tuple[np.ndarray, np.ndarray]:
        """``(ppm, intensity)`` inside ``[low, high]``, sorted by ascending ppm."""
        m = self.mask(low, high)
        x, y = self.ppm[m], self.intensity[m]
        order = np.argsort(x)
        return x[order], y[order]

    def with_intensity(self, intensity: np.ndarray) -> "Spectrum":
        """A copy carrying a different intensity array (same axis and metadata)."""
        return Spectrum(self.ppm.copy(), np.asarray(intensity, dtype=float), dict(self.meta))


def fourier_transform(acq: Acquisition, zero_fill: bool = True) -> np.ndarray:
    """FID -> complex spectrum, including the axis reversal.

    **Why the reversal step exists.**  ``nmrglue.proc_base.fft`` applies the
    plain numpy FFT, which orders the output by increasing frequency index.
    Bruker digitises with the opposite sense, so the transform comes out as a
    *mirror image* of the spectrum: the aromatic region lands where the
    aliphatic region belongs.  It is a genuinely nasty trap, because a mirrored
    polyolefin spectrum still looks entirely plausible -- one tall aliphatic
    peak with small satellites -- and every integral you take from it is wrong
    while nothing raises an error.  ``proc_base.rev`` puts the axis back.
    """
    if acq.data is None:
        raise ValueError("acquisition carries no 1D data to transform")
    if acq.is_frequency_domain:
        return np.asarray(acq.data)

    fid = np.asarray(acq.data)
    if acq.vendor == "bruker":
        fid = ng.bruker.remove_digital_filter(acq.dic, fid)
    if zero_fill:
        n = fid.shape[-1]
        fid = ng.proc_base.zf_size(fid, 2 ** int(np.ceil(np.log2(n))) * 2)
    spec = ng.proc_base.fft(fid)
    spec = ng.proc_base.rev(spec)
    return spec


def autophase(spec: np.ndarray, algorithm: str = "acme") -> np.ndarray:
    """Automatic zero- and first-order phase correction, then take the real part."""
    phased = ng.proc_autophase.autops(spec, algorithm, disp=False)
    return ng.proc_base.di(phased)


def correct_baseline(intensity: np.ndarray, window: int = 20) -> np.ndarray:
    """Polynomial baseline correction on the real spectrum."""
    return ng.proc_bl.baseline_corrector(np.asarray(intensity), wd=window)


def ppm_axis(acq: Acquisition, intensity: np.ndarray) -> np.ndarray:
    """Chemical-shift axis for ``intensity``, from the vendor parameters."""
    if acq.vendor == "bruker":
        udic = ng.bruker.guess_udic(acq.dic, intensity)
    else:
        udic = ng.jcampdx.guess_udic(acq.dic, intensity)
    uc = ng.fileiobase.uc_from_udic(udic, dim=0)
    return uc.ppm_scale()


def calibrate_ppm(
    ppm: np.ndarray,
    intensity: np.ndarray,
    reference_ppm: float,
    search_halfwidth: float = 0.6,
) -> tuple[np.ndarray, float | None]:
    """Shift the axis so the tallest line near ``reference_ppm`` sits on it.

    The reference is normally the residual solvent line, whose shift is
    tabulated.  Returns ``(shifted_ppm, applied_shift)``; ``applied_shift`` is
    ``None`` when no point of the axis falls inside the search window, in
    which case the axis is returned untouched rather than silently moved.
    """
    ppm = np.asarray(ppm, dtype=float)
    intensity = np.asarray(intensity, dtype=float)
    window = (ppm > reference_ppm - search_halfwidth) & (ppm < reference_ppm + search_halfwidth)
    if not window.any():
        return ppm, None
    found = float(ppm[window][int(np.argmax(np.abs(intensity[window])))])
    offset = found - reference_ppm
    return ppm - offset, offset


def process(
    acq: Acquisition,
    reference_ppm: float | None = None,
    search_halfwidth: float = 0.6,
    baseline_window: int = 20,
    phase_algorithm: str = "acme",
) -> Spectrum:
    """Run the whole workup and return a :class:`Spectrum`.

    ``reference_ppm=None`` skips calibration and keeps the vendor's own
    referencing.
    """
    spec = fourier_transform(acq)
    real = autophase(spec, phase_algorithm) if np.iscomplexobj(spec) else np.asarray(spec).real
    real = correct_baseline(real, window=baseline_window)
    ppm = ppm_axis(acq, real)

    applied = None
    if reference_ppm is not None:
        ppm, applied = calibrate_ppm(ppm, real, reference_ppm, search_halfwidth)

    meta = dict(acq.meta)
    meta["calibration_shift_ppm"] = applied
    meta["vendor"] = acq.vendor
    return Spectrum(ppm, real, meta)
