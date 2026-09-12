"""Draw docs/demo.png from the synthetic spectrum in synthetic_demo.py.

Two panels: the CH3 window sitting beside the main-chain CH2 line at two shim
qualities, and the vinyl -CH= window before and after the residual solvent
line is subtracted.

    pip install -r examples/requirements-plot.txt
    python examples/plot_demo.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "examples"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from endgroup_nmr import load_regions, subtract_solvent_lines  # noqa: E402
from synthetic_demo import build_spectrum  # noqa: E402

OUT = REPO_ROOT / "docs" / "demo.png"

GOOD = "#1f4e8c"
BAD = "#d1662a"
WINDOW = "#c8cdd6"


def slice_ppm(ppm, y, lo, hi):
    sel = (ppm >= lo) & (ppm <= hi)
    return ppm[sel], y[sel]


def panel_shim(ax, config):
    for label, hwhm, skirt, colour in (
        ("well shimmed", 0.004, 0.030, GOOD),
        ("broadened", 0.090, 0.200, BAD),
    ):
        spectrum = build_spectrum(hwhm, skirt, label)
        x, y = slice_ppm(spectrum.ppm, spectrum.intensity, 0.55, 1.75)
        y = y / y.max()
        ax.plot(x, np.clip(y, 1e-6, None), color=colour, lw=1.0, label=label)

    ch3 = config.region("CH3")
    ax.axvspan(ch3.low, ch3.high, color=WINDOW, alpha=0.55, zorder=0)
    ax.set_xlim(1.75, 0.55)
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 3.0)
    ax.set_xlabel("$\\delta$ / ppm")
    ax.set_ylabel("intensity / main-chain maximum")
    ax.set_title("(a)", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="upper left")


def panel_solvent(ax, config):
    spectrum = build_spectrum(0.004, 0.030, "well shimmed")
    corrected, _ = subtract_solvent_lines(spectrum, config)

    x, raw = slice_ppm(spectrum.ppm, spectrum.intensity, 5.55, 6.15)
    _, fixed = slice_ppm(corrected.ppm, corrected.intensity, 5.55, 6.15)

    vinyl = config.region("vinyl_CH")
    ax.axvspan(vinyl.low, vinyl.high, color=WINDOW, alpha=0.55, zorder=0)
    ax.axhline(0.0, color="#888888", lw=0.6, zorder=1)
    ax.plot(x, raw, color=BAD, lw=1.0, label="raw")
    ax.plot(x, fixed, color=GOOD, lw=1.0, label="solvent line subtracted")

    ax.set_xlim(6.15, 5.55)
    ax.set_ylim(-120, 260)
    ax.set_xlabel("$\\delta$ / ppm")
    ax.set_ylabel("intensity / a.u.")
    ax.set_title("(b)", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="upper left")


def main() -> int:
    config = load_regions(REPO_ROOT / "regions.example.json")

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6))
    panel_shim(axes[0], config)
    panel_solvent(axes[1], config)
    for ax in axes:
        ax.tick_params(labelsize=8)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.tight_layout()

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=160)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
