"""End-to-end demo on a synthetic spectrum, with no instrument data needed.

Builds a 1H spectrum of a linear polyolefin from Lorentzians at the textbook
shifts, with the features that make real end-group analysis awkward:

* a main-chain CH2 line ~5000x taller than any end-group signal, sitting on a
  broad skirt,
* a residual solvent line overlapping the vinyl -CH= window, phased so its tail
  is partly dispersive (the reason naive integration of that window can go
  negative),
* noise.

It then runs the library on two copies, one well shimmed and one with a
broadened main-chain line, and prints what the valley_ratio gate, the solvent
subtraction and the quantification make of each.  The spectrum was built from
known amounts, so the last table compares the recovered numbers against them.

Run it from the repository root:

    python examples/synthetic_demo.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from endgroup_nmr import (  # noqa: E402
    Spectrum,
    check_valley_ratio,
    load_regions,
    quantify,
    subtract_solvent_lines,
)

# ---------------------------------------------------------------- truth ----
# Relative molar amounts the synthetic sample is built from.  One "chain"
# carries one unsaturated terminus, so chains = vinyl + vinylidene = 1.0.
TRUE_AMOUNTS = {
    "backbone": 200.0,        # CH2 units per chain
    "methyl_end": 1.30,
    "vinyl_end": 0.60,
    "vinylidene_end": 0.40,
    "internal_olefin": 0.20,
}
SOLVENT_AREA = 60.0           # residual solvent line, in the same units
SOLVENT_PPM = 5.91
SOLVENT_PHASE_RAD = 0.35      # nonzero -> the tail is partly dispersive
NOISE_SIGMA = 0.05
RNG_SEED = 20260912


def lorentzian(x, center, hwhm, area, phase=0.0):
    """Area-normalised Lorentzian, optionally phase-mixed with its dispersion."""
    delta = x - center
    denom = delta * delta + hwhm * hwhm
    absorption = hwhm / denom
    dispersion = delta / denom
    return (area / np.pi) * (np.cos(phase) * absorption + np.sin(phase) * dispersion)


def build_spectrum(main_chain_hwhm: float, skirt_hwhm: float, label: str) -> Spectrum:
    """Synthesise one 1H spectrum on a descending ppm axis."""
    rng = np.random.default_rng(RNG_SEED)
    ppm = np.linspace(10.0, -1.0, 32768)  # descending, like the vendors write it
    y = np.zeros_like(ppm)

    backbone_area = 2.0 * TRUE_AMOUNTS["backbone"]      # 2 protons per CH2
    y += lorentzian(ppm, 1.30, main_chain_hwhm, backbone_area * 0.94)
    y += lorentzian(ppm, 1.30, skirt_hwhm, backbone_area * 0.06)

    y += lorentzian(ppm, 0.88, 0.004, 3.0 * TRUE_AMOUNTS["methyl_end"])
    y += lorentzian(ppm, 4.70, 0.005, 2.0 * TRUE_AMOUNTS["vinylidene_end"])
    y += lorentzian(ppm, 4.94, 0.005, 2.0 * TRUE_AMOUNTS["vinyl_end"])
    y += lorentzian(ppm, 5.79, 0.006, 1.0 * TRUE_AMOUNTS["vinyl_end"])
    y += lorentzian(ppm, 5.40, 0.006, 2.0 * TRUE_AMOUNTS["internal_olefin"])

    y += lorentzian(ppm, SOLVENT_PPM, 0.004, SOLVENT_AREA, phase=SOLVENT_PHASE_RAD)
    y += rng.normal(0.0, NOISE_SIGMA, ppm.size)

    return Spectrum(ppm, y, {"label": label})


def fmt(value, digits=3):
    if value is None:
        return "-"
    if isinstance(value, bool):
        return "pass" if value else "REJECT"
    if isinstance(value, float) and not np.isfinite(value):
        return "-"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def print_table(title, columns, rows):
    print(f"\n{title}")
    widths = [max(len(c), *(len(r[i]) for r in rows)) for i, c in enumerate(columns)]
    print("  ".join(c.ljust(w) if i == 0 else c.rjust(w)
                    for i, (c, w) in enumerate(zip(columns, widths))))
    print("  ".join("-" * w for w in widths))
    for row in rows:
        print("  ".join(cell.ljust(w) if i == 0 else cell.rjust(w)
                        for i, (cell, w) in enumerate(zip(row, widths))))


def main() -> int:
    config = load_regions(REPO_ROOT / "regions.example.json")

    cases = [
        ("well shimmed", build_spectrum(0.004, 0.030, "well shimmed")),
        ("poorly shimmed", build_spectrum(0.090, 0.200, "poorly shimmed")),
    ]

    qc_rows = []
    quant_rows = []
    solvent_rows = []
    recovered: dict[str, object] = {}

    for name, spectrum in cases:
        qc = check_valley_ratio(spectrum, config)
        qc_rows.append([
            name,
            fmt(qc.ratio, 4),
            fmt(qc.threshold, 3),
            fmt(qc.passed),
            fmt(qc.valley_ppm, 2),
            f"{qc.peak_height:.4g}",
        ])

        raw = quantify(spectrum, config)
        corrected, fits = subtract_solvent_lines(spectrum, config)
        fitted = quantify(corrected, config)
        fit = fits[0]

        solvent_rows.append([
            name,
            fmt(fit.center_ppm, 3),
            fmt(fit.hwhm_ppm, 4),
            fmt(fit.phase_rad, 3),
            fmt(raw.per_1000_backbone_protons["vinyl_CH"], 3),
            fmt(fitted.per_1000_backbone_protons["vinyl_CH"], 3),
        ])

        quant_rows.append([
            name,
            fmt(fitted.per_1000_backbone_protons["CH3"], 2),
            fmt(fitted.amounts["vinyl_end"] / fitted.chains, 3),
            fmt(fitted.amounts["vinylidene_end"] / fitted.chains, 3),
            fmt(fitted.ends_per_chain, 3),
            fmt(fitted.backbone_units_per_chain, 1),
            f"{fitted.mn:,.0f}",
        ])
        recovered[name] = fitted

    print("Synthetic 1H spectrum of a linear polyolefin (no instrument data involved).")
    print(f"Axis 10.0 to -1.0 ppm, 32768 points, noise sigma {NOISE_SIGMA}, seed {RNG_SEED}.")

    print_table(
        "1. Shim-quality gate (valley_ratio between the CH3 and main-chain maxima)",
        ["case", "valley_ratio", "max", "verdict", "valley/ppm", "peak height"],
        qc_rows,
    )

    print_table(
        "2. Residual solvent line at 5.91 ppm, fitted and subtracted",
        ["case", "fit centre", "fit hwhm", "fit phase", "-CH= raw", "-CH= corrected"],
        solvent_rows,
    )
    print("   (-CH= columns are the 5.70-5.88 window per 1000 main-chain protons.)")

    print_table(
        "3. Quantification after solvent subtraction",
        ["case", "CH3/1000H", "vinyl/chain", "vinylidene/chain",
         "ends/chain", "CH2 units/chain", "Mn"],
        quant_rows,
    )

    truth_carbons = (
        TRUE_AMOUNTS["backbone"] * 1
        + TRUE_AMOUNTS["methyl_end"] * 1
        + TRUE_AMOUNTS["vinyl_end"] * 2
        + TRUE_AMOUNTS["vinylidene_end"] * 2
        + TRUE_AMOUNTS["internal_olefin"] * 2
    )
    true_chains = TRUE_AMOUNTS["vinyl_end"] + TRUE_AMOUNTS["vinylidene_end"]
    true_ends = (
        TRUE_AMOUNTS["methyl_end"] + TRUE_AMOUNTS["vinyl_end"] + TRUE_AMOUNTS["vinylidene_end"]
    )
    backbone_protons = 2.0 * TRUE_AMOUNTS["backbone"]
    truth = {
        "CH3 / 1000 main-chain H": 3.0 * TRUE_AMOUNTS["methyl_end"] / backbone_protons * 1000,
        "-CH= / 1000 main-chain H": 1.0 * TRUE_AMOUNTS["vinyl_end"] / backbone_protons * 1000,
        "vinyl ends / chain": TRUE_AMOUNTS["vinyl_end"] / true_chains,
        "ends / chain": true_ends / true_chains,
        "CH2 units / chain": TRUE_AMOUNTS["backbone"] / true_chains,
        "Mn": 14.0266 * truth_carbons / true_chains,
    }

    def recovered_values(result):
        return {
            "CH3 / 1000 main-chain H": result.per_1000_backbone_protons["CH3"],
            "-CH= / 1000 main-chain H": result.per_1000_backbone_protons["vinyl_CH"],
            "vinyl ends / chain": result.amounts["vinyl_end"] / result.chains,
            "ends / chain": result.ends_per_chain,
            "CH2 units / chain": result.backbone_units_per_chain,
            "Mn": result.mn,
        }

    good = recovered_values(recovered["well shimmed"])
    bad = recovered_values(recovered["poorly shimmed"])
    truth_rows = [
        [label, fmt(truth[label], 2), fmt(good[label], 2), fmt(bad[label], 2)]
        for label in truth
    ]
    print_table(
        "4. Recovered vs. the amounts the spectrum was built from",
        ["quantity", "truth", "well shimmed", "poorly shimmed"],
        truth_rows,
    )

    print(
        "\n   Read it this way. The solvent-overlapped -CH= window and the two unsaturated\n"
        "   end counts come back on the nose once the solvent line is subtracted. The CH3\n"
        "   window does not: even in the well-shimmed case it reads high, because a\n"
        "   Lorentzian main-chain line 5000x taller than the CH3 signal still has\n"
        "   measurable tail under 0.75-1.00 ppm, and that tail is integrated as if it were\n"
        "   chain ends. ends/chain inherits the same bias, always upwards.\n"
        "   That is the point of the gate rather than an argument against it: valley_ratio\n"
        "   measures the size of exactly this leak, and the two cases differ by more than\n"
        "   two orders of magnitude in it. Below the threshold the CH3 window reads high by\n"
        "   tens of per cent and ends/chain with it, while Mn and the chain length (which\n"
        "   rest on the main-chain and unsaturation windows, not on CH3) stay within a\n"
        "   couple of per cent. Above the threshold the CH3 window is mostly main-chain\n"
        "   tail and every quantity is out by roughly a factor of three, chain length and\n"
        "   Mn included. So the gate does not certify CH3 as unbiased; it separates a bias\n"
        "   you can state and live with from one that swamps the measurement.\n"
        "   A real line shape is closer to Gaussian in the far wings than the Lorentzian\n"
        "   used here, so this demo is the pessimistic case for CH3 leakage."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
