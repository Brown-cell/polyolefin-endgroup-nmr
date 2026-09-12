"""Unit tests on synthetic spectra with known answers.

Everything here is built from Lorentzians whose areas are chosen in advance,
so each assertion compares a measured number against arithmetic rather than
against a previous run of the same code.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from endgroup_nmr import (  # noqa: E402
    Spectrum,
    calibrate_ppm,
    check_valley_ratio,
    fit_solvent_line,
    group_amounts,
    integrate_region,
    load_regions,
    quantify,
    subtract_solvent_lines,
    valley_ratio,
)

REGIONS_FILE = REPO_ROOT / "regions.example.json"

# Relative molar amounts for the reference sample used by most tests.
AMOUNTS = {
    "backbone": 50.0,
    "methyl_end": 1.30,
    "vinyl_end": 0.60,
    "vinylidene_end": 0.40,
    "internal_olefin": 0.20,
}


def lorentzian(x, center, hwhm, area, phase=0.0):
    delta = np.asarray(x) - center
    denom = delta * delta + hwhm * hwhm
    return (area / np.pi) * (
        np.cos(phase) * hwhm / denom + np.sin(phase) * delta / denom
    )


def axis(n=120001, high=8.0, low=-1.0):
    """Descending ppm axis, as the vendors write it."""
    return np.linspace(high, low, n)


def clean_spectrum(main_chain_hwhm=0.001, solvent_area=0.0, solvent_phase=0.0):
    """Well-resolved synthetic polyolefin spectrum built from AMOUNTS."""
    ppm = axis()
    y = np.zeros_like(ppm)
    y += lorentzian(ppm, 1.30, main_chain_hwhm, 2.0 * AMOUNTS["backbone"])
    y += lorentzian(ppm, 0.88, 0.001, 3.0 * AMOUNTS["methyl_end"])
    y += lorentzian(ppm, 4.70, 0.001, 2.0 * AMOUNTS["vinylidene_end"])
    y += lorentzian(ppm, 4.94, 0.001, 2.0 * AMOUNTS["vinyl_end"])
    y += lorentzian(ppm, 5.79, 0.001, 1.0 * AMOUNTS["vinyl_end"])
    y += lorentzian(ppm, 5.40, 0.001, 2.0 * AMOUNTS["internal_olefin"])
    if solvent_area:
        y += lorentzian(ppm, 5.91, 0.004, solvent_area, phase=solvent_phase)
    return Spectrum(ppm, y)


@pytest.fixture(scope="module")
def config():
    return load_regions(REGIONS_FILE)


# --------------------------------------------------------------- windows ---

def test_regions_file_declares_the_expected_windows(config):
    names = {r.name for r in config.regions}
    assert names == {
        "main_chain_CH2",
        "CH3",
        "vinylidene_CH2",
        "vinyl_CH2",
        "vinyl_CH",
        "internal_olefin",
    }
    assert config.backbone_group == "backbone"
    assert config.region("main_chain_CH2").protons_per_group == 2


# ------------------------------------------------------------ integration ---

def test_integrate_region_recovers_a_known_area():
    ppm = axis()
    y = lorentzian(ppm, 1.30, 0.002, 12.0)
    spectrum = Spectrum(ppm, y)
    # +/- 0.15 ppm is 75 half-widths, so >99.5% of a Lorentzian's area
    assert integrate_region(spectrum, 1.15, 1.45) == pytest.approx(12.0, rel=0.01)


def test_integration_does_not_depend_on_axis_direction():
    ppm = axis()
    y = lorentzian(ppm, 1.30, 0.002, 12.0)
    descending = Spectrum(ppm, y)
    ascending = Spectrum(ppm[::-1], y[::-1])
    assert integrate_region(descending, 1.15, 1.45) == pytest.approx(
        integrate_region(ascending, 1.15, 1.45), rel=1e-9
    )


def test_window_outside_the_spectrum_integrates_to_zero():
    spectrum = clean_spectrum()
    assert integrate_region(spectrum, 20.0, 21.0) == 0.0


def test_group_amounts_takes_the_median_over_a_group_two_windows(config):
    # vinyl_end is seen through vinyl_CH (1H) and vinyl_CH2 (2H); both must
    # report the same molar amount.
    integrals = {"vinyl_CH": 0.6, "vinyl_CH2": 1.2, "main_chain_CH2": 100.0,
                 "CH3": 3.9, "vinylidene_CH2": 0.8, "internal_olefin": 0.4}
    amounts = group_amounts(integrals, config)
    assert amounts["vinyl_end"] == pytest.approx(0.6)
    assert amounts["backbone"] == pytest.approx(50.0)
    assert amounts["methyl_end"] == pytest.approx(1.3)


# ------------------------------------------------------------ valley_ratio ---

def test_valley_ratio_is_small_when_the_peaks_are_resolved(config):
    qc = check_valley_ratio(clean_spectrum(main_chain_hwhm=0.001), config)
    assert qc is not None
    assert qc.ratio < 1e-3
    assert qc.passed is True
    # the trough must sit between the two maxima
    assert 0.88 < qc.valley_ppm < 1.30


def test_valley_ratio_grows_when_the_main_chain_line_broadens(config):
    narrow = check_valley_ratio(clean_spectrum(main_chain_hwhm=0.001), config)
    broad = check_valley_ratio(clean_spectrum(main_chain_hwhm=0.09), config)
    assert broad.ratio > 20 * narrow.ratio
    assert broad.passed is False


def test_valley_ratio_without_a_threshold_reports_no_verdict():
    spectrum = clean_spectrum()
    result = valley_ratio(spectrum, (1.10, 1.45), (0.75, 1.00))
    assert result.threshold is None
    assert result.passed is None
    assert np.isfinite(result.ratio)


def test_broadening_inflates_the_CH3_window_which_is_what_the_gate_guards(config):
    narrow = quantify(clean_spectrum(main_chain_hwhm=0.001), config)
    broad = quantify(clean_spectrum(main_chain_hwhm=0.09), config)
    assert (
        broad.per_1000_backbone_protons["CH3"]
        > 2 * narrow.per_1000_backbone_protons["CH3"]
    )
    assert broad.ends_per_chain > narrow.ends_per_chain


# ------------------------------------------------------------- quantities ---

def test_ends_per_chain_matches_the_amounts_the_spectrum_was_built_from(config):
    result = quantify(clean_spectrum(), config)
    chains = AMOUNTS["vinyl_end"] + AMOUNTS["vinylidene_end"]
    expected = (
        AMOUNTS["methyl_end"] + AMOUNTS["vinyl_end"] + AMOUNTS["vinylidene_end"]
    ) / chains
    assert result.ends_per_chain == pytest.approx(expected, rel=0.02)


def test_backbone_units_and_mn_match_the_construction(config):
    result = quantify(clean_spectrum(), config)
    chains = AMOUNTS["vinyl_end"] + AMOUNTS["vinylidene_end"]
    carbons = (
        AMOUNTS["backbone"]
        + AMOUNTS["methyl_end"]
        + 2 * AMOUNTS["vinyl_end"]
        + 2 * AMOUNTS["vinylidene_end"]
        + 2 * AMOUNTS["internal_olefin"]
    )
    assert result.backbone_units_per_chain == pytest.approx(
        AMOUNTS["backbone"] / chains, rel=0.02
    )
    assert result.mn == pytest.approx(14.0266 * carbons / chains, rel=0.02)


def test_internal_olefin_counts_in_mn_but_not_in_ends_per_chain(config):
    result = quantify(clean_spectrum(), config)
    assert "internal_olefin" not in result.assumptions["chain_end_groups"]
    assert result.amounts["internal_olefin"] == pytest.approx(
        AMOUNTS["internal_olefin"], rel=0.05
    )


def test_tail_only_end_group_windows_give_a_large_number_not_nan(config):
    """Documented behaviour, and a trap worth knowing about.

    NaN is returned only when the chain-counting windows integrate to zero or
    less.  A single huge main-chain line puts a little Lorentzian tail into
    every window, so ``chains`` comes out tiny-but-positive and ends/chain
    comes out large and finite rather than undefined.  There is no detection
    floor in this package: judge the result on ``valley_ratio`` and on the
    ``chains`` figure itself.
    """
    ppm = axis()
    y = lorentzian(ppm, 1.30, 0.002, 100.0)
    result = quantify(Spectrum(ppm, y), config)
    assert result.chains > 0
    assert np.isfinite(result.ends_per_chain)
    assert result.ends_per_chain > 10  # nonsense, and visibly so


def test_a_spectrum_with_no_end_groups_gives_nan_rather_than_a_number(config):
    ppm = axis()
    y = lorentzian(ppm, 1.30, 0.002, 100.0)
    y[np.abs(ppm - 1.30) > 0.1] = 0.0  # main chain only, nothing anywhere else
    result = quantify(Spectrum(ppm, y), config)
    assert np.isnan(result.ends_per_chain)
    assert np.isnan(result.mn)
    assert result.notes


# ------------------------------------------------------------ calibration ---

def test_calibrate_ppm_moves_a_known_line_onto_the_reference():
    ppm = axis()
    y = lorentzian(ppm, 5.85, 0.004, 100.0)  # line is 0.06 ppm off
    shifted, offset = calibrate_ppm(ppm, y, reference_ppm=5.91, search_halfwidth=0.3)
    assert offset == pytest.approx(-0.06, abs=1e-3)
    assert shifted[int(np.argmax(y))] == pytest.approx(5.91, abs=1e-3)


def test_calibration_is_skipped_rather_than_guessed_when_nothing_is_in_range():
    ppm = axis()
    y = lorentzian(ppm, 1.30, 0.004, 100.0)
    shifted, offset = calibrate_ppm(ppm, y, reference_ppm=40.0, search_halfwidth=0.3)
    assert offset is None
    assert np.array_equal(shifted, ppm)


# --------------------------------------------------------- solvent fitting ---

def test_solvent_fit_finds_the_line_it_was_pointed_at():
    spectrum = clean_spectrum(solvent_area=60.0, solvent_phase=0.35)
    fit = fit_solvent_line(
        spectrum,
        [(5.90, 6.40), (5.55, 5.68)],
        center_bounds=(5.85, 6.00),
    )
    assert fit.success
    assert fit.center_ppm == pytest.approx(5.91, abs=0.005)
    assert fit.hwhm_ppm == pytest.approx(0.004, rel=0.15)
    assert fit.phase_rad == pytest.approx(0.35, abs=0.05)


def test_subtracting_the_solvent_line_restores_the_overlapped_window(config):
    without = quantify(clean_spectrum(), config)
    with_solvent = clean_spectrum(solvent_area=60.0, solvent_phase=0.35)

    uncorrected = quantify(with_solvent, config)
    corrected_spectrum, fits = subtract_solvent_lines(with_solvent, config)
    corrected = quantify(corrected_spectrum, config)

    # the dispersive tail drives the raw window negative ...
    assert uncorrected.integrals["vinyl_CH"] < 0
    # ... and subtraction brings it back to the clean value
    assert corrected.integrals["vinyl_CH"] == pytest.approx(
        without.integrals["vinyl_CH"], rel=0.05
    )
    assert fits[0].applies_to == "vinyl_CH"


def test_a_config_without_solvent_lines_returns_the_spectrum_untouched():
    bare = load_regions(REGIONS_FILE)  # a fresh copy: never mutate the shared fixture
    bare.solvent_lines = []
    spectrum = clean_spectrum()
    corrected, fits = subtract_solvent_lines(spectrum, bare)
    assert fits == []
    assert corrected is spectrum
