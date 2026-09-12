"""Quantify polyolefin chain-end groups from 1H NMR.

Read a raw dataset (a Bruker directory, or a JCAMP-DX export from JEOL Delta),
work it up, gate it on shim quality, subtract an overlapping residual solvent
line, and turn region integrals into ends per chain and a number-average molar
mass.

Typical use::

    from endgroup_nmr import load_regions, read_any, process, quantify
    from endgroup_nmr import check_valley_ratio, subtract_solvent_lines

    config = load_regions("regions.example.json")
    spectrum = process(read_any(path), reference_ppm=config.reference_ppm)

    qc = check_valley_ratio(spectrum, config)
    if qc.passed is False:
        raise SystemExit(f"shim quality too poor: valley_ratio={qc.ratio:.3f}")

    corrected, fits = subtract_solvent_lines(spectrum, config)
    result = quantify(corrected, config)
    print(result.ends_per_chain, result.mn)

Every chemical shift lives in the regions file; there is no hard-coded ppm
value in this package.
"""

from .processing import (
    Spectrum,
    autophase,
    calibrate_ppm,
    correct_baseline,
    fourier_transform,
    ppm_axis,
    process,
)
from .qc import ValleyRatio, check_valley_ratio, valley_ratio
from .quantify import (
    EndGroupResult,
    Region,
    RegionConfig,
    group_amounts,
    integrate_region,
    integrate_regions,
    load_regions,
    quantify,
)
from .readers import (
    Acquisition,
    find_bruker_experiments,
    read_any,
    read_bruker,
    read_jcampdx,
)
from .solvent_fit import (
    SolventLineFit,
    fit_solvent_line,
    mixed_lorentzian,
    subtract_solvent_lines,
)

__version__ = "0.1.0"

__all__ = [
    "Acquisition",
    "EndGroupResult",
    "Region",
    "RegionConfig",
    "SolventLineFit",
    "Spectrum",
    "ValleyRatio",
    "autophase",
    "calibrate_ppm",
    "check_valley_ratio",
    "correct_baseline",
    "find_bruker_experiments",
    "fit_solvent_line",
    "fourier_transform",
    "group_amounts",
    "integrate_region",
    "integrate_regions",
    "load_regions",
    "mixed_lorentzian",
    "ppm_axis",
    "process",
    "quantify",
    "read_any",
    "read_bruker",
    "read_jcampdx",
    "subtract_solvent_lines",
    "valley_ratio",
    "__version__",
]
