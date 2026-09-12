"""Readers for raw 1H NMR datasets.

Two vendor paths are supported, both through ``nmrglue``:

* Bruker / TopSpin: a dataset is a *directory* holding an ``acqus`` file plus a
  raw ``fid`` (1D) or ``ser`` (nD).  A notebook folder usually holds several
  numbered sub-experiments (``1/``, ``2/``, ``3/`` ...), each its own dataset.
  Read natively with :func:`nmrglue.bruker.read`.
* JEOL: exported as JCAMP-DX (``.jdx`` / ``.dx``).  Delta writes JCAMP-DX from
  *File > Export*; the binary ``.jdf`` container is not read here (see the
  README).

Both readers return an :class:`Acquisition`: the vendor dictionary, the data
array, and a small normalised metadata dict.  Everything downstream
(:mod:`endgroup_nmr.processing`) consumes that one shape, so adding a third
vendor means adding a reader and nothing else.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import nmrglue as ng

__all__ = [
    "Acquisition",
    "find_bruker_experiments",
    "read_bruker",
    "read_jcampdx",
    "read_any",
]


@dataclass
class Acquisition:
    """One raw dataset as it came off the spectrometer.

    Attributes
    ----------
    dic
        The vendor parameter dictionary exactly as ``nmrglue`` returned it.
    data
        Complex FID (time domain) or real/complex spectrum (frequency domain).
    vendor
        ``"bruker"`` or ``"jcampdx"``.
    is_frequency_domain
        ``True`` when ``data`` has already been Fourier transformed by the
        vendor software, so :mod:`endgroup_nmr.processing` must not FFT again.
    meta
        Normalised acquisition metadata; keys are the same for every vendor:
        ``nucleus``, ``solvent``, ``temperature_c``, ``num_scans``,
        ``pulse_program``, ``source``.
    """

    dic: Any
    data: np.ndarray | None
    vendor: str
    is_frequency_domain: bool = False
    meta: dict = field(default_factory=dict)

    @property
    def ndim(self) -> int:
        return 0 if self.data is None else int(self.data.ndim)


def find_bruker_experiments(root: str) -> list[tuple[str, str]]:
    """Walk ``root`` and return every Bruker dataset below it.

    Returns a sorted list of ``(absolute_dir, key)`` where ``key`` is the path
    relative to ``root`` with forward slashes, a stable identifier to use as a
    row label in output tables.
    """
    found: list[tuple[str, str]] = []
    for dirpath, _dirnames, filenames in os.walk(root):
        if "acqus" in filenames and ("fid" in filenames or "ser" in filenames):
            rel = os.path.relpath(dirpath, root).replace("\\", "/")
            found.append((dirpath, rel))
    found.sort(key=lambda pair: pair[1])
    return found


def _acqus(dic, key, default=""):
    try:
        return dic["acqus"][key]
    except (KeyError, TypeError):
        return default


def read_bruker(expdir: str) -> Acquisition:
    """Read one Bruker experiment directory (the one holding ``acqus``)."""
    dic, data = ng.bruker.read(expdir)

    nucleus = str(_acqus(dic, "NUC1", "")).strip() or "unknown"
    solvent = str(_acqus(dic, "SOLVENT", "")).strip()
    pulse_program = str(_acqus(dic, "PULPROG", "")).strip()
    num_scans = _acqus(dic, "NS", "")

    temperature_c = None
    te_kelvin = _acqus(dic, "TE", "")
    try:
        if te_kelvin != "":
            temperature_c = float(te_kelvin) - 273.15
    except (TypeError, ValueError):
        temperature_c = None

    # nD datasets are inventoried but not handed to the 1D pipeline.
    if data is not None and data.ndim > 1:
        data = None

    return Acquisition(
        dic=dic,
        data=data,
        vendor="bruker",
        is_frequency_domain=False,  # raw fid/ser is always time domain
        meta={
            "nucleus": nucleus,
            "solvent": solvent,
            "temperature_c": temperature_c,
            "num_scans": num_scans,
            "pulse_program": pulse_program,
            "source": expdir,
        },
    )


def _jcamp_first(dic, *keys, default=""):
    for key in keys:
        value = dic.get(key)
        if isinstance(value, list) and value:
            value = value[0]
        if value not in (None, ""):
            return value
    return default


def read_jcampdx(path: str) -> Acquisition:
    """Read a JCAMP-DX file (the JEOL Delta export path).

    JCAMP-DX files most often carry the *processed* spectrum rather than the
    FID, so :attr:`Acquisition.is_frequency_domain` is set from the file's own
    ``DATATYPE`` record instead of being assumed.
    """
    dic, data = ng.jcampdx.read(path)
    if isinstance(data, (list, tuple)):
        data = data[0]
    data = np.asarray(data)
    if data.ndim > 1:
        data = data[0] if data.shape[0] == 1 else None

    datatype = str(_jcamp_first(dic, "DATATYPE", "DATA TYPE")).upper()
    is_fid = "FID" in datatype or "FREE INDUCTION" in datatype

    nucleus = str(_jcamp_first(dic, ".OBSERVENUCLEUS", "OBSERVENUCLEUS")).strip()
    nucleus = nucleus.lstrip("^") or "unknown"

    temperature_c = None
    raw_temp = _jcamp_first(dic, ".TEMPERATURE", "TEMPERATURE", default="")
    try:
        if raw_temp != "":
            temperature_c = float(str(raw_temp).split()[0])
    except (TypeError, ValueError, IndexError):
        temperature_c = None

    return Acquisition(
        dic=dic,
        data=data,
        vendor="jcampdx",
        is_frequency_domain=not is_fid,
        meta={
            "nucleus": nucleus,
            "solvent": str(_jcamp_first(dic, ".SOLVENTNAME", "SOLVENTNAME")).strip(),
            "temperature_c": temperature_c,
            "num_scans": _jcamp_first(dic, ".AVERAGES", "AVERAGES"),
            "pulse_program": str(_jcamp_first(dic, ".PULSESEQUENCE")).strip(),
            "source": path,
        },
    )


def read_any(path: str) -> Acquisition:
    """Dispatch on what ``path`` actually is.

    A directory holding ``acqus`` is read as Bruker; a file is read as
    JCAMP-DX.  Raises ``FileNotFoundError`` / ``ValueError`` rather than
    guessing when neither applies.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    if os.path.isdir(path):
        if os.path.exists(os.path.join(path, "acqus")):
            return read_bruker(path)
        raise ValueError(
            f"{path!r} is a directory but holds no 'acqus' file; point at the "
            "numbered sub-experiment directory of a Bruker dataset"
        )
    return read_jcampdx(path)
