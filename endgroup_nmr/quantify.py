"""Region integration, end-group bookkeeping, ends/chain and Mn.

Everything spectral is behind :func:`integrate_region`; everything chemical is
driven by the regions file, which the code reads and never second-guesses.  A
region declares

* the ppm window,
* which chemical **group** it reports on,
* how many protons of that group fall inside the window,
* optionally how many carbons the group contributes to the chain,
* and its **role**: ``backbone``, ``chain_end`` or ``in_chain``.

From there the arithmetic is short:

    amount(group)  = integral / protons_per_group      (relative molar amount)
    chains         = sum of amount over chain_count_groups
    ends_per_chain = sum of amount over chain_end_groups / chains
    Mn             = repeat_unit_mass * carbons_per_chain + end_group_mass

``ends_per_chain`` needs an independent count of chains, and there is no way
to get one from the spectrum alone -- you have to assume something.  The
default assumption in the shipped regions file is the textbook one for a
polyolefin terminated by beta-hydride elimination: **exactly one unsaturated
terminus per chain**.  A perfectly linear chain then gives 2.0 (one saturated
end plus one unsaturated end), and the departure from 2.0 is the diagnostic.
Change ``chain_count_groups`` if your chemistry counts chains some other way;
the number means nothing without knowing which assumption produced it.
"""

from __future__ import annotations

import json
import statistics
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .processing import Spectrum

__all__ = [
    "Region",
    "RegionConfig",
    "load_regions",
    "integrate_region",
    "integrate_regions",
    "group_amounts",
    "EndGroupResult",
    "quantify",
]


@dataclass(frozen=True)
class Region:
    """One integration window declared in the regions file."""

    name: str
    low: float
    high: float
    group: str
    protons_per_group: float = 1.0
    carbons_per_group: float = 0.0
    role: str = "chain_end"
    note: str = ""


@dataclass
class RegionConfig:
    """Parsed regions file."""

    nucleus: str = "1H"
    reference: dict = field(default_factory=dict)
    regions: list[Region] = field(default_factory=list)
    quantification: dict = field(default_factory=dict)
    qc: dict = field(default_factory=dict)
    solvent_lines: list = field(default_factory=list)
    source: str = ""

    def region(self, name: str) -> Region:
        for region in self.regions:
            if region.name == name:
                return region
        raise KeyError(f"no region named {name!r} in {self.source or 'the regions config'}")

    def regions_for_group(self, group: str) -> list[Region]:
        return [r for r in self.regions if r.group == group]

    def groups_with_role(self, role: str) -> list[str]:
        seen: list[str] = []
        for region in self.regions:
            if region.role == role and region.group not in seen:
                seen.append(region.group)
        return seen

    @property
    def backbone_group(self) -> str:
        declared = self.quantification.get("backbone_group")
        if declared:
            return declared
        backbone = self.groups_with_role("backbone")
        if not backbone:
            raise ValueError("the regions config declares no region with role 'backbone'")
        return backbone[0]

    @property
    def reference_ppm(self) -> float | None:
        value = self.reference.get("ppm")
        return None if value is None else float(value)


def load_regions(path: str | Path) -> RegionConfig:
    """Read a regions JSON file.

    The file is the only place ppm windows live.  Nothing in this package has
    a hard-coded chemical shift.
    """
    path = Path(path)
    with path.open("r", encoding="utf-8") as handle:
        raw = json.load(handle)

    regions = []
    for item in raw.get("regions", []):
        regions.append(
            Region(
                name=item["name"],
                low=float(item["low"]),
                high=float(item["high"]),
                group=item.get("group", item["name"]),
                protons_per_group=float(item.get("protons_per_group", 1.0)),
                carbons_per_group=float(item.get("carbons_per_group", 0.0)),
                role=item.get("role", "chain_end"),
                note=item.get("note", ""),
            )
        )
    if not regions:
        raise ValueError(f"{path} declares no regions")

    return RegionConfig(
        nucleus=raw.get("nucleus", "1H"),
        reference=raw.get("reference", {}),
        regions=regions,
        quantification=raw.get("quantification", {}),
        qc=raw.get("qc", {}),
        solvent_lines=raw.get("solvent_lines", []),
        source=str(path),
    )


def integrate_region(spectrum: Spectrum, low: float, high: float) -> float:
    """Trapezoidal integral of the spectrum over ``[low, high]`` ppm.

    The ppm axis usually descends, so the points are sorted before
    integrating: the result is positive for a positive peak whichever way the
    axis runs.  An empty window integrates to ``0.0`` rather than raising --
    a window that falls off the edge of the spectrum contains no signal.
    """
    x, y = spectrum.slice(low, high)
    if x.size < 2:
        return 0.0
    return float(np.trapezoid(y, x))


def integrate_regions(spectrum: Spectrum, config: RegionConfig) -> dict[str, float]:
    """Integrate every declared region.  Keys are region names."""
    return {r.name: integrate_region(spectrum, r.low, r.high) for r in config.regions}


def group_amounts(integrals: dict[str, float], config: RegionConfig) -> dict[str, float]:
    """Relative molar amount of each chemical group.

    A group observed through more than one window (a vinyl end shows up as
    both ``-CH=`` and ``=CH2``) is reduced with the **median** of the
    per-window amounts, which is insensitive to one window being spoiled by an
    overlap.
    """
    amounts: dict[str, float] = {}
    buckets: dict[str, list[float]] = {}
    for region in config.regions:
        value = integrals.get(region.name, 0.0) / region.protons_per_group
        buckets.setdefault(region.group, []).append(value)
    for group, values in buckets.items():
        amounts[group] = float(statistics.median(values))
    return amounts


@dataclass
class EndGroupResult:
    """Everything one spectrum yields."""

    integrals: dict[str, float]
    amounts: dict[str, float]
    per_1000_backbone_protons: dict[str, float]
    chains: float
    ends_per_chain: float
    backbone_units_per_chain: float
    mn: float
    assumptions: dict = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    def as_row(self) -> dict:
        """Flat dict for a results table."""
        row: dict = {}
        for name, value in self.per_1000_backbone_protons.items():
            row[f"{name}_per1000H"] = value
        row["ends_per_chain"] = self.ends_per_chain
        row["backbone_units_per_chain"] = self.backbone_units_per_chain
        row["Mn_g_per_mol"] = self.mn
        return row


def quantify(spectrum: Spectrum, config: RegionConfig) -> EndGroupResult:
    """Integrate, normalise, and turn the integrals into ends/chain and Mn."""
    integrals = integrate_regions(spectrum, config)
    amounts = group_amounts(integrals, config)

    backbone_group = config.backbone_group
    backbone_regions = config.regions_for_group(backbone_group)
    backbone_integral = sum(integrals[r.name] for r in backbone_regions)

    per_1000 = {
        name: (1000.0 * value / backbone_integral if backbone_integral > 0 else float("nan"))
        for name, value in integrals.items()
    }

    quant = config.quantification
    chain_count_groups = quant.get("chain_count_groups") or config.groups_with_role("chain_end")
    end_groups = quant.get("chain_end_groups") or config.groups_with_role("chain_end")
    repeat_mass = float(quant.get("repeat_unit_mass", 14.0266))
    extra_end_mass = float(quant.get("end_group_mass", 0.0))

    notes: list[str] = []
    chains = float(sum(max(amounts.get(g, 0.0), 0.0) for g in chain_count_groups))
    ends_total = float(sum(max(amounts.get(g, 0.0), 0.0) for g in end_groups))

    backbone_units = amounts.get(backbone_group, 0.0)
    # one carbon count per group, not per window: a group seen through two
    # windows must not contribute its carbons twice.
    carbons_per_group: dict[str, float] = {}
    for region in config.regions:
        carbons_per_group[region.group] = max(
            carbons_per_group.get(region.group, 0.0), region.carbons_per_group
        )
    carbons_total = float(
        sum(amounts.get(group, 0.0) * carbons for group, carbons in carbons_per_group.items())
    )

    if chains <= 0:
        notes.append(
            "no chain-counting group was detected, so ends/chain and Mn are undefined; "
            "either the spectrum has no measurable end groups or the windows are wrong"
        )
        ends_per_chain = float("nan")
        units_per_chain = float("nan")
        mn = float("nan")
    else:
        ends_per_chain = ends_total / chains
        units_per_chain = backbone_units / chains
        mn = repeat_mass * carbons_total / chains + extra_end_mass

    if backbone_integral <= 0:
        notes.append("the backbone window integrates to zero or less; check phase and baseline")

    return EndGroupResult(
        integrals=integrals,
        amounts=amounts,
        per_1000_backbone_protons=per_1000,
        chains=chains,
        ends_per_chain=ends_per_chain,
        backbone_units_per_chain=units_per_chain,
        mn=mn,
        assumptions={
            "backbone_group": backbone_group,
            "chain_count_groups": list(chain_count_groups),
            "chain_end_groups": list(end_groups),
            "repeat_unit_mass": repeat_mass,
            "end_group_mass": extra_end_mass,
        },
        notes=notes,
    )
