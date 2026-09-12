"""Command-line front end: point it at datasets, get a table.

    python -m endgroup_nmr.cli --regions regions.example.json DATA_DIR
    python -m endgroup_nmr.cli --regions regions.example.json spectrum.jdx --csv out.csv

A path that is a directory is searched recursively for Bruker datasets (any
directory holding ``acqus`` next to ``fid`` or ``ser``), so you can hand it a
whole day's folder.  A path that is a file is read as JCAMP-DX.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from .processing import process
from .qc import check_valley_ratio
from .quantify import load_regions, quantify
from .readers import find_bruker_experiments, read_any, read_bruker
from .solvent_fit import subtract_solvent_lines

__all__ = ["analyze", "collect_datasets", "main"]


def collect_datasets(path: str) -> list[tuple[str, str, str]]:
    """Expand one command-line path into ``(key, path, kind)`` triples."""
    p = Path(path)
    if p.is_file():
        return [(p.name, str(p), "file")]
    if not p.is_dir():
        raise FileNotFoundError(path)
    if (p / "acqus").exists():
        return [(p.name, str(p), "bruker")]
    found = find_bruker_experiments(str(p))
    return [(key, expdir, "bruker") for expdir, key in found]


def analyze(dataset_path: str, config, kind: str = "auto") -> dict:
    """Read, process, QC, solvent-correct and quantify one dataset.

    Returns a flat row: metadata, QC, then the quantities.  Quantification
    runs even when the QC gate fails (the row carries ``qc_passed`` so the
    caller can drop it), because seeing the number you are rejecting is more
    useful than a blank.
    """
    acq = read_bruker(dataset_path) if kind == "bruker" else read_any(dataset_path)
    spectrum = process(acq, reference_ppm=config.reference_ppm,
                       search_halfwidth=float(config.reference.get("search_halfwidth_ppm", 0.6)))

    row: dict = {
        "nucleus": acq.meta.get("nucleus", ""),
        "solvent": acq.meta.get("solvent", ""),
        "temperature_c": acq.meta.get("temperature_c"),
        "num_scans": acq.meta.get("num_scans"),
        "calibration_shift_ppm": spectrum.meta.get("calibration_shift_ppm"),
    }

    qc = check_valley_ratio(spectrum, config)
    if qc is not None:
        row.update(qc.as_row())
        row["qc_passed"] = qc.passed

    corrected, fits = subtract_solvent_lines(spectrum, config)
    for fit in fits:
        row.update({f"{k}[{fit.applies_to}]" if fit.applies_to else k: v
                    for k, v in fit.as_row().items()})

    result = quantify(corrected, config)
    row.update(result.as_row())
    if result.notes:
        row["notes"] = "; ".join(result.notes)
    return row


def _print_table(rows: list[dict], columns: list[str]) -> None:
    widths = {c: max(len(c), *(len(_fmt(r.get(c))) for r in rows)) for c in columns}
    print("  ".join(c.rjust(widths[c]) for c in columns))
    print("  ".join("-" * widths[c] for c in columns))
    for row in rows:
        print("  ".join(_fmt(row.get(c)).rjust(widths[c]) for c in columns))


def _fmt(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "yes" if value else "NO"
    if isinstance(value, float):
        if value != value:  # NaN
            return "-"
        return f"{value:.4g}"
    return str(value)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="endgroup_nmr",
        description="Quantify polyolefin chain-end groups from 1H NMR.",
    )
    parser.add_argument("paths", nargs="+", help="Bruker dataset directories, or JCAMP-DX files")
    parser.add_argument("--regions", required=True, help="path to the regions JSON file")
    parser.add_argument("--csv", help="also write the table to this CSV file")
    parser.add_argument(
        "--drop-failed-qc",
        action="store_true",
        help="omit spectra that fail the valley_ratio gate instead of flagging them",
    )
    args = parser.parse_args(argv)

    config = load_regions(args.regions)

    rows: list[dict] = []
    for path in args.paths:
        for key, dataset_path, kind in collect_datasets(path):
            try:
                row = analyze(dataset_path, config, kind=kind)
            except Exception as exc:  # one bad dataset must not kill the batch
                rows.append({"key": key, "notes": f"FAILED: {exc}"})
                print(f"failed: {key}: {exc}", file=sys.stderr)
                continue
            row = {"key": key, **row}
            if args.drop_failed_qc and row.get("qc_passed") is False:
                print(f"rejected by valley_ratio: {key}", file=sys.stderr)
                continue
            rows.append(row)

    if not rows:
        print("no datasets found", file=sys.stderr)
        return 1

    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)

    summary = [c for c in ("key", "solvent", "valley_ratio", "qc_passed",
                           "ends_per_chain", "Mn_g_per_mol", "notes") if c in columns]
    _print_table(rows, summary)

    if args.csv:
        with open(args.csv, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nwrote {args.csv} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
