# polyolefin-endgroup-nmr

Quantify polyolefin chain-end groups from ¹H NMR: region integration against
the main-chain CH₂ reference, a shim-quality gate, fit-and-subtract for a
residual solvent line that overlaps an end-group window, and ends per chain
plus a number-average molar mass.

Raw data is read with [nmrglue](https://www.nmrglue.com/): Bruker/TopSpin
dataset directories natively, JEOL through a JCAMP-DX export.

```
pip install -r requirements.txt
python examples/synthetic_demo.py          # runs on a synthetic spectrum, no data needed
python -m endgroup_nmr.cli --regions regions.example.json /path/to/data
```

## What it measures

A linear polyethylene chain is a very long run of CH₂ with something different
at each end. Integrate the ends against the backbone and you get the chain
length, and therefore Mₙ, without a calibration curve, plus a count of which
kind of end you have.

The windows below are the textbook assignments and live in
`regions.example.json`, not in the code:

| window | ppm | what it reports |
|---|---|---|
| main-chain CH₂ | 1.10–1.45 | the `(CH₂)ₙ` envelope; the intensity reference |
| CH₃ | 0.75–1.00 | saturated ends (and short-chain branches) |
| vinylidene =CH₂ | 4.60–4.80 | 1,1-disubstituted olefin end |
| vinyl =CH₂ | 4.90–5.00 | terminal vinyl, two protons |
| vinyl –CH= | 5.70–5.88 | the same vinyl end, one proton |
| internal –CH=CH– | 5.30–5.50 | in-chain unsaturation, not an end |

From the integrals:

```
amount(group)  = integral / protons_per_group
chains         = sum of amount over the chain-counting groups
ends_per_chain = sum of amount over the chain-end groups / chains
Mn             = repeat_unit_mass * (carbons per chain) + end_group_mass
```

`ends_per_chain` needs a count of chains, and no spectrum gives you one for
free — you have to assume something. The shipped default is the textbook
assumption for a chain terminated by β-hydride elimination: **exactly one
unsaturated terminus per chain**. A perfectly linear chain then reads 2.0 (one
saturated end, one unsaturated end) and the departure from 2.0 is the finding.
Change `chain_count_groups` if your chemistry counts chains differently; the
number is meaningless without knowing which assumption produced it.

## The reversal step, and why leaving it out is so dangerous

`nmrglue.proc_base.fft` is the plain FFT, which orders its output by
increasing frequency index. Bruker digitises with the opposite sense, so the
transform comes out as a **mirror image**: the aromatic end of the spectrum
lands where the aliphatic end belongs. `proc_base.rev` puts it back, and
`processing.fourier_transform` always calls it.

The reason to be careful about this one is that a mirrored polyolefin spectrum
still looks completely plausible — one enormous aliphatic line with small
satellites around it, exactly what you expected to see. Nothing raises an
error, no peak is missing, and every integral you take is from the wrong part
of the spectrum. Check the residual solvent line lands where it should before
believing anything else.

The rest of the workup is the standard recipe: remove the Bruker digital
filter (the `grpdly` group delay, which otherwise wrecks the baseline),
zero-fill, FFT, reverse, ACME automatic phasing, discard the imaginary
channel, polynomial baseline, then calibrate the axis on a line of known
shift.

## The `valley_ratio` shim-quality gate

The main-chain CH₂ line is three or four orders of magnitude taller than every
signal you actually want. A Lorentzian has heavy tails and a badly shimmed
line has heavier ones, so the *foot* of the main-chain peak spills into the
neighbouring CH₃ window and gets integrated as chain ends. The spectrum looks
fine and the number is inflated.

```
valley_ratio = (lowest point between the CH₃ maximum and the main-chain maximum)
               / (main-chain maximum)
```

On a well-shimmed spectrum the two peaks come down to the noise between them
and the ratio is tiny. As the main-chain line broadens, the trough lifts off
the baseline and the ratio grows. Rejecting spectra above a threshold discards
exactly the ones whose CH₃ integral is contaminated, which is a failure no
later arithmetic can undo.

The threshold is empirical and belongs in the regions file. It depends on
field, temperature, solvent and how much CH₃ you expect. Set it by measuring
several spectra you trust and several you know are bad, then put the number
where it separates them.

The gate does not certify the CH₃ window as unbiased — see *Limitations*. It
separates a bias you can state from one that swamps the measurement.

## Solvent-line subtraction

Useful end-group windows sit next to the residual solvent line, and that line
is enormous. Its tail leaks into the window, and because the leak is partly
dispersive (a phasing artefact, negative on one side of the line), naive
integration of such a window can come out **negative**. The demo shows this:
the raw –CH= window integrates to about −26 per 1000 main-chain protons, and
to +1.49 after subtraction, against a true value of 1.50.

`solvent_fit` fits a phase-mixed Lorentzian — `cos(φ)·absorption +
sin(φ)·dispersion` plus a local linear baseline — and subtracts the line
shape. The fit windows are declared in the regions file and must **exclude**
the window being protected:

```json
{
  "applies_to": "vinyl_CH",
  "center_bounds": [5.85, 6.00],
  "fit_windows": [[5.90, 6.40], [5.55, 5.68]],
  "max_hwhm_ppm": 0.2
}
```

Keeping the fit off the peak of interest is what makes this honest. Fit across
the window and the model will happily absorb the signal and report that there
was nothing there. Only the line shape is subtracted, not the fitted local
baseline: that term soaks up the offset inside the fit windows and would tilt
the rest of the spectrum if it were extrapolated.

## Configuring the windows

Copy `regions.example.json` and edit it. Every shift in this package comes
from that file; there is no hard-coded ppm anywhere in `endgroup_nmr/`.

Each region declares its window, the chemical `group` it reports on, how many
protons of that group fall inside the window, optionally how many carbons the
group contributes to the chain (used for Mₙ), and a `role`:

* `backbone` — the intensity reference,
* `chain_end` — counted in `ends_per_chain`,
* `in_chain` — counted in Mₙ but not as an end.

A group may be declared through more than one window. The vinyl end appears as
both `vinyl_CH` (1H) and `vinyl_CH2` (2H); both should report the same molar
amount, and the code reduces them with the **median**, so one window spoiled by
an overlap does not drag the answer.

Calibrate the axis by setting `reference.ppm` to the shift of your residual
solvent line — 5.91 for 1,1,2,2-tetrachloroethane-d₂ (the usual choice for
polyolefins, which need heat to dissolve), 7.26 for CDCl₃, 5.32 for CD₂Cl₂,
7.16 for C₆D₆. Set it to `null` to keep the spectrometer's own referencing.

## Running it on real folders

**Bruker.** A dataset is a directory holding `acqus` plus `fid` (1D) or `ser`
(nD); a day's folder usually contains numbered sub-experiments `1/`, `2/`,
`3/`. Point the CLI at any level and it walks down:

```
python -m endgroup_nmr.cli --regions my_regions.json /data/2026-09-12 --csv results.csv
```

Each row is keyed by the dataset's path relative to the directory you gave.
`--drop-failed-qc` omits spectra that fail the `valley_ratio` gate instead of
flagging them; by default they are kept with `qc_passed = NO`, because seeing
the number you are rejecting is more useful than a blank.

**JEOL datasets.** Export to JCAMP-DX from Delta (*File → Export*) and pass
the `.jdx` file:

```
python -m endgroup_nmr.cli --regions my_regions.json spectrum.jdx
```

The binary `.jdf` container is deliberately not parsed here. The free parsers
for it are GPL-licensed, and vendoring one into an MIT project would relicense
this package; JCAMP-DX is a documented interchange format that nmrglue reads
directly. If the export is already Fourier-transformed — most are — the reader
detects that from the file's `DATATYPE` record and the pipeline skips straight
to baseline correction and calibration rather than transforming twice.

## Library use

```python
from endgroup_nmr import (load_regions, read_any, process, check_valley_ratio,
                          subtract_solvent_lines, quantify)

config = load_regions("regions.example.json")
spectrum = process(read_any(path), reference_ppm=config.reference_ppm)

qc = check_valley_ratio(spectrum, config)
if qc.passed is False:
    raise SystemExit(f"shim quality too poor: valley_ratio={qc.ratio:.4f}")

corrected, fits = subtract_solvent_lines(spectrum, config)
result = quantify(corrected, config)
print(result.ends_per_chain, result.mn, result.notes)
```

`process` returns a `Spectrum` (a ppm array, an intensity array and metadata),
so you can also build one from arrays you already have and use the
quantification half on its own — that is what the demo and the tests do.

## The demo

`examples/synthetic_demo.py` builds a ¹H spectrum from Lorentzians at the
textbook shifts, with a main-chain line thousands of times taller than the end
groups sitting on a broad skirt, a phase-mixed residual solvent line
overlapping the –CH= window, and noise. It then runs the whole pipeline on two
copies — one well shimmed, one with a deliberately broadened main-chain line —
and prints the QC verdicts, the solvent fit, the quantities, and a comparison
against the amounts the spectrum was built from.

It needs no instrument data, so it doubles as a smoke test:

```
python examples/synthetic_demo.py
python -m pytest tests/
```

## Limitations

* **The CH₃ window is biased upwards, always.** Even on a well-shimmed
  spectrum, a main-chain line thousands of times taller has measurable
  Lorentzian tail under 0.75–1.00 ppm, and it is integrated as if it were
  chain ends. In the demo the well-shimmed case reads about 28 % high on CH₃
  and 16 % high on ends/chain, while Mₙ and the chain length — which rest on
  the main-chain and unsaturation windows — stay within a couple of per cent.
  `valley_ratio` measures the size of that leak; it does not remove it. Modelling
  the main-chain tail and subtracting it, the way the solvent line is handled
  here, would be the way to do that, and this package does not do it.
* **CH₃ is not only chain ends.** Short-chain branches land in the same
  window. On a branched sample the CH₃ count is ends plus branches, and
  nothing in a 1D ¹H spectrum separates them.
* **Mₙ is a first-order estimate** — carbons per chain times the repeat-unit
  mass, from a chain count that rests on the assumption you configured. It is
  an end-group Mₙ: it says nothing about the distribution, and it degrades as
  chains get longer and end groups scarcer.
* **There is no detection floor.** `ends_per_chain` and `Mₙ` come back as
  `NaN` only when the chain-counting windows integrate to zero or less. A
  spectrum whose end-group windows hold nothing but main-chain tail still has
  a tiny positive chain count, so it returns a large finite `ends_per_chain`
  and a small `Mₙ` rather than refusing to answer. Judge that case on
  `valley_ratio` and on `result.chains`, which is reported for exactly this
  reason; a threshold would have to be calibrated per instrument and this
  package does not invent one for you.
* **Windows are windows, not assignments.** Anything nonzero deserves a look
  at the spectrum before it is believed. Overlaps between end-group windows
  and other species are common and the code cannot see them.
* **Automatic phasing can fail** on spectra with a huge dynamic range. If the
  numbers look strange, plot the trace before blaming the arithmetic.
* **2D datasets are skipped.** Bruker `ser` files are found but not handed to
  the 1D pipeline.

## Requirements

Python 3.10+, with `numpy`, `scipy` and `nmrglue`; `pytest` to run the tests.

## License

MIT — see [LICENSE](LICENSE).
