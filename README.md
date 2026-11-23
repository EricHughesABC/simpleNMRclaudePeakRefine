# simpleNMRclaudePeakRefine

Automatic refinement and harmonisation of NMR peak picking for
[simpleNMR](https://github.com/EricHughesABC/simpleNMRtools).

The project is directed by Eric Hughes (Durham University). The code is written by
Claude (Anthropic) under his direction.

Mnova exports the 1D and 2D peak lists (1H, pureshift 1H, 13C, HSQC, HMBC, COSY,
HSQC-CLIP-COSY, DEPT-HSQC) as a simpleNMR JSON file. This program builds **one proton
shift list and one carbon shift list** from all the evidence. It then rewrites every
spectrum so that its peaks sit exactly on those lists, and the refined file is read back
into Mnova.

The program never moves a peak further than the spectrum's resolution allows. A pick that
doesn't fit is typed as an artifact, not pushed onto the nearest resonance. The full logic
is in [`docs/refinement_logic.md`](docs/refinement_logic.md), also available as a
[PDF](docs/refinement_logic.pdf).

## Installation

```bash
pip install "simpleNMRclaudePeakRefine[formula] @ git+https://github.com/EricHughesABC/simpleNMRclaudePeakRefine.git"
```

The `[formula]` extra installs RDKit, which is used for the molecular-formula checks. The
program works without it.

For development:

```bash
git clone https://github.com/EricHughesABC/simpleNMRclaudePeakRefine.git
cd simpleNMRclaudePeakRefine
pip install -e ".[formula,test]"
```

## Use

```bash
simplenmr-claude-peak-refine export.json                       # writes export_refined.json
simplenmr-claude-peak-refine export.json refined.json --report report.txt
python -m simplenmr_claude_peak_refine export.json             # same thing
simplenmr-claude-peak-check refined.json                       # consistency check
```

The exit code is 0 on success. It is 2 when the spectrum labels chosen in Mnova contradict
the data; in that case no output is written and the error report says which spectrum to fix.

From Python:

```python
from simplenmr_claude_peak_refine import refine, Params
doc, report = refine("export.json", "refined.json")
```

## The Mnova round trip

1. In Mnova, run `collectPeaks.qs`. Its dialog (`peakRefineDialog.qs`) asks only which
   spectra to use and what each one is (H1_1D, H1_pureshift, C13_1D, HSQC, HMBC, …, or SKIP).
   It writes `<document>_claudePeakPicked_mresnova.json` next to the Mnova document. No
   server is contacted.
2. Run `simplenmr-claude-peak-refine` on the exported JSON.
3. Read the refined JSON back into Mnova with `readRevisedPeaks.qs` (clear-and-replace of
   the peaks). Before changing anything, it saves the current peaks of every spectrum it
   replaces to `<document>_peaksBeforeRefine.json` next to the Mnova document.
4. To go back to the unrefined peaks (for testing, or to try a different choice of
   spectra), run `undoRefinedPeaks.qs`. It restores the saved peaks: positions,
   intensities and, where Mnova allows it, type, kind, flags and annotation.

The export dialog remembers its choices, by spectrum, in `peakRefine_exptIdentifiers.json`
in the document's folder. This file is separate from the simpleNMR verify step's
`exptIdentifiers.json`, because the two steps usually use different spectra. The first
time, the verify step's choices are used as a starting point.

The two Mnova scripts stay in simpleNMRtools, because Mnova loads its scripts and toolbar
(`simpleNMRtools.mui`) from that one folder.

## Examples and tests

| Folder | Tracked in git | Contents |
|---|---|---|
| `examplesPublic/` | yes | exports that may be published: `inputs/<name>.json` and `expected/<name>_refined.json` |
| `examples/` | **no** (in `.gitignore`) | local examples, including proprietary data, in the same layout |

`pytest` re-runs every example in both folders, compares the refined spectra with the
expected output, and runs the consistency check. Folders that are missing are skipped, so
the tests also pass on a fresh clone. An empty file `expected/<name>.labelerror` marks an
example that must be rejected as mislabelled.

Mnova documents (`*.mnova`) are never tracked.

## Licence

MIT. See [LICENSE](LICENSE).
