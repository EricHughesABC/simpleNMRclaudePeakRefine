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

1. In Mnova, choose the spectra in the simpleNMR dialog and export the peaks. This uses the
   export script in `simpleNMRtools/MNOVAscripts`.
2. Run `simplenmr-claude-peak-refine` on the exported JSON.
3. Read the refined JSON back into Mnova with the read-back script in
   `simpleNMRtools/MNOVAscripts` (clear-and-replace of the 2D peaks, the 1H multiplets and
   integrals, and the 13C peaks).

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
