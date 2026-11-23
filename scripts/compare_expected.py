"""Re-run every example and show how the output differs from expected/.

    python scripts/compare_expected.py              # examplesPublic/ and examples/
    python scripts/compare_expected.py --update     # also overwrite expected/ with the new output

For each example it prints which spectra changed. For 2D spectra it lists the compound
peaks (f2, f1) removed and added, with their annotations, so a code change can be checked
by eye before the expected files are updated.
"""
import argparse
import json
import shutil
import tempfile
from pathlib import Path

from simplenmr_claude_peak_refine import LabelError, refine

ROOT = Path(__file__).resolve().parents[1]


def compound_peaks(spec):
    out = {}
    for p in spec.get("peaks", {}).get("data", {}).values():
        if p.get("type", 0) == 0:
            out[(round(p["delta2"], 4), round(p["delta1"], 4))] = p.get("annotation", "")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--update", action="store_true", help="overwrite expected/ with the new output")
    a = ap.parse_args()
    cases = sorted(p for d in ("examplesPublic", "examples") for p in (ROOT / d / "inputs").glob("*.json"))
    if not cases:
        print("no example inputs found in examplesPublic/inputs or examples/inputs")
        return
    tmp = Path(tempfile.mkdtemp())
    n_same = n_diff = 0
    for inp in cases:
        name = f"{inp.parents[1].name}/{inp.stem}"
        expected = inp.parents[1] / "expected" / f"{inp.stem}_refined.json"
        out = tmp / f"{inp.stem}_refined.json"
        try:
            refine(str(inp), str(out))
        except LabelError as e:
            print(f"{name}: label error ({e})")
            continue
        if not expected.exists():
            print(f"{name}: no expected output")
            continue
        got, want = json.load(open(out)), json.load(open(expected))
        changed = [k for k in want if k != "peakRefinement" and isinstance(want[k], dict)
                   and want[k].get("datatype") == "nmrspectrum" and got.get(k) != want[k]]
        if not changed:
            n_same += 1
            print(f"{name}: identical")
            continue
        n_diff += 1
        print(f"{name}: changed {', '.join(changed)}")
        for k in changed:
            if want[k].get("type") == "2D":
                old, new = compound_peaks(want[k]), compound_peaks(got[k])
                for pk in sorted(set(old) - set(new), reverse=True):
                    print(f"    - {k} {pk[0]:.4f} / {pk[1]:.4f}  {old[pk][:90]}")
                for pk in sorted(set(new) - set(old), reverse=True):
                    print(f"    + {k} {pk[0]:.4f} / {pk[1]:.4f}  {new[pk][:90]}")
                same = set(old) & set(new)
                n_ann = sum(old[pk] != new[pk] for pk in same)
                if n_ann:
                    print(f"    ~ {k}: {n_ann} peak(s) at the same position with a new annotation")
                if set(old) == set(new) and not n_ann:
                    print(f"    ~ {k}: same compound peaks; other fields differ (intensity, width or non-compound peaks)")
        if a.update:
            shutil.copy(out, expected)
            print(f"    expected/{expected.name} updated")
    print(f"\n{n_same} identical, {n_diff} changed")


if __name__ == "__main__":
    main()
