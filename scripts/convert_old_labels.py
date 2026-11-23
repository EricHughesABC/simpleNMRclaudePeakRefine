"""Convert old exports to the current experiment names.

Older Mnova exports used names that simpleNMR does not recognise:
    PureShift       -> H1_pureshift
    DDEPT_CH3_ONLY  -> DDEPTCH3ONLY
This renames the spectrum keys (PureShift_0 -> H1_pureshift_0), the labels in
exptIdentifiers and the label at the end of each chosenSpectra entry.

    python scripts/convert_old_labels.py old.json [more.json ...]      # converts in place
    python scripts/convert_old_labels.py old.json -o new.json
"""
import argparse
import json
import re

RENAME = {"PureShift": "H1_pureshift", "DDEPT_CH3_ONLY": "DDEPTCH3ONLY"}


def convert(doc: dict) -> tuple[dict, int]:
    changes = 0
    out = {}
    for key, value in doc.items():
        m = re.fullmatch(r"(.+)_(\d+)", key)
        if m and m.group(1) in RENAME:
            key = f"{RENAME[m.group(1)]}_{m.group(2)}"
            changes += 1
        out[key] = value
    ids = out.get("exptIdentifiers", {}).get("data", {})
    for k, v in ids.items():
        if v in RENAME:
            ids[k] = RENAME[v]
            changes += 1
    chosen = out.get("chosenSpectra", {}).get("data", {})
    for k, v in chosen.items():
        head, _, label = str(v).rpartition(" ")
        if label in RENAME:
            chosen[k] = f"{head} {RENAME[label]}"
            changes += 1
    return out, changes


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("-o", "--output", help="output file (only with a single input; default: in place)")
    a = ap.parse_args()
    if a.output and len(a.files) != 1:
        ap.error("-o needs exactly one input file")
    for f in a.files:
        with open(f) as fh:
            doc = json.load(fh)
        new, n = convert(doc)
        target = a.output or f
        if n or a.output:
            with open(target, "w") as fh:
                json.dump(new, fh, indent=4)
        print(f"{f}: {n} name(s) changed" + (f" -> {target}" if a.output else ""))


if __name__ == "__main__":
    main()
