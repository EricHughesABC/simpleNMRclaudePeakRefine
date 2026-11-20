"""Consistency check of a refined JSON file.

Every compound peak in every 2D spectrum must lie exactly on the refined proton list
(f2, and f1 for 1H-1H spectra) and the refined carbon list (f1 for 1H-13C spectra).
"""
import json
import sys


def consistency(path_or_doc):
    """Return (n_protons, n_carbons, n_compound_2d_peaks, off_list) for a refined file."""
    d = json.load(open(path_or_doc)) if isinstance(path_or_doc, str) else path_or_doc
    H, C = [], []
    for k, v in d.items():
        if not isinstance(v, dict) or v.get("datatype") != "nmrspectrum":
            continue
        if k.startswith("H1_1D"):
            H = [m["delta1"] for m in v["multiplets"]["data"].values() if m.get("type", 0) == 0]
        if k.startswith("C13_1D"):
            C = [p["delta1"] for p in v["peaks"]["data"].values() if p.get("type", 0) == 0]
    if not C:
        C = [c["delta"] for c in d["peakRefinement"]["data"]["0"]["carbons"]]
    bad, n = [], 0
    for k, v in d.items():
        if not isinstance(v, dict) or v.get("type") != "2D" or k.upper().startswith("SKIP"):
            continue
        hh = "[1H, 1H]" in str(v.get("nucleus")) or "['1H', '1H']" in str(v.get("nucleus"))
        for p in v["peaks"]["data"].values():
            if p.get("type", 0) != 0:
                continue
            n += 1
            ok_h = any(abs(p["delta2"] - x) < 1e-9 for x in H)
            ok_1 = any(abs(p["delta1"] - x) < 1e-9 for x in (H if hh else C))
            if not (ok_h and ok_1):
                bad.append((k, round(p["delta2"], 4), round(p["delta1"], 3), p.get("annotation", "")[:40]))
    return len(H), len(C), n, bad


def main():
    status = 0
    for f in sys.argv[1:]:
        h, c, n, bad = consistency(f)
        print(f"{f}: H {h}  C {c}  2D compound peaks {n}  off-list {len(bad)} {bad[:6]}")
        status |= bool(bad)
    sys.exit(status)


if __name__ == "__main__":
    main()
