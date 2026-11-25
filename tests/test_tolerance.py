"""The shared matching tolerance (tolerance.py) and its use in the refinement.

tolerance.py is a copy of simpleNMRtools core/tolerance.py and must stay identical to it.
If simpleNMRtools is checked out next to this repo (../simpleNMRtools), or its folder is
given in the environment variable SIMPLENMRTOOLS_DIR, the two files are compared byte for
byte; otherwise that test is skipped.
"""
import os
from pathlib import Path

import pytest

from simplenmr_claude_peak_refine import tolerance
from simplenmr_claude_peak_refine.refine import Spectrum, proton_tol_ppm

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / "src" / "simplenmr_claude_peak_refine" / "tolerance.py"


def _simplenmr_copy():
    cands = [os.environ.get("SIMPLENMRTOOLS_DIR"), ROOT.parent / "simpleNMRtools", ROOT.parent / "simplenmrtools"]
    for c in cands:
        if c and (Path(c) / "core" / "tolerance.py").exists():
            return Path(c) / "core" / "tolerance.py"
    return None


def test_same_as_simplenmr_copy():
    other = _simplenmr_copy()
    if other is None:
        pytest.skip("simpleNMRtools not found (set SIMPLENMRTOOLS_DIR to compare tolerance.py)")
    assert HERE.read_bytes().replace(b"\r\n", b"\n") == other.read_bytes().replace(b"\r\n", b"\n"), \
        f"tolerance.py differs from {other} - copy the simpleNMR file across unchanged"


def _spec(sw, si, sf, nucleus="[13C, 1H]"):
    raw = {"nucleus": nucleus, "sw_Hz": f"[{sw[0]}, {sw[1]}]", "spectral_size": f"[{si[0]}, {si[1]}]",
           "specfrequency": f"[{sf[0]}, {sf[1]}]"}
    return Spectrum("X", 2, ["13C", "1H"], list(sf), list(sw), [256, 1024], list(si), [], raw)


def test_hsqc_proton_window():
    # sugar HSQC, 700 MHz: f2 point 11363.6/2048 = 5.55 Hz -> N_POINTS points
    hs = _spec((35193.1, 11363.64), (4096, 2048), (175.97, 699.743))
    want = max(tolerance.N_POINTS * 11363.64 / 2048, tolerance.FLOOR_HZ) / 699.743
    assert proton_tol_ppm(hs) == pytest.approx(want)


def test_other_spectrum_uses_coarser_point():
    # a finer DDEPT f2 grid is matched against the coarser HSQC grid
    hs = _spec((35193.1, 11363.64), (4096, 2048), (175.97, 699.743))
    dd = _spec((15836.0, 4854.37), (1024, 1024), (175.96, 699.742))
    want = max(tolerance.N_POINTS * max(4854.37 / 1024, 11363.64 / 2048), tolerance.FLOOR_HZ) / 699.742
    assert proton_tol_ppm(dd, hs) == pytest.approx(want)


def test_floor():
    hs = _spec((35000.0, 2000.0), (4096, 4096), (150.9, 600.0))      # 0.49 Hz points
    assert proton_tol_ppm(hs) * 600.0 == pytest.approx(tolerance.FLOOR_HZ)
