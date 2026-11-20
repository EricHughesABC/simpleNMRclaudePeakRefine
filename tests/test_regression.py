"""Regression tests: re-run every example and compare with its expected output.

Examples are read from examplesPublic/ (tracked) and examples/ (local only, git-ignored,
proprietary data).  Each folder holds inputs/<name>.json and expected/<name>_refined.json.
Folders or expected files that are missing are skipped, so the tests also run on a fresh clone.
"""
import json
from pathlib import Path

import pytest

from simplenmr_claude_peak_refine import LabelError, consistency, refine

ROOT = Path(__file__).resolve().parents[1]
CASES = sorted(p for d in ("examplesPublic", "examples") for p in (ROOT / d / "inputs").glob("*.json"))


def _strip(doc):
    doc = dict(doc)
    doc.pop("peakRefinement", None)
    return doc


@pytest.mark.parametrize("inp", CASES, ids=[f"{p.parents[1].name}/{p.stem}" for p in CASES])
def test_example(inp, tmp_path):
    expected = inp.parents[1] / "expected" / f"{inp.stem}_refined.json"
    out = tmp_path / f"{inp.stem}_refined.json"
    if not expected.exists():
        mislabel = inp.parents[1] / "expected" / f"{inp.stem}.labelerror"
        if mislabel.exists():                # this example must be rejected
            with pytest.raises(LabelError):
                refine(str(inp), str(out))
            return
        pytest.skip("no expected output")
    refine(str(inp), str(out))
    got, want = json.load(open(out)), json.load(open(expected))
    assert _strip(got) == _strip(want), "refined spectra differ from the expected output"
    _, _, _, off = consistency(str(out))
    assert not off, f"compound 2D peaks not on the proton/carbon lists: {off[:5]}"


def test_found_examples():
    if not CASES:
        pytest.skip("no example inputs in examplesPublic/ or examples/")
