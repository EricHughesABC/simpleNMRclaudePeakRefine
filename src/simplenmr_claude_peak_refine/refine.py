"""
simpleNMRclaudePeakRefine (refine.py)  -  v0.28

Refine Mnova HSQC / HMBC peak picking using 1D 1H, pureshift and 13C peak lists.

Input : simpleNMR-style JSON exported from Mnova (keys from 'exptIdentifiers',
        e.g. H1_1D_0, H1_pureshift_0, C13_1D_0, HSQC_0, HMBC_0, COSY_0).
Output: the same JSON with the peak lists of every correlation spectrum replaced
        (clear-and-replace write-back): HSQC, HMBC, COSY and any other 1H-13C 2D
        (HSQC-CLIP-COSY, DEPT-edited HSQC ...), all snapped to one proton and one carbon
        shift list.  A top-level 'peakRefinement' block holds the report, the proton/carbon
        shift lists and the original peak lists.
        The 1D spectra are made consistent with the same lists:
          1H  - multiplets and integrals rebuilt, one compound multiplet per proton at its
                refined shift (Mnova category / J values kept where the multiplet was not split,
                nH from Mnova, the integrals, or the HSQC CH2 rule); 1H line peaks unchanged;
          13C - one compound peak per carbon; compound-typed peaks that are not carbons are
                retyped (satellite/noise and line-fit splits -> artifact, minor -> impurity).

Conventions (as exported by the Mnova QS script)
  1D peaks : delta1 = shift (ppm)
  2D peaks : delta1 = f1 (indirect, e.g. 13C), delta2 = f2 (1H)
  width1/width2 in ppm; 2D metadata strings are "[f1, f2]"
  acquired_size = complex points
  type  : 0 compound, 1 artifact, 2 impurity, 3 solvent, 4 CS ref, 5 quant ref
  kind  : 1 maximum (positive), 2 minimum (negative)
  flags : integer; documented Peak.Flags live in the low bits and are extracted
          with FLAG_MASK.  Higher bits are undocumented Mnova internals and are
          preserved but ignored.

1H proton list: Mnova multiplets + integrals when present (pureshift peaks take priority
inside each multiplet); otherwise multiplets are built from the 1H line picks (symmetry
centres, lines hidden under solvent allowed) and split by HSQC carbons.

Spectrum labels (exptIdentifiers, chosen by the user in the Mnova dialog) are authoritative:
H1_1D, H1_pureshift, C13_1D, COSY, HSQC, HMBC, HSQC_CLIPCOSY, DDEPTCH3ONLY are used as labelled
(dataset <label>_0; a second dataset <label>_1 is left unchanged for now); NOESY and DEPT135
are checked but not refined yet; SKIP (and unknown labels, predicted spectra, spectra without
acquisition metadata) are not analysed.  Everything not refined is left unchanged in the output.  Each label is checked against the nucleus,
dimension, pulse programme and Mnova experiment type; a contradiction stops the program with
an error report (exit code 2, no output written) so the labels can be fixed in Mnova.
H1_1D and HSQC are required.

Molecular formula (v0.24): if the export carries the proposed structure ('smiles' or 'molfile',
count > 0), only its molecular formula is used (computed with RDKit) - never the connectivity
or symmetry, so the peak picking stays unbiased between isomers.  At present the formula is
used for checks in the report (carbon count, total H, CH3 counts, NH/OH); no peak decision
depends on it.  Without a structure, or without RDKit, the program runs blind as before.

Requirements: Python >= 3.10, numpy (scipy optional, rdkit optional for the molecular formula).

Usage
  simplenmr-claude-peak-refine input.json [output.json] [--report report.txt]
  python -m simplenmr_claude_peak_refine input.json [output.json] [--report report.txt]
  (output defaults to <input>_refined.json; the report is also printed to the screen)
"""
from __future__ import annotations

import argparse
import copy
import json
import sys
import re
import math
import statistics as st
from dataclasses import dataclass, field

import numpy as np

try:                                    # package import (normal use)
    from . import tolerance
except ImportError:                     # refine.py run on its own, tolerance.py beside it
    import tolerance

try:
    from scipy.optimize import linear_sum_assignment
except ImportError:  # pragma: no cover
    linear_sum_assignment = None

VERSION = "0.28"

TYPE_NAMES = {0: "compound", 1: "artifact", 2: "impurity", 3: "solvent",
              4: "cs_ref", 5: "quant_ref"}
COMPOUND, ARTIFACT, IMPURITY, SOLVENT = 0, 1, 2, 3
FLAG_BITS = {"hidden": 1, "weak": 2, "c13_satellite": 4, "rotational": 8,
             "labile": 16, "solvent_2ary": 256}
FLAG_MASK = sum(FLAG_BITS.values())


# --------------------------------------------------------------------------- #
# Parameters (all tolerances in Hz or digital points, converted per spectrum)
# --------------------------------------------------------------------------- #
@dataclass
class Params:
    ps_area_frac: float = 0.25      # pureshift peak kept as resonance if area >= frac * largest in multiplet range
    ps_range_pad_hz: float = 2.0
    ps_dispersive_frac: float = 0.25  # pureshift unreliable in a multiplet when |negative| >= this x max positive
    ps_merge_hz: float = 3.0        # pureshift peaks closer than this inside one multiplet are one resonance    # padding around Mnova 1H integral ranges
    h1_cluster_gap_hz: float = 16.0  # no Mnova multiplets: 1H lines closer than this form one multiplet
    h1_min_cluster_frac: float = 0.1 # ... a line cluster below this x the median cluster area needs 2D support
    split_f2_gap_hz: float = 16.0    # ... HSQC picks on one carbon group further apart than this in f2 are different protons
    h1_line_assign_hz: float = 25.0 # 1H lines assigned to nearest resonance within this
    hsqc_tol_pts: float = 2.0       # HSQC carbon (f1) match tolerance, in f1 digital points.  The proton (f2)
                                    # tolerance of the HSQC and other H-C 2D spectra comes from tolerance.py,
                                    # shared with simpleNMR (tolerance.N_POINTS, tolerance.FLOOR_HZ)
    split_tol_pts: float = 2.0      # HSQC proton (f2) resolution test - two signals are different when further apart
    split_min_tol_hz: float = 4.0   # than this many HSQC f2 points (at least this in Hz): HSQC row vs pureshift peak,
                                    # and pooling HSQC picks that have no proton.  A resolution test, not a matching
                                    # window, so it stays tighter than the shared matching tolerance
    hmbc_f1_tol_pts: float = 2.0    # HMBC f1 -> 13C tolerance in f1 points (picks sit on the grid)
    hmbc_f1_tol_res: float = 0.25   # ... or this fraction of the f1 resolution (1/AQ), whichever is larger
    hmbc_gap_hz: float = 12.0       # lines closer than this in a carbon row form one cluster
    hmbc_pair_max_hz: float = 25.0  # max separation of two clusters merged as a symmetric pair
    hmbc_sym_tol_hz: float = 3.0    # midpoint of a symmetric pair must be this close to a resonance
    hmbc_match_hz: float = 6.0      # cluster centroid -> resonance tolerance
    hmbc_presence_hz: float = 15.0  # single line within this of a resonance = presence-only evidence
    hmbc_single_hz: float = 3.0
    hmbc_single_res: float = 0.6
    hmbc_weak_frac: float = 0.2     # correlation weaker than this x the proton's strongest is marked weak (2J/4J?)    # ... or this x the HMBC f2 resolution, whichever is larger
    wiggle_ratio: float = 5.0       # a pick this many times weaker than a peak in the same f2 column ...
    wiggle_res: float = 2.5         # ... within this x the f1 resolution is a truncation wiggle     # a single line this close to a resonance is a (weak) position estimate
    onebond_min_hz: float = 115.0   # 1J(CH) range used to recognise breakthrough
    onebond_max_hz: float = 180.0
    sideband_tol_hz: float = 2.0
    cosy_f2_tol_hz: float = 6.0
    cosy_f1_tol_pts: float = 1.0
    c13_min_frac: float = 0.02       # 13C compound peaks below this fraction of the reference intensity are ignored
    c13_minor_frac: float = 0.1      # 13C peaks below this x reference, separated by a gap ...
    c13_minor_gap: float = 4.0       # ... of at least this ratio from the rest, need 2D support
    split_weak_frac: float = 0.1     # HSQC splitting: a weak group is still a proton if on its own distant 13C carbon
    c13_merge_hz: float = 3.0        # 13C compound peaks closer than this are one carbon (GSD-split lines)
    c13_ref_n: int = 5               # reference = median of the N strongest compound peaks (solvent excluded)
    c13_solvent_tol_ppm: float = 0.3
    solvent_h_clear_hz: float = 10.0  # a pick this far from the solvent 1H position can belong to the compound # compound-typed 13C peaks this close to a solvent position are not used
    hsqc_keep_duplicates: bool = False
    split_min_frac: float = 0.25     # HSQC-row splitting: a carbon group needs this share of the strongest in the multiplet
    ch2_merge_frac: float = 0.25     # no 13C: two CH2 HSQC groups closer than this x f1 resolution = one carbon
    hmbc_q_min_peaks: int = 2        # no 13C: an HMBC f1 row needs this many peaks to become a quaternary carbon  # False: duplicate HSQC picks are dropped from the output
    strong_coupling_hz: float = 20.0  # COSY-connected pair closer than this flagged
    cosy_group_share: float = 0.35       # unresolved group: 2nd member also kept if the f2 centroid lies between the two, at least this fraction of the way from the nearer one
    cosy_weak_frac: float = 0.2          # COSY pair weaker than this x the strongest pair of either proton is marked weak
    cosy_keep_diagonal: bool = False     # COSY output: one peak per coupled pair, f2 = higher shift
    clipcosy_keep_direct: bool = False   # HSQC-CLIP-COSY output: relay peaks only (direct ones duplicate the HSQC)
    hc_tol_pts: float = 2.0              # other H-C 2D spectra: strict f1 tolerance, digital points
    hc_expected_res: float = 2.0         # ... wider f1 tolerance (x f1 resolution) for expected direct/relay positions
    reclassify_ps: bool = True        # test Mnova solvent/impurity-typed pureshift peaks inside compound multiplets
    reclassify_hmbc: bool = True      # absorb impurity/artifact-typed HMBC lines that join compound clusters


# --------------------------------------------------------------------------- #
# Data model
# --------------------------------------------------------------------------- #
@dataclass
class Peak:
    key: str
    d1: float
    d2: float
    intensity: float
    width1: float
    width2: float
    type: int
    kind: int
    flags_raw: int
    annotation: str
    raw: dict

    @property
    def flags(self) -> int:
        return self.flags_raw & FLAG_MASK

    def has(self, name: str) -> bool:
        return bool(self.flags & FLAG_BITS[name])

    @property
    def area(self) -> float:
        return abs(self.intensity) * max(self.width1, 1e-9)


@dataclass
class Spectrum:
    name: str
    ndim: int
    nuclei: list
    sf: list          # MHz, [f1, f2] for 2D, [f] for 1D
    sw: list          # Hz
    td: list          # complex points acquired
    si: list          # points after zero filling
    peaks: list
    raw: dict
    key: str = ""     # key of this spectrum in the JSON document
    kind_fixed: int = 0   # picks whose 'kind' disagreed with the sign of the intensity
    nh_scale: int = 1     # 1H only: integrals found normalised to nh_scale x 1H
    excluded: dict = field(default_factory=dict)   # 13C only: peak key -> (type, reason)

    def res_hz(self, dim: int = -1) -> float:
        """Approximate resolution limit 1/AQ in Hz."""
        return self.sw[dim] / self.td[dim]

    def pt_hz(self, dim: int = -1) -> float:
        """Digital point spacing in Hz."""
        return self.sw[dim] / self.si[dim]

    def compound(self):
        return [p for p in self.peaks if p.type == COMPOUND]


def proton_tol_ppm(spec: "Spectrum", ref: "Spectrum | None" = None) -> float:
    """Tolerance in ppm for giving the proton (f2) of a 2D pick to a proton.

    The same rule simpleNMR uses when it matches HMBC/COSY protons to the HSQC rows
    (tolerance.py, copied unchanged from simpleNMRtools core/tolerance.py):
        max(N_POINTS x the coarser f2 point spacing of this spectrum and the HSQC, FLOOR_HZ)
    ref is the HSQC (None when spec is the HSQC itself)."""
    ref = ref if ref is not None else spec
    doc = {"SPEC_0": spec.raw, "HSQC_0": ref.raw}
    tol = tolerance.match_tolerance_ppm(doc, "SPEC", -1, ref_label="HSQC", ref_nucleus="1H")
    if tol is None:                      # metadata not readable by tolerance.py: same rule from the parsed values
        tol = max(tolerance.N_POINTS * max(spec.pt_hz(-1), ref.pt_hz(-1)), tolerance.FLOOR_HZ) / spec.sf[-1]
    return tol


def _as_list(v):
    if isinstance(v, str):
        v = v.strip().strip("[]")
        out = []
        for x in v.split(","):
            x = x.strip()
            try:
                out.append(float(x))
            except ValueError:
                out.append(x)
        return out
    return v if isinstance(v, list) else [v]


# 13C solvent positions (ppm), Fulmer et al. 2010; residual-solvent lines only
SOLVENT_C13 = {
    "CDCL3": [77.16], "DMSO": [39.52], "CD3OD": [49.00], "MEOD": [49.00],
    "ACETONE": [29.84, 206.26], "C6D6": [128.06], "CD3CN": [1.32, 118.26], "ACN": [1.32, 118.26],
    "CD2CL2": [53.84], "THF": [25.31, 67.21], "TOLUENE": [20.43, 125.13, 127.96, 128.87, 137.48],
    "D2O": [],
}


# residual-solvent 1H positions (ppm), Fulmer et al. 2010
SOLVENT_H1 = {
    "CDCL3": [7.26], "DMSO": [2.50], "CD3OD": [3.31], "MEOD": [3.31], "ACETONE": [2.05],
    "C6D6": [7.16], "CD3CN": [1.94], "ACN": [1.94], "CD2CL2": [5.32], "THF": [1.72, 3.58],
    "TOLUENE": [2.08, 6.97, 7.01, 7.09], "D2O": [4.79],
}


def solvent_h1_positions(solvent: str) -> list:
    key = str(solvent).upper().replace("-", "").replace("_", "")
    for k, v in SOLVENT_H1.items():
        if k in key:
            return v
    return []


def solvent_c13_positions(solvent: str) -> list:
    key = str(solvent).upper().replace("-", "").replace("_", "")
    for k, v in SOLVENT_C13.items():
        if k in key:
            return v
    return []


# Experiment names - the same as simpleNMR's NMREXPERIMENTS (simpleNMRtools config/globals.py),
# simpleNMRbuilder and SIMPLENMR_EXPERIMENTS in the Mnova scripts' simpleUtils.qs.  Spectra in
# the JSON are keyed <name>_<n> (HSQC_0, H1_pureshift_0, ...).  Use these names, never the strings.
H1_1D, C13_1D, DEPT135, H1_PURESHIFT = "H1_1D", "C13_1D", "DEPT135", "H1_pureshift"
HSQC, HMBC, COSY, NOESY = "HSQC", "HMBC", "COSY", "NOESY"
HSQC_CLIPCOSY, DDEPTCH3ONLY = "HSQC_CLIPCOSY", "DDEPTCH3ONLY"
SKIP = "SKIP"
NMREXPERIMENTS = (H1_1D, C13_1D, DEPT135, H1_PURESHIFT, HSQC, HMBC, COSY, NOESY, HSQC_CLIPCOSY, DDEPTCH3ONLY)

PURESHIFT_PULSE = ("psyche", "pure", "reset_", "zangger", "bird", "pshift")


# the spectra keys of the Mnova/simpleNMR spectrum-selection dialog
KNOWN_LABELS = (H1_1D, H1_PURESHIFT, C13_1D, COSY, HSQC, HMBC, HSQC_CLIPCOSY, DDEPTCH3ONLY,
                NOESY, DEPT135)
REFINED_LABELS = (H1_1D, H1_PURESHIFT, C13_1D, COSY, HSQC, HMBC, HSQC_CLIPCOSY, DDEPTCH3ONLY)


class LabelError(Exception):
    """The spectrum labels chosen in Mnova contradict the data; the user must fix them."""
    def __init__(self, errors, warnings=()):
        super().__init__("; ".join(errors))
        self.errors, self.warnings = list(errors), list(warnings)


def check_label(s: dict, label: str):
    """Check one spectrum against the role the user chose for it in Mnova.

    The label is authoritative - the program never changes it.  Returns (errors, warnings):
    an error is a clear contradiction (wrong nucleus or dimension, a pulse programme or
    Mnova experiment type of another kind); a warning is a doubt the data cannot settle."""
    err, warn = [], []
    nuc = str(s.get("nucleus", "")).replace(" ", "")
    pp = str(s.get("pulsesequence", "")).lower()
    exp = str(s.get("experimenttype", "")).upper()
    dim = s.get("type")
    is_ps = any(t in pp for t in PURESHIFT_PULSE)
    desc = f"{dim} {nuc}, pulse programme '{s.get('pulsesequence', '')}', Mnova type '{s.get('experimenttype', '')}'"

    def need(cond, what):
        if not cond:
            err.append(f"labelled {label} but the data is {desc} - {what}")
    hc = nuc in ("[13C,1H]", "[1H,13C]")
    if label == H1_1D:
        need(dim == "1D" and nuc == "1H", "expected a 1D 1H spectrum")
        if dim == "1D" and nuc == "1H" and is_ps:
            err.append(f"labelled H1_1D but the pulse programme '{s.get('pulsesequence')}' is a pureshift "
                       f"experiment - should it be H1_pureshift?")
    elif label == H1_PURESHIFT:
        need(dim == "1D" and nuc == "1H", "expected a 1D 1H pureshift spectrum")
        if dim == "1D" and nuc == "1H" and not is_ps:
            if re.match(r"zg\d*$", pp):
                err.append(f"labelled H1_pureshift but the pulse programme '{s.get('pulsesequence')}' is a standard "
                           f"1H experiment - should it be H1_1D?")
            else:
                warn.append(f"labelled H1_pureshift; pulse programme '{s.get('pulsesequence')}' not recognised as "
                            f"pureshift - used as labelled")
    elif label == C13_1D:
        need(dim == "1D" and nuc == "13C", "expected a 1D 13C spectrum")
    elif label == DEPT135:
        need(dim == "1D" and nuc == "13C", "expected a 1D 13C DEPT-135 spectrum")
        if dim == "1D" and nuc == "13C" and "dept" not in pp:
            err.append(f"labelled DEPT135 but the pulse programme '{s.get('pulsesequence')}' is not a DEPT "
                       f"experiment - should it be C13_1D?")
    elif label == NOESY:
        need(dim == "2D" and nuc == "[1H,1H]", "expected a 2D 1H-1H NOESY/ROESY")
        if dim == "2D" and nuc == "[1H,1H]" and (COSY in exp or "cosy" in pp or "TOCSY" in exp) \
                and not any(t in exp or t.lower() in pp for t in (NOESY, "ROESY")):
            err.append(f"labelled NOESY but the data is a COSY/TOCSY ({desc})")
    elif label == COSY:
        need(dim == "2D" and nuc == "[1H,1H]", "expected a 2D 1H-1H COSY")
        other = next((t for t in (NOESY, "ROESY", "TOCSY") if t in exp or t.lower() in pp), None)
        if other:
            err.append(f"labelled COSY but the data is a {other} ({desc})")
    elif label in (HSQC, HMBC, HSQC_CLIPCOSY, DDEPTCH3ONLY):
        need(dim == "2D" and hc, "expected a 2D 13C-1H spectrum")
        kind = (HMBC if HMBC in exp or "hmbc" in pp else
                "CLIP-COSY" if "clip" in pp else
                "DEPT-HSQC" if "dept" in pp else
                HSQC if HSQC in exp or "hsqc" in pp else None)
        expect = {HSQC: (HSQC,), HMBC: (HMBC,), HSQC_CLIPCOSY: ("CLIP-COSY",),
                  DDEPTCH3ONLY: ("DEPT-HSQC",)}[label]
        if kind and kind not in expect and not (label == HSQC_CLIPCOSY and kind == HSQC and "clip" in pp):
            err.append(f"labelled {label} but the data looks like a {kind} ({desc})")
        if kind is None and dim == "2D" and hc:
            warn.append(f"note: pulse programme '{s.get('pulsesequence')}' and Mnova type do not identify the "
                        f"experiment (e.g. NOAH) - used as labelled {label}")
    return err, warn


def load(path: str):
    with open(path) as fh:
        doc = json.load(fh)
    spectra, seen = {}, {}
    load.skipped, load.errors, load.warnings = [], [], []
    for ident in doc["exptIdentifiers"]["data"].values():
        n = seen.get(ident, 0)
        seen[ident] = n + 1
        key = f"{ident}_{n}"
        s = doc.get(key)
        if not s:
            continue
        why = ("labelled SKIP" if str(ident).strip().upper().startswith(SKIP) else
               "predicted spectrum" if str(s.get("origin", "")).lower() == "predictor"
               or "predicted" in str(s.get("subtype", "")).lower() else
               "no acquisition metadata" if not str(s.get("acquired_size", "")).strip("[] ,") else None)
        if why:               # ignored for the analysis and left unchanged in the output
            load.skipped.append(f"{key} ({s.get('subtype', '') or s.get('experiment', '')}, {why})")
            continue
        if ident not in KNOWN_LABELS:
            load.warnings.append(f"{key}: label '{ident}' is not one of the experiment names "
                                 f"({', '.join(NMREXPERIMENTS)}) - left unchanged; an export made before "
                                 f"v0.26 can be updated with scripts/convert_old_labels.py")
            load.skipped.append(f"{key} ({s.get('subtype', '') or s.get('experiment', '')}, unknown label)")
            continue
        e, w = check_label(s, ident)          # the label is authoritative; this only catches mislabels
        load.errors += [f"{key}: {x}" for x in e]
        load.warnings += [f"{key}: {x}" for x in w]
        if ident not in REFINED_LABELS:
            load.skipped.append(f"{key} ({ident}, not refined yet)")
            continue
        name = ident
        if name in spectra:                   # e.g. HSQC_1: a second, higher-resolution or narrowed dataset
            load.warnings.append(f"{key}: only {spectra[name].key} is refined at present (combining several "
                                 f"{ident} datasets is not supported yet) - {key} left unchanged")
            load.skipped.append(f"{key} (second {ident} dataset)")
            continue
        peaks, nfix = [], 0
        for k, p in s["peaks"]["data"].items():
            kind = p.get("kind", 1)
            if p["intensity"] < 0 and kind != 2 or p["intensity"] > 0 and kind == 2:
                kind = 2 if p["intensity"] < 0 else 1      # kind disagrees with the sign: trust the sign
                nfix += 1
            peaks.append(Peak(k, p["delta1"], p.get("delta2", 0.0), p["intensity"],
                              p.get("width1", 0.0), p.get("width2", 0.0),
                              p.get("type", 0), kind, p.get("flags", 0),
                              p.get("annotation", ""), p))
        spectra[name] = Spectrum(name, 2 if s["type"] == "2D" else 1,
                                  _as_list(s["nucleus"]),
                                  [float(x) for x in _as_list(s["specfrequency"])],
                                  [float(x) for x in _as_list(s["sw_Hz"])],
                                  [float(x) for x in _as_list(s["acquired_size"])],
                                  [float(x) for x in _as_list(s["spectral_size"])],
                                  peaks, s, key)
        spectra[name].kind_fixed = nfix
    for req in (H1_1D, HSQC):
        if req not in spectra:
            load.errors.append(f"no spectrum labelled {req} - it is required")
    if load.errors:
        raise LabelError(load.errors, load.warnings)
    return doc, spectra


# --------------------------------------------------------------------------- #
# Systematic line splitting (distorted lineshape, e.g. poor shimming)
# --------------------------------------------------------------------------- #
def split_tolerance(pos_hz, ints, base_hz, cap_hz=6.0, min_peaks=5, min_frac=0.6):
    """Merge tolerance for one spectrum.  Normally base_hz.  If most strong peaks have a partner
    of comparable height within cap_hz (every resonance picked as 2-3 lines by the line fit -
    a distorted lineshape, typically poor shimming), the tolerance is raised to cover the
    observed splitting.  Returns (tolerance_hz, fraction_split, n_peaks)."""
    pk = sorted(zip(pos_hz, ints))
    n = len(pk)
    if n < min_peaks:
        return base_hz, 0.0, n
    gaps = []
    for i, (x, a) in enumerate(pk):
        nb = [abs(y - x) for j, (y, b) in enumerate(pk) if j != i and b >= 0.25 * a]
        g = min(nb) if nb else None
        if g is not None and g <= cap_hz:
            gaps.append(g)
    frac = len(gaps) / n
    if frac < min_frac:
        return base_hz, frac, n
    gaps.sort()
    q = gaps[min(len(gaps) - 1, int(0.9 * len(gaps)))]
    return max(base_hz, min(cap_hz, q + 0.5)), frac, n


def merge_clusters(peaks, pos_hz, tol_hz, max_span_hz):
    """Single-linkage clusters of peaks (gap to the neighbour <= tol_hz, total span <= max_span_hz).
    Returns a list of clusters, each a list of peaks, strongest first."""
    order = sorted(range(len(peaks)), key=lambda i: pos_hz[i])
    clusters, cur = [], []
    for i in order:
        if cur and (pos_hz[i] - pos_hz[cur[-1]] > tol_hz or pos_hz[i] - pos_hz[cur[0]] > max_span_hz):
            clusters.append(cur)
            cur = []
        cur.append(i)
    if cur:
        clusters.append(cur)
    return [sorted((peaks[i] for i in c), key=lambda p: -abs(p.intensity)) for c in clusters]


def _centroid_copy(cl):
    """The strongest peak of a split cluster, moved to the intensity-weighted centre."""
    rep = copy.copy(cl[0])
    w = [abs(p.intensity) for p in cl]
    rep.d1 = sum(p.d1 * x for p, x in zip(cl, w)) / sum(w)
    return rep


# --------------------------------------------------------------------------- #
# Molecular formula from the proposed structure (formula only - no connectivity)
# --------------------------------------------------------------------------- #
def molecular_formula(doc):
    """Return ({'C': n, 'H': n, ...}, formula string, source) or (None, reason, None).
    Uses RDKit on the SMILES (preferred) or the molfile; implicit hydrogens are counted."""
    def first(key):
        b = doc.get(key)
        if isinstance(b, dict) and b.get("count") and isinstance(b.get("data"), dict):
            v = next(iter(b["data"].values()), "")
            return v if isinstance(v, str) and v.strip() else None
        return None
    smi, mol_block = first("smiles"), first("molfile")
    if not smi and not mol_block:
        return None, "no structure in the export - formula checks not run", None
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem.rdMolDescriptors import CalcMolFormula
        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        return None, "structure present but RDKit is not installed (pip install rdkit) - formula checks not run", None
    for src, make in (("SMILES", lambda: Chem.MolFromSmiles(smi) if smi else None),
                      ("molfile", lambda: Chem.MolFromMolBlock(mol_block.replace("\r\n", "\n")) if mol_block else None)):
        try:
            m = make()
        except Exception:      # pragma: no cover - malformed input
            m = None
        if m is None:
            continue
        counts = {}
        for a in Chem.AddHs(m).GetAtoms():
            counts[a.GetSymbol()] = counts.get(a.GetSymbol(), 0) + 1
        if "D" in counts:
            counts["H"] = counts.get("H", 0) + counts.pop("D")
        return counts, CalcMolFormula(m), src
    return None, "structure present but RDKit could not read the SMILES or the molfile - formula checks not run", None


def formula_checks(formula, protons, carbons, new, sp, extra, rep):
    """Report-only checks of the refined lists against the molecular formula."""
    counts, text, src = formula
    if counts is None:
        rep.append(f"Molecular formula: {text}")
        return {}
    nC, nH = counts.get("C", 0), counts.get("H", 0)
    nX = sum(counts.get(x, 0) for x in ("F", "Cl", "Br", "I"))
    nN = counts.get("N", 0) + counts.get("P", 0)
    dbe = nC - (nH + nX) / 2 + nN / 2 + 1
    het = {x: counts[x] for x in ("N", "O", "S", "P") if counts.get(x)}
    rep.append(f"Molecular formula: {text} (from the {src}; formula only - structure not used), "
               f"degrees of unsaturation {dbe:g}")
    out = {"formula": text, "source": src, "counts": counts, "dbe": dbe}
    # carbons: the formula count is an upper limit (symmetry can make carbons equivalent)
    nq = sum(1 for c in carbons if not c.protons)
    if len(carbons) > nC:
        rep.append(f"  Formula check - carbons: {len(carbons)} found but the formula has {nC} - at least "
                   f"{len(carbons) - nC} must be impurity/artefact carbons (or the proposed structure is wrong)")
    elif len(carbons) < nC:
        rep.append(f"  Formula check - carbons: {len(carbons)} found ({len(carbons) - nq} protonated, {nq} without "
                   f"HSQC) of {nC} in the formula - {nC - len(carbons)} not found: symmetry-equivalent, overlapped "
                   f"or quaternaries without enough HMBC evidence")
    else:
        rep.append(f"  Formula check - carbons: {len(carbons)} found, as in the formula")
    # hydrogens
    h1 = sp.get(H1_1D)
    if h1 is not None:
        mults = [m for m in new[h1.key]["multiplets"]["data"].values() if m.get("type", 0) == COMPOUND]
        tot = sum(m["nH"] for m in mults)
        by_delta = {round(m["delta1"], 9): m["nH"] for m in mults}
        nh_of = {h.id: by_delta.get(round(h.delta, 9)) for h in protons}
        on_c = sum(n or 0 for h in protons if h.carbons for n in [nh_of[h.id]])
        off_c = [h for h in protons if not h.carbons]
        msg = f"  Formula check - hydrogens: {tot} H in the refined 1H multiplets, formula {nH}"
        msg += " - agrees" if tot == nH else f" - {'short' if tot < nH else 'over'} by {abs(nH - tot)}"
        rep.append(msg)
        rep.append(f"  Formula check - {on_c} H on carbon-bearing protons; {len(off_c)} proton(s) without a carbon "
                   f"({', '.join(f'{h.delta:.4f}' for h in off_c) or 'none'}); heteroatoms "
                   f"{', '.join(f'{k}{v}' for k, v in het.items()) or 'none'}"
                   + (" - no N/O/S/P in the formula, so protons without a carbon cannot be NH/OH/SH"
                      if off_c and not het else ""))
        # CH3 carbons seen in the DEPT-HSQC (CH3 only) should carry 3H
        dd = extra.get(DDEPTCH3ONLY) if extra else None
        if dd:
            ch3 = {p["delta1"] for p in dd if p.get("type") == COMPOUND}
            cmap = {c.id: c for c in carbons}
            bad = [(h, nh_of[h.id]) for h in protons for c in h.carbons
                   if c in cmap and any(abs(cmap[c].delta - x) < 1e-6 for x in ch3) and nh_of[h.id] not in (None, 3)]
            for h, n in bad:
                rep.append(f"  Formula check - {h.id} {h.delta:.4f} is on a CH3 carbon (DDEPTCH3ONLY) but its 1H "
                           f"multiplet has nH {n} - expected 3; check the 1H integral")
            out["ch3_carbons"] = sorted(ch3)
        out.update({"H_written": tot, "H_on_carbon": on_c, "protons_without_carbon": [h.delta for h in off_c]})
    out.update({"carbons_found": len(carbons)})
    return out


# --------------------------------------------------------------------------- #
# Proton resonances (pureshift + 1H)
# --------------------------------------------------------------------------- #
@dataclass
class Proton:
    id: str
    delta: float
    source: str                     # 'pureshift' | 'H1-multiplet' | 'consensus'
    mult_delta: float | None = None # Mnova multiplet delta
    mult_category: str = ""
    range: tuple = (0.0, 0.0)
    h1_centroid: float | None = None
    nH_est: float | None = None
    evidence: dict = field(default_factory=dict)   # source -> list of (delta, weight)
    flags: list = field(default_factory=list)
    carbons: list = field(default_factory=list)    # HSQC partners
    spin_system: int | None = None
    mnova_type: int = COMPOUND      # Mnova type of the seeding pureshift peak
    split_c: float | None = None    # f1 of the HSQC carbon group this proton was split out for
    strong: bool = False            # in a strongly coupled region (pureshift unreliable)
    mult_nh: float | None = None    # Mnova nH of the (merged) multiplet the proton came from


def split_by_hsqc(sp, protons, P: Params, rep: list, strong_pos=None):
    """Use HSQC rows to find protons hidden inside one 1H multiplet.

    Each HSQC carbon row is a carbon-filtered 1H spectrum, so a 1H multiplet whose range
    contains cross-peaks at two or more different carbons holds at least that many protons.
      - no pureshift peak in the multiplet: replace the multiplet proton by one proton per
        carbon, positioned at the intensity-weighted f2 centroid of that carbon's HSQC peaks;
      - pureshift present but a carbon group has no pureshift proton within tolerance
        (e.g. a peak lost to strong coupling): add a proton for that group."""
    hs = sp.get(HSQC)
    if hs is None or not hs.compound():
        return protons
    sfC, sfH = hs.sf
    tolH = max(P.split_tol_pts * hs.pt_hz(1), P.split_min_tol_hz) / sfH   # resolution test (see Params)
    tolC = P.hsqc_tol_pts * hs.pt_hz(0) / sfC
    pad = P.ps_range_pad_hz / sfH
    by_range = {}
    for pr in protons:
        by_range.setdefault(tuple(pr.range), []).append(pr)
    out = []
    for rng, members in by_range.items():
        if rng[1] <= rng[0]:
            out += members
            continue
        pks = sorted((p for p in hs.compound() if rng[0] - pad <= p.d2 <= rng[1] + pad),
                     key=lambda p: (p.kind, p.d1))
        groups = []
        for p in pks:
            if groups and groups[-1][-1].kind == p.kind and p.d1 - groups[-1][-1].d1 <= tolC:
                groups[-1].append(p)
            else:
                groups.append([p])
        from_lines = all(m.source == "H1-lines" for m in members)
        if from_lines:
            lns = [p for p in _h1_lines(sp[H1_1D]) if rng[0] <= p.d1 <= rng[1]]
            symmetric = multiplet_centre(lns, _masked_positions(sp[H1_1D]), sfH)[1] is not None if lns else True
        if from_lines and not symmetric:
            # crude line cluster with no single symmetry centre: one carbon group may still hold
            # protons far apart in f2 (a symmetric cluster is one multiplet, e.g. a dd picked on both halves)
            sub = []
            for g in groups:
                g = sorted(g, key=lambda p: -p.d2)
                cur = [g[0]]
                for p in g[1:]:
                    if (cur[-1].d2 - p.d2) * sfH > P.split_f2_gap_hz:
                        sub.append(cur)
                        cur = []
                    cur.append(p)
                sub.append(cur)
            groups = sub
        if groups:
            gmax = max(sum(abs(p.intensity) for p in g) for g in groups)
            weak = [g for g in groups if sum(abs(p.intensity) for p in g) < P.split_min_frac * gmax]
            c13 = sp.get(C13_1D)
            if weak:
                strong_f1 = [st.mean(p.d1 for p in g) for g in groups if g not in weak]
                cpk = [p.d1 for p in c13.compound()] if c13 is not None else None   # no 13C: HSQC row is the carbon
                res1 = 3 * hs.res_hz(0) / sfC
                def explained_elsewhere(g):
                    f1g, sg = st.mean(p.d1 for p in g), sum(abs(p.intensity) for p in g)
                    return any(abs(q.d1 - f1g) <= tolC and not (rng[0] - pad <= q.d2 <= rng[1] + pad)
                               and abs(q.intensity) > sg for q in hs.compound())
                rescued = [g for g in weak if sum(abs(p.intensity) for p in g) >= P.split_weak_frac * gmax
                           and all(abs(st.mean(p.d1 for p in g) - x) > res1 for x in strong_f1)
                           and (cpk is None or any(abs(st.mean(p.d1 for p in g) - x) <= tolC for x in cpk))
                           and not explained_elsewhere(g)]
                for g in rescued:
                    rep.append(f"Protons: weak HSQC group {st.mean(p.d2 for p in g):.4f}/{st.mean(p.d1 for p in g):.2f} "
                               f"({sum(abs(p.intensity) for p in g) / gmax:.0%} of the strongest) kept - its own carbon, "
                               f"far from the others in f1")
                weak = [g for g in weak if g not in rescued]
            groups = [g for g in groups if g not in weak]
        if len(groups) <= 1:
            out += members
            continue
        cen = [(_wmean([p.d2 for p in g], [abs(p.intensity) for p in g]),
                _wmean([p.d1 for p in g], [abs(p.intensity) for p in g]), g) for g in groups]
        ref = members[0]
        if not any(m.source == "pureshift" for m in members):
            new = []
            cen_s = sorted(cen, key=lambda c: -c[0])
            centres = [None] * len(cen_s)
            strong = any(m.strong for m in members)
            if not from_lines:
                lns = [p for p in _h1_lines(sp[H1_1D]) if rng[0] <= p.d1 <= rng[1]]
            if lns and not strong:  # one symmetry centre per proton, chosen jointly, near its HSQC column
                win = max(hs.res_hz(1), 6.0)    # a line centre must agree with its HSQC row within the f2 resolution
                centres = split_centres(lns, [c[0] for c in cen_s], _masked_positions(sp[H1_1D]), sfH, win)
            for i, (f2, f1, g) in enumerate(cen_s):
                pos, src = f2, "HSQC-row"
                if centres[i] is not None:
                    pos, src = centres[i], "H1-lines"
                if strong and strong_pos and round(f1, 4) in strong_pos:
                    pos, src = strong_pos[round(f1, 4)][0], "HSQC+HMBC"
                pr = Proton("", pos, src, ref.mult_delta, ref.mult_category, ref.range,
                            evidence={"HSQC-row": [f2]} | ({"H1-lines": [pos]} if src == "H1-lines" else {}))
                pr.flags.append(f"split from 1H multiplet {ref.mult_delta:.4f} using HSQC row at {f1:.2f}")
                pr.split_c = f1
                pr.mult_nh = ref.mult_nh
                if strong:
                    pr.strong = True
                    pr.flags.append("strongly coupled region: pureshift unreliable, position from HSQC/HMBC"
                                    + (f" ({strong_pos[round(f1, 4)][1]} HSQC + {strong_pos[round(f1, 4)][2]} HMBC picks)"
                                       if src == "HSQC+HMBC" else ""))
                new.append(pr)
            rep.append(f"Protons: 1H multiplet {ref.mult_delta:.4f} ({ref.mult_category}) contains HSQC peaks at "
                       f"{len(groups)} carbons - split into " +
                       ", ".join(f"{f2:.4f} (C {f1:.2f})" for f2, f1, _ in sorted(cen, key=lambda c: -c[0])))
            out += new
        else:
            out += members
            # degenerate protons: one pureshift proton owning HSQC groups on two distinct carbons
            # is two protons at the same shift (one per carbon)
            ps_m = [m for m in members if m.source == "pureshift"]
            owned = {}
            for f2, f1, g in cen:
                near = min(ps_m, key=lambda m: abs(m.delta - f2))
                if abs(near.delta - f2) <= tolH:
                    owned.setdefault(id(near), (near, []))[1].append((f2, f1, g))
            c13 = sp.get(C13_1D)
            cpk = [p.d1 for p in c13.compound()] if c13 is not None else None
            for m, gs in owned.values():
                if len(gs) < 2:
                    continue
                gs.sort(key=lambda c: -sum(abs(p.intensity) for p in c[2]))
                f1s = [c[1] for c in gs]
                if any(abs(a - b) <= 2 * tolC for i, a in enumerate(f1s) for b in f1s[i + 1:]):
                    continue                     # not clearly different carbons
                if cpk is not None and not all(any(abs(f1 - x) <= tolC for x in cpk) for f1 in f1s):
                    continue                     # each carbon must be a 13C peak
                m.split_c = f1s[0]
                m.flags.append(f"degenerate: proton at the same shift as {len(gs) - 1} other(s) "
                               f"(HSQC carbons {', '.join(f'{x:.2f}' for x in f1s)})")
                for f2, f1, g in gs[1:]:
                    pr = Proton("", m.delta, m.source, m.mult_delta, m.mult_category, m.range,
                                evidence={k: list(v) for k, v in m.evidence.items()})
                    pr.ps_key, pr.mnova_type, pr.mult_nh = m.ps_key, m.mnova_type, m.mult_nh
                    pr.split_c = f1
                    pr.flags.append(f"degenerate: second proton at the pureshift shift {m.delta:.4f} "
                                    f"(HSQC carbons {', '.join(f'{x:.2f}' for x in f1s)}); HMBC/COSY "
                                    f"correlations at this shift are listed under the first proton")
                    out.append(pr)
                rep.append(f"Protons: pureshift {m.delta:.4f} carries HSQC cross-peaks at {len(gs)} distinct carbons "
                           f"({', '.join(f'{x:.2f}' for x in f1s)}; weakest "
                           f"{sum(abs(p.intensity) for p in gs[-1][2]) / sum(abs(p.intensity) for p in gs[0][2]):.0%}"
                           f" of the strongest) - degenerate protons, {len(gs)} protons at the same shift")
            for f2, f1, g in cen:
                if not any(abs(m.delta - f2) <= tolH for m in members):
                    pr = Proton("", f2, "HSQC-row", ref.mult_delta, ref.mult_category, ref.range,
                                evidence={"HSQC-row": [f2]})
                    pr.flags.append(f"no pureshift peak for HSQC carbon {f1:.2f}; position from HSQC row")
                    rep.append(f"Protons: HSQC carbon {f1:.2f} in 1H multiplet {ref.mult_delta:.4f} has no pureshift "
                               f"peak within {tolH * sfH:.1f} Hz - proton added at {f2:.4f} from the HSQC row")
                    out.append(pr)
    return out


def _h1_lines(h1):
    return sorted((p for p in h1.peaks if p.type == COMPOUND and p.intensity > 0
                   and not (p.has("weak") or p.has("c13_satellite") or p.has("hidden"))),
                  key=lambda p: -p.d1)


def _masked_positions(h1, span_hz=8.0):
    """Positions where a compound 1H line can be hidden: solvent/impurity/artifact-typed picks, and
    the residual-solvent table positions +/- span_hz (in case the solvent itself was not picked)."""
    sf = h1.sf[-1]
    out = [(p.d1, abs(p.intensity)) for p in h1.peaks if p.type != COMPOUND]
    for d in solvent_h1_positions(h1.raw.get("solvent", "")):
        out += [(d + k / sf, math.inf) for k in range(-int(span_hz), int(span_hz) + 1)]
    return out


def _hidden_at(m, x, masked, tol):
    """The mirror m of line x can be hidden by a masked line at least half as intense as x."""
    return any(abs(y - m) <= 2 * tol and iy >= 0.5 * abs(x.intensity) for y, iy in masked)


def _centre_candidates(lines, sfH, sub_hz=3.0):
    """Line-pair midpoints, plus midpoints between a line and the centroid of two adjacent lines
    (one half of a multiplet resolved, the mirror half a single broad line)."""
    xs = [p.d1 for p in lines]
    cands = {(xs[i] + xs[j]) / 2 for i in range(len(xs)) for j in range(i, len(xs))}
    subs = [_wmean([a.d1, b.d1], [a.area, b.area]) for a, b in zip(lines, lines[1:])
            if abs(a.d1 - b.d1) * sfH <= sub_hz]
    cands |= {(x + c) / 2 for x in xs for c in subs}
    return sorted(cands)


def _has_mirror(m, lines, tol, x=None, broad=False):
    """A line of comparable intensity (ratio >= 0.2) lies at m, within tol - or, with broad=True,
    within a quarter of the summed widths (broad, partly unresolved lines)."""
    wx = x.width1 if x is not None else 0.0
    ix = abs(x.intensity) if x is not None else 0.0
    for y in lines:
        t = max(tol, 0.25 * (y.width1 + wx)) if broad else tol
        if abs(y.d1 - m) <= t and (x is None or y is x or
                                   min(ix, abs(y.intensity)) >= 0.2 * max(ix, abs(y.intensity))):
            return True
    return False


def _main_lines(lines, frac=0.05):
    """Lines carrying at least frac of the strongest line (tiny lines cannot define a centre)."""
    if not lines:
        return lines
    imax = max(abs(p.intensity) for p in lines)
    return [p for p in lines if abs(p.intensity) >= frac * imax]


def multiplet_centre(lines, masked, sfH, tol_hz=0.6):
    """Strict line tolerance first; broad (width-based) tolerance only if no centre is found."""
    lines = _main_lines(lines)
    c, h = _multiplet_centre(lines, masked, sfH, tol_hz, broad=False)
    if h is None and len(lines) > 1:
        c2, h2 = _multiplet_centre(lines, masked, sfH, tol_hz, broad=True)
        if h2 is not None:
            return c2, h2
    return c, h


def _multiplet_centre(lines, masked, sfH, tol_hz=0.6, broad=False):
    """Centre of a first-order multiplet from its line positions.

    Line intensities are distorted by roofing, so the centre is the point about which the line
    positions are symmetric.  Candidates are midpoints of line pairs; a candidate is valid when
    every line has a mirror line, or its mirror falls on a solvent/impurity line (a line hidden
    under another signal).  Most mirrored lines wins, then fewest hidden lines.
    Returns (centre, n_hidden) or (area-weighted centroid, None) when no symmetric centre exists."""
    xs = [p.d1 for p in lines]
    if len(xs) == 1:
        return xs[0], 0
    tol = tol_hz / sfH
    best = None
    for c in _centre_candidates(lines, sfH):
        matched = hidden = 0
        ok = True
        for ln in lines:
            m = 2 * c - ln.d1
            if _has_mirror(m, lines, tol, ln, broad):
                matched += 1
            elif _hidden_at(m, ln, masked, tol):
                hidden += 1
            else:
                ok = False
                break
        if ok and hidden < matched:
            key = (matched, -hidden)
            if best is None or key > best[0]:
                best = (key, c, hidden)
    if best is None:
        return _wmean(xs, [p.area for p in lines]), None
    return best[1], best[2]


def split_centres(lines, ests, masked, sfH, win_hz, tol_hz=0.6, half_hz=25.0, ncand=12):
    """Strict line tolerance first; broad (width-based) tolerance only if it explains more lines.
    (No small-line filter here: the protons sharing a cluster can differ greatly in intensity.)"""
    a, na = _split_centres(lines, ests, masked, sfH, win_hz, tol_hz, half_hz, ncand, broad=False)
    if na < len(lines):
        b, nb = _split_centres(lines, ests, masked, sfH, win_hz, tol_hz, half_hz, ncand, broad=True)
        if nb > na:
            return b
    return a


def _split_centres(lines, ests, masked, sfH, win_hz, tol_hz=0.6, half_hz=25.0, ncand=12, broad=False):
    """Centres for k protons sharing one 1H line cluster.

    ests are the protons' f2 positions from the HSQC.  Each proton's candidate centres are line-pair
    midpoints within win_hz of its estimate; the combination that explains the most lines (a line
    is explained by a centre when its mirror line, or a solvent/impurity line, lies at 2c - x,
    within half_hz of c) wins, ties broken by the smallest total distance to the estimates.
    Returns a list of centres (None where no candidate)."""
    import itertools
    xs = [p.d1 for p in lines]
    tol, half, win = tol_hz / sfH, half_hz / sfH, win_hz / sfH
    mids = _centre_candidates(lines, sfH)

    def covers(c):
        out = set()
        for i, x in enumerate(xs):
            if abs(x - c) > half:
                continue
            m = 2 * c - x
            if _has_mirror(m, lines, tol, lines[i], broad) or _hidden_at(m, lines[i], masked, tol):
                out.add(i)
        return frozenset(out)

    cands = []
    ncand = ncand if len(ests) > 3 else 3 * ncand
    for e in ests:
        cs = sorted((c for c in mids if abs(c - e) <= win), key=lambda c: abs(c - e))[:ncand]
        cands.append([(c, covers(c)) for c in cs] or [(None, frozenset())])
    def score(combo):
        cov = frozenset().union(*(cv for _, cv in combo))
        dist = sum(abs(c - e) for (c, _), e in zip(combo, ests) if c is not None)
        cs_ = [(c, e) for (c, _), e in zip(combo, ests) if c is not None]
        if any(abs(a[0] - b[0]) < 1e-9 and abs(a[1] - b[1]) > win / 2
               for i, a in enumerate(cs_) for b in cs_[i + 1:]):
            return None    # two protons on the same centre, though the HSQC puts them apart
        cs = [c for c, _ in combo]
        own = 0          # lines explained by the centre they are nearest to
        for i, x in enumerate(xs):
            near = min((j for j, c in enumerate(cs) if c is not None), key=lambda j: abs(cs[j] - x), default=None)
            if near is not None and i in combo[near][1]:
                own += 1
        return (len(cov), own, -dist)

    best = None
    n_combo = 1
    for c in cands:
        n_combo *= len(c)
    if n_combo <= 5000:                  # exhaustive
        for combo in itertools.product(*cands):
            key = score(combo)
            if key is not None and (best is None or key > best[0]):
                best = (key, [c for c, _ in combo])
    else:                                # large clusters: coordinate ascent from the nearest candidates
        cur = [c[0] for c in cands]
        cur_key = score(cur)
        for _ in range(4):
            improved = False
            for i in range(len(cands)):
                for cand in cands[i]:
                    trial = cur[:i] + [cand] + cur[i + 1:]
                    k = score(trial)
                    if k is not None and (cur_key is None or k > cur_key):
                        cur, cur_key, improved = trial, k, True
            if not improved:
                break
        if cur_key is not None:
            best = (cur_key, [c for c, _ in cur])
    return (best[1], best[0][0]) if best else ([None] * len(ests), 0)


def h1_line_multiplets(sp, P: Params, rep: list, covered=None):
    """No Mnova multiplets: group the 1H line picks into pseudo-multiplets.

    With covered (a list of (lo, hi) ranges of the Mnova multiplets) only the lines outside
    those ranges are used ("gap filling"), and a cluster is kept only when the H-C 2D spectra
    support it (an HSQC pick, or HMBC picks from >= 2 carbons, inside its range).

    Compound lines (not weak, not 13C satellites, positive) closer than h1_cluster_gap_hz form one
    multiplet; its position is the area-weighted centroid of the lines and its range spans the
    lines.  Small clusters (< h1_min_cluster_frac x median area) are kept only when an H-C 2D spectrum
    has a compound pick inside them.  Clusters holding several protons are split later by
    split_by_hsqc."""
    h1 = sp[H1_1D]
    sfH = h1.sf[-1]
    lines = _h1_lines(h1)
    if covered is not None:
        pad = P.ps_range_pad_hz / sfH
        lines = [p for p in lines if not any(lo - pad <= p.d1 <= hi + pad for lo, hi in covered)]
    clusters = []
    for ln in lines:
        if clusters and (clusters[-1][-1].d1 - ln.d1) * sfH <= P.h1_cluster_gap_hz:
            clusters[-1].append(ln)
        else:
            clusters.append([ln])
    if not clusters:
        return []
    areas = [sum(p.area for p in c) for c in clusters]
    ref = st.median(areas)
    cols = [pk.d2 for s in sp.values() if s.ndim == 2 and s.nuclei[0] != "1H" for pk in s.compound()]
    out, dropped, notes, asym = [], [], [], []
    masked = _masked_positions(h1)
    for i, (c, a) in enumerate(zip(clusters, areas)):
        # range: line extent + half width (>= 2 Hz), never more than 40 % of the gap to a neighbour
        up = (clusters[i - 1][-1].d1 - c[0].d1) * 0.4 if i else 1e9
        dn = (c[-1].d1 - clusters[i + 1][0].d1) * 0.4 if i + 1 < len(clusters) else 1e9
        lo = c[-1].d1 - min(max(c[-1].width1 / 2, 2.0 / sfH), dn)
        hi = c[0].d1 + min(max(c[0].width1 / 2, 2.0 / sfH), up)
        if covered is not None:
            n_hsqc = sum(lo <= pk.d2 <= hi for pk in sp[HSQC].compound())
            hm = sp.get(HMBC)
            n_hmbc = len({round(pk.d1, 1) for pk in hm.compound() if lo <= pk.d2 <= hi}) if hm else 0
            if not (n_hsqc or n_hmbc >= 2):
                dropped.append(f"{_wmean([p.d1 for p in c], [p.area for p in c]):.4f}")
                continue
        elif a < P.h1_min_cluster_frac * ref and not any(lo <= x <= hi for x in cols):
            dropped.append(f"{_wmean([p.d1 for p in c], [p.area for p in c]):.4f}")
            continue
        centre, hidden = multiplet_centre(c, masked, sfH)
        cat = "m" if hidden else {1: "s", 2: "d", 3: "t", 4: "q"}.get(len(c), "m")
        if hidden:      # the multiplet is symmetric about its centre: extend the range over the hidden lines
            span = max(abs(p.d1 - centre) for p in c) + 2.0 / sfH
            lo, hi = min(lo, centre - span), max(hi, centre + span)
            notes.append(f"{centre:.4f} ({hidden} line(s) under solvent/impurity)")
        elif hidden is None and len(c) > 1:
            asym.append(f"{centre:.4f}")
        out.append({"delta1": centre, "category": cat, "type": COMPOUND, "_range": (lo, hi), "_lines": True,
                    "_nlines": len(c), "_hidden": hidden, "_gap": covered is not None})
        if covered is not None:
            rep.append(f"Protons: 1H signal {centre:.4f} ({len(c)} lines) is not in the Mnova multiplets but has 2D "
                       f"support (HSQC {n_hsqc}, HMBC {n_hmbc} carbons) - added"
                       + (f"; {hidden} line(s) under solvent/impurity" if hidden else ""))
    if covered is not None:
        if dropped:
            rep.append("Protons: compound 1H lines outside the Mnova multiplets without 2D support ignored: "
                       + ", ".join(dropped))
        return out
    rep.append(f"Protons: no Mnova 1H multiplets - {len(out)} multiplets built from {len(lines)} compound 1H lines "
               f"(gap > {P.h1_cluster_gap_hz:g} Hz separates multiplets; weak/13C-satellite lines ignored)"
               + (f"; small clusters without 2D support ignored: {', '.join(dropped)}" if dropped else ""))
    rep.append("Protons: multiplet centres are the symmetry centres of the line positions (intensities are "
               "roofed)" + (f"; completed with hidden lines: {'; '.join(notes)}" if notes else "")
               + (f"; no symmetric centre, area centroid used: {', '.join(asym)}" if asym else ""))
    return out


def merge_strong_multiplets(mults, ints, ps, sfH, P: Params, rep: list):
    """Mark multiplets where the pureshift is unreliable (large negative/dispersive peaks: strong
    coupling) and merge touching unreliable multiplets into one region, split later by HSQC."""
    pad = P.ps_range_pad_hz / sfH
    info = []
    for m in mults:
        rng = m.get("_range") or next(((i["rangeMin1"], i["rangeMax1"]) for i in ints
                                       if i["rangeMin1"] <= m["delta1"] <= i["rangeMax1"]), None)
        bad = False
        if rng is not None and m.get("type", 0) == COMPOUND:
            pk = [p for p in ps.peaks if p.type == COMPOUND and not p.has("weak") and rng[0] - pad <= p.d1 <= rng[1] + pad]
            pos = max((p.intensity for p in pk if p.intensity > 0), default=0)
            neg = max((-p.intensity for p in pk if p.intensity < 0), default=0)
            bad = pos > 0 and neg >= P.ps_dispersive_frac * pos
        info.append((m, rng, bad))
    info.sort(key=lambda x: -x[0]["delta1"])
    out, i = [], 0
    while i < len(info):
        m, rng, bad = info[i]
        if not bad:
            out.append(m)
            i += 1
            continue
        grp = [(m, rng)]
        j = i + 1
        while j < len(info) and info[j][2] and abs(grp[-1][1][0] - info[j][1][1]) * sfH <= 2.0:
            grp.append((info[j][0], info[j][1]))
            j += 1
        nh = sum(g[0].get("nH", 0) or 0 for g in grp) or None
        lo, hi = min(g[1][0] for g in grp), max(g[1][1] for g in grp)
        out.append({"delta1": _wmean([g[0]["delta1"] for g in grp], [g[0].get("nH", 1) or 1 for g in grp]),
                    "category": "m" if len(grp) > 1 else m.get("category", "m"), "nH": nh, "type": COMPOUND,
                    "_range": (lo, hi), "_strong": True})
        rep.append(f"Protons: pureshift unreliable in {lo:.4f}-{hi:.4f} (dispersive negative peaks - strong "
                   f"coupling)" + (f"; multiplets " + ", ".join(f"{g[0]['delta1']:.4f}" for g in grp)
                                   + " merged" + (f" ({nh}H)" if nh else "") if len(grp) > 1 else "")
                   + "; protons from HSQC carbons, positions from HSQC/HMBC")
        i = j
    return out


def impurity_carbons(sp, protons, carbons, hmbc_corr, P):
    """Weak carbons (< c13_minor_frac of the reference) with no compound proton and no HMBC
    correlation, whose HSQC cross-peaks all lie under a compound proton that already sits on
    its own (stronger) carbon: a protonated impurity carbon hidden under a compound signal.
    (A real quaternary carbon has no HSQC cross-peak at all, so it is never caught here.)"""
    c13, hs = sp.get(C13_1D), sp.get(HSQC)
    if c13 is None or hs is None or not getattr(c13, "ref", None):
        return {}
    sfC, sfH = hs.sf
    tolC = P.hsqc_tol_pts * hs.pt_hz(0) / sfC
    pad = P.ps_range_pad_hz / sfH
    inten = {p.key: abs(p.intensity) for p in c13.peaks}
    out = {}
    for c in carbons:
        if c.protons or c.peak_key is None or inten.get(c.peak_key, 0) >= P.c13_minor_frac * c13.ref:
            continue
        if any(x["c"] is c for x in hmbc_corr):
            continue
        picks = [p for p in hs.compound() if abs(p.d1 - c.delta) <= tolC
                 and min(carbons, key=lambda q: abs(q.delta - p.d1)) is c]
        if not picks:
            continue                    # no HSQC cross-peak: quaternary candidate, keep
        owners = [next((h for h in protons if h.range[0] - pad <= p.d2 <= h.range[1] + pad), None) for p in picks]
        if all(o is not None and o.carbons for o in owners):
            out[c.peak_key] = (f"weak ({inten[c.peak_key] / c13.ref:.0%} of reference), no HMBC; its HSQC cross-peak "
                               f"{picks[0].d2:.4f}/{picks[0].d1:.2f} lies under {owners[0].id} {owners[0].delta:.4f}, which "
                               f"sits on its own carbon - impurity hidden under a compound signal")
    return out


def strong_positions(protons, hmbc_corr, sfH):
    """Strongly coupled protons: position = mean of their HSQC rows and HMBC columns (every
    assigned HMBC correlation, single lines included)."""
    out = {}
    for h in protons:
        if not h.strong or h.split_c is None:
            continue
        hs = h.evidence.get(HSQC, [])
        hb = [h.delta + c["dH"] / sfH for c in hmbc_corr if c["h"] is h]
        if hs or hb:
            out[round(h.split_c, 4)] = (st.mean(hs + hb), len(hs), len(hb))
    return out


def build_protons(sp, P: Params, rep: list, exclude=frozenset(), strong_pos=None):
    h1, ps = sp.get(H1_1D), sp.get(H1_PURESHIFT)
    if h1 is None:
        raise ValueError("H1_1D spectrum required")
    sfH = h1.sf[-1]
    ints = list(h1.raw["integrals"]["data"].values())
    norm = h1.raw["integrals"].get("normValue", 1.0) or 1.0
    mults = [m for m in h1.raw["multiplets"]["data"].values() if m.get("type", 0) == COMPOUND]
    from_lines = not mults
    if from_lines:
        mults = h1_line_multiplets(sp, P, rep)
    else:                       # gap filling: signals Mnova has no multiplet for, but the 2D supports
        covered = []
        for m in mults:
            r = next(((i["rangeMin1"], i["rangeMax1"]) for i in ints if i.get("type", 0) == COMPOUND
                      and i["rangeMin1"] <= m["delta1"] <= i["rangeMax1"]), None)
            covered.append(r or (m["delta1"] - 2.0 / sfH, m["delta1"] + 2.0 / sfH))
        mults = mults + h1_line_multiplets(sp, P, rep, covered=covered)
    if ps is not None:
        mults = merge_strong_multiplets(mults, ints, ps, sfH, P, rep)

    ps_cands = []
    if ps is not None:
        ok_types = (COMPOUND, SOLVENT, IMPURITY) if P.reclassify_ps else (COMPOUND,)
        ps_cands = [p for p in ps.peaks if p.type in ok_types and not p.has("weak")
                    and p.key not in exclude and p.intensity > 0]   # negative = dispersive residue

    ps_tol, ps_split = P.ps_merge_hz, False
    if ps_cands:
        top = max(p.intensity for p in ps_cands)
        strong = [p for p in ps_cands if p.type == COMPOUND and p.intensity >= 0.05 * top]
        ps_tol, frac, nst = split_tolerance([p.d1 * sfH for p in strong], [p.intensity for p in strong],
                                            P.ps_merge_hz)
        ps_split = ps_tol > P.ps_merge_hz
        if ps_split:
            rep.append(f"Pureshift WARNING: {frac:.0%} of the {nst} strong peaks are picked as 2-3 lines within "
                       f"{ps_tol:.1f} Hz - distorted lineshape (check the shimming / line fit); lines of one "
                       f"resonance merged at their intensity-weighted centre")
    protons, used_ps = [], set()
    for m in mults:
        if m.get("type", 0) != COMPOUND:
            continue
        rng = m.get("_range") or next(((i["rangeMin1"], i["rangeMax1"]) for i in ints
                    if i["rangeMin1"] <= m["delta1"] <= i["rangeMax1"]), None)
        if rng is None:
            w = 0.02
            rng = (m["delta1"] - w, m["delta1"] + w)
        pad = P.ps_range_pad_hz / sfH
        inside = [p for p in ps_cands if rng[0] - pad <= p.d1 <= rng[1] + pad]
        if m.get("_strong"):
            used_ps |= {p.key for p in inside}
            inside = []
        found = []
        if inside:
            amax = max(p.area for p in inside)
            found = [p for p in inside if p.area >= P.ps_area_frac * amax]
            kept_ps = []
            for cl in merge_clusters(found, [p.d1 * sfH for p in found], ps_tol, 2.5 * ps_tol):
                c = _centroid_copy(cl) if ps_split and len(cl) > 1 else cl[0]
                kept_ps.append(c)
                for p in cl[1:]:
                    used_ps.add(p.key)
                if len(cl) > 1:
                    rep.append(f"Pureshift: {', '.join(f'{p.d1:.4f}' for p in sorted(cl, key=lambda p: -p.d1))} "
                               f"merged into {c.d1:.4f} (span {(max(p.d1 for p in cl) - min(p.d1 for p in cl)) * sfH:.1f} Hz "
                               f"in one multiplet)")
            found = kept_ps
        if found:
            for p in sorted(found, key=lambda p: -p.d1):
                if p.key in used_ps:
                    continue
                used_ps.add(p.key)
                pr = Proton("", p.d1, "pureshift", m["delta1"], m["category"], rng,
                            evidence={"pureshift": [p.d1]}, mnova_type=p.type)
                pr.ps_key = p.key
                protons.append(pr)
        else:
            src = "H1-lines" if m.get("_lines") else "H1-multiplet"
            pr = Proton("", m["delta1"], src, m["delta1"], m["category"], rng,
                        evidence={src: [m["delta1"]]})
            pr.hidden_lines = m.get("_hidden") or 0
            if m.get("_gap"):
                pr.flags.append("not in the Mnova multiplets - added from the 1H lines (2D support)")
            if m.get("_strong"):
                pr.strong = True
                pr.flags.append("strongly coupled region: pureshift unreliable, position from HSQC/HMBC")
            elif ps is not None:
                pr.flags.append("no pureshift peak (broad/exchangeable or missing)")
            pr.mult_nh = m.get("nH")
            protons.append(pr)
    for pr in protons:
        if pr.mult_nh is None:
            pr.mult_nh = next((m.get("nH") for m in mults if abs(m["delta1"] - pr.mult_delta) < 1e-9), None)

    protons = split_by_hsqc(sp, protons, P, rep, strong_pos)
    protons.sort(key=lambda p: -p.delta)
    for i, p in enumerate(protons):
        p.id = f"H{i + 1}"

    # --- 1H line centroids and proton counts --------------------------------
    # solvent-typed lines count only inside Mnova multiplet ranges (never for multiplets built from lines)
    line_ranges = {tuple(m["_range"]) for m in mults if m.get("_lines")}
    lines = [p for p in h1.peaks if p.type == COMPOUND or
             (p.type == SOLVENT and not from_lines and not p.has("c13_satellite") and
              any(pr.range[0] <= p.d1 <= pr.range[1] for pr in protons if tuple(pr.range) not in line_ranges))]
    tol = P.h1_line_assign_hz / sfH
    groups = {pr.id: [] for pr in protons}
    for ln in lines:
        near = min(protons, key=lambda pr: abs(pr.delta - ln.d1))
        if abs(near.delta - ln.d1) <= tol:
            groups[near.id].append(ln)
    # integrals dominated by a solvent line (e.g. one drawn over DMSO) cannot calibrate the counts
    def solvent_dominated(i):
        inside = [p for p in h1.peaks if i["rangeMin1"] <= p.d1 <= i["rangeMax1"]]
        sv = max((abs(p.intensity) for p in inside if p.type == SOLVENT), default=0)
        return sv > max((abs(p.intensity) for p in inside if p.type == COMPOUND), default=0)
    cal = [i for i in ints if i.get("type", 0) == COMPOUND and not solvent_dominated(i)]
    tot_area = sum(ln.area for g in groups.values() for ln in g
                   if any(i["rangeMin1"] <= ln.d1 <= i["rangeMax1"] for i in cal))
    tot_H = sum(i["intensity"] / norm for i in cal)
    k = tot_area / tot_H if tot_H else None
    if not k:      # no integrals: relative counts, the smallest resolved single-proton signals = 1H
        solo = [sum(ln.area for ln in groups[pr.id]) for pr in protons
                if groups[pr.id] and pr.source != "HSQC-row" and "split" not in " ".join(pr.flags)]
        solo = [a for a in solo if a >= 0.25 * st.median(solo)] if solo else solo
        if solo:
            amin = min(solo)
            k = st.median([a for a in solo if a < 1.5 * amin])
            rep.append("Protons: no 1H integrals - nH is relative (smallest single-proton signals set to 1H)")
    for pr in protons:
        g = groups[pr.id]
        if g:
            a = sum(ln.area for ln in g)
            pr.h1_centroid = sum(ln.d1 * ln.area for ln in g) / a
            pr.nH_est = a / k if k else None
            pr.evidence.setdefault("H1-centroid", []).append(pr.h1_centroid)
    # degenerate protons (same shift): the 1H lines cannot be divided - share the area equally
    for pr in protons:
        same = [q for q in protons if abs(q.delta - pr.delta) < 1e-9]
        if len(same) > 1 and pr is same[0] and pr.nH_est:
            share = pr.nH_est / len(same)
            for q in same:
                q.nH_est, q.h1_centroid = share, pr.h1_centroid
    # integrals normalised to n x 1H (e.g. a 3H methyl set to 9): rescale if all counts fit
    est = [pr.nH_est for pr in protons if pr.nH_est and pr.nH_est >= 0.5]
    if len(est) >= 3 and min(est) >= 1.6:
        r = round(min(est))
        fit = sum(abs(e / r - round(e / r)) <= 0.3 for e in est) / len(est)
        if fit >= 0.75:
            for pr in protons:
                if pr.nH_est:
                    pr.nH_est /= r
            h1.nh_scale = r
            rep.append(f"Protons: 1H integrals appear normalised to {r} x 1H (smallest signals ~{min(est):.1f}) - "
                       f"nH divided by {r}; check the integral normalisation in Mnova")

    # --- pureshift peaks that look real but are not used ---------------------
    if ps is not None:
        res_areas = [p.area for p in ps_cands if p.key in used_ps]
        ref = st.median(res_areas) if res_areas else 0
        pad = P.ps_range_pad_hz / sfH
        for p in ps.compound():
            if p.key in used_ps or p.has("weak") or p.area < 0.1 * ref:
                continue
            if any(pr.range[0] - pad <= p.d1 <= pr.range[1] + pad for pr in protons):
                continue          # residue inside a compound multiplet
            h1hit = min(h1.peaks, key=lambda q: abs(q.d1 - p.d1))
            dh = (h1hit.d1 - p.d1) * sfH
            rep.append(f"Classification: pureshift {p.d1:.4f} typed compound is not inside any compound "
                       f"1H multiplet (nearest 1H peak {h1hit.d1:.4f}, {dh:+.1f} Hz, typed "
                       f"{TYPE_NAMES[h1hit.type]}) - excluded from compound protons")
    return protons


def pureshift_sidebands(sp, protons, P: Params, rep: list):
    ps, h1 = sp.get(H1_PURESHIFT), sp.get(H1_1D)
    if ps is None:
        return None
    sf = ps.sf[-1]
    res = [p.delta for p in protons if p.source == "pureshift"]
    offs = []
    for p in ps.peaks:
        if p.type in (ARTIFACT, IMPURITY):
            for r in res:
                o = (p.d1 - r) * sf
                if 35 <= abs(o) <= 75:
                    offs.append(abs(o))
    if not offs:
        rep.append("Pureshift sidebands: none typed by Mnova; spacing not measured")
        return None
    spacing = st.median(offs)
    rep.append(f"Pureshift sideband spacing measured from Mnova-typed sidebands: {spacing:.1f} Hz "
               f"({len(offs)} offsets)")
    for p in ps.compound():
        for r in res:
            for kk in (1, 2):
                o = (p.d1 - r) * sf
                if abs(abs(o) - kk * spacing) <= P.sideband_tol_hz:
                    h1line = [q for q in h1.peaks if abs(q.d1 - p.d1) * sf < 3] if h1 else []
                    if any(q.type == COMPOUND for q in h1line):
                        continue  # real compound 1H line there: multiplet residue, not a sideband
                    extra = (" - but a 1H line exists there (" +
                             ", ".join(f"{q.d1:.4f} {TYPE_NAMES[q.type]}" for q in h1line) + ")") if h1line else ""
                    rep.append(f"  compound-typed pureshift {p.d1:.4f} = {r:.4f} {o:+.1f} Hz "
                               f"(order {kk} sideband){extra}")
    return spacing


# --------------------------------------------------------------------------- #
# Carbons
# --------------------------------------------------------------------------- #
@dataclass
class Carbon:
    id: str
    delta: float
    mult: str = ""            # CH / CH2 / Cq? from HSQC
    protons: list = field(default_factory=list)
    source: str = "C13"       # where delta came from: C13 (1D), HSQC or HMBC (f1)
    peak_key: str | None = None   # key of the 13C 1D peak it came from


def _wmean(vals, weights):
    w = sum(weights)
    return sum(v * x for v, x in zip(vals, weights)) / w if w else st.mean(vals)


def carbons_from_2d(sp, protons, P, rep):
    """No 13C 1D: carbon list from HSQC f1 (protonated) and HMBC f1 rows (quaternary).

    HSQC peaks are grouped by the proton they belong to, so carbons closer than the f1
    resolution stay separate when their protons differ (e.g. 124.45/8.64 vs 124.51/7.69).
    A positive (CH/CH3) carbon cannot carry two different protons, so those are never
    merged; two negative (CH2) groups on different protons at the same f1 are merged as a
    diastereotopic pair.  Quaternaries come from HMBC rows with no HSQC carbon."""
    hs, hm, h1 = sp.get(HSQC), sp.get(HMBC), sp.get(H1_1D)
    carbons = []
    sfH1 = h1.sf[-1] if h1 else 600.0
    solvH = [p.d1 for p in h1.peaks if p.type == SOLVENT] if h1 else []
    solvC = solvent_c13_positions((h1.raw if h1 else {}).get("solvent", ""))
    tolC_hs = 0.0
    if hs is not None and hs.compound():
        sfC, sfH = hs.sf
        tolC_hs = P.hsqc_tol_pts * hs.pt_hz(0) / sfC
        hgroups, _, _ = hsqc_groups(hs, protons, P)
        best = {}
        for g in hgroups:
            if g["h"] is None:
                continue
            if solvH and min(abs(g["f2"] - x) for x in solvH) * sfH < 5 and \
                    solvC and min(abs(g["f1"] - x) for x in solvC) <= tolC_hs:
                continue
            # a proton sits on one carbon: keep its strongest HSQC group, the rest are wiggles/noise
            if g["h"].id not in best or g["strength"] > best[g["h"].id]["strength"]:
                best[g["h"].id] = g
        dropped = [g for g in hgroups if g["h"] is not None and best.get(g["h"].id) is not g]
        if dropped:
            rep.append("13C from HSQC: weaker HSQC groups of already-placed protons not used as carbons: " +
                       ", ".join(f"{g['h'].id} {g['f1']:.2f}" for g in dropped))
        cands = [dict(f1=g["f1"], kind=g["kind"], hids={g["h"].id}, n=len(g["peaks"])) for g in best.values()]
        # merge diastereotopic CH2 pairs
        merge_tol = P.ch2_merge_frac * hs.res_hz(0) / sfC
        cands.sort(key=lambda c: c["f1"])
        merged = []
        for c in cands:
            m = next((o for o in merged if o["kind"] == 2 and c["kind"] == 2
                      and abs(o["f1"] - c["f1"]) <= merge_tol and not (o["hids"] & c["hids"])), None)
            if m:
                m["f1"] = (m["f1"] * m["n"] + c["f1"] * c["n"]) / (m["n"] + c["n"])
                m["n"] += c["n"]; m["hids"] |= c["hids"]
                rep.append(f"13C from HSQC: CH2 at {m['f1']:.2f} carries protons {'/'.join(sorted(m['hids']))} "
                           f"(diastereotopic pair)")
            else:
                merged.append(c)
        for c in merged:
            carbons.append(Carbon("", c["f1"], source=HSQC))
        close = [(a, b) for i, a in enumerate(merged) for b in merged[i + 1:]
                 if abs(a["f1"] - b["f1"]) * sfC < hs.res_hz(0)]
        for a, b in close:
            rep.append(f"13C from HSQC: carbons {a['f1']:.2f} ({'/'.join(sorted(a['hids']))}) and "
                       f"{b['f1']:.2f} ({'/'.join(sorted(b['hids']))}) are closer than the f1 resolution - "
                       f"kept separate because their protons differ")
    if hm is not None and hm.compound():
        sfC = hm.sf[0]
        tolC_hm = max(P.hmbc_f1_tol_pts * hm.pt_hz(0), P.hmbc_f1_tol_res * hm.res_hz(0)) / sfC
        pks = sorted((p for p in hm.compound()
                      if not (solvC and min(abs(p.d1 - x) for x in solvC) <= tolC_hm)), key=lambda p: p.d1)
        rows = []
        for p in pks:
            if rows and p.d1 - rows[-1][-1].d1 <= tolC_hm:
                rows[-1].append(p)
            else:
                rows.append([p])
        known = [c.delta for c in carbons]
        tol = max(tolC_hs, tolC_hm)
        sfHm = hm.sf[1]
        res1, res2 = hm.res_hz(0), hm.res_hz(1)
        hmap = {h.id: h for h in protons}
        hs_pairs = [(c.delta, hid) for c, cd in zip(carbons, merged) for hid in cd["hids"] if hid in hmap] \
            if hs is not None and hs.compound() else []
        pad_fp = (P.ps_range_pad_hz + res2 / 2) / sfHm
        ptol = max(P.hmbc_presence_hz, 1.5 * res2) / sfHm

        def at_proton(p):
            """proton a pick belongs to, or None (also None for 1J-breakthrough positions)"""
            for cf1, hid in hs_pairs:
                if abs(cf1 - p.d1) * sfC <= res1:
                    o = abs(p.d2 - hmap[hid].delta) * sfHm
                    if P.onebond_min_hz / 2 <= o <= P.onebond_max_hz / 2:
                        return None
            h = min(protons, key=lambda h: abs(h.delta - p.d2)) if protons else None
            if h is not None and abs(h.delta - p.d2) <= ptol:
                return h
            return proton_owner(protons, p.d2, pad_fp, 0.0)

        for r in rows:
            f1 = _wmean([p.d1 for p in r], [abs(p.intensity) for p in r])
            if known and min(abs(f1 - x) for x in known) <= tol:
                continue
            at = [h for h in (at_proton(p) for p in r) if h is not None]
            if len(at) < P.hmbc_q_min_peaks:
                rep.append(f"13C from HMBC: row at f1 {f1:.2f} ({len(r)} pick(s), {len(at)} at a proton and not "
                           f"1J breakthrough) not used as a carbon")
                continue
            carbons.append(Carbon("", f1, source=HMBC))
    carbons.sort(key=lambda c: -c.delta)
    for i, c in enumerate(carbons):
        c.id = f"C{i + 1}"
    if carbons:
        rep.append(f"13C: no 1D carbon spectrum - {sum(c.source == HSQC for c in carbons)} carbons from HSQC f1, "
                   f"{sum(c.source == HMBC for c in carbons)} quaternary from HMBC f1 "
                   f"(precision limited by 2D f1 resolution)")
    return carbons


def solvent_carbon_evidence(sp, protons, f1, P, others=()):
    """Reason to keep a compound-typed 13C peak that sits on the solvent, or None: an HSQC
    cross-peak from a compound proton clear of the solvent 1H, or HMBC picks from at least
    two such protons."""
    for name in (HSQC, HMBC):
        spec = sp.get(name)
        if spec is None:
            continue
        sfC, sfH = spec.sf
        solvH = solvent_h1_positions(spec.raw.get("solvent", ""))
        if name == HSQC:
            tolC = P.hsqc_tol_pts * spec.pt_hz(0) / sfC
        else:
            tolC = max(P.hmbc_f1_tol_pts * spec.pt_hz(0), P.hmbc_f1_tol_res * spec.res_hz(0)) / sfC
        if any(abs(o - f1) <= tolC for o in others):
            continue           # a stronger carbon is within tolerance: the 2D cannot tell them apart
        pad = P.ps_range_pad_hz / sfH
        hs = set()
        for p in spec.peaks:
            if p.type not in (COMPOUND, SOLVENT) or abs(p.d1 - f1) > tolC:
                continue
            if any(abs(p.d1 - o) < abs(p.d1 - f1) for o in others):
                continue           # this pick is nearer another carbon
            if solvH and min(abs(p.d2 - x) for x in solvH) * sfH < P.solvent_h_clear_hz:
                continue
            h = proton_owner(protons, p.d2, pad, 0.0)
            if h is not None:
                hs.add(h.id)
        if name == HSQC and hs:
            return "HSQC cross-peak from " + "/".join(sorted(hs))
        if name == HMBC and len(hs) >= 2:
            return "HMBC correlations from " + "/".join(sorted(hs))
    return None


def build_carbons(sp, P: "Params" = None, rep: list | None = None, protons=None, exclude_c=None):
    c13 = sp.get(C13_1D)
    if c13 is None:
        return carbons_from_2d(sp, protons or [], P or Params(), rep if rep is not None else [])
    P = P or Params()
    rep = rep if rep is not None else []
    sf = c13.sf[-1]
    c13.excluded = {}       # peak key -> (new type, reason), used to retype peaks in the written 13C list
    # solvent positions: Mnova-typed solvent peaks plus the table value for the stated solvent
    solv = [p.d1 for p in c13.peaks if p.type == SOLVENT]
    solv += solvent_c13_positions(c13.raw.get("solvent", ""))
    cmp, near_solvent = [], []
    for p in c13.compound():
        if solv and min(abs(p.d1 - x) for x in solv) <= P.c13_solvent_tol_ppm:
            near_solvent.append(p)
        else:
            cmp.append(p)
    if near_solvent:
        kept_ns = []
        for p in near_solvent:
            why = solvent_carbon_evidence(sp, protons or [], p.d1, P)
            if why:
                kept_ns.append(p)
                rep.append(f"13C: {p.d1:.3f} lies at the solvent position but is kept - {why}")
        near_solvent = [p for p in near_solvent if p not in kept_ns]
        cmp += kept_ns
    for p in near_solvent:
        c13.excluded[p.key] = (SOLVENT, "at the solvent position, no 2D evidence for a compound carbon")
    if near_solvent:
        rep.append("13C: compound-typed peaks at solvent positions not used: " +
                   ", ".join(f"{p.d1:.3f}" for p in near_solvent) + " - check by hand if the compound has a carbon there")
    if not cmp:
        return []
    # reference intensity: median of the strongest few compound peaks, so neither the solvent
    # nor a single very tall peak sets the threshold (important for dilute samples)
    top = sorted((abs(p.intensity) for p in cmp), reverse=True)[:max(1, P.c13_ref_n)]
    ref = st.median(top)
    c13.ref = ref
    thr = P.c13_min_frac * ref
    keep = [p for p in cmp if abs(p.intensity) >= thr]
    drop = [p for p in cmp if abs(p.intensity) < thr]
    for p in drop:
        c13.excluded[p.key] = (ARTIFACT, "below the 13C threshold (13C-13C satellite / noise)")
    weakest = min(abs(p.intensity) for p in keep)
    rep.append(f"13C: reference intensity {ref:.3g} (median of {len(top)} strongest compound peaks), "
               f"threshold {thr:.3g}; {len(keep)} carbons kept, weakest kept at {weakest / ref:.1%} of reference")
    if drop:
        notes = []
        for p in sorted(drop, key=lambda p: -p.d1):
            par = min(keep, key=lambda q: abs(q.d1 - p.d1))
            notes.append(f"{p.d1:.3f} ({(p.d1 - par.d1) * sf:+.0f} Hz from {par.d1:.3f}, "
                         f"{abs(p.intensity) / abs(par.intensity):.1%})")
        rep.append(f"13C: {len(drop)} weak compound-typed peaks ignored "
                   f"(13C-13C satellites / noise): " + ", ".join(notes))
    # gap check: warn if the kept/dropped boundary is not a clear gap
    if drop:
        strongest_drop = max(abs(p.intensity) for p in drop)
        if weakest < 5 * strongest_drop:
            rep.append(f"13C WARNING: weakest kept peak only {weakest / strongest_drop:.1f}x the strongest "
                       f"ignored one - threshold may be cutting real (e.g. quaternary) carbons")
    strong = [p for p in keep if abs(p.intensity) >= 0.1 * ref]
    tol, frac, nst = split_tolerance([p.d1 * sf for p in strong], [abs(p.intensity) for p in strong],
                                     P.c13_merge_hz)
    c13.split_tol = tol if tol > P.c13_merge_hz else None
    if c13.split_tol:
        rep.append(f"13C WARNING: {frac:.0%} of the {nst} strong peaks are picked as 2-3 lines within {tol:.1f} Hz - "
                   f"distorted lineshape (check the shimming / line fit); each cluster merged into one carbon at "
                   f"its intensity-weighted centre")
    merged_c = []
    for cl in merge_clusters(keep, [p.d1 * sf for p in keep], tol, 2.5 * tol):
        if len(cl) == 1:
            merged_c.append(cl[0])
            continue
        c = _centroid_copy(cl) if c13.split_tol else cl[0]
        merged_c.append(c)
        for p in cl[1:]:
            c13.excluded[p.key] = (ARTIFACT, f"line-fit split of the carbon at {c.d1:.4f}")
        rep.append(f"13C: {', '.join(f'{p.d1:.4f}' for p in sorted(cl, key=lambda p: -p.d1))} merged into "
                   f"{c.d1:.4f} (span {(max(p.d1 for p in cl) - min(p.d1 for p in cl)) * sf:.1f} Hz - "
                   f"one carbon split by the line fit)")
    keep = merged_c
    if exclude_c:                     # impurity carbons found after a first HSQC/HMBC pass
        for p in [p for p in keep if p.key in exclude_c]:
            c13.excluded[p.key] = (IMPURITY, exclude_c[p.key])
        keep = [p for p in keep if p.key not in exclude_c]
    # minor peaks: far weaker than every other carbon (clear intensity gap) -> impurity unless a
    # 2D spectrum shows a compound proton on/near it and no stronger carbon is within its tolerance
    minor = sorted((p for p in keep if abs(p.intensity) < P.c13_minor_frac * ref), key=lambda p: -abs(p.intensity))
    major = [p for p in keep if p not in minor]
    if minor and major and min(abs(p.intensity) for p in major) >= P.c13_minor_gap * abs(minor[0].intensity):
        excl = []
        for p in minor:
            why = solvent_carbon_evidence(sp, protons or [], p.d1, P, others=[q.d1 for q in major])
            if why:
                rep.append(f"13C: minor peak {p.d1:.3f} ({abs(p.intensity) / ref:.1%}) kept - {why}")
            else:
                excl.append(p)
        for p in excl:
            c13.excluded[p.key] = (IMPURITY, "minor peak with no 2D support")
        if excl:
            keep = [p for p in keep if p not in excl]
            rep.append(f"13C: {len(excl)} minor peaks (<= {abs(minor[0].intensity) / min(abs(q.intensity) for q in major):.0%} "
                       f"of the weakest other carbon, no 2D support) treated as impurity: "
                       + ", ".join(f"{p.d1:.3f}" for p in sorted(excl, key=lambda p: -p.d1)))
    keep = sorted(keep, key=lambda p: -p.d1)
    out = [Carbon(f"C{i + 1}", p.d1) for i, p in enumerate(keep)]
    for c, p in zip(out, keep):
        c.peak_key = p.key
    return out


# --------------------------------------------------------------------------- #
# HSQC
# --------------------------------------------------------------------------- #
def new_peak(d1, d2, intensity, kind, width1, width2, annotation, type_=COMPOUND, flags=0):
    return {"delta1": d1, "delta2": d2, "intensity": intensity, "type": type_, "kind": kind,
            "width1": width1, "width2": width2, "annotation": annotation, "flags": flags}


def hsqc_extra_picks(hs, protons, carbons, P: Params, rep=None):
    """Solvent-typed HSQC picks that belong to the compound: owned by a compound proton (1H
    footprint), clear of the residual-solvent 1H position and on a compound carbon row."""
    if not carbons:
        return []
    sfC, sfH = hs.sf
    tolC = P.hsqc_tol_pts * hs.pt_hz(0) / sfC
    pad = P.ps_range_pad_hz / sfH
    solvH = solvent_h1_positions(hs.raw.get("solvent", ""))
    extra = []
    for p in hs.peaks:
        if p.type != SOLVENT:
            continue
        if solvH and min(abs(p.d2 - x) for x in solvH) * sfH < P.solvent_h_clear_hz:
            continue
        h = proton_owner(protons, p.d2, pad, 0.0)
        if h is None or not any(abs(c.delta - p.d1) <= tolC for c in carbons):
            continue
        extra.append(p)
        if rep is not None:
            rep.append(f"HSQC {p.d2:.4f}/{p.d1:.3f}: typed solvent by Mnova but inside {h.id} {h.delta:.4f} and "
                       f"{min(abs(p.d2 - x) for x in solvH) * sfH:.0f} Hz from the solvent 1H - used as compound")
    return extra


def hsqc_groups(hs, protons, P: Params, extra=()):
    """Group HSQC compound picks by owner proton (1H footprint), sign and f1.

    Returns a list of dicts: h (Proton or None), kind, f1, f2 (intensity-weighted), strength
    (sum |I|), intensity (signed sum), peaks."""
    sfC, sfH = hs.sf
    tolH = proton_tol_ppm(hs)
    tolH_res = max(P.split_tol_pts * hs.pt_hz(1), P.split_min_tol_hz) / hs.sf[1]   # resolution test (see Params)
    tolC = P.hsqc_tol_pts * hs.pt_hz(0) / sfC
    pad = P.ps_range_pad_hz / sfH
    byowner = {}
    for pk in list(hs.compound()) + list(extra):
        h = proton_owner(protons, pk.d2, pad, tolH)
        if h is not None:   # protons split out of one multiplet by HSQC carbon: the carbon decides
            sibs = [q for q in protons if q.split_c is not None and tuple(q.range) == tuple(h.range)
                    and abs(q.split_c - pk.d1) <= tolC]
            if sibs:
                h = min(sibs, key=lambda q: abs(q.split_c - pk.d1))
        key = (h.id if h else None, pk.kind)
        byowner.setdefault(key, (h, []))[1].append(pk)
    groups = []
    for (hid, kind), (h, pks) in byowner.items():
        pks.sort(key=lambda p: p.d1)
        sub = [[pks[0]]]
        for p in pks[1:]:
            same_f2 = h is not None or abs(p.d2 - sub[-1][-1].d2) <= 2 * tolH_res
            if p.d1 - sub[-1][-1].d1 <= tolC and same_f2:
                sub[-1].append(p)
            else:
                sub.append([p])
        for g in sub:
            w = [abs(p.intensity) for p in g]
            groups.append(dict(h=h, kind=kind, peaks=g, strength=sum(w), intensity=sum(p.intensity for p in g),
                               f1=_wmean([p.d1 for p in g], w), f2=_wmean([p.d2 for p in g], w)))
    return groups, tolH, tolC


def refine_hsqc(sp, protons, carbons, P: Params, rep: list):
    """Snap the HSQC to the proton and carbon lists.

    Picks are grouped by owner proton, sign and f1 (multiplet line picks become one group).
    Groups are accepted strongest first under the structural rules
      - a proton is attached to exactly one carbon,
      - a CH/CH3 (positive) carbon carries one proton, a CH2 (negative) carbon at most two,
      - all protons on one carbon have the same sign;
    a carbon is chosen from those within the f1 tolerance (nearest free one).  Groups that
    break a rule are retyped artifact with the reason (truncation wiggles, t1 noise, ...)."""
    hs = sp.get(HSQC)
    if hs is None:
        return None, []
    sfC, sfH = hs.sf
    extra = hsqc_extra_picks(hs, protons, carbons, P, rep)
    groups, tolH, tolC = hsqc_groups(hs, protons, P, extra)
    rep.append(f"HSQC tolerances: H {tolH * sfH:.1f} Hz ({tolH:.4f} ppm), "
               f"C {tolC * sfC:.1f} Hz ({tolC:.3f} ppm); f1 resolution {hs.res_hz(0):.0f} Hz")
    w1 = [p.width1 * sfC for p in hs.compound() if p.width1 > 0]
    if w1 and st.median(w1) < 0.5 * hs.res_hz(0):
        rep.append(f"HSQC: f1 linewidths (median {st.median(w1):.0f} Hz) are much narrower than 1/AQ "
                   f"({hs.res_hz(0):.0f} Hz) - linear prediction or similar applied in f1; carbons are separated "
                   f"by the picks down to {0.5 * hs.pt_hz(0):.0f} Hz")
    out, matches = [], []
    on_carbon = {}            # carbon id -> list of accepted groups
    accepted = {}             # proton id -> group
    n_dup = 0
    for g in sorted(groups, key=lambda g: -g["strength"]):
        h = g["h"]
        cands = sorted((c for c in carbons if abs(c.delta - g["f1"]) <= tolC), key=lambda c: abs(c.delta - g["f1"]))
        reason = None
        chosen = None
        if h is not None and h.id in accepted:
            a = accepted[h.id]
            reason = (f"{h.id} already on {a['c'].id} ({g['strength'] / a['strength']:.0%} of that peak)"
                      + (" - f1 truncation wiggle?" if abs(g["f1"] - a["c"].delta) * sfC <= P.wiggle_res * hs.res_hz(0)
                         else ""))
        elif h is not None and cands:
            for c in cands:
                have = on_carbon.get(c.id, [])
                cap = 2 if g["kind"] == 2 else 1
                if len(have) < cap and all(x["kind"] == g["kind"] for x in have):
                    chosen = c
                    break
            if chosen is None:
                c = cands[0]
                reason = f"{c.id} already carries " + "/".join(x["h"].id for x in on_carbon[c.id]) + \
                         f" ({'CH2' if on_carbon[c.id][0]['kind'] == 2 else 'CH/CH3'})"
        if chosen is not None:
            g["c"] = chosen
            # carbons the HSQC cannot tell apart (closer than half its f1 resolution)
            g["near"] = [c for c in cands if c is not chosen
                         and abs(c.delta - chosen.delta) * sfC < 0.5 * hs.res_hz(0)]
            accepted[h.id] = g
            on_carbon.setdefault(chosen.id, []).append(g)
            continue
        g["reason"] = reason
    # Carbons closer than the f1 resolution: the HSQC only decides between them if the picks
    # do - the competing proton's picks must be at least half a digital point apart, in the same
    # order as the carbons; for a free (unassigned) neighbour, the chosen carbon must be nearer
    # to the picks by at least half a point.
    half_pt = 0.5 * hs.pt_hz(0) / sfC
    owner_of = {g["c"].id: g for gs in on_carbon.values() for g in gs}
    for gs in on_carbon.values():
        for g in gs:
            g["ambig"] = []
            for c2 in g.get("near", []):
                g2 = owner_of.get(c2.id)
                if g2 is not None and g2 is not g:
                    df, dc = g["f1"] - g2["f1"], g["c"].delta - c2.delta
                    if abs(df) < half_pt or df * dc < 0 or abs(dc) < half_pt:
                        g["ambig"].append(c2)
                elif g2 is None and abs(g["f1"] - c2.delta) - abs(g["f1"] - g["c"].delta) < half_pt:
                    g["ambig"].append(c2)
    # outputs for accepted groups
    for cid, gs in on_carbon.items():
        for g in gs:
            h, c = g["h"], g["c"]
            dH, dC = (g["f2"] - h.delta) * sfH, (g["f1"] - c.delta) * sfC
            h.carbons.append(c.id)
            c.protons.append(h.id)
            c.mult = "CH2" if g["kind"] == 2 else "CH/CH3"
            h.evidence.setdefault(HSQC, []).append(g["f2"])
            rep_pk = max(g["peaks"], key=lambda p: abs(p.intensity))
            matches.append((rep_pk, h, c, dH, dC))
            n_dup += len(g["peaks"]) - 1
            out.append(new_peak(c.delta, h.delta, g["intensity"], g["kind"],
                                st.mean(p.width1 for p in g["peaks"]), st.mean(p.width2 for p in g["peaks"]),
                                f"{h.id}/{c.id} snapped dH {dH:+.1f}Hz dC {dC:+.1f}Hz"
                                + (f", {len(g['peaks'])} picks merged" if len(g["peaks"]) > 1 else "")
                                + (" - ambiguous C: " + "/".join([c.id] + [x.id for x in g["ambig"]])
                                   + " (closer than HSQC f1 resolution) - review" if g["ambig"] else ""),
                                flags=rep_pk.flags_raw))
    for cid, gs in on_carbon.items():
        for g in gs:
            if g["ambig"]:
                others = ", ".join(f"{x.id} {x.delta:.3f}" for x in g["ambig"])
                rep.append(f"HSQC: {g['h'].id} {g['h'].delta:.4f} placed on {g['c'].id} {g['c'].delta:.3f}, but {others} "
                           f"is within {max(abs(x.delta - g['c'].delta) for x in g['ambig']) * sfC:.0f} Hz "
                           f"(f1 resolution {hs.res_hz(0):.0f} Hz) - assignment not determined by the HSQC, review")
                g["h"].flags.append("carbon ambiguous in HSQC: " + "/".join([g["c"].id] + [x.id for x in g["ambig"]]))
    if n_dup:
        rep.append(f"HSQC: {n_dup} multiplet-line picks merged into their cross-peaks"
                   + ("" if not P.hsqc_keep_duplicates else ""))
    # everything else
    for g in groups:
        if "c" in g:
            continue
        h = g["h"]
        cn = min(carbons, key=lambda c: abs(c.delta - g["f1"])) if carbons else None
        c_ok = cn is not None and abs(cn.delta - g["f1"]) <= tolC
        where = f"{g['f2']:.4f}/{g['f1']:.3f} ({len(g['peaks'])} pick(s))"
        reason = g.get("reason")
        strongest_row = max((x["strength"] for x in on_carbon.get(cn.id, [])), default=0) if c_ok else 0
        if reason is None and h is None and c_ok and on_carbon.get(cn.id):
            if g["strength"] < strongest_row:
                reason = f"no proton; {cn.id} already assigned ({g['strength'] / strongest_row:.0%} of its peak)"
        if reason is None and h is None and not c_ok:
            col = [x for x in on_carbon.values() for x in x if abs(x["f2"] - g["f2"]) * sfH <= max(10, hs.res_hz(1))]
            if col and max(x["strength"] for x in col) >= P.wiggle_ratio * g["strength"]:
                reason = "weak pick with no proton and no carbon - t1 noise?"
        if reason is not None:
            for p in g["peaks"]:
                out.append(new_peak(p.d1, p.d2, p.intensity, p.kind, p.width1, p.width2,
                                    "artifact: " + reason, type_=ARTIFACT, flags=p.flags_raw))
            rep.append(f"HSQC {where}: {reason} -> artifact")
        elif h is not None and not c_ok:
            h.evidence.setdefault(HSQC, []).append(g["f2"])
            out.append(new_peak(g["f1"], h.delta, g["intensity"], g["kind"], g["peaks"][0].width1,
                                g["peaks"][0].width2, f"{h.id}/C? f2 snapped, C from HSQC - review"))
            rep.append(f"HSQC {where}: no 13C peak within tolerance - f1 kept, review")
        elif h is None and c_ok:
            out.append(new_peak(cn.delta, g["f2"], g["intensity"], g["kind"], g["peaks"][0].width1,
                                g["peaks"][0].width2, f"H?/{cn.id} f1 snapped, H from HSQC - review"))
            rep.append(f"HSQC {where}: no proton resonance within tolerance - f2 kept, review")
        else:
            for p in g["peaks"]:
                out.append(new_peak(p.d1, p.d2, p.intensity, p.kind, p.width1, p.width2,
                                    "unmatched - review", flags=p.flags_raw))
            rep.append(f"HSQC {where}: unmatched in both dimensions - review")
    used = {id(p) for p in extra}
    for pk in hs.peaks:
        if pk.type != COMPOUND and id(pk) not in used:
            out.append(dict(pk.raw))
    return out, matches


# --------------------------------------------------------------------------- #
# HMBC
# --------------------------------------------------------------------------- #
def refine_hmbc(sp, protons, carbons, hsqc_matches, P: Params, rep: list):
    hm = sp.get(HMBC)
    if hm is None:
        return None, []
    sfC, sfH = hm.sf
    tolC = max(P.hmbc_f1_tol_pts * hm.pt_hz(0), P.hmbc_f1_tol_res * hm.res_hz(0)) / sfC
    rep.append(f"HMBC: f1 point {hm.pt_hz(0):.1f} Hz, f1 tolerance {tolC * sfC:.1f} Hz ({tolC:.3f} ppm), "
               f"f1 resolution {hm.res_hz(0):.0f} Hz; f2 point {hm.pt_hz(1):.2f} Hz, f2 resolution {hm.res_hz(1):.1f} Hz")

    # solvent recognition from the 1D lists
    h1, c13 = sp.get(H1_1D), sp.get(C13_1D)
    solvH = [p.d1 for p in h1.peaks if p.type == SOLVENT] if h1 else []
    solvH += solvent_h1_positions(hm.raw.get("solvent", ""))
    solvC = [p.d1 for p in c13.peaks if p.type == SOLVENT] if c13 else []
    if not solvC and h1 is not None:
        solvC = solvent_c13_positions(h1.raw.get("solvent", ""))
    h1_imp = [p for p in h1.peaks if p.type == IMPURITY and not p.has("weak")] if h1 else []

    res1, res2 = hm.res_hz(0), hm.res_hz(1)
    noise = {}                # proton id -> intensities of t1-noise picks in its column
    single_tol = max(P.hmbc_single_hz, P.hmbc_single_res * res2)
    presence_tol = max(P.hmbc_presence_hz, 1.5 * res2)
    bt_c_tol = max(tolC, res1 / sfC)
    mirror_tol = max(8.0, 1.5 * res2)
    pad_fp = (P.ps_range_pad_hz + res2 / 2) / sfH
    out, correlations = [], []
    rows = {}
    cand = [p for p in hm.peaks if p.type == COMPOUND or
            (P.reclassify_hmbc and p.type in (ARTIFACT, IMPURITY))]
    for pk in cand:
        if solvH and solvC and min(abs(x - pk.d2) for x in solvH) * sfH < 5 and \
                min(abs(x - pk.d1) for x in solvC) <= max(tolC, res1 / sfC):
            out.append(new_peak(pk.d1, pk.d2, pk.intensity, pk.kind, pk.width1, pk.width2,
                                "solvent (matched 1H/13C solvent peaks)", type_=SOLVENT, flags=pk.flags_raw))
            rep.append(f"HMBC {pk.d2:.4f}/{pk.d1:.3f}: compound-typed peak on solvent positions -> solvent")
            continue
        cands = sorted((c for c in carbons if abs(c.delta - pk.d1) <= tolC), key=lambda c: abs(c.delta - pk.d1))
        wide = False
        if not cands:
            # nothing within the tolerance: carbons within one f1 resolution width are still
            # candidates (low-resolution HMBC), decided after the wiggle test
            cands = sorted((c for c in carbons if abs(c.delta - pk.d1) * sfC <= res1),
                           key=lambda c: abs(c.delta - pk.d1))
            wide = bool(cands)
        ckey = (cands[0].id + ("~" if wide else "")) if cands else f"f1={pk.d1:.3f}"
        rows.setdefault(ckey, {"carbons": cands, "peaks": [], "wide": wide})["peaks"].append(pk)

    cmap = {c.id: c for c in carbons}
    cmap.update({c.id + "~": c for c in carbons})
    for ckey, row in rows.items():
        pks = sorted(row["peaks"], key=lambda p: p.d2)
        clusters = [[pks[0]]]
        for p in pks[1:]:
            if (p.d2 - clusters[-1][-1].d2) * sfH <= P.hmbc_gap_hz:
                clusters[-1].append(p)
            else:
                clusters.append([p])

        def centroid(cl):
            w = sum(abs(p.intensity) for p in cl)
            return sum(p.d2 * abs(p.intensity) for p in cl) / w

        def nearest(x):
            h = min(protons, key=lambda h: abs(h.delta - x))
            return h, (x - h.delta) * sfH

        # merge symmetric pairs about a resonance
        merged, i = [], 0
        while i < len(clusters):
            if i + 1 < len(clusters):
                a, b = clusters[i], clusters[i + 1]
                ca, cb = centroid(a), centroid(b)
                mid = (ca + cb) / 2
                h, dm = nearest(mid)
                _, da = nearest(ca)
                _, db = nearest(cb)
                if (cb - ca) * sfH <= P.hmbc_pair_max_hz and abs(dm) <= P.hmbc_sym_tol_hz \
                        and abs(da) > P.hmbc_sym_tol_hz and abs(db) > P.hmbc_sym_tol_hz:
                    merged.append((a + b, "pair", mid))
                    i += 2
                    continue
            merged.append((clusters[i], "lines", centroid(clusters[i])))
            i += 1

        # split clusters that are really two overlapping protons: if the whole cluster does not
        # centre on a resonance, but partitioning its lines by nearest resonance gives parts
        # that each centre on their own resonance, use the partition.
        split = []
        for cl, how, cen in merged:
            _, d_all = nearest(cen)
            parts = {}
            for p in cl:
                parts.setdefault(nearest(p.d2)[0].id, []).append(p)
            if len(parts) >= 2 and abs(d_all) > P.hmbc_sym_tol_hz and \
                    all(abs(nearest(centroid(v))[1]) <= P.hmbc_match_hz for v in parts.values()):
                for v in parts.values():
                    split.append((v, "split", centroid(v)))
                rep.append(f"HMBC row {ckey}: cluster {min(p.d2 for p in cl):.4f}-{max(p.d2 for p in cl):.4f} "
                           f"split between overlapping protons " + "/".join(parts))
            else:
                split.append((cl, how, cen))
        merged = split

        for cl, how, cen in merged:
            if not any(p.type == COMPOUND for p in cl):
                for p in cl:          # only Mnova non-compound lines: leave them as they were
                    out.append(dict(p.raw))
                continue
            n_re = sum(1 for p in cl if p.type != COMPOUND)
            f1 = sum(p.d1 * abs(p.intensity) for p in cl) / sum(abs(p.intensity) for p in cl)
            inten = sum(p.intensity for p in cl)
            h, dh = nearest(cen)
            footprint = False
            if abs(dh) > presence_tol:
                fp = proton_owner(protons, cen, pad_fp, 0.0)
                if fp is not None:
                    h, dh, footprint = fp, (cen - fp.delta) * sfH, True
            carbon = cmap.get(ckey)
            # 1J breakthrough: at H +/- 1J/2 in the row of that proton's own HSQC carbon.  Accepted when
            # no proton explains the position, or when the mirror line at H -/+ 1J/2 is also present.
            bt = None
            for _pk, hh, cc, *_ in hsqc_matches:
                if abs(cc.delta - f1) <= bt_c_tol:
                    o = (cen - hh.delta) * sfH
                    if P.onebond_min_hz / 2 <= abs(o) <= P.onebond_max_hz / 2:
                        mir = [q for q in cand if abs((q.d2 - hh.delta) * sfH + o) <= mirror_tol
                               and abs(q.d1 - cc.delta) <= bt_c_tol]
                        mirror = bool(mir)
                        if mirror and abs(dh) <= single_tol and \
                                sum(abs(p.intensity) for p in cl) >= 2 * max(abs(q.intensity) for q in mir):
                            # sits on a proton and is much stronger than its would-be mirror:
                            # a real correlation that happens to coincide with a breakthrough line
                            continue
                        if (abs(dh) > presence_tol and not footprint) or mirror:
                            bt = (hh, cc, o, mirror)
                            break
            if bt is not None:
                hh, cc, o, mirror = bt
                for p in cl:
                    out.append(new_peak(p.d1, p.d2, p.intensity, p.kind, p.width1, p.width2,
                                        f"1J breakthrough of {hh.id}/{cc.id} (1J~{2 * abs(o):.0f}Hz)"
                                        + (", mirror line present" if mirror else ""),
                                        type_=ARTIFACT, flags=p.flags_raw))
                rep.append(f"HMBC {cen:.4f}/{f1:.3f} ({len(cl)} line(s)): 1J breakthrough of "
                           f"{hh.id} {hh.delta:.4f}/{cc.id} {cc.delta:.3f}, offset {o:+.1f} Hz (1J~{2 * abs(o):.0f} Hz)"
                           + (", mirror line present" if mirror else ""))
                continue
            if carbon is not None and row["wide"]:
                imax = max(abs(p.intensity) for p in cl)
                strong = [q for q in cand if abs(q.d2 - cen) * sfH <= max(10.0, 1.5 * res2)
                          and abs(q.intensity) >= P.wiggle_ratio * imax
                          and 0 < abs(q.d1 - f1) * sfC <= P.wiggle_res * res1]
                if strong or abs(dh) > presence_tol and not footprint:
                    carbon = None          # fall through to the artefact / unassigned handling
            if carbon is None or (abs(dh) > presence_tol and not footprint):
                imax = max(abs(p.intensity) for p in cl)
                strong = [q for q in cand if abs(q.d2 - cen) * sfH <= max(10.0, 1.5 * res2)
                          and abs(q.intensity) >= P.wiggle_ratio * imax
                          and 0 < abs(q.d1 - f1) * sfC <= P.wiggle_res * res1]
                if strong:
                    q = max(strong, key=lambda q: abs(q.intensity))
                    why = (f"f1 truncation wiggle of {q.d2:.4f}/{q.d1:.2f} "
                           f"({imax / abs(q.intensity):.0%}, {(f1 - q.d1) * sfC:+.0f} Hz in f1)")
                elif carbon is None and abs(dh) <= presence_tol:
                    why = f"no carbon at f1 {f1:.2f} in the column of {h.id} - t1 noise?"
                    noise.setdefault(h.id, []).append(imax)
                else:
                    why = None
                if why:
                    for p in cl:
                        out.append(new_peak(p.d1, p.d2, p.intensity, p.kind, p.width1, p.width2,
                                            "artifact: " + why, type_=ARTIFACT, flags=p.flags_raw))
                    rep.append(f"HMBC {cen:.4f}/{f1:.3f} ({len(cl)} line(s)): {why} -> artifact")
                    continue
                if h1_imp and abs(dh) > presence_tol:
                    hit = min(h1_imp, key=lambda q: abs(q.d1 - cen))
                    if abs(hit.d1 - cen) * sfH <= max(P.hmbc_match_hz, single_tol):
                        why = f"on 1H {TYPE_NAMES[hit.type]} line {hit.d1:.4f}, no compound proton"
                        for p in cl:
                            out.append(new_peak(p.d1, p.d2, p.intensity, p.kind, p.width1, p.width2,
                                                why, type_=IMPURITY, flags=p.flags_raw))
                        rep.append(f"HMBC {cen:.4f}/{f1:.3f} ({len(cl)} line(s)): {why} -> impurity")
                        continue
                for p in cl:
                    out.append(new_peak(p.d1, p.d2, p.intensity, p.kind, p.width1, p.width2,
                                        "unassigned - review", flags=p.flags_raw))
                rep.append(f"HMBC {cen:.4f}/{f1:.3f} ({len(cl)} line(s)): unassigned "
                           f"(nearest proton {h.id} {dh:+.1f} Hz, carbon {'none' if carbon is None else carbon.id})")
                continue
            quality = "position" if (len(cl) >= 2 and abs(dh) <= max(P.hmbc_match_hz, single_tol)) else \
                      ("position-1line" if abs(dh) <= single_tol else "presence")
            if how == "pair":
                quality = "position"
            amb = [c.id for c in row["carbons"][1:]]
            if row["wide"]:
                quality = "presence"
            correlations.append(dict(h=h, c=carbon, n=len(cl), n_reclassified=n_re, centroid=cen, dH=dh, f1=f1,
                                     wide=row["wide"], cands=row["carbons"],
                                     dC=(f1 - carbon.delta) * sfC, intensity=inten, quality=quality,
                                     ambiguous=amb, how=how,
                                     w1=st.mean(p.width1 for p in cl), w2=st.mean(p.width2 for p in cl)))
            if quality == "position":
                h.evidence.setdefault(HMBC, []).append(cen)

    # combine duplicates of the same (proton, carbon)
    comb = {}
    taken = {(c["h"].id, c["c"].id) for c in correlations if not c.get("wide")}
    for c in correlations:
        if c.get("wide") and (c["h"].id, c["c"].id) in taken:
            # this proton already correlates with that carbon within tolerance: the off-tolerance
            # pick is more likely its other candidate carbon
            alt = next((x for x in c["cands"] if (c["h"].id, x.id) not in taken), None)
            if alt is not None:
                c["ambiguous"] = [x.id for x in c["cands"] if x is not alt]
                c["c"], c["dC"] = alt, (c["f1"] - alt.delta) * sfC
        k = (c["h"].id, c["c"].id)
        if k in comb:
            o = comb[k]
            o["n"] += c["n"]; o["intensity"] += c["intensity"]; o["n_reclassified"] += c["n_reclassified"]
            rank = {"position": 0, "position-1line": 1, "presence": 2}
            if rank[c["quality"]] < rank[o["quality"]]:
                o["quality"], o["dH"] = c["quality"], c["dH"]
        else:
            comb[k] = c
    strongest = {}
    for c in comb.values():
        strongest[c["h"].id] = max(strongest.get(c["h"].id, 0), abs(c["intensity"]))
    for c in comb.values():
        ann = f"{c['h'].id}->{c['c'].id} HMBC {c['n']} line(s) {c['quality']} dH {c['dH']:+.1f}Hz"
        rel = abs(c["intensity"]) / strongest[c["h"].id]
        if rel < P.hmbc_weak_frac:
            c["weak"] = rel
            ann += f" - weak ({rel:.0%} of strongest for {c['h'].id}): 2J/4J?"
        nz = noise.get(c["h"].id, [])
        if len(nz) >= 2 and abs(c["intensity"]) <= 1.5 * max(nz):
            c["noise_level"] = True
            ann += f" - weak, comparable to t1 noise in this column ({len(nz)} noise picks) - review"
        if c["n_reclassified"]:
            ann += f" (incl. {c['n_reclassified']} Mnova non-compound line(s))"
        if c.get("wide"):
            ann += (f" - f1 {c['dC']:+.0f} Hz from {c['c'].id}, beyond the tolerance but within the f1 "
                    f"resolution ({res1:.0f} Hz) - review")
        if c["ambiguous"]:
            ann += " ambiguous C: " + "/".join([c["c"].id] + c["ambiguous"])
        out.append(new_peak(c["c"].delta, c["h"].delta, c["intensity"], 1, c["w1"], c["w2"], ann))
    used = {id(p) for p in cand}
    for pk in hm.peaks:
        if pk.type != COMPOUND and id(pk) not in used:
            out.append(dict(pk.raw))
    return out, list(comb.values())


# --------------------------------------------------------------------------- #
# COSY (clustering / diagnostics only)
# --------------------------------------------------------------------------- #
def cosy_spin_systems(sp, protons, P: Params, rep: list):
    cs = sp.get(COSY)
    if cs is None or not cs.compound():
        rep.append("COSY: not available - spin systems not built")
        return []
    sf = cs.sf[-1]
    t2 = P.cosy_f2_tol_hz / sf
    t1 = max(P.cosy_f1_tol_pts * cs.res_hz(0), P.cosy_f2_tol_hz) / sf
    parent = {h.id: h.id for h in protons}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    links = {}
    for pk in cs.compound():
        a = min(protons, key=lambda h: abs(h.delta - pk.d2))
        b = min(protons, key=lambda h: abs(h.delta - pk.d1))
        if abs(a.delta - pk.d2) > t2 or abs(b.delta - pk.d1) > t1 or a.id == b.id:
            continue
        k = tuple(sorted((a.id, b.id)))
        links[k] = links.get(k, 0) + 1
        parent[find(a.id)] = find(b.id)
    systems = {}
    for h in protons:
        systems.setdefault(find(h.id), []).append(h)
    out = [sorted(v, key=lambda h: -h.delta) for v in systems.values()]
    out.sort(key=lambda v: -v[0].delta)
    hmap = {h.id: h for h in protons}
    for i, v in enumerate(out):
        for h in v:
            h.spin_system = i + 1
    rep.append(f"COSY: f2 tol {P.cosy_f2_tol_hz:.0f} Hz, f1 tol {t1 * sf:.0f} Hz; {len(links)} proton pairs linked")
    for i, v in enumerate(out):
        rep.append(f"  spin system {i + 1}: " + ", ".join(f"{h.id} {h.delta:.4f}" for h in v))
    for (a, b), n in sorted(links.items()):
        dv = abs(hmap[a].delta - hmap[b].delta) * sf
        if dv < P.strong_coupling_hz:
            msg = f"COSY-linked {a}/{b} only {dv:.1f} Hz apart: possible strong coupling (or near-diagonal artefact)"
            hmap[a].flags.append(msg); hmap[b].flags.append(msg)
            rep.append("  " + msg)
    return out


# --------------------------------------------------------------------------- #
# Referencing check
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# COSY and other H-C 2D spectra (CLIP-COSY, DEPT-HSQC ...): snap to the shift lists
# --------------------------------------------------------------------------- #
def proton_owner(protons, x, pad, tol):
    """Proton a line at x belongs to: inside its 1H multiplet footprint (padded), else nearest within tol."""
    inside = [h for h in protons if h.range[1] > h.range[0] and h.range[0] - pad <= x <= h.range[1] + pad]
    if inside:
        return min(inside, key=lambda h: abs(h.delta - x))
    h = min(protons, key=lambda h: abs(h.delta - x)) if protons else None
    return h if h is not None and abs(h.delta - x) <= tol else None


def refine_cosy(sp, protons, P: Params, rep: list, confirm=None):
    """One peak per coupled proton pair at (f2 = higher shift, f1 = lower shift).

    Protons closer than the COSY f1 resolution form an unresolved group (their f1 positions
    cannot be told apart, and their f2 multiplet lines interleave).  Picks are first collected
    per (group, group) cross-peak; the member of a multi-proton group is then chosen once per
    cross-peak from the intensity-weighted centroid of all its picks - the f2 centroid (high
    resolution) when the group appears in f2, otherwise the f1 centroid - and annotated."""
    cs = sp.get(COSY)
    if cs is None:
        return None, {}
    sf = cs.sf[-1]
    pad2 = max(P.ps_range_pad_hz, 2 * cs.pt_hz(1)) / sf
    pad1 = max(P.ps_range_pad_hz, cs.res_hz(0)) / sf
    t2 = P.cosy_f2_tol_hz / sf
    t1 = max(P.cosy_f1_tol_pts * cs.res_hz(0), P.cosy_f2_tol_hz) / sf
    res1 = max(cs.res_hz(0), P.cosy_f2_tol_hz) / sf
    # unresolved groups
    gid = {h.id: h.id for h in protons}

    def find(x):
        while gid[x] != x:
            x = gid[x]
        return x
    for i, a in enumerate(protons):
        for b in protons[i + 1:]:
            if abs(a.delta - b.delta) < res1:
                gid[find(a.id)] = find(b.id)
    groups = {}
    for h in protons:
        groups.setdefault(find(h.id), []).append(h)
    for g in groups.values():
        if len(g) > 1:
            rep.append("COSY: protons " + "/".join(f"{h.id} {h.delta:.4f}" for h in g) +
                       f" closer than f1 resolution ({cs.res_hz(0):.0f} Hz) - treated as one group, "
                       f"couplings within the group cannot be read from the COSY")

    def wc(vals, w):
        return _wmean(vals, [abs(x) for x in w])

    out, cross, n_diag, n_un, n_imp = [], {}, 0, 0, 0
    h1 = sp.get(H1_1D)
    h1_imp = [p for p in h1.peaks if p.type == IMPURITY and not p.has("weak")] if h1 else []
    for pk in cs.compound():
        a = proton_owner(protons, pk.d2, pad2, t2)
        b = proton_owner(protons, pk.d1, pad1, t1)
        if a is None or b is None:
            if abs(pk.d2 - pk.d1) * sf <= max(cs.res_hz(0), P.cosy_f2_tol_hz):
                n_diag += 1
                sat = [h for h in protons
                       if P.onebond_min_hz / 2 <= abs(pk.d2 - h.delta) * sf <= P.onebond_max_hz / 2]
                why = (f"13C satellite of {sat[0].id}? (1J~{2 * abs(pk.d2 - sat[0].delta) * sf:.0f} Hz)"
                       if sat else "no compound proton")
                out.append(new_peak(pk.d1, pk.d2, pk.intensity, pk.kind, pk.width1, pk.width2,
                                    f"diagonal, {why}", type_=ARTIFACT, flags=pk.flags_raw))
            elif a is None and b is None:
                n_imp += 1
                out.append(new_peak(pk.d1, pk.d2, pk.intensity, pk.kind, pk.width1, pk.width2,
                                    "no compound proton in either dimension - impurity?", type_=IMPURITY,
                                    flags=pk.flags_raw))
            else:
                free = pk.d2 if a is None else pk.d1      # the dimension with no compound proton
                tolf = t2 if a is None else t1
                hit = min(h1_imp, key=lambda q: abs(q.d1 - free)) if h1_imp else None
                if hit is not None and abs(hit.d1 - free) <= max(tolf, P.cosy_f2_tol_hz / sf):
                    out.append(new_peak(pk.d1, pk.d2, pk.intensity, pk.kind, pk.width1, pk.width2,
                                        f"on 1H {TYPE_NAMES[hit.type]} line {hit.d1:.4f}, no compound proton",
                                        type_=IMPURITY, flags=pk.flags_raw))
                    n_imp += 1
                    continue
                n_un += 1
                out.append(new_peak(pk.d1, pk.d2, pk.intensity, pk.kind, pk.width1, pk.width2,
                                    "unassigned - review", flags=pk.flags_raw))
            continue
        ga, gb = find(a.id), find(b.id)
        if ga == gb:
            n_diag += 1
            continue
        key = tuple(sorted((ga, gb)))
        cross.setdefault(key, []).append((pk, ga, gb))

    def choose(group_key, picks):
        """members of a group taking part in one cross-peak: [(proton, note)]"""
        members = groups[group_key]
        if len(members) == 1:
            return [(members[0], "")]
        as_f2 = [pk for pk, ga, gb in picks if ga == group_key]
        if as_f2:
            # The f2 lines of one multiplet spread over its whole width, so giving each line to
            # its nearest member splits them evenly between close members.  Use the
            # intensity-weighted f2 centroid of the cross-peak instead; keep a second member
            # only when the centroid lies near the middle between the two (genuinely shared).
            c = wc([p.d2 for p in as_f2], [p.intensity for p in as_f2])
            ranked = sorted(members, key=lambda h: abs(h.delta - c))
            best, nxt = ranked[0], ranked[1]
            gap = abs(nxt.delta - best.delta) or 1e-9
            kept = [best]
            between = (c - best.delta) * (nxt.delta - best.delta) > 0
            # one line is a single multiplet line, not a centre: never enough to keep two members
            if len(as_f2) > 1 and between and abs(c - best.delta) / gap > P.cosy_group_share:
                kept.append(nxt)
            note = f"f2 centroid {c:.4f} -> " + "/".join(h.id for h in kept)
            return [(h, note) for h in kept]
        as_f1 = [pk for pk, ga, gb in picks if gb == group_key]
        c = wc([p.d1 for p in as_f1], [p.intensity for p in as_f1])
        h = min(members, key=lambda h: abs(h.delta - c))
        return [(h, f"{h.id} chosen over " + "/".join(x.id for x in members if x is not h) + " by f1 centroid")]

    links = {}
    inten = {k: sum(abs(p.intensity) for p, _, _ in v) for k, v in cross.items()}
    strongest = {}
    for (g1, g2), v in inten.items():
        for g in (g1, g2):
            strongest[g] = max(strongest.get(g, 0), v)
    group_links = {}
    n_pairs = 0
    for (g1, g2), picks in cross.items():
        for x in groups[g1]:
            for y in groups[g2]:
                group_links.setdefault(x.id, set()).add(y.id)
                group_links.setdefault(y.id, set()).add(x.id)
        c1, c2 = choose(g1, picks), choose(g2, picks)
        combos = [(a, na, b, nb) for a, na in c1 for b, nb in c2]
        multi = len(combos) > 1
        how = ""
        if confirm is not None and (len(groups[g1]) > 1 or len(groups[g2]) > 1):
            ok = [c for c in combos if frozenset((c[0].id, c[2].id)) in confirm]
            if ok:
                combos, how = ok, "confirmed by CLIP-COSY relay"
                multi = False
        pks = [p for p, _, _ in picks]
        total = sum(abs(p.intensity) for p in pks)
        for a, na, b, nb in combos:
            hi, lo = (a, b) if a.delta >= b.delta else (b, a)
            notes = [n for n in (na, nb, how) if n]
            ann = f"{hi.id}-{lo.id} COSY {len(pks)} line(s)" + "".join(f"; {n}" for n in dict.fromkeys(notes))
            if multi:
                ann += "; unresolved group, coupling may belong to only one member - review"
            dv = abs(hi.delta - lo.delta) * sf
            if dv < P.strong_coupling_hz:
                ann += f"; near diagonal ({dv:.1f} Hz) - review"
            rel = inten[(g1, g2)] / max(strongest[g1], strongest[g2])
            if rel < P.cosy_weak_frac:
                ann += f"; weak ({rel:.0%} of strongest coupling of these protons) - long-range?"
            out.append(new_peak(lo.delta, hi.delta, total / len(combos), 1,
                                st.mean(p.width1 for p in pks), st.mean(p.width2 for p in pks), ann))
            links.setdefault(hi.id, set()).add(lo.id)
            links.setdefault(lo.id, set()).add(hi.id)
            n_pairs += 1
    for pk in cs.peaks:
        if pk.type != COMPOUND:
            out.append(dict(pk.raw))
    rep.append(f"COSY output: {len(cs.compound())} compound picks -> {n_pairs} coupled pairs; "
               f"{n_diag} diagonal picks dropped" + (f"; {n_imp} impurity (no compound proton in one or both dimensions)" if n_imp else "")
               + (f"; {n_un} unassigned" if n_un else ""))
    return out, group_links


def _is_clipcosy(spec):
    """CLIP-COSY is the dialog label HSQC_CLIPCOSY (the label is authoritative)."""
    return spec.name == HSQC_CLIPCOSY


def refine_hc_other(spec, protons, carbons, links, P: Params, rep: list, cosy_present=True, hsqc=None):
    """Snap an additional H-C correlation spectrum (CLIP-COSY, DEPT-edited HSQC ...).

    Each compound pick is given to a proton (1H footprint) and a carbon.  Expected positions
    get a wider f1 tolerance: direct = the proton's own HSQC carbon; relay (CLIP-COSY only) =
    carbons carrying a proton COSY-coupled to this one, or - when no COSY was supplied -
    any protonated carbon.  Picks of the same (proton, carbon, sign) are merged.  Artefacts:
      - weak opposite-sign picks at an assigned (proton, carbon): truncation wiggles,
      - side lines in the row of an assigned peak (within 1J/2 in f2, no own proton),
      - picks at an f1 with no carbon, much weaker than the strongest pick in their column."""
    sfC, sfH = spec.sf
    clip = _is_clipcosy(spec)
    pad = max(P.ps_range_pad_hz, 2 * spec.pt_hz(1)) / sfH
    tolH = proton_tol_ppm(spec, hsqc)
    strict = max(P.hc_tol_pts * spec.pt_hz(0), 0.5 * spec.res_hz(0)) / sfC
    wide = max(strict, P.hc_expected_res * spec.res_hz(0) / sfC)
    cmap, hmap = {c.id: c for c in carbons}, {h.id: h for h in protons}
    on_c = {c.id: [h.id for h in protons if c.id in h.carbons] for c in carbons}
    groups, leftovers, best_partner = {}, [], {}
    for pk in spec.compound():
        h = proton_owner(protons, pk.d2, pad, tolH)
        if h is None:
            leftovers.append((pk, None))
            continue
        expected = {cid: "direct" for cid in h.carbons}
        partner_of = {}
        if clip:
            if cosy_present:
                rel = {cid: oid for oid in links.get(h.id, ()) for cid in hmap[oid].carbons}
            else:
                rel = {cid: "/".join(on_c[cid]) for cid in on_c if on_c[cid] and cid not in h.carbons}
            for cid, oid in rel.items():
                if cid not in expected:
                    expected[cid] = f"relay {h.id}->{cid} (coupled to {oid})"
                    partner_of[cid] = oid
        best = None
        tol = wide if cosy_present else strict
        for cid, how in expected.items():
            d = abs(cmap[cid].delta - pk.d1)
            if d <= (wide if how == "direct" else tol) and (best is None or d < best[2]):
                best = (cmap[cid], how, d)
        if best is None and carbons:
            c = min(carbons, key=lambda c: abs(c.delta - pk.d1))
            if abs(c.delta - pk.d1) <= strict and c.protons:
                best = (c, "unexpected (not direct, no COSY link) - review", 0)
        if best is None:
            leftovers.append((pk, h))
            continue
        c, how, _ = best
        if c.id in partner_of and how.startswith("relay"):
            best_partner[id(pk)] = partner_of[c.id]
        sign = -1 if pk.intensity < 0 else 1
        groups.setdefault((h.id, c.id, sign), [h, c, how, [], best_partner.get(id(pk))])[3].append(pk)

    out, n_art = [], 0

    def artifact(pks, why):
        nonlocal n_art
        for p in pks:
            out.append(new_peak(p.d1, p.d2, p.intensity, p.kind, p.width1, p.width2,
                                "artifact: " + why, type_=ARTIFACT, flags=p.flags_raw))
        n_art += len(pks)

    # opposite-sign wiggles at the same (proton, carbon)
    strength = {k: sum(abs(p.intensity) for p in v[3]) for k, v in groups.items()}
    for (hid, cid, sign) in list(groups):
        other = (hid, cid, -sign)
        if other in groups and strength[(hid, cid, sign)] < 0.25 * strength[other]:
            h, c, how, pks, _ = groups.pop((hid, cid, sign))
            artifact(pks, f"opposite-sign wiggle of {hid}/{cid} "
                          f"({strength[(hid, cid, sign)] / strength[other]:.0%})")
    colw = max(10.0, 1.5 * spec.res_hz(1)) / sfH
    allpk = spec.compound()
    for pk, h in leftovers:
        c = min(carbons, key=lambda c: abs(c.delta - pk.d1)) if carbons else None
        c_ok = c is not None and abs(c.delta - pk.d1) <= strict
        side = None
        if c_ok:
            for (hid, cid, sign), g in groups.items():
                if cid == c.id and 0 < abs(pk.d2 - g[0].delta) * sfH <= P.onebond_max_hz / 2:
                    side = g
                    break
        cm = max((abs(q.intensity) for q in allpk if abs(q.d2 - pk.d2) <= colw), default=0)
        if side is not None:
            artifact([pk], f"side line in the row of {side[0].id}/{side[1].id} "
                           f"({(pk.d2 - side[0].delta) * sfH:+.0f} Hz)")
        elif not c_ok and cm >= P.wiggle_ratio * abs(pk.intensity):
            artifact([pk], f"no carbon at f1 {pk.d1:.2f}, {abs(pk.intensity) / cm:.0%} of the strongest pick "
                           f"in its column - wiggle/t1 noise?")
        elif h is not None:
            out.append(new_peak(pk.d1, h.delta, pk.intensity, pk.kind, pk.width1, pk.width2,
                                f"{h.id}/C? f2 snapped, no carbon - review", flags=pk.flags_raw))
        else:
            out.append(new_peak(pk.d1, pk.d2, pk.intensity, pk.kind, pk.width1, pk.width2,
                                "unassigned - review", flags=pk.flags_raw))
    n_un = sum(1 for p in out if p["type"] == COMPOUND and "review" in p["annotation"])
    kept = dropped = 0
    relay_pairs = set()
    for (hid, cid, sign), (h, c, how, pks, partner) in groups.items():
        if partner:
            for pid in partner.split("/"):
                relay_pairs.add(frozenset((h.id, pid)))
        if clip and how == "direct" and not P.clipcosy_keep_direct:
            dropped += 1
            continue
        kept += 1
        inten = sum(p.intensity for p in pks)
        out.append(new_peak(c.delta, h.delta, inten, 2 if inten < 0 else 1,
                            st.mean(p.width1 for p in pks), st.mean(p.width2 for p in pks),
                            f"{h.id}/{c.id} {how}, {len(pks)} pick(s)"))
    for pk in spec.peaks:
        if pk.type != COMPOUND:
            out.append(dict(pk.raw))
    kind = "CLIP-COSY" if clip else "H-C 2D"
    rep.append(f"{spec.name} ({kind}, {spec.raw.get('pulsesequence', '')}): {len(spec.compound())} compound picks -> "
               f"{kept} peaks" + (f"; {dropped} direct peaks dropped (in HSQC)" if dropped else "")
               + (f"; {n_art} artifact picks" if n_art else "")
               + (f"; {n_un} for review" if n_un else "")
               + f"; f1 tolerance {strict * sfC:.0f} Hz"
               + (f" (expected positions {wide * sfC:.0f} Hz)" if cosy_present else
                  " (no COSY: relays matched to any protonated carbon)" if clip else ""))
    signs = {}
    for (hid, cid, sign), (h, c, how, pks, partner) in groups.items():
        signs.setdefault("direct" if how == "direct" else "relay" if how.startswith("relay") else "other", set()).add(sign)
    if clip and len(signs.get("direct", set())) == 1 and signs.get("relay") == signs.get("direct"):
        rep.append(f"  note: direct and relay peaks have the same sign in {spec.name}")
    return out, relay_pairs


def referencing(hsqc_matches, hmbc_corr, rep: list):
    def summarise(label, dHs, dCs, Cs):
        if not dHs:
            return
        line = f"{label}: median dH {st.median(dHs):+.2f} Hz"
        if dCs:
            line += f", median dC {st.median(dCs):+.1f} Hz"
            if len(Cs) >= 4 and (max(Cs) - min(Cs)) > 20:
                slope, icpt = np.polyfit(Cs, dCs, 1)
                line += f", dC trend {slope * 100:+.1f} Hz per 100 ppm"
        rep.append(line + f" (n={len(dHs)})")
    summarise("HSQC vs 1D", [m[3] for m in hsqc_matches], [m[4] for m in hsqc_matches],
              [m[2].delta for m in hsqc_matches])
    pos = [c for c in hmbc_corr if c["quality"] == "position"]
    summarise("HMBC vs 1D (multi-line correlations)", [c["dH"] for c in pos], [c["dC"] for c in pos],
              [c["c"].delta for c in pos])


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def refine(path_in: str, path_out: str | None = None, params: Params | None = None):
    P = params or Params()
    doc, sp = load(path_in)
    rep = [f"simpleNMRclaudePeakRefine v{VERSION}: {path_in}"]
    if getattr(load, "skipped", None):
        rep.append("Skipped (not analysed, left unchanged): " + "; ".join(load.skipped))
    for w in getattr(load, "warnings", []):
        rep.append(("Label " if w.split(": ", 1)[-1].startswith("note") else "Label warning: ") + w)
    missing = [n for n in (H1_PURESHIFT, C13_1D) if n not in sp]
    if missing:
        rep.append("Not supplied: " + ", ".join(missing) + " - " +
                   ("proton shifts from 1H multiplets, split by HSQC rows where needed; " if H1_PURESHIFT in missing else "") +
                   ("carbon shifts from HSQC/HMBC f1" if C13_1D in missing else ""))
    for s in sp.values():
        if s.kind_fixed:
            rep.append(f"  {s.name}: {s.kind_fixed} picks had 'kind' inconsistent with the sign of the intensity "
                       f"- sign used")
        rep.append(f"  {s.name:10s} {s.raw.get('pulsesequence', ''):24s} peaks {len(s.peaks):4d} "
                   f"(compound {len(s.compound())}) res {', '.join(f'{s.res_hz(i):.1f}' for i in range(s.ndim))} Hz")

    # Pass 1 may include pureshift peaks Mnova typed solvent/impurity that sit inside compound
    # 1H multiplets.  They are kept only if HSQC or HMBC confirm them; otherwise they are
    # excluded and the analysis is repeated.
    exclude, head = set(), list(rep)
    strong_pos = None
    excl_c = {}
    for _ in range(8):
        rep = list(head)
        protons = build_protons(sp, P, rep, frozenset(exclude), strong_pos)
        pureshift_sidebands(sp, protons, P, rep)
        carbons = build_carbons(sp, P, rep, protons, excl_c)
        hsqc_out, hsqc_matches = refine_hsqc(sp, protons, carbons, P, rep)
        hmbc_out, hmbc_corr = refine_hmbc(sp, protons, carbons, hsqc_matches, P, rep)
        new_c = {k: v for k, v in impurity_carbons(sp, protons, carbons, hmbc_corr, P).items() if k not in excl_c}
        if new_c:
            excl_c.update(new_c)
            c13 = sp[C13_1D]
            head.append("13C: " + "; ".join(f"{next(p.d1 for p in c13.peaks if p.key == k):.3f} treated as impurity - {v}"
                                            for k, v in new_c.items()))
            continue
        # surplus pureshift peaks: a Mnova multiplet whose proton count is already filled by protons
        # with HSQC carbons cannot hold further carbon-less pureshift peaks (sidebands/artefacts)
        h1 = sp[H1_1D]
        surplus = []
        groups = {}
        for h in protons:
            if h.mult_nh and h.source == "pureshift":
                groups.setdefault((h.mult_delta, tuple(h.range)), []).append(h)
        sfH_ = h1.sf[-1]
        for hs_ in groups.values():
            n_mult = round(hs_[0].mult_nh / (h1.nh_scale or 1))
            with_c = [h for h in hs_ if h.carbons]
            extra_ = [h for h in hs_ if not h.carbons and h.ps_key is not None]
            if extra_ and with_c and len(with_c) >= n_mult:
                if n_mult == 1 and len(with_c) == 1:
                    # one proton, several pureshift peaks: keep the peak that best matches the
                    # independent evidence (Mnova multiplet centre, HSQC row, HMBC columns)
                    hc_ = with_c[0]
                    ev = [hc_.mult_delta] + hc_.evidence.get(HSQC, []) + \
                         [c["h"].delta + c["dH"] / sfH_ for c in hmbc_corr if c["h"] in hs_]
                    target = st.median(ev)
                    best = min([hc_] + extra_, key=lambda h: abs(h.delta - target))
                    if best is not hc_ and abs(hc_.delta - target) - abs(best.delta - target) > 1.0 / sfH_:
                        surplus += [hc_] + [h for h in extra_ if h is not best]
                        continue
                surplus += extra_
        if surplus:
            exclude |= {h.ps_key for h in surplus}
            head.append("Protons: surplus pureshift peaks in Mnova multiplets already filled by protons with HSQC "
                        "carbons (the peak best matching the multiplet centre, HSQC row and HMBC columns is kept) - "
                        "treated as pureshift artefacts: " + ", ".join(
                            f"{h.delta:.4f} (multiplet {h.mult_delta:.4f}, {round(h.mult_nh)}H)" for h in surplus))
            continue
        rejected = []
        for h in protons:
            if h.mnova_type == COMPOUND:
                continue
            n_hmbc = sum(1 for c in hmbc_corr if c["h"] is h and c["quality"] == "position")
            if h.carbons or n_hmbc:
                msg = (f"pureshift {h.delta:.4f} typed {TYPE_NAMES[h.mnova_type]} by Mnova reclassified as "
                       f"compound (HSQC partner: {','.join(h.carbons) or 'none'}; HMBC correlations: {n_hmbc})")
                h.flags.append(f"Mnova typed {TYPE_NAMES[h.mnova_type]}; reclassified")
                rep.append("Classification: " + msg)
            else:
                rejected.append(h)
        if not rejected:
            if strong_pos is None and any(h.strong and h.split_c is not None for h in protons):
                strong_pos = strong_positions(protons, hmbc_corr, sp[H1_1D].sf[-1])   # rerun with HSQC/HMBC positions
                if strong_pos:
                    continue
            break
        exclude |= {h.ps_key for h in rejected}
        head.append("Classification: pureshift " + ", ".join(f"{h.delta:.4f}" for h in rejected) +
                    " typed non-compound by Mnova and not confirmed by HSQC/HMBC - kept as non-compound")
    systems = cosy_spin_systems(sp, protons, P, rep)
    # pass 1: COSY links (unresolved groups linked as a whole) to set CLIP-COSY relay expectations;
    # pass 2: COSY output, with couplings inside unresolved groups confirmed by CLIP-COSY relays
    _, links = refine_cosy(sp, protons, P, [])
    # a COSY with almost no cross-peaks picked cannot say which relays to expect
    n_cosy = len(sp[COSY].compound()) if COSY in sp else 0
    cosy_usable = COSY in sp and n_cosy >= max(3, len(protons) // 3)
    if COSY in sp and not cosy_usable:
        rep.append(f"COSY: only {n_cosy} compound picks - not used to set CLIP-COSY/relay expectations")
    extra, relay_pairs = {}, set()
    for name, spec in sp.items():
        if spec.ndim == 2 and name not in (HSQC, HMBC, COSY) and \
                [str(n).strip() for n in spec.nuclei] == ["13C", "1H"]:
            extra[name], rp = refine_hc_other(spec, protons, carbons, links, P, rep,
                                              cosy_present=cosy_usable, hsqc=sp.get(HSQC))
            if _is_clipcosy(spec):
                relay_pairs |= rp
    cosy_out, _ = refine_cosy(sp, protons, P, rep, confirm=relay_pairs or None)
    if cosy_out is not None:
        extra[COSY] = cosy_out
    referencing(hsqc_matches, hmbc_corr, rep)

    sfH = sp[H1_1D].sf[-1]
    # per-proton evidence summary
    tbl = ["", "Protons (delta used for snapping; other estimates as offsets in Hz)",
           f"{'id':4s} {'delta':>8s} {'source':12s} {'H1 mult':10s} {'nH':>4s} {'H1cen':>6s} "
           f"{HSQC:>6s} {'HMBC(n)':>14s} {'C (HSQC)':10s} flags"]
    for h in protons:
        def off(v):
            return f"{(v - h.delta) * sfH:+.1f}" if v is not None else "-"
        hs = h.evidence.get(HSQC, [])
        hb = h.evidence.get(HMBC, [])
        hb_s = f"{st.mean([(x - h.delta) * sfH for x in hb]):+.1f}({len(hb)})" if hb else "-"
        tbl.append(f"{h.id:4s} {h.delta:8.4f} {h.source:12s} {h.mult_category:3s}{(h.mult_delta or 0):7.4f} "
                   f"{(f'{h.nH_est:4.1f}' if h.nH_est else '   -')} {off(h.h1_centroid):>6s} {off(hs[0]) if hs else '-':>6s} {hb_s:>14s} "
                   f"{','.join(h.carbons) or '-':10s} {'; '.join(h.flags)}")
    tbl += ["", "Carbons", f"{'id':4s} {'delta':>9s} {HSQC:8s} {'source':6s} protons"]
    for c in carbons:
        tbl.append(f"{c.id:4s} {c.delta:9.4f} {c.mult or 'no HSQC':8s} {c.source:6s} {','.join(c.protons) or '(quaternary?)'}")
    tbl += ["", "HMBC correlations", f"{'H':4s} {'C':4s} {'n':>2s} {'dH Hz':>6s} {'dC Hz':>6s} quality"]
    for c in sorted(hmbc_corr, key=lambda c: (-c["h"].delta, -c["c"].delta)):
        tbl.append(f"{c['h'].id:4s} {c['c'].id:4s} {c['n']:2d} {c['dH']:+6.1f} {c['dC']:+6.1f} {c['quality']}"
                   + ("  wide-f1" if c.get("wide") else "")
                   + (f"  ambiguous with {'/'.join(c['ambiguous'])}" if c["ambiguous"] else "")
                   + ("  t1-noise level" if c.get("noise_level") else "")
                   + (f"  weak {c['weak']:.0%}" if "weak" in c else ""))
    rep += tbl

    # ---- write-back JSON -----------------------------------------------------
    new = copy.deepcopy(doc)
    originals = {}
    for name, peaks in [(HSQC, hsqc_out), (HMBC, hmbc_out)] + list(extra.items()):
        if peaks is None:
            continue
        key = sp[name].key
        originals[name] = doc[key]["peaks"]
        new[key]["peaks"] = {"datatype": "peaks", "count": len(peaks),
                             "data": {str(i): p for i, p in enumerate(peaks)}}
    originals.update(write_1d(new, doc, sp, protons, carbons, rep))
    fcheck = formula_checks(molecular_formula(doc), protons, carbons, new, sp, extra, rep)
    # safety net: compound 2D peaks must lie on the proton/carbon lists (what simpleNMR checks);
    # anything left over (e.g. "unassigned - review") is written as artifact, the note is kept
    hset = [h.delta for h in protons]
    cset = [c.delta for c in carbons]
    n_off = 0
    for name in [n for n, pk in [(HSQC, hsqc_out), (HMBC, hmbc_out)] + list(extra.items()) if pk is not None]:
        hh = sp[name].nuclei[0] == "1H"
        for p in new[sp[name].key]["peaks"]["data"].values():
            if p["type"] != COMPOUND:
                continue
            on_h = any(abs(p["delta2"] - x) < 1e-9 for x in hset)
            on_1 = any(abs(p["delta1"] - x) < 1e-9 for x in (hset if hh else cset))
            if not (on_h and on_1):
                p["type"] = ARTIFACT
                p["annotation"] = (p.get("annotation") or "") + " (not on the proton/carbon lists - typed artifact)"
                n_off += 1
    if n_off:
        rep.append(f"Write-back: {n_off} unresolved 2D picks typed artifact so that every compound peak lies on the "
                   f"proton/carbon lists (see their annotations)")
    new["peakRefinement"] = {"datatype": "peakRefinement", "count": 1, "data": {"0": {
        "version": VERSION,
        "parameters": P.__dict__,
        "tolerance": {"version": tolerance.TOLERANCE_VERSION, "n_points": tolerance.N_POINTS,
                      "floor_hz": tolerance.FLOOR_HZ},
        "molecularFormula": fcheck,
        "report": rep,
        "protons": [{"id": h.id, "delta": h.delta, "source": h.source, "multiplet": h.mult_category,
                     "nH_est": h.nH_est, "carbons": h.carbons, "spin_system": h.spin_system,
                     "flags": h.flags,
                     "evidence": {k: v for k, v in h.evidence.items()}} for h in protons],
        "carbons": [{"id": c.id, "delta": c.delta, "hsqc": c.mult, "protons": c.protons, "source": c.source} for c in carbons],
        "original_peaks": originals}}}
    if path_out:
        with open(path_out, "w") as fh:
            json.dump(new, fh, indent=1)
    return new, rep

# --------------------------------------------------------------------------- #
# 1D write-back: 1H multiplets/integrals and 13C peaks consistent with the 2D lists
# --------------------------------------------------------------------------- #
def _split_counts(total, ests):
    """Distribute an integer proton count over several protons in proportion to ests (each >= 1)."""
    n = len(ests)
    if not total or total < n:
        return [max(1, round(e)) for e in ests]
    w = [max(e, 0.3) for e in ests]
    raw = [total * x / sum(w) for x in w]
    out = [max(1, math.floor(x)) for x in raw]
    while sum(out) < total:
        i = max(range(n), key=lambda i: raw[i] - out[i])
        out[i] += 1
    while sum(out) > total:
        i = max((i for i in range(n) if out[i] > 1), key=lambda i: out[i] - raw[i])
        out[i] -= 1
    return out


def write_1d(new, doc, sp, protons, carbons, rep):
    """Rewrite the 1H multiplets/integrals and the 13C peak list so that the 1D spectra carry
    exactly the proton and carbon shifts used for the 2D peaks.  Returns the originals."""
    originals = {}
    h1 = sp.get(H1_1D)
    if h1 is not None:
        raw = doc[h1.key]
        r = h1.nh_scale or 1
        norm0 = raw["multiplets"].get("normValue") or raw["integrals"].get("normValue") or 1.0
        norm = norm0 * r
        orig_m = [m for m in raw["multiplets"]["data"].values()]
        orig_i = [i for i in raw["integrals"]["data"].values()]
        groups = {}
        for h in protons:
            groups.setdefault((h.mult_delta, tuple(h.range)), []).append(h)
        # HSQC rule: two separate protons on one CH2 carbon are 1H each (diastereotopic pair)
        one_h = {pid for c in carbons if c.mult == "CH2" and len(c.protons) == 2 for pid in c.protons}
        nh_notes = []
        mults, ints, used_i = [], [], set()
        for (md, rng), hs in groups.items():
            hs = sorted(hs, key=lambda h: -h.delta)
            m0 = next((m for m in orig_m if m.get("type", 0) == COMPOUND and md is not None
                       and abs(m["delta1"] - md) < 1e-9), None)
            i0 = None
            if m0 is not None:
                i0 = next((k for k, i in enumerate(orig_i) if i.get("type", 0) == COMPOUND
                           and i["rangeMin1"] <= m0["delta1"] <= i["rangeMax1"]), None)
            area = (m0.get("integralValue") if m0 else None) or \
                   (orig_i[i0]["intensity"] if i0 is not None else None)
            if m0 is None and i0 is None and any(h.strong for h in hs):   # merged region: the integrals inside it
                inside_i = [k for k, i in enumerate(orig_i) if i.get("type", 0) == COMPOUND
                            and rng[0] - 1e-6 <= (i["rangeMin1"] + i["rangeMax1"]) / 2 <= rng[1] + 1e-6]
                if inside_i:
                    used_i |= set(inside_i)
                    area = sum(orig_i[k]["intensity"] for k in inside_i)
            if len(hs) == 1:
                h = hs[0]
                if m0 is not None and m0.get("nH"):
                    nhs = [max(1, round(m0["nH"] / r))]
                else:
                    nhs = [max(1, round(h.nH_est or 1))]
            else:
                mn = m0.get("nH") if m0 is not None else hs[0].mult_nh
                total = round(mn / r) if mn else None
                cmap = {c.id: c for c in carbons}

                def fixed(h):          # HSQC CH2 rule: 2H for a lone CH2 proton, 1H each for a pair
                    cs = [cmap[c] for c in h.carbons if c in cmap]
                    if cs and cs[0].mult == "CH2":
                        return 2 if len(cs[0].protons) == 1 else 1
                    return None

                def hsqc_strength(h):  # this proton's HSQC cross-peak volume (CH3 ~ 3 x CH)
                    hq = sp.get(HSQC)
                    cs = [cmap[c] for c in h.carbons if c in cmap]
                    if hq is None or not cs:
                        return None
                    tolC = 2 * hq.pt_hz(0) / hq.sf[0]
                    v = max((abs(p.intensity) for p in hq.compound()   # peak height: a doublet picked on
                             if abs(p.d1 - cs[0].delta) <= tolC          # both lines must not count twice
                             and h.range[0] <= p.d2 <= h.range[1]), default=0)
                    return v or None
                fx = [fixed(h) for h in hs]
                rest = [k for k, x in enumerate(fx) if x is None]
                strengths = [hsqc_strength(hs[k]) for k in rest]
                if total and rest and all(strengths) and total - sum(x for x in fx if x) >= len(rest):
                    rem = total - sum(x for x in fx if x)
                    part = _split_counts(rem, [rem * v / sum(strengths) for v in strengths])
                    nhs = [x for x in fx]
                    for k, n in zip(rest, part):
                        nhs[k] = n
                else:
                    nhs = _split_counts(total, [fixed(h) or h.nH_est or 1 for h in hs])
            for k, h in enumerate(hs):
                if h.id in one_h and nhs[k] != 1:
                    nh_notes.append(f"{h.id} {h.delta:.4f} {nhs[k]}H -> 1H")
                    nhs[k] = 1
            tot = sum(nhs)
            lo, hi = (orig_i[i0]["rangeMin1"], orig_i[i0]["rangeMax1"]) if i0 is not None else rng
            if i0 is not None:
                used_i.add(i0)
            cuts = [hi] + [(a.delta + b.delta) / 2 for a, b in zip(hs, hs[1:])] + [lo]
            for k, (h, nh) in enumerate(zip(hs, nhs)):
                val = area * nh / tot if area else nh * norm
                if len(hs) == 1 and m0 is not None:
                    m = copy.deepcopy(m0)
                else:
                    m = {"delta1": 0.0, "nH": 0, "realH": 0, "integralValue": 0.0,
                         "category": h.mult_category if len(hs) == 1 and h.mult_category else "m",
                         "type": COMPOUND, "jlistcount": 0, "jvals": []}
                m.update({"delta1": h.delta, "nH": nh, "realH": nh, "integralValue": val, "type": COMPOUND})
                mults.append(m)
                ints.append({"intensity": val, "rangeMin1": cuts[k + 1], "rangeMin2": 0,
                             "rangeMax1": cuts[k], "rangeMax2": 0,
                             "delta1": (cuts[k] + cuts[k + 1]) / 2, "delta2": 0, "type": COMPOUND})
        mults.sort(key=lambda m: -m["delta1"])
        ints.sort(key=lambda i: -i["rangeMax1"])
        keep_m = [m for m in orig_m if m.get("type", 0) != COMPOUND]
        keep_i = [i for k, i in enumerate(orig_i) if i.get("type", 0) != COMPOUND]
        dropped_i = [i for k, i in enumerate(orig_i) if i.get("type", 0) == COMPOUND and k not in used_i]
        originals[H1_1D] = {"multiplets": raw["multiplets"], "integrals": raw["integrals"]}
        new[h1.key]["multiplets"] = {"datatype": "multiplets", "normValue": norm, "count": len(mults) + len(keep_m),
                                     "data": {str(k): m for k, m in enumerate(mults + keep_m)}}
        new[h1.key]["integrals"] = {"datatype": "integrals", "normValue": norm, "count": len(ints) + len(keep_i),
                                    "data": {str(k): i for k, i in enumerate(ints + keep_i)}}
        if nh_notes:
            rep.append("1D write-back: nH corrected - one of two protons on a CH2 carbon (HSQC): " + "; ".join(nh_notes))
        odd = [h for h in protons if h.nH_est and not h.strong and not getattr(h, "hidden_lines", 0)
               and h.split_c is None                      # split protons: counts set from the HSQC
               and abs(h.nH_est - round(h.nH_est)) > 0.3]
        rep.append(f"1D write-back: total {sum(m['nH'] for m in mults)} H in the compound multiplets"
                   + (" - nH doubtful for " + ", ".join(f"{h.id} {h.delta:.3f} (integral {h.nH_est:.2f})" for h in odd)
                      + "; check the 1H integrals / normalisation" if odd else ""))
        rep.append(f"1D write-back: 1H multiplets {len([m for m in orig_m if m.get('type', 0) == COMPOUND])} -> "
                   f"{len(mults)} (one per proton, delta = refined shift), integrals rebuilt to match"
                   + (f"; integral normalisation corrected (normValue x {r})" if r > 1 else "")
                   + (f"; compound integrals without a proton dropped: "
                      + ", ".join(f"{i['rangeMin1']:.3f}-{i['rangeMax1']:.3f}" for i in dropped_i) if dropped_i else ""))
    c13 = sp.get(C13_1D)
    if c13 is not None and carbons:
        raw = doc[c13.key]
        by_key = {c.peak_key: c for c in carbons if c.peak_key is not None}
        out, n_re = [], 0
        for k, p in raw["peaks"]["data"].items():
            q = copy.deepcopy(p)
            if p.get("type", 0) == COMPOUND:
                if k in by_key:
                    q["delta1"] = by_key[k].delta
                    q["annotation"] = by_key[k].id
                elif k in c13.excluded:
                    t, why = c13.excluded[k]
                    q["type"], q["annotation"] = t, why
                    n_re += 1
            out.append(q)
        originals[C13_1D] = raw["peaks"]
        new[c13.key]["peaks"] = {"datatype": "peaks", "count": len(out),
                                 "data": {str(i): p for i, p in enumerate(out)}}
        rep.append(f"1D write-back: 13C {len(by_key)} compound peaks (one per carbon); {n_re} compound-typed peaks "
                   f"retyped (satellite/noise -> artifact, minor -> impurity, line-fit splits -> artifact, "
                   f"solvent position -> solvent)")
    return originals


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input")
    ap.add_argument("output", nargs="?")
    ap.add_argument("--report")
    a = ap.parse_args()
    if not a.output:
        a.output = (a.input[:-5] if a.input.lower().endswith(".json") else a.input) + "_refined.json"
    try:
        _, rep = refine(a.input, a.output)
    except LabelError as e:
        text = "\n".join([f"simpleNMRclaudePeakRefine v{VERSION}: {a.input}",
                          "ERROR: the spectrum labels chosen in Mnova do not match the data.",
                          "Nothing was refined and no output file was written.",
                          "Fix the labels in the Mnova spectrum-selection dialog (or set the spectrum to SKIP) "
                          "and export again.", ""]
                         + [f"  - {x}" for x in e.errors]
                         + ([""] + [f"  warning: {x}" for x in e.warnings] if e.warnings else []))
        print(text, file=sys.stderr)
        if a.report:
            with open(a.report, "w") as fh:
                fh.write(text + "\n")
        sys.exit(2)
    text = "\n".join(rep)
    print(text)
    print(f"\nRefined JSON written to {a.output}")
    if a.report:
        with open(a.report, "w") as fh:
            fh.write(text + "\n")

