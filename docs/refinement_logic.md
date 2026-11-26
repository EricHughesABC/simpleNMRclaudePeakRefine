# simpleNMRclaudePeakRefine: how peak positions are refined and harmonised

Automatic peak-picking refinement for Mnova exports used by simpleNMR. Program version 0.29, November 2026. Also available as [refinement_logic.pdf](refinement_logic.pdf) (written for 0.25; the tolerance note in section 2 is new in 0.28, the DEPT-135 notes are new in 0.29).

> The program never trusts a single spectrum. From all the available evidence it builds one list of proton shifts and one list of carbon shifts. It then rewrites every spectrum so that its peaks sit exactly on those lists. Each shift is decided once, and every spectrum is made to agree with it.

## 1. Read and check the input
- The input is a JSON file exported from Mnova. The role of each spectrum comes from the label the user chose in the Mnova spectrum-selection dialog (`HSQC_0`, `HMBC_0`, …). Spectra labelled SKIP, unknown labels, NOESY and second datasets are left untouched.
- Each label is checked against the spectrum's own metadata (nucleus, dimension, pulse programme). A mislabelled spectrum stops the run with an error report and exit code 2, so nothing wrong is written.
- The 1H 1D and the HSQC are required. The 13C 1D, DEPT-135, pureshift 1H, HMBC, COSY, HSQC-CLIP-COSY and DEPT-HSQC (CH₃ only) spectra are used when present.
- If the export includes the proposed structure, only its molecular formula is used (computed with RDKit). It is used for checks, never for peak decisions, so the picking stays unbiased between isomers.

## 2. Set tolerances from each spectrum's own resolution
Every matching window is a number of digital points in that spectrum's dimension (`sw_Hz / spectral_size`), with small minimum values in Hz, converted to ppm with that dimension's spectrometer frequency. A 1D shift is precise to well under 1 Hz, while a 2D cross-peak is only as good as its point grid. So the 1D spectra decide the positions, and the 2D spectra decide which positions are real.

The proton window used to give an HSQC (or other H-C 2D) pick to a proton when the pick lies outside every 1H multiplet is the same rule simpleNMR uses to match HMBC/COSY protons to the HSQC rows: `tolerance.py`, copied unchanged from simpleNMRtools `core/tolerance.py` (3 points of the coarser f2 grid of the spectrum and the HSQC, at least 4 Hz). The tests that decide whether two HSQC signals are one or two (an HSQC row against a pureshift peak, and pooling HSQC picks that have no proton) are resolution tests, not matching windows, so they keep their own tighter value (`split_tol_pts`, 2 points, at least 4 Hz): a wider value merges real protons, as in EWK18 (1.1118 on C12, 15 Hz from 1.1333 on C14), and pools separate weak picks, as in santonin (1.3711 and 1.3940). The HMBC and COSY windows work on clusters of multiplet lines and keep their own values.

## 3. Build the carbon list
**With a 13C 1D spectrum:**
- Peaks above a threshold (2% of the median of the five strongest) become candidates. Solvent peaks and 13C satellites are excluded.
- Lines of one carbon that the line fit split are merged into one carbon. When most peaks are split (a distorted lineshape, typically from poor shimming), the merge tolerance widens and the report gives a warning.
- Two kinds of carbon are typed as impurities: weak minor peaks with no 2D support, and weak carbons whose only HSQC peaks sit under protons that already belong to another carbon.

**Without a 13C 1D spectrum:** protonated carbons come from the HSQC f1 positions (moved to the DEPT-135 positions when a DEPT-135 is supplied). Quaternary carbons come from HMBC f1 rows, but only when a row has at least two HMBC picks.

**DEPT-135 (v0.29):**
- DEPT picks above the same threshold as the 13C are matched to the 13C carbons (within 3 Hz or 3 DEPT points). Each carbon takes the sign of its strongest DEPT peak: negative = CH₂, positive = CH or CH₃. Further picks on the same carbon are line-fit splits.
- A DEPT peak with no 13C carbon becomes a carbon when an HSQC pick confirms it and no carbon lies within the HSQC f1 tolerance. The 13C pick sometimes misses protonated carbons: in the menthyl anthranilate export 116.85, 116.44 and 74.14 came only from the DEPT. Such carbons are added to the written 13C list, and the report asks for them to be checked.
- **Non-edited HSQC** (no negative compound picks) with a DEPT-135: each HSQC pick takes the sign of the nearest DEPT peak within the HSQC f1 tolerance. The HSQC rules (a CH₂ carbon carries up to two protons, a CH or CH₃ carbon one) then work as for an edited HSQC. The written HSQC peaks keep their measured sign; the CH₂ peaks are noted "CH2 from DEPT-135".
- **Edited HSQC** with a DEPT-135: the HSQC sign is used, and the DEPT is a cross-check. Disagreements, and protonated carbons without a DEPT peak, are reported.
- A non-edited HSQC without a DEPT-135 gives a warning: CH₂ carbons cannot be recognised, so the second proton of a diastereotopic CH₂ is typed artifact.
- Whether the HSQC is edited comes from its name (Mnova experiment `HSQC-EDITED`, Bruker `hsqced…`). In a non-edited HSQC, a negative pick next to a stronger positive peak in the same f2 column (within 2.5 × the f1 resolution) is an f1 truncation wiggle and is typed artifact. With two or more negative picks that are not wiggles, the HSQC is treated as edited and the report warns. (Colchicine: two wiggles of the OCH₃ peaks had made v0.29 treat the HSQC as edited and create two false protons.)
- A DEPT peak at a 13C position that the 13C step rejected (solvent, impurity, noise) follows the 13C decision; it does not bring the carbon back.
- When several carbons lie within the HSQC f1 tolerance, the one whose DEPT sign fits the HSQC group is chosen first; a carbon with no DEPT peak (quaternary) comes after. Colchicine: 134.90 (DEPT CH) rather than 134.70 (quaternary).

## 4. Build the proton list
Candidate protons come from these sources, in order of precision:
1. pureshift peaks inside Mnova multiplets;
2. Mnova multiplet centres;
3. 1H signals with no Mnova multiplet but with H–C 2D support;
4. multiplets built from the 1H line picks, if there are no multiplets at all.

The candidates are then split, merged or removed according to the HSQC:
- **One proton per HSQC carbon group.** A multiplet covering several carbons is split, and the protons are centred jointly using the HSQC rows.
- **Degenerate protons.** One shift with HSQC groups on clearly different carbons becomes one proton per carbon at the same shift.
- **Strongly coupled regions.** The pureshift is not trusted where its negative artefacts are large. The protons are taken from the HSQC carbons and positioned at the mean of the HSQC rows and HMBC columns.
- **Surplus pureshift peaks.** In a multiplet already filled by protons that have carbons, extra carbon-less peaks are treated as artefacts. Where several peaks compete, the one kept is nearest the consensus of the Mnova centre, the HSQC rows and the HMBC columns.
- **Weak HSQC groups.** A weak group can still become its own proton if it is clearly separate and not explained by a stronger peak elsewhere.
- **Carbons closer than the HSQC f1 tolerance (v0.29).** One HSQC carbon group whose picks form f2 clusters more than 16 Hz apart, each nearest a different 13C carbon (same f1 order, at least half an f1 point apart; with a DEPT-135 both protonated), holds one proton per carbon. Colchicine: OCH₃ 56.51/3.871 and 56.30/3.833, which Mnova had picked as one 6H "doublet".
- **Line centres for split protons (v0.29).** When protons are centred jointly on the 1H line picks, the choice that explains the most line intensity wins (in steps of 10%), then the most lines: two small lines no longer outvote a strong singlet.
- **Carbon-less protons** (OH, NH) are kept if they are real 1H signals.

## 5. Pair and validate with the 2D spectra
- Each HSQC peak is matched to a proton and a carbon within the resolution tolerance.
- HMBC peaks are matched to a proton in f2 and to a carbon in f1. COSY, HSQC-CLIP-COSY and DEPT-HSQC peaks are matched the same way.
- Picks that fit no pair, or are explained as artefacts, are typed artifact. These include t1 noise, one-bond breakthrough and peaks too weak relative to their row or column. Such picks are never moved to a far-away resonance.

## 6. Write the results back
- **2D spectra:** the peaks are cleared and replaced. Every compound peak sits exactly at its proton/carbon or proton/proton shifts. The original picks are kept in a `peakRefinement` block in the file. As a safety net, any compound peak that is not on the lists is typed artifact.
- **1H 1D:** one multiplet and one integral per proton. The proton count comes from Mnova, then from the integrals. When a multiplet is split, a CH₂ carbon fixes its count at 2, and the rest of the integral is shared by HSQC peak height.
- **13C 1D:** one compound peak per carbon. The remaining peaks are retyped as impurity, solvent or artifact.
- **DEPT-135:** one compound peak per protonated carbon, at the carbon shift, with its measured sign. Split lines and noise are retyped artifact (simpleNMR counts the DEPT peaks against the protonated carbons of the molecule).
- **Report:** every decision with its reason, and the offset of each proton from each source, in Hz.

## 7. Final checks (reported, not applied)
- A consistency check that every compound 2D peak lies on the proton and carbon lists.
- If a molecular formula is available, three checks: the carbon count against the formula (an upper limit, since symmetry can make carbons equivalent), the total number of hydrogens, and CH₃ carbons whose proton count is not 3.

> The principle throughout is that a peak is never moved further than the spectrum's resolution allows. A peak that does not fit is typed as an artifact, not pushed onto the nearest resonance.
