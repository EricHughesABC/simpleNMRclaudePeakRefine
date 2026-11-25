# simpleNMRclaudePeakRefine: how peak positions are refined and harmonised

Automatic peak-picking refinement for Mnova exports used by simpleNMR. Program version 0.28, November 2026. Also available as [refinement_logic.pdf](refinement_logic.pdf) (written for 0.25; the tolerance note in section 2 is new in 0.28).

> The program never trusts a single spectrum. From all the available evidence it builds one list of proton shifts and one list of carbon shifts. It then rewrites every spectrum so that its peaks sit exactly on those lists. Each shift is decided once, and every spectrum is made to agree with it.

## 1. Read and check the input
- The input is a JSON file exported from Mnova. The role of each spectrum comes from the label the user chose in the Mnova spectrum-selection dialog (`HSQC_0`, `HMBC_0`, …). Spectra labelled SKIP, unknown labels, NOESY, DEPT135 and second datasets are left untouched.
- Each label is checked against the spectrum's own metadata (nucleus, dimension, pulse programme). A mislabelled spectrum stops the run with an error report and exit code 2, so nothing wrong is written.
- The 1H 1D and the HSQC are required. The 13C 1D, pureshift 1H, HMBC, COSY, HSQC-CLIP-COSY and DEPT-HSQC (CH₃ only) spectra are used when present.
- If the export includes the proposed structure, only its molecular formula is used (computed with RDKit). It is used for checks, never for peak decisions, so the picking stays unbiased between isomers.

## 2. Set tolerances from each spectrum's own resolution
Every matching window is a number of digital points in that spectrum's dimension (`sw_Hz / spectral_size`), with small minimum values in Hz, converted to ppm with that dimension's spectrometer frequency. A 1D shift is precise to well under 1 Hz, while a 2D cross-peak is only as good as its point grid. So the 1D spectra decide the positions, and the 2D spectra decide which positions are real.

The proton window used to give an HSQC (or other H-C 2D) pick to a proton when the pick lies outside every 1H multiplet is the same rule simpleNMR uses to match HMBC/COSY protons to the HSQC rows: `tolerance.py`, copied unchanged from simpleNMRtools `core/tolerance.py` (3 points of the coarser f2 grid of the spectrum and the HSQC, at least 4 Hz). The tests that decide whether two HSQC signals are one or two (an HSQC row against a pureshift peak, and pooling HSQC picks that have no proton) are resolution tests, not matching windows, so they keep their own tighter value (`split_tol_pts`, 2 points, at least 4 Hz): a wider value merges real protons, as in EWK18 (1.1118 on C12, 15 Hz from 1.1333 on C14), and pools separate weak picks, as in santonin (1.3711 and 1.3940). The HMBC and COSY windows work on clusters of multiplet lines and keep their own values.

## 3. Build the carbon list
**With a 13C 1D spectrum:**
- Peaks above a threshold (2% of the median of the five strongest) become candidates. Solvent peaks and 13C satellites are excluded.
- Lines of one carbon that the line fit split are merged into one carbon. When most peaks are split (a distorted lineshape, typically from poor shimming), the merge tolerance widens and the report gives a warning.
- Two kinds of carbon are typed as impurities: weak minor peaks with no 2D support, and weak carbons whose only HSQC peaks sit under protons that already belong to another carbon.

**Without a 13C 1D spectrum:** protonated carbons come from the HSQC f1 positions. Quaternary carbons come from HMBC f1 rows, but only when a row has at least two HMBC picks.

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
- **Carbon-less protons** (OH, NH) are kept if they are real 1H signals.

## 5. Pair and validate with the 2D spectra
- Each HSQC peak is matched to a proton and a carbon within the resolution tolerance.
- HMBC peaks are matched to a proton in f2 and to a carbon in f1. COSY, HSQC-CLIP-COSY and DEPT-HSQC peaks are matched the same way.
- Picks that fit no pair, or are explained as artefacts, are typed artifact. These include t1 noise, one-bond breakthrough and peaks too weak relative to their row or column. Such picks are never moved to a far-away resonance.

## 6. Write the results back
- **2D spectra:** the peaks are cleared and replaced. Every compound peak sits exactly at its proton/carbon or proton/proton shifts. The original picks are kept in a `peakRefinement` block in the file. As a safety net, any compound peak that is not on the lists is typed artifact.
- **1H 1D:** one multiplet and one integral per proton. The proton count comes from Mnova, then from the integrals. When a multiplet is split, a CH₂ carbon fixes its count at 2, and the rest of the integral is shared by HSQC peak height.
- **13C 1D:** one compound peak per carbon. The remaining peaks are retyped as impurity, solvent or artifact.
- **Report:** every decision with its reason, and the offset of each proton from each source, in Hz.

## 7. Final checks (reported, not applied)
- A consistency check that every compound 2D peak lies on the proton and carbon lists.
- If a molecular formula is available, three checks: the carbon count against the formula (an upper limit, since symmetry can make carbons equivalent), the total number of hydrogens, and CH₃ carbons whose proton count is not 3.

> The principle throughout is that a peak is never moved further than the spectrum's resolution allows. A peak that does not fit is typed as an artifact, not pushed onto the nearest resonance.
