# AMSOL -> MOPAC AM1 + own CM1A: equivalence evidence (2026-10-18)

CM1A charges are now computed by `src/opls_ff_gen/cm1a.py` from a MOPAC 23.2.5
AM1 single point (density matrix from the `.aux` file). CM1A follows Storer
et al., *J. Comput.-Aided Mol. Des.* 9, 87 (1995).

- Branch: `mopac-cm1a`.
- Reference: tag `amsol-reference`, the last commit using AMSOL 7.1.
- Check script: `src/scripts/check_cm1a_equivalence.py`. It needs the AMSOL 7.1
  executable and carries its own copy of the pre-replacement AMSOL code.

## 1. Parameters

All 17 published CM1A parameters in `cm1a.py` equal the values in AMSOL 7.1's
`PARAMS.i`. These were extracted programmatically from the DATA statements,
not copied by hand:
- c (scale): N, F, S, Cl, Br, I;
- d (offset): O, F, S, Cl, Br, I;
- pair offsets: H-N, H-O, H-Si, O-S, N-O;
- the N-C triple-bond switch.

## 2. Charges, molecule by molecule

Each molecule gets one RDKit 3D geometry: the SDF's if it is 3D with explicit
H, otherwise ETKDG + MMFF94. Both programs run a single AM1 SCF on it with
the formal net charge, and the raw CM1A charges are compared.

The acceptance criterion, set before the run, was every atom within
0.002 e.

| | |
|---|---|
| Molecules (RDKit NCI first_5K, ChEMBL2321810, PubChem 200, BZR, CDK2, EGFR, NCI first_200; project test sets, rare-type set, project SDFs) | 7,208 |
| Compared (skipped: 168 unsupported elements, 137 multi-fragment, 30 RDKit embedding failures, 8 unreadable) | **6,803 molecules, 228,190 atoms** |
| Atoms over 0.002 e | **0** |
| Largest per-atom difference | **0.0009 e** (P in a cyclic phosphazene; the largest differences are all on P) |
| Per-molecule maximum: median, 90th, 99th, 99.9th percentile | 0.00005, 0.00025, 0.00033, 0.00054 e |
| AMSOL SCF failed, MOPAC converged | 62 |
| MOPAC failed, AMSOL succeeded | **0** |

Result file: `equivalence_all.json`.

Two harness problems were found and fixed before this final run; neither
affects the charge code:
- **2D input records.** `pubchem.200` and `first_200.props` are 2D drawings
  without hydrogens. Hydrogens added onto flat coordinates gave absurd
  structures, with heats of formation of 3,000-12,000 kcal/mol in both
  programs and differences up to 1.4 e. Such records now get an ETKDG 3D
  conformer.
- **SCF convergence.** MOPAC's default SCF convergence left up to 0.001 e of
  noise in large conjugated molecules (protoporphyrin IX: 0.0024 e). The
  workflow now uses `SCFCRT=1.D-10` (porphyrin: 0.0007 e).

## 3. End-to-end: the 100-molecule xTB workflow

The batch `benchmarks/opls2020_xtb_mopac` (MOPAC charges) was compared with
`benchmarks/opls2020_xtb_foyerfree` (AMSOL charges, otherwise the same code)
using `src/scripts/compare_batches_numeric.py`. Result file: `end_to_end_100.json`.

- **All 100 generated:** pass `grompp` (0 warnings) and a zero-step `mdrun`.
- **Identical:** atom types, classes, rules, every table-derived bonded term,
  and the QM-fitted bonds and angles (these come from the Hessian, which does
  not depend on charges).
- **Final charges:** largest difference 0.0003 e; median per molecule
  0.00003 e.
- **Fitted torsions:** unchanged within 0.01 kJ/mol in 99 molecules.
  Carbamazepine's urea-aryl fit moved by up to 54 kJ/mol in RB coefficients,
  but its scan RMSE improved (11.0 -> 9.9 and 4.0 -> 3.4 kJ/mol). That fit is
  poorly determined by its scans, so a 1e-4 e charge change sends the damped
  fit to a different, equally good solution. This is a property of the torsion
  fit (open item), not of the charges.
- **Scan RMSE:** 113 scans, mean change -0.01 kJ/mol, largest 1.1 kJ/mol
  (carbamazepine).

## 4. Package without AMSOL

All 64 tests pass with a guard that raises if anything tries to run AMSOL.
The `--amsol` option and the AMSOL code were removed from the package.
