# GROMACS vs OpenMM as the MM engine of the torsion fit (2026-10-19)

The restrained MM relaxations of the on-the-fly torsion fit can run with
OpenMM (default) or GROMACS (`--mm-engine gromacs`, `src/opls_ff_gen/mm.py`).
The checks follow the verification plan agreed on 2026-10-18; all were made
with `src/scripts/check_mm_engines.py`.

GROMACS settings for a gas-phase molecule: periodic box at least twice a
cut-off that exceeds the molecule by 1 nm (every pair interacts once, no
periodic image does); plain cut-off Coulomb (reaction field, epsilon_rf = 1)
and LJ with no potential modifiers; L-BFGS, emtol 1e-3 kJ/mol/nm;
`[ dihedral_restraints ]` for the driven and locked dihedrals. GROMACS 2026.3
(mixed precision), OpenMM Reference platform.

## Data

- `scans_100.json`: the 46 molecules of the 100-molecule xTB batch of
  2026-10-19 that have torsion scans (113 scans, 2,427 scan geometries),
  their final fitted topologies.
- `ions_scans.json`: 6 ions from `ions.csv` that have scans (10 scans, 240
  geometries).
- `ions_perturbed.json`: all 8 ions of `ions.csv`, 20 random displacements
  (0.05 A per coordinate) of the input geometry each (160 geometries).

## Single points (energies per component and forces)

Largest absolute difference over all geometries of a set (kJ/mol; forces
kJ/mol/nm):

| Component | 100-molecule scans | Ion scans | Perturbed ions |
|---|---|---|---|
| Bonds | 1.2e-3 | 2.8e-4 | 1.6e-3 |
| Angles | 7.2e-4 | 2.0e-4 | 5.6e-4 |
| Proper torsions (RB) | 9.6e-4 | 9.0e-5 | 1.5e-4 |
| Impropers | 8.7e-4 | 8.5e-5 | 8.7e-5 |
| LJ 1-4 | 4.0e-4 | 1.1e-4 | 1.4e-4 |
| Coulomb 1-4 | 6.1e-4 | 3.1e-4 | 3.1e-4 |
| LJ, longer range | 2.4e-3 | 2.0e-4 | 5.6e-5 |
| Coulomb, longer range | 3.6e-4 | 3.1e-4 | 2.8e-4 |
| Total potential | 2.8e-3 | 6.8e-4 | 1.7e-3 |
| Forces (largest force) | 0.21 (12,539) | 0.10 (3,948) | 0.15 (15,122) |

Every component agrees to the precision of GROMACS's mixed-precision build
(about 1e-5 relative for energies and 2e-5 for forces), charged molecules
included.

## Relaxed scan profiles

The restrained relaxed profile of every scan with both engines (the fit's
settings, hydrogen rotors locked), compared relative to each profile's
minimum:

- 123 scans (2,667 points). Median of the per-scan largest difference
  0.003 kJ/mol; 122 scans below 0.1 kJ/mol; ions below 0.007 kJ/mol.
- One point differs: alpha-terpineol, scan 6-7, phi = 35 deg, 22.6 kJ/mol
  (relaxed geometries 1.07 A RMSD apart). The ring changes conformation
  between 30 and 40 deg (the profile drops from -54 to -74 kJ/mol); from the
  same QM start OpenMM's minimiser crossed into the other ring conformer,
  GROMACS's did not. Both are local minima and all other points of the scan
  agree to 0.003 kJ/mol: a minimiser-path difference at a bistable point,
  not an energy difference.

## Cost

GROMACS runs one `gmx grompp` + `gmx mdrun` per scan point (the fit runs the
points in parallel). Methyl thiocyanate (one scan, 8 threads): torsion fit
3.9 s with GROMACS, 0.5 s with OpenMM, identical parameters to
0.003 kJ/mol. The 100-molecule timings are below.

## End-to-end: 100-molecule batch with each engine

`src/scripts/run_batch_opls2020.py` on `batch_2026_10_17` + `batch_2026_10_18`
(xTB, SMILES optimized with GFN2-xTB, the default since 2026-10-19), once with
`--mm-engine openmm` and once with `--mm-engine gromacs`, 8 threads:

| | OpenMM | GROMACS |
|---|---|---|
| generate / grompp (0 warnings) / mdrun | 100 / 100 / 100 | 100 / 100 / 100 |
| Bonded rule counts | identical for all 100 molecules | |
| Scans; fitted RMSE mean / median (kJ/mol) | 113; 2.448 / 1.45 | 113; 2.457 / 1.45 |
| Wall time, whole batch | 14.5 min | 30.5 min |
| of which torsion fits | 3.4 min | 19.6 min |

Per scan, the fitted RMSE agrees within 0.01 kJ/mol for 97 of 113 scans and
within 0.43 kJ/mol for all. The four scans above 0.1 kJ/mol are carbamazepine
(the known non-unique fit) and acetazolamide, where the thiadiazole ring
torsion is split differently between CW-SA-CW-N and CW-SA-CW-SY
(V3 9.7/9.0 vs 13.3/6.6 kJ/mol, similar sum and objective): a ring torsion
the scans barely determine. The non-ring terms of both molecules agree to
about 0.05 kJ/mol.

While running this batch, the SMILES route was found to depend on the
optimizer: RDKit's MMFF set-up applies the MMFF94 aromaticity model in place,
so with the xTB default the 4-hydroxycoumarin ring of warfarin became
aromatic and two atoms fell back to GAFF2. SMILES inputs now always get the
MMFF94 aromaticity model (`structure.molecule_from_smiles_optimized`); the
explicit model reproduces the earlier flags for all 100 molecules, and
warfarin was rerun with both engines (numbers above include the rerun).
Structure files (SDF/MOL2/PDB) now get the same model, so all 100 molecules
type identically from SMILES and from SDF.
