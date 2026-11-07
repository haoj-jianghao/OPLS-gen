# GROMACS OPLS-AA bonded-parameter comparison

Parameter source:
`/home/jiang/software/gromacs-2026.3/share/top/oplsaa.ff/ffbonded.itp`

Molecular connectivity is enumerated from the explicit-hydrogen SDF records
with RDKit. Atom classes are the Foyer OPLS-AA classes already validated
against LigParGen.

## Topology enumeration

| Molecule | Bonds | Angles | Proper torsions | 1-4 pairs |
|---|---:|---:|---:|---:|
| Ethanol | 8 | 13 | 12 | 12 |
| Phenoxyethanol | 20 | 32 | 41 | 38 |

These sets agree exactly with the interactions in the corresponding
LigParGen ITP files. Each pair is a unique atom pair with graph distance three
and is also the terminal pair of at least one enumerated proper torsion.

## Numerical comparison

All 28 bond assignments and all 45 angle assignments reproduce the LigParGen
values exactly at the precision printed in the ITP files.

Most proper torsions also agree. Six individual torsions, belonging to two
parameter families, differ:

| OPLS classes | Occurrences | GROMACS `C0..C5` (kJ/mol) | LigParGen `C0..C5` (kJ/mol) |
|---|---:|---|---|
| `HC-CT-OH-HO` | 4 | `0.94140 2.82420 0 -3.76560 0 0` | `0.736 2.209 0 -2.946 0 0` |
| `CA-CA-OS-CT` | 2 | `12.552 0 -12.552 0 0 0` | `14.100 0 -9.079 0 -5.021 0` |

The GROMACS file defines aromatic planarity through
`improper_Z_CA_X_Y = 180.0 4.60240 2`. LigParGen uses
`180.0 10.460 2` for each of the six aromatic carbons in phenoxyethanol.

These are source-force-field differences rather than graph-enumeration or unit
conversion errors. The module retains the GROMACS values and their source
lines. Exact LigParGen compatibility would require a documented BOSS-compatible
override table for these parameter families.
