# Assignment-module validation

The combined `opls_ff_gen` module was run end to end on both LigParGen
benchmarks. It reproduced all canonical Foyer atom types and reconstructed the
LigParGen 1.14*CM1A-LBCC charges within the established tolerance.

| Molecule | Atoms | Canonical type matches | LBCC applications | Charge MAE (e) | Maximum charge difference (e) |
|---|---:|---:|---:|---:|---:|
| Ethanol | 9 | 9 | 1 | 0.00005754 | 0.00011851 |
| Phenoxyethanol | 20 | 20 | 6 | 0.00005889 | 0.00023076 |

The signed comparison used for the absolute statistics is
`module charge - LigParGen ITP charge`. Both final molecular charge sums are
zero to floating-point precision.

The benchmark outputs are under `validation/assignments`. A separate
SMILES-only acetone test is under `validation/new_molecule_smoke/acetone`; it
assigned the expected carbonyl types and applied one `C-O` LBCC.
