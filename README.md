# Open-source OPLS/2020 force-field generation

`opls_ff_gen` parameterizes organic molecules with the OPLS-AA force field
and writes GROMACS-ready topologies. It uses the charge model CM1A with the 
1.14 scaling and LBCC correction (1.14*CM1A-LBCC) for neutral molecule and 
the unscaled charge for charged molecules.

Bonded parameter, whenever availiable comes from the published OPLS/2020 tables. 
Anything the tables do not cover is filled per molecule:
Lennard-Jones parameters for untypable atoms come from GAFF2. And if QM fit is not
enabled, any missing bonded terms come from GAFF2 or recorded QM values or substitution/empirical
fallbacks. Default option is to not run QM fit, but with QM enabled with `--qm xtb|dft`, bonded 
terms are fitted on the fly to QM using xTB (only for bonded terms that are missing from the 
original OPLS/2020 publication).

**Layout.** All source code is under `src/`: the pipeline package
`src/opls_ff_gen/` (data tables in `src/opls_ff_gen/data/`) and the build,
batch, check and benchmark tools in `src/scripts/`. `src/README.md` documents
every Python file. Tests are in `tests/`, validation notes in `validation/`,
the handoff and project status in `PROGRESS.md`.

## Complete workflow

`opls_ff_gen` combines:

1. explicit-hydrogen 3D structures from a structure file, or a deterministic
   RDKit conformer from a SMILES string, optimized with GFN2-xTB by default;
2. atom typing: OPLS-AA SMARTS rule plus the atom-type archive, then
   refinement to OPLS/2020 type numbers (new alkane, alkyl-arene,
   cyclobutane, alkyne-adjacent and allene types);
3. classes and Lennard-Jones parameters from the OPLS/2020 table
   (`oplsaa2020.par`, installed with `opls-ff-gen install-si`); atoms that no OPLS type
   covers get a GAFF2 type from a built-in rule table (no AmberTools needed), then
   the OPLS/2020 type analogous to that GAFF2 type where there is one, else
   GAFF2 LJ;
4. fixed-geometry AM1 single points with MOPAC and CM1A charges computed by
   the project (`src/opls_ff_gen/cm1a.py`; see "CM1A charges" below), 1.14 scaling for
   neutral molecules, graph-symmetry averaging and the 19 published
   localized bond-charge corrections;
5. RDKit enumeration of bonds, angles, proper torsions, planar impropers and
   1-4 pairs;
6. bonds, angles, torsions and impropers from the OPLS/2020 tables
   (`oplsaa2020.sb` and the torsion section of `oplsaa2020.par`);
7. for proper torsions missing from those tables: a GAFF2-derived general
   X-B-C-X term, then AmberTools parmchk2-style substitution
   (`src/opls_ff_gen/gaff2_torsions.py`, `src/opls_ff_gen/parmchk.py`; see "Missing
   torsions" below);
8. for every bond, angle and proper torsion still missing: by default
   (`--qm none`) the gaff2.dat and recorded-QM tables, then the fallbacks of
   `src/opls_ff_gen/bonded_substitution.py` (torsions: zero as a last
   resort), so every term gets a value; with `--qm xtb|dft` instead, an
   on-the-fly QM fit of every term the OPLS/2020 table lacks, skipping step 7
   (`src/opls_ff_gen/qm_fit.py`, see below); impropers are never fitted; and
9. a self-contained `.itp` (atom types with LJ, explicit parameters on every
   bonded line), matching `.gro`, a `.top` with the OPLS `[ defaults ]`, and
   an auditable JSON record. GROMACS's own `oplsaa.ff` is not used.


## OPLS/2020 parameters

Source: Supporting Information of W. L. Jorgensen, M. M. Ghahremanpour,
A. Saar and J. Tirado-Rives, "OPLS/2020 Force Field for Unsaturated
Hydrocarbons, Alcohols, and Ethers", *J. Phys. Chem. B* (2024),
https://doi.org/10.1021/acs.jpcb.3c06602. The two text files (BOSS
`oplsaa.par` and `oplsaa.sb` formats, dated Oct 1, 2023) are installed by
each user (below). They hold about 990 atom types, 465 bond, 1417
angle and 1225 torsion entries, far more than GROMACS's `oplsaa.ff`.

**The parameter files are not distributed** (ACS Supporting Information; no
permission to redistribute). Download the SI from the article page and
install the two text files once:

```bash
opls-ff-gen install-si jp3c06602_si_002.txt jp3c06602_si_003.txt
```

`install-si` identifies each file (`_002`: non-bonded and torsions, `_003`:
bonds and angles), checks it against the SHA-256 of the Oct 1, 2023 release
the tool was validated with (a different version needs `--force`, and runs
then warn), and copies it to `~/.local/share/opls-ff-gen/opls2020/`
(`$XDG_DATA_HOME`; `%APPDATA%` on Windows), where it survives upgrades.
`OPLS_FF_GEN_SI_DIR` selects another directory (`install-si --dest DIR`).
In a development checkout, copies in `src/opls_ff_gen/data/opls2020/`
(`oplsaa2020.par`, `oplsaa2020.sb`) are used if no user installation
exists. `opls-ff-gen doctor` shows which files are read.

## Optional run-time QM fitting of missing bonded terms

The OPLS/2020 table is always the first choice for every bonded term. With
`--qm xtb` or `--qm dft` (default `none`), every term the OPLS/2020 table
lacks goes straight to the QM fit.

## Command-line use

Installing from a release (conda environments, MOPAC, the SI files,
programs you already have, `doctor`), or running the source without
installing (`python run_opls_ff_gen.py ...` with any Python that has RDKit
and NumPy): see `INSTALL.md`.


```bash
python -m opls_ff_gen \
  --sdf molecules_to_test/ethanol.sdf \
  --name ethanol \
  --output-dir generated/ethanol
```

From a SMILES string:

```bash
python -m opls_ff_gen \
  --smiles 'c1ccc(cc1)OCCO' \
  --name phenoxyethanol \
  --output-dir generated/phenoxyethanol
```

`--charge-model auto` is the default. It selects 1.14*CM1A-LBCC for neutral
molecules and unscaled CM1A without LBCC for ions. `cm1a` and `cm1a-lbcc` can
also be selected explicitly.

## Python API

```python
from opls_ff_gen import generate_gromacs_topology, molecule_from_smiles

molecule = molecule_from_smiles("CCO", name="ethanol")
result = generate_gromacs_topology(
    molecule,
    "generated/ethanol",
    name="ethanol",
    charge_model="auto",
)

for atom in result.assignment.atoms:
    print(atom.index, atom.opls_type, atom.final_charge_e)
print(result.itp_path)
```

## Outputs

Each complete run creates:

- `NAME.itp`: molecule atom types, partial charges, bonds, pairs, angles,
  proper torsions, and improper torsions;
- `NAME.gro`: coordinates in a padded orthorhombic box;
- `NAME.top`: a minimal one-molecule topology with the OPLS `[ defaults ]`
  (nbfunc 1, comb-rule 3, gen-pairs yes, fudgeLJ/fudgeQQ 0.5) that includes
  `NAME.itp`;
- `NAME.topology.json`: every interaction, parameter value, OPLS class match,
  rule and source line in the OPLS/2020 files;
- `NAME.assignment.json`: the stable machine-readable handoff for the bonded
  parameter stage;
- `NAME.assignment.csv`: the same atom-level data in tabular form;
- `NAME.assigned.sdf`: the 3D structure with OPLS types, classes, and final
  charges stored as SDF properties;
- `mopac/NAME.mop`, `NAME.out`, `NAME.aux`: the complete AM1 calculation
  record (the CM1A charges are recomputed from the density matrix in the
  `.aux` file).

Every atom record contains its preserved index, element, coordinates,
OPLS type and class, mass, sigma, epsilon, raw CM1A charge, scaled charge,
symmetry-averaged charge, LBCC increment, and final charge.



