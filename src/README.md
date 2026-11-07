# Source code

All source code of the pipeline lives under `src/`:

- `src/opls_ff_gen/` is the pipeline itself, an installable Python package.
  It is run as `opls-ff-gen` or `python -m opls_ff_gen`. Its data tables are
  in `src/opls_ff_gen/data/`.
- `src/scripts/` holds the tools around it. They build the data tables,
  run batches, and check or benchmark the pipeline. The pipeline never
  imports them at runtime.

The tests are in `tests/` at the top level. Install the package in editable
mode once, so that the tests and scripts import it from `src/`:

    pip install --no-deps -e .
    python -m unittest discover -s tests

## Pipeline package: `src/opls_ff_gen/`

The order below follows a run, from input to topology.

### Interfaces

- **`__init__.py`**: the package's public API: `generate_gromacs_topology`,
  `assign_molecule`, `assign_bonded_interactions`, `load_sdf`,
  `molecule_from_smiles`, their result and interaction record classes, and
  `AssignmentError`.
- **`__main__.py`**: the command-line interface (`opls-ff-gen`). It reads a
  SMILES string or a structure file and the options (`--qm none|xtb|dft`,
  `--optimize`, `--torsion-substitution`, `--qm-bonded-table`, DFT level,
  MM engine). It writes the GROMACS topology and prints a count of
  parameters per rule.
- **`workflow.py`**: the end-to-end driver. It runs the stages in order: atom
  typing and charges, then bonded assignment, then the optional run-time QM
  fit. It writes the ITP, TOP, GRO and JSON outputs; the ITP comments name
  every parameter that did not come straight from the OPLS/2020 table. With
  QM, the pre-fit topology goes to `qm_pending/`, and the final-named files
  are written only after the fit succeeds.

### Input structure

- **`structure.py`**:
  - Reads the input: SMILES, SDF/MOL, MOL2 or PDB. For PDB and MOL2, bond
    orders and charges can come from a template SMILES.
  - Optimizes the geometry if asked (`none`, `mmff` or `xtb`).
  - Applies the MMFF94 aromaticity model, so every input route types the
    same way.

### Atom typing

- **`opls_typer.py`**: OPLS-AA SMARTS atom typing in the Foyer dialect,
  without Foyer. It is a line-by-line port of Foyer 1.2.0's algorithm on top
  of RDKit: rule iteration, white- and blacklists, overrides.
- **`assignment.py`**: the atom-typing stage, and the entry point for
  charges.
  - Typing: OPLS-AA SMARTS rules plus the atom-type archive (definitions,
    precedence, borrowed types), then refinement to OPLS/2020 type numbers.
    The archive's refinement rules cover Si, the halogens, allenes and the
    newer alkane types.
  - Collisions: if several types match an atom and no precedence entry
    decides, the highest-numbered type wins.
  - Parameters: classes and LJ come from the OPLS/2020 table. An atom no
    OPLS rule can type takes the OPLS/2020 analog of its GAFF2 type if there
    is one, else GAFF2 LJ.
  - Charges: CM1A charges (scaled 1.14 and LBCC-corrected for neutral
    molecules). Writes the assignment JSON, CSV and SDF.
- **`smarts_typing.py`**: a generic table of ordered RDKit SMARTS rules,
  where the first rule that matches an atom types it. Used by the GAFF2
  typers.
- **`amber.py`**: GAFF2 atom types and Lennard-Jones parameters for atoms
  that OPLS cannot type. It uses the built-in rule table
  (`data/gaff2_lj.json`); no AmberTools is needed at runtime. Also the
  cosmetic `opls_8NN` names under which GAFF2 types are written to the ITP
  (`data/gaff2_display_types.json`; labels, not OPLS/2020 types).
- **`gaff2_analogs.py`**: OPLS/2020 analog types for GAFF2-typed atoms
  (`data/gaff2_opls_analogs.json`): a GAFF2 type plus SMARTS context mapped to
  an OPLS/2020 type of a class with the same element, hybridization and bond
  order. Atoms without an analog keep GAFF2. For bonded terms the OPLS/2020
  table lacks, an analog atom keeps its GAFF2 type (`analog_gaff2_classes`,
  passed to `assign_bonded_interactions` as `gaff2_classes`).

### Charges

- **`cm1a.py`**: CM1A partial charges, computed from a MOPAC AM1
  calculation with our own implementation of the CM1A charge model. This
  replaces AMSOL.
- **`lbcc.py`**: the published localized bond-charge corrections (LBCC) for
  1.14×CM1A-LBCC, and finding where they apply in a molecule.

### Parameter tables

- **`opls2020.py`**: parser for the OPLS/2020 SI tables in BOSS format
  (`oplsaa2020.par`, `oplsaa2020.sb`). Provides atom types, bonds, angles,
  torsion lookup (exact, `?` patterns, X-B-C-X), impropers, and the
  conversion from Fourier to Ryckaert-Bellemans coefficients. The SI values
  are never changed.
- **`gaff2_bonded.py`**: types every atom with its full GAFF2 type (a second
  SMARTS table; conjugated type pairs are assigned by alternation along
  Kekulé bonds). Looks up gaff2.dat bonds and angles for terms with a
  GAFF2-typed atom, and gaff2.dat general X-b-c-X torsions for the torsion
  fallback. The data are in `data/gaff2_bonded.json`.
- **`gaff2_torsions.py`**: GAFF2-derived general X-B-C-X torsions for OPLS
  central class pairs that OPLS/2020 has no general term for
  (`data/gaff2_general_torsions.par`). It also flags general terms that need
  a QM check.
- **`parmchk.py`**: parmchk2-style substitution of proper torsions missing
  from OPLS/2020. It uses AmberTools' `PARMCHK.DAT` scores for equal and
  corresponding types, and gives each substitution a penalty.
- **`qm_bonded.py`**: bonds and angles recorded from earlier QM fits
  (`data/qm_bonded_terms.par`). Angles are keyed by classes and small-ring
  context. It is used only when no run-time QM fit is requested.

### Bonded assignment

- **`bonded.py`**: enumerates bonds, angles, propers, impropers and 1-4
  pairs, and assigns every bonded parameter with a `rule` and a source line.
  The OPLS/2020 table is always the first choice.
  - With `--qm xtb|dft`, every term the table lacks becomes a `qm-pending`
    placeholder for the QM fit, and nothing else is consulted.
  - Without QM (the default), a missing term goes through the GAFF2 and
    recorded-QM tables, then the substitutions and fallbacks, so every term
    gets a value.
  - Impropers use the OPLS improper families.
- **`bonded_substitution.py`**: the fallback used when no table has a term
  (without QM).
  - Bonds: b0 from a parmchk2-style substitute (GAFF2 types if a GAFF2-typed
    atom is involved, then OPLS classes) or a default-class substitute, with
    the empirical k; otherwise both b0 and k are empirical.
  - Angles: empirical; θ0 is the mean angle about the centre's class, and k
    follows the GAFF empirical formula.
  - Torsions, in order:
    1. `default-general`: X-B-C-X with default central classes.
    2. `gaff2-dat-general`: gaff2.dat X-b-c-X for the full GAFF2 types.
    3. `central-pair`: the closest OPLS/2020 entry about the same central
       classes.
    4. `zero-torsion`: zero, as a last resort.

### Run-time QM fit (optional)

- **`qm_fit.py`**: the on-the-fly QM fit of `qm-pending` terms.
  - Optimizes the geometry and computes the Hessian with GFN2-xTB, or at the
    chosen DFT level.
  - Bonds and angles come from the modified Seminario method.
  - Torsions come from a joint damped ridge fit to relaxed scans.
  - Results are cached per molecule.
- **`qm.py`**: the QM side of the fit: DFT levels (PySCF), single-point
  energies, and relaxed GFN2-xTB dihedral scans (rotatable and ring bonds).
- **`xtb.py`**: runs the GFN2-xTB executable and reads its geometries and
  outputs.
- **`mm.py`**: the restrained MM relaxations along each scan for the
  torsion fit, run with OpenMM (default) or GROMACS.
- **`gromacs.py`**: finds the GROMACS executable (through `tools.py`).
- **`datafiles.py`**: integrity of the package data: `data/MANIFEST.json`
  (size and SHA-256 of every shipped data file, the user-installed SI
  excluded), checked before each run and by `opls-ff-gen doctor`.
- **`si.py`**: the user-installed OPLS/2020 SI files: where they are read
  from (`OPLS_FF_GEN_SI_DIR`, the per-user directory, the package's
  `data/opls2020/`), their checksums, and `opls-ff-gen install-si`.
- **`tools.py`**: external programs (MOPAC, xtb, GROMACS) and optional
  Python packages: lookup (environment variable, user config file, `PATH`),
  `opls-ff-gen config`, `opls-ff-gen doctor`, the pre-flight check of a run,
  and the default MM engine (OpenMM if importable, else GROMACS).
- **`context.py`**: chemical-context keys (classes plus SMARTS context) for
  bonded terms. The torsion fit uses them to share parameters among
  equivalent quartets.

## Tools: `src/scripts/`

Run them from the repository root, for example
`python src/scripts/run_batch_opls2020.py --help`.

### Running batches

- **`run_batch_opls2020.py`**: runs the pipeline on CSVs of molecules (with
  `--qm none|xtb|dft` and `--optimize`), then checks each topology with
  `gmx grompp` and a zero-step `gmx mdrun`.
- **`summarize_opls2020_batch.py`**: summarizes a batch: counts, bonded
  rules and QM fit quality.
- **`coverage_opls2020.py`**: for a molecule CSV, reports OPLS/2020 typing
  coverage and the bonded terms missing from the table, with the scans a QM
  fit would need.
- **`estimate_qm_cost.py`**: estimates the compute time of the run-time QM
  fit for a batch, at the xTB or a DFT level.

### Building pipeline data

- **`build_gaff2_table.py`**: builds `data/gaff2_lj.json`, the GAFF2 typing
  rules and LJ parameters for atoms OPLS cannot type.
- **`build_gaff2_bonded.py`**: builds `data/gaff2_bonded.json`: the full
  GAFF2 typing rules, plus gaff2.dat bonds, angles and general torsions
  (`--gaff2-dat PATH`).
- **`gaff2_general_torsions.py`**: derives the GAFF2 general X-B-C-X terms
  for OPLS central class pairs and writes `data/gaff2_general_torsions.par`.
  The `KEPT` entries override this rebuild.
- **`qm_bonded_terms.py`**: collects QM-fitted bonds and angles from batch
  outputs into `data/qm_bonded_terms.par`, with their provenance in
  `data/qm_bonded_terms.json`.
- **`export_foyer_rules.py`**: one-time export of Foyer 1.2.0's OPLS-AA
  typing rules to `data/opls_rules.json`. Needs Foyer.
- **`gaff2_reference.py`**: antechamber GAFF2 types for a molecule set,
  used as the reference when building and checking the GAFF2 rules. Needs
  AmberTools; never run by the pipeline.

### Equivalence checks against replaced tools

- **`check_foyer_equivalence.py`**: checks `opls_typer.py` against Foyer
  1.2.0, pattern graphs and per-atom typing.
- **`check_cm1a_equivalence.py`**: checks MOPAC AM1 + `cm1a.py` charges
  against AMSOL 7.1.
- **`check_mm_engines.py`**: checks the GROMACS MM engine against OpenMM
  (single-point energies and forces, relaxed scan profiles).
- **`compare_gaff2_rules.py`**: compares the GAFF2 LJ rule table with the
  antechamber reference types, by LJ.
- **`compare_gaff2_bonded_types.py`**: compares the full GAFF2 bonded types,
  and the resulting gaff2.dat terms, with antechamber.

### Benchmarks and batch comparisons

- **`benchmark_ligpargen.py`**: compares generated ITP files with LigParGen
  ITP files, atom by atom and term by term.
- **`benchmark_table_torsions.py`**: compares OPLS/2020 table torsions with
  relaxed xTB scans; this is the reference error for fitted or substituted
  torsions.
- **`benchmark_geometry_optimization.py`**: compares MMFF94 and GFN2-xTB
  geometry optimization (time, RMSD, effect on CM1A charges).
- **`compare_general_torsions.py`**: compares GAFF2 general and parmchk2
  torsions with the QM-fitted torsions on the same scans.
- **`compare_batches.py`**: file-by-file identity check of two batch output
  directories.
- **`compare_batches_numeric.py`**: numerical comparison of two batches
  (types, charges, bonded parameters, fit RMSE).
- **`check_dist.py`**: checks a built wheel and sdist before release
  (package data and checksums, licences, no SI or local-only files,
  version = pyproject = tag); used by the release CI.
- **`release_smoke_test.py`**: end-to-end check of an installed
  `opls-ff-gen` (SMILES, SDF, cation, zwitterion; `--qm` adds an xTB fit);
  used by the release CI.
- **`build_data_manifest.py`**: rewrites `src/opls_ff_gen/data/MANIFEST.json`
  after a data file is edited (the scripts that rebuild data files call it
  themselves; `tests/test_datafiles.py` fails while it is stale).
- **`reorder_to_itp.py`**: writes an SDF of a SMILES in the atom order of a
  LigParGen ITP, by graph isomorphism.

## Tests: `tests/`

- **`test_assignment.py`**: regression tests for SMARTS atom typing,
  borrowed types and CM1A-LBCC charges.
- **`test_opls2020.py`**: covers:
  - the OPLS/2020 tables and typing (precedence, collisions, Si, halogens);
  - the bonded lookup order, the GAFF2 and parmchk2 lookups, and the
    bond/angle/torsion fallbacks;
  - QM mode using the table only;
  - the recorded QM table, the Seminario method and the on-the-fly fit;
  - DFT levels.
- **`test_gaff2.py`**: the GAFF2 LJ and bonded rule tables against stored
  antechamber types.
- **`test_structure.py`**: structure input (SMILES, SDF, MOL2, PDB) and
  geometry optimization.
- **`test_context.py`**: chemical-context keys of bonded terms.
- **`test_fit_math.py`**: torsion conversions and the ridge solve of the
  torsion fit.
- **`test_workflow.py`**: end-to-end topology generation, including GROMACS
  checks.
