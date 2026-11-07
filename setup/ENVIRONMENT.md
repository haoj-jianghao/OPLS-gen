# Environment

The workflow runs in one conda environment (here `openMM`).

## Python packages

| Package | Needed for | Install |
|---|---|---|
| RDKit | molecules, SMARTS atom typing, GAFF2 fallback | conda-forge |
| NumPy | fitting | conda-forge |
| OpenMM | MM relaxations of the on-the-fly torsion fit with `--mm-engine openmm` (default); not needed with `--mm-engine gromacs` | conda-forge |
| PySCF, geomeTRIC | `--qm dft` only | pip |

Foyer, ParmEd, AmberTools and AMSOL are **not** needed:
- atom typing uses the project's own engine (`src/opls_ff_gen/opls_typer.py`) and rule data
  (`src/opls_ff_gen/data/opls_rules.json`, exported from Foyer 1.2.0);
- GAFF2 types for atoms OPLS cannot type come from `src/opls_ff_gen/data/gaff2_lj.json`;
- CM1A charges are computed by `src/opls_ff_gen/cm1a.py` from a MOPAC AM1 calculation.

Foyer (with ParmEd), AmberTools and AMSOL 7.1 are only needed to re-run the
equivalence checks against them:
- `src/scripts/check_foyer_equivalence.py` (Foyer 1.2.0);
- `src/scripts/gaff2_reference.py` (antechamber);
- `src/scripts/check_cm1a_equivalence.py` (AMSOL 7.1).

## External programs

| Program | Needed for | Path setting |
|---|---|---|
| MOPAC 23 (conda-forge `mopac`) | AM1 wavefunction for the CM1A charges (always) | `OPLS_FF_GEN_MOPAC`, config `mopac`, `PATH` |
| xtb 6.7 (GFN2-xTB) | SMILES geometry (default `--optimize xtb`) and the on-the-fly QM fit (`--qm xtb`, and scan geometries for `--qm dft`) | `OPLS_FF_GEN_XTB`, config `xtb`, `PATH` (parameter files: `OPLS_FF_GEN_XTBPATH` or config `xtbpath`) |
| GROMACS | `--mm-engine gromacs`, and checking topologies (`gmx grompp`/`mdrun` in `src/scripts/run_batch_opls2020.py`, the tests) | `OPLS_FF_GEN_GMX`, config `gmx`, `gmx`/`gmx_mpi` on `PATH` |

Every program is looked up in its environment variable, then the user config
file (`~/.config/opls-ff-gen/config.toml`, or `OPLS_FF_GEN_CONFIG`), then
`PATH` (`src/opls_ff_gen/tools.py`); no path is hard-coded. A program you
already have is registered with `opls-ff-gen config set KEY PATH`.
`opls-ff-gen doctor` lists every dependency (found or missing, version,
path, install command), and each run checks what its options need before it
starts. Optional Python packages (OpenMM, PySCF, geomeTRIC) must be
installed in the environment that runs opls-ff-gen; a copy in another
environment cannot be used. Without OpenMM, the `--qm` fits use GROMACS
if it is found.

## OPLS/2020 parameter files

These are not in git. See the README section "OPLS/2020 parameters".

## History

- `openMM-with-foyer-*.txt` and `chem-after-openmm-removal-explicit.txt` are
  environment snapshots from 2026-10-16, when Foyer 1.2.0 was installed from
  locally built wheels (the conda-forge Foyer 1.1.0 needed older NumPy and OpenMM).
- The Foyer installation notes were removed with the switch to the native
  typer (2026-10-18). The wheels stay in `setup/wheels/` (not in git) for the
  reference checks.
