# Release and distribution plan (2026-11-07)

Goal: users download a GitHub release, install it in their own environment
and get the same results as here, with as few external dependencies as
possible. The install must fail early and clearly, never halfway through a
run. Status: plan only; nothing below is implemented yet, except where
marked.

## Decisions

- **OPLS/2020 SI is not bundled (user, 2026-11-07).** The user has no
  permission from ACS to redistribute the Supporting Information of
  doi:10.1021/acs.jpcb.3c06602. Users download it themselves and install
  it with a command that checks its checksums (step 1). Re-encoding the SI
  values as a Python "library" shipped with the code was considered and
  rejected (user, 2026-11-07): it distributes the same content in another
  format (compilation copyright, EU database right; the Crossref licence
  record lists only ACS's standard policies, no CC licence). The user keeps
  the user-install plan.
- **MOPAC is a required, user-installed program (user, 2026-11-07).**
  It is not bundled in the wheel; users install it themselves (conda-forge
  `mopac`, pinned to 23.2.5, the version behind all validation data), and
  the release documents and checks it. MOPAC is Apache-2.0, so bundling
  would be allowed, but the binary needs its GCC runtime and BLAS libraries
  and would make the wheels platform-specific.
- **The GNN charge model is dropped (user, 2026-11-07).** Charges stay
  MOPAC AM1 + CM1A; no neural network replaces MOPAC (step 9).
- **A missing xtb falls back to MMFF94** for the geometry optimization
  (user, 2026-11-07; done, `structure.optimize_geometry`, warning and
  `Geometry_Fallback`). `--qm xtb` without xtb cannot fall back and stays an
  error (to be reported by the pre-flight check, step 4).
- Still open (my recommendations in brackets):
  - platforms [Linux, then macOS; Windows depends on MOPAC/xtb];
  - Python version: the environment files pin python=3.12 (the only tested
    version; conda would otherwise pick 3.15); `requires-python >= 3.12`.
- **MIT license** for the code (user, 2026-11-07): `LICENSE`; `NOTICE` gives
  the source and terms of every data file.

## Distribution channels (2026-11-07)

- **GitHub release** (the plan): the wheel, the sdist, `environment.yml` and
  `environment-full.yml`, SHA-256 sums, and GitHub's source archive. Users
  install the wheel with pip into a conda or pip environment, or run the
  source directly. opls-ff-gen itself is **not** published on PyPI or
  conda-forge; only its dependencies come from there. **Not published on
  PyPI or conda-forge for now (user, 2026-11-07).** Either remains possible
  later: PyPI is simple because the wheel already exists; conda-forge needs a
  feedstock and a review (MOPAC and xtb are already there).
- **Run from source without installing** (user, 2026-11-07):
  `run_opls_ff_gen.py` at the repository root. It runs with any Python 3.12+
  that has RDKit and NumPy, uses the source tree even if a copy is
  installed, and takes the same arguments and the same program lookup as the
  installed command. Verified with an existing environment (`chem`: RDKit
  2025.09.6, no opls-ff-gen installed): config set mopac/xtb, doctor, then a
  run. The bonded terms were identical to the installed run; the charges
  differed by about 1e-6 e (another RDKit version gives a slightly
  different MMFF94 geometry). A Python without RDKit gets a clear message.

## Dependency tiers

| tier | what | needed for | source |
|---|---|---|---|
| required, Python | Python >= 3.12, RDKit >= 2025.9, NumPy | everything | pip / conda |
| required, program | MOPAC 23.2.5 (user install) | CM1A charges | conda-forge `mopac` (Apache-2.0) |
| required, data | package data; OPLS/2020 SI (user-installed) | everything | wheel; `opls-ff-gen install-si` |
| optional | xtb | SMILES geometry, `--qm xtb`, DFT scans | conda-forge `xtb` |
| optional | OpenMM (extra `[qm]`) | torsion fits (default MM engine) | pip / conda |
| optional | GROMACS | `--mm-engine gromacs`, grompp checks | user install |
| optional | PySCF + geomeTRIC (extra `[dft]`) | `--qm dft` | pip |

## State found (2026-11-07; the first three points fixed, steps 2-3)

- The wheel built from `pyproject.toml` has none of the 11 files in
  `src/opls_ff_gen/data/` (no `package-data` entry). It works here only
  because of the editable install.
- `openmm` is a required dependency, but it is imported only in `mm.py` for
  the torsion fits; `--qm none` (the default) does not need it.
- xtb and GROMACS fall back to hard-coded paths under `/home/jiang/software/`
  (`xtb.py`, `gromacs.py`).
- The SMILES default is the xtb optimization, so an install without xtb
  fails on the most common input.
- There is no `LICENSE` file.

## Steps

1. **SI and data licensing.** `opls-ff-gen install-si FILE FILE` (done
   2026-11-07, `si.py`): identifies the two files by checksum or header (any
   order), checks the SHA-256 of `jp3c06602_si_002.txt` (.par) and
   `jp3c06602_si_003.txt` (.sb) against the Oct 1, 2023 release (other
   versions only with `--force`; runs then warn), and copies them to the
   per-user directory `$XDG_DATA_HOME/opls-ff-gen/opls2020` (default
   `~/.local/share/...`), which survives upgrades. Lookup:
   `OPLS_FF_GEN_SI_DIR` if set (authoritative), else the user directory,
   else the package's `data/opls2020/` (development checkouts). Verified in
   the clean venv: install-si, then an ITP identical to the checkout.
   `LICENSE` (MIT) and `NOTICE` added 2026-11-07. **AmberTools terms
   (checked 2026-11-07** in AmberTools/LICENSE of ambertools26_rc7.tar.bz2,
   SHA-256 as conda-forge): "the force field parameter files in the dat/leap
   subdirectory have been placed in the public domain", so the GAFF2-derived
   files (gaff2_lj.json, gaff2_bonded.json, gaff2_general_torsions.par;
   gaff2.dat is dat/leap/parm/) are public domain. `PARMCHK.DAT` and
   `PARM_BLBA_GAFF2.DAT` are in dat/antechamber/, which no exception covers,
   so they are GPL v3 (the AmberTools default). They are shipped unchanged
   with the GPL text (data/parmchk/AmberTools-LICENSE.txt, README.txt); the
   package licence expression is `MIT AND GPL-3.0-only`. Open (user):
   parmchk.py/bonded_substitution.py reimplement parmchk2's steps "following
   parmchk2.c". A behavioural reimplementation is normally independent work;
   a line-by-line translation of the C code would be a derived work of GPL
   code. Only the author of that port can say which it is. The alternative
   to shipping the GPL files is to read them from the user's AmberTools
   installation. Clean-up done: the unused GROMACS sigma/epsilon values
   (LGPL-2.1) were removed from opls_rules.json (now Foyer MIT only).
2. **Packaging (done 2026-11-07).** `[tool.setuptools.package-data]` plus `MANIFEST.in`, so
   the wheel and sdist contain `data/**` (SI excluded). Leave out the
   local-only material (benchmarks, skills, validation outputs). Move OpenMM
   to an extra `[qm]` and keep `[dft]`. Add `data/MANIFEST.json` (file,
   SHA-256), checked at start-up: a missing or damaged file is a clear error.
   Result: the wheel holds all 11 data files + the manifest and no SI; the
   sdist holds the package, its data, README and pyproject (no SI, scripts,
   tests or benchmarks). The wheel installed into a clean venv (RDKit and
   NumPy from PyPI, MOPAC set with `config set`) gives an ITP identical to
   the development checkout; a deleted or edited data file stops the run.
   In a development checkout an edited data file is only a warning; run
   `python src/scripts/build_data_manifest.py` (the data-building scripts
   do it themselves; `tests/test_datafiles.py` fails while it is stale).
3. **Finding programs (done 2026-11-07, `tools.py`).** No hard-coded
   paths. Each program is looked up in its environment variable, then the
   config file `~/.config/opls-ff-gen/config.toml` (`opls-ff-gen config
   set/unset/show`), then `PATH`. OpenMM moved to the extra `[qm]`; without
   OpenMM the fits default to GROMACS if found. Version notes (done
   2026-11-07): MOPAC other than 23.x, xtb other than 6.7, or RDKit outside
   2025.09-2026.09 gives a NOTE (run) or `note:` (doctor), never an error
   (user: other versions will likely work). An xtb that is present but fails
   now falls back to MMFF94, like a missing one.
4. **Pre-flight checks (done 2026-11-07).** Before any work starts, the
   CLI checks what the chosen options need (SI files and MOPAC always; xtb
   for `--qm`; OpenMM or GROMACS for the fits; PySCF and geomeTRIC for DFT)
   and reports every missing item at once, with "if you have it" (config
   set / environment variable) and "to install it" commands; exit code 2.
   `opls-ff-gen doctor` lists each dependency (found/missing, version, path,
   source) and exits 1 if a required one is missing. Python packages must be
   in the running environment (another environment's copy cannot be loaded).
5. **Defaults when a program is missing.** No xtb: the geometry optimization
   uses MMFF94 with a warning and `Geometry_Fallback` in the JSON, the same
   mechanism as the bonding fallback (user decision and done, 2026-11-07).
   `--qm xtb` without xtb is an error.
6. **Install routes (done 2026-11-07).** Release assets `environment.yml` (python, rdkit,
   numpy, mopac) and `environment-full.yml` (adds xtb, openmm, pyscf,
   geometric), plus the wheel. The pip-only route is the wheel plus a MOPAC
   binary. `INSTALL.md` covers both routes, the SI step and `doctor`.
   Verified by creating both environments from conda-forge. The minimal one
   (RDKit 2026.9.2, MOPAC 23.2.5) plus the wheel, run as a new user (empty
   config and data directories), went through install-si, doctor (exit 0)
   and an MMFF94 run, with an ITP identical to the checkout. The full one
   adds xtb 6.7.1, OpenMM 8.6.1, PySCF 2.15 and geomeTRIC 1.1.1; there a
   `--qm xtb` run fitted 11 missing terms in 31 s.
7. **Release CI (written 2026-11-07, `.github/workflows/release.yml`; first
   GitHub run pending).** Triggers: tag `v*`, pull requests to main, manual.
   Jobs:
   - `build`: `python -m build` plus `src/scripts/check_dist.py` (all
     manifest data present with matching SHA-256, licences present, no SI or
     local-only files, version = pyproject = tag).
   - `test-conda`: minimal env; full env + conda-forge GROMACS. Each runs
     doctor, the unit tests against the installed wheel, and
     `src/scripts/release_smoke_test.py`: SMILES, SDF, cation, zwitterion,
     plus a `--qm xtb` fit in the full env.
   - `test-pip`: PyPI RDKit/NumPy, MOPAC from a separate conda prefix via
     `config set`, smoke test with MMFF94.
   - `release` (tags only): a DRAFT GitHub release with the wheel, sdist,
     environment files, INSTALL.md and SHA256SUMS.

   The SI comes from two repository secrets (xz + base64; 31 KB and 17 KB,
   under GitHub's 48 KB limit), verified by install-si and never printed.
   The user must add them before the first run (workflow header). Replayed
   locally with the same steps (scratch script): see PROGRESS.md.
8. **Reproducibility record.** The assignment JSON gets the package version,
   the data-manifest hash, the SI hashes, and the MOPAC/xtb/RDKit versions.
9. **MOPAC as a user install.** `environment.yml` pins `mopac=23.2.5`;
   `INSTALL.md` covers the conda-forge route and `OPLS_FF_GEN_MOPAC`;
   `doctor` and the pre-flight check report a missing or other-version
   MOPAC; the release notes carry the MOPAC Apache-2.0 attribution.
   (Replaces the GNN MOPAC bypass planned earlier the same day, dropped by
   the user.)

Order: steps 2-4 first (independent of licensing), then 1, 5, 6, then 7 and
8; step 9 goes with step 6.
