# Installing opls-ff-gen

opls-ff-gen is pure Python, so nothing needs to be compiled. A run needs four
things:

| what | required? | how you get it |
|---|---|---|
| Python 3.12, RDKit (>= 2025.9), NumPy | required | conda (below) or pip |
| MOPAC 23.2.5 | required: AM1 for the CM1A charges | conda-forge `mopac`, or an existing installation |
| the OPLS/2020 parameter files | required | you download them (they may not be redistributed), then `opls-ff-gen install-si` |
| xtb 6.7, OpenMM, PySCF + geomeTRIC, GROMACS | optional | only for the features listed under [Optional features](#optional-features) |

There are two ways to use it:
- **install the wheel** from the GitHub release into a conda or pip
  environment (sections 1-5);
- **run the source directly**, without installing anything
  ([section 6](#6-run-from-the-source-without-installing)).

opls-ff-gen itself is not on PyPI or conda-forge. Its dependencies come from
conda-forge or PyPI; the package comes from the GitHub release.

Linux is tested. macOS should work through conda-forge but is untested so
far. On Windows, use WSL.

## 1. Install (conda, recommended)

Download from the GitHub release:
- the wheel `opls_ff_gen-<version>-py3-none-any.whl`;
- one of the environment files:
  - `environment.yml`: required dependencies only;
  - `environment-full.yml`: adds every optional feature except GROMACS.

```bash
conda env create -f environment.yml          # or environment-full.yml
conda activate opls-ff-gen                   # or opls-ff-gen-full
pip install opls_ff_gen-<version>-py3-none-any.whl
```

To install into an environment you already have, add the dependencies
yourself, then the wheel:

```bash
conda install -c conda-forge "rdkit>=2025.9" numpy mopac=23.2.5
pip install opls_ff_gen-<version>-py3-none-any.whl
```

### pip only

`pip install opls_ff_gen-<version>-py3-none-any.whl` brings RDKit and NumPy
from PyPI. MOPAC is not a Python package, so you install it separately:
- conda-forge `mopac`; or
- the binaries from https://github.com/openmopac/mopac/releases. With these,
  point opls-ff-gen at the program: `opls-ff-gen config set mopac /path/to/mopac`.

## 2. Install the OPLS/2020 parameter files (once per user)

The OPLS/2020 parameters are the Supporting Information of W. L. Jorgensen
et al., *J. Phys. Chem. B* (2024), https://doi.org/10.1021/acs.jpcb.3c06602.
Download the two text files `jp3c06602_si_002.txt` and
`jp3c06602_si_003.txt` from the article page, then:

```bash
opls-ff-gen install-si jp3c06602_si_002.txt jp3c06602_si_003.txt
```

- Each file is checked against the release opls-ff-gen was validated with
  (Oct 1, 2023).
- The files go to `~/.local/share/opls-ff-gen/opls2020/`, where they survive
  upgrades and reinstalls of the package.
- To keep them elsewhere, use `--dest DIR` and set `OPLS_FF_GEN_SI_DIR=DIR`.

## 3. Use programs you already have

opls-ff-gen looks for each program, in order:
1. the environment variable;
2. its config file `~/.config/opls-ff-gen/config.toml`;
3. `PATH`.

If you already have MOPAC, xtb or GROMACS somewhere else, register it:

```bash
opls-ff-gen config set xtb  /opt/xtb-6.7.1/bin/xtb     # keys: mopac, xtb, gmx, xtbpath
opls-ff-gen config set gmx  /opt/gromacs/bin/gmx
opls-ff-gen config show
```

The environment variables `OPLS_FF_GEN_MOPAC`, `OPLS_FF_GEN_XTB`,
`OPLS_FF_GEN_GMX` and `OPLS_FF_GEN_XTBPATH` take precedence over the config file.

Python packages (OpenMM, PySCF, geomeTRIC) are different: they have to be
installed in the environment that runs opls-ff-gen. A copy in another
environment cannot be used. If you already have an environment with OpenMM,
install opls-ff-gen into that environment instead.

## 4. Check the installation

```bash
opls-ff-gen doctor
```

`doctor` lists every dependency: found or missing, version, path, and where
the path came from (environment variable, config file or `PATH`). It also
checks the package data checksums, and the SI files and whether they match
the validated release. It gives the install command for anything missing,
and exits with 1 if a required item is missing.

Each run also checks what its options need before it starts. If something
is missing, it reports everything missing at once and stops (exit code 2),
before writing any output.

## 5. First run

```bash
opls-ff-gen --smiles 'CC(=O)Nc1ccc(O)cc1' --name paracetamol --output-dir out/paracetamol
```

This writes a GROMACS topology (`.itp`, `.top`, `.gro`) and JSON records to
`out/paracetamol/`. `opls-ff-gen --help` lists the options.

## 6. Run from the source without installing

Download the source code from the GitHub release (or clone the repository),
and run `run_opls_ff_gen.py` with any Python 3.12+ that can import RDKit and
NumPy. For example, use an environment you already have:

```bash
PY=/path/to/env/bin/python            # a Python with RDKit and NumPy
$PY run_opls_ff_gen.py doctor
$PY run_opls_ff_gen.py config set mopac /path/to/mopac    # programs you already have
$PY run_opls_ff_gen.py config set xtb   /path/to/xtb
$PY run_opls_ff_gen.py install-si jp3c06602_si_002.txt jp3c06602_si_003.txt
$PY run_opls_ff_gen.py --smiles 'CCO' --name ethanol --output-dir out/ethanol
```

- `run_opls_ff_gen.py` takes the same arguments as the `opls-ff-gen` command,
  and finds programs the same way (config file, environment variables,
  `PATH`).
- It uses the source tree even if a copy of the package is installed.
- The Python you run it with must itself provide RDKit and NumPy, plus
  OpenMM or PySCF for the features that need them. Python packages cannot be
  borrowed from another environment by a path. If the Python lacks them, the
  script says so and stops.
- Instead of `install-si`, you can copy the SI files into the source tree as
  `src/opls_ff_gen/data/opls2020/oplsaa2020.par` and `.sb`.

## Optional features

| feature | needs | conda-forge | notes |
|---|---|---|---|
| xTB geometry for SMILES (default) | xtb | `xtb=6.7` | without xtb, MMFF94 is used, with a warning; the JSON records it |
| `--qm xtb` on-the-fly fit of missing bonded terms | xtb, plus OpenMM or GROMACS | `xtb=6.7 openmm` | without OpenMM, GROMACS is used if found |
| `--qm dft` | the above, plus PySCF and geomeTRIC | `pyscf geometric` | |
| `--mm-engine gromacs` | GROMACS | `gromacs` (CPU build) or your own build | `opls-ff-gen config set gmx PATH` |

To add one to an existing environment: `conda install -c conda-forge xtb=6.7 openmm`
(or `pip install openmm pyscf geometric`). Then run `opls-ff-gen doctor`.

## Upgrade or remove

- **Upgrade:** `pip install --force-reinstall opls_ff_gen-<new>-py3-none-any.whl`.
  The SI files and the config file are kept.
- **Remove:** `pip uninstall opls-ff-gen`. Then delete
  `~/.local/share/opls-ff-gen/` (the SI files) and `~/.config/opls-ff-gen/`
  (the config) if you no longer need them.

## Troubleshooting

| message | fix |
|---|---|
| `OPLS/2020 SI not installed` | step 2 |
| `MOPAC: not found` | `conda install -c conda-forge mopac=23.2.5`, or `opls-ff-gen config set mopac PATH` |
| `package data incomplete` / `differs from its manifest` | the installation is damaged; reinstall the wheel |
| `SI files ... differ from the Oct 1, 2023 release` | warning only: files installed with `--force` (another SI version); results may differ from the validated ones |
| `xtb not found: the geometry is optimized with MMFF94 instead` | warning only; install xtb for the default geometry |
| `geometry fallback: xtb optimization failed` | warning only: xtb ran but failed (log in the output directory's `geometry/`); MMFF94 was used |
| `NOTE: ... differs from the tested ...` | information only: another MOPAC, xtb or RDKit version than the tested ones (MOPAC 23.2.5, xtb 6.7.1, RDKit 2025.09-2026.09); expected to work |

## Development install (from the git repository)

```bash
conda env create -f environment-full.yml && conda activate opls-ff-gen-full
pip install --no-deps --no-build-isolation -e .
python -m unittest discover -s tests
```

In a checkout, copies of the SI in `src/opls_ff_gen/data/opls2020/`
(`oplsaa2020.par`, `oplsaa2020.sb`) are used when no user installation
exists. After you edit a data file, run
`python src/scripts/build_data_manifest.py` to update the checksum manifest.
