"""External programs and optional Python packages: lookup, status, pre-flight.

Programs (``PROGRAMS``) are found, in this order:

1. the environment variable (e.g. ``OPLS_FF_GEN_XTB``);
2. the user config file, ``[programs]`` table (``opls-ff-gen config set xtb PATH``);
   the file is ``OPLS_FF_GEN_CONFIG``, else
   ``$XDG_CONFIG_HOME/opls-ff-gen/config.toml`` (default ``~/.config/...``);
3. ``PATH``.

No path is hard-coded. A program the user already has is used wherever it is
installed; one that is missing is reported with its install command, before a
run starts (``preflight``) and by ``opls-ff-gen doctor``.

Optional Python packages (OpenMM, PySCF, geomeTRIC) must be importable by the
Python that runs opls-ff-gen: a package in another environment cannot be
loaded safely (compiled extensions are tied to their interpreter), so the
fix is to install it into this environment, or to install opls-ff-gen into
the environment that already has it.
"""

from __future__ import annotations

import importlib.util
import os
import re
import shutil
import subprocess
import tomllib
from dataclasses import dataclass
from pathlib import Path

CONFIG_ENV = "OPLS_FF_GEN_CONFIG"


@dataclass(frozen=True)
class Program:
    key: str  # config key
    label: str
    env: str  # environment variable
    names: tuple[str, ...]  # executable names on PATH
    required: bool
    used_for: str
    install: str
    version_args: tuple[str, ...]
    version_pattern: str


@dataclass(frozen=True)
class Package:
    module: str
    label: str
    used_for: str
    install: str


PROGRAMS = {
    "mopac": Program(
        "mopac", "MOPAC", "OPLS_FF_GEN_MOPAC", ("mopac",), True,
        "the CM1A charges (AM1; every run)",
        "conda install -c conda-forge mopac=23.2.5", ("--version",),
        r"MOPAC version (\S+)"),
    "xtb": Program(
        "xtb", "xtb (GFN2-xTB)", "OPLS_FF_GEN_XTB", ("xtb",), False,
        "SMILES geometry (MMFF94 if missing), --qm xtb, scans for --qm dft",
        "conda install -c conda-forge xtb=6.7", ("--version",),
        r"xtb version (\S+)"),
    "gmx": Program(
        "gmx", "GROMACS", "OPLS_FF_GEN_GMX", ("gmx", "gmx_mpi"), False,
        "--mm-engine gromacs; GROMACS checks of the batch scripts",
        "conda install -c conda-forge gromacs (or an existing GROMACS build)",
        ("--version",), r"GROMACS version:\s*(\S+)"),
}

PACKAGES = {
    "openmm": Package("openmm", "OpenMM", "--mm-engine openmm (torsion fits of --qm)",
                      "conda install -c conda-forge openmm  (or: pip install openmm)"),
    "pyscf": Package("pyscf", "PySCF", "--qm dft", "pip install pyscf"),
    "geometric": Package("geometric", "geomeTRIC", "--qm dft geometry optimization",
                         "pip install geometric"),
}

# extra settings the config file may hold besides the programs
SETTINGS = {
    "xtbpath": ("OPLS_FF_GEN_XTBPATH", "directory of the xtb parameter files "
                "(default: share/xtb next to the executable)"),
}


# -- config file --------------------------------------------------------------
def config_path() -> Path:
    if os.environ.get(CONFIG_ENV):
        return Path(os.environ[CONFIG_ENV]).expanduser()
    base = os.environ.get("XDG_CONFIG_HOME") or Path.home() / ".config"
    return Path(base) / "opls-ff-gen" / "config.toml"


def read_config() -> dict[str, str]:
    """The ``[programs]`` table of the config file ({} if there is none)."""
    path = config_path()
    if not path.is_file():
        return {}
    try:
        data = tomllib.loads(path.read_text())
    except (OSError, tomllib.TOMLDecodeError) as exc:
        raise RuntimeError(f"cannot read the opls-ff-gen config file {path}: {exc}") from exc
    table = data.get("programs", {})
    return {str(key): str(value) for key, value in table.items()}


def write_config(values: dict[str, str]) -> Path:
    path = config_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["# opls-ff-gen settings (edit, or use: opls-ff-gen config set KEY PATH)",
             "[programs]"]
    for key in sorted(values):
        escaped = values[key].replace("\\", "\\\\").replace('"', '\\"')
        lines.append(f'{key} = "{escaped}"')
    path.write_text("\n".join(lines) + "\n")
    return path


def set_config(key: str, value: str | None) -> Path:
    """Set (``value``) or remove (None) one config entry, after validating it."""
    if key not in PROGRAMS and key not in SETTINGS:
        raise ValueError(f"unknown key {key!r}; choose from "
                         f"{', '.join([*PROGRAMS, *SETTINGS])}")
    values = read_config()
    if value is None:
        values.pop(key, None)
        return write_config(values)
    path = Path(value).expanduser().resolve()
    if key in PROGRAMS and not _executable(path):
        raise ValueError(f"{path} is not an executable file")
    if key in SETTINGS and not path.is_dir():
        raise ValueError(f"{path} is not a directory")
    values[key] = str(path)
    return write_config(values)


# -- lookup -------------------------------------------------------------------
def _executable(path: Path) -> bool:
    return path.is_file() and os.access(path, os.X_OK)


@dataclass(frozen=True)
class Location:
    path: Path | None
    source: str  # "env OPLS_FF_GEN_XTB", "config", "PATH", or "not found"

    @property
    def ok(self) -> bool:
        return self.path is not None and _executable(self.path)


def locate(key: str) -> Location:
    program = PROGRAMS[key]
    if os.environ.get(program.env):
        return Location(Path(os.environ[program.env]).expanduser(), f"env {program.env}")
    configured = read_config().get(key)
    if configured:
        return Location(Path(configured).expanduser(), f"config {config_path()}")
    for name in program.names:
        found = shutil.which(name)
        if found:
            return Location(Path(found), "PATH")
    return Location(None, "not found")


def program_path(key: str) -> Path:
    """Path of a program, or its bare name if it is not found (running it then
    fails with a 'not found' error)."""
    return locate(key).path or Path(PROGRAMS[key].names[0])


def setting(key: str) -> str | None:
    env, _ = SETTINGS[key]
    return os.environ.get(env) or read_config().get(key)


def program_version(key: str, path: Path | None = None) -> str | None:
    program = PROGRAMS[key]
    path = path or locate(key).path
    if path is None or not _executable(path):
        return None
    env = dict(os.environ)
    if key == "xtb":
        env["XTBPATH"] = setting("xtbpath") or str(path.resolve().parent.parent / "share" / "xtb")
    try:
        result = subprocess.run([str(path), *program.version_args], capture_output=True,
                                text=True, timeout=30, check=False, env=env,
                                stdin=subprocess.DEVNULL)
    except (OSError, subprocess.TimeoutExpired):
        return None
    match = re.search(program.version_pattern, result.stdout + result.stderr)
    return match.group(1) if match else "unknown"


def package_available(module: str) -> bool:
    return importlib.util.find_spec(module) is not None


def package_version(module: str) -> str | None:
    from importlib import metadata

    if not package_available(module):
        return None
    try:
        return metadata.version(module)
    except metadata.PackageNotFoundError:
        return "unknown"


# -- tested versions ------------------------------------------------------------
# Versions opls-ff-gen was validated with. Others are expected to work (the
# interfaces used are stable), so a different version only gives a note.
TESTED = {
    "mopac": ("23",),  # 23.2.5
    "xtb": ("6", "7"),  # 6.7.1
}
TESTED_TEXT = {"mopac": "23.2.5", "xtb": "6.7.1", "rdkit": "2025.09-2026.09"}


def _parts(version: str) -> tuple[str, ...]:
    return tuple(re.findall(r"\d+", version))


def version_note(key: str, version: str | None) -> str | None:
    """None for a tested version (or a program without a tested version, such
    as GROMACS), else a note (never an error)."""
    if version is None or key not in TESTED_TEXT:
        return None
    label = PROGRAMS[key].label if key in PROGRAMS else "RDKit"
    if version == "unknown":
        return f"{label}: the version could not be read; tested with {TESTED_TEXT[key]}"
    parts = _parts(version)
    if key == "rdkit":
        tested = len(parts) >= 2 and (2025, 9) <= (int(parts[0]), int(parts[1])) <= (2026, 9)
    else:
        tested = parts[:len(TESTED[key])] == TESTED[key]
    if tested:
        return None
    return (f"{label} {version} differs from the tested {TESTED_TEXT[key]}; it is expected "
            "to work, but results are not validated with it")


def run_version_notes(*, xtb_used: bool) -> list[str]:
    """Version notes for the programs and packages a run uses."""
    notes = []
    for key in ("mopac", "xtb") if xtb_used else ("mopac",):
        location = locate(key)
        if location.ok:
            note = version_note(key, program_version(key, location.path))
            if note:
                notes.append(note)
    note = version_note("rdkit", package_version("rdkit"))
    return notes + ([note] if note else [])


# -- pre-flight ---------------------------------------------------------------
def _program_hints(key: str, indent: str) -> str:
    program = PROGRAMS[key]
    return (f"{indent}If you have it: opls-ff-gen config set {key} /path/to/{program.names[0]}"
            f"  (or set {program.env})\n"
            f"{indent}To install it:  {program.install}")


def _program_problem(key: str) -> str:
    program, location = PROGRAMS[key], locate(key)
    where = (f"{location.path} ({location.source}) is not an executable file"
             if location.path is not None else "not found")
    return (f"{program.label}: {where}; needed for {program.used_for}.\n"
            + _program_hints(key, "    "))


def _package_problem(module: str) -> str:
    package = PACKAGES[module]
    return (f"{package.label}: not importable by this Python; needed for {package.used_for}.\n"
            f"    Install it into this environment: {package.install}\n"
            f"    (a copy in another environment cannot be used; alternatively install "
            f"opls-ff-gen into that environment)")


def preflight(*, qm: str | None, mm_engine: str, optimize_xtb: bool) -> tuple[list[str], list[str]]:
    """(errors, warnings) for a run with these options, before any work starts."""
    from . import si
    from .datafiles import verify

    errors, warnings = verify()
    si_errors, si_warnings = si.check()
    errors += si_errors
    warnings += si_warnings
    if not locate("mopac").ok:
        errors.append(_program_problem("mopac"))
    xtb_ok = locate("xtb").ok
    if qm in ("xtb", "dft") and not xtb_ok:
        errors.append(_program_problem("xtb"))
    elif optimize_xtb and not xtb_ok:
        warnings.append("xtb not found: the geometry is optimized with MMFF94 instead.\n"
                        + _program_problem("xtb"))
    if qm in ("xtb", "dft"):
        if mm_engine == "openmm" and not package_available("openmm"):
            errors.append(_package_problem("openmm"))
        if mm_engine == "gromacs" and not locate("gmx").ok:
            errors.append(_program_problem("gmx"))
    if qm == "dft":
        errors += [_package_problem(m) for m in ("pyscf", "geometric")
                   if not package_available(m)]
    return errors, warnings


def default_mm_engine() -> str:
    """``OPLS_FF_GEN_MM_ENGINE``, else OpenMM if importable, else GROMACS if
    found, else openmm (reported as missing by the pre-flight check)."""
    configured = os.environ.get("OPLS_FF_GEN_MM_ENGINE")
    if configured:
        return configured
    if package_available("openmm"):
        return "openmm"
    if locate("gmx").ok:
        return "gromacs"
    return "openmm"


# -- doctor ---------------------------------------------------------------------
def _tag(ok: bool, required: bool = False) -> str:
    return f"{'[ok]' if ok else '[MISSING]' if required else '[missing]':<10}"


def doctor() -> tuple[str, bool]:
    """Status report of every dependency; True if the required ones are present."""
    from . import si

    lines = [f"config file: {config_path()}"
             + ("" if config_path().is_file() else " (none)"), "", "Programs:"]
    required_ok = True
    for key, program in PROGRAMS.items():
        location = locate(key)
        tag = "required" if program.required else "optional"
        if location.ok:
            version = program_version(key, location.path)
            lines.append(f"  {_tag(True)}{program.label} {version or ''} ({tag})\n"
                         f"            {location.path}  [{location.source}]")
            note = version_note(key, version)
            if note:
                lines.append(f"            note: {note}")
        else:
            required_ok &= not program.required
            broken = (f"; {location.path} [{location.source}] is not executable"
                      if location.path is not None else "")
            lines.append(f"  {_tag(False, program.required)}{program.label} ({tag}; "
                         f"{program.used_for}){broken}\n" + _program_hints(key, " " * 12))
    lines += ["", "Python packages (this interpreter):"]
    for module, package in PACKAGES.items():
        version = package_version(module)
        if version:
            lines.append(f"  {_tag(True)}{package.label} {version} (optional)")
        else:
            lines.append(f"  {_tag(False)}{package.label} (optional; {package.used_for}): "
                         f"{package.install}")
    for module, label in (("rdkit", "RDKit"), ("numpy", "NumPy")):
        version = package_version(module)
        required_ok &= version is not None
        lines.append(f"  {_tag(version is not None, True)}{label} {version or ''} (required)")
        note = version_note(module, version) if module == "rdkit" else None
        if note:
            lines.append(f"            note: {note}")
    from .datafiles import DATA_DIR, MANIFEST, verify

    data_errors, data_warnings = verify()
    required_ok &= not data_errors
    lines += ["", "Package data:"]
    if not data_errors and not data_warnings:
        import json

        count = len(json.loads(MANIFEST.read_text())["files"])
        lines.append(f"  {_tag(True)}{count} files match {MANIFEST.name} ({DATA_DIR})")
    lines += [f"  {_tag(False, True)}{message}" for message in data_errors]
    lines += [f"  [warning] {message}" for message in data_warnings]
    directory, source = si.si_dir()
    si_errors, si_warnings = si.check()
    required_ok &= not si_errors
    lines += ["", f"OPLS/2020 SI (user-installed): {directory} [{source}]"]
    if si_errors:
        lines += [f"  {_tag(False, True)}" + message.replace("\n", "\n  ")
                  for message in si_errors]
    else:
        for si_file in si.SI_FILES:
            lines.append(f"  {_tag(True)}{si_file.name} ({si_file.original})")
        lines += [f"  [warning] {message}" for message in si_warnings]
        if not si_warnings:
            lines.append("            checksums match the Oct 1, 2023 release")
    engine = default_mm_engine()
    lines += ["", f"default MM engine for --qm fits: {engine}"]
    lines.append("" if required_ok else
                 "\nRequired dependencies are missing: runs will stop until they are installed.")
    return "\n".join(lines).rstrip() + "\n", required_ok
