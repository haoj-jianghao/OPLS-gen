"""Restrained MM relaxations for the torsion fit, with OpenMM or GROMACS.

``mm_profile`` relaxes the molecule's topology at every point of a QM scan:
starting from the QM geometry, the driven dihedral is restrained to its QM
value (harmonic, k = 5e4 kJ/mol/rad^2) and every other rotatable bond is held
near its QM value through one dihedral (k = 2e3); the returned energy is the
MM potential without these restraints. Both engines compute the same gas-phase
energy (all pairs, nrexcl 3, 1-4 pairs scaled by the topology's fudge factors):

``openmm``
    ``app.GromacsTopFile`` with ``NoCutoff``, restraints as
    ``CustomTorsionForce``, ``LocalEnergyMinimizer`` (tolerance 1e-3).
``gromacs``
    ``gmx grompp``/``mdrun`` per point: L-BFGS minimisation (emtol 1e-3) with
    ``[ dihedral_restraints ]``; vacuum through a periodic box at least twice
    a cut-off that exceeds the molecule by 1 nm, so every pair interacts once
    and no periodic image does (plain cut-off, i.e. reaction field with
    epsilon_rf = 1, and no potential modifiers: the bare Coulomb and LJ sums).

The engine is chosen with ``engine`` (default ``OPLS_FF_GEN_MM_ENGINE``, else
openmm if it is importable, else gromacs if it is found; ``tools.py``).
OpenMM's platform is ``OPLS_FF_GEN_MM_PLATFORM`` (default Reference).
"""

from __future__ import annotations

import math
import os
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from rdkit import Chem

from .qm import dihedral_deg
from .tools import default_mm_engine

ENGINES = ("openmm", "gromacs")
# OPLS_FF_GEN_MM_ENGINE, else OpenMM if importable, else GROMACS if found
DEFAULT_ENGINE = default_mm_engine()
# "Reference" reproduces earlier runs exactly and is as fast as "CPU" for one
# small molecule.
MM_PLATFORM = os.environ.get("OPLS_FF_GEN_MM_PLATFORM", "Reference")
K_DRIVE = 5.0e4  # kJ/mol/rad^2
K_LOCK = 2.0e3
EM_TOLERANCE = 1e-3  # kJ/mol/nm


def locked_dihedrals(molecule: Chem.Mol, bond,
                     lock_hydrogens: bool = False) -> list[tuple[int, int, int, int]]:
    """One heavy-atom dihedral per rotatable bond, except the driven bond.

    With ``lock_hydrogens`` a bond whose outer neighbours on one side are only
    hydrogens (OH, NH2, CH3 rotors) is locked through such a hydrogen, so that
    e.g. a hydroxyl cannot swing into a hydrogen bond the QM scan does not have.
    """
    pattern = Chem.MolFromSmarts("[!#1&!D1]-&!@[!#1&!D1]")
    result = []
    for j, k in molecule.GetSubstructMatches(pattern):
        if {j, k} == set(bond):
            continue

        def outer(center, other):
            neighbors = [a for a in molecule.GetAtomWithIdx(center).GetNeighbors()
                         if a.GetIdx() != other]
            heavy = [a.GetIdx() for a in neighbors if a.GetAtomicNum() > 1]
            if heavy:
                return heavy[0]
            return neighbors[0].GetIdx() if lock_hydrogens and neighbors else None

        i, l = outer(j, k), outer(k, j)
        if i is not None and l is not None:
            result.append((i, j, k, l))
    return result


def mm_profile(top: Path, points, drive, molecule: Chem.Mol, return_geometries: bool = False,
               lock_hydrogens: bool = False, engine: str | None = None):
    """MM energies (kJ/mol) at the QM scan points, relaxed with restraints.

    ``points`` carry ``phi_deg`` and ``coordinates`` (angstrom). With
    ``return_geometries`` also the relaxed coordinates (angstrom).
    ``lock_hydrogens``: see ``locked_dihedrals``.
    """
    engine = engine or DEFAULT_ENGINE
    if engine not in ENGINES:
        raise ValueError(f"MM engine must be one of {ENGINES}, not {engine!r}")
    locks = locked_dihedrals(molecule, drive[1:3], lock_hydrogens)
    relax = _openmm_profile if engine == "openmm" else _gromacs_profile
    energies, geometries = relax(Path(top), points, tuple(drive), locks)
    return (energies, geometries) if return_geometries else energies


# --------------------------------------------------------------- OpenMM ----
def _openmm_profile(top: Path, points, drive, locks):
    import openmm as mm
    from openmm import app, unit

    system = app.GromacsTopFile(str(top)).createSystem(nonbondedMethod=app.NoCutoff)
    for force in system.getForces():
        force.setForceGroup(0)
    expression = "0.5*k*min(dt, 2*pi-dt)^2; dt=abs(theta-t0); pi=3.141592653589793"
    restraint = mm.CustomTorsionForce(expression)
    restraint.addPerTorsionParameter("t0")
    restraint.addGlobalParameter("k", K_DRIVE)
    restraint.addTorsion(*drive, [0.0])
    restraint.setForceGroup(1)
    system.addForce(restraint)
    lock = mm.CustomTorsionForce(expression.replace("k*", "klock*"))
    lock.addPerTorsionParameter("t0")
    lock.addGlobalParameter("klock", K_LOCK)
    for quartet in locks:
        lock.addTorsion(*quartet, [0.0])
    lock.setForceGroup(2)
    system.addForce(lock)
    platform = mm.Platform.getPlatformByName(MM_PLATFORM)
    properties = {"Threads": "1"} if MM_PLATFORM == "CPU" else {}
    context = mm.Context(system, mm.VerletIntegrator(1e-3), platform, properties)
    energies, geometries = [], []
    for point in points:
        restraint.setTorsionParameters(0, *drive, [math.radians(point["phi_deg"])])
        restraint.updateParametersInContext(context)
        for index, quartet in enumerate(locks):
            value = math.radians(dihedral_deg(point["coordinates"], quartet))
            lock.setTorsionParameters(index, *quartet, [value])
        lock.updateParametersInContext(context)
        context.setPositions(point["coordinates"] * 0.1 * unit.nanometer)
        mm.LocalEnergyMinimizer.minimize(context, EM_TOLERANCE, 0)
        state = context.getState(getEnergy=True, getPositions=True, groups={0})
        energies.append(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
        positions = state.getPositions(asNumpy=True)
        geometries.append(np.asarray(positions.value_in_unit(unit.angstrom)))
    return energies, geometries


# -------------------------------------------------------------- GROMACS ----
def _gromacs_mdp(cutoff_nm: float) -> str:
    return "\n".join([
        "integrator = l-bfgs", "nsteps = 100000", f"emtol = {EM_TOLERANCE}",
        "cutoff-scheme = Verlet", "pbc = xyz", "nstlist = 10", "verlet-buffer-tolerance = -1",
        f"rlist = {cutoff_nm:.3f}", f"rcoulomb = {cutoff_nm:.3f}", f"rvdw = {cutoff_nm:.3f}",
        "coulombtype = cut-off", "coulomb-modifier = None",
        "vdwtype = cut-off", "vdw-modifier = None", "DispCorr = no",
        "nstenergy = 0", "nstlog = 0", "nstcalcenergy = 1", "",
    ])


def _atom_names(itp: Path) -> list[str]:
    names, section = [], None
    for line in itp.read_text().splitlines():
        stripped = line.split(";")[0].strip()
        if stripped.startswith("["):
            section = stripped.strip("[] ")
        elif section == "atoms" and stripped:
            names.append(stripped.split()[4])
    return names


def _write_g96(path: Path, coordinates_nm: np.ndarray, box_nm: float, names: list[str]) -> None:
    lines = ["TITLE", "scan point", "END", "POSITION"]
    for index, ((x, y, z), name) in enumerate(zip(coordinates_nm, names), 1):
        lines.append(f"{1:5d} {'MOL':<5s} {name:<5s}{index:7d}{x:15.9f}{y:15.9f}{z:15.9f}")
    lines += ["END", "BOX", f"{box_nm:15.9f}{box_nm:15.9f}{box_nm:15.9f}", "END"]
    path.write_text("\n".join(lines) + "\n")


def _read_g96(path: Path) -> np.ndarray:
    lines = path.read_text().splitlines()
    start = lines.index("POSITION") + 1
    rows = []
    for line in lines[start:]:
        if line.strip() == "END":
            break
        rows.append([float(line[24:39]), float(line[39:54]), float(line[54:69])])
    return np.array(rows)


def _restraint_energy(coordinates_angstrom, restraints) -> float:
    total = 0.0
    for quartet, target_deg, k in restraints:
        delta = math.radians((dihedral_deg(coordinates_angstrom, quartet) - target_deg + 180.0)
                             % 360.0 - 180.0)
        total += 0.5 * k * delta ** 2
    return total


def _gromacs_profile(top: Path, points, drive, locks):
    from .gromacs import require_gmx

    gmx = str(require_gmx())
    top_text = top.read_text()
    include = re.search(r'^#include "([^"]+\.itp)"\s*$', top_text, re.M)
    if include is None:
        raise RuntimeError(f"{top}: no molecule #include to attach the restraints to")
    first = np.asarray(points[0]["coordinates"]) / 10.0
    size = float(np.max(np.linalg.norm(first[:, None] - first[None], axis=-1)))
    cutoff = size + 1.0
    box = 2.0 * cutoff + 0.2
    energies, geometries = [], []
    with tempfile.TemporaryDirectory(prefix="mm_gmx_") as scratch:
        work = Path(scratch)
        (work / "em.mdp").write_text(_gromacs_mdp(cutoff))
        itp = (top.parent / include.group(1)).resolve()
        names = _atom_names(itp)
        for index, point in enumerate(points):
            coordinates = np.asarray(point["coordinates"])
            restraints = [(tuple(drive), float(point["phi_deg"]), K_DRIVE)]
            restraints += [(q, dihedral_deg(coordinates, q), K_LOCK) for q in locks]
            block = ["[ dihedral_restraints ]", "; ai aj ak al type phi dphi kfac"]
            block += [" ".join(str(a + 1) for a in q) + f" 1 {phi:.6f} 0 {k:.1f}"
                      for q, phi, k in restraints]
            text = top_text.replace(include.group(0),
                                    f'#include "{itp}"\n' + "\n".join(block) + "\n", 1)
            (work / "point.top").write_text(text)
            nm = coordinates / 10.0
            _write_g96(work / "start.g96", nm - nm.mean(axis=0) + box / 2.0, box, names)
            for command in (
                [gmx, "grompp", "-f", "em.mdp", "-c", "start.g96", "-p", "point.top",
                 # the one warning: grompp advises switched potentials for L-BFGS
                 "-o", "em.tpr", "-po", "em.out.mdp", "-maxwarn", "1"],
                [gmx, "mdrun", "-s", "em.tpr", "-deffnm", "em", "-c", "em.g96",
                 "-nt", "1", "-nb", "cpu"],
            ):
                completed = subprocess.run(command, cwd=work, capture_output=True, text=True)
                if completed.returncode != 0:
                    raise RuntimeError(f"{command[1]} failed for scan point {index}:\n"
                                       f"{(completed.stdout + completed.stderr)[-1500:]}")
            log = (work / "em.log").read_text()
            match = re.findall(r"Potential Energy\s+=\s+(\S+)", log)
            if not match:
                raise RuntimeError(f"no final energy in the GROMACS log of scan point {index}")
            relaxed = (_read_g96(work / "em.g96") - box / 2.0 + nm.mean(axis=0)) * 10.0
            energies.append(float(match[-1]) - _restraint_energy(relaxed, restraints))
            geometries.append(relaxed)
    return energies, geometries
