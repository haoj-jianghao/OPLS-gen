"""Check the GROMACS MM engine against OpenMM on a batch of generated topologies.

Two checks per molecule directory (``run_batch_opls2020.py`` output with a
``qm_fit/`` that has scans, or any directory given with ``--perturb``):

``single``
    Every scan geometry (or, with ``--perturb N``, N random displacements of
    the input geometry): single-point energies per component and forces.
    GROMACS: ``mdrun -rerun`` with the engine's vacuum settings (``mm.py``);
    OpenMM: ``NoCutoff``. Components: bonds, angles, proper (RB) and improper
    torsions, and LJ and Coulomb split into 1-4 pairs and longer-range pairs.
``relaxed``
    Every scan: the restrained relaxed profile (``mm.mm_profile``, the fit's
    settings) with both engines; profiles are compared relative to their own
    minimum, and the relaxed geometries by RMSD.

    python src/scripts/check_mm_engines.py BATCH_DIR --out report.json [--workers 10]
        [--perturb 20 DIR ...]
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from opls_ff_gen import mm  # noqa: E402
from opls_ff_gen.gromacs import require_gmx  # noqa: E402
from opls_ff_gen.xtb import read_xyz  # noqa: E402

# gmx energy selection names (spaces shown as "-")
GMX_TERMS = {"bond": "Bond", "angle": "Angle", "proper": "Ryckaert-Bell.",
             "improper": "Per.-Imp.-Dih.", "lj14": "LJ-14", "coulomb14": "Coulomb-14",
             "lj_sr": "LJ-(SR)", "coulomb_sr": "Coulomb-(SR)", "potential": "Potential"}


# ----------------------------------------------------------------- OpenMM ----
def openmm_components(top: Path, frames_angstrom: list[np.ndarray]):
    """Per-component energies and forces (kJ/mol, kJ/mol/nm) per frame.

    ``GromacsTopFile`` (comb-rule 3) builds Coulomb as a ``NonbondedForce``
    (1-4 pairs as its exceptions), LJ as a ``CustomNonbondedForce`` and the 1-4
    LJ pairs as a ``CustomBondForce``; any LJ left in the ``NonbondedForce`` is
    split the same way as Coulomb.
    """
    import openmm as omm
    from openmm import app, unit

    def build(charges: bool, lj: bool, pairs: bool):
        """Context whose NonbondedForce keeps only the selected parts."""
        system = app.GromacsTopFile(str(top)).createSystem(nonbondedMethod=app.NoCutoff)
        for index, force in enumerate(system.getForces()):
            force.setForceGroup(index)
            if isinstance(force, omm.NonbondedForce):
                for i in range(force.getNumParticles()):
                    q, sigma, epsilon = force.getParticleParameters(i)
                    force.setParticleParameters(i, q if charges and pairs else 0.0, sigma,
                                                epsilon if lj and pairs else 0.0)
                for i in range(force.getNumExceptions()):
                    a, b, qq, sigma, epsilon = force.getExceptionParameters(i)
                    force.setExceptionParameters(i, a, b, qq if charges else 0.0, sigma,
                                                 epsilon if lj else 0.0)
        context = omm.Context(system, omm.VerletIntegrator(1e-3),
                              omm.Platform.getPlatformByName("Reference"))
        return system, context

    system, full = build(True, True, True)
    groups: dict[str, set[int]] = {}
    for index, force in enumerate(system.getForces()):
        key = {"HarmonicBondForce": "bond", "HarmonicAngleForce": "angle",
               "RBTorsionForce": "proper", "PeriodicTorsionForce": "improper",
               "CustomNonbondedForce": "lj_sr", "CustomBondForce": "lj14",
               "NonbondedForce": "nonbonded"}.get(type(force).__name__)
        if key:
            groups.setdefault(key, set()).add(index)
    nonbonded = groups.pop("nonbonded", set())
    variants = {"coulomb_all": build(True, False, True)[1], "coulomb14": build(True, False, False)[1],
                "lj_all": build(False, True, True)[1], "lj14": build(False, True, False)[1]}

    def energy(context, selected):
        if not selected:
            return 0.0
        return (context.getState(getEnergy=True, groups=selected).getPotentialEnergy()
                .value_in_unit(unit.kilojoule_per_mole))

    rows, forces = [], []
    for frame in frames_angstrom:
        positions = frame * 0.1 * unit.nanometer
        full.setPositions(positions)
        row = {key: energy(full, selected) for key, selected in groups.items()}
        for key in ("bond", "angle", "proper", "improper", "lj_sr", "lj14"):
            row.setdefault(key, 0.0)
        part = {}
        for key, context in variants.items():
            context.setPositions(positions)
            part[key] = energy(context, nonbonded)
        row["coulomb14"] = part["coulomb14"]
        row["coulomb_sr"] = part["coulomb_all"] - part["coulomb14"]
        row["lj14"] += part["lj14"]
        row["lj_sr"] += part["lj_all"] - part["lj14"]
        state = full.getState(getEnergy=True, getForces=True)
        row["potential"] = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        rows.append(row)
        forces.append(np.asarray(state.getForces(asNumpy=True)
                                 .value_in_unit(unit.kilojoule_per_mole / unit.nanometer)))
    return rows, forces


# ---------------------------------------------------------------- GROMACS ----
def _g96_trajectory(path: Path, frames_nm, box, names) -> None:
    lines = ["TITLE", "rerun frames", "END"]
    for step, frame in enumerate(frames_nm):
        lines += ["TIMESTEP", f"{step:15d}{float(step):15.6f}", "END", "POSITION"]
        for index, ((x, y, z), name) in enumerate(zip(frame, names), 1):
            lines.append(f"{1:5d} {'MOL':<5s} {name:<5s}{index:7d}{x:15.9f}{y:15.9f}{z:15.9f}")
        lines += ["END", "BOX", f"{box:15.9f}{box:15.9f}{box:15.9f}", "END"]
    path.write_text("\n".join(lines) + "\n")


def gromacs_components(top: Path, frames_angstrom: list[np.ndarray]):
    gmx = str(require_gmx())
    frames_nm = [np.asarray(f) / 10.0 for f in frames_angstrom]
    size = max(float(np.max(np.linalg.norm(f[:, None] - f[None], axis=-1))) for f in frames_nm)
    cutoff = size + 1.0
    box = 2.0 * cutoff + 0.2
    itp = top.parent / re.search(r'^#include "([^"]+\.itp)"', top.read_text(), re.M).group(1)
    names = mm._atom_names(itp)
    shift = box / 2.0
    with tempfile.TemporaryDirectory(prefix="mmcheck_") as scratch:
        work = Path(scratch)
        (work / "sp.mdp").write_text(mm._gromacs_mdp(cutoff).replace(
            "integrator = l-bfgs", "integrator = md").replace("nstenergy = 0", "nstenergy = 1")
            + "nstfout = 1\nnstxout = 0\n")
        centred = [f - f.mean(axis=0) + shift for f in frames_nm]
        _g96_trajectory(work / "frames.g96", centred, box, names)
        mm._write_g96(work / "start.g96", centred[0], box, names)
        (work / "sp.top").write_text(top.read_text().replace(
            f'"{itp.name}"', f'"{itp.resolve()}"'))
        for command in (
            [gmx, "grompp", "-f", "sp.mdp", "-c", "start.g96", "-p", "sp.top", "-o", "sp.tpr",
             "-maxwarn", "1"],
            [gmx, "mdrun", "-s", "sp.tpr", "-rerun", "frames.g96", "-deffnm", "sp", "-nt", "1",
             "-nb", "cpu"],
        ):
            done = subprocess.run(command, cwd=work, capture_output=True, text=True)
            if done.returncode:
                raise RuntimeError((done.stdout + done.stderr)[-1500:])
        listed = subprocess.run([gmx, "energy", "-f", "sp.edr"], cwd=work, input="\n",
                                capture_output=True, text=True)
        listing = listed.stdout + listed.stderr
        wanted = set(GMX_TERMS.values())
        available = {label: int(number) for number, label in re.findall(r"(\d+)\s+(\S+)", listing)
                     if label in wanted}
        # xvg columns follow the energy-file order, not the selection order
        selection = sorted(((key, available[term]) for key, term in GMX_TERMS.items()
                            if term in available), key=lambda item: item[1])
        done = subprocess.run([gmx, "energy", "-f", "sp.edr", "-o", "energy.xvg", "-xvg", "none"],
                              cwd=work, input="\n".join(str(n) for _, n in selection) + "\n\n",
                              capture_output=True, text=True)
        if done.returncode:
            raise RuntimeError(f"gmx energy {selection}: {(done.stdout + done.stderr)[-800:]}")
        table = np.loadtxt(work / "energy.xvg", ndmin=2)
        rows = [{key: float(line[1 + i]) for i, (key, _) in enumerate(selection)} for line in table]
        for row in rows:
            for key in GMX_TERMS:
                row.setdefault(key, 0.0)
        subprocess.run([gmx, "traj", "-f", "sp.trr", "-s", "sp.tpr", "-of", "forces.xvg",
                        "-xvg", "none"], cwd=work, input="0\n", capture_output=True, text=True,
                       check=True)
        raw = np.loadtxt(work / "forces.xvg", ndmin=2)[:, 1:]
        forces = [line.reshape(-1, 3) for line in raw]
    return rows, forces


# ------------------------------------------------------------------ checks ----
def _scan_frames(directory: Path):
    scans = []
    for scan_json in sorted((directory / "qm_fit").glob("scan_*/scan.json")):
        scan = json.loads(scan_json.read_text())
        points = [{"phi_deg": p["phi_deg"],
                   "coordinates": read_xyz((Path(p["dir"]) / "xtbopt.xyz").read_text())}
                  for p in scan["points"]]
        scans.append({"name": scan_json.parent.name, "drive": tuple(scan["drive"]),
                      "points": points})
    return scans


def check_directory(job) -> dict:
    directory, perturb, relaxed = job
    from rdkit import Chem

    directory = Path(directory)
    topology = json.loads(next(directory.glob("*.topology.json")).read_text())
    top = Path(topology["files"]["top"])
    molecule = Chem.MolFromMolFile(str(next(directory.glob("*.assigned.sdf"))), removeHs=False)
    scans = _scan_frames(directory)
    if perturb:
        rng = np.random.default_rng(7)
        start = molecule.GetConformer().GetPositions()
        frames = [start + rng.normal(0.0, 0.05, start.shape) for _ in range(perturb)]
    else:
        frames = [p["coordinates"] for scan in scans for p in scan["points"]]
    result = {"name": directory.name, "charge": Chem.GetFormalCharge(molecule),
              "atoms": molecule.GetNumAtoms(), "frames": len(frames)}
    if not frames:
        return result
    omm_rows, omm_forces = openmm_components(top, frames)
    gmx_rows, gmx_forces = gromacs_components(top, frames)
    diffs = {key: float(max(abs(a[key] - b[key]) for a, b in zip(omm_rows, gmx_rows)))
             for key in GMX_TERMS}
    result["max_abs_energy_diff_kj_mol"] = diffs
    result["max_abs_force_diff_kj_mol_nm"] = float(max(
        np.abs(a - b).max() for a, b in zip(omm_forces, gmx_forces)))
    result["max_abs_force_kj_mol_nm"] = float(max(np.abs(a).max() for a in omm_forces))
    if relaxed and not perturb:
        result["scans"] = []
        for scan in scans:
            profiles, geometries = {}, {}
            for engine in mm.ENGINES:
                energies, relaxed_geometries = mm.mm_profile(
                    top, scan["points"], scan["drive"], molecule, return_geometries=True,
                    lock_hydrogens=True, engine=engine)
                profiles[engine] = np.asarray(energies) - min(energies)
                geometries[engine] = relaxed_geometries
            rmsd = max(float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))
                       for a, b in zip(geometries["openmm"], geometries["gromacs"]))
            result["scans"].append({
                "scan": scan["name"], "points": len(scan["points"]),
                "max_abs_profile_diff_kj_mol": float(np.abs(profiles["openmm"]
                                                            - profiles["gromacs"]).max()),
                "max_rmsd_angstrom": rmsd})
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("batch", nargs="?", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--no-relaxed", action="store_true")
    parser.add_argument("--perturb", type=int, default=0,
                        help="single points on N random displacements (0.05 A) of the input "
                        "geometry of each directory in --dirs")
    parser.add_argument("--dirs", nargs="*", type=Path, default=[])
    args = parser.parse_args()
    jobs = [(d, args.perturb, False) for d in args.dirs]
    if args.batch:
        jobs += [(d, 0, not args.no_relaxed) for d in sorted(args.batch.iterdir())
                 if (d / "qm_fit").is_dir() and list((d / "qm_fit").glob("scan_*/scan.json"))]
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for result in pool.map(check_directory, jobs):
            results.append(result)
            worst = max(result.get("max_abs_energy_diff_kj_mol", {"-": 0}).values())
            scans = result.get("scans", [])
            profile = max((s["max_abs_profile_diff_kj_mol"] for s in scans), default=0.0)
            print(f"{result['name']}: {result['frames']} frames, max |dE| {worst:.2e}, "
                  f"max |dF| {result.get('max_abs_force_diff_kj_mol_nm', 0):.2e}, "
                  f"{len(scans)} relaxed scans max |dprofile| {profile:.3f} kJ/mol", flush=True)
    args.out.write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
