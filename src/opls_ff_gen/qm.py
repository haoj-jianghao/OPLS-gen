"""QM side of the on-the-fly fit: DFT levels, single points and relaxed xTB scans.

``DFTLevel`` selects the PySCF level for ``--qm dft``. ``relaxed_scan`` drives
one dihedral through constrained GFN2-xTB optimizations:

- rotatable bond: 24 points over 360 degrees (15 degree steps), each point
  started from the previous one;
- ring bond: the ring-internal dihedral, +-30 degrees in 5 degree steps when
  the ring is planar (|theta0| < 15), otherwise from theta0 through 0 to
  -theta0 and 20 degrees beyond both, scanned outwards in two branches.
"""

from __future__ import annotations

import json
import math
import re
import shutil
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rdkit import Chem

from .xtb import read_xyz, run_xtb, to_xyz

HARTREE_KJ_MOL = 2625.499639
FIT_CAP_KJ_MOL = 50.0  # scan points higher than this above the QM minimum are not fitted
ROTATABLE_STEP, RING_STEP = 15.0, 5.0
PLANAR_LIMIT, PLANAR_WINDOW, PUCKER_MARGIN = 15.0, 30.0, 20.0


@dataclass(frozen=True)
class DFTLevel:
    """DFT level for ``--qm dft`` (PySCF, restricted Kohn-Sham, density fitting).

    ``functional`` is a PySCF/libxc name (e.g. b3lyp, pbe0, m06-2x, wb97x-v),
    ``basis`` a PySCF basis name (def2-svp, def2-tzvp, aug-cc-pvdz, ...),
    ``dispersion`` a PySCF ``disp`` model (d3bj, d3zero, d4) or None.
    """

    functional: str = "b3lyp"
    basis: str = "def2-svp"
    dispersion: str | None = "d3bj"

    @property
    def label(self) -> str:
        name = self.functional.upper() + (f"-{self.dispersion.upper()}" if self.dispersion else "")
        return f"{name}/{self.basis}"

    @property
    def tag(self) -> str:
        """File-name-safe identifier, used for per-level cache files."""
        return re.sub(r"[^a-z0-9]+", "-", f"dft-{self.functional}-{self.dispersion or 'nodisp'}"
                      f"-{self.basis}".lower()).strip("-")

    def scf(self, elements, coordinates, charge: int, threads: int, dispersion: bool = True):
        """Unconverged PySCF RKS object (density fitting) for this level."""
        from pyscf import dft, gto, lib

        lib.num_threads(threads)
        molecule = gto.M(atom=list(zip(elements, map(tuple, coordinates))), basis=self.basis,
                         charge=charge, verbose=0)
        solver = dft.RKS(molecule).density_fit()
        solver.xc = self.functional
        if dispersion and self.dispersion:
            solver.disp = self.dispersion
        return solver

    def validate(self) -> None:
        """Fail fast (seconds) if PySCF cannot run this level: a water single
        point and DF Hessian."""
        try:
            solver = self.scf(["O", "H", "H"], [(0, 0, 0), (0.757, 0.586, 0), (-0.757, 0.586, 0)],
                              0, 1)
            solver.kernel()
            if not solver.converged:
                raise RuntimeError("test SCF did not converge")
            self.scf(["O", "H", "H"], [(0, 0, 0), (0.757, 0.586, 0), (-0.757, 0.586, 0)], 0, 1,
                     dispersion=False).run().Hessian().kernel()
        except Exception as exc:
            raise ValueError(f"DFT level {self.label} cannot be run with PySCF: "
                             f"{type(exc).__name__}: {exc}") from exc


DEFAULT_DFT = DFTLevel()


def xtb_energy_eh(output: str) -> float:
    return float(re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh", output)[-1])


def dft_energy_eh(elements, coordinates, charge: int, threads: int,
                  level: DFTLevel = DEFAULT_DFT) -> float:
    """DFT single point (default B3LYP-D3(BJ)/def2-SVP, density fitting)."""
    solver = level.scf(elements, coordinates, charge, threads)
    energy = solver.kernel()
    if not solver.converged:
        raise RuntimeError("DFT single point did not converge")
    return float(energy)


def dihedral_deg(coordinates: np.ndarray, atoms) -> float:
    """Dihedral (degrees) of 0-based atom indices."""
    p0, p1, p2, p3 = (coordinates[i] for i in atoms)
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return math.degrees(math.atan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def fourier_to_rb(v1, v2, v3, v4) -> tuple[float, ...]:
    """OPLS Fourier (phi, IUPAC) to GROMACS Ryckaert-Bellemans (psi = phi - 180)."""
    return (v2 + 0.5 * (v1 + v3), 0.5 * (-v1 + 3 * v3), -v2 + 4 * v4, -2 * v3, -4 * v4, 0.0)


def heavy_count(molecule: Chem.Mol, quartet) -> int:
    return sum(molecule.GetAtomWithIdx(a).GetAtomicNum() > 1 for a in quartet)


def ring_drive(molecule: Chem.Mol, j: int, k: int) -> tuple[int, int, int, int]:
    """Ring-internal dihedral about ring bond j-k with the most heavy atoms."""
    rings = [set(r) for r in molecule.GetRingInfo().AtomRings()]
    best = None
    for i in (a.GetIdx() for a in molecule.GetAtomWithIdx(j).GetNeighbors() if a.GetIdx() != k):
        for l in (a.GetIdx() for a in molecule.GetAtomWithIdx(k).GetNeighbors() if a.GetIdx() != j):
            if any({i, j, k, l} <= ring for ring in rings):
                quartet = (i, j, k, l)
                if best is None or heavy_count(molecule, quartet) > heavy_count(molecule, best):
                    best = quartet
    if best is None:
        raise ValueError(f"no ring-internal dihedral about {j + 1}-{k + 1}")
    return best


def relaxed_scan(task: dict) -> dict:
    """Relaxed xTB scan for one bond (runs in a worker process), cached in
    ``scan_dir/scan.json``; a cached scan about another dihedral is redone.

    ``task``: elements, start (angstrom), drive (0-based quartet), ring,
    charge, scan_dir.
    """
    scan_dir = Path(task["scan_dir"])
    done = scan_dir / "scan.json"
    drive = tuple(task["drive"])
    if done.exists():
        result = json.loads(done.read_text())
        if tuple(result.get("drive", ())) == drive:
            return result
        shutil.rmtree(scan_dir)
    elements = task["elements"]
    start = np.array(task["start"])
    charge_args = ["--chrg", str(task["charge"])] if task.get("charge") else []
    clock = time.perf_counter()
    run_xtb(scan_dir / "start", to_xyz(elements, start), ["--opt", "tight", *charge_args])
    optimized = read_xyz((scan_dir / "start" / "xtbopt.xyz").read_text())
    theta0 = dihedral_deg(optimized, drive)
    if task["ring"]:
        if abs(theta0) < PLANAR_LIMIT:
            low, high = theta0 - PLANAR_WINDOW, theta0 + PLANAR_WINDOW
        else:
            a = theta0 + math.copysign(PUCKER_MARGIN, theta0)
            b = -theta0 - math.copysign(PUCKER_MARGIN, theta0)
            low, high = min(a, b), max(a, b)
        steps = range(int(math.floor((low - theta0) / RING_STEP)),
                      int(math.ceil((high - theta0) / RING_STEP)) + 1)
        offsets = [m * RING_STEP for m in steps]
        branches = [sorted(o for o in offsets if o >= 0), sorted((o for o in offsets if o < 0), reverse=True)]
    else:
        branches = [[m * ROTATABLE_STEP for m in range(int(360 / ROTATABLE_STEP))]]
    atoms = ",".join(str(i + 1) for i in drive)
    points = []
    for branch in branches:
        seed = optimized
        for offset in branch:
            target = theta0 + offset
            value = (target + 180.0) % 360.0 - 180.0
            point_dir = scan_dir / f"p{len(points):03d}"
            constrain = f"$constrain\n  force constant=5.0\n  dihedral: {atoms},{value:.4f}\n$end\n"
            run_xtb(point_dir, to_xyz(elements, seed), ["--opt", "tight", *charge_args], constrain)
            geometry = read_xyz((point_dir / "xtbopt.xyz").read_text())
            seed = geometry
            points.append({"dir": str(point_dir), "target_offset": offset,
                           "phi_deg": dihedral_deg(geometry, drive)})
    result = {"drive": list(drive), "theta0": theta0, "points": points,
              "xtb_seconds": round(time.perf_counter() - clock, 1)}
    done.write_text(json.dumps(result, indent=2) + "\n")
    return result
