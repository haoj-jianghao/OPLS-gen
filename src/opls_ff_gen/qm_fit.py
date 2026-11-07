"""On-the-fly QM fitting of bonded terms missing from the OPLS/2020 tables.

``generate_gromacs_topology(..., qm_method=...)`` writes every bond, angle and
proper torsion that the parameter table lacks as a ``qm-pending`` placeholder
and calls :func:`fit_pending_terms`, which replaces them with values fitted to
QM for this molecule (rule ``qm-fit``). Two levels are available:

``xtb``
    GFN2-xTB for everything: optimized geometry, Hessian, relaxed torsion
    scans and scan energies.
``dft``
    A DFT level (``qm.DFTLevel``; default B3LYP-D3(BJ)/def2-SVP;
    PySCF, density fitting): geometry optimized with geomeTRIC from the xTB
    geometry, analytic Hessian (same functional and basis, no dispersion
    term), and single points on the xTB relaxed scan geometries
    (DFT//GFN2-xTB). Cache files carry the level's tag, so runs at different
    levels in one directory never reuse each other's results.

Bonds and angles
    Equilibrium values are the QM-optimized bond lengths and angles; force
    constants come from the QM Hessian by the modified Seminario method
    (A. E. A. Allen, M. C. Payne and D. J. Cole, J. Chem. Theory Comput. 14,
    274 (2018)). Terms related by molecular graph symmetry are averaged.

Torsions
    Every central bond with missing quartets is scanned (one bond per set of
    symmetry-equivalent bonds): rotatable bonds over 360 degrees in 15 degree
    steps, ring bonds over the ring-internal dihedral (``qm.relaxed_scan``),
    each point a constrained GFN2-xTB optimization. The MM baseline
    is the molecule's topology with the fitted bonds/angles, every table
    torsion kept and the missing torsions zeroed.

    Chemical context is respected when several torsions share one bond:
    - torsions already in the table are never changed; the fit only supplies
      what they miss (QM - MM residual);
    - each missing quartet is keyed by its class quartet plus the SMARTS
      context of its four atoms (``context.context_key``), and every quartet
      with the same key, about any bond of the molecule, shares one set of
      OPLS Fourier V1-V4 (V1-V3 for ring bonds); all scans are fitted jointly;
    - per bond, the missing quartets with the most heavy atoms carry the term
      (a key that carries the term on any bond carries it everywhere); the
      other missing quartets are zero. This removes the degeneracy between
      quartets that differ only by a phase about the same bond.
    The fit is a damped self-consistent ridge fit: relax MM (``mm.mm_profile``,
    OpenMM or GROMACS) at every scan point with the current parameters (driven dihedral
    restrained; every other rotatable bond, including OH/NH2/CH3 rotors, held
    at its QM value), linearize, solve min |A x - r|^2 + ridge |x|^2, and take
    the largest of the steps 1, 1/2, 1/4, 1/8 that lowers the relaxed
    objective; stop when no step does or the parameters stop changing.

Everything is written to ``<output>/qm_fit/``: geometries, Hessian, scans,
``profile_*.csv`` per scan and ``report.json`` with fitted values, RMSEs and
timings.

Caches are reused only for the same inputs. ``qm_fit/manifest.json`` holds a
fingerprint of the atoms, bond graph and orders, formal and net charges,
starting coordinates, xtb version and scan protocol; a directory with a
different (or no) fingerprint is moved aside to ``qm_fit.stale-<time>`` and
everything is recomputed. Level-specific files (DFT geometry, Hessian, scan
energies) carry the level's tag; the manifest also records the PySCF version
per DFT level and drops that level's files when it changes. A scan whose
driven dihedral differs from the one requested is redone (``qm.relaxed_scan``).
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from rdkit import Chem

from .bonded import RULE_QM_PENDING, BondedAssignment
from .context import context_key
from .mm import DEFAULT_ENGINE, ENGINES, mm_profile
from .qm import (
    DEFAULT_DFT,
    FIT_CAP_KJ_MOL,
    HARTREE_KJ_MOL,
    PLANAR_LIMIT,
    PLANAR_WINDOW,
    PUCKER_MARGIN,
    RING_STEP,
    ROTATABLE_STEP,
    DFTLevel,
    dft_energy_eh,
    fourier_to_rb,
    heavy_count,
    relaxed_scan,
    ring_drive,
    xtb_energy_eh,
)
from .xtb import read_xyz, run_xtb, to_xyz

RULE_QM_FIT = "qm-fit"
# Missing torsion about a linear (sp) atom: its dihedral is ill-defined, so
# it is neither scanned nor fitted, and is written as zero with a warning.
RULE_ZERO_LINEAR = "zero-linear"
METHODS = ("xtb", "dft")
RIDGE = 0.1
ITERATIONS = 8
CONVERGED_V = 0.05
BOHR_NM = 0.0529177210903
HESSIAN_TO_KJ_MOL_NM2 = HARTREE_KJ_MOL / BOHR_NM ** 2


def fourier_basis(phi_deg: float, order: int) -> np.ndarray:
    phi = math.radians(phi_deg)
    return np.array([0.5 * (1 + (-1) ** (n + 1) * math.cos(n * phi)) for n in range(1, order + 1)])


def solve(a: np.ndarray, r: np.ndarray, x0: np.ndarray, ridge: float):
    """Ridge correction to the prior: min |A x - r|^2 + ridge |x - x0|^2.

    Returns x and the diagonal of the resolution matrix (1 = set by the scans,
    0 = left at the prior).
    """
    u, s, vt = np.linalg.svd(a, full_matrices=False)
    gain = s / (s ** 2 + ridge)
    x = x0 + vt.T @ (gain * (u.T @ (r - a @ x0)))
    determined = np.clip(np.einsum("ki,k,ki->i", vt, s ** 2 / (s ** 2 + ridge), vt), 0.0, 1.0)
    return x, determined


def method_label(method: str, level: DFTLevel = DEFAULT_DFT) -> str:
    if method == "dft":
        dft = level.label
        return (f"{dft} (PySCF, density fitting): geometry and Hessian; torsion scans "
                f"{dft}//GFN2-xTB")
    return "GFN2-xTB (xtb 6.7): geometry, Hessian, torsion scans"


# ---------------------------------------------------------------- cache ----
CACHE_VERSION = 1  # bump when cached files change meaning


def cache_inputs(molecule: Chem.Mol, charge: int) -> dict:
    """Everything the cached QM results depend on, apart from the level."""
    from .xtb import xtb_version

    return {
        "cache_version": CACHE_VERSION,
        "elements": [atom.GetSymbol() for atom in molecule.GetAtoms()],
        "formal_charges": [atom.GetFormalCharge() for atom in molecule.GetAtoms()],
        "net_charge": charge,
        "bonds": sorted([bond.GetBeginAtomIdx() + 1, bond.GetEndAtomIdx() + 1,
                         bond.GetBondTypeAsDouble()] for bond in molecule.GetBonds()),
        "start_angstrom": np.round(molecule.GetConformer().GetPositions(), 4).tolist(),
        "xtb": xtb_version(),
        "scan_protocol": {"rotatable_step": ROTATABLE_STEP, "ring_step": RING_STEP,
                          "planar_limit": PLANAR_LIMIT, "planar_window": PLANAR_WINDOW,
                          "pucker_margin": PUCKER_MARGIN},
    }


def _fingerprint(inputs: dict) -> str:
    return hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest()


def prepare_cache(workdir: Path, inputs: dict, method: str, level: DFTLevel) -> list[str]:
    """Make ``workdir`` a cache for ``inputs`` only; returns warnings about
    caches that were set aside or dropped."""
    warnings = []
    manifest_path = workdir / "manifest.json"
    fingerprint = _fingerprint(inputs)
    manifest = {}
    if manifest_path.exists():
        try:
            manifest = json.loads(manifest_path.read_text())
        except json.JSONDecodeError:
            manifest = {}
    if workdir.exists() and any(workdir.iterdir()) and manifest.get("fingerprint") != fingerprint:
        base = f"{workdir.name}.stale-{time.strftime('%Y%m%d-%H%M%S')}"
        stale, suffix = workdir.with_name(base), 1
        while stale.exists():
            stale, suffix = workdir.with_name(f"{base}-{suffix}"), suffix + 1
        workdir.rename(stale)
        warnings.append(f"QM cache {workdir} was made for a different molecule, charge, "
                        f"start geometry or xtb version; moved to {stale.name}, recomputed")
        manifest = {}
    workdir.mkdir(parents=True, exist_ok=True)
    levels = manifest.get("levels", {})
    if method == "dft":
        import pyscf

        tag = _tag(method, level)
        if tag in levels and levels[tag].get("pyscf") != pyscf.__version__:
            for path in workdir.rglob(f"*_{tag}.*"):
                path.unlink()
            warnings.append(f"QM cache: {tag} results from PySCF {levels[tag].get('pyscf')} "
                            f"dropped (now {pyscf.__version__}), recomputed")
        levels[tag] = {"pyscf": pyscf.__version__}
    manifest_path.write_text(json.dumps({"fingerprint": fingerprint, "levels": levels,
                                         "inputs": inputs}) + "\n")
    return warnings


# ------------------------------------------------------------- geometry ----
def _tag(method: str, level: DFTLevel) -> str:
    """Cache-file tag: "xtb", or the DFT level (so levels never share caches)."""
    return level.tag if method == "dft" else "xtb"


def optimized_geometry(elements, start, method: str, charge: int, threads: int,
                       workdir: Path, level: DFTLevel = DEFAULT_DFT) -> np.ndarray:
    """QM-optimized coordinates (angstrom); cached in ``workdir``."""
    cache = workdir / f"opt_{_tag(method, level)}.xyz"
    if cache.exists():
        return read_xyz(cache.read_text())
    run_xtb(workdir / "opt_xtb", to_xyz(elements, start), ["--opt", "tight", "--chrg", str(charge)])
    geometry = read_xyz((workdir / "opt_xtb" / "xtbopt.xyz").read_text())
    if method == "dft":
        from pyscf.geomopt.geometric_solver import optimize

        mf = level.scf(elements, geometry, charge, threads)
        optimized = optimize(mf, maxsteps=200)
        geometry = optimized.atom_coords(unit="Angstrom")
    cache.write_text(to_xyz(elements, geometry))
    return geometry


def hessian(elements, geometry, method: str, charge: int, threads: int,
            workdir: Path, level: DFTLevel = DEFAULT_DFT) -> np.ndarray:
    """Cartesian Hessian (3N x 3N) in kJ/mol/nm^2; cached in ``workdir``."""
    cache = workdir / f"hessian_{_tag(method, level)}.npy"
    if cache.exists():
        return np.load(cache)
    n = len(elements)
    if method == "xtb":
        directory = workdir / "hess_xtb"
        run_xtb(directory, to_xyz(elements, geometry), ["--hess", "--chrg", str(charge)])
        text = (directory / "hessian").read_text().replace("$hessian", "")
        values = np.array([float(v) for v in text.split()])
        matrix = values.reshape(3 * n, 3 * n)
    else:
        mf = level.scf(elements, geometry, charge, threads, dispersion=False)
        mf.kernel()
        blocks = mf.Hessian().kernel()  # (n, n, 3, 3), Eh/bohr^2
        matrix = blocks.transpose(0, 2, 1, 3).reshape(3 * n, 3 * n)
    matrix = 0.5 * (matrix + matrix.T) * HESSIAN_TO_KJ_MOL_NM2
    np.save(cache, matrix)
    return matrix


# ------------------------------------------------------ modified Seminario ----
def _block(matrix: np.ndarray, a: int, b: int) -> tuple[np.ndarray, np.ndarray]:
    """Eigenvalues/vectors of -H_ab (the interatomic force-constant block)."""
    block = -matrix[3 * a:3 * a + 3, 3 * b:3 * b + 3]
    values, vectors = np.linalg.eig(block)
    return np.real(values), np.real(vectors)


def _projected(matrix, a, b, direction) -> float:
    values, vectors = _block(matrix, a, b)
    return float(sum(values[i] * abs(np.dot(direction, vectors[:, i])) for i in range(3)))


def seminario_bond(matrix, coordinates_nm, a: int, b: int) -> tuple[float, float]:
    """(r0 nm, k kJ/mol/nm^2) for bond a-b."""
    vector = coordinates_nm[b] - coordinates_nm[a]
    length = float(np.linalg.norm(vector))
    unit_vector = vector / length
    k = 0.5 * (_projected(matrix, a, b, unit_vector) + _projected(matrix, b, a, -unit_vector))
    return length, k


def _perpendicular(u_ab, u_cb):
    normal = np.cross(u_cb, u_ab)
    norm = np.linalg.norm(normal)
    if norm < 1e-6:  # linear angle: any normal to the bond
        trial = np.array([1.0, 0.0, 0.0]) if abs(u_ab[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        normal = np.cross(u_ab, trial)
        norm = np.linalg.norm(normal)
    normal /= norm
    return normal, np.cross(normal, u_ab), np.cross(u_cb, normal)


def seminario_angle(matrix, coordinates_nm, molecule: Chem.Mol,
                    a: int, b: int, c: int) -> tuple[float, float]:
    """(theta0 deg, k kJ/mol/rad^2) for angle a-b-c (modified Seminario)."""
    r_ab = coordinates_nm[a] - coordinates_nm[b]
    r_cb = coordinates_nm[c] - coordinates_nm[b]
    len_ab, len_cb = np.linalg.norm(r_ab), np.linalg.norm(r_cb)
    u_ab, u_cb = r_ab / len_ab, r_cb / len_cb
    theta = math.degrees(math.acos(np.clip(np.dot(u_ab, u_cb), -1.0, 1.0)))
    _, u_pa, u_pc = _perpendicular(u_ab, u_cb)

    def scale(outer: int, u_p: np.ndarray, u_bond: np.ndarray) -> float:
        """Modified-Seminario factor for the other angles sharing bond outer-b."""
        others = []
        for neighbor in molecule.GetAtomWithIdx(b).GetNeighbors():
            d = neighbor.GetIdx()
            if d in (a, c):
                continue
            r_db = coordinates_nm[d] - coordinates_nm[b]
            _, u_p_other, _ = _perpendicular(u_bond, r_db / np.linalg.norm(r_db))
            others.append(abs(np.dot(u_p, u_p_other)) ** 2)
        return 1.0 + (sum(others) / len(others) if others else 0.0)

    k_a = len_ab ** 2 * _projected(matrix, a, b, u_pa) / scale(a, u_pa, u_ab)
    k_c = len_cb ** 2 * _projected(matrix, c, b, u_pc) / scale(c, u_pc, u_cb)
    k = 1.0 / (1.0 / k_a + 1.0 / k_c)
    return theta, k


# ------------------------------------------------------------- torsions ----
def _drive(molecule: Chem.Mol, j: int, k: int) -> tuple[int, int, int, int]:
    """Heavy-atom dihedral to drive about j-k (ring-internal for ring bonds)."""
    if molecule.GetBondBetweenAtoms(j, k).IsInRing():
        return ring_drive(molecule, j, k)
    quartets = [(i.GetIdx(), j, k, l.GetIdx())
                for i in molecule.GetAtomWithIdx(j).GetNeighbors() if i.GetIdx() != k
                for l in molecule.GetAtomWithIdx(k).GetNeighbors()
                if l.GetIdx() != j and l.GetIdx() != i.GetIdx()]
    return max(quartets, key=lambda q: (heavy_count(molecule, q), [-a for a in q]))


def _linear(atom: Chem.Atom) -> bool:
    return atom.GetHybridization() == Chem.HybridizationType.SP


def _rewrite_itp(itp: Path, top: Path, bonds: dict, angles: dict, torsions: dict,
                 suffix: str) -> tuple[Path, Path]:
    """Copy of the ITP/TOP with the given 1-based bonds (b0 nm, k), angles
    (theta deg, k) and RB torsions (C0-C5) replaced."""
    lines = itp.read_text().splitlines()
    section = None
    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("["):
            section = stripped
            continue
        fields = stripped.split(";")[0].split()
        if section == "[ bonds ]" and len(fields) >= 5:
            atoms = tuple(int(x) for x in fields[:2])
            if atoms in bonds:
                b0, k = bonds[atoms]
                lines[index] = f"{atoms[0]:5d} {atoms[1]:5d} {1:5d} {b0:13.7f} {k:18.6f} ; {suffix}"
        elif section == "[ angles ]" and len(fields) >= 6:
            atoms = tuple(int(x) for x in fields[:3])
            if atoms in angles:
                theta, k = angles[atoms]
                lines[index] = (f"{atoms[0]:5d} {atoms[1]:5d} {atoms[2]:5d} {1:5d} "
                                f"{theta:14.7f} {k:20.6f} ; {suffix}")
        elif section == "[ dihedrals ]" and len(fields) >= 11 and fields[4] == "3":
            atoms = tuple(int(x) for x in fields[:4])
            if atoms in torsions:
                lines[index] = ("".join(f"{a:5d} " for a in atoms) + f"{3:5d} "
                                + " ".join(f"{c:11.6f}" for c in torsions[atoms])
                                + f" ; {suffix}")
    new_itp = itp.with_name(f"{itp.stem}.{suffix}.itp")
    new_itp.write_text("\n".join(lines) + "\n")
    new_top = top.with_name(f"{top.stem}.{suffix}.top")
    new_top.write_text(top.read_text().replace(f'"{itp.name}"', f'"{new_itp.name}"'))
    return new_itp, new_top


def _relative(values) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return values - values.min()


DFT_THREADS_PER_POINT = 2  # PySCF scales poorly beyond ~2 threads for 15-40 atoms


def _point_energy(job: tuple) -> None:
    """QM single point of one scan geometry, cached in its directory."""
    directory, elements, method, charge, threads, level = job
    directory = Path(directory)
    cache = directory / f"energy_{_tag(method, level)}.json"
    if cache.exists():
        return
    geometry = read_xyz((directory / "xtbopt.xyz").read_text())
    clock = time.perf_counter()
    if method == "dft":
        energy = dft_energy_eh(elements, geometry, charge, threads, level)
    else:
        energy = xtb_energy_eh(run_xtb(directory / "sp", to_xyz(elements, geometry),
                                         ["--chrg", str(charge)]))
    cache.write_text(json.dumps({"energy_eh": energy, "level": _tag(method, level), "threads": threads,
                                 "seconds": round(time.perf_counter() - clock, 2)}) + "\n")


def _scan_energies(tasks: list[dict], elements, method: str, charge: int,
                   threads: int, level: DFTLevel = DEFAULT_DFT) -> list[list[dict]]:
    """QM energies of every scan point (several points at a time), then each
    scan's points sorted by angle with energies in kJ/mol."""
    scans = [json.loads((Path(task["scan_dir"]) / "scan.json").read_text()) for task in tasks]
    per_point = DFT_THREADS_PER_POINT if method == "dft" else 1
    tag = _tag(method, level)
    jobs = [(point["dir"], elements, method, charge, per_point, level)
            for scan in scans for point in scan["points"]
            if not (Path(point["dir"]) / f"energy_{tag}.json").exists()]
    workers = max(1, min(len(jobs), threads // per_point))
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            list(pool.map(_point_energy, jobs))
    else:
        for job in jobs:
            _point_energy(job)
    result = []
    for scan in scans:
        points = []
        for point in scan["points"]:
            directory = Path(point["dir"])
            energy = json.loads((directory / f"energy_{tag}.json").read_text())["energy_eh"]
            points.append({"phi_deg": point["phi_deg"],
                           "coordinates": read_xyz((directory / "xtbopt.xyz").read_text()),
                           "qm_kj_mol": energy * HARTREE_KJ_MOL})
        result.append(sorted(points, key=lambda p: p["phi_deg"]))
    return result


class _TorsionModel:
    """Unknowns (Fourier blocks per carrier key) and their quartets."""

    def __init__(self, molecule: Chem.Mol, pending: dict[tuple[int, ...], tuple[str, ...]]):
        self.keys = {q: context_key(molecule, q, classes) for q, classes in pending.items()}
        carriers = set()
        for bond in {frozenset(q[1:3]) for q in pending}:
            about = [q for q in pending if frozenset(q[1:3]) == bond]
            most = max(heavy_count(molecule, q) for q in about)
            carriers |= {self.keys[q] for q in about if heavy_count(molecule, q) == most}
        self.carrier_quartets = [q for q in pending if self.keys[q] in carriers]
        orders = {}
        for q in self.carrier_quartets:
            ring = molecule.GetBondBetweenAtoms(q[1], q[2]).IsInRing()
            orders[self.keys[q]] = min(orders.get(self.keys[q], 4), 3 if ring else 4)
        self.columns, offset = {}, 0
        for key in sorted(orders, key=lambda k: (k[0], k[1])):
            self.columns[key] = (offset, orders[key])
            offset += orders[key]
        self.parameters = offset
        self.pending = pending

    def rb(self, x: np.ndarray) -> dict[tuple[int, ...], tuple[float, ...]]:
        """1-based quartet -> RB coefficients for every pending quartet."""
        terms = {}
        for q in self.pending:
            key = self.keys[q]
            if key in self.columns:
                start, order = self.columns[key]
                v = np.pad(x[start:start + order], (0, 4 - order))
                terms[tuple(a + 1 for a in q)] = fourier_to_rb(*v)
            else:
                terms[tuple(a + 1 for a in q)] = (0.0,) * 6
        return terms

    def design_row(self, geometry: np.ndarray) -> np.ndarray:
        from .qm import dihedral_deg

        row = np.zeros(self.parameters)
        for q in self.carrier_quartets:
            start, order = self.columns[self.keys[q]]
            row[start:start + order] += fourier_basis(dihedral_deg(geometry, q), order)
        return row


def fit_pending_terms(
    molecule: Chem.Mol,
    bonded: BondedAssignment,
    topology_json: Path,
    *,
    method: str = "xtb",
    threads: int = 8,
    dft_level: DFTLevel = DEFAULT_DFT,
    ridge: float = RIDGE,
    iterations: int = ITERATIONS,
    mm_engine: str | None = None,
    workdir: Path | None = None,
) -> BondedAssignment:
    """Replace the ``qm-pending`` terms of ``bonded`` by QM-fitted values.

    ``mm_engine`` ("openmm" or "gromacs", default ``mm.DEFAULT_ENGINE``) runs
    the restrained MM relaxations of the torsion fit. ``workdir`` (default
    ``qm_fit`` next to ``topology_json``) holds the QM cache and the report."""
    if method not in METHODS:
        raise ValueError(f"QM method must be one of {METHODS}, not {method!r}")
    mm_engine = mm_engine or DEFAULT_ENGINE
    if mm_engine not in ENGINES:
        raise ValueError(f"MM engine must be one of {ENGINES}, not {mm_engine!r}")
    clock_start = time.perf_counter()
    timings: dict[str, float] = {}
    topology = json.loads(Path(topology_json).read_text())
    itp, top = Path(topology["files"]["itp"]), Path(topology["files"]["top"])
    workdir = Path(workdir) if workdir is not None else Path(topology_json).parent / "qm_fit"
    elements = [atom.GetSymbol() for atom in molecule.GetAtoms()]
    charge = Chem.GetFormalCharge(molecule)
    cache_warnings = prepare_cache(workdir, cache_inputs(molecule, charge), method, dft_level)
    ranks = list(Chem.CanonicalRankAtoms(molecule, breakTies=False))
    report: dict = {"method": method_label(method, dft_level), "level": method,
                    "dft_level": ({"functional": dft_level.functional, "basis": dft_level.basis,
                                   "dispersion": dft_level.dispersion}
                                  if method == "dft" else None), "bonds": [], "angles": [],
                    "scans": [], "torsion_keys": []}

    pending_bonds = [b for b in bonded.bonds if b.rule == RULE_QM_PENDING]
    pending_angles = [a for a in bonded.angles if a.rule == RULE_QM_PENDING]
    pending_torsions, linear_torsions = [], []
    for value in bonded.proper_torsions:
        if value.rule == RULE_QM_PENDING:
            j, k = value.atoms[1] - 1, value.atoms[2] - 1
            linear = _linear(molecule.GetAtomWithIdx(j)) or _linear(molecule.GetAtomWithIdx(k))
            (linear_torsions if linear else pending_torsions).append(value)

    clock = time.perf_counter()
    # The DFT geometry only feeds the bond/angle fit (torsion scans are xTB
    # constrained optimizations), so it is skipped when only torsions are missing.
    geometry_method = method if (pending_bonds or pending_angles) else "xtb"
    if pending_bonds or pending_angles or pending_torsions:  # not for zero-linear alone
        geometry = optimized_geometry(elements, molecule.GetConformer().GetPositions(),
                                      geometry_method, charge, threads, workdir, dft_level)
    timings["optimization_s"] = time.perf_counter() - clock

    # ---- bonds and angles: modified Seminario on the QM Hessian ----------
    fitted_bonds, fitted_angles = {}, {}
    if pending_bonds or pending_angles:
        clock = time.perf_counter()
        matrix = hessian(elements, geometry, method, charge, threads, workdir, dft_level)
        timings["hessian_s"] = time.perf_counter() - clock
        coordinates_nm = geometry / 10.0
        groups: dict[tuple, list] = {}
        for value in pending_bonds:
            a, b = (i - 1 for i in value.atoms)
            r0, k = seminario_bond(matrix, coordinates_nm, a, b)
            key = ("bond",) + tuple(sorted((ranks[a], ranks[b])))
            groups.setdefault(key, []).append((value.atoms, r0, k))
        for value in pending_angles:
            a, b, c = (i - 1 for i in value.atoms)
            theta, k = seminario_angle(matrix, coordinates_nm, molecule, a, b, c)
            key = ("angle", ranks[b]) + tuple(sorted((ranks[a], ranks[c])))
            groups.setdefault(key, []).append((value.atoms, theta, k))
        for key, members in groups.items():
            mean0 = float(np.mean([m[1] for m in members]))
            mean_k = float(np.mean([m[2] for m in members]))
            for atoms, _, _ in members:
                target = fitted_bonds if key[0] == "bond" else fitted_angles
                target[atoms] = (round(mean0, 5 if key[0] == "bond" else 3), round(mean_k, 1))
        report["bonds"] = [{"atoms": list(a), "b0_nm": v[0], "k_kj_mol_nm2": v[1]}
                           for a, v in fitted_bonds.items()]
        report["angles"] = [{"atoms": list(a), "theta0_deg": v[0], "k_kj_mol_rad2": v[1]}
                            for a, v in fitted_angles.items()]

    # ---- torsions: relaxed scans and joint self-consistent fit ----------
    fitted_torsions: dict[tuple[int, ...], tuple[float, ...]] = {}
    if pending_torsions:
        pending = {tuple(a - 1 for a in t.atoms): tuple(t.classes) for t in pending_torsions}
        model = _TorsionModel(molecule, pending)
        bonds_to_scan, seen = [], set()
        for quartet in model.carrier_quartets:
            j, k = sorted(quartet[1:3])
            key = tuple(sorted((ranks[j], ranks[k])))
            if key not in seen:
                seen.add(key)
                bonds_to_scan.append((j, k))
        tasks = []
        for j, k in bonds_to_scan:
            tasks.append({"bond": [j + 1, k + 1],
                          "ring": molecule.GetBondBetweenAtoms(j, k).IsInRing(),
                          "drive": list(_drive(molecule, j, k)), "elements": elements,
                          "start": geometry.tolist(), "charge": charge,
                          "scan_dir": str(workdir / f"scan_{j + 1}_{k + 1}")})
        clock = time.perf_counter()
        workers = max(1, min(threads, len(tasks)))
        if workers == 1:
            for task in tasks:
                relaxed_scan(task)
        else:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                list(pool.map(relaxed_scan, tasks))
        timings["scan_geometries_s"] = time.perf_counter() - clock
        clock = time.perf_counter()
        scans = []
        for task, points in zip(tasks, _scan_energies(tasks, elements, method, charge, threads, dft_level)):
            qm = _relative([p["qm_kj_mol"] for p in points])
            scans.append({"task": task, "points": points, "qm": qm,
                          "mask": qm <= FIT_CAP_KJ_MOL})
        timings["scan_energies_s"] = time.perf_counter() - clock

        clock = time.perf_counter()
        base_bonds = {**fitted_bonds}
        base_angles = {**fitted_angles}
        zero = {tuple(a + 1 for a in q): (0.0,) * 6 for q in pending}
        _, zero_top = _rewrite_itp(itp, top, base_bonds, base_angles, zero, "qmfit_zero")

        workers = max(1, min(threads, sum(len(scan["points"]) for scan in scans)))
        pool = ProcessPoolExecutor(max_workers=workers)
        # Points relax independently: each scan is split into chunks so that
        # all workers are busy even when there are fewer scans than workers.
        chunks = max(1, -(-workers // max(1, len(scans))))

        def relax(trial_top, with_geometries):
            jobs = []
            for scan in scans:
                points = scan["points"]
                size = -(-len(points) // chunks)
                jobs.append([pool.submit(mm_profile, trial_top, points[i:i + size],
                                         tuple(scan["task"]["drive"]), molecule, True,
                                         lock_hydrogens=True, engine=mm_engine)
                             for i in range(0, len(points), size)])
            results = []
            for parts in jobs:
                energies, geometries = [], []
                for part in parts:
                    e, g = part.result()
                    energies += e
                    geometries += g
                results.append((energies, geometries) if with_geometries else energies)
            return results

        def linearize(x):
            _, trial_top = _rewrite_itp(itp, top, base_bonds, base_angles, model.rb(x),
                                        "qmfit_trial")
            rows, rhs, profiles = [], [], []
            for scan, (energies, geometries) in zip(scans, relax(trial_top, True)):
                energies = np.asarray(energies)
                design = np.array([model.design_row(g) for g in geometries])
                baseline = energies - design @ x
                mask = scan["mask"]
                a, r = design[mask], (scan["qm"] - baseline)[mask]
                rows.append(a - a.mean(axis=0))
                rhs.append(r - r.mean())
                profiles.append(_relative(energies))
            return np.vstack(rows), np.concatenate(rhs), profiles

        def objective(profiles, x):
            total = 0.0
            for scan, profile in zip(scans, profiles):
                residual = (profile - scan["qm"])[scan["mask"]]
                total += float(np.sum((residual - residual.mean()) ** 2))
            return total + ridge * float(np.sum(x ** 2))

        for scan, energies in zip(scans, relax(zero_top, False)):
            scan["zero"] = _relative(energies)
        x = np.zeros(model.parameters)
        history = []
        if model.parameters and scans:
            # Damped self-consistent iteration: the full linearized step is
            # tried first and halved (backtracking) until the relaxed objective
            # decreases, so MM conformer jumps cannot make the fit oscillate.
            a, r, profiles = linearize(x)
            value = objective(profiles, x)
            history.append({"iteration": 0, "objective": value, "step": 0.0,
                            "max_abs_v": 0.0})
            for iteration in range(1, iterations + 1):
                x_new, _ = solve(a, r, np.zeros_like(x), ridge)
                if float(np.abs(x_new - x).max()) < CONVERGED_V:
                    break
                for alpha in (1.0, 0.5, 0.25, 0.125):
                    trial = x + alpha * (x_new - x)
                    a_t, r_t, profiles_t = linearize(trial)
                    value_t = objective(profiles_t, trial)
                    if value_t < value:
                        break
                else:
                    break  # no step decreases the relaxed objective
                x, a, r, profiles, value = trial, a_t, r_t, profiles_t, value_t
                history.append({"iteration": iteration, "objective": value, "step": alpha,
                                "max_abs_v": float(np.abs(x).max())})
            _, determined = solve(a, r, np.zeros_like(x), ridge)
        else:
            profiles = [scan["zero"] for scan in scans]
            determined = np.zeros(0)
        pool.shutdown()
        timings["torsion_fit_s"] = time.perf_counter() - clock
        fitted_torsions = model.rb(x)
        report["fit_history"] = history
        for scan, profile in zip(scans, profiles):
            mask = scan["mask"]

            def rmse(model_profile):
                return float(np.sqrt(np.mean((model_profile[mask] - scan["qm"][mask]) ** 2)))

            task = scan["task"]
            report["scans"].append({
                "bond": task["bond"], "ring": task["ring"],
                "drive": [d + 1 for d in task["drive"]], "points": len(scan["points"]),
                "fitted_points": int(mask.sum()), "qm_range_kj_mol": float(scan["qm"].max()),
                "rmse_kj_mol": {"missing_zero": rmse(scan["zero"]), "fitted": rmse(profile)},
            })
            with (workdir / f"profile_{task['bond'][0]}_{task['bond'][1]}.csv").open("w") as out:
                out.write("phi_deg,qm,missing_zero,fitted\n")
                for i, point in enumerate(scan["points"]):
                    out.write(f"{point['phi_deg']:.2f},{scan['qm'][i]:.4f},"
                              f"{scan['zero'][i]:.4f},{profile[i]:.4f}\n")
        for key, (start, order) in model.columns.items():
            report["torsion_keys"].append({
                "classes": "-".join(key[0]), "context": key[1],
                "fourier_v_kj_mol": [float(v) for v in x[start:start + order]],
                "determined": [float(v) for v in determined[start:start + order]]
                if len(determined) else [],
                "quartets": [[a + 1 for a in q] for q in model.carrier_quartets
                             if model.keys[q] == key],
            })
        report["zero_by_convention"] = [[a + 1 for a in q] for q in pending
                                        if model.keys[q] not in model.columns]
    report["zero_linear"] = [list(value.atoms) for value in linear_torsions]

    for stale in itp.parent.glob("*.qmfit_*.*"):
        stale.unlink()
    timings["total_s"] = time.perf_counter() - clock_start
    report["timings_s"] = {k: round(v, 1) for k, v in timings.items()}
    report["threads"] = threads
    report["mm_engine"] = mm_engine
    (workdir / "report.json").write_text(json.dumps(report, indent=2) + "\n")

    def replace_bond(value):
        if value.rule != RULE_QM_PENDING:
            return value
        b0, k = fitted_bonds[value.atoms]
        return dataclasses.replace(value, length_nm=b0, force_constant_kj_mol_nm2=k,
                                   rule=RULE_QM_FIT, source_line=0)

    def replace_angle(value):
        if value.rule != RULE_QM_PENDING:
            return value
        theta, k = fitted_angles[value.atoms]
        return dataclasses.replace(value, angle_deg=theta, force_constant_kj_mol_rad2=k,
                                   rule=RULE_QM_FIT, source_line=0)

    def replace_torsion(value):
        if value.rule != RULE_QM_PENDING:
            return value
        if value in linear_torsions:
            return dataclasses.replace(value, coefficients_kj_mol=(0.0,) * 6,
                                       rule=RULE_ZERO_LINEAR, source_line=0)
        coefficients = fitted_torsions.get(value.atoms) or fitted_torsions.get(value.atoms[::-1])
        return dataclasses.replace(value, coefficients_kj_mol=tuple(coefficients),
                                   rule=RULE_QM_FIT, source_line=0)

    warnings = list(bonded.warnings) + cache_warnings
    by_bond: dict[tuple[int, int], int] = {}
    for value in linear_torsions:
        central = tuple(sorted(value.atoms[1:3]))
        by_bond[central] = by_bond.get(central, 0) + 1
    for (j, k), count in sorted(by_bond.items()):
        warnings.append(f"bond {j}-{k}: {count} torsion(s) about a linear atom are missing "
                        f"from the table; set to zero ({RULE_ZERO_LINEAR}), not fitted")
    for scan in report["scans"]:
        if scan["rmse_kj_mol"]["fitted"] > 4.0:
            warnings.append(f"qm-fit scan about bond {scan['bond'][0]}-{scan['bond'][1]}: "
                            f"RMSE {scan['rmse_kj_mol']['fitted']:.1f} kJ/mol after fitting")
    return dataclasses.replace(
        bonded,
        bonds=tuple(replace_bond(v) for v in bonded.bonds),
        angles=tuple(replace_angle(v) for v in bonded.angles),
        proper_torsions=tuple(replace_torsion(v) for v in bonded.proper_torsions),
        warnings=tuple(warnings),
        qm_method=method_label(method, dft_level)
        + (f"; MM relaxations: {mm_engine}" if pending_torsions else ""),
    )
