"""Benchmark: OPLS/2020 table torsions against relaxed xTB scans.

Reference for the QM comparisons of fitted or borrowed torsions: how close do
torsions taken entirely from the OPLS/2020 table come to GFN2-xTB? Molecules
of a batch run (``src/scripts/run_batch_opls2020.py``) whose bonded terms all come
from the table (no QM fit) are scanned with the protocol of the on-the-fly
fit (``qm_fit``): the same drive dihedral (``_drive``), relaxed xTB scan
(``qm.relaxed_scan``: 15 degree steps over 360 degrees, ring bonds over the
ring-internal dihedral), MM relaxed with the driven dihedral restrained and
every other rotor locked (``mm.mm_profile``), RMSE over points within
``FIT_CAP_KJ_MOL`` of the QM minimum, profiles relative to their minima.

Scanned: every bond between two non-terminal, non-linear atoms, one per set of
symmetry-equivalent bonds; ``kind`` tells rotors with heavy atoms on both
sides (``heavy``), rotors with only hydrogens on one side (``h-rotor``, e.g.
CH3, OH) and ring bonds (``ring``).

    python src/scripts/benchmark_table_torsions.py --batch DIR --out DIR \\
        [--threads 8] [--dry-run]
"""

from __future__ import annotations

import argparse
import json
import re
import statistics
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rdkit import Chem, RDLogger  # noqa: E402

from opls_ff_gen.mm import mm_profile  # noqa: E402
from opls_ff_gen.qm import FIT_CAP_KJ_MOL, relaxed_scan  # noqa: E402
from opls_ff_gen.qm_fit import _drive, _linear, _scan_energies  # noqa: E402

RDLogger.DisableLog("rdApp.*")
TABLE_RULES = {"exact", "wildcard", "improper-family"}


def _bonds(molecule: Chem.Mol) -> list[tuple[int, int, str]]:
    ranks = list(Chem.CanonicalRankAtoms(molecule, breakTies=False))
    seen, bonds = set(), []
    for bond in molecule.GetBonds():
        j, k = sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        a, b = molecule.GetAtomWithIdx(j), molecule.GetAtomWithIdx(k)
        if a.GetDegree() < 2 or b.GetDegree() < 2 or _linear(a) or _linear(b):
            continue
        key = tuple(sorted((ranks[j], ranks[k])))
        if key in seen:
            continue
        seen.add(key)
        if bond.IsInRing():
            kind = "ring"
        elif all(any(n.GetAtomicNum() > 1 for n in atom.GetNeighbors() if n.GetIdx() != other)
                 for atom, other in ((a, k), (b, j))):
            kind = "heavy"
        else:
            kind = "h-rotor"
        bonds.append((j, k, kind))
    return bonds


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--dry-run", action="store_true", help="list molecules and bonds only")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    summary = json.loads((args.batch / "batch_summary.json").read_text())

    molecules = []
    for name, entry in summary.items():
        rules = entry.get("rule_counts", {})
        if not entry.get("mdrun_ok") or set(rules) - TABLE_RULES:
            continue  # QM-fitted, GAFF2 impropers or failed: not table-only
        stem = re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_").lower()
        directory = args.batch / stem
        topology = json.loads(next(directory.glob("*.topology.json")).read_text())
        molecule = Chem.MolFromMolFile(str(next(directory.glob("*.assigned.sdf"))), removeHs=False)
        molecules.append((name, directory, topology, molecule, _bonds(molecule)))
    kinds = defaultdict(int)
    for *_, bonds in molecules:
        for *_, kind in bonds:
            kinds[kind] += 1
    print(f"{len(molecules)} table-only molecules, {sum(kinds.values())} bonds {dict(kinds)}")
    if args.dry_run:
        for name, _, _, molecule, bonds in molecules:
            print(f"  {name}: {molecule.GetNumAtoms()} atoms, {len(bonds)} bonds")
        return

    results = []
    for name, directory, topology, molecule, bonds in molecules:
        clock = time.perf_counter()
        elements = [atom.GetSymbol() for atom in molecule.GetAtoms()]
        charge = Chem.GetFormalCharge(molecule)
        start = molecule.GetConformer().GetPositions()
        workdir = args.out / directory.name
        tasks = [{"bond": [j + 1, k + 1], "ring": kind == "ring", "kind": kind,
                  "drive": list(_drive(molecule, j, k)), "elements": elements,
                  "start": start.tolist(), "charge": charge,
                  "scan_dir": str(workdir / f"scan_{j + 1}_{k + 1}")} for j, k, kind in bonds]
        with ProcessPoolExecutor(max_workers=max(1, min(args.threads, len(tasks)))) as pool:
            list(pool.map(relaxed_scan, tasks))
        scans = _scan_energies(tasks, elements, "xtb", charge, args.threads, None)
        top = Path(topology["files"]["top"])
        with ProcessPoolExecutor(max_workers=max(1, min(args.threads, len(tasks)))) as pool:
            profiles = list(pool.map(
                mm_profile, [top] * len(tasks), scans, [tuple(t["drive"]) for t in tasks],
                [molecule] * len(tasks), [False] * len(tasks), [True] * len(tasks)))
        assignment = json.loads(Path(topology["files"]["assignment"]).read_text())
        atom_classes = [a["opls_class"] for a in assignment["atoms"]]
        for task, points, energies in zip(tasks, scans, profiles):
            qm = np.array([p["qm_kj_mol"] for p in points])
            qm -= qm.min()
            mm = np.asarray(energies) - min(energies)
            mask = qm <= FIT_CAP_KJ_MOL
            j, k = (a - 1 for a in task["bond"])
            results.append({
                "molecule": name, "bond": task["bond"], "kind": task["kind"],
                "central": "-".join(sorted((atom_classes[j], atom_classes[k]))),
                "qm_range": round(float(qm.max()), 2), "points": len(points),
                "rmse": round(float(np.sqrt(np.mean((mm[mask] - qm[mask]) ** 2))), 3),
                "profile": {"phi": [round(p["phi_deg"], 2) for p in points],
                            "qm": qm.round(3).tolist(), "mm": mm.round(3).tolist()},
            })
        print(f"{name}: {len(tasks)} scans, {time.perf_counter() - clock:.0f} s", flush=True)

    def stats(subset):
        values = [r["rmse"] for r in subset]
        return {"scans": len(values), "mean": round(statistics.mean(values), 2),
                "median": round(statistics.median(values), 2)} if values else None

    by_central = defaultdict(list)
    for r in results:
        by_central[r["central"]].append(r)
    report = {
        "batch": str(args.batch), "molecules": len(molecules),
        "summary": {"all": stats(results),
                    **{kind: stats([r for r in results if r["kind"] == kind])
                       for kind in ("heavy", "h-rotor", "ring")}},
        "by_central": {c: stats(rs) for c, rs in sorted(by_central.items(), key=lambda kv: -len(kv[1]))},
        "results": sorted(results, key=lambda r: -r["rmse"]),
    }
    (args.out / "benchmark.json").write_text(json.dumps(report, indent=1) + "\n")
    print(json.dumps(report["summary"], indent=1))


if __name__ == "__main__":
    main()
