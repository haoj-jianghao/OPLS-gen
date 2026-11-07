"""Estimate the compute time of on-the-fly QM fitting for a batch.

Inputs:
- a batch directory from ``run_batch_opls2020.py`` (xTB run): per molecule the
  wall time, the QM report (scan count, points, which steps ran) and the atom
  count;
- a timing JSON of QM building blocks (``timings``: per molecule atom count and
  seconds for DFT single points, gradients and Hessians at given basis sets).

For the xTB level the batch's measured wall times are summed. For a DFT level,
the measured xTB-specific steps are replaced by DFT costs from power-law fits
t = a * N_atoms^b to the timing data:
  scan energies:  points x t_SP(N)
  geometry:       OPT_CYCLES x (t_SP + t_grad)   (geomeTRIC from the xTB minimum;
                  only when bonds/angles are fitted)
  Hessian:        t_hess(N) when bonds/angles were fitted
The scan geometries (xTB) and the MM torsion fit are the same at every level.

    python src/scripts/estimate_qm_cost.py BATCH_DIR TIMINGS.json [--cores 8]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

OPT_CYCLES = 15  # typical geomeTRIC cycles from an xTB-optimized start (checked on test molecules)


def power_law(sizes, seconds):
    b, log_a = np.polyfit(np.log(sizes), np.log(seconds), 1)
    return lambda n: float(np.exp(log_a) * n ** b), (float(np.exp(log_a)), float(b))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("batch", type=Path)
    parser.add_argument("timings", type=Path)
    args = parser.parse_args()
    timing = json.loads(args.timings.read_text())
    fits = {}
    for key in ("b3lyp_def2-svp_sp_s", "b3lyp_def2-svp_grad_s", "b3lyp_def2-tzvp_sp_s",
                "b3lyp_def2-tzvp_grad_s", "b3lyp_def2-svp_dfhess_s"):
        rows = [(v["atoms"], v[key]) for v in timing.values() if key in v]
        if len(rows) >= 2:
            fits[key] = power_law(*zip(*rows))
    for key, (_, (a, b)) in fits.items():
        print(f"fit {key}: t = {a:.3g} * N^{b:.2f} s")

    summary = json.loads((args.batch / "batch_summary.json").read_text())
    totals = {"xtb": 0.0, "dft_svp": 0.0, "dft_tzvp_sp": 0.0}
    rows = []
    for name, entry in summary.items():
        if not entry.get("generate_ok"):
            continue
        wall = entry["wall_s"]
        qm = entry.get("qm")
        totals["xtb"] += wall
        if not qm:
            totals["dft_svp"] += wall
            totals["dft_tzvp_sp"] += wall
            continue
        atoms = _atoms(args.batch / _slug(name))
        points = sum(24 if not s["ring"] else _ring_points(args.batch / _slug(name), s)
                     for s in qm["scans"])
        t = qm["timings_s"]
        xtb_specific = t.get("optimization_s", 0) + t.get("hessian_s", 0) + t.get("scan_energies_s", 0)
        common = wall - xtb_specific
        hess = bool(qm["bonds"] or qm["angles"])
        estimates = {}
        for label, sp_key, grad_key in (("dft_svp", "b3lyp_def2-svp_sp_s", "b3lyp_def2-svp_grad_s"),
                                        ("dft_tzvp_sp", "b3lyp_def2-tzvp_sp_s", "b3lyp_def2-tzvp_grad_s")):
            sp = fits[sp_key][0](atoms)
            grad = fits[grad_key][0](atoms)
            h = fits["b3lyp_def2-svp_dfhess_s"][0](atoms) if hess and "b3lyp_def2-svp_dfhess_s" in fits else 0.0
            optimization = OPT_CYCLES * (sp + grad) if hess else 0.0
            estimates[label] = common + points * sp + optimization + h
            totals[label] += estimates[label]
        rows.append((name, atoms, len(qm["scans"]), points, hess, wall, estimates["dft_svp"],
                     estimates["dft_tzvp_sp"]))
    rows.sort(key=lambda r: -r[6])
    print(f"\n{'molecule':24s} atoms scans points hess   xTB(s)  DFT/SVP(h) DFT/TZVP(h)")
    for r in rows:
        print(f"{r[0]:24s} {r[1]:5d} {r[2]:5d} {r[3]:6d} {str(r[4]):5s} {r[5]:8.0f} {r[6] / 3600:10.2f} {r[7] / 3600:10.2f}")
    print(f"\nmolecules: {sum(1 for e in summary.values() if e.get('generate_ok'))}, with QM fit: {len(rows)}")
    for label, seconds in totals.items():
        print(f"total {label}: {seconds / 3600:.2f} h")
    (args.batch / "qm_cost_estimate.json").write_text(json.dumps({
        "fits": {k: {"a": v[1][0], "b": v[1][1]} for k, v in fits.items()},
        "opt_cycles": OPT_CYCLES,
        "totals_h": {k: v / 3600 for k, v in totals.items()},
        "molecules": [dict(zip(("name", "atoms", "scans", "points", "hessian", "xtb_s",
                                "dft_svp_s", "dft_tzvp_sp_s"), r)) for r in rows],
    }, indent=2) + "\n")


def _slug(name: str) -> str:
    import re
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_").lower()


def _atoms(directory: Path) -> int:
    assignment = json.loads(next(directory.glob("*.assignment.json")).read_text())
    return assignment["atom_count"]


def _ring_points(directory: Path, scan: dict) -> int:
    path = directory / "qm_fit" / f"scan_{scan['bond'][0]}_{scan['bond'][1]}" / "scan.json"
    return len(json.loads(path.read_text())["points"])


if __name__ == "__main__":
    main()
