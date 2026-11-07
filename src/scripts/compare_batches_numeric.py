"""Numerical comparison of two run_batch_opls2020.py batches (e.g. two charge engines).

Per molecule (from the topology/assignment JSON and the QM-fit report):
- atom types, classes and type rules must be identical;
- final charges: largest per-atom difference;
- every bonded term must have the same atoms and rule; table terms must have
  identical parameters; for QM-fitted terms the largest parameter change is
  reported (bond b0 nm / k, angle theta0 deg / k, torsion RB coefficients);
- torsion scans: fitted RMSE before/after.

    python src/scripts/compare_batches_numeric.py REFERENCE_DIR NEW_DIR [--out report.json]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def _load(directory: Path):
    topology = json.loads(next(directory.glob("*.topology.json")).read_text())
    assignment = json.loads(next(directory.glob("*.assignment.json")).read_text())
    report_path = directory / "qm_fit" / "report.json"
    report = json.loads(report_path.read_text()) if report_path.exists() else None
    return topology, assignment, report


def compare(ref_dir: Path, new_dir: Path) -> dict:
    t1, a1, r1 = _load(ref_dir)
    t2, a2, r2 = _load(new_dir)
    out = {"problems": []}
    key = lambda a: (a["opls_type"], a["opls_class"], a["type_rule"])  # noqa: E731
    if [key(a) for a in a1["atoms"]] != [key(a) for a in a2["atoms"]]:
        out["problems"].append("atom types differ")
    q1 = np.array([a["final_charge_e"] for a in a1["atoms"]])
    q2 = np.array([a["final_charge_e"] for a in a2["atoms"]])
    out["max_dq_final"] = float(np.abs(q1 - q2).max())
    raw1 = np.array([a["cm1a_raw_e"] for a in a1["atoms"]])
    raw2 = np.array([a["cm1a_raw_e"] for a in a2["atoms"]])
    out["max_dq_raw"] = float(np.abs(raw1 - raw2).max())
    fitted_change = {"bond_b0": 0.0, "bond_k_rel": 0.0, "angle_theta": 0.0, "angle_k_rel": 0.0,
                     "torsion_rb": 0.0}
    fields = {"bonds": ("length_nm", "force_constant_kj_mol_nm2"),
              "angles": ("angle_deg", "force_constant_kj_mol_rad2"),
              "proper_torsions": ("coefficients_kj_mol",),
              "improper_torsions": ("phase_deg", "force_constant_kj_mol", "multiplicity")}
    for group, names in fields.items():
        g1, g2 = t1[group], t2[group]
        if [(v["atoms"], v["rule"]) for v in g1] != [(v["atoms"], v["rule"]) for v in g2]:
            out["problems"].append(f"{group}: different terms or rules")
            continue
        for v1, v2 in zip(g1, g2):
            values1 = np.ravel([v1[n] for n in names]).astype(float)
            values2 = np.ravel([v2[n] for n in names]).astype(float)
            if v1["rule"] != "qm-fit":
                if not np.array_equal(values1, values2):
                    out["problems"].append(f"{group} {v1['atoms']}: table parameters differ")
                continue
            if group == "bonds":
                fitted_change["bond_b0"] = max(fitted_change["bond_b0"], abs(values1[0] - values2[0]))
                fitted_change["bond_k_rel"] = max(fitted_change["bond_k_rel"],
                                                  abs(values1[1] - values2[1]) / values1[1])
            elif group == "angles":
                fitted_change["angle_theta"] = max(fitted_change["angle_theta"],
                                                   abs(values1[0] - values2[0]))
                fitted_change["angle_k_rel"] = max(fitted_change["angle_k_rel"],
                                                   abs(values1[1] - values2[1]) / values1[1])
            else:
                fitted_change["torsion_rb"] = max(fitted_change["torsion_rb"],
                                                  float(np.abs(values1 - values2).max()))
    out["fitted_change"] = fitted_change
    if (r1 is None) != (r2 is None):
        out["problems"].append("QM fit in only one batch")
    elif r1 is not None:
        out["scan_rmse"] = [(s1["bond"], round(s1["rmse_kj_mol"]["fitted"], 3),
                             round(s2["rmse_kj_mol"]["fitted"], 3))
                            for s1, s2 in zip(r1["scans"], r2["scans"])]
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("new", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    names = sorted({p.name for p in args.reference.iterdir() if p.is_dir()}
                   & {p.name for p in args.new.iterdir() if p.is_dir()})
    results = {name: compare(args.reference / name, args.new / name) for name in names}
    with_problems = {n: r["problems"] for n, r in results.items() if r["problems"]}
    dq = np.array([r["max_dq_final"] for r in results.values()])
    dq_raw = np.array([r["max_dq_raw"] for r in results.values()])
    torsion = np.array([r["fitted_change"]["torsion_rb"] for r in results.values()])
    rmse = [(n, b, a, c) for n, r in results.items() for b, a, c in r.get("scan_rmse", [])]
    rmse_delta = np.array([c - a for _, _, a, c in rmse]) if rmse else np.zeros(0)
    summary = {
        "molecules": len(names), "molecules_with_problems": len(with_problems),
        "max_dq_raw_e": float(dq_raw.max()), "max_dq_final_e": float(dq.max()),
        "median_molecule_max_dq_final_e": float(np.median(dq)),
        "max_fitted_torsion_rb_change_kj_mol": float(torsion.max()),
        "max_fitted_bond_b0_change_nm": max(r["fitted_change"]["bond_b0"] for r in results.values()),
        "max_fitted_bond_k_relative_change": max(r["fitted_change"]["bond_k_rel"] for r in results.values()),
        "max_fitted_angle_theta_change_deg": max(r["fitted_change"]["angle_theta"] for r in results.values()),
        "max_fitted_angle_k_relative_change": max(r["fitted_change"]["angle_k_rel"] for r in results.values()),
        "scans": len(rmse),
        "scan_rmse_change_kj_mol": {"mean": float(rmse_delta.mean()) if len(rmse) else None,
                                    "max_abs": float(np.abs(rmse_delta).max()) if len(rmse) else None},
    }
    print(json.dumps(summary, indent=1))
    for name, problems in with_problems.items():
        print("  PROBLEM", name, problems)
    if args.out:
        args.out.write_text(json.dumps({"summary": summary, "problems": with_problems,
                                        "molecules": results}, indent=1) + "\n")


if __name__ == "__main__":
    main()
