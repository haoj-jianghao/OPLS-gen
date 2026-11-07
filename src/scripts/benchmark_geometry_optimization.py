"""Compare MMFF94 and GFN2-xTB geometry optimization for SMILES inputs.

For each molecule: one ETKDGv3 conformer (fixed seed), optimized separately
with MMFF94 and with GFN2-xTB (--opt tight, 1 thread); reports the wall time
of each optimization, the heavy-atom RMSD between the two results, and the
largest difference in raw CM1A charges (MOPAC AM1) computed on each geometry.
Also checks that the mmff path reproduces assignment.molecule_from_smiles.

    python src/scripts/benchmark_geometry_optimization.py CSV [CSV ...] --out report.json
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import rdMolAlign  # noqa: E402

from opls_ff_gen.assignment import molecule_from_smiles  # noqa: E402
from opls_ff_gen.cm1a import run_cm1a  # noqa: E402
from opls_ff_gen.structure import molecule_from_smiles_optimized  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", nargs="+", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for path in args.csv:
        with path.open() as handle:
            rows += list(csv.DictReader(handle))
    results = []
    for row in rows:
        entry = {"name": row["name"]}
        try:
            old = molecule_from_smiles(row["smiles"], name=row["name"])
            mmff = molecule_from_smiles_optimized(row["smiles"], row["name"], "mmff")
            xtb = molecule_from_smiles_optimized(row["smiles"], row["name"], "xtb")
        except Exception as exc:  # noqa: BLE001
            entry["error"] = f"{type(exc).__name__}: {exc}"[:300]
            results.append(entry)
            print(entry, flush=True)
            continue
        entry["atoms"] = mmff.GetNumAtoms()
        entry["mmff_reproduces_old"] = bool(np.array_equal(
            old.GetConformer().GetPositions(), mmff.GetConformer().GetPositions()))
        entry["mmff_s"] = float(mmff.GetProp("Geometry_Optimization_Seconds"))
        entry["xtb_s"] = float(xtb.GetProp("Geometry_Optimization_Seconds"))
        heavy_m, heavy_x = Chem.RemoveHs(mmff), Chem.RemoveHs(xtb)
        entry["heavy_rmsd_A"] = float(rdMolAlign.GetBestRMS(heavy_m, heavy_x))
        with tempfile.TemporaryDirectory() as scratch:
            q_m = run_cm1a(mmff, Path(scratch) / "m", "m", Chem.GetFormalCharge(mmff)).charges
            q_x = run_cm1a(xtb, Path(scratch) / "x", "x", Chem.GetFormalCharge(xtb)).charges
        entry["max_dq_cm1a"] = float(np.abs(q_m - q_x).max())
        entry["mean_dq_cm1a"] = float(np.abs(q_m - q_x).mean())
        results.append(entry)
        print(f"{row['name']:24s} atoms {entry['atoms']:3d}  MMFF {entry['mmff_s']:6.3f} s  "
              f"xTB {entry['xtb_s']:6.2f} s  RMSD {entry['heavy_rmsd_A']:.3f} A  "
              f"max|dq| {entry['max_dq_cm1a']:.3f} e", flush=True)
    ok = [r for r in results if "error" not in r]
    mm, xt = np.array([r["mmff_s"] for r in ok]), np.array([r["xtb_s"] for r in ok])
    summary = {
        "molecules": len(results), "optimized": len(ok),
        "errors": [r for r in results if "error" in r],
        "mmff_reproduces_old_path": sum(r["mmff_reproduces_old"] for r in ok),
        "mmff_seconds": {"total": float(mm.sum()), "mean": float(mm.mean()), "max": float(mm.max())},
        "xtb_seconds": {"total": float(xt.sum()), "mean": float(xt.mean()), "max": float(xt.max())},
        "xtb_over_mmff_mean_ratio": float(xt.mean() / mm.mean()),
        "heavy_rmsd_A": {"median": float(np.median([r["heavy_rmsd_A"] for r in ok])),
                         "max": float(max(r["heavy_rmsd_A"] for r in ok))},
        "cm1a_max_dq_e": {"median": float(np.median([r["max_dq_cm1a"] for r in ok])),
                          "max": float(max(r["max_dq_cm1a"] for r in ok))},
    }
    args.out.write_text(json.dumps({"summary": summary, "molecules": results}, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
