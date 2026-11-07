"""Compare torsions from the GAFF2 general terms / parmchk2 substitution with QM.

Two batch runs of the same molecules (``src/scripts/run_batch_opls2020.py``):

- ``--qm-batch``: ``--torsion-substitution none``, so every torsion missing
  from OPLS/2020 is scanned (relaxed xTB) and fitted;
- ``--general-batch``: the default, where those torsions take the GAFF2
  general term (or a parmchk2 substitution).

For every scan of the QM batch, the QM-fit topology (fitted bonds and angles,
table torsions) is rewritten with the formerly missing torsions set to the
values of the general batch, and relaxed along the scan exactly as the fit
does (``mm.mm_profile``, driven dihedral restrained, other rotors locked).
RMSE against QM over points within ``FIT_CAP_KJ_MOL`` of the QM minimum,
profiles relative to their minima, as in ``qm_fit/report.json``:

- ``zero``: missing torsions zero (from the QM batch report);
- ``general``: missing torsions from the general batch;
- ``fitted``: the on-the-fly fit (from the QM batch report).

    python src/scripts/compare_general_torsions.py --qm-batch A --general-batch B \\
        --out report.json [--workers 8] [--tag general]

``kept_fit`` counts the missing quartets about the scanned bond that the
general batch still fits (e.g. no parmchk2 substitution found); they keep
their fitted values, so such a scan is only partly a test of the fallback.
``--tag`` names the rewritten topology in the QM batch directory.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rdkit import Chem, RDLogger  # noqa: E402

from opls_ff_gen.mm import mm_profile  # noqa: E402
from opls_ff_gen.qm import FIT_CAP_KJ_MOL  # noqa: E402
from opls_ff_gen.qm_fit import _rewrite_itp  # noqa: E402
from opls_ff_gen.xtb import read_xyz  # noqa: E402

RDLogger.DisableLog("rdApp.*")
TABLE_RULES = {"exact", "wildcard"}


def _topology(directory: Path) -> dict:
    return json.loads(next(directory.glob("*.topology.json")).read_text())


def _scan_job(job: dict) -> dict:
    scan_dir = Path(job["scan_dir"])
    scan = json.loads((scan_dir / "scan.json").read_text())
    # qm_fit sorts the points by angle (profile CSV order); scan.json is in scan order
    points = [{"phi_deg": p["phi_deg"],
               "coordinates": read_xyz((Path(p["dir"]) / "xtbopt.xyz").read_text())}
              for p in sorted(scan["points"], key=lambda p: p["phi_deg"])]
    molecule = Chem.AddHs(Chem.MolFromSmiles(job["smiles"]))
    energies = np.asarray(mm_profile(Path(job["top"]), points, tuple(scan["drive"]), molecule,
                                     lock_hydrogens=True))
    with Path(job["profile_csv"]).open() as handle:
        rows = list(csv.DictReader(handle))
    assert all(abs(float(r["phi_deg"]) - p["phi_deg"]) < 0.01 for r, p in zip(rows, points))
    qm = np.array([float(r["qm"]) for r in rows])
    zero = np.array([float(r["missing_zero"]) for r in rows])
    fitted = np.array([float(r["fitted"]) for r in rows])
    general = energies - energies.min()
    mask = qm <= FIT_CAP_KJ_MOL

    def rmse(profile):
        return round(float(np.sqrt(np.mean((profile[mask] - qm[mask]) ** 2))), 3)

    return {**{k: job[k] for k in ("molecule", "bond", "central", "terms", "qm_check",
                                   "kept_fit")},
            "qm_range": round(float(qm.max()), 2),
            "rmse": {"zero": rmse(zero), "general": rmse(general), "fitted": rmse(fitted)},
            "profile": {"phi": [round(float(r["phi_deg"]), 2) for r in rows],
                        "qm": qm.round(3).tolist(), "general": general.round(3).tolist()}}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qm-batch", type=Path, required=True)
    parser.add_argument("--general-batch", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--tag", default="general")
    args = parser.parse_args()
    summary = json.loads((args.qm_batch / "batch_summary.json").read_text())

    jobs = []
    for name, entry in summary.items():
        stem = re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_").lower()  # run_batch slug
        qm_dir, general_dir = args.qm_batch / stem, args.general_batch / stem
        report_path = qm_dir / "qm_fit" / "report.json"
        if not report_path.is_file() or not general_dir.is_dir():
            continue
        report = json.loads(report_path.read_text())
        if not report["scans"]:
            continue
        qm_top, general_top = _topology(qm_dir), _topology(general_dir)
        general_terms = {tuple(t["atoms"]): t for t in general_top["proper_torsions"]}
        swap, terms_by_bond = {}, defaultdict(lambda: defaultdict(int))
        flagged_by_bond, kept_by_bond = defaultdict(bool), defaultdict(int)
        for torsion in qm_top["proper_torsions"]:
            if torsion["rule"] in TABLE_RULES:
                continue
            atoms = tuple(torsion["atoms"])
            other = general_terms[atoms]
            if other["rule"] in TABLE_RULES or other["rule"].startswith("qm"):
                kept_by_bond[frozenset(atoms[1:3])] += 1
                continue  # still fitted in the general batch: keep the fit
            swap[atoms] = tuple(other["coefficients_kj_mol"])
            bond = frozenset(atoms[1:3])
            terms_by_bond[bond][f"{other['rule']} {'-'.join(other['parameter_classes'])}"] += 1
            flagged_by_bond[bond] |= bool(other.get("qm_check"))
        if not swap:
            continue
        itp, top = Path(qm_top["files"]["itp"]), Path(qm_top["files"]["top"])
        _, new_top = _rewrite_itp(itp, top, {}, {}, swap, args.tag)
        assignment = json.loads(Path(qm_top["files"]["assignment"]).read_text())
        atom_classes = {a["index"]: a["opls_class"] for a in assignment["atoms"]}
        for scan in report["scans"]:
            j, k = scan["bond"]
            bond = frozenset((j, k))
            if bond not in terms_by_bond:
                continue
            jobs.append({
                "molecule": name, "smiles": entry["smiles"], "bond": [j, k],
                "central": "-".join(sorted((atom_classes[j], atom_classes[k]))),
                "terms": dict(terms_by_bond[bond]), "qm_check": flagged_by_bond[bond],
                "kept_fit": kept_by_bond[bond],
                "scan_dir": str(qm_dir / "qm_fit" / f"scan_{j}_{k}"),
                "profile_csv": str(qm_dir / "qm_fit" / f"profile_{j}_{k}.csv"),
                "top": str(new_top),
            })

    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(_scan_job, jobs))

    def stats(key, subset):
        values = [r["rmse"][key] for r in subset]
        return {"mean": round(statistics.mean(values), 2),
                "median": round(statistics.median(values), 2)} if values else None

    groups = {"all": results,
              "qm_check": [r for r in results if r["qm_check"]],
              "unflagged": [r for r in results if not r["qm_check"]]}
    by_term = defaultdict(list)
    for r in results:
        for term in r["terms"]:
            by_term[term].append(r)
    out = {
        "scans": len(results),
        "summary": {name: {"scans": len(subset),
                           **{key: stats(key, subset) for key in ("zero", "general", "fitted")},
                           "general_worse_than_zero": sum(
                               r["rmse"]["general"] > r["rmse"]["zero"] + 0.5 for r in subset)}
                    for name, subset in groups.items()},
        "by_term": {term: {"scans": len(rs),
                           **{key: stats(key, rs) for key in ("zero", "general", "fitted")}}
                    for term, rs in sorted(by_term.items(), key=lambda kv: -len(kv[1]))},
        "results": sorted(results, key=lambda r: -r["rmse"]["general"]),
    }
    args.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out["summary"], indent=1))


if __name__ == "__main__":
    main()
