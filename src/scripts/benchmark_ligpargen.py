#!/usr/bin/env python3
"""Benchmark opls_ff_gen ITP files against LigParGen ITP files.

Usage:
    python src/scripts/benchmark_ligpargen.py NAME [NAME ...]

For each NAME, compares ``benchmarks/NAME/NAME.itp`` (opls_ff_gen) with
``molecules_to_test/NAME.itp`` (LigParGen). Both files must use the same atom
order. Per-atom and per-interaction tables are written to ``benchmarks/NAME/``
and a summary table is printed in Markdown.
"""

from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

from rdkit import Chem


PROJECT_ROOT = Path(__file__).resolve().parents[2]
REFERENCE_DIR = PROJECT_ROOT / "molecules_to_test"
BENCHMARK_DIR = PROJECT_ROOT / "benchmarks"

# Printed precision of the LigParGen ITP files.
TOLERANCE = {
    "bond_b0": 5e-5,
    "bond_kb": 0.05,
    "angle_theta": 1e-3,
    "angle_k": 1e-3,
    "torsion": 1e-3,
    "improper": 1e-3,
}


def parse_itp(path: Path) -> dict:
    atomtypes, atoms = {}, []
    bonds, angles, propers, impropers, pairs = {}, {}, {}, {}, set()
    section = None
    for raw in path.read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("["):
            section = line.strip("[] ").lower()
            continue
        f = line.split()
        if section == "atomtypes" and len(f) >= 7:
            # LigParGen: name class mass charge ptype sigma eps
            # opls_ff_gen: name class at.num mass charge ptype sigma eps
            if f[4] == "A":
                atomtypes[f[0]] = (f[1], float(f[2]), float(f[5]), float(f[6]))
            else:
                atomtypes[f[0]] = (f[1], float(f[3]), float(f[6]), float(f[7]))
        elif section == "atoms" and len(f) >= 8:
            atoms.append({"type": f[1], "charge": float(f[6]), "mass": float(f[7])})
        elif section == "bonds" and len(f) >= 5:
            i, j = int(f[0]), int(f[1])
            bonds[tuple(sorted((i, j)))] = (float(f[3]), float(f[4]))
        elif section == "angles" and len(f) >= 6:
            i, j, k = int(f[0]), int(f[1]), int(f[2])
            angles[min((i, j, k), (k, j, i))] = (float(f[4]), float(f[5]))
        elif section == "pairs" and len(f) >= 2:
            pairs.add(tuple(sorted((int(f[0]), int(f[1])))))
        elif section == "dihedrals" and len(f) >= 5:
            idx = tuple(int(v) for v in f[:4])
            funct = int(f[4])
            if funct == 3:
                propers[min(idx, idx[::-1])] = tuple(float(v) for v in f[5:11])
            elif funct in (1, 4):
                center = idx[1]
                key = (center, frozenset(idx) - {center})
                impropers[key] = (float(f[5]), float(f[6]), int(float(f[7])))
    return {
        "atomtypes": atomtypes,
        "atoms": atoms,
        "bonds": bonds,
        "angles": angles,
        "propers": propers,
        "impropers": impropers,
        "pairs": pairs,
    }


def atom_classes(itp: dict) -> list[str]:
    return [itp["atomtypes"][atom["type"]][0] for atom in itp["atoms"]]


def compare(name: str) -> dict:
    out_dir = BENCHMARK_DIR / name
    ref = parse_itp(REFERENCE_DIR / f"{name}.itp")
    assignment = json.loads((out_dir / f"{name}.assignment.json").read_text())
    if len(ref["atoms"]) != len(assignment["atoms"]):
        raise RuntimeError(f"{name}: atom-count mismatch")
    molecule = Chem.SDMolSupplier(
        str(out_dir / f"{name}.assigned.sdf"), removeHs=False
    )[0]
    sdf_bonds = {
        tuple(sorted((b.GetBeginAtomIdx() + 1, b.GetEndAtomIdx() + 1)))
        for b in molecule.GetBonds()
    }
    if sdf_bonds != set(ref["bonds"]):
        raise RuntimeError(f"{name}: connectivity or atom order differs")
    ref_cls = atom_classes(ref)

    # Atoms: workflow values come from the assignment record, which exists
    # even when bonded assignment fails and no ITP is written.
    atom_rows = []
    for i, (r, atom) in enumerate(zip(ref["atoms"], assignment["atoms"]), start=1):
        rt = ref["atomtypes"][r["type"]]
        rs, re_ = rt[2], rt[3]
        gs, ge = atom["sigma_nm"], atom["epsilon_kj_mol"]
        lj_equal = (re_ == 0 and ge == 0) or (
            abs(rs - gs) < 5e-6 and abs(re_ - ge) < 5e-6
        )
        charge = atom["final_charge_e"]
        atom_rows.append(
            {
                "atom": i,
                "element": atom["element"],
                "foyer_type": atom["opls_type"],
                "ligpargen_class": rt[0],
                "workflow_class": atom["opls_class"],
                "class_match": rt[0] == atom["opls_class"],
                "ligpargen_sigma_nm": rs,
                "workflow_sigma_nm": gs,
                "ligpargen_epsilon_kj_mol": re_,
                "workflow_epsilon_kj_mol": ge,
                "lj_match": lj_equal,
                "ligpargen_charge_e": r["charge"],
                "workflow_charge_e": round(charge, 6),
                "charge_diff_e": round(charge - r["charge"], 6),
                "lbcc_delta_e": round(atom["lbcc_delta_e"], 6),
            }
        )
    write_csv(out_dir / "benchmark_atoms.csv", atom_rows)

    charges = [abs(r["charge_diff_e"]) for r in atom_rows]
    summary = {
        "molecule": name,
        "atoms": len(atom_rows),
        "class_matches": sum(r["class_match"] for r in atom_rows),
        "lj_matches": sum(r["lj_match"] for r in atom_rows),
        "charge_mae_e": sum(charges) / len(charges),
        "charge_max_e": max(charges),
        "lbcc_applications": len(assignment["lbcc_applications"]),
        "itp_generated": False,
    }
    # The workflow writes a sanitized, lower-case stem (e.g. m_1_octanol.itp).
    itp_candidates = sorted(out_dir.glob("*.itp"))
    itp_path = itp_candidates[0] if len(itp_candidates) == 1 else out_dir / f"{name}.itp"
    if not itp_path.is_file():
        (out_dir / "benchmark_summary.json").write_text(
            json.dumps(summary, indent=2) + "\n"
        )
        return summary
    summary["itp_generated"] = True
    gen = parse_itp(itp_path)

    # Interactions
    term_rows = []

    def classes(atoms):
        return "-".join(ref_cls[a - 1] for a in atoms)

    def record(kind, atoms, ref_value, gen_value, tolerances):
        if ref_value is None or gen_value is None:
            status = "only_ligpargen" if gen_value is None else "only_workflow"
        else:
            status = "match"
            for rv, gv, tol in zip(ref_value, gen_value, tolerances):
                if abs(rv - gv) > tol:
                    status = "differs"
        term_rows.append(
            {
                "kind": kind,
                "atoms": " ".join(str(a) for a in atoms),
                "ligpargen_classes": classes(atoms),
                "status": status,
                "ligpargen": " ".join(f"{v:g}" for v in ref_value) if ref_value else "",
                "workflow": " ".join(f"{v:g}" for v in gen_value) if gen_value else "",
            }
        )

    for key in sorted(set(ref["bonds"]) | set(gen["bonds"])):
        record("bond", key, ref["bonds"].get(key), gen["bonds"].get(key),
               (TOLERANCE["bond_b0"], TOLERANCE["bond_kb"]))
    for key in sorted(set(ref["angles"]) | set(gen["angles"])):
        record("angle", key, ref["angles"].get(key), gen["angles"].get(key),
               (TOLERANCE["angle_theta"], TOLERANCE["angle_k"]))
    for key in sorted(set(ref["propers"]) | set(gen["propers"])):
        record("proper", key, ref["propers"].get(key), gen["propers"].get(key),
               (TOLERANCE["torsion"],) * 6)
    for key in sorted(set(ref["impropers"]) | set(gen["impropers"]),
                      key=lambda k: k[0]):
        center, others = key
        atoms = (center, *sorted(others))
        record("improper", atoms, ref["impropers"].get(key),
               gen["impropers"].get(key), (TOLERANCE["improper"],) * 3)
    for key in sorted(ref["pairs"] | gen["pairs"]):
        record("pair", key, (1.0,) if key in ref["pairs"] else None,
               (1.0,) if key in gen["pairs"] else None, (0.0,))
    write_csv(out_dir / "benchmark_terms.csv", term_rows)

    # Grouped parameter differences by LigParGen class pattern
    families = defaultdict(lambda: {"count": 0, "ligpargen": set(), "workflow": set()})
    for row in term_rows:
        if row["status"] == "match" or row["kind"] == "pair":
            continue
        pattern = row["ligpargen_classes"]
        parts = pattern.split("-")
        if row["kind"] != "improper" and parts[::-1] < parts:
            pattern = "-".join(parts[::-1])
        fam = families[(row["kind"], row["status"], pattern)]
        fam["count"] += 1
        fam["ligpargen"].add(row["ligpargen"])
        fam["workflow"].add(row["workflow"])
    family_rows = [
        {
            "kind": kind,
            "status": status,
            "classes": pattern,
            "count": fam["count"],
            "ligpargen": " | ".join(sorted(fam["ligpargen"])),
            "workflow": " | ".join(sorted(fam["workflow"])),
        }
        for (kind, status, pattern), fam in sorted(families.items())
    ]
    write_csv(out_dir / "benchmark_term_families.csv", family_rows)

    for kind in ("bond", "angle", "proper", "improper", "pair"):
        rows = [r for r in term_rows if r["kind"] == kind]
        summary[f"{kind}_ligpargen"] = sum(r["status"] != "only_workflow" for r in rows)
        summary[f"{kind}_workflow"] = sum(r["status"] != "only_ligpargen" for r in rows)
        summary[f"{kind}_match"] = sum(r["status"] == "match" for r in rows)
        summary[f"{kind}_differs"] = sum(r["status"] == "differs" for r in rows)
    (out_dir / "benchmark_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main(names: list[str]) -> None:
    summaries = [compare(name) for name in names]
    print("| Molecule | Atoms | Class | LJ | Charge MAE (e) | Max |dq| (e) |"
          " Bonds | Angles | Propers | Impropers | Pairs |")
    print("|---|---:|---:|---:|---:|---:|---|---|---|---|---|")
    for s in summaries:
        if not s["itp_generated"]:
            print(
                f"| {s['molecule']} | {s['atoms']} | {s['class_matches']} |"
                f" {s['lj_matches']} | {s['charge_mae_e']:.4f} |"
                f" {s['charge_max_e']:.4f} | no ITP | | | | |"
            )
            continue
        terms = []
        for kind in ("bond", "angle", "proper", "improper", "pair"):
            terms.append(
                f"{s[f'{kind}_match']}/{s[f'{kind}_ligpargen']}"
                f" (wf {s[f'{kind}_workflow']})"
            )
        print(
            f"| {s['molecule']} | {s['atoms']} | {s['class_matches']} |"
            f" {s['lj_matches']} | {s['charge_mae_e']:.4f} |"
            f" {s['charge_max_e']:.4f} | " + " | ".join(terms) + " |"
        )


if __name__ == "__main__":
    main(sys.argv[1:])
