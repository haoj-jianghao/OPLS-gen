"""OPLS/2020 coverage of a molecule CSV: typing and missing bonded terms.

For every molecule: atom typing (Foyer + archive + OPLS/2020 refinement, GAFF2
for atoms OPLS cannot type) and the bonded lookup in the OPLS/2020 tables with
missing terms deferred. Reports per molecule the GAFF2 atoms, the missing
bonds/angles/torsions, and the scans the on-the-fly QM fit would need
(unique central bonds after graph symmetry, rotatable vs ring).

    python src/scripts/coverage_opls2020.py CSV [CSV ...] --out OUT.json
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rdkit import Chem  # noqa: E402

from opls_ff_gen.assignment import _assign_atom_types, molecule_from_smiles  # noqa: E402
from opls_ff_gen.bonded import RULE_QM_PENDING, assign_bonded_interactions  # noqa: E402


def scan_bonds(molecule: Chem.Mol, torsions) -> tuple[int, int]:
    """Unique (symmetry-reduced) central bonds with pending torsions: rotatable, ring."""
    ranks = list(Chem.CanonicalRankAtoms(molecule, breakTies=False))
    seen, rotatable, ring = set(), 0, 0
    for torsion in torsions:
        j, k = (a - 1 for a in torsion.atoms[1:3])
        key = tuple(sorted((ranks[j], ranks[k])))
        if key in seen:
            continue
        seen.add(key)
        if molecule.GetBondBetweenAtoms(j, k).IsInRing():
            ring += 1
        else:
            rotatable += 1
    return rotatable, ring


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", nargs="+", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for path in args.csv:
        with path.open() as handle:
            for record in csv.DictReader(handle):
                rows.append((path.parent.name, record["name"], record["smiles"]))
    results, missing_classes = [], Counter()
    for batch, name, smiles in rows:
        entry = {"batch": batch, "name": name, "smiles": smiles}
        try:
            molecule = molecule_from_smiles(smiles, name=name)
            types, parameters, rules = _assign_atom_types(molecule)
            classes = [parameters[t]["class"] for t in types]
            bonded = assign_bonded_interactions(molecule, classes, defer_missing=True)
        except Exception as exc:  # noqa: BLE001 - reported per molecule
            entry["error"] = f"{type(exc).__name__}: {exc}"
            results.append(entry)
            print(f"{name}: ERROR {entry['error'][:200]}")
            continue
        pending = {kind: [v for v in group if v.rule == RULE_QM_PENDING]
                   for kind, group in (("bond", bonded.bonds), ("angle", bonded.angles),
                                       ("proper", bonded.proper_torsions))}
        for kind, values in pending.items():
            for value in values:
                missing_classes[(kind, "-".join(min(value.classes, value.classes[::-1])))] += 1
        rotatable, ring = scan_bonds(molecule, pending["proper"])
        entry.update(
            atoms=molecule.GetNumAtoms(), heavy=molecule.GetNumHeavyAtoms(),
            gaff2_atoms=[f"{i + 1}:{r.split(':')[-1]}" for i, r in enumerate(rules)
                         if r.startswith("amber")],
            borrowed_atoms=sum(r.startswith("borrowed") for r in rules),
            pending={kind: len(values) for kind, values in pending.items()},
            total={"bond": len(bonded.bonds), "angle": len(bonded.angles),
                   "proper": len(bonded.proper_torsions)},
            scan_bonds={"rotatable": rotatable, "ring": ring},
            pending_classes=sorted({f"{k}:{'-'.join(v.classes)}" for k, vs in pending.items()
                                    for v in vs}),
        )
        results.append(entry)
        print(f"{name}: {entry['atoms']} atoms, gaff2 {len(entry['gaff2_atoms'])}, pending "
              f"{entry['pending']}, scans {entry['scan_bonds']}")
    args.out.write_text(json.dumps({
        "molecules": results,
        "missing_class_counts": [{"kind": k, "classes": c, "count": n}
                                 for (k, c), n in missing_classes.most_common()],
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
