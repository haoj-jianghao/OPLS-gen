"""Compare the GAFF2 rule table with antechamber reference types (by LJ).

Reads a reference JSON from ``gaff2_reference.py``. For every atom, the rule
table's type and antechamber's type are compared on their Lennard-Jones
parameters (types with identical LJ are equivalent). Prints the agreement and
the mismatches grouped by (antechamber type -> rule type) with examples.

    python src/scripts/compare_gaff2_rules.py REF.json [--examples 3]
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rdkit import Chem, RDLogger  # noqa: E402

from opls_ff_gen.amber import gaff2_lennard_jones, gaff2_types  # noqa: E402

RDLogger.DisableLog("rdApp.*")


def lj_key(table, gaff_type):
    sigma, epsilon = table[gaff_type]
    return round(sigma, 6), round(epsilon, 6)


def compare(reference: dict, examples: int = 3) -> tuple[int, int, dict]:
    table = gaff2_lennard_jones()
    total = agree = 0
    mismatches: dict[tuple[str, str], list] = defaultdict(list)
    for entry in reference["molecules"]:
        molecule = Chem.AddHs(Chem.MolFromSmiles(entry["smiles"]))
        if molecule.GetNumAtoms() != len(entry["types"]):
            continue
        try:
            ours = gaff2_types(molecule)
        except Exception as exc:  # noqa: BLE001
            for i, ref in enumerate(entry["types"]):
                mismatches[(ref, "ERROR")].append((entry["name"], entry["smiles"], i + 1))
            total += len(entry["types"])
            continue
        for i, ref in enumerate(entry["types"]):
            if ref not in table:
                continue
            total += 1
            if lj_key(table, ours[i]) == lj_key(table, ref):
                agree += 1
            else:
                mismatches[(ref, ours[i])].append((entry["name"], entry["smiles"], i + 1))
    return total, agree, mismatches


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("--examples", type=int, default=3)
    args = parser.parse_args()
    reference = json.loads(args.reference.read_text())
    total, agree, mismatches = compare(reference)
    print(f"atoms {total}; LJ agreement {agree} ({100 * agree / total:.3f}%); "
          f"mismatched atoms {total - agree} in {len(mismatches)} type pairs")
    for (ref, ours), cases in sorted(mismatches.items(), key=lambda kv: -len(kv[1])):
        print(f"  antechamber {ref:3s} -> rules {ours:5s}: {len(cases)} atoms")
        for name, smiles, atom in cases[: args.examples]:
            print(f"      {name}: atom {atom} in {smiles}")


if __name__ == "__main__":
    main()
