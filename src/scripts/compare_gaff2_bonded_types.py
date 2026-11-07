"""Compare the full GAFF2 bonded types (``gaff2_bonded.py``) with antechamber.

Reads a reference JSON from ``gaff2_reference.py`` (SMILES and antechamber
GAFF2 types in RDKit atom order). Each molecule gets the workflow's MMFF94
aromaticity. Reports

- atoms: our type vs antechamber's, with the two names of a conjugated pair
  merged (which one an atom gets depends only on the alternation);
- bonds and angles: the gaff2.dat term looked up with our types vs with
  antechamber's, compared by value (so a pair swap that gives the same term
  counts as agreement); "absent" = gaff2.dat has no term for our types.

    python src/scripts/compare_gaff2_bonded_types.py REF.json [--examples 2] [--json OUT]
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rdkit import Chem, RDLogger  # noqa: E402

from opls_ff_gen.gaff2_bonded import load_gaff2_bonded  # noqa: E402
from opls_ff_gen.structure import apply_typing_aromaticity  # noqa: E402

RDLogger.DisableLog("rdApp.*")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("--examples", type=int, default=2)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    table = load_gaff2_bonded()
    merge = {second: first for first, second in table.pairs.items()}
    norm = lambda t: merge.get(t, t)  # noqa: E731
    counts = Counter()
    atom_diff, term_diff = defaultdict(list), defaultdict(list)
    for entry in json.loads(args.reference.read_text())["molecules"]:
        types = entry["types"].split() if isinstance(entry["types"], str) else entry["types"]
        molecule = Chem.MolFromSmiles(entry["smiles"])
        if molecule is None:
            continue
        molecule = apply_typing_aromaticity(Chem.AddHs(molecule))
        if molecule.GetNumAtoms() != len(types):
            counts["skipped"] += 1
            continue
        ours = table.types(molecule)
        if None in ours:
            counts["untypable"] += 1
            continue
        counts["molecules"] += 1
        for i, (ref, our) in enumerate(zip(types, ours)):
            counts["atoms"] += 1
            if norm(ref) != norm(our):
                atom_diff[(norm(ref), norm(our))].append(f"{entry['name']}#{i + 1}")
        terms = [("bond", (b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in molecule.GetBonds()]
        for center in range(molecule.GetNumAtoms()):
            neighbors = [a.GetIdx() for a in molecule.GetAtomWithIdx(center).GetNeighbors()]
            terms += [("angle", (a, center, c)) for a, c in combinations(neighbors, 2)]
        for kind, atoms in terms:
            lookup = table.bond if kind == "bond" else table.angle
            expected = lookup(tuple(types[i] for i in atoms))
            got = lookup(tuple(ours[i] for i in atoms))
            counts[kind + "s"] += 1
            if expected is None:
                counts[kind + "s_absent_for_antechamber"] += 1
                continue
            if got is not None and got.values == expected.values:
                continue
            key = (kind, "-".join(expected.classes),
                   "absent " + "-".join(ours[i] for i in atoms) if got is None
                   else "-".join(got.classes))
            term_diff[key].append(entry["name"])
    counts["atoms_differ"] = sum(len(v) for v in atom_diff.values())
    for kind in ("bond", "angle"):
        counts[kind + "s_differ"] = sum(len(v) for k, v in term_diff.items() if k[0] == kind)
    print(dict(counts))
    for (ref, our), where in sorted(atom_diff.items(), key=lambda kv: -len(kv[1])):
        print(f"atom {ref} -> {our}: {len(where)}  {where[:args.examples]}")
    for key, where in sorted(term_diff.items(), key=lambda kv: -len(kv[1]))[:40]:
        print(f"{key[0]} {key[1]} -> {key[2]}: {len(where)}  {where[:args.examples]}")
    if args.json:
        args.json.write_text(json.dumps({
            "counts": counts,
            "atoms": [{"antechamber": r, "ours": o, "count": len(w), "examples": w[:5]}
                      for (r, o), w in sorted(atom_diff.items(), key=lambda kv: -len(kv[1]))],
            "terms": [{"kind": k[0], "antechamber": k[1], "ours": k[2], "count": len(w),
                       "examples": w[:5]}
                      for k, w in sorted(term_diff.items(), key=lambda kv: -len(kv[1]))],
        }, indent=1) + "\n")


if __name__ == "__main__":
    main()
