"""Check opls_ff_gen.opls_typer against Foyer 1.2.0 (needs Foyer installed).

Two checks, both required before Foyer is dropped:

``patterns``
    For every SMARTS definition (Foyer's oplsaa.xml + the atom-type archive),
    our pattern graph must equal Foyer's SMARTSGraph: same nodes, same
    expression trees, same edges with the same bond orders, same rule-filter
    element.
``molecules``
    For every molecule of the given sets, the white- and blacklists of every
    atom after Foyer's iterative typing must be identical, and so must the
    final outcome of ``assignment._assign_atom_types`` (types, parameters and
    rules, or the error).

    python src/scripts/check_foyer_equivalence.py patterns
    python src/scripts/check_foyer_equivalence.py molecules --smi A.smi --csv B.csv \\
        --sdf C.sdf --out report.json [--workers 8]
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402

RDLogger.DisableLog("rdApp.*")
_FORCEFIELD_CACHE = None


def _foyer_definitions():
    global _FORCEFIELD_CACHE
    if _FORCEFIELD_CACHE is not None:
        return _FORCEFIELD_CACHE
    import foyer

    from opls_ff_gen.assignment import load_atomtype_archive

    forcefield = foyer.Forcefield(name="oplsaa")
    archive = load_atomtype_archive()
    forcefield.atomTypeDefinitions.update(archive.definitions)
    for preferred, shadowed in archive.precedence.items():
        forcefield.atomTypeOverrides.setdefault(preferred, set()).update(shadowed)
    _FORCEFIELD_CACHE = forcefield
    return forcefield


def _lark_tree(tree):
    import lark

    if isinstance(tree, lark.Tree):
        return (str(tree.data), tuple(_lark_tree(c) for c in tree.children))
    return str(tree)


def _our_tree(tree):
    from opls_ff_gen.opls_typer import Tree

    if isinstance(tree, Tree):
        return (tree.data, tuple(_our_tree(c) for c in tree.children))
    return str(tree)


def check_patterns() -> int:
    from foyer.smarts import SMARTS
    from foyer.smarts_graph import SMARTSGraph

    from opls_ff_gen.opls_typer import build_pattern

    forcefield = _foyer_definitions()
    parser = SMARTS()
    bad = 0
    for name, smarts in forcefield.atomTypeDefinitions.items():
        if not smarts:
            continue
        theirs = SMARTSGraph(smarts, parser=parser, name=name)
        ours = build_pattern(name, smarts)
        their_nodes = [_lark_tree(theirs.nodes[n]["atom"].children[0]) for n in sorted(theirs.nodes)]
        our_nodes = [_our_tree(node) for node in ours.nodes]
        their_edges = {(min(a, b), max(a, b)): float(d.get("bond_order", 0.0) or 0.0)
                       for a, b, d in theirs.edges(data=True)}
        atom0 = theirs.nodes[0]["atom"]
        symbols = list(atom0.find_data("atom_symbol"))
        their_filter = (str(symbols[0].children[0]) if len(symbols) == 1
                        and not list(atom0.find_data("not_expression")) else None)
        problems = []
        if their_nodes != our_nodes:
            problems.append("nodes")
        if their_edges != ours.edges:
            problems.append(f"edges {their_edges} vs {ours.edges}")
        if their_filter != ours.prefilter:
            problems.append(f"filter {their_filter} vs {ours.prefilter}")
        if problems:
            bad += 1
            print(f"{name} {smarts}: {'; '.join(problems)}")
    print(f"definitions checked: {sum(bool(s) for s in forcefield.atomTypeDefinitions.values())}, "
          f"mismatches: {bad}")
    return bad


# ------------------------------------------------------------------ molecules --
_FORCEFIELD = None


def _molecule(smiles: str | None, molblock: str | None):
    molecule = Chem.MolFromSmiles(smiles) if smiles else Chem.MolFromMolBlock(molblock, removeHs=False)
    if molecule is None:
        return None
    molecule = Chem.AddHs(molecule, addCoords=molblock is not None)
    if molecule.GetNumConformers() == 0:
        AllChem.Compute2DCoords(molecule)
    return molecule


def _rdkit_to_parmed(molecule):
    """ParmEd structure as the workflow built it for Foyer (bond orders as
    RDKit's GetBondTypeAsDouble, atoms and bonds in RDKit order)."""
    import parmed as pmd

    structure = pmd.Structure()
    atoms = []
    periodic_table = Chem.GetPeriodicTable()
    for index, atom in enumerate(molecule.GetAtoms(), start=1):
        parmed_atom = pmd.Atom(name=f"{atom.GetSymbol()}{index}", atomic_number=atom.GetAtomicNum(),
                               mass=periodic_table.GetAtomicWeight(atom.GetAtomicNum()))
        structure.add_atom(parmed_atom, resname="MOL", resnum=1)
        atoms.append(parmed_atom)
    for bond in molecule.GetBonds():
        structure.bonds.append(pmd.Bond(atoms[bond.GetBeginAtomIdx()], atoms[bond.GetEndAtomIdx()],
                                        order=bond.GetBondTypeAsDouble()))
    structure.coordinates = molecule.GetConformer().GetPositions().tolist()
    return structure


def _foyer_smarts_atom_types(molecule, archive):
    """Drop-in for assignment._smarts_atom_types that runs Foyer itself (the
    pre-removal implementation, tag foyer-reference; collisions stay errors)."""
    from foyer import atomtyper

    from opls_ff_gen.assignment import AssignmentError

    def resolve(graph, typemap):
        for atom_id, atom in typemap.items():
            candidates = atom["whitelist"] - atom["blacklist"]
            if len(candidates) == 1:
                atom["atomtype"] = next(iter(candidates))
            elif len(candidates) > 1:
                raise ValueError(f"Found multiple types for atom {atom_id}: {sorted(candidates)}.")

    try:
        forcefield = _foyer_definitions()
        with mock.patch.object(atomtyper, "_resolve_atomtypes", resolve):
            typemap = forcefield.run_atomtyping(_rdkit_to_parmed(molecule), use_residue_map=False)
    except Exception as exc:
        raise AssignmentError(f"OPLS-AA atom typing failed: {exc}") from exc
    return [typemap[index]["atomtype"] for index in range(molecule.GetNumAtoms())], {}


def _foyer_typemap(molecule):
    from foyer import atomtyper

    global _FORCEFIELD
    if _FORCEFIELD is None:
        _FORCEFIELD = _foyer_definitions()
    # raw white-/blacklists: skip Foyer's resolution step (compared separately)
    with mock.patch.object(atomtyper, "_resolve_atomtypes", lambda graph, typemap: None):
        typemap = _FORCEFIELD.run_atomtyping(_rdkit_to_parmed(molecule), use_residue_map=False)
    return [(sorted(typemap[i]["whitelist"]), sorted(typemap[i]["blacklist"]))
            for i in range(molecule.GetNumAtoms())]


def _our_typemap(molecule):
    from opls_ff_gen.assignment import load_atomtype_archive
    from opls_ff_gen.opls_typer import find_atomtypes
    from opls_ff_gen.assignment import typing_rules

    definitions, overrides = typing_rules(load_atomtype_archive())
    return [(sorted(a["whitelist"]), sorted(a["blacklist"]))
            for a in find_atomtypes(molecule, definitions, overrides)]


def _outcome(function, molecule):
    try:
        types, parameters, rules = function(molecule)
        return {"types": types, "rules": rules,
                "parameters": {t: parameters[t] for t in sorted(set(types))}}
    except Exception as exc:  # noqa: BLE001 - the error itself is compared
        return {"error": f"{type(exc).__name__}: {exc}"}


def _check_one(item):
    name, smiles, molblock = item
    molecule = _molecule(smiles, molblock)
    if molecule is None:
        return {"name": name, "skip": "rdkit"}
    result = {"name": name, "atoms": molecule.GetNumAtoms()}
    try:
        theirs = _foyer_typemap(molecule)
    except Exception as exc:  # noqa: BLE001
        theirs = f"error {type(exc).__name__}: {exc}"
    try:
        ours = _our_typemap(molecule)
    except Exception as exc:  # noqa: BLE001
        ours = f"error {type(exc).__name__}: {exc}"
    result["typemap_equal"] = theirs == ours
    if theirs != ours:
        result["typemap_diff"] = (
            [(i + 1, t, o) for i, (t, o) in enumerate(zip(theirs, ours)) if t != o][:5]
            if isinstance(theirs, list) and isinstance(ours, list) else [theirs, ours])
    from opls_ff_gen import assignment

    new = _outcome(assignment._assign_atom_types, molecule)
    with mock.patch.object(assignment, "_smarts_atom_types", _foyer_smarts_atom_types):
        old = _outcome(assignment._assign_atom_types, molecule)
    result["opls2020_equal"] = new == old
    if new != old:
        result["opls2020_diff"] = {"foyer": str(old)[:300], "native": str(new)[:300]}
    return result


def check_molecules(args) -> int:
    items = []
    for path in args.smi or []:
        for index, line in enumerate(path.read_text().splitlines()):
            fields = line.split()
            if fields and not fields[0].upper().startswith("SMILES"):
                items.append((f"{path.stem}:{index}", fields[0], None))
    for path in args.csv or []:
        with path.open() as handle:
            items += [(f"{path.stem}:{row['name']}", row["smiles"], None)
                      for row in csv.DictReader(handle)]
    for path in args.sdf or []:
        supplier = Chem.SDMolSupplier(str(path), removeHs=False)
        for index, molecule in enumerate(supplier):
            if molecule is not None:
                items.append((f"{path.stem}:{index}", None, Chem.MolToMolBlock(molecule)))
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for result in pool.map(_check_one, items, chunksize=4):
            results.append(result)
            if len(results) % 250 == 0 or len(results) == len(items):
                bad = sum(1 for r in results if "skip" not in r and not all(
                    r[k] for k in ("typemap_equal", "opls2020_equal")))
                print(f"progress {len(results)}/{len(items)} molecules, mismatches so far {bad}",
                      flush=True)
    checked = [r for r in results if "skip" not in r]
    keys = ("typemap_equal", "opls2020_equal")
    summary = {key: sum(r[key] for r in checked) for key in keys}
    bad = [r for r in checked if not all(r[key] for key in keys)]
    report = {"molecules": len(items), "checked": len(checked),
              "atoms": sum(r["atoms"] for r in checked), "equal": summary,
              "mismatches": bad}
    args.out.write_text(json.dumps(report, indent=1) + "\n")
    print(f"molecules {len(items)}, checked {len(checked)}, atoms {report['atoms']}; equal: "
          f"{summary}; mismatching molecules {len(bad)}")
    for r in bad[:20]:
        print("  ", r["name"], {k: v for k, v in r.items() if k.endswith("diff")})
    return len(bad)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("patterns", "molecules"))
    parser.add_argument("--smi", type=Path, nargs="*")
    parser.add_argument("--csv", type=Path, nargs="*")
    parser.add_argument("--sdf", type=Path, nargs="*")
    parser.add_argument("--out", type=Path, default=Path("foyer_equivalence.json"))
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    bad = check_patterns() if args.mode == "patterns" else check_molecules(args)
    raise SystemExit(1 if bad else 0)


if __name__ == "__main__":
    main()
