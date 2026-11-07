"""Antechamber GAFF2 atom types for a molecule set (reference for the lookup).

Needs AmberTools (antechamber); used only to build and check the rule table in
``src/opls_ff_gen/data/gaff2_lj.json``. Molecules are read as SMILES, given
explicit hydrogens and 2D coordinates (typing does not depend on geometry),
and typed with ``antechamber -at gaff2``; the RDKit atom order is kept.
Molecules with elements outside H C N O F P S Cl Br I, radicals, or that
antechamber rejects are skipped (and counted).

    python src/scripts/gaff2_reference.py --out REF.json [--workers 8] \\
        [--smi FILE.smi ...] [--csv FILE.csv ...]
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402

RDLogger.DisableLog("rdApp.*")
ELEMENTS = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53}
NCI = Path(sys.prefix) / "share" / "RDKit" / "Data" / "NCI" / "first_5K.smi"


AMBERHOME = Path(os.environ.get("AMBERHOME", "/home/jiang/miniforge3/envs/chem"))


def antechamber_types(molecule: Chem.Mol) -> list[str]:
    """GAFF2 types from ``antechamber -at gaff2`` (input atom order kept)."""
    with tempfile.TemporaryDirectory() as scratch:
        scratch = Path(scratch)
        writer = Chem.SDWriter(str(scratch / "in.sdf"))
        writer.write(molecule)
        writer.close()
        result = subprocess.run(
            [str(AMBERHOME / "bin" / "antechamber"), "-i", "in.sdf", "-fi", "sdf", "-o",
             "out.mol2", "-fo", "mol2", "-at", "gaff2", "-dr", "no", "-pf", "y"],
            cwd=scratch, env=dict(os.environ, AMBERHOME=str(AMBERHOME)),
            capture_output=True, text=True, check=False,
        )
        mol2 = scratch / "out.mol2"
        if result.returncode != 0 or not mol2.is_file():
            raise RuntimeError(result.stdout[-300:])
        block = mol2.read_text().split("@<TRIPOS>ATOM", 1)[1].split("@<TRIPOS>", 1)[0]
    types = [line.split()[5] for line in block.strip().splitlines()]
    if len(types) != molecule.GetNumAtoms():
        raise RuntimeError("atom count mismatch")
    return types


def _type(item: tuple[str, str]) -> dict:
    name, smiles = item
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return {"name": name, "smiles": smiles, "skip": "rdkit"}
    if any(a.GetAtomicNum() not in ELEMENTS or a.GetNumRadicalElectrons()
           for a in molecule.GetAtoms()):
        return {"name": name, "smiles": smiles, "skip": "elements"}
    if len(Chem.GetMolFrags(molecule)) > 1:
        return {"name": name, "smiles": smiles, "skip": "fragments"}
    molecule = Chem.AddHs(molecule)
    AllChem.Compute2DCoords(molecule)
    try:
        types = antechamber_types(molecule)
    except Exception as exc:  # noqa: BLE001 - antechamber failure is a skip
        return {"name": name, "smiles": smiles, "skip": f"antechamber: {str(exc)[:80]}"}
    # input SMILES: MolFromSmiles + AddHs rebuilds the same atom order as the types
    return {"name": name, "smiles": smiles, "types": types}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--smi", type=Path, nargs="*", default=[NCI])
    parser.add_argument("--csv", type=Path, nargs="*", default=[])
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    items = []
    for path in args.smi:
        for index, line in enumerate(path.read_text().splitlines()):
            fields = line.split()
            if fields and not fields[0].upper().startswith("SMILES"):
                items.append((f"{path.stem}:{index}", fields[0]))
    for path in args.csv:
        with path.open() as handle:
            items += [(row["name"], row["smiles"]) for row in csv.DictReader(handle)]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(_type, items, chunksize=20))
    typed = [r for r in results if "types" in r]
    skipped = {}
    for r in results:
        if "skip" in r:
            key = r["skip"].split(":")[0]
            skipped[key] = skipped.get(key, 0) + 1
    args.out.write_text(json.dumps({"source": "antechamber -at gaff2 (AmberTools 26)",
                                    "molecules": typed}, indent=0) + "\n")
    print(f"typed {len(typed)} molecules, {sum(len(r['types']) for r in typed)} atoms; "
          f"skipped {skipped}")


if __name__ == "__main__":
    main()
