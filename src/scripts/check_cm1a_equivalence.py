"""Compare CM1A charges from MOPAC AM1 + opls_ff_gen.cm1a with AMSOL 7.1.

Needs the AMSOL 7.1 executable (``--amsol``); it carries its own copy of the
pre-replacement AMSOL code, so it does not depend on the package using AMSOL.
Each molecule gets one RDKit 3D geometry (the SDF's if it is 3D with explicit
hydrogens, otherwise ETKDG + MMFF94 as for SMILES); both engines run a single AM1 SCF on that geometry with the formal
net charge. Reported per molecule: the largest per-atom difference of the raw
CM1A charges, and whether either engine failed. Acceptance (decided before
the run): every atom within ``--tolerance`` (default 0.002 e).

    python src/scripts/check_cm1a_equivalence.py --smi A.smi --csv B.csv --sdf C.sdf \\
        --out report.json [--workers 9] [--tolerance 0.002]
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402

RDLogger.DisableLog("rdApp.*")
ELEMENTS = {1, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}  # CM1A-parameterized (AM1) elements


def _molecule(smiles: str | None, molblock: str | None):
    if molblock is not None:
        molecule = Chem.MolFromMolBlock(molblock, removeHs=False)
        if molecule is None:
            return None, "rdkit"
        flat = abs(molecule.GetConformer().GetPositions()[:, 2]).max() < 1e-3
        if flat or not any(a.GetAtomicNum() == 1 for a in molecule.GetAtoms()):
            # 2D drawings or heavy-atom-only records: build a 3D conformer as for SMILES
            smiles, molblock = Chem.MolToSmiles(Chem.RemoveHs(molecule)), None
    if smiles:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return None, "rdkit"
        molecule = Chem.AddHs(molecule)
        try:
            if AllChem.EmbedMolecule(molecule, randomSeed=20261018) != 0:
                return None, "embed"
        except RuntimeError:  # RDKit distance-geometry invariant failures on odd inputs
            return None, "embed"
        if AllChem.MMFFHasAllMoleculeParams(molecule):
            AllChem.MMFFOptimizeMolecule(molecule, maxIters=2000)
    else:
        molecule = Chem.MolFromMolBlock(molblock, removeHs=False)
        if molecule is None:
            return None, "rdkit"
        molecule = Chem.AddHs(molecule, addCoords=True)
    if len(Chem.GetMolFrags(molecule)) > 1:
        return None, "fragments"
    if any(a.GetAtomicNum() not in ELEMENTS or a.GetNumRadicalElectrons()
           for a in molecule.GetAtoms()):
        return None, "elements"
    return molecule, None


# ---- AMSOL 7.1 reference (the pre-replacement code, tag amsol-reference) ----
def _write_amsol_input(
    path: Path, name: str, molecule: Chem.Mol, net_charge: int
) -> None:
    keywords = ["AM1", "CM1", "1SCF", "CART"]
    if net_charge:
        keywords.append(f"CHARGE={net_charge}")
    conformer = molecule.GetConformer()
    lines = [
        " ".join(keywords),
        f"{name}: gas-phase CM1A single point on input coordinates",
        "",
    ]
    for index, atom in enumerate(molecule.GetAtoms()):
        position = conformer.GetAtomPosition(index)
        lines.append(
            f"  {atom.GetSymbol():<2s} {position.x:14.8f} 1"
            f" {position.y:14.8f} 1 {position.z:14.8f} 1"
        )
    lines.append("  0      0.00000000 0     0.00000000 0     0.00000000 0")
    path.write_text("\n".join(lines) + "\n")


def _run_amsol(
    executable: Path, input_path: Path, output_path: Path, timeout_seconds: int
) -> None:
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise RuntimeError(f"AMSOL executable is unavailable: {executable}")
    try:
        with input_path.open("rb") as input_file, output_path.open("wb") as output_file:
            result = subprocess.run(
                [str(executable)],
                stdin=input_file,
                stdout=output_file,
                stderr=subprocess.PIPE,
                cwd=output_path.parent,
                check=False,
                timeout=timeout_seconds,
            )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"AMSOL exceeded the {timeout_seconds}-second timeout"
        ) from exc
    if result.returncode != 0:
        stderr = result.stderr.decode(errors="replace")
        raise RuntimeError(
            f"AMSOL failed with exit code {result.returncode}: {stderr.strip()}"
        )


def _parse_amsol_output(path: Path) -> tuple[list[tuple[int, str, float]], float]:
    text = path.read_text(errors="replace")
    if "SCF FIELD WAS ACHIEVED" not in text:
        raise RuntimeError(f"{path}: AMSOL did not report SCF convergence")
    heat_match = re.search(
        r"Final heat of formation\s*=\s*([-+]?\d+\.\d+)\s+kcal", text
    )
    if heat_match is None:
        raise RuntimeError(f"{path}: AM1 heat of formation was not found")

    section_start = text.find("calculated with CM1")
    section_end = text.find("For the charges calculated by CM1:", section_start)
    if section_start < 0 or section_end < 0:
        raise RuntimeError(f"{path}: CM1A charge table was not found")
    pattern = re.compile(
        r"^\s*(\d+)\s+([A-Za-z]{1,2})\s+([-+]?\d+\.\d+)", re.MULTILINE
    )
    charges = [
        (int(match.group(1)), match.group(2), float(match.group(3)))
        for match in pattern.finditer(text[section_start:section_end])
    ]
    if not charges:
        raise RuntimeError(f"{path}: CM1A charge table was empty")
    return charges, float(heat_match.group(1))


def _amsol_charges(molecule, workdir: Path, amsol: Path):
    net_charge = Chem.GetFormalCharge(molecule)
    _write_amsol_input(workdir / "m.amsol.in", "m", molecule, net_charge)
    _run_amsol(amsol, workdir / "m.amsol.in", workdir / "m.amsol.out", 600)
    records, heat = _parse_amsol_output(workdir / "m.amsol.out")
    if [(i, e) for i, e, _ in records] != [(a.GetIdx() + 1, a.GetSymbol())
                                          for a in molecule.GetAtoms()]:
        raise RuntimeError("AMSOL atom order differs from the input")
    return np.array([q for _, _, q in records]), heat


def _charges(engine: str, molecule, workdir: Path, amsol: Path):
    workdir.mkdir(parents=True, exist_ok=True)
    if engine == "amsol":
        return _amsol_charges(molecule, workdir, amsol)
    from opls_ff_gen.cm1a import run_cm1a

    result = run_cm1a(molecule, workdir, "m", Chem.GetFormalCharge(molecule), 600)
    return result.charges, result.heat_of_formation_kcal_mol


def _check_one(item):
    name, smiles, molblock, amsol = item
    molecule, skip = _molecule(smiles, molblock)
    if molecule is None:
        return {"name": name, "skip": skip}
    result = {"name": name, "atoms": molecule.GetNumAtoms(),
              "net_charge": Chem.GetFormalCharge(molecule)}
    with tempfile.TemporaryDirectory() as scratch:
        scratch = Path(scratch)
        outcome = {}
        for engine in ("amsol", "mopac"):
            try:
                outcome[engine] = _charges(engine, molecule, scratch / engine, Path(amsol))
            except Exception as exc:  # noqa: BLE001 - failures are compared
                outcome[engine] = f"{type(exc).__name__}: {str(exc)[:200]}"
    if isinstance(outcome["amsol"], str) or isinstance(outcome["mopac"], str):
        result["failed"] = {e: o for e, o in outcome.items() if isinstance(o, str)}
        return result
    (ref, ref_heat), (new, new_heat) = outcome["amsol"], outcome["mopac"]
    diff = np.abs(new - ref)
    worst = int(diff.argmax())
    result.update(max_abs_dq=float(diff.max()), mean_abs_dq=float(diff.mean()),
                  worst_atom=worst + 1,
                  worst_element=molecule.GetAtomWithIdx(worst).GetSymbol(),
                  heat_amsol=ref_heat, heat_mopac=new_heat, smiles=Chem.MolToSmiles(
                      Chem.RemoveHs(molecule)))
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smi", type=Path, nargs="*", default=[])
    parser.add_argument("--csv", type=Path, nargs="*", default=[])
    parser.add_argument("--sdf", type=Path, nargs="*", default=[])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--amsol", type=Path, default=Path("/home/jiang/software/amsol7.1/amsol7.1"))
    parser.add_argument("--workers", type=int, default=9)
    parser.add_argument("--tolerance", type=float, default=0.002)
    args = parser.parse_args()
    items = []
    for path in args.smi:
        for index, line in enumerate(path.read_text().splitlines()):
            fields = line.split()
            if fields and not fields[0].upper().startswith("SMILES"):
                items.append((f"{path.stem}:{index}", fields[0], None, str(args.amsol)))
    for path in args.csv:
        with path.open() as handle:
            items += [(f"{path.stem}:{row['name']}", row["smiles"], None, str(args.amsol))
                      for row in csv.DictReader(handle)]
    for path in args.sdf:
        for index, molecule in enumerate(Chem.SDMolSupplier(str(path), removeHs=False)):
            if molecule is not None:
                items.append((f"{path.stem}:{index}", None, Chem.MolToMolBlock(molecule),
                              str(args.amsol)))
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for result in pool.map(_check_one, items, chunksize=8):
            results.append(result)
            if len(results) % 500 == 0 or len(results) == len(items):
                compared = [r for r in results if "max_abs_dq" in r]
                worst = max((r["max_abs_dq"] for r in compared), default=0.0)
                print(f"progress {len(results)}/{len(items)}: compared {len(compared)}, "
                      f"worst |dq| so far {worst:.4f} e", flush=True)
    compared = [r for r in results if "max_abs_dq" in r]
    failed = [r for r in results if "failed" in r]
    skipped = {}
    for r in results:
        if "skip" in r:
            skipped[r["skip"]] = skipped.get(r["skip"], 0) + 1
    over = [r for r in compared if r["max_abs_dq"] > args.tolerance]
    maxima = np.array([r["max_abs_dq"] for r in compared])
    summary = {
        "molecules": len(items), "compared": len(compared),
        "atoms": int(sum(r["atoms"] for r in compared)), "skipped": skipped,
        "tolerance_e": args.tolerance, "over_tolerance": len(over),
        "max_abs_dq": float(maxima.max()) if len(maxima) else None,
        "molecule_max_abs_dq_percentiles": {p: float(np.percentile(maxima, p))
                                            for p in (50, 90, 99, 99.9)} if len(maxima) else {},
        "failed_both": sum(1 for r in failed if len(r["failed"]) == 2),
        "failed_amsol_only": sum(1 for r in failed if list(r["failed"]) == ["amsol"]),
        "failed_mopac_only": sum(1 for r in failed if list(r["failed"]) == ["mopac"]),
    }
    args.out.write_text(json.dumps({"summary": summary,
                                    "over_tolerance": sorted(over, key=lambda r: -r["max_abs_dq"]),
                                    "failed": failed,
                                    "worst_20": sorted(compared, key=lambda r: -r["max_abs_dq"])[:20]},
                                   indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    # AMSOL-only failures (MOPAC succeeds) are reported but are not a regression
    raise SystemExit(1 if over or summary["failed_mopac_only"] else 0)


if __name__ == "__main__":
    main()
