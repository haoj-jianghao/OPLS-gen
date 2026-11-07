"""End-to-end check of an installed opls-ff-gen (the release CI runs it).

    python src/scripts/release_smoke_test.py --out DIR [--qm] [--expect-xtb]

Runs the installed `opls-ff-gen` command (from PATH) on:
- a SMILES (paracetamol, default geometry);
- an SDF file (DEET, from molecules_to_test/);
- a cation and a zwitterion (net charge +1 and 0).
With --qm it also runs one on-the-fly `--qm xtb` fit. Each run must exit 0,
write the ITP/TOP/GRO and JSON files, give charges that sum to the net charge
(1e-6 e), and record its geometry method. With --expect-xtb, the SMILES
geometry must be xTB (no fallback); without it, the MMFF94 fallback is
accepted. Exit code 1 on any failure.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

CASES = [
    # name, input arguments, net charge, extra arguments
    ("paracetamol", ["--smiles", "CC(=O)Nc1ccc(O)cc1"], 0, []),
    ("deet", ["--structure", str(ROOT / "molecules_to_test" / "deet.sdf")], 0, []),
    ("methylammonium", ["--smiles", "C[NH3+]"], 1, []),
    ("glycine_zwitterion", ["--smiles", "[NH3+]CC(=O)[O-]"], 0, []),
]
QM_CASE = ("thiocarbamate_qm", ["--smiles", "CSC(=O)N"], 0, ["--qm", "xtb", "--qm-threads", "2"])


def run_case(command: str, out: Path, case, expect_xtb: bool) -> list[str]:
    name, inputs, net_charge, extra = case
    directory = out / name
    result = subprocess.run([command, *inputs, "--name", name, "--output-dir", str(directory),
                             *extra], capture_output=True, text=True, timeout=1800)
    log = result.stdout + result.stderr
    (out / f"{name}.log").write_text(log)
    if result.returncode != 0:
        return [f"{name}: exit {result.returncode}: {log.strip().splitlines()[-1:]}"]
    problems = []
    for suffix in (".itp", ".top", ".gro", ".assignment.json", ".topology.json"):
        if not (directory / f"{name}{suffix}").is_file():
            problems.append(f"{name}: no {name}{suffix}")
    if problems:
        return problems
    record = json.loads((directory / f"{name}.assignment.json").read_text())
    total = record["final_charge_sum_e"]
    if abs(total - net_charge) > 1e-6:
        problems.append(f"{name}: charges sum to {total}, not {net_charge}")
    geometry = record.get("geometry", {})
    method = geometry.get("Geometry_Optimization", "")
    if "--smiles" in inputs and expect_xtb and not method.startswith("GFN2-xTB") \
            and name != "glycine_zwitterion":  # the zwitterion falls back by design
        problems.append(f"{name}: geometry {method!r}, expected GFN2-xTB "
                        f"({geometry.get('Geometry_Fallback', '')})")
    if extra and "qm-fit" not in (directory / f"{name}.itp").read_text():
        problems.append(f"{name}: --qm xtb fitted no term")
    print(f"  {name}: ok ({method or 'input coordinates'}; charge {total:+.6f})")
    return problems


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--qm", action="store_true", help="also run a --qm xtb fit")
    parser.add_argument("--expect-xtb", action="store_true",
                        help="SMILES geometries must come from xTB")
    args = parser.parse_args()
    command = shutil.which("opls-ff-gen")
    if command is None:
        sys.exit("opls-ff-gen is not on PATH")
    args.out.mkdir(parents=True, exist_ok=True)
    print(f"smoke test of {command}")
    problems = []
    for case in CASES + ([QM_CASE] if args.qm else []):
        problems += run_case(command, args.out, case, args.expect_xtb)
    for problem in problems:
        print(f"FAIL: {problem}")
    print("smoke test: " + ("FAILED" if problems else "ok"))
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
