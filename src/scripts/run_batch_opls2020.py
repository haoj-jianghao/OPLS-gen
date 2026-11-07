"""Run the OPLS/2020 workflow (on-the-fly QM fitting) on CSVs of molecules.

Each molecule is generated once with ``--qm METHOD`` (default none; OPLS/2020 tables,
GAFF2 LJ for atoms OPLS cannot type, missing bonded terms fitted to QM with --qm), then
checked with ``gmx grompp`` (0 warnings allowed) and a zero-step ``gmx mdrun``
(finite potential energy). Molecules run one at a time; the QM fit inside a
molecule uses ``--qm-threads`` workers. GROMACS is found as in
``opls_ff_gen.gromacs`` (``OPLS_FF_GEN_GMX``, else ``gmx`` on PATH).

    python src/scripts/run_batch_opls2020.py CSV [CSV ...] --out DIR [--qm xtb|dft]
        [--mm-engine openmm|gromacs] [--torsion-substitution parmchk2|none]
        [--qm-bonded-table use|none] [--optimize none|mmff|xtb]
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from opls_ff_gen.gromacs import require_gmx  # noqa: E402

MDP = """integrator = md
nsteps = 0
dt = 0.001
continuation = yes
cutoff-scheme = Verlet
nstlist = 10
rlist = 0.8
coulombtype = Cut-off
rcoulomb = 0.8
vdwtype = Cut-off
rvdw = 0.8
pbc = xyz
constraints = none
nstenergy = 1
"""


def slug(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_").lower()


def run(command: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        command, cwd=cwd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )


def gromacs_check(directory: Path, stem: str) -> dict:
    """``gmx grompp`` (0 warnings allowed) and a zero-step ``gmx mdrun``."""
    gmx = require_gmx()
    (directory / "check.mdp").write_text(MDP)
    grompp = run(
        [str(gmx), "grompp", "-f", "check.mdp", "-c", f"{stem}.gro", "-p", f"{stem}.top",
         "-o", "check.tpr", "-po", "check.out.mdp", "-maxwarn", "0"],
        cwd=directory,
    )
    (directory / "grompp.log").write_text(grompp.stdout)
    result = {
        "grompp_ok": grompp.returncode == 0,
        "grompp_warnings": len(re.findall(r"^WARNING \d+", grompp.stdout, re.M)),
        "grompp_notes": len(re.findall(r"^NOTE \d+", grompp.stdout, re.M)),
    }
    if not result["grompp_ok"]:
        errors = re.findall(r"ERROR \d+.*?\n\n(.*?)\n\n", grompp.stdout, re.S)
        result["grompp_error"] = (errors[0] if errors else grompp.stdout[-800:]).strip()
        return result
    mdrun = run(
        [str(gmx), "mdrun", "-s", "check.tpr", "-deffnm", "check", "-nt", "1",
         "-nb", "cpu", "-pme", "cpu"],
        cwd=directory,
    )
    result["mdrun_ok"] = mdrun.returncode == 0
    if result["mdrun_ok"]:
        energy = subprocess.run(
            [str(gmx), "energy", "-f", "check.edr", "-o", "check_energy.xvg"],
            cwd=directory, input="Potential\n\n", text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )
        match = re.search(r"^Potential\s+(\S+)", energy.stdout, re.M)
        result["potential_kj_mol"] = float(match.group(1)) if match else None
    else:
        result["mdrun_error"] = mdrun.stdout[-800:]
    return result




def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", nargs="+", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--qm", default="none", choices=("xtb", "dft", "none"))
    parser.add_argument("--qm-threads", type=int, default=8)
    parser.add_argument("--optimize", choices=("none", "mmff", "xtb"),
                        help="geometry optimization of the SMILES conformer (default xtb)")
    parser.add_argument("--mm-engine", choices=("openmm", "gromacs"))
    parser.add_argument("--torsion-substitution", choices=("parmchk2", "none"))
    parser.add_argument("--qm-bonded-table", choices=("use", "none"))
    parser.add_argument("--dft-functional")
    parser.add_argument("--dft-basis")
    parser.add_argument("--dft-dispersion")
    parser.add_argument("--only", nargs="*", help="molecule names to run (default all)")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    summary_path = args.out / "batch_summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    molecules = []
    for path in args.csv:
        with path.open() as handle:
            molecules += [(path.parent.name, row["name"], row["smiles"])
                          for row in csv.DictReader(handle)]
    for batch, name, smiles in molecules:
        if args.only and name not in args.only:
            continue
        stem = slug(name)
        directory = args.out / stem
        directory.mkdir(exist_ok=True)
        clock = time.perf_counter()
        generate = subprocess.run(
            [sys.executable, "-m", "opls_ff_gen", "--smiles", smiles, "--name", stem,
             "--output-dir", str(directory), "--qm", args.qm,
             "--qm-threads", str(args.qm_threads),
             *(["--optimize", args.optimize] if args.optimize else []),
             *(["--mm-engine", args.mm_engine] if args.mm_engine else []),
             *(["--torsion-substitution", args.torsion_substitution]
               if args.torsion_substitution else []),
             *(["--qm-bonded-table", args.qm_bonded_table]
               if args.qm_bonded_table else []),
             *[arg for option, value in (("--dft-functional", args.dft_functional),
                                         ("--dft-basis", args.dft_basis),
                                         ("--dft-dispersion", args.dft_dispersion))
               if value for arg in (option, value)]],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )
        seconds = time.perf_counter() - clock
        (directory / "generate.log").write_text(generate.stdout)
        entry = {"batch": batch, "name": name, "smiles": smiles, "generate_ok":
                 generate.returncode == 0, "wall_s": round(seconds, 1)}
        if generate.returncode == 0:
            topology = json.loads(next(directory.glob("*.topology.json")).read_text())
            entry["rule_counts"] = topology["bonded_rule_counts"]
            entry["warnings"] = topology["bonded_warnings"]
            entry["amber_atoms"] = {k: v for k, v in topology["borrowed_atoms"].items()
                                    if v.startswith("amber")}
            entry["borrowed_atoms"] = {k: v for k, v in topology["borrowed_atoms"].items()
                                       if v.startswith("borrowed")}
            report = directory / "qm_fit" / "report.json"
            if report.exists():
                qm = json.loads(report.read_text())
                entry["qm"] = {
                    "timings_s": qm["timings_s"],
                    "bonds": len(qm["bonds"]), "angles": len(qm["angles"]),
                    "scans": [{"bond": s["bond"], "ring": s["ring"],
                               "rmse_zero": round(s["rmse_kj_mol"]["missing_zero"], 2),
                               "rmse_fit": round(s["rmse_kj_mol"]["fitted"], 2)}
                              for s in qm["scans"]],
                    "torsion_keys": len(qm["torsion_keys"]),
                    "zero_by_convention": len(qm.get("zero_by_convention", [])),
                    "zero_linear": len(qm.get("zero_linear", [])),
                }
            entry.update(gromacs_check(directory, Path(topology["files"]["top"]).stem))
        else:
            entry["error"] = generate.stdout.strip().splitlines()[-1][:400] if generate.stdout else ""
        summary[name] = entry
        summary_path.write_text(json.dumps(summary, indent=2) + "\n")
        status = ("OK" if entry.get("mdrun_ok") else "FAIL")
        print(f"{name}: {status} {entry['wall_s']} s {entry.get('rule_counts', entry.get('error'))}",
              flush=True)


if __name__ == "__main__":
    main()
