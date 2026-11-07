"""End-to-end tests for reusable GROMACS topology generation."""

from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

from opls_ff_gen import generate_gromacs_topology, load_sdf, molecule_from_smiles
from opls_ff_gen.gromacs import gmx_executable


ROOT = Path(__file__).resolve().parents[1]
MOLECULES = ROOT / "molecules_to_test"
GMX = gmx_executable() or Path("gmx-not-found")


def _section_counts(path: Path) -> dict[str, int]:
    counts: dict[str, int] = {}
    section = ""
    for raw_line in path.read_text().splitlines():
        line = raw_line.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("["):
            section = line.strip("[] ").lower()
            counts.setdefault(section, 0)
        elif section:
            counts[section] += 1
    return counts


class CompleteWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp_directory = tempfile.TemporaryDirectory()
        root = Path(cls.temp_directory.name)
        cls.results = {
            "sdf": generate_gromacs_topology(
                load_sdf(MOLECULES / "ethanol.sdf"),
                root / "from_sdf",
                name="ethanol_sdf",
            ),
            "smiles": generate_gromacs_topology(
                molecule_from_smiles("COc1ccccc1", name="anisole_smiles"),
                root / "from_smiles",
                name="anisole_smiles",
            ),
        }

    @classmethod
    def tearDownClass(cls):
        cls.temp_directory.cleanup()

    def test_both_input_routes_write_complete_artifact_sets(self):
        for result in self.results.values():
            self.assertTrue(result.itp_path.is_file())
            self.assertTrue(result.gro_path.is_file())
            self.assertTrue(result.top_path.is_file())
            self.assertTrue(result.topology_json_path.is_file())
            self.assertTrue(result.assignment.sdf_path.is_file())
            self.assertTrue(result.assignment.charge_output_path.is_file())

    def test_itp_contains_complete_topologies(self):
        expected = {
            "sdf": {
                "atomtypes": 9,
                "moleculetype": 1,
                "atoms": 9,
                "bonds": 8,
                "pairs": 12,
                "angles": 13,
                "dihedrals": 12,
            },
            "smiles": {
                "atomtypes": 16,
                "moleculetype": 1,
                "atoms": 16,
                "bonds": 16,
                "pairs": 26,
                "angles": 25,
                "dihedrals": 35,
            },
        }
        for route, result in self.results.items():
            counts = _section_counts(result.itp_path)
            self.assertEqual(counts, expected[route])
            self.assertNotIn("defaults", counts)
            # self-contained: the TOP carries the OPLS [ defaults ], no oplsaa.ff include
            top = result.top_path.read_text()
            self.assertIn("[ defaults ]", top)
            self.assertNotIn("oplsaa.ff", top)

    def test_manifest_is_auditable(self):
        expected = {
            "sdf": (8, 12, 0),
            "smiles": (16, 29, 6),
        }
        for route, result in self.results.items():
            payload = json.loads(result.topology_json_path.read_text())
            self.assertEqual(payload["charge_model"], "1.14*CM1A-LBCC")
            self.assertEqual(payload["defaults"]["nrexcl"], 3)
            self.assertEqual(payload["defaults"]["fudge_lj"], 0.5)
            self.assertEqual(payload["defaults"]["fudge_qq"], 0.5)
            bonds, propers, impropers = expected[route]
            self.assertEqual(len(payload["bonds"]), bonds)
            self.assertEqual(len(payload["proper_torsions"]), propers)
            self.assertEqual(len(payload["improper_torsions"]), impropers)
            self.assertTrue(all(item["source_line"] > 0 for item in payload["bonds"]))

    @unittest.skipUnless(GMX.is_file(), "local GROMACS executable is unavailable")
    def test_generated_topologies_pass_grompp(self):
        mdp = Path(self.temp_directory.name) / "check.mdp"
        mdp.write_text(
            "\n".join(
                [
                    "integrator = md",
                    "nsteps = 0",
                    "dt = 0.001",
                    "continuation = yes",
                    "cutoff-scheme = Verlet",
                    "nstlist = 10",
                    "rlist = 0.8",
                    "coulombtype = Cut-off",
                    "rcoulomb = 0.8",
                    "vdwtype = Cut-off",
                    "rvdw = 0.8",
                    "pbc = xyz",
                    "constraints = none",
                    "",
                ]
            )
        )
        for route, result in self.results.items():
            completed = subprocess.run(
                [
                    str(GMX),
                    "grompp",
                    "-f",
                    str(mdp),
                    "-c",
                    str(result.gro_path),
                    "-p",
                    str(result.top_path),
                    "-o",
                    str(result.top_path.with_suffix(".tpr")),
                    "-po",
                    str(result.top_path.with_suffix(".mdp")),
                ],
                cwd=result.top_path.parent,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, f"{route}:\n{completed.stdout}")


if __name__ == "__main__":
    unittest.main()
