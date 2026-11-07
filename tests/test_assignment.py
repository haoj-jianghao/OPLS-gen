"""Regression tests for atom typing (SMARTS stage) and CM1A-LBCC charges."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from rdkit import Chem

from opls_ff_gen import assign_molecule, load_sdf
from opls_ff_gen.assignment import smarts_atom_types
from opls_ff_gen.lbcc import LBCC_RULES, apply_lbcc, find_lbcc_applications


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MOLECULE_DIR = PROJECT_ROOT / "molecules_to_test"


def read_itp_charges(path: Path) -> list[float]:
    charges = []
    in_atoms = False
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if line.startswith("[") and line.endswith("]"):
            in_atoms = line.strip("[] ").lower() == "atoms"
            continue
        if in_atoms and line and not line.startswith(";"):
            fields = line.split(";", 1)[0].split()
            if len(fields) >= 8:
                charges.append(float(fields[6]))
    return charges


class AssignmentRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp_directory = tempfile.TemporaryDirectory()
        cls.results = {}
        for name in ("ethanol", "phenoxyethanol"):
            molecule = load_sdf(MOLECULE_DIR / f"{name}.sdf")
            cls.results[name] = assign_molecule(
                molecule,
                Path(cls.temp_directory.name) / name,
                name=name,
                charge_model="cm1a-lbcc",
            )

    @classmethod
    def tearDownClass(cls):
        cls.temp_directory.cleanup()

    def test_ligpargen_charge_regressions(self):
        for name, result in self.results.items():
            reference = read_itp_charges(MOLECULE_DIR / f"{name}.itp")
            candidate = [atom.final_charge_e for atom in result.atoms]
            self.assertEqual(len(reference), len(candidate))
            self.assertLessEqual(
                max(abs(ref - value) for ref, value in zip(reference, candidate)),
                3e-4,
            )
            self.assertAlmostEqual(sum(candidate), 0.0, places=10)

    def test_smarts_type_regressions(self):
        expected = {
            "ethanol": [
                "opls_135",
                "opls_157",
                "opls_154",
                "opls_140",
                "opls_140",
                "opls_140",
                "opls_140",
                "opls_140",
                "opls_155",
            ],
            "phenoxyethanol": [
                *("opls_145" for _ in range(6)),
                "opls_179",
                "opls_182",
                "opls_157",
                "opls_154",
                *("opls_146" for _ in range(5)),
                "opls_185",
                "opls_185",
                "opls_140",
                "opls_140",
                "opls_155",
            ],
        }
        for name in self.results:
            types, _ = smarts_atom_types(load_sdf(MOLECULE_DIR / f"{name}.sdf"))
            self.assertEqual(types, expected[name])

    def test_expected_lbcc_applications(self):
        ethanol_rules = [app.rule for app in self.results["ethanol"].lbcc_applications]
        phenoxy_rules = [
            app.rule for app in self.results["phenoxyethanol"].lbcc_applications
        ]
        self.assertEqual(ethanol_rules, ["CT-OH"])
        self.assertEqual(phenoxy_rules.count("CT-OH"), 1)
        self.assertEqual(phenoxy_rules.count("HA-CA"), 5)

    def test_machine_readable_outputs_exist(self):
        for result in self.results.values():
            self.assertTrue(result.json_path.is_file())
            self.assertTrue(result.csv_path.is_file())
            self.assertTrue(result.sdf_path.is_file())
            self.assertTrue(result.charge_output_path.is_file())

    def test_known_rule_overlap_precedence_is_deterministic(self):
        expected = {
            "icaridin": (4, "opls_467"),
            "transfluthrin": (7, "opls_490"),
        }
        for name, (zero_index, atom_type) in expected.items():
            molecule = load_sdf(MOLECULE_DIR / f"{name}.sdf")
            assigned, _ = smarts_atom_types(molecule)
            self.assertEqual(assigned[zero_index], atom_type)

    def test_alkene_r2c_carbon_types_in_and_out_of_rings(self):
        from rdkit.Chem import AllChem

        expected = {
            "CC1=CCC(CC1)C(=C)C": {1: "opls_141", 7: "opls_141"},  # limonene
            "CC(C)=CCCC(C)(O)C=C": {1: "opls_141"},  # linalool
            "Cc1ccccc1": {1: "opls_145"},  # toluene ipso C stays aromatic
        }
        for smiles, atoms in expected.items():
            molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
            AllChem.EmbedMolecule(molecule, randomSeed=1)
            assigned, _ = smarts_atom_types(molecule)
            for zero_index, atom_type in atoms.items():
                self.assertEqual(assigned[zero_index], atom_type, smiles)

    def test_typing_fixes_from_batch_2026_10_17(self):
        from rdkit.Chem import AllChem

        expected = {
            "CCN(CC)CC(=O)Nc1c(C)cccc1C": {3: "opls_908", 5: "opls_908"},  # lidocaine
            "CCN(CC)CC": {1: "opls_908"},  # triethylamine keeps its type
            # Coumarin lactone is typed as an ester, not a phenyl ether.
            "O=c1ccc2ccccc2o1": {0: "opls_466", 1: "opls_465", 10: "opls_467"},
            "COc1ccccc1Oc1ccccc1": {8: "opls_473"},  # diaryl ether O unchanged
            # OPLS-AA sulfonamide, imidazole and isoxazole types (no Foyer SMARTS).
            "NS(=O)(=O)c1ccccc1": {
                0: "opls_478", 1: "opls_474", 2: "opls_475", 4: "opls_488",
                10: "opls_479",
            },
            "c1cnc[nH]1": {
                0: "opls_561", 1: "opls_560", 2: "opls_559", 3: "opls_558",
                4: "opls_557", 5: "opls_565", 6: "opls_564", 7: "opls_563",
                8: "opls_562",
            },
            "Cc1ccno1": {1: "opls_583", 2: "opls_582", 3: "opls_581", 4: "opls_580",
                         5: "opls_579"},
        }
        for smiles, atoms in expected.items():
            molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
            AllChem.EmbedMolecule(molecule, randomSeed=1)
            assigned, _ = smarts_atom_types(molecule)
            for zero_index, atom_type in atoms.items():
                self.assertEqual(assigned[zero_index], atom_type, smiles)


class BorrowedTypeTests(unittest.TestCase):
    """The atom-type archive borrows existing OPLS-AA types for gaps."""

    EXPECTED = {
        # SMILES: {zero-based atom index: (SMARTS-stage type, rule)}
        "CCCC[N+](C)(C)[O-]": {
            4: ("opls_288", "borrowed:amine-oxide-N"),
            7: ("opls_272", "borrowed:amine-oxide-O"),
            5: ("opls_291", "borrowed:ammonium-C"),
        },
        "NNC(=O)c1ccncc1": {
            0: ("opls_900", "borrowed:hydrazine-NH2"),
            1: ("opls_238", "borrowed:hydrazide-NH"),
        },
        "Cn1cncn1": {
            0: ("opls_135", "borrowed:ring-N-CH3"),
            1: ("opls_557", "borrowed:triazole-N1"),
            2: ("opls_558", "borrowed:triazole-C"),
            3: ("opls_559", "borrowed:triazole-N4"),
            5: ("opls_559", "borrowed:triazole-N2"),
        },
        "Cn1c(=O)c2c(ncn2C)n(C)c1=O": {
            4: ("opls_624", "borrowed:purine-C5"),
            5: ("opls_623", "borrowed:purine-C4"),
            6: ("opls_626", "borrowed:purine-N"),
            7: ("opls_627", "borrowed:purine-C8"),
            8: ("opls_628", "borrowed:purine-N-subst"),
            9: ("opls_135", "borrowed:ring-N-CH3"),
        },
    }

    def test_borrowed_types_and_provenance(self):
        from rdkit.Chem import AllChem

        from opls_ff_gen.assignment import AssignmentError, _assign_atom_types

        for smiles, atoms in self.EXPECTED.items():
            molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
            AllChem.EmbedMolecule(molecule, randomSeed=1)
            with self.assertRaisesRegex(AssignmentError, "no type for atoms"):
                _assign_atom_types(molecule, None, amber_fallback=False)  # rules only
            types, rules = smarts_atom_types(molecule)
            for index, (opls_type, rule) in atoms.items():
                self.assertEqual((types[index], rules[index]), (opls_type, rule), smiles)
            _assign_atom_types(molecule)  # every borrowed type is in OPLS/2020

    def test_borrowing_leaves_rule_typed_molecules_unchanged(self):
        for name in ("deet", "icaridin", "transfluthrin"):
            molecule = load_sdf(MOLECULE_DIR / f"{name}.sdf")
            _, rules = smarts_atom_types(molecule)
            self.assertEqual(set(rules), {"smarts"})



    def test_all_19_published_rules_are_reachable_and_charge_conserving(self):
        representative_smiles = (
            "c1ccccc1",
            "n1ccccc1",
            "O=[N+]([O-])c1ccccc1",
            "Nc1ccccc1",
            "Fc1ccccc1",
            "Brc1ccccc1",
            "N#Cc1ccccc1",
            "Oc1ccccc1",
            "CS(C)=O",
            "CCO",
            "CC(=O)C",
            "CSC",
            "CCS",
            "CC[N+](=O)[O-]",
            "CCCl",
            "CCBr",
            "CCN",
            "CC#N",
        )
        observed = set()
        for smiles in representative_smiles:
            molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
            applications = find_lbcc_applications(molecule)
            observed.update(application.rule for application in applications)
            corrected, deltas = apply_lbcc(
                [0.0] * molecule.GetNumAtoms(), applications
            )
            self.assertAlmostEqual(sum(deltas), 0.0, places=12)
            self.assertAlmostEqual(sum(corrected), 0.0, places=12)

        self.assertEqual(len(LBCC_RULES), 19)
        self.assertEqual(observed, {rule.name for rule in LBCC_RULES})


if __name__ == "__main__":
    unittest.main()
