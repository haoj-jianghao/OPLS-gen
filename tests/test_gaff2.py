"""GAFF2 rule table (replaces antechamber) against stored antechamber types."""

from __future__ import annotations

import json
import unittest
from itertools import combinations
from pathlib import Path

from rdkit import Chem

from opls_ff_gen import AssignmentError
from opls_ff_gen.amber import gaff2_lennard_jones, gaff2_types, load_gaff2_table
from opls_ff_gen.gaff2_bonded import load_gaff2_bonded
from opls_ff_gen.structure import apply_typing_aromaticity

REFERENCE = Path(__file__).resolve().parent / "data" / "gaff2_antechamber_reference.json"
# Atoms where antechamber is known to differ (see src/scripts/build_gaff2_table.py):
# cyclic N-nitrosamine N (antechamber n3, acyclic nitrosamines nh).
KNOWN_DIFFERENCES = {("first_5K:338", 3)}


def _lj(table, gaff_type):
    sigma, epsilon = table[gaff_type]
    return round(sigma, 6), round(epsilon, 6)


class Gaff2RuleTests(unittest.TestCase):
    def test_rules_reproduce_antechamber_lennard_jones(self):
        reference = json.loads(REFERENCE.read_text())["molecules"]
        table = gaff2_lennard_jones()
        atoms, wrong = 0, []
        for entry in reference:
            molecule = Chem.AddHs(Chem.MolFromSmiles(entry["smiles"]))
            expected = entry["types"].split()
            self.assertEqual(molecule.GetNumAtoms(), len(expected), entry["name"])
            ours = gaff2_types(molecule)
            for index, ref in enumerate(expected):
                atoms += 1
                if _lj(table, ours[index]) != _lj(table, ref) and (
                    (entry["name"], index + 1) not in KNOWN_DIFFERENCES
                ):
                    wrong.append((entry["name"], index + 1, ref, ours[index]))
        self.assertGreater(atoms, 15000)
        self.assertEqual(wrong, [])

    def test_every_rule_type_has_lennard_jones_parameters(self):
        table = load_gaff2_table()
        self.assertTrue({rule.label for rule in table.rules} <= table.lennard_jones.keys())
        sigma, epsilon = table.lennard_jones["c3"]  # gaff2.dat 2.2.30: c3 1.9069 0.1078
        self.assertAlmostEqual(sigma, 2 * 1.9069 / 2 ** (1 / 6) / 10)
        self.assertAlmostEqual(epsilon, 0.1078 * 4.184)

    def test_selected_atoms_only_and_unknown_elements(self):
        molecule = Chem.AddHs(Chem.MolFromSmiles("CC(=O)NC"))
        self.assertEqual(gaff2_types(molecule, [3]), {3: "ns"})
        with self.assertRaisesRegex(AssignmentError, "No GAFF2 rule"):
            gaff2_types(Chem.AddHs(Chem.MolFromSmiles("C[Si](C)(C)C")))


class Gaff2BondedTypeTests(unittest.TestCase):
    # antechamber re-perceives the resonance of cationic dyes; cyclic
    # N-nitrosamine N; a tetrazole imine and a diazo ester (see build_gaff2_bonded.py)
    KNOWN = {"first_5K:4192", "first_5K:4382", "first_5K:3574", "first_5K:338",
             "diazomethane_ester"}

    def test_full_types_give_antechamber_bond_and_angle_terms(self):
        table = load_gaff2_bonded()
        wrong, terms = [], 0
        for entry in json.loads(REFERENCE.read_text())["molecules"]:
            expected = entry["types"].split()
            molecule = apply_typing_aromaticity(Chem.AddHs(Chem.MolFromSmiles(entry["smiles"])))
            ours = table.types(molecule)
            self.assertNotIn(None, ours, entry["name"])
            atoms = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in molecule.GetBonds()]
            atoms += [(a.GetIdx(), center.GetIdx(), c.GetIdx()) for center in molecule.GetAtoms()
                      for a, c in combinations(center.GetNeighbors(), 2)]
            for term in atoms:
                lookup = table.bond if len(term) == 2 else table.angle
                reference = lookup(tuple(expected[i] for i in term))
                if reference is None:
                    continue
                terms += 1
                got = lookup(tuple(ours[i] for i in term))
                if (got is None or got.values != reference.values) and entry["name"] not in self.KNOWN:
                    wrong.append((entry["name"], term))
        self.assertGreater(terms, 45000)
        self.assertEqual(wrong, [])

    def test_conjugated_pairs_alternate(self):
        table = load_gaff2_bonded()

        def heavy(smiles):
            molecule = apply_typing_aromaticity(Chem.AddHs(Chem.MolFromSmiles(smiles)))
            return [t for t, a in zip(table.types(molecule), molecule.GetAtoms())
                    if a.GetAtomicNum() > 1]

        # which atom of a system gets the first name is arbitrary (same terms)
        furan = heavy("c1ccoc1")  # C0-C1 single, C1=C2, C0=C4
        self.assertEqual(furan[0], furan[1])
        self.assertEqual({furan[1], furan[2]}, {"cc", "cd"})
        self.assertEqual(furan[2], furan[4])
        hexatriene = heavy("C=CC=CC=C")
        self.assertEqual(hexatriene[0], "c2")
        self.assertEqual(hexatriene[1], hexatriene[2])
        self.assertEqual({hexatriene[2], hexatriene[3]}, {"ce", "cf"})
        self.assertEqual(hexatriene[3], hexatriene[4])
        self.assertEqual(heavy("c1ccccc1-c1ccccc1").count("cp"), 2)
        self.assertEqual(heavy("c1ccc2ccccc2c1"), ["ca"] * 10)  # fused: not cp

    def test_units_and_lookup(self):
        table = load_gaff2_bonded()
        bond = table.bond(("c3", "c3"))  # gaff2.dat: 228.89 1.5354
        self.assertAlmostEqual(bond.values[0], 0.15354)
        self.assertAlmostEqual(bond.values[1], 2 * 228.89 * 4.184 * 100)
        angle = table.angle(("c3", "c3", "c3"))  # 59.87 112.63
        self.assertEqual(angle.values, (112.63, 2 * 59.87 * 4.184))
        self.assertEqual(table.bond(("cd", "cc")).classes, ("cc", "cd"))


if __name__ == "__main__":
    unittest.main()
