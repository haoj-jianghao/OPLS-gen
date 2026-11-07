"""Structure input (SMILES, SDF, MOL2, PDB) and geometry optimization."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np
from rdkit import Chem

from opls_ff_gen import AssignmentError, load_sdf, molecule_from_smiles
from opls_ff_gen.structure import apply_typing_aromaticity, load_structure, prepare_structure
from opls_ff_gen.xtb import XTB

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data"


def _smiles(molecule):
    return Chem.MolToSmiles(Chem.RemoveHs(molecule))


class StructureInputTests(unittest.TestCase):
    def setUp(self):
        self.deet = load_sdf(ROOT / "molecules_to_test" / "deet.sdf")
        self.metronidazole = load_sdf(ROOT / "molecules_to_test" / "Metronidazole.sdf")
        self.scratch = tempfile.TemporaryDirectory()
        self.tmp = Path(self.scratch.name)

    def tearDown(self):
        self.scratch.cleanup()

    def test_sdf_and_mol2_give_the_same_molecule(self):
        sdf = prepare_structure(path=ROOT / "molecules_to_test" / "deet.sdf")
        mol2 = prepare_structure(path=DATA / "deet.mol2")
        self.assertEqual(_smiles(sdf), _smiles(mol2))
        np.testing.assert_allclose(sdf.GetConformer().GetPositions(),
                                   mol2.GetConformer().GetPositions(), atol=1e-3)
        self.assertEqual(mol2.GetProp("Geometry_Optimization"), "none (input coordinates)")

    def test_mol2_nitro_needs_template_and_net_charge_catches_it(self):
        with self.assertRaisesRegex(AssignmentError, "formal charges sum to -2"):
            load_structure(DATA / "Metronidazole.mol2", net_charge=0)
        fixed = load_structure(DATA / "Metronidazole.mol2", net_charge=0,
                               template_smiles=_smiles(self.metronidazole))
        self.assertEqual(_smiles(fixed), _smiles(self.metronidazole))

    def test_pdb_bond_orders_perceived_or_from_template(self):
        pdb = self.tmp / "deet.pdb"
        Chem.MolToPDBFile(self.deet, str(pdb))
        self.assertEqual(_smiles(load_structure(pdb, net_charge=0)), _smiles(self.deet))
        templated = load_structure(pdb, template_smiles=_smiles(self.deet))
        self.assertEqual(templated.GetProp("Input_Bond_Orders"), "template SMILES")
        self.assertEqual(_smiles(templated), _smiles(self.deet))

    def test_missing_hydrogens_need_optimization(self):
        pdb = self.tmp / "deet_noH.pdb"
        Chem.MolToPDBFile(Chem.RemoveHs(self.deet), str(pdb))
        with self.assertRaisesRegex(AssignmentError, "hydrogen atom"):
            prepare_structure(path=pdb, template_smiles=_smiles(self.deet))
        molecule = prepare_structure(path=pdb, template_smiles=_smiles(self.deet), optimize="mmff")
        self.assertEqual(molecule.GetNumAtoms(), self.deet.GetNumAtoms())

    def test_smiles_always_optimized_and_mmff_matches_previous_path(self):
        with self.assertRaisesRegex(AssignmentError, "always optimized"):
            prepare_structure(smiles="CCO", optimize="none")
        new = prepare_structure(smiles="CCN(CC)C(=O)c1cccc(C)c1", name="deet", optimize="mmff")
        old = molecule_from_smiles("CCN(CC)C(=O)c1cccc(C)c1", name="deet")
        np.testing.assert_array_equal(new.GetConformer().GetPositions(),
                                      old.GetConformer().GetPositions())

    @unittest.skipUnless(XTB.is_file(), "xtb unavailable")
    def test_smiles_default_to_xtb(self):
        molecule = prepare_structure(smiles="CCO", name="ethanol")
        self.assertTrue(molecule.GetProp("Geometry_Optimization").startswith("GFN2-xTB"))

    def test_file_and_smiles_routes_share_the_aromaticity_model(self):
        # caffeine's pyrimidinedione ring is aromatic in RDKit's default model
        # (as read from the SDF) but not in MMFF94; both routes use MMFF94.
        from opls_ff_gen.assignment import smarts_atom_types

        from_file = prepare_structure(path=ROOT / "molecules_to_test" / "caffeine.sdf")
        from_smiles = prepare_structure(smiles=_smiles(from_file), name="caffeine",
                                        optimize="mmff")
        self.assertEqual(sorted(smarts_atom_types(from_file)[0], key=str),
                         sorted(smarts_atom_types(from_smiles)[0], key=str))
        self.assertEqual([a.GetIsAromatic() for a in from_file.GetAtoms()],
                         [a.GetIsAromatic() for a in apply_typing_aromaticity(Chem.Mol(from_file))
                          .GetAtoms()])

    @unittest.skipUnless(XTB.is_file(), "xtb unavailable")
    def test_smiles_typing_does_not_depend_on_the_optimizer(self):
        # RDKit's default model makes the 4-hydroxycoumarin ring aromatic; the
        # SMILES route uses MMFF94 aromaticity with either optimizer.
        from opls_ff_gen.assignment import smarts_atom_types

        smiles = "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O"  # warfarin
        types = [smarts_atom_types(prepare_structure(smiles=smiles, name="w", optimize=method))[0]
                 for method in ("mmff", "xtb")]
        self.assertEqual(types[0], types[1])
        self.assertEqual(types[1][11:13], ["opls_141", "opls_518"])

    @unittest.skipUnless(XTB.is_file(), "xtb unavailable")
    def test_xtb_optimization_keeps_bonding(self):
        molecule = prepare_structure(path=ROOT / "molecules_to_test" / "deet.sdf", optimize="xtb")
        self.assertTrue(molecule.GetProp("Geometry_Optimization").startswith("GFN2-xTB"))
        self.assertEqual(_smiles(molecule), _smiles(self.deet))
        moved = np.abs(molecule.GetConformer().GetPositions()
                       - self.deet.GetConformer().GetPositions()).max()
        self.assertGreater(moved, 0.001)


    @unittest.skipUnless(XTB.is_file(), "xtb unavailable")
    def test_xtb_bonding_change_falls_back_to_mmff(self):
        # gas-phase GFN2-xTB moves an ammonium H to the carboxylate of zwitterionic glycine
        molecule = prepare_structure(smiles="[NH3+]CC(=O)[O-]", name="glycine")
        self.assertEqual(molecule.GetProp("Geometry_Optimization"), "MMFF94")
        self.assertIn("changed the bonding", molecule.GetProp("Geometry_Fallback"))
        self.assertEqual(Chem.MolToSmiles(Chem.RemoveHs(molecule)), "[NH3+]CC(=O)[O-]")


    def test_missing_xtb_falls_back_to_mmff(self):
        from unittest import mock

        with mock.patch("opls_ff_gen.xtb.XTB", Path("/nonexistent/xtb")):
            molecule = prepare_structure(smiles="CCO", name="ethanol")
        self.assertEqual(molecule.GetProp("Geometry_Optimization"), "MMFF94")
        self.assertIn("xtb not found", molecule.GetProp("Geometry_Fallback"))


    def test_failing_xtb_falls_back_to_mmff(self):
        from unittest import mock

        with mock.patch("opls_ff_gen.xtb.xtb_available", return_value=True), \
             mock.patch("opls_ff_gen.structure._optimize_xtb",
                        side_effect=RuntimeError("xtb failed in /tmp/x")):
            molecule = prepare_structure(smiles="CCO", name="ethanol")
        self.assertEqual(molecule.GetProp("Geometry_Optimization"), "MMFF94")
        self.assertIn("xtb optimization failed", molecule.GetProp("Geometry_Fallback"))


if __name__ == "__main__":
    unittest.main()
