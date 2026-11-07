"""Chemical-context keys of bonded terms (shared torsion parameters in qm_fit)."""

from __future__ import annotations

import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from opls_ff_gen import load_sdf
from opls_ff_gen.assignment import _assign_atom_types
from opls_ff_gen.context import context_key, context_smarts

MOLECULES = Path(__file__).resolve().parents[1] / "molecules_to_test"


def _molecule(smiles):
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(molecule, randomSeed=1)
    return molecule


class ContextTests(unittest.TestCase):
    def test_context_separates_triazole_from_imidazole(self):
        triazole, imidazole = _molecule("CCn1cncn1"), _molecule("CCn1ccnc1")
        # quartet C-C-N1-ring C in both
        t = context_smarts(triazole, (0, 1, 2, 3))
        i = context_smarts(imidazole, (0, 1, 2, 3))
        self.assertNotEqual(t, i)
        self.assertFalse(imidazole.HasSubstructMatch(Chem.MolFromSmarts(t)))

    def test_key_is_orientation_independent(self):
        isoniazid = load_sdf(MOLECULES / "isoniazid.sdf")
        types, parameters, _ = _assign_atom_types(isoniazid)
        quartet = (1, 8, 9, 10)  # C2(=O) - N9(H) - N10(H2) - H11
        classes = tuple(parameters[types[i]]["class"] for i in quartet)
        forward = context_key(isoniazid, quartet, classes)
        reverse = context_key(isoniazid, quartet[::-1], classes[::-1])
        self.assertEqual(forward, reverse)
        # classes and SMARTS map numbers run the same way
        oriented = quartet if forward[0] == classes else quartet[::-1]
        self.assertEqual(forward[1], context_smarts(isoniazid, oriented))

    def test_nitro_oxygens_share_context(self):
        nitro = _molecule("Cn1cncc1[N+](=O)[O-]")
        ring_c, n, o1, o2 = 5, 6, 7, 8
        neighbor = next(a.GetIdx() for a in nitro.GetAtomWithIdx(ring_c).GetNeighbors()
                        if a.GetIdx() != n and a.GetAtomicNum() > 1)
        self.assertEqual(context_smarts(nitro, (neighbor, ring_c, n, o1)),
                         context_smarts(nitro, (neighbor, ring_c, n, o2)))

    def test_context_independent_of_input_route(self):
        # caffeine: the 6-ring is non-aromatic from SMILES, aromatic via SDF
        molecule = _molecule("Cn1c(=O)c2c(ncn2C)n(C)c1=O")
        reread = Chem.MolFromMolBlock(Chem.MolToMolBlock(molecule), removeHs=False)
        quartet = (0, 1, 2, 4)
        self.assertEqual(context_smarts(molecule, quartet), context_smarts(reread, quartet))


if __name__ == "__main__":
    unittest.main()
