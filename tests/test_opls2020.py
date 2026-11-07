"""OPLS/2020 parameter set (default), GAFF2 fallback and on-the-fly QM fit."""

from __future__ import annotations

import json
import subprocess
import tempfile
import importlib.util
import unittest
from dataclasses import replace
from pathlib import Path
from unittest import mock

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from opls_ff_gen import AssignmentError, load_sdf, molecule_from_smiles
from opls_ff_gen.assignment import (
    _assign_atom_types, _smarts_atom_types, load_atomtype_archive, smarts_atom_types)
from opls_ff_gen.bonded import (
    RULE_GAFF2_IMPROPER,
    RULE_IMPROPER,
    RULE_QM_PENDING,
    assign_bonded_interactions,
)
from opls_ff_gen.opls2020 import load_atom_types, load_bonded_table
from opls_ff_gen.parmchk import RULE_CORR, RULE_CORR_GENERAL, load_parmchk_scores
from opls_ff_gen.parmchk import RULES as PARMCHK_RULES
from opls_ff_gen.gaff2_torsions import RULE_GAFF2_GENERAL, load_gaff2_general
from opls_ff_gen.gaff2_bonded import RULE_GAFF2_BONDED, load_gaff2_bonded
from opls_ff_gen.qm_bonded import RULE_QM_TABLE, angle_ring_context, load_qm_bonded
from opls_ff_gen.bonded_substitution import (
    RULE_CENTRAL_PAIR, RULE_DEFAULT_GENERAL, RULE_EMPIRICAL, RULE_GAFF2_DAT_GENERAL,
    RULE_ZERO_TORSION, default_class, load_bonded_substitution)
from opls_ff_gen.qm_fit import (
    RULE_QM_FIT,
    RULE_ZERO_LINEAR,
    cache_inputs,
    prepare_cache,
    seminario_bond,
)
from opls_ff_gen.gromacs import gmx_executable
from opls_ff_gen.qm import DEFAULT_DFT
from opls_ff_gen.xtb import XTB
from opls_ff_gen.workflow import generate_gromacs_topology

ROOT = Path(__file__).resolve().parents[1]
MOLECULES = ROOT / "molecules_to_test"
GMX = gmx_executable() or Path("gmx-not-found")


def _types(smiles: str):
    molecule = molecule_from_smiles(smiles)
    types, parameters, rules = _assign_atom_types(molecule)
    heavy = [i for i, atom in enumerate(molecule.GetAtoms()) if atom.GetAtomicNum() > 1]
    return [types[i] for i in heavy], parameters, [rules[i] for i in heavy]


class ParameterTableTests(unittest.TestCase):
    def test_values_are_converted_to_gromacs_units(self):
        table = load_bonded_table()
        # SI Table S3: CT-CT 1.529 A, 268 kcal/mol/A^2; CT-CT-CT 112.7 deg, 58.35
        np.testing.assert_allclose(table.bond(("CT", "CT")).values, (0.1529, 268 * 836.8))
        theta, k = table.angle(("CT", "CT", "CT")).values
        self.assertEqual(theta, 112.7)
        self.assertAlmostEqual(k, 58.35 * 8.368)
        benzene_c = load_atom_types()[145]
        self.assertEqual((benzene_c.opls_class, benzene_c.sigma_angstrom,
                          benzene_c.epsilon_kcal_mol), ("CA", 3.55, 0.068))
        # O-C-X-Y V2 = 21 kcal/mol -> k = 43.932 kJ/mol, phase 180, n = 2
        self.assertEqual(table.improper("improper_O_C_X_Y").values, (180.0, 43.932, 2.0))
        # allene twist V2 = -8 kcal/mol -> phase 0, k = 16.736 kJ/mol, n = 2
        self.assertEqual(table.allene_twist("HC", "CT").values, (0.0, 16.736, 2.0))
        self.assertEqual(table.allene_twist("CT", "HC").values, table.allene_twist("HC", "CT").values)
        self.assertIsNone(table.allene_twist("CA", "HC"))

    def test_torsion_lookup_follows_boss_priority(self):
        table = load_bonded_table()
        # first entry wins: general OPLS/2020 CT-CT-CT-CT (003), not "butane only" (004)
        general = table.proper(("CT", "CT", "CT", "CT"))
        self.assertEqual(table.torsion_type(general.source_line), "003")
        # reversed quartet; first entry is the trans-favouring tertiary-amide term
        amide = table.proper(("CT", "CT", "N", "CT"))
        self.assertEqual(table.torsion_type(amide.source_line), "108")
        # one-character wildcard: CA-CA-CA-N? matches CA-CA-CA-N
        found, rule = table.lookup("proper", ("CA", "CA", "CA", "N"))
        self.assertEqual((found.classes[3], rule), ("N?", "wildcard"))
        # an exact entry beats the ??-CA-CA-?? wildcard
        _, rule = table.lookup("proper", ("HA", "CA", "CA", "HA"))
        self.assertEqual(rule, "exact")


class TypingTests(unittest.TestCase):
    def test_alkane_types_are_opls2020(self):
        types, parameters, rules = _types("CCCC(C)C")
        self.assertEqual(types, ["opls_54", "opls_57", "opls_57", "opls_58", "opls_55", "opls_55"])
        self.assertEqual(rules[0], "smarts:opls_135")
        self.assertAlmostEqual(parameters["opls_55"]["sigma"], 0.340)

    def test_cyclobutane_and_alkyne_neighbours(self):
        types, parameters, _ = _types("C1CCC1")
        self.assertEqual(set(types), {"opls_714"})
        self.assertEqual(parameters["opls_714"]["class"], "CY")
        types, _, _ = _types("CC#CC")
        self.assertEqual(types, ["opls_927", "opls_939", "opls_939", "opls_927"])

    def test_allene_types(self):
        types, parameters, _ = _types("C=C=C")
        self.assertEqual(types, ["opls_1154", "opls_1157", "opls_1154"])
        self.assertEqual(parameters["opls_1157"]["class"], "C:")
        types, _, rules = _types("CC(C)=C=C(C)C")  # Foyer leaves this centre untyped
        self.assertEqual(types[3], "opls_1157")
        self.assertEqual(rules[3], "archive:opls2020")
        types, _, _ = _types("CC=C=C")
        self.assertEqual(types[1:], ["opls_1155", "opls_1157", "opls_1154"])
        self.assertEqual(_types("C=C=O")[0][1], "gaff2_c1")  # ketene is not an allene

    def test_untypable_atoms_get_gaff2(self):
        types, parameters, rules = _types("Cc1nnc(s1)N")
        gaff = [t for t in types if t.startswith("gaff2_")]
        self.assertEqual(gaff, ["gaff2_nb", "gaff2_nb"])
        self.assertEqual(parameters["gaff2_nb"]["class"], "nb")
        self.assertIn("amber:gaff2:nb", rules)
        with self.assertRaises(AssignmentError):
            _assign_atom_types(molecule_from_smiles("Cc1nnc(s1)N"), amber_fallback=False)

    def test_opls2020_types_without_foyer_smarts_are_borrowed(self):
        # heavy-atom types; these parents had GAFF2 atoms before the 2026-10-29 rules
        expected = {
            "c1cn[nH]c1": ["opls_551", "opls_550", "opls_549", "opls_548", "opls_552"],
            "c1ccnnc1": ["opls_539", "opls_539", "opls_538", "opls_537", "opls_537", "opls_538"],
            "c1cnccn1": ["opls_528", "opls_528", "opls_527", "opls_528", "opls_528", "opls_527"],
            "c1ncncn1": ["opls_642", "opls_641"] * 3,
            "CC1CO1": ["opls_55", "opls_1027", "opls_1026", "opls_1025"],
            "C[N+](C)(C)C": ["opls_291", "opls_288", "opls_291", "opls_291", "opls_291"],
            "CC(Br)C": ["opls_55", "opls_978", "opls_975", "opls_55"],
            "CN(C)c1ccccc1": ["opls_905", "opls_902", "opls_905", "opls_918"] + ["opls_145"] * 5,
            "CS(N)(=O)=O": ["opls_476", "opls_474", "opls_478", "opls_475", "opls_475"],
            "CC1(C)OCCO1": ["opls_56", "opls_197", "opls_56", "opls_186", "opls_182",
                            "opls_182", "opls_186"],
            # guanidinium cation ("*" refinement rules: the rule set types its N as amines)
            "CNC(=[NH2+])N": ["opls_305", "opls_303", "opls_302", "opls_300", "opls_300"],
            "NC(=[NH2+])Nc1ccccc1": ["opls_300", "opls_302", "opls_300", "opls_303", "opls_793"]
                                    + ["opls_145"] * 5,
        }
        for smiles, heavy in expected.items():
            types, _, rules = _types(smiles)
            self.assertEqual(types, heavy, smiles)
            self.assertFalse(any(r.startswith(("amber:", "analog:")) for r in rules), smiles)
        molecule = molecule_from_smiles("c1cn[nH]c1")  # the ring H get the pyrazole H types
        types, _, _ = _assign_atom_types(molecule)
        self.assertEqual(sorted(types[5:]), ["opls_553", "opls_554", "opls_555", "opls_556"])

    def test_gaff2_types_with_an_opls2020_analog(self):
        molecule = molecule_from_smiles("CS(=O)(=O)O")  # methanesulfonic acid
        types, parameters, rules = _assign_atom_types(molecule)
        self.assertEqual(types[:5], ["opls_135", "opls_493", "opls_494", "opls_494", "opls_268"])
        self.assertEqual(types[8], "opls_270")
        self.assertEqual(rules[1], "analog:sulfonyl-S")
        self.assertEqual(parameters["opls_493"]["class"], "SY")
        off, _, _ = _assign_atom_types(molecule, gaff2_analogs=False)
        self.assertEqual(off[:5], ["gaff2_c3", "gaff2_s6", "gaff2_o", "gaff2_o", "gaff2_oh"])
        # thiosemicarbazone: C=S, hydrazide N and C=N keep GAFF2 (with the N-H)
        molecule = molecule_from_smiles("NC(=S)NN=Cc1ccccc1")
        types, _, _ = _assign_atom_types(molecule)
        self.assertEqual(types[1:6], ["gaff2_c", "gaff2_s", "gaff2_ns", "gaff2_n2", "gaff2_c2"])
        self.assertEqual(types[14], "gaff2_hn")
        # dithiocarbamate: amide-type N analog next to the GAFF2 C=S
        types, _, rules = _types("CCN(CC)C(=S)SC")
        self.assertEqual((types[2], rules[2]), ("opls_239", "analog:amide-N-H0"))
        self.assertEqual(types[5:7], ["gaff2_c", "gaff2_s"])
        # 5-ring heteroaromatics have no analog
        self.assertEqual(_types("Cc1nnc(s1)N")[0].count("gaff2_nb"), 2)

    def test_analog_atoms_keep_their_gaff2_type_after_the_table(self):
        from opls_ff_gen.gaff2_analogs import analog_gaff2_classes

        molecule = molecule_from_smiles("CS(=O)(=O)OC")  # methyl methanesulfonate
        types, parameters, rules = _assign_atom_types(molecule)
        classes = [parameters[t]["class"] for t in types]
        gaff2 = analog_gaff2_classes(molecule, rules)
        self.assertEqual(gaff2[:6], ["c3", "s6", "o", "o", "os", "c3"])
        rules_of = {}
        for label, value in (("opls", None), ("gaff2", gaff2)):
            bonded = assign_bonded_interactions(molecule, classes, gaff2_classes=value)
            rules_of[label] = {x.atoms: x.rule for x in (*bonded.angles, *bonded.proper_torsions)}
            # the OPLS/2020 table is still first: S=O and C-S terms are table terms
            self.assertEqual(rules_of[label][(3, 2, 4)], "exact")
        # CT-SY-OS and SY-OS-CT are not in the table: gaff2.dat with the GAFF2 types
        self.assertEqual([rules_of["opls"][a] for a in ((1, 2, 5), (2, 5, 6))], ["empirical"] * 2)
        self.assertEqual([rules_of["gaff2"][a] for a in ((1, 2, 5), (2, 5, 6))], ["gaff2"] * 2)
        self.assertEqual(rules_of["gaff2"][(1, 2, 5, 6)], RULE_GAFF2_GENERAL)

    def test_ring_only_general_torsion(self):
        # X-CA-CU-X (xTB refit on fused pyrazolo-pyrimidines) is for ring bonds only
        general = load_gaff2_general()
        ring = general.lookup(("CA", "CA", "CU", "NB"), ring=True)
        self.assertTrue(ring.ring_only)
        self.assertEqual(ring.classes, ("X", "CA", "CU", "X"))
        self.assertIsNone(general.lookup(("CA", "CA", "CU", "NB")))
        for smiles, rule in (("Oc1ncnc2[nH]ncc12", RULE_GAFF2_GENERAL),  # fused ring bond
                             ("c1ccc(cc1)-c1cc[nH]n1", RULE_CORR_GENERAL)):  # aryl-pyrazole
            molecule = molecule_from_smiles(smiles)
            types, parameters, _ = _assign_atom_types(molecule)
            bonded = assign_bonded_interactions(molecule, [parameters[t]["class"] for t in types])
            about = {t.rule for t in bonded.proper_torsions if sorted(t.classes[1:3]) == ["CA", "CU"]}
            self.assertEqual(about, {rule}, smiles)

    def test_gaff2_atoms_get_cosmetic_opls_names_in_the_itp(self):
        from opls_ff_gen.amber import gaff2_display_types, load_gaff2_table

        names = gaff2_display_types()
        self.assertEqual(set(names), set(load_gaff2_table().lennard_jones))
        self.assertEqual(len(set(names.values())), len(names))
        numbers = {int(name.removeprefix("opls_")) for name in names.values()}
        self.assertTrue(all(800 < n < 900 for n in numbers))
        self.assertFalse(numbers & set(load_atom_types()))  # never an OPLS/2020 number
        with tempfile.TemporaryDirectory() as directory:
            molecule = molecule_from_smiles("CC(C)=NO", name="oxime")  # C=N-OH stays GAFF2
            generate_gromacs_topology(molecule, directory, name="oxime")
            itp = (Path(directory) / "oxime.itp").read_text()
        section = itp.split("[ atomtypes ]")[1].split("[ moleculetype ]")[0]
        bond_types = [line.split()[1] for line in section.splitlines()
                      if line.strip() and not line.startswith(";")]
        self.assertEqual([bond_types[i] for i in (1, 3, 4)], [names["c2"], names["n2"], names["oh"]])
        self.assertFalse(set(bond_types) & set(names))  # no GAFF2 type names
        self.assertIn(f"; gaff2-lj {names['n2']}", itp)
        self.assertIn(f"{names['n2']} = GAFF2 n2", itp)
        self.assertNotIn("amber:gaff2", itp)

    def test_rule_set_ambiguities_are_resolved(self):
        types, _, _ = _types("O=c1[nH]c2cc(Cl)ccc2o1")  # chlorzoxazone
        self.assertIn("opls_467", types)  # ring O as ester/carbamate OS
        types, parameters, _ = _types("Cn1c(=O)c2[nH]cnc2n(C)c1=O")  # theophylline
        self.assertIn("opls_628", types)  # N7-H is purine NA, not an amine
        self.assertEqual(parameters["opls_628"]["class"], "NA")

    def test_fused_heteroaromatic_and_alpha_precedence(self):
        # Fused rings: the pyridine/pyrimidine rule over naphthalene fusion C (opls_147)
        types, parameters, _ = _types("c1ccc2ncccc2c1")  # quinoline
        self.assertEqual(types[3], "opls_521")  # C8a, bonded to N
        self.assertEqual(types[8], "opls_522")  # C4a
        self.assertNotIn("opls_147", types)
        types, _, _ = _types("Nc1ncnc2[nH]cnc12")  # adenine: C6 is pyrimidine C4, not 917
        self.assertEqual(types[1], "opls_532")
        self.assertEqual(types[9], "opls_533")
        # Carbon groups: the nitrile alpha-H over the carbonyl alpha-H
        molecule = molecule_from_smiles("COC(=O)CC#N")
        types, _, _ = _assign_atom_types(molecule)
        self.assertEqual(types[4], "opls_756")
        self.assertEqual(types[10:], ["opls_759", "opls_759"])
        types, _, _ = _types("O=C1C=CC(=O)C=C1")  # benzoquinone C=O is a carbonyl C
        self.assertEqual(types[1], "opls_233")

    def test_fused_cycloalkene_is_not_ring_fusion_c(self):
        types, _, _ = _types("C1=CC2CCC1CC2")  # bicyclo[2.2.2]octene
        self.assertEqual(types[:2], ["opls_142", "opls_142"])
        types, _, _ = _types("C1CCC2=C(C1)CCCC2")  # octalin
        self.assertEqual(types[3:5], ["opls_141", "opls_141"])
        types, _, _ = _types("c1ccc2ccccc2c1")  # naphthalene keeps opls_147
        self.assertEqual(types[3], "opls_147")
        types, _, _ = _types("CC1=C2C=CC=CC2=CCN1")  # 145 also matches: keeps 147
        self.assertEqual([types[2], types[7]], ["opls_147", "opls_147"])

    def test_collision_without_precedence_takes_highest_number(self):
        ergoline = "CC1Cc2c[nH]c3cccc(c23)C1"  # indole C7a/C3a match 145 and 543/544
        types, _, rules = _types(ergoline)
        self.assertEqual(types[6] + types[11], "opls_543opls_544")  # precedence entries
        self.assertEqual({rules[6], rules[11]}, {"smarts"})
        archive = load_atomtype_archive()
        precedence = {k: v - {"opls_145"} if k in ("opls_543", "opls_544") else v
                      for k, v in archive.precedence.items()}
        molecule = molecule_from_smiles(ergoline)
        types, collisions = _smarts_atom_types(molecule, replace(archive, precedence=precedence))
        self.assertEqual((types[6], types[11]), ("opls_543", "opls_544"))
        self.assertEqual(collisions, {6: ("opls_145", "opls_543"), 11: ("opls_145", "opls_544")})

    def test_silicon_types(self):
        def all_types(smiles):
            return _assign_atom_types(molecule_from_smiles(smiles))[0]

        types = all_types("C[Si](C)(C)OC")  # methoxytrimethylsilane
        self.assertEqual(types[:6], ["opls_1065", "opls_1075", "opls_1065", "opls_1065",
                                     "opls_1078", "opls_181"])
        self.assertEqual(set(types[15:]), {"opls_185"})  # H alpha to the silyl ether O
        types = all_types("C[SiH](C)C")
        self.assertEqual((types[1], types[7], types[4]), ("opls_1061", "opls_1064", "opls_60"))
        types = all_types("C[Si](C)(C)O")
        self.assertEqual((types[1], types[4], types[-1]), ("opls_1070", "opls_1073", "opls_1074"))
        types = all_types("C[Si](C)(C)c1ccccc1")
        self.assertEqual(types[4], "opls_1069")

    def test_halogen_types(self):
        cases = {  # SMILES -> type of the (first) halogen, heavy-atom order
            "ClC=C": "opls_226", "ClC(Cl)=C": "opls_398", "CC(=O)Cl": "opls_981",
            "Clc1ncccn1": "opls_264", "ClC(Cl)(Cl)Cl": "opls_123", "ClC(Cl)Cl": "opls_121",
            "ClCCl": "opls_119", "CCCCCl": "opls_151", "CNCl": "gaff2_cl",
            "Brc1ccccc1": "opls_730", "CCBr": "opls_975", "CC(=O)Br": "opls_982",
            "CCCI": "opls_1014", "Ic1ccccc1": "opls_732",
            "CCCF": "opls_956", "FC(F)F": "opls_965", "CC(=O)F": "opls_980",
            "Fc1ccncc1": "opls_719", "Fc1ccccc1": "opls_719",
        }
        for smiles, expected in cases.items():
            types, parameters, _ = _types(smiles)
            molecule = Chem.MolFromSmiles(smiles)
            index = next(a.GetIdx() for a in molecule.GetAtoms()
                         if a.GetSymbol() in ("F", "Cl", "Br", "I"))
            with self.subTest(smiles=smiles):
                self.assertEqual(types[index], expected)
        types, _ = smarts_atom_types(Chem.MolFromSmiles("[Cl-]"))
        self.assertEqual(types, ["opls_401"])  # the ion keeps the ion type


class BondedTests(unittest.TestCase):
    def test_missing_terms_are_listed_or_deferred(self):
        molecule = load_sdf(MOLECULES / "deet.sdf")
        types, parameters, _ = _assign_atom_types(molecule)
        classes = [parameters[t]["class"] for t in types]
        with self.assertRaisesRegex(AssignmentError, "proper CA-C-N-CT: 2x"):
            assign_bonded_interactions(molecule, classes, torsion_substitution=False)
        bonded = assign_bonded_interactions(molecule, classes, defer_missing=True,
                                            torsion_substitution=False)
        pending = [t for t in bonded.proper_torsions if t.rule == RULE_QM_PENDING]
        self.assertEqual(len(pending), 2)
        self.assertEqual({t.classes for t in pending}, {("CA", "C", "N", "CT")})
        self.assertTrue(all(t.coefficients_kj_mol == (0.0,) * 6 for t in pending))


    def test_allene_twist_terms(self):
        molecule = molecule_from_smiles("C=C=C")
        types, parameters, _ = _assign_atom_types(molecule)
        bonded = assign_bonded_interactions(molecule, [parameters[t]["class"] for t in types])
        twists = [v for v in bonded.improper_torsions if v.parameter_name == "allene_twist"]
        self.assertEqual(sorted(v.atoms for v in twists),
                         [(4, 1, 3, 6), (4, 1, 3, 7), (5, 1, 3, 6), (5, 1, 3, 7)])
        self.assertTrue(all((v.phase_deg, v.multiplicity) == (0.0, 2) for v in twists))
        # a substituent class without a twist entry (CA) is reported
        molecule = molecule_from_smiles("C=C=Cc1ccccc1")
        types, parameters, _ = _assign_atom_types(molecule)
        bonded = assign_bonded_interactions(molecule, [parameters[t]["class"] for t in types],
                                            defer_missing=True)
        self.assertEqual(sum(v.parameter_name == "allene_twist" for v in bonded.improper_torsions), 2)
        self.assertTrue(any("no twist term for substituent classes" in w for w in bonded.warnings))

    def test_gaff2_centres_get_opls_impropers(self):
        # benzaldehyde thiosemicarbazone: C=S C (c), hydrazide N (ns), imine C (c2) stay GAFF2
        molecule = molecule_from_smiles("NC(=S)NN=Cc1ccccc1")
        types, parameters, _ = _assign_atom_types(molecule)
        bonded = assign_bonded_interactions(molecule, [parameters[t]["class"] for t in types],
                                            defer_missing=True)
        gaff = {v.atoms[1]: v.parameter_name for v in bonded.improper_torsions
                if v.rule == RULE_GAFF2_IMPROPER}
        self.assertEqual(gaff, {2: "improper_Z_CM_X_Y", 4: "improper_Z_N_X_Y",
                                6: "improper_Z_CM_X_Y"})

    def test_gaff2_aromatic_c_takes_the_opls_analog(self):
        # warfarin with an aromatic coumarin ring (as from an aromatic SDF): C3/C4 have
        # GAFF2 type ca, whose analog is benzene C, so they get the OPLS CA improper
        molecule = Chem.AddHs(Chem.MolFromSmiles("CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O"))
        AllChem.EmbedMolecule(molecule, randomSeed=7)
        types, parameters, rules = _assign_atom_types(molecule)
        self.assertEqual([(types[i], rules[i]) for i in (11, 12)],
                         [("opls_145", "analog:aromatic-C-6ring")] * 2)
        bonded = assign_bonded_interactions(molecule, [parameters[t]["class"] for t in types],
                                            defer_missing=True)
        self.assertFalse(any(v.rule == RULE_GAFF2_IMPROPER for v in bonded.improper_torsions))


class ParmchkSubstitutionTests(unittest.TestCase):
    """parmchk2 steps 5 and 6 with the uppercase PARMCHK.DAT types."""

    def test_penalty_scores_follow_parmchk2(self):
        table, scores = load_bonded_table(), load_parmchk_scores()
        cases = {
            # outer OS -> OH: tor 0.1
            ("CM", "CM", "C", "OS"): (RULE_CORR, ("CM", "CM", "C", "OH"), 0.1),
            # central CB -> CA: 10 x (0.5 x 19.9 + 0.5 x 19.8)
            ("C", "OS", "CB", "CA"): (RULE_CORR, ("CA", "CA", "OS", "C"), 198.5),
            # no specific entry: general X-C-CB-X for C-CA
            ("CA", "C", "CA", "CA"): (RULE_CORR_GENERAL, ("X", "C", "CB", "X"), 198.5),
        }
        for classes, (rule, source, penalty) in cases.items():
            hit = scores.substitute(table, classes)
            self.assertEqual(hit.rule, rule, classes)
            self.assertIn(hit.parameter.classes, (source, source[::-1]), classes)
            self.assertAlmostEqual(hit.penalty, penalty, places=6, msg=classes)

    def test_classes_without_correspondences_stay_missing(self):
        # C has no uppercase correspondence and OPLS/2020 no X-C-N-X term.
        self.assertIsNone(load_parmchk_scores().substitute(
            load_bonded_table(), ("CA", "C", "N", "CT")))

    def test_gaff2_general_terms_come_before_substitution(self):
        molecule = molecule_from_smiles("O=c1ccc2ccccc2o1", name="coumarin")
        types, parameters, _ = _assign_atom_types(molecule)
        classes = [parameters[t]["class"] for t in types]
        bonded = assign_bonded_interactions(molecule, classes)
        general = [t for t in bonded.proper_torsions if t.rule == RULE_GAFF2_GENERAL]
        self.assertTrue(general)
        for t in general:
            self.assertEqual(t.parameter_classes[0], "X")
            self.assertEqual(set(t.parameter_classes[1:3]), set(t.classes[1:3]))
        # CM-CM-C-OS: X-C-CM-X now, the parmchk2 substitution (CM-CM-C-OH) before
        lactone = [t for t in bonded.proper_torsions
                   if min(t.classes, t.classes[::-1]) == ("CM", "CM", "C", "OS")]
        self.assertEqual({t.rule for t in lactone}, {RULE_GAFF2_GENERAL})
        # with the run-time QM fit these torsions are fitted instead
        strict = assign_bonded_interactions(molecule, classes, defer_missing=True)
        self.assertFalse(any(t.rule in PARMCHK_RULES or t.rule == RULE_GAFF2_GENERAL
                             for t in strict.proper_torsions))
        self.assertEqual({tuple(t.atoms) for t in strict.proper_torsions
                          if t.rule == RULE_QM_PENDING},
                         {tuple(t.atoms) for t in bonded.proper_torsions
                          if t.rule not in ("exact", "wildcard")})

    def test_qm_check_flags(self):
        cases = {
            # GAFF2 has X-c2-no-X only by type substitution
            "Cc1ncc([N+](=O)[O-])n1CCO": "gaff2-general X-CW-NO-X",
            # OPLS/2020 X-CB-N-X kept, flagged (aryl-N)
            "Cn1c(=O)c2c(ncn2C)n(C)c1=O": "OPLS/2020 X-N-CB-X",
        }
        for smiles, entry in cases.items():
            molecule = molecule_from_smiles(smiles, name="m")
            types, parameters, _ = _assign_atom_types(molecule)
            bonded = assign_bonded_interactions(
                molecule, [parameters[t]["class"] for t in types])
            self.assertTrue(any(t.qm_check for t in bonded.proper_torsions), smiles)
            self.assertTrue(any(entry in w and "QM check" in w for w in bonded.warnings),
                            (smiles, bonded.warnings))


class BondedFallbackTests(unittest.TestCase):
    def test_missing_angle_takes_empirical_rule_without_qm(self):
        molecule = molecule_from_smiles("CSCC(=O)O")  # S-CT-C is not in OPLS/2020
        types, parameters, _ = _assign_atom_types(molecule)
        classes = [parameters[t]["class"] for t in types]
        bonded = assign_bonded_interactions(molecule, classes, qm_bonded_table=False)
        (angle,) = [a for a in bonded.angles if a.rule == RULE_EMPIRICAL]
        self.assertEqual(angle.classes, ("S", "CT", "C"))
        self.assertAlmostEqual(angle.angle_deg, 110.2, places=1)  # mean OPLS angle about CT
        self.assertGreater(angle.force_constant_kj_mol_rad2, 300)
        # the run-time QM fit is an option: with it, the term is deferred instead
        deferred = assign_bonded_interactions(molecule, classes, defer_missing=True)
        (angle,) = [a for a in deferred.angles if a.classes == ("S", "CT", "C")]
        self.assertEqual(angle.rule, RULE_QM_PENDING)
        self.assertLessEqual(set(deferred.rule_counts()), {RULE_QM_PENDING, "exact", "wildcard", RULE_IMPROPER})

    def test_default_classes_and_empirical_bond(self):
        molecule = molecule_from_smiles("CC(=O)Nc1ccccc1S(=O)(=O)N")
        names = [default_class(molecule, a.GetIdx()) for a in molecule.GetAtoms()
                 if a.GetAtomicNum() > 1]
        self.assertEqual(names, ["CT", "C", "O", "N", "CA", "CA", "CA", "CA", "CA", "CA",
                                 "SY", "O", "O", "N"])  # sulfonamide N is N, as in OPLS/2020
        rules = load_bonded_substitution(False).empirical
        b0, k = rules.bond(6, 6, 1.5)  # aromatic C-C: 1.526 - 0.71 log10(1.5) A
        self.assertAlmostEqual(b0, 0.1401, places=4)
        self.assertAlmostEqual(rules.bond_k(6, 6, b0), k)

    def test_missing_torsion_takes_default_general_term(self):
        # tolazamide: the azepane N-amino N is GAFF2-typed (n3); no table or
        # parmchk2 term about N-n3, so X-N-n3-X becomes X-N-NT-X (NT default class)
        molecule = molecule_from_smiles("Cc1ccc(cc1)S(=O)(=O)NC(=O)NN1CCCCCC1")
        types, parameters, _ = _assign_atom_types(molecule)
        classes = [parameters[t]["class"] for t in types]
        bonded = assign_bonded_interactions(molecule, classes)
        substituted = [t for t in bonded.proper_torsions if t.rule == RULE_DEFAULT_GENERAL]
        self.assertEqual(len(substituted), 4)
        for torsion in substituted:
            self.assertEqual({torsion.classes[1], torsion.classes[2]}, {"N", "n3"})
            self.assertEqual(set(torsion.parameter_classes[1:3]), {"N", "NT"})
            self.assertEqual(torsion.substitution_penalty, 1.0)
        deferred = assign_bonded_interactions(molecule, classes, defer_missing=True)
        self.assertNotIn(RULE_DEFAULT_GENERAL, deferred.rule_counts())

    def test_torsion_fallback_gaff2_dat_and_central_pair(self):
        def propers(smiles):
            molecule = molecule_from_smiles(smiles)
            types, parameters, _ = _assign_atom_types(molecule)
            classes = [parameters[t]["class"] for t in types]
            return assign_bonded_interactions(molecule, classes).proper_torsions

        # phorate: P(=S)-S about a GAFF2-typed p5 takes gaff2.dat X-p5-ss-X
        # (IDIVF 3, PK 11.4, n 2, 180: V2 = 2*11.4/3 kcal/mol per quartet)
        (term, *_) = [t for t in propers("CCOP(=S)(OCC)SCSCC") if t.rule == RULE_GAFF2_DAT_GENERAL]
        self.assertEqual(term.parameter_classes, ("X", "p5", "ss", "X"))
        self.assertAlmostEqual(term.coefficients_kj_mol[0], 2 * 11.4 / 3 * 4.184, places=6)
        # phenyltrimethylsilane: OPLS/2020 has only specific Si torsions
        hits = {(t.classes, t.parameter_classes, t.substitution_penalty)
                for t in propers("C[Si](C)(C)c1ccccc1") if t.rule == RULE_CENTRAL_PAIR}
        self.assertIn((("CA", "Si", "CT", "HC"), ("CT", "Si", "CT", "HC"), 1.0),
                      {(c if c[1] == "Si" else c[::-1], e, p) for c, e, p in hits})

    def test_last_resort_zero_torsion(self):
        # TMS-acetylene: Si-C(sp) with a GAFF2-typed alkyne C, which no rule covers
        molecule = molecule_from_smiles("C[Si](C)(C)C#C")
        types, parameters, _ = _assign_atom_types(molecule)
        classes = [parameters[t]["class"] for t in types]
        bonded = assign_bonded_interactions(molecule, classes)
        zeros = [t for t in bonded.proper_torsions if t.rule == RULE_ZERO_TORSION]
        self.assertEqual(len(zeros), 3)
        for torsion in zeros:
            self.assertEqual(torsion.coefficients_kj_mol, (0.0,) * 6)
            self.assertIsNone(torsion.substitution_penalty)
        self.assertTrue(any("set to zero" in w for w in bonded.warnings))
        deferred = assign_bonded_interactions(molecule, classes, defer_missing=True)
        self.assertNotIn(RULE_ZERO_TORSION, deferred.rule_counts())

    def test_untypable_gaff2_atom_only_affects_its_own_terms(self):
        # Si has no GAFF2 type; the GAFF2 alkyne C-C bond still comes from gaff2.dat
        molecule = molecule_from_smiles("C[Si](C)(C)C#C")
        types, parameters, _ = _assign_atom_types(molecule)
        classes = [parameters[t]["class"] for t in types]
        self.assertIsNone(load_gaff2_bonded().types(molecule)[1])
        bonded = assign_bonded_interactions(molecule, classes)
        rules = {tuple(b.atoms): b.rule for b in bonded.bonds}
        self.assertEqual(rules[(5, 6)], RULE_GAFF2_BONDED)  # C#C
        self.assertNotEqual(rules[(2, 5)], RULE_GAFF2_BONDED)  # Si-C(sp)

    def test_qm_mode_uses_only_the_opls2020_table(self):
        # phorate (GAFF2 p5, P-S torsions) and phenyltrimethylsilane (Si torsions)
        for smiles in ("CCOP(=S)(OCC)SCSCC", "C[Si](C)(C)c1ccccc1"):
            molecule = molecule_from_smiles(smiles)
            types, parameters, _ = _assign_atom_types(molecule)
            classes = [parameters[t]["class"] for t in types]
            default = assign_bonded_interactions(molecule, classes)
            deferred = assign_bonded_interactions(molecule, classes, defer_missing=True)
            with self.subTest(smiles=smiles):
                self.assertLessEqual(set(deferred.rule_counts()), {RULE_QM_PENDING, "exact", "wildcard", RULE_IMPROPER} | {RULE_GAFF2_IMPROPER})
                for kind in ("bonds", "angles", "proper_torsions"):
                    for a, b in zip(getattr(default, kind), getattr(deferred, kind)):
                        self.assertEqual(a.atoms, b.atoms)
                        if a.rule in ("exact", "wildcard"):
                            self.assertEqual(b.rule, a.rule)
                        else:
                            self.assertEqual(b.rule, RULE_QM_PENDING)

    def test_empirical_angle_without_a_bond_value_is_missing_not_keyerror(self):
        molecule = molecule_from_smiles("CCO")
        types, parameters, _ = _assign_atom_types(molecule)
        classes = tuple(parameters[t]["class"] for t in types[:3])
        fallback = load_bonded_substitution(False)
        self.assertIsNone(fallback.angle(molecule, (0, 1, 2), classes, {}))


class QmFailureTests(unittest.TestCase):
    def test_failed_fit_leaves_no_final_topology(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            stale = output / "tms_acetylene.itp"
            stale.write_text("; from an earlier run\n")
            molecule = molecule_from_smiles("C[Si](C)(C)C#C", name="tms_acetylene")
            with mock.patch("opls_ff_gen.qm_fit.fit_pending_terms",
                            side_effect=RuntimeError("xtb failed")):
                with self.assertRaisesRegex(RuntimeError, "xtb failed"):
                    generate_gromacs_topology(molecule, output, name="tms_acetylene",
                                              qm_method="xtb")
            for suffix in (".itp", ".gro", ".top", ".topology.json"):
                self.assertFalse((output / f"tms_acetylene{suffix}").exists(), suffix)
            pending = (output / "qm_pending" / "tms_acetylene.itp").read_text()
            self.assertIn(RULE_QM_PENDING, pending)


class QmBondedTableTests(unittest.TestCase):
    def test_angle_ring_context(self):
        carbonate = Chem.AddHs(Chem.MolFromSmiles("CC1COC(=O)O1"))  # O5 C4(=O6) O7
        self.assertEqual(angle_ring_context(carbonate, 3, 4, 6), "in5")
        self.assertEqual(angle_ring_context(carbonate, 3, 4, 5), "exo5")
        self.assertEqual(angle_ring_context(carbonate, 0, 1, 2), "exo5")
        chain = Chem.AddHs(Chem.MolFromSmiles("CCOC(=O)OCC"))
        self.assertEqual(angle_ring_context(chain, 2, 3, 5), "-")
        benzene = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
        self.assertEqual(angle_ring_context(benzene, 0, 1, 2), "-")

    def test_recorded_terms_fill_table_gaps_without_qm(self):
        caffeine = molecule_from_smiles("Cn1c(=O)c2c(ncn2C)n(C)c1=O", name="caffeine")
        types, parameters, _ = _assign_atom_types(caffeine)
        classes = [parameters[t]["class"] for t in types]
        # with the run-time QM fit, only the OPLS/2020 table comes before the fit
        strict = assign_bonded_interactions(caffeine, classes, defer_missing=True)
        self.assertNotIn(RULE_QM_TABLE, strict.rule_counts())
        pending = [v for v in strict.bonds + strict.angles if v.rule == RULE_QM_PENDING]
        self.assertTrue(pending)
        bonded = assign_bonded_interactions(caffeine, classes)
        recorded = load_qm_bonded()
        for value in bonded.bonds:
            if value.rule == RULE_QM_TABLE:
                term = recorded.bond(value.classes)
                self.assertEqual((value.length_nm, value.force_constant_kj_mol_nm2), term.values)
        self.assertEqual({tuple(v.atoms) for v in bonded.bonds + bonded.angles
                          if v.rule == RULE_QM_TABLE}, {tuple(v.atoms) for v in pending})


class Gaff2BondedLookupTests(unittest.TestCase):
    def test_gaff2_typed_atoms_take_gaff2_bonds_and_angles(self):
        # acetazolamide: the thiadiazole N=N atoms are GAFF2-typed (nb)
        molecule = molecule_from_smiles("CC(=O)Nc1nnc(s1)S(N)(=O)=O", name="acetazolamide")
        types, parameters, _ = _assign_atom_types(molecule)
        classes = [parameters[t]["class"] for t in types]
        self.assertTrue(any(c[:1].islower() for c in classes))
        bonded = assign_bonded_interactions(molecule, classes)
        gaff2 = [v for v in bonded.bonds + bonded.angles if v.rule == RULE_GAFF2_BONDED]
        self.assertTrue(gaff2)
        for value in gaff2:
            self.assertTrue(any(c[:1].islower() for c in value.classes), value)
        # with the run-time QM fit these terms are fitted instead
        deferred = assign_bonded_interactions(molecule, classes, defer_missing=True)
        self.assertNotIn(RULE_GAFF2_BONDED, deferred.rule_counts())
        self.assertEqual({tuple(v.atoms) for v in deferred.bonds + deferred.angles
                          if v.rule == RULE_QM_PENDING},
                         {tuple(v.atoms) for v in bonded.bonds + bonded.angles
                          if v.rule != "exact"})


class CacheTests(unittest.TestCase):
    def test_cache_is_reused_only_for_the_same_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            workdir = Path(directory) / "qm_fit"
            workdir.mkdir()
            (workdir / "opt_xtb.xyz").write_text("old\n")  # cache without a manifest
            ethanol = cache_inputs(molecule_from_smiles("CCO"), 0)
            warnings = prepare_cache(workdir, ethanol, "xtb", DEFAULT_DFT)
            self.assertEqual(len(warnings), 1)
            self.assertFalse((workdir / "opt_xtb.xyz").exists())
            self.assertEqual(len(list(Path(directory).glob("qm_fit.stale-*"))), 1)
            (workdir / "opt_xtb.xyz").write_text("new\n")
            self.assertEqual(prepare_cache(workdir, ethanol, "xtb", DEFAULT_DFT), [])
            self.assertTrue((workdir / "opt_xtb.xyz").exists())
            methanol = cache_inputs(molecule_from_smiles("CO"), 0)
            self.assertEqual(len(prepare_cache(workdir, methanol, "xtb", DEFAULT_DFT)), 1)
            self.assertEqual(len(list(Path(directory).glob("qm_fit.stale-*"))), 2)


    def test_every_term_records_its_rule(self):
        molecule = molecule_from_smiles("COc1ccccc1", name="anisole")  # fully tabulated
        types, parameters, _ = _assign_atom_types(molecule)
        bonded = assign_bonded_interactions(molecule, [parameters[t]["class"] for t in types])
        terms = bonded.bonds + bonded.angles + bonded.proper_torsions + bonded.improper_torsions
        self.assertTrue(all(value.rule for value in terms))
        self.assertEqual(sum(bonded.rule_counts().values()), len(terms))
        self.assertEqual({v.rule for v in bonded.improper_torsions}, {RULE_IMPROPER})
        self.assertTrue(all(v.source_line > 0 for v in bonded.bonds + bonded.angles))

    def test_no_torsion_spans_a_three_membered_ring(self):
        molecule = molecule_from_smiles("CC1CC1", name="methylcyclopropane")
        types, parameters, _ = _assign_atom_types(molecule)
        bonded = assign_bonded_interactions(molecule, [parameters[t]["class"] for t in types],
                                            defer_missing=True)
        for value in bonded.proper_torsions:
            outer1, center1, center2, outer2 = (index - 1 for index in value.atoms)
            self.assertIsNone(molecule.GetBondBetweenAtoms(outer1, center2))
            self.assertIsNone(molecule.GetBondBetweenAtoms(outer2, center1))


class SeminarioTests(unittest.TestCase):
    def test_diatomic_force_constant(self):
        k = 3.0e5
        coordinates = np.array([[0.0, 0.0, 0.0], [0.12, 0.0, 0.0]])
        hessian = np.zeros((6, 6))
        hessian[0, 0] = hessian[3, 3] = k
        hessian[0, 3] = hessian[3, 0] = -k
        r0, fitted = seminario_bond(hessian, coordinates, 0, 1)
        self.assertAlmostEqual(r0, 0.12)
        self.assertAlmostEqual(fitted, k)


@unittest.skipUnless(GMX.is_file() and XTB.is_file()
                     and importlib.util.find_spec("openmm") is not None,
                     "GROMACS, xtb or OpenMM unavailable")
class OnTheFlyFitTests(unittest.TestCase):
    def test_xtb_fit_writes_a_grompp_ready_topology(self):
        with tempfile.TemporaryDirectory() as directory:
            molecule = molecule_from_smiles("CN(C)C(=O)c1ccccc1", name="dmb")
            result = generate_gromacs_topology(molecule, directory, name="dmb",
                                               qm_method="xtb", qm_threads=4,
                                               torsion_substitution=False)
            fitted = [t for t in result.bonded.proper_torsions if t.rule == RULE_QM_FIT]
            self.assertEqual({min(t.classes, t.classes[::-1]) for t in fitted},
                             {("CA", "C", "N", "CT")})
            self.assertNotIn(RULE_QM_PENDING, result.bonded.rule_counts())
            self.assertFalse((Path(directory) / "qm_pending").exists())
            self.assertFalse(list(Path(directory).rglob("*.qmfit_*.*")))
            report = json.loads((Path(directory) / "qm_fit" / "report.json").read_text())
            self.assertEqual(len(report["scans"]), 1)
            scan = report["scans"][0]
            self.assertLess(scan["rmse_kj_mol"]["fitted"], scan["rmse_kj_mol"]["missing_zero"])
            self.assertNotIn("oplsaa.ff", result.top_path.read_text())
            self.assertIn("tertiary amide N uses OPLS type 239", result.itp_path.read_text())
            (Path(directory) / "check.mdp").write_text("integrator = md\nnsteps = 0\n")
            grompp = subprocess.run(
                [str(GMX), "grompp", "-f", "check.mdp", "-c", result.gro_path.name,
                 "-p", result.top_path.name, "-o", "check.tpr", "-maxwarn", "0"],
                cwd=directory, capture_output=True, text=True,
            )
            self.assertEqual(grompp.returncode, 0, grompp.stderr[-2000:])

    def test_torsions_about_linear_atoms_are_zero_not_fitted(self):
        with tempfile.TemporaryDirectory() as directory:
            molecule = molecule_from_smiles("N#CSC", name="mescn")  # methyl thiocyanate
            # Without substitution: parmchk2 would fill NZ-CZ-S-CT from X-CY-S-X.
            result = generate_gromacs_topology(molecule, directory, name="mescn",
                                               qm_method="xtb", qm_threads=4,
                                               torsion_substitution=False)
            linear = [t for t in result.bonded.proper_torsions if t.rule == RULE_ZERO_LINEAR]
            self.assertEqual([t.atoms for t in linear], [(1, 2, 3, 4)])
            self.assertEqual(linear[0].coefficients_kj_mol, (0.0,) * 6)
            self.assertTrue(any(RULE_ZERO_LINEAR in w for w in result.bonded.warnings))
            report = json.loads((Path(directory) / "qm_fit" / "report.json").read_text())
            self.assertEqual(report["zero_linear"], [[1, 2, 3, 4]])

    def test_gromacs_and_openmm_engines_agree(self):
        from opls_ff_gen.mm import mm_profile
        from opls_ff_gen.xtb import read_xyz

        with tempfile.TemporaryDirectory() as directory:
            molecule = molecule_from_smiles("CN(C)C(=O)c1ccccc1", name="dmb")  # one scan
            result = generate_gromacs_topology(molecule, directory, name="dmb",
                                               qm_method="xtb", qm_threads=4,
                                               torsion_substitution=False,
                                               mm_engine="gromacs")
            self.assertIn("MM relaxations: gromacs", result.bonded.qm_method)
            scan = json.loads(next((Path(directory) / "qm_fit").glob("scan_*/scan.json"))
                              .read_text())
            points = [{"phi_deg": p["phi_deg"],
                       "coordinates": read_xyz((Path(p["dir"]) / "xtbopt.xyz").read_text())}
                      for p in scan["points"][::4]]
            profiles = [np.asarray(mm_profile(result.top_path, points, tuple(scan["drive"]),
                                              molecule, lock_hydrogens=True, engine=engine))
                        for engine in ("openmm", "gromacs")]
            np.testing.assert_allclose(profiles[0] - profiles[0].min(),
                                       profiles[1] - profiles[1].min(), atol=0.05)


if __name__ == "__main__":
    unittest.main()


class DFTLevelTests(unittest.TestCase):
    def test_default_level_and_cache_tags(self):
        from opls_ff_gen.qm import DEFAULT_DFT, DFTLevel
        from opls_ff_gen.qm_fit import method_label

        self.assertEqual(DEFAULT_DFT.label, "B3LYP-D3BJ/def2-svp")
        self.assertIn("B3LYP-D3BJ/def2-svp", method_label("dft"))
        tz = DFTLevel("b3lyp", "def2-tzvp", "d3bj")
        self.assertNotEqual(tz.tag, DEFAULT_DFT.tag)  # levels never share cache files
        self.assertEqual(DFTLevel("wb97x-v", "def2-svp", None).label, "WB97X-V/def2-svp")
        self.assertRegex(tz.tag, r"^[a-z0-9-]+$")

    def test_unsupported_level_fails_fast(self):
        from opls_ff_gen.qm import DFTLevel

        with self.assertRaisesRegex(ValueError, "cannot be run with PySCF"):
            DFTLevel("wb97x-d", "def2-svp", None).validate()  # not implemented in PySCF
        with self.assertRaisesRegex(ValueError, "cannot be run with PySCF"):
            DFTLevel("b3lyp", "no-such-basis", "d3bj").validate()
