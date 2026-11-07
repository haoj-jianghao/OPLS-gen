"""Build ``src/opls_ff_gen/data/gaff2_bonded.json``: full GAFF2 types + bonds/angles.

The LJ rule table (``build_gaff2_table.py``, ``gaff2_lj.json``) gives one
representative type where GAFF2 types share Lennard-Jones parameters. Bonded
parameters need the full types, so this second rule table follows
antechamber's ``ATOMTYPE_GFF2.DEF`` type by type, in its order, with RDKit
SMARTS on the workflow's MMFF94 aromaticity. The conjugated pair types
(cc/cd, ce/cf, cg/ch, cp/cq, nc/nd, ne/nf, pc/pd, pe/pf) are written with the
first name; ``gaff2_bonded.py`` assigns the second by alternation along the
Kekulé bonds (single bond: same name, double bond: the other), as antechamber
does. Antechamber's ring classes are approximated: AR1 (pure aromatic) = a
6-ring of aromatic C and 2-connected aromatic N; AR2/AR3 (other planar rings)
= a 5- or 6-ring of sp2 C, N, 2-connected O and S. Checked against antechamber
with ``src/scripts/compare_gaff2_bonded_types.py``.

Also extracts the BOND and ANGL sections of gaff2.dat (K, r0 / theta0, as
written; ``gaff2_bonded.py`` converts to GROMACS units) and the general
X-b-c-X entries of its DIHE section (IDIVF, PK, phase, periodicity per term).

    python src/scripts/build_gaff2_bonded.py --gaff2-dat PATH/gaff2.dat

Without ``--gaff2-dat`` the parameter sections of the existing JSON are kept,
so rules can be edited without AmberTools.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from build_gaff2_table import AMIDE, CONJ, EW, H1, PL, exact_h, n_h  # noqa: E402

TABLE = Path(__file__).resolve().parents[1] / "opls_ff_gen" / "data" / "gaff2_bonded.json"

# AR1 (pure aromatic): the atom is in a 6-ring of aromatic C, 2-connected aromatic N
# and aromatic N+ (pyridinium)
_AR = "[c,$([n;X2]),$([n+;X3])]"
AR1 = f"$(*1:{_AR}:{_AR}:{_AR}:{_AR}:{_AR}:1)"
# AR2/AR3 (other planar rings): the atom is in a 5- to 7-ring of planar atoms
PLANAR = ",".join("$(*1" + f"~{PL}" * (size - 1) + "1)" for size in (5, 6, 7))
SB_DB = "$([a]),$(*=*)"  # has a single and a double bond (any 3-connected sp2 atom)
# single-bond partners that make a chain atom "conjugated" (ce, ne, pe, cg)
C3, C2 = "[#6;X3]", "[#6;X2]"
XB2 = "[#7,#15;X2]"
XD_DB = "[#16,#15;X3,X4;$(*=*)]"
XA1 = "[#8,#16;X1]"


def single_to(*partners: str) -> str:
    return ",".join(f"$(*-{p})" for p in partners)


RULES: list[tuple[str, str, str]] = [
    # hydrogen (as in the LJ table; hw (water) not used)
    ("hn", "[#1:1]~[#7]", ""),
    ("ho", "[#1:1]~[#8]", ""),
    ("hs", "[#1:1]~[#16]", ""),
    ("hp", "[#1:1]~[#15]", ""),
    ("hx", "[#1:1]~[#6]~[#7;X4]", ""),
    ("h3", f"[#1:1]~[#6;X4;$(*(~{EW})(~{EW})~{EW})]", ""),
    ("h2", f"[#1:1]~[#6;X4;$(*(~{EW})~{EW})]", ""),
    ("h1", f"[#1:1]~[#6;X4;$(*~{EW})]", ""),
    ("hc", "[#1:1]~[#6;X4]", ""),
    ("h5", f"[#1:1]~[#6;X3;$(*(~{EW})~{EW})]", ""),
    ("h4", f"[#1:1]~[#6;X3;$(*~{EW})]", ""),
    ("ha", "[#1:1]", ""),
    # carbon
    ("cx", "[#6;X4;r3:1]", "sp3 C in a 3-ring"),
    ("cy", "[#6;X4;r4:1]", "sp3 C in a 4-ring"),
    ("c5", "[#6;X4;r5:1]", "sp3 C in a 5-ring"),
    ("c6", "[#6;X4;r6:1]", "sp3 C in a 6-ring"),
    ("c3", "[#6;X4:1]", ""),
    ("cs", "[#6;X3:1]~[#16;X1]", "C=S"),
    ("c", "[#6;X3:1]~[#8,#16;X1]", "C=O"),
    ("cz", "[#6;X3:1](~[#7;X3])(~[#7;X3])~[#7;X3]", "guanidine C"),
    # cp: AR1 C in one 6-ring whose three neighbours are all AR1 atoms (biphenyl,
    # carbazole bridgeheads); gaff2_bonded.py demotes one in two 6-rings to ca
    ("cp", f"[#6;X3;{AR1};$(*(~[{AR1}])(~[{AR1}])~[{AR1}]):1]", "aromatic bridgehead"),
    ("ca", f"[#6;X3;{AR1}:1]", "pure aromatic C"),
    ("cc", f"[#6;X3;{PLANAR};{SB_DB}:1]", "sp2 C of conjugated ring systems (cc/cd)"),
    ("ce", f"[#6;X3;$(*=*);{single_to(C3, C2, XB2, XD_DB)}:1]",
     "sp2 C of conjugated chains (ce/cf)"),
    ("cu", "[#6;X3;r3:1]", ""),
    ("cv", "[#6;X3;r4:1]", ""),
    ("c2", "[#6;X3:1]", ""),
    ("cg", f"[#6;X2;$(*#*);{single_to(C2, C3, '[#7;X1]', XB2)}:1]",
     "sp C of conjugated systems (cg/ch)"),
    ("c1", "[#6;X2,X1:1]", ""),
    # nitrogen (ATOMTYPE_GFF2.DEF order)
    ("ns", f"[#7;X3;{exact_h(1)};{AMIDE}:1]", ""),
    ("nt", f"[#7;X3;{exact_h(2)};{AMIDE}:1]", ""),
    ("ni", f"[#7;X3;r3;{AMIDE}:1]", ""),
    ("nj", f"[#7;X3;r4;{AMIDE}:1]", ""),
    ("n", f"[#7;X3;{AMIDE}:1]", ""),
    ("nk", "[#7;X4;r3:1]", ""),
    ("nl", "[#7;X4;r4:1]", ""),
    ("nx", f"[#7;X4;{exact_h(1)}:1]", ""),
    ("ny", f"[#7;X4;{exact_h(2)}:1]", ""),
    ("nz", f"[#7;X4;{exact_h(3)}:1]", ""),
    ("n+", f"[#7;X4;{n_h(4)}:1]", ""),
    ("n4", "[#7;X4:1]", ""),
    ("no", "[#7;X3;$(*(~[#8;X1])~[#8;X1]):1]", ""),
    ("na", "[#7;X3;a:1]", ""),
    ("na", f"[#7;X3;$(*1~{PL}~{PL}~{PL}~{PL}1):1]", "planar non-aromatic 5-ring (AR3)"),
    ("na", f"[#7;X3;$(*1~{PL}~{PL}~{PL}~{PL}~{PL}1):1]", "planar non-aromatic 6-ring (AR3)"),
    ("nu", f"[#7;X3;{exact_h(1)};{CONJ}:1]", ""),
    ("nv", f"[#7;X3;{exact_h(2)};{CONJ}:1]", ""),
    ("nm", f"[#7;X3;r3;{CONJ}:1]", ""),
    ("nn", f"[#7;X3;r4;{CONJ}:1]", ""),
    ("nh", f"[#7;X3;{CONJ}:1]", ""),
    ("np", f"[#7;X3;r3;{exact_h(0)}:1]", ""),
    ("nq", f"[#7;X3;r4;{exact_h(0)}:1]", ""),
    ("n5", f"[#7;X3;r3;{exact_h(1)}:1]", ""),
    ("n6", f"[#7;X3;r4;{exact_h(1)}:1]", ""),
    ("n7", f"[#7;X3;{exact_h(1)}:1]", ""),
    ("n8", f"[#7;X3;{exact_h(2)}:1]", ""),
    ("n9", f"[#7;X3;{exact_h(3)}:1]", ""),
    ("n3", "[#7;X3:1]", ""),
    ("nb", f"[#7;X2;{AR1}:1]", "pure aromatic N"),
    ("nc", f"[#7;X2;{PLANAR};{SB_DB}:1]", "sp2 N of conjugated ring systems (nc/nd)"),
    ("ne", f"[#7;X2;$(*=*);{single_to(C2, C3, XA1, XB2, XD_DB)}:1]",
     "sp2 N of conjugated chains (ne/nf)"),
    ("n1", "[#7;X2;$(*#*),$(*(=*)=*):1]", ""),
    ("n2", "[#7;X2:1]", ""),
    ("n1", "[#7;X1:1]", ""),
    # oxygen
    ("o", "[#8;X1:1]", ""),
    ("oh", f"[#8;{H1}:1]", ""),
    ("op", "[#8;X2;r3:1]", ""),
    ("oq", "[#8;X2;r4:1]", ""),
    ("os", "[#8:1]", ""),
    # sulfur
    ("s", "[#16;X1:1]", ""),
    ("s2", "[#16;X2;$(*=*),$(*#*):1]", ""),
    ("sh", f"[#16;X2;{H1}:1]", ""),
    ("sp", "[#16;X2;r3:1]", ""),
    ("sq", "[#16;X2;r4:1]", ""),
    ("ss", "[#16;X2:1]", ""),
    ("sx", f"[#16;X3;$(*=*);{single_to(XB2, C3, XD_DB)}:1]", "conjugated s4"),
    ("s4", "[#16;X3:1]", ""),
    ("sy", f"[#16;X4;$(*=*);{single_to(XB2, C3, XD_DB)}:1]", "conjugated s6"),
    ("s6", "[#16:1]", ""),
    # phosphorus
    ("pb", f"[#15;X2;{AR1}:1]", ""),
    ("pc", f"[#15;X2;{PLANAR};{SB_DB}:1]", "sp2 P of conjugated ring systems (pc/pd)"),
    ("pe", f"[#15;X2;$(*=*);{single_to(C2, C3, XA1, XB2, XD_DB)}:1]",
     "sp2 P of conjugated chains (pe/pf)"),
    ("p2", "[#15;X1,X2:1]", ""),
    ("px", f"[#15;X3;$(*=*);{single_to(XB2, C3, XD_DB)}:1]", ""),
    ("p4", f"[#15;X3;$(*={XA1}):1]", ""),
    ("p3", "[#15;X3:1]", ""),
    ("py", f"[#15;X4;$(*=*);{single_to(XB2, C3, XD_DB)}:1]", ""),
    ("p5", "[#15:1]", ""),
    # halogens
    ("f", "[#9:1]", ""),
    ("cl", "[#17:1]", ""),
    ("br", "[#35:1]", ""),
    ("i", "[#53:1]", ""),
]

# second name of each conjugated pair type
PAIRS = {"cc": "cd", "ce": "cf", "cg": "ch", "cp": "cq",
         "nc": "nd", "ne": "nf", "pc": "pd", "pe": "pf"}

DESCRIPTION = (
    "Full GAFF2 atom types for bonded parameters (RDKit SMARTS, first match, mapped :1; "
    "conjugated pair types written with the first name, the second assigned by alternation "
    "along Kekule bonds), the BOND/ANGL sections of gaff2.dat (K in kcal/mol/A^2 or "
    "kcal/mol/rad^2, E = K (x - x0)^2) and its general X-b-c-X torsions ([IDIVF, PK, phase, "
    "periodicity] per term, E = PK/IDIVF (1 + cos(n phi - phase)); source line last). "
    "Generated by src/scripts/build_gaff2_bonded.py."
)


def read_sections(path: Path) -> tuple[str, dict, dict, dict]:
    lines = path.read_text().splitlines()
    end_mass = next(i for i, line in enumerate(lines) if i > 0 and not line.strip())
    i = end_mass + 2  # skip the hydrophilic-types line
    bonds, angles = {}, {}
    while lines[i].strip():
        name = lines[i][:5]
        k, r0 = lines[i][5:].split()[:2]
        bonds["-".join(t.strip() for t in name.split("-"))] = [float(k), float(r0), i + 1]
        i += 1
    i += 1
    while lines[i].strip():
        name = lines[i][:8]
        k, theta = lines[i][8:].split()[:2]
        angles["-".join(t.strip() for t in name.split("-"))] = [float(k), float(theta), i + 1]
        i += 1
    i += 1
    torsions: dict[str, list] = {}
    while lines[i].strip():
        quartet = [t.strip() for t in lines[i][:11].split("-")]
        if quartet[0] == quartet[3] == "X":  # general X-b-c-X entries only
            idivf, pk, phase, periodicity = lines[i][11:].split()[:4]
            entry = torsions.setdefault("-".join(quartet), [[], i + 1])
            entry[0].append([int(float(idivf)), float(pk), float(phase), float(periodicity)])
        i += 1
    return f"{path.name}: {lines[0].strip()}", bonds, angles, torsions


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gaff2-dat", type=Path)
    args = parser.parse_args()
    if args.gaff2_dat:
        source, bonds, angles, torsions = read_sections(args.gaff2_dat)
    else:
        existing = json.loads(TABLE.read_text())
        source, bonds, angles = existing["source"], existing["bonds"], existing["angles"]
        torsions = existing.get("torsions", {})
    rules = [{"type": label, "smarts": smarts, "note": note} for label, smarts, note in RULES]
    TABLE.write_text(json.dumps({"description": DESCRIPTION, "source": source, "pairs": PAIRS,
                                 "rules": rules, "bonds": bonds, "angles": angles,
                                 "torsions": torsions},
                                indent=1) + "\n")
    print(f"{TABLE}: {len(rules)} rules, {len(bonds)} bonds, {len(angles)} angles, "
          f"{len(torsions)} general torsions")


if __name__ == "__main__":
    main()
    from opls_ff_gen.datafiles import write_manifest  # keep data/MANIFEST.json current

    write_manifest()
