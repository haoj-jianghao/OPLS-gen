"""Build ``src/opls_ff_gen/data/gaff2_lj.json``: GAFF2 typing rules + LJ parameters.

The rules are written here with readable macros and expanded to plain RDKit
SMARTS in the JSON (which is what the package reads). They follow the order of
antechamber's ``ATOMTYPE_GFF2.DEF`` but only separate GAFF2 types whose
Lennard-Jones parameters differ; types with identical LJ get one
representative name (noted per rule). Checked against antechamber with
``src/scripts/gaff2_reference.py`` + ``src/scripts/compare_gaff2_rules.py``: 99.990%
LJ agreement on 145,971 atoms of RDKit's NCI first_5K set and 100% on a
targeted set of rare types (tests/data/gaff2_antechamber_reference.json).
Known antechamber differences (not reproduced):
- 3-connected N of cationic triarylmethane dyes (crystal violet type):
  antechamber re-perceives charges and bonds across the resonance system;
- cyclic N-nitrosamines: antechamber types the amine N n3 (acyclic: nh).

    python src/scripts/build_gaff2_table.py [--gaff2-dat PATH/gaff2.dat]

Without ``--gaff2-dat`` the LJ section of the existing JSON is kept, so rules
can be edited without AmberTools.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

TABLE = Path(__file__).resolve().parents[1] / "opls_ff_gen" / "data" / "gaff2_lj.json"

EW = "[#7,#8,#9,#16,#17,#35,#53]"  # electron-withdrawing neighbours of an H-bearing C
H1 = "$(*~[#1])"  # has at least one H


def n_h(n: int) -> str:
    return "$(*" + "(~[#1])" * (n - 1) + "~[#1])"


def exact_h(n: int) -> str:
    return f"{n_h(n)};!{n_h(n + 1)}" if n else f"!{H1}"


# N bonded to an aromatic atom or to C3/N2/P2 carrying a double bond
CONJ = "$(*-[a,$([#6;X3]=*),$([#7;X2]=*),$([#15;X2]=*)])"
# N bonded to a 3-connected C that carries a terminal O or S (amide, urea, ...)
AMIDE = "$(*~[#6;X3]~[#8,#16;X1])"
# planar ring members for antechamber's AR3 rings (no sp3 C)
PL = "[#6&X3,#7,#8&X2,#16&X2]"

RULES: list[tuple[str, str, str]] = [
    # hydrogen
    ("hn", "[#1:1]~[#7]", ""),
    ("ho", "[#1:1]~[#8]", ""),
    ("hs", "[#1:1]~[#16]", ""),
    ("hp", "[#1:1]~[#15]", ""),
    ("hx", "[#1:1]~[#6]~[#7;X4]", "H on C next to a 4-connected (cationic) N"),
    ("h3", f"[#1:1]~[#6;X4;$(*(~{EW})(~{EW})~{EW})]", ""),
    ("h2", f"[#1:1]~[#6;X4;$(*(~{EW})~{EW})]", ""),
    ("h1", f"[#1:1]~[#6;X4;$(*~{EW})]", ""),
    ("hc", "[#1:1]~[#6;X4]", ""),
    ("h5", f"[#1:1]~[#6;X3;$(*(~{EW})~{EW})]", ""),
    ("h4", f"[#1:1]~[#6;X3;$(*~{EW})]", ""),
    ("ha", "[#1:1]", "any other H (aromatic/sp2/sp C)"),
    # carbon: LJ depends only on connectivity
    ("c3", "[#6;X4:1]", "sp3 C (c3/c5/c6/cx/cy share LJ)"),
    ("c", "[#6;X3:1]~[#8,#16;X1]", "C=O/C=S (sp2 group)"),
    ("ca", "[#6;X3;a:1]", "aromatic C (sp2 group)"),
    ("c2", "[#6;X3:1]", "other sp2 C (sp2 group)"),
    ("c1", "[#6;X2,X1:1]", "sp C (c1/cg/ch share LJ)"),
    # nitrogen, in the order of antechamber's ATOMTYPE_GFF2.DEF
    ("nk", "[#7;X4;r3:1]", ""),
    ("nl", "[#7;X4;r4:1]", ""),
    ("nx", f"[#7;X4;{exact_h(1)}:1]", ""),
    ("ny", f"[#7;X4;{exact_h(2)}:1]", ""),
    ("nz", f"[#7;X4;{exact_h(3)}:1]", ""),
    ("n+", f"[#7;X4;{n_h(4)}:1]", ""),
    ("n4", "[#7;X4:1]", ""),
    ("ns", f"[#7;X3;{exact_h(1)};{AMIDE}:1]", "amide N-H"),
    ("nt", f"[#7;X3;{exact_h(2)};{AMIDE}:1]", "amide NH2"),
    ("ni", f"[#7;X3;r3;{AMIDE}:1]", ""),
    ("nj", f"[#7;X3;r4;{AMIDE}:1]", ""),
    ("n", f"[#7;X3;{AMIDE}:1]", "amide N"),
    ("no", "[#7;X3;$(*(~[#8;X1])~[#8;X1]):1]", "nitro N"),
    ("na", "[#7;X3;a:1]", "aromatic 3-connected N"),
    ("na", f"[#7;X3;$(*1~{PL}~{PL}~{PL}~{PL}1):1]", "N in a planar non-aromatic 5-ring (antechamber AR3)"),
    ("na", f"[#7;X3;$(*1~{PL}~{PL}~{PL}~{PL}~{PL}1):1]",
     "N in a planar non-aromatic 6-ring, e.g. phenothiazine (AR3)"),
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
    ("n1", "[#7;X2;$(*#*),$(*(=*)=*):1]", "sp N"),
    ("nb", "[#7;X2;a:1]", "aromatic sp2 N (nb/nc/nd share LJ with n2)"),
    ("n2", "[#7;X2:1]", "sp2 N (nb/nc/nd/ne/nf/n2 share LJ)"),
    ("n1", "[#7;X1:1]", "terminal N"),
    # oxygen, sulfur, phosphorus, halogens
    ("o", "[#8;X1:1]", ""),
    ("oh", f"[#8;{H1}:1]", ""),
    ("os", "[#8:1]", "os/op/oq share LJ"),
    ("s6", "[#16;X4,X5,X6:1]", "s6/sy share LJ"),
    ("s", "[#16;X1:1]", ""),
    ("sh", f"[#16;{H1}:1]", ""),
    ("ss", "[#16:1]", "s2/s4/sp/sq/ss/sx share LJ with s/sh"),
    ("p5", "[#15:1]", "all P types share LJ"),
    ("f", "[#9:1]", ""),
    ("cl", "[#17:1]", ""),
    ("br", "[#35:1]", ""),
    ("i", "[#53:1]", ""),
]

DESCRIPTION = (
    "GAFF2 Lennard-Jones parameters and an RDKit SMARTS rule table assigning GAFF2 types, used "
    "for atoms no OPLS type covers (replaces antechamber). Rules are tried in order; the first "
    "match types the atom (mapped :1). Where GAFF2 distinguishes types with identical LJ (e.g. "
    "ca/cc/cd/ce/cp, nb/nc/nd/ne/n2, os/op/oq, s/s2/sh/ss/s4/sx) the rules assign one "
    "representative type. Generated by src/scripts/build_gaff2_table.py; edit the rules there."
)


def read_gaff2_dat(path: Path) -> tuple[str, dict]:
    lines = path.read_text().splitlines()
    descriptions = {}
    for line in lines[1:120]:
        match = re.match(r"^(\S{1,2})\s+[0-9.]+\s+[0-9.]*\s+(.*)$", line)
        if match:
            descriptions.setdefault(match[1], match[2].strip())
    start = next(i for i, line in enumerate(lines) if line.startswith("MOD4"))
    lennard_jones = {}
    for number, line in enumerate(lines[start + 1:], start + 2):
        if not line.strip() or line.strip() == "END":
            break
        match = re.match(r"^\s*(\S+)\s+([0-9.]+)\s+([0-9.]+)", line)
        if match:
            lennard_jones[match[1]] = {
                "rmin_half_angstrom": float(match[2]),
                "epsilon_kcal_mol": float(match[3]),
                "line": number,
                "description": descriptions.get(match[1], ""),
            }
    source = (f"{path.name}: {lines[0].strip()}; J. Wang et al., J. Comput. Chem. 25, 1157 "
              "(2004) (GAFF); GAFF2 (AmberTools)")
    return source, lennard_jones


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gaff2-dat", type=Path, help="take LJ parameters from this gaff2.dat")
    args = parser.parse_args()
    if args.gaff2_dat:
        source, lennard_jones = read_gaff2_dat(args.gaff2_dat)
    else:
        existing = json.loads(TABLE.read_text())
        source, lennard_jones = existing["source"], existing["lennard_jones"]
    missing = sorted({label for label, _, _ in RULES} - lennard_jones.keys())
    if missing:
        raise SystemExit(f"rule types without LJ parameters: {missing}")
    rules = [{"type": label, "smarts": smarts, "note": note} for label, smarts, note in RULES]
    TABLE.write_text(json.dumps({"description": DESCRIPTION, "source": source,
                                 "lennard_jones": lennard_jones, "rules": rules}, indent=1) + "\n")
    print(f"{TABLE}: {len(lennard_jones)} LJ types, {len(rules)} rules")


if __name__ == "__main__":
    main()
    from opls_ff_gen.datafiles import write_manifest  # keep data/MANIFEST.json current

    write_manifest()
