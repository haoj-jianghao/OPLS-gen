"""Export Foyer 1.2.0's OPLS-AA atom-typing data to ``src/opls_ff_gen/data/opls_rules.json``.

One-time migration (needs Foyer installed). Copies, in the XML's order (rule
order matters for Foyer's iterative typing), every atom type of Foyer's
``forcefields/xml/oplsaa.xml``: name, class, element, mass, SMARTS definition
(Foyer dialect), overrides and description. The XML's NonbondedForce
sigma/epsilon values (from GROMACS' oplsaa.ff, LGPL-2.1) are not exported
(2026-11-07): the workflow never used them (Lennard-Jones comes from
OPLS/2020 or GAFF2), and leaving them out keeps the file under Foyer's
MIT License alone (Copyright (c) 2015 Vanderbilt University).

    python src/scripts/export_foyer_rules.py
"""

from __future__ import annotations

import importlib.resources as resources
import json
from pathlib import Path
from xml.etree import ElementTree as ET

import foyer

OUT = Path(__file__).resolve().parents[1] / "opls_ff_gen" / "data" / "opls_rules.json"


def main() -> None:
    xml_path = resources.files("foyer.forcefields").joinpath("xml/oplsaa.xml")
    root = ET.parse(xml_path).getroot()
    types = []
    for atom_type in root.find("AtomTypes"):
        a = atom_type.attrib
        types.append({
            "name": a["name"],
            "class": a.get("class", ""),
            "element": a.get("element", ""),
            "mass": float(a["mass"]),
            "smarts": a.get("def", ""),
            "overrides": [o.strip() for o in a.get("overrides", "").split(",") if o.strip()],
            "description": a.get("desc", ""),
        })
    OUT.write_text(json.dumps({
        "description": (
            "OPLS-AA atom types with SMARTS typing rules (Foyer dialect), in Foyer's order. "
            "Rules are applied by src/opls_ff_gen/opls_typer.py; the atom-type archive adds to them."),
        "source": (f"Foyer {foyer.__version__} forcefields/xml/oplsaa.xml (MIT License, "
                   "Copyright (c) 2015 Vanderbilt University, mosdef-hub/foyer; C. Klein "
                   "et al., Comput. Phys. Commun. 2019, doi:10.1016/j.cpc.2018.10.023); "
                   "atom types, SMARTS rules and overrides only (no parameter values)"),
        "types": types,
    }, indent=1) + "\n")
    print(f"{OUT}: {len(types)} types, {sum(bool(t['smarts']) for t in types)} with SMARTS")


if __name__ == "__main__":
    main()
    from opls_ff_gen.datafiles import write_manifest  # keep data/MANIFEST.json current

    write_manifest()
