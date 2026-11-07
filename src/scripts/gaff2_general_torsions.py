"""GAFF2 general torsions (X-b-c-X) mapped to OPLS central class pairs.

Candidate X-B-C-X terms for OPLS/2020 from GAFF2, for review only (nothing
here feeds the workflow):

1. Every molecule of a GAFF2 reference set (``src/scripts/gaff2_reference.py``:
   SMILES and antechamber GAFF2 types in RDKit atom order) is also typed with
   the project's OPLS/2020 typer. Each bond that can carry torsions gives one
   observation: OPLS classes (B, C) <-> GAFF2 types (b, c) of the same atoms.
2. The general term GAFF2 uses about b-c is what tleap/parmchk2 would use
   when no specific quartet exists: X-b-c-X in gaff2.dat (parmchk2 step 3),
   else X-b'-c'-X with equal types (step 4, penalty 0), else with
   corresponding types (step 6, with parmchk2's penalty).
3. AMBER's E = PK/IDIVF [1 + cos(n phi - gamma)] per quartet is converted to
   OPLS Fourier V1-V4 per quartet (kcal/mol): V_n = +-2 PK/IDIVF, + when
   (n odd, gamma 0) or (n even, gamma 180), - otherwise; constants dropped.
4. Per OPLS pair (B, C): the converted terms of all GAFF2 pairs seen on it,
   weighted by bond count; ``unanimous`` if they agree, ``conflict`` if not,
   ``opls-general`` if OPLS/2020 already has X-B-C-X (kept as a check of the
   conversion), ``no-gaff`` if GAFF2 has no general term either.

    python src/scripts/gaff2_general_torsions.py --types REF.json --out DIR \\
        [--missing census.json] [--workers 8] [--par src/opls_ff_gen/data/gaff2_general_torsions.par]

Writes DIR/gaff2_general_torsions.json and DIR/summary.md; with ``--par``
also the supplementary torsion file the workflow reads
(``src/opls_ff_gen/gaff2_torsions.py``): one X-B-C-X line per OPLS pair without
an OPLS/2020 general term, the majority converted term by bond count (user
decision 2026-10-27). Terms that GAFF2 itself only has by type substitution
(parmchk2 step 6) and the OPLS/2020 general terms listed in ``QM_CHECK_OPLS``
are flagged ``qm-check``.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rdkit import Chem, RDLogger  # noqa: E402

from opls_ff_gen.assignment import _assign_atom_types  # noqa: E402
from opls_ff_gen.opls2020 import _TORSION_LINE, load_bonded_table  # noqa: E402
from opls_ff_gen.parmchk import ParmchkScores  # noqa: E402
from opls_ff_gen.structure import apply_typing_aromaticity  # noqa: E402

RDLogger.DisableLog("rdApp.*")
# OPLS/2020 general terms kept as they are but flagged for a QM check: the
# converted GAFF2 term differs strongly (aryl-N and aryl-carbonyl bonds).
QM_CHECK_OPLS = ("CB-N", "CB-NA", "CQ-N", "CA-NC", "C-CB")
# Central pairs whose GAFF2 term is replaced by the xTB fragment refit
# (validation/torsion_fallbacks/README.md, section 8; user decision 2026-10-28):
# V1-V3 (kcal/mol, per quartet) and the number of training fragments.
XTB_REFIT = {
    "CT-OS": ((-4.275, 1.027, -0.007), 19),
    "C-N": ((-0.076, 1.178, -2.762), 28),
    "CT-NT": ((-1.472, -0.594, 0.207), 9),
    "CT-OH": ((0.934, -0.730, 0.333), 18),
    "CW-SA": ((-1.352, 9.746, -4.053), 10),
    "CW-nb": ((-2.637, 5.688, -4.891), 5),
}
# Pairs without a GAFF2 term in the mapping, appended as xtb-refit lines
# (user decision 2026-10-29; benchmarks/borrowed_fit_2026_10_29, fit
# validation/torsion_fallbacks/data/refit_1029): phenyl-1,3,5-triazine CA-CQ;
# CA-CU on fused pyrazolo-pyrimidine rings, for ring bonds only (it is wrong on
# a rotatable aryl-pyrazole CA-CU bond). Value: V1-V3 (kcal/mol, per quartet),
# training fragments, ring bonds only.
XTB_REFIT_NEW = {
    "CA-CQ": ((0.133, 1.685, 0.183), 2, False),
    "CA-CU": ((-0.908, 14.059, -6.145), 7, True),
}
# Central pairs kept at their earlier term when a rebuild moves the majority
# (user decision 2026-10-28): V1-V4 (kcal/mol, per quartet), source, reason.
KEPT = {
    "CA-ca": ((0.0, 2.0, 0.0, 0.0), "X-cp-cp-X",
              "kept 2026-10-28: the new majority X-c2-ca-X is a parmchk2 substitution "
              "(penalty 232) on ring bonds"),
}
AMBERHOME = Path(os.environ.get("AMBERHOME", "/home/jiang/miniforge3/envs/chem"))
GAFF2_DAT = AMBERHOME / "dat" / "leap" / "parm" / "gaff2.dat"
_GAFF_TORSION = re.compile(
    r"^(.{2})-(.{2})-(.{2})-(.{2})\s+(\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)")


# -- GAFF2 general terms ------------------------------------------------------
class _Entry:
    def __init__(self, classes):
        self.classes = classes


class Gaff2Torsions:
    """gaff2.dat torsions; only the general X-b-c-X entries are offered to the
    parmchk2 steps, so steps 2 and 5 (specific entries) never match."""

    def __init__(self, path: Path = GAFF2_DAT):
        self.terms: dict[tuple[str, ...], list[tuple[int, float, float, float]]] = {}
        for line in path.read_text().splitlines():
            match = _GAFF_TORSION.match(line)
            if match:
                quartet = tuple(x.strip() for x in match.groups()[:4])
                self.terms.setdefault(quartet, []).append(
                    (int(match[5]), float(match[6]), float(match[7]), float(match[8])))

    def proper_entry(self, classes, general):
        if not general:
            return None
        for quartet in (tuple(classes), tuple(classes)[::-1]):
            if quartet in self.terms:
                return _Entry(quartet)
        return None

    def has_specific_central(self, b, c):
        return False

    def general(self, b: str, c: str, scores: ParmchkScores):
        """(GAFF2 entry, rule, penalty) used about b-c, or None."""
        direct = self.proper_entry(("X", b, c, "X"), general=True)
        if direct is not None:
            return direct.classes, "gaff2-general", 0.0
        hit = scores.substitute(self, ("X", b, c, "X"))
        if hit is None:
            return None
        return hit.parameter.classes, hit.rule, round(hit.penalty, 1)


def to_fourier(terms) -> tuple[float, float, float, float]:
    """AMBER (IDIVF, PK, phase, periodicity) lines -> OPLS V1-V4 per quartet."""
    v = [0.0] * 4
    for idivf, pk, phase, periodicity in terms:
        n = abs(int(periodicity))
        k = pk / idivf
        if k == 0.0:
            continue
        if n > 4 or phase not in (0.0, 180.0):
            raise ValueError(f"cannot convert n={n} phase={phase}")
        positive = (n % 2 == 1 and phase == 0.0) or (n % 2 == 0 and phase == 180.0)
        v[n - 1] += (2.0 if positive else -2.0) * k
    return tuple(round(x, 4) + 0.0 for x in v)


# -- observations -------------------------------------------------------------
def _observe(item: dict) -> dict:
    molecule = Chem.MolFromSmiles(item["smiles"])
    molecule = apply_typing_aromaticity(Chem.AddHs(molecule))
    if molecule.GetNumAtoms() != len(item["types"]):
        return {"skip": "atom count"}
    try:
        types, parameters, _ = _assign_atom_types(molecule)
    except Exception as exc:  # noqa: BLE001 - untypable molecule is a skip
        return {"skip": f"opls typing: {str(exc)[:60]}"}
    classes = [parameters[t]["class"] for t in types]
    pairs = []
    for bond in molecule.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if bond.GetBeginAtom().GetDegree() < 2 or bond.GetEndAtom().GetDegree() < 2:
            continue
        (b, gb), (c, gc) = sorted([(classes[i], item["types"][i]), (classes[j], item["types"][j])])
        pairs.append((b, c, gb, gc, bond.IsInRing()))
    return {"name": item["name"], "pairs": pairs}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--types", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--missing", type=Path,
                        help="census JSON; marks central pairs of the still-missing torsions")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--par", type=Path, help="write the supplementary torsion file")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    molecules = json.loads(args.types.read_text())["molecules"]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        observed = list(pool.map(_observe, molecules, chunksize=20))
    skipped = Counter(r["skip"].split(":")[0] for r in observed if "skip" in r)
    observed = [r for r in observed if "skip" not in r]

    bonds = defaultdict(Counter)  # (B, C) -> Counter((b, c))
    molecules_per_pair = defaultdict(set)
    ring = defaultdict(Counter)
    for record in observed:
        for b, c, gb, gc, in_ring in record["pairs"]:
            bonds[(b, c)][(gb, gc)] += 1
            molecules_per_pair[(b, c)].add(record["name"])
            ring[(b, c)]["ring" if in_ring else "chain"] += 1

    gaff = Gaff2Torsions()
    scores = ParmchkScores(uppercase_only=False)
    table = load_bonded_table()
    par_lines = table.path.read_text().splitlines()
    still_missing = set()
    if args.missing:
        from opls_ff_gen.parmchk import load_parmchk_scores
        upper = load_parmchk_scores()
        for record in json.loads(args.missing.read_text())["records"]:
            quartet = tuple(record["classes"].split("-"))
            if upper.substitute(table, quartet) is None:
                still_missing.add(tuple(sorted(quartet[1:3])))

    entries = []
    for (b, c), gaff_pairs in sorted(bonds.items()):
        opls = table.proper_entry(("X", b, c, "X"), general=True)
        opls_v = None
        if opls is not None:
            match = _TORSION_LINE.match(par_lines[opls.source_line - 1])
            opls_v = [float(match[i]) for i in (2, 3, 4, 5)]
        sources = []
        for (gb, gc), count in gaff_pairs.most_common():
            found = gaff.general(gb, gc, scores)
            if found is None:
                sources.append({"gaff2": f"{gb}-{gc}", "bonds": count, "term": None})
                continue
            quartet, rule, penalty = found
            sources.append({
                "gaff2": f"{gb}-{gc}", "bonds": count, "term": "-".join(quartet),
                "rule": rule, "penalty": penalty,
                "amber": gaff.terms[quartet], "fourier_kcal": to_fourier(gaff.terms[quartet]),
            })
        values = Counter()
        for source in sources:
            if source["term"] is not None:
                values[tuple(source["fourier_kcal"])] += source["bonds"]
        if opls is not None:
            status = "opls-general"
        elif not values:
            status = "no-gaff"
        elif len(values) == 1:
            status = "unanimous"
        else:
            status = "conflict"
        total = sum(gaff_pairs.values())
        majority = values.most_common(1)[0] if values else None
        entries.append({
            "central": f"{b}-{c}", "status": status, "bonds": total,
            "molecules": len(molecules_per_pair[(b, c)]), "ring_chain": dict(ring[(b, c)]),
            "opls_specific_entries": sum(
                1 for p in table._central_entries(b, c) if not p.classes[0] == p.classes[3] == "X"),
            "opls_general": None if opls is None else {
                "entry": "-".join(opls.classes), "fourier_kcal": opls_v,
                "source_line": opls.source_line},
            "majority_fourier_kcal": None if majority is None else list(majority[0]),
            "majority_share": None if majority is None else round(
                majority[1] / sum(values.values()), 3),
            "still_missing_in_100": (b, c) in still_missing,
            "gaff2_sources": sources,
        })

    report = {
        "types": str(args.types), "gaff2_dat": str(GAFF2_DAT),
        "molecules": len(observed), "skipped": dict(skipped),
        "status_counts": dict(Counter(e["status"] for e in entries)),
        "pairs": entries,
    }
    (args.out / "gaff2_general_torsions.json").write_text(json.dumps(report, indent=1) + "\n")
    _write_summary(args.out / "summary.md", report)
    if args.par:
        _write_par(args.par, report)
    print(f"{len(observed)} molecules, {len(entries)} OPLS central pairs: "
          f"{report['status_counts']}; skipped {dict(skipped)}")


def _write_par(path: Path, report: dict) -> None:
    """BOSS-style torsion lines (same columns as oplsaa2020.par)."""
    from datetime import date

    def pad(name: str) -> str:
        return f"{name:<2s}"

    lines = [
        "# GAFF2 general torsions (X-B-C-X) converted to OPLS Fourier form.",
        "# Supplementary to OPLS/2020 (not part of the ACS SI); read by",
        "# src/opls_ff_gen/gaff2_torsions.py after the OPLS/2020 lookup and before the",
        "# parmchk2 type substitution. Generated by src/scripts/gaff2_general_torsions.py",
        f"# on {date.today().isoformat()} from {report['molecules']} molecules (NCI first_5K +",
        "# 100 test molecules; antechamber -at gaff2 and the OPLS/2020 typer on the",
        "# same atoms) and gaff2.dat of AmberTools 26.",
        "# Per OPLS central pair B-C: the GAFF2 general term used about the GAFF2",
        "# types seen on the same bonds (X-b-c-X, else parmchk2 equal/corresponding",
        "# types), majority by bond count; V_n = +-2 PK/IDIVF per quartet (kcal/mol).",
        "# Columns: id V1 V2 V3 V4 X-B-C-X, then the GAFF2 source, the bond count",
        "# and majority share; 'qm-check' marks terms flagged for a QM check.",
        "# 'xtb-refit' lines replace the GAFF2 term with one V1-V3 fitted jointly",
        "# to relaxed GFN2-xTB scans of NCI fragments (2026-10-28; see",
        "# validation/torsion_fallbacks/README.md, section 8); the last lines add",
        "# pairs without a GAFF2 term (2026-10-29, section 11).",
        "# 'kept' lines keep an earlier term where a rebuild moved the majority.",
        "# OPLS/2020 general terms flagged for a QM check (kept as they are):",
        *(f"# qm-check-opls: X -{pad(b)}-{pad(c)}-X  aryl-N/aryl-carbonyl, GAFF2 differs"
          for b, c in (pair.split("-") for pair in QM_CHECK_OPLS)),
        "",
    ]
    number = 9000
    for e in sorted(report["pairs"], key=lambda e: (-e["bonds"], e["central"])):
        if e["status"] == "opls-general" or e["majority_fourier_kcal"] is None:
            continue
        values = tuple(e["majority_fourier_kcal"])
        sources = [s for s in e["gaff2_sources"]
                   if s["term"] is not None and tuple(s["fourier_kcal"]) == values]
        top = sources[0]
        # GAFF2 itself has no general term for the main source of this value
        flag = top["rule"] == "parmchk2-corr-general"
        b, c = e["central"].split("-")
        number += 1
        if e["central"] in XTB_REFIT:
            (v1, v2, v3), fragments = XTB_REFIT[e["central"]]
            lines.append(
                f"{number:4d} {v1:8.3f} {v2:8.3f} {v3:8.3f} {0.0:8.3f}   X -{pad(b)}-{pad(c)}-X "
                f"   xtb-refit, {fragments} fragments (was gaff2 {top['term']}"
                f" {' '.join(f'{v:g}' for v in values)}); {e['bonds']} bonds")
            continue
        if e["central"] in KEPT:
            (v1, v2, v3, v4), term, reason = KEPT[e["central"]]
            lines.append(
                f"{number:4d} {v1:8.3f} {v2:8.3f} {v3:8.3f} {v4:8.3f}   X -{pad(b)}-{pad(c)}-X "
                f"   gaff2 {term}; {reason} (majority {top['term']}"
                f" {' '.join(f'{v:g}' for v in values)}); {e['bonds']} bonds")
            continue
        v1, v2, v3, v4 = values
        lines.append(
            f"{number:4d} {v1:8.3f} {v2:8.3f} {v3:8.3f} {v4:8.3f}   X -{pad(b)}-{pad(c)}-X "
            f"   gaff2 {top['term']}"
            + ("" if top["rule"] == "gaff2-general" else f" ({top['rule']} {top['penalty']})")
            + f"; {e['bonds']} bonds, majority {e['majority_share']}"
            + ("; qm-check" if flag else ""))
    seen = {e["central"] for e in report["pairs"]}
    for central, ((v1, v2, v3), fragments, ring_only) in XTB_REFIT_NEW.items():
        if central in seen and not ring_only:
            raise SystemExit(f"{central} now has a GAFF2 term: move it to XTB_REFIT")
        b, c = central.split("-")
        number += 1
        lines.append(
            f"{number:4d} {v1:8.3f} {v2:8.3f} {v3:8.3f} {0.0:8.3f}   X -{pad(b)}-{pad(c)}-X "
            f"   xtb-refit, {fragments} fragments"
            + ("; ring bonds only" if ring_only else " (no GAFF2 term)"))
    path.write_text("\n".join(lines) + "\n")
    print(f"wrote {number - 9000} torsions to {path}")


def _fmt(v) -> str:
    return "—" if v is None else " ".join(f"{x:g}" for x in v)


def _write_summary(path: Path, report: dict) -> None:
    lines = [
        "# GAFF2 general torsions mapped to OPLS central pairs",
        "",
        f"{report['molecules']} molecules; status counts {report['status_counts']}.",
        "V1-V4 are OPLS Fourier coefficients per quartet (kcal/mol).",
        "",
        "| B-C | status | bonds | mols | OPLS specific | still missing | OPLS X-B-C-X V | GAFF2 (bonds: term, V) |",
        "|---|---|---:|---:|---:|---|---|---|",
    ]
    def source(s: dict) -> str:
        if s["term"] is None:
            return f"{s['gaff2']} ({s['bonds']}): none"
        rule = "" if s["rule"] == "gaff2-general" else f" [{s['rule']} {s['penalty']}]"
        return f"{s['gaff2']} ({s['bonds']}): {s['term']}{rule} {_fmt(s['fourier_kcal'])}"

    for e in sorted(report["pairs"], key=lambda e: -e["bonds"]):
        gaff = "; ".join(source(s) for s in e["gaff2_sources"][:6])
        if len(e["gaff2_sources"]) > 6:
            gaff += f"; +{len(e['gaff2_sources']) - 6} more"
        lines.append(
            f"| {e['central']} | {e['status']} | {e['bonds']} | {e['molecules']} | "
            f"{e['opls_specific_entries']} | {'yes' if e['still_missing_in_100'] else ''} | "
            f"{_fmt(e['opls_general']['fourier_kcal']) if e['opls_general'] else ''} | {gaff} |")
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
    from opls_ff_gen.datafiles import write_manifest  # keep data/MANIFEST.json current

    write_manifest()
