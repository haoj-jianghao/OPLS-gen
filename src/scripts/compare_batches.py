"""Compare two batch output directories (run_batch_opls2020.py) file by file.

For every molecule present in both: the ITP, GRO and TOP files must be
byte-identical; the topology and assignment JSON must be identical after the
output directory is normalized; the QM-fit report must be identical apart
from timings. Prints one line per differing molecule and a summary.

    python src/scripts/compare_batches.py REFERENCE_DIR NEW_DIR [--allow NAME ...]
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

IGNORE: list[re.Pattern] = []


def _content(path: Path) -> bytes:
    """File bytes without lines matching an --ignore-line pattern."""
    data = path.read_bytes()
    if not IGNORE:
        return data
    lines = data.decode(errors="replace").splitlines(keepends=True)
    return "".join(l for l in lines if not any(r.search(l) for r in IGNORE)).encode()


def _normalized_json(path: Path, directory: Path, drop: tuple[str, ...] = ()):
    text = path.read_text().replace(str(directory.resolve()), "<DIR>")
    data = json.loads(text)
    for key in drop:
        data.pop(key, None)
    return data


def compare_molecule(reference: Path, new: Path) -> list[str]:
    problems = []
    for pattern in ("*.itp", "*.gro", "*.top"):
        ref_files = sorted(p.name for p in reference.glob(pattern) if ".qmfit_" not in p.name)
        new_files = sorted(p.name for p in new.glob(pattern) if ".qmfit_" not in p.name)
        if ref_files != new_files:
            problems.append(f"{pattern}: files {ref_files} vs {new_files}")
            continue
        for name in ref_files:
            if _content(reference / name) != _content(new / name):
                problems.append(f"{name} differs")
    for pattern in ("*.topology.json", "*.assignment.json"):
        for ref_path in reference.glob(pattern):
            new_path = new / ref_path.name
            if not new_path.exists():
                problems.append(f"{ref_path.name} missing")
            elif (_normalized_json(ref_path, reference)
                  != _normalized_json(new_path, new)):
                problems.append(f"{ref_path.name} differs")
    ref_report, new_report = reference / "qm_fit" / "report.json", new / "qm_fit" / "report.json"
    if ref_report.exists() != new_report.exists():
        problems.append("qm_fit/report.json present in only one")
    elif ref_report.exists():
        if (_normalized_json(ref_report, reference, ("timings_s",))
                != _normalized_json(new_report, new, ("timings_s",))):
            problems.append("qm_fit/report.json differs (excluding timings)")
    return problems


def first_difference(a: Path, b: Path) -> str:
    for number, (x, y) in enumerate(zip(a.read_text().splitlines(), b.read_text().splitlines()), 1):
        if x != y:
            if not any(r.search(x) for r in IGNORE):
                return f"line {number}:\n      - {x}\n      + {y}"
    return "length differs"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("new", type=Path)
    parser.add_argument("--ignore-line", nargs="*", default=[],
                        help="regex: lines of ITP/GRO/TOP files to ignore (e.g. a changed header comment)")
    parser.add_argument("--allow", nargs="*", default=[],
                        help="molecule directories with expected differences (reported, not counted)")
    args = parser.parse_args()
    IGNORE.extend(re.compile(pattern) for pattern in args.ignore_line)
    ref_dirs = {p.name for p in args.reference.iterdir() if p.is_dir()}
    new_dirs = {p.name for p in args.new.iterdir() if p.is_dir()}
    common = sorted(ref_dirs & new_dirs)
    print(f"molecules: reference {len(ref_dirs)}, new {len(new_dirs)}, compared {len(common)}")
    for name in sorted(ref_dirs ^ new_dirs):
        print(f"  only in {'reference' if name in ref_dirs else 'new'}: {name}")
    identical, unexpected = 0, 0
    for name in common:
        problems = compare_molecule(args.reference / name, args.new / name)
        if not problems:
            identical += 1
            continue
        tag = "expected" if name in args.allow else "DIFFERENT"
        unexpected += tag == "DIFFERENT"
        print(f"  {tag} {name}: {'; '.join(problems)}")
        for itp in (args.reference / name).glob("*.itp"):
            other = args.new / name / itp.name
            if ".qmfit_" not in itp.name and other.exists() and _content(itp) != _content(other):
                print(f"    first ITP difference, {first_difference(itp, other)}")
    print(f"identical {identical}/{len(common)}; unexpected differences {unexpected}")
    raise SystemExit(1 if unexpected or ref_dirs ^ new_dirs else 0)


if __name__ == "__main__":
    main()
