"""Check built release files (wheel and sdist) before they are published.

    python src/scripts/check_dist.py dist/

- the wheel contains every file listed in data/MANIFEST.json with the listed
  SHA-256, the manifest itself, LICENSE, NOTICE and the AmberTools GPL text;
- neither the wheel nor the sdist contains the OPLS/2020 SI (data/opls2020/,
  jp3c06602_si_*, oplsaa2020.*) or local-only material (benchmarks/, skills/,
  validation/, src/scripts/);
- the wheel's version matches pyproject.toml (and the tag, with --tag).
Exit code 1 on any failure.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import tarfile
import tomllib
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FORBIDDEN = re.compile(r"(^|/)(opls2020/|jp3c06602_si_|oplsaa2020\.|benchmarks/|skills/|"
                       r"validation/|src/scripts/|\.claude/)")


def check(dist: Path, tag: str | None) -> list[str]:
    problems = []
    version = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["version"]
    if tag is not None and tag.lstrip("v") != version:
        problems.append(f"tag {tag} does not match pyproject version {version}")
    wheels = sorted(dist.glob("opls_ff_gen-*.whl"))
    sdists = sorted(dist.glob("opls_ff_gen-*.tar.gz"))
    if len(wheels) != 1 or len(sdists) != 1:
        return problems + [f"expected one wheel and one sdist in {dist}, found "
                           f"{[w.name for w in wheels]} {[s.name for s in sdists]}"]
    wheel, sdist = wheels[0], sdists[0]
    if f"-{version}-" not in wheel.name:
        problems.append(f"{wheel.name} is not version {version}")

    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        problems += [f"wheel contains {n}" for n in names if FORBIDDEN.search(n)]
        manifest_name = "opls_ff_gen/data/MANIFEST.json"
        if manifest_name not in names:
            return problems + ["wheel has no data/MANIFEST.json"]
        listed = json.loads(archive.read(manifest_name))["files"]
        for name, entry in listed.items():
            path = f"opls_ff_gen/data/{name}"
            if path not in names:
                problems.append(f"wheel lacks package data {name}")
            elif hashlib.sha256(archive.read(path)).hexdigest() != entry["sha256"]:
                problems.append(f"wheel data {name} differs from MANIFEST.json")
        for license_name in ("LICENSE", "NOTICE", "AmberTools-LICENSE.txt"):
            if not any(n.endswith("/" + license_name) for n in names):
                problems.append(f"wheel lacks {license_name}")
        metadata = archive.read(next(n for n in names if n.endswith(".dist-info/METADATA")))
        if b"License-Expression: MIT AND GPL-3.0-only" not in metadata:
            problems.append("wheel metadata lacks License-Expression: MIT AND GPL-3.0-only")
        print(f"{wheel.name}: {len(names)} files, {len(listed)} data files checked")

    with tarfile.open(sdist) as archive:
        names = [m.name.split("/", 1)[1] for m in archive.getmembers() if "/" in m.name]
        problems += [f"sdist contains {n}" for n in names if FORBIDDEN.search(n)]
        for required in ("pyproject.toml", "LICENSE", "NOTICE", "INSTALL.md",
                         "src/opls_ff_gen/data/MANIFEST.json"):
            if required not in names:
                problems.append(f"sdist lacks {required}")
        print(f"{sdist.name}: {len(names)} files")
    return problems


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("dist", type=Path)
    parser.add_argument("--tag", help="release tag, e.g. v0.3.0")
    args = parser.parse_args()
    problems = check(args.dist, args.tag)
    for problem in problems:
        print(f"FAIL: {problem}")
    print("dist check: " + ("FAILED" if problems else "ok"))
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
