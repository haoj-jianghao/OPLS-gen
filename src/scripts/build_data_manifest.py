"""Write ``src/opls_ff_gen/data/MANIFEST.json`` (size and SHA-256 of every
shipped data file; ``opls_ff_gen.datafiles``). Run after editing a data file;
the scripts that rebuild data files call it themselves.

    python src/scripts/build_data_manifest.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from opls_ff_gen.datafiles import build_manifest, write_manifest  # noqa: E402


def main() -> None:
    path = write_manifest()
    print(f"{path}: {len(build_manifest()['files'])} files")


if __name__ == "__main__":
    main()
