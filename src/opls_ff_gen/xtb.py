"""Running GFN2-xTB (xtb 6.x) and reading its geometries.

The executable is found by ``tools.locate("xtb")``: ``OPLS_FF_GEN_XTB``, the
user config file (``opls-ff-gen config set xtb PATH``), else ``xtb`` on
``PATH``; ``XTBPATH`` (parameter files) from ``OPLS_FF_GEN_XTBPATH``, the config
key ``xtbpath``, or ``share/xtb`` next to the executable.
"""

from __future__ import annotations

import os
import re
import subprocess
from functools import lru_cache
from pathlib import Path

import numpy as np

from .tools import program_path, setting

XTB = program_path("xtb")
XTBPATH = Path(setting("xtbpath") or XTB.resolve().parent.parent / "share" / "xtb")


def xtb_available() -> bool:
    """True if the xtb executable exists and can be run."""
    return XTB.is_file() and os.access(XTB, os.X_OK)


@lru_cache(maxsize=1)
def xtb_version() -> str:
    """Version line of the xtb executable ("6.7.1 (edcfbbe)"), "unknown" if it cannot run."""
    try:
        result = subprocess.run([str(XTB), "--version"], capture_output=True, text=True,
                                check=False, env=dict(os.environ, XTBPATH=str(XTBPATH)))
    except OSError:
        return "unknown"
    match = re.search(r"xtb version (\S+(?: \([0-9a-f]+\))?)", result.stdout)
    return match.group(1) if match else "unknown"


def run_xtb(workdir: Path, xyz: str, args: list[str], constrain: str | None = None,
            threads: int = 1) -> str:
    """Run ``xtb in.xyz --gfn 2 ARGS`` in ``workdir``; return stdout (raises on failure)."""
    workdir.mkdir(parents=True, exist_ok=True)
    (workdir / "in.xyz").write_text(xyz)
    command = [str(XTB), "in.xyz", "--gfn", "2", *args]
    if constrain is not None:
        (workdir / "constrain.inp").write_text(constrain)
        command += ["--input", "constrain.inp"]
    env = dict(os.environ, XTBPATH=str(XTBPATH), OMP_NUM_THREADS=str(threads),
               OMP_STACKSIZE="1G", MKL_NUM_THREADS=str(threads))
    result = subprocess.run(command, cwd=workdir, env=env, capture_output=True,
                            text=True, check=False)
    output = result.stdout + result.stderr
    (workdir / "xtb.log").write_text(output)
    if "normal termination" not in output:
        raise RuntimeError(f"xtb failed in {workdir}")
    return result.stdout


def read_xyz(text: str) -> np.ndarray:
    lines = text.strip().splitlines()
    count = int(lines[0])
    return np.array([[float(v) for v in line.split()[1:4]] for line in lines[2:2 + count]])


def to_xyz(elements: list[str], coordinates: np.ndarray) -> str:
    body = "\n".join(f"{e} {x:.8f} {y:.8f} {z:.8f}" for e, (x, y, z) in zip(elements, coordinates))
    return f"{len(elements)}\n\n{body}\n"
