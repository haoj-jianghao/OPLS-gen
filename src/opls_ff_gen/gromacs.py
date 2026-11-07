"""Locating the GROMACS executable.

Found by ``tools.locate("gmx")``: ``OPLS_FF_GEN_GMX``, the user config file
(``opls-ff-gen config set gmx PATH``), else ``gmx`` (or ``gmx_mpi``) on ``PATH``.
Used by the GROMACS MM engine (``mm.py``), the batch checks and the tests.
"""

from __future__ import annotations

from pathlib import Path

from .tools import locate


def gmx_executable() -> Path | None:
    """The GROMACS executable, or None if none is found."""
    return locate("gmx").path


def require_gmx() -> Path:
    location = locate("gmx")
    if not location.ok:
        raise RuntimeError("GROMACS not found: run 'opls-ff-gen config set gmx /path/to/gmx', "
                           "set OPLS_FF_GEN_GMX, or put gmx on PATH (opls-ff-gen doctor)")
    return location.path
