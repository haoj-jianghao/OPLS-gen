"""External-program lookup, config file and pre-flight checks (tools.py)."""

from __future__ import annotations

import os
import stat
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from opls_ff_gen import tools


def _fake_program(directory: Path, name: str) -> Path:
    path = directory / name
    path.write_text("#!/bin/sh\necho 'xtb version 6.7.1 (fake)'\n")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


class ToolsTest(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory()
        self.dir = Path(self.scratch.name)
        clean = {k: v for k, v in os.environ.items()
                 if not k.startswith("OPLS_FF_GEN_")}
        clean.update(OPLS_FF_GEN_CONFIG=str(self.dir / "config.toml"), PATH=str(self.dir / "bin"))
        (self.dir / "bin").mkdir()
        patcher = mock.patch.dict(os.environ, clean, clear=True)
        patcher.start()
        self.addCleanup(patcher.stop)
        self.addCleanup(self.scratch.cleanup)

    def test_lookup_order_env_config_path(self):
        self.assertFalse(tools.locate("xtb").ok)
        on_path = _fake_program(self.dir / "bin", "xtb")
        self.assertEqual(tools.locate("xtb"), tools.Location(on_path, "PATH"))
        (self.dir / "opt").mkdir()
        configured = _fake_program(self.dir / "opt", "xtb")
        tools.set_config("xtb", str(configured))
        self.assertEqual(tools.locate("xtb").path, configured.resolve())
        self.assertTrue(tools.locate("xtb").source.startswith("config"))
        os.environ["OPLS_FF_GEN_XTB"] = str(on_path)
        self.assertEqual(tools.locate("xtb").source, "env OPLS_FF_GEN_XTB")
        self.assertEqual(tools.program_version("xtb"), "6.7.1")
        tools.set_config("xtb", None)
        self.assertNotIn("xtb", tools.read_config())

    def test_config_rejects_bad_entries(self):
        with self.assertRaises(ValueError):
            tools.set_config("xtb", str(self.dir / "missing"))
        with self.assertRaises(ValueError):
            tools.set_config("python", "/usr/bin/python3")

    def test_preflight(self):
        # the SI files are user-installed and may be absent: leave them out of the counts
        real = tools.preflight

        def preflight(**options):
            errors, warnings = real(**options)
            return [e for e in errors if not e.startswith("OPLS/2020 SI")], warnings

        errors, warnings = preflight(qm=None, mm_engine="openmm", optimize_xtb=True)
        self.assertEqual(len(errors), 1)
        self.assertIn("MOPAC", errors[0])
        self.assertIn("config set mopac", errors[0])
        self.assertIn("MMFF94", warnings[0])
        errors, _ = preflight(qm="xtb", mm_engine="gromacs", optimize_xtb=False)
        self.assertEqual([e.split(":")[0] for e in errors], ["MOPAC", "xtb (GFN2-xTB)", "GROMACS"])
        with mock.patch.object(tools, "package_available", return_value=False):
            errors, _ = preflight(qm="dft", mm_engine="openmm", optimize_xtb=False)
        self.assertEqual([e.split(":")[0] for e in errors],
                         ["MOPAC", "xtb (GFN2-xTB)", "OpenMM", "PySCF", "geomeTRIC"])

    def test_default_mm_engine_uses_what_is_available(self):
        with mock.patch.object(tools, "package_available", return_value=False):
            self.assertEqual(tools.default_mm_engine(), "openmm")  # nothing: reported later
            _fake_program(self.dir / "bin", "gmx")
            self.assertEqual(tools.default_mm_engine(), "gromacs")
        os.environ["OPLS_FF_GEN_MM_ENGINE"] = "openmm"
        self.assertEqual(tools.default_mm_engine(), "openmm")



class VersionNoteTest(unittest.TestCase):
    def test_notes_only_outside_the_tested_versions(self):
        for key, version in (("mopac", "23.2.5"), ("mopac", "23.1.0"), ("xtb", "6.7.1"),
                             ("rdkit", "2025.09.6"), ("rdkit", "2026.09.2"), ("gmx", "2026.3"),
                             ("xtb", None)):
            self.assertIsNone(tools.version_note(key, version), (key, version))
        for key, version in (("mopac", "22.1.1"), ("xtb", "6.6.1"), ("rdkit", "2024.03.5"),
                             ("rdkit", "2027.03.1"), ("xtb", "unknown")):
            note = tools.version_note(key, version)
            self.assertIsNotNone(note, (key, version))
            self.assertNotIn("error", note.lower())


class SourceLauncherTest(unittest.TestCase):
    def test_run_from_source_tree(self):
        import subprocess
        import sys

        launcher = Path(__file__).resolve().parents[1] / "run_opls_ff_gen.py"
        result = subprocess.run([sys.executable, str(launcher), "config", "--help"],
                                capture_output=True, text=True, timeout=120)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("opls-ff-gen config set KEY PATH", result.stdout)


if __name__ == "__main__":
    unittest.main()
