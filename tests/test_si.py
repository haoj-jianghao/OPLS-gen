"""User-installed OPLS/2020 SI files (si.py): lookup order, install-si checks."""

from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from opls_ff_gen import si

SOURCE = [si.PACKAGE_SI_DIR / f.name for f in si.SI_FILES]


@unittest.skipUnless(all(p.is_file() for p in SOURCE), "OPLS/2020 SI not installed here")
class SITest(unittest.TestCase):
    def setUp(self):
        scratch = tempfile.TemporaryDirectory()
        self.addCleanup(scratch.cleanup)
        self.dir = Path(scratch.name)
        clean = {k: v for k, v in os.environ.items() if k != si.SI_DIR_ENV}
        clean["XDG_DATA_HOME"] = str(self.dir / "xdg")
        patcher = mock.patch.dict(os.environ, clean, clear=True)
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_lookup_order(self):
        self.assertEqual(si.si_dir(), (si.PACKAGE_SI_DIR, "package directory"))
        report = si.install(list(reversed(SOURCE)))  # any order
        self.assertEqual(len(report), 2)
        self.assertTrue(all("verified" in line for line in report))
        self.assertEqual(si.si_dir(), (si.user_si_dir(), "user directory"))
        self.assertEqual(si.check(), ([], []))
        other = self.dir / "other"
        si.install(SOURCE, dest=other)
        os.environ[si.SI_DIR_ENV] = str(other)
        self.assertEqual(si.si_dir(), (other, f"env {si.SI_DIR_ENV}"))
        os.environ[si.SI_DIR_ENV] = str(self.dir / "empty")
        self.assertEqual(si.si_dir()[0], self.dir / "empty")
        errors, _ = si.check()
        self.assertIn("install-si", errors[0])

    def test_install_rejects_wrong_or_modified_files(self):
        junk = self.dir / "junk.txt"
        junk.write_text("not a parameter file\n")
        with self.assertRaisesRegex(si.SIInstallError, "not an OPLS/2020 SI"):
            si.install([junk, SOURCE[1]])
        with self.assertRaisesRegex(si.SIInstallError, "both files are needed"):
            si.install([SOURCE[1]])
        with self.assertRaisesRegex(si.SIInstallError, "are both"):
            si.install([SOURCE[0], SOURCE[0]])
        lines = SOURCE[0].read_text().splitlines(keepends=True)
        modified = self.dir / "modified.par"
        modified.write_text("".join(lines[:2] + [lines[2].replace("0.030", "0.031")] + lines[3:]))
        with self.assertRaisesRegex(si.SIInstallError, "--force"):
            si.install([modified, SOURCE[1]])
        report = si.install([modified, SOURCE[1]], force=True)
        self.assertIn("UNVERIFIED", report[0])
        errors, warnings = si.check()
        self.assertEqual(errors, [])
        self.assertIn("differ", warnings[0])


if __name__ == "__main__":
    unittest.main()
