"""Package data manifest (datafiles.py)."""

from __future__ import annotations

import json
import unittest

from opls_ff_gen import datafiles


class DataManifestTest(unittest.TestCase):
    def test_manifest_is_current(self):
        stored = json.loads(datafiles.MANIFEST.read_text())
        current = datafiles.build_manifest()
        self.assertEqual(stored["files"], current["files"],
                         "data/MANIFEST.json is stale: run python src/scripts/build_data_manifest.py")
        self.assertEqual(datafiles.verify(), ([], []))

    def test_si_is_never_listed(self):
        names = json.loads(datafiles.MANIFEST.read_text())["files"]
        self.assertFalse([name for name in names if name.startswith("opls2020/")])
        self.assertTrue({"opls_rules.json", "atomtype_archive.json",
                         "parmchk/PARMCHK.DAT"} <= set(names))


if __name__ == "__main__":
    unittest.main()
