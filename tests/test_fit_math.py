"""Torsion conversions and the ridge solve of the on-the-fly fit."""

from __future__ import annotations

import math
import unittest

import numpy as np

from opls_ff_gen.opls2020 import load_bonded_table
from opls_ff_gen.qm import fourier_to_rb
from opls_ff_gen.qm_fit import fourier_basis, solve


def fourier_energy(v, phi):
    return sum(vn / 2 * (1 + (-1) ** (n + 1) * math.cos(n * phi)) for n, vn in enumerate(v, 1))


def rb_energy(c, phi):
    psi = phi - math.pi
    return sum(cn * math.cos(psi) ** n for n, cn in enumerate(c))


class TorsionMathTests(unittest.TestCase):
    def test_rb_conversion_reproduces_the_fourier_energy(self):
        v = (5.4392, -0.2092, 0.8368, 0.3)
        c = fourier_to_rb(*v)
        for phi in np.linspace(-math.pi, math.pi, 13):
            self.assertAlmostEqual(rb_energy(c, phi), fourier_energy(v, phi), places=9)

    def test_rb_conversion_matches_the_opls2020_table(self):
        # OPLS/2020 CT-CT-CT-CT (type 003) read through the table parser
        table = load_bonded_table()
        expected = table.proper(("CT", "CT", "CT", "CT")).values
        for phi in np.linspace(-math.pi, math.pi, 13):
            self.assertAlmostEqual(rb_energy(expected, phi),
                                   rb_energy(fourier_to_rb(*_fourier_kj(table)), phi), places=9)

    def test_prior_keeps_unidentifiable_split(self):
        # Two carriers 180 degrees apart: only the sum of their V2 terms is
        # determined by the scan; the ridge correction keeps the prior
        # difference and fits the sum.
        phi = np.linspace(-180, 165, 24)
        design = np.array([np.concatenate([fourier_basis(p, 2), fourier_basis(p + 180, 2)])
                           for p in phi])
        truth = np.array([1.0, 6.0, -1.0, 4.0])
        prior = np.array([0.0, 3.0, 0.0, 1.0])
        r = design @ truth
        a, r = design - design.mean(axis=0), r - r.mean()
        x, determined = solve(a, r, prior, 1e-6)
        np.testing.assert_allclose(a @ x, r, atol=1e-6)          # scan reproduced
        self.assertAlmostEqual(x[1] + x[3], truth[1] + truth[3], places=5)  # V2 sum determined
        self.assertAlmostEqual(x[1] - x[3], prior[1] - prior[3], places=5)  # V2 split from prior
        np.testing.assert_allclose(determined[[1, 3]], 0.5, atol=1e-6)

    def test_ridge_keeps_flat_column_at_prior(self):
        # A torsion that barely moves (ring) has a nearly flat column; the
        # ridge leaves it near the prior instead of amplifying noise.
        rng = np.random.default_rng(0)
        a = np.column_stack([np.cos(np.radians(np.linspace(-180, 165, 24))),
                             1e-4 * rng.standard_normal(24)])
        a -= a.mean(axis=0)
        r = a @ np.array([5.0, 0.0]) + 0.1 * rng.standard_normal(24)
        x, determined = solve(a, r, np.array([0.0, 2.0]), 0.1)
        self.assertAlmostEqual(x[0], 5.0, delta=0.2)
        self.assertAlmostEqual(x[1], 2.0, delta=0.01)
        self.assertLess(determined[1], 1e-3)


def _fourier_kj(table):
    """V1-V4 (kJ/mol) of the CT-CT-CT-CT row, read from the .par text."""
    line = table.path.read_text().splitlines()[table.proper(("CT", "CT", "CT", "CT")).source_line - 1]
    return [4.184 * float(value) for value in line.split()[1:5]]


if __name__ == "__main__":
    unittest.main()
