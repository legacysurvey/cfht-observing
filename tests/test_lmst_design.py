"""Circular balancing and survey-progress checks."""

import unittest

import numpy as np
from astropy.table import Table

from update_lmst_design import design_values, uniform_lmst


class DesignTests(unittest.TestCase):
    def check_balanced(self, ra):
        result = uniform_lmst(ra)
        self.assertTrue(np.all((result >= 0) & (result < 360)))
        ha = (result - ra + 180) % 360 - 180
        self.assertAlmostEqual(float(np.mean(ha)), 0, places=9)
        ordered = np.sort(result)
        spacing = np.diff(np.r_[ordered, ordered[0] + 360])
        np.testing.assert_allclose(spacing, 360 / len(ra), atol=1e-9)
        return result

    def test_gapped_sky_wraparound_and_duplicate_ra(self):
        for ra in ([350, 359, 0, 1, 30, 60, 130, 180, 240, 300],
                   [359, 1], [0, 0, 0, 90, 90, 180, 270], [42]):
            with self.subTest(ra=ra):
                self.check_balanced(np.array(ra))

    def test_uniform_input_unchanged(self):
        ra = np.arange(0, 360, 10.0)
        np.testing.assert_allclose(self.check_balanced(ra), ra, atol=1e-10)

    def test_empty_and_nonfinite(self):
        self.assertEqual(len(uniform_lmst([])), 0)
        with self.assertRaises(ValueError):
            uniform_lmst([np.nan])

    def test_remaining_tiles_per_filter_and_no_program_cut(self):
        t = Table({'RA': [350., 10., 150., 250., 0., 90., 180.],
                   'FILTER': ['a', 'a', 'a', 'a', 'b', 'b', 'b'],
                   'IN_IBIS': [1, 1, 1, 0, 1, 1, 1],
                   'IN_HSC': [0, 0, 1, 0, 0, 0, 0],
                   'DONE': [0, 0, 0, 0, 0, 0, 1],
                   'PROGRAM': ['LBNL', 'NAOC', '', '', '', '', '']})
        result = design_values(t)
        self.assertTrue(np.all(np.isnan(result[[2, 3, 6]])))
        np.testing.assert_allclose(result[:2], uniform_lmst(t['RA'][:2]))
        np.testing.assert_allclose(result[4:6], uniform_lmst(t['RA'][4:6]))
        t['DONE'][0] = 1
        updated = design_values(t)
        self.assertTrue(np.isnan(updated[0]))
        self.assertAlmostEqual(updated[1], 10.)


if __name__ == '__main__':
    unittest.main()
