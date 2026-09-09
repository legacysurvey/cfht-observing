"""Boundary checks for the HSC exclusion region."""

import unittest

import numpy as np
from astropy.table import Table

from update_tile_columns import hsc_mask, labels


def stripe(dec, ra_min, ra_max):
    return {'DEC': dec, 'RA_MIN': ra_min, 'RA_MAX': ra_max}


class TileColumnTests(unittest.TestCase):
    def test_isolated_stripe_inset(self):
        stripes = [stripe(2.08, 126.60, 249.93)]
        ra = [127.59, 127.60, 248.93, 248.94, 180, 180, 180, 180]
        dec = [2.08] * 4 + [1.679, 1.68, 2.48, 2.481]
        np.testing.assert_array_equal(hsc_mask(ra, dec, stripes),
                                      [0, 1, 1, 0, 0, 1, 1, 0])

    def test_joined_stripes_and_exposed_step(self):
        stripes = [stripe(2.08, 126.60, 249.93),
                   stripe(-0.72, 128.16, 228.30)]
        # The shared edge at DEC=0.68 has no inset. At RA=240 the lower
        # stripe is absent, so the upper stripe's southern inset remains.
        ra = [180] * 5 + [240, 240, 240]
        dec = [-1.121, -1.12, 0.68, 2.48, 2.481, 0.68, 1.679, 1.68]
        np.testing.assert_array_equal(hsc_mask(ra, dec, stripes),
                                      [0, 1, 1, 1, 0, 0, 0, 1])

    def test_wraparound_and_independent_hemispheres(self):
        stripes = [stripe(3.29, 315.57, 44.37),
                   stripe(0.49, 313.82, 47.15),
                   stripe(2.08, 126.60, 249.93),
                   stripe(-0.72, 128.16, 228.30)]
        ra = [0, 360, 359.9, 180, 316.56, 316.57, 43.37, 43.38, 0, 0]
        dec = [1.89, 1.89, 1.89, 0.68, 3.29, 3.29, 3.29, 3.29,
               3.69, 3.691]
        np.testing.assert_array_equal(hsc_mask(ra, dec, stripes),
                                      [1, 1, 1, 1, 0, 1, 1, 0, 1, 0])

    def test_gap_between_stripes(self):
        stripes = [stripe(0, 100, 200), stripe(4, 100, 200)]
        np.testing.assert_array_equal(hsc_mask([150] * 3, [0, 2, 4], stripes),
                                      [1, 0, 1])

    def test_program_plan_precedence(self):
        table = Table({'OBJECT': ['a', 'b', 'c', 'd'],
                       'IN_IBIS': [1, 1, 0, 0],
                       'RA': [180] * 4, 'DEC': [0] * 4})
        program, _ = labels(table, {'a', 'c'}, [])
        np.testing.assert_array_equal(program, ['NAOC', 'LBNL', 'NAOC', ''])

    def test_temporary_northern_extension(self):
        stripes = [dict(stripe(2.08, 126.60, 249.93), FIELD='NGC-5'),
                   dict(stripe(-0.72, 128.16, 228.30), FIELD='NGC-6'),
                   dict(stripe(3.29, 315.57, 44.37), FIELD='SGC-1'),
                   dict(stripe(0.49, 313.82, 47.15), FIELD='SGC-2')]
        ra = [180, 180, 127.59, 0, 360, 359, 43.38, 180, 0]
        dec = [4, 4.001, 3, 6.121, 6.121, 6.122, 5, -1.121, 0.089]
        original = hsc_mask(ra, dec, stripes)
        extended = hsc_mask(ra, dec, stripes,
                            {'NGC-5': 4.0, 'SGC-1': 6.121})
        np.testing.assert_array_equal(extended, [1, 0, 0, 1, 1, 0, 0, 0, 0])
        self.assertTrue(np.all(extended >= original))
        np.testing.assert_array_equal(hsc_mask(ra, dec, stripes), original)

    def test_sgc_limit_uses_only_sgc_ibis_tiles(self):
        stripes = [dict(stripe(2.08, 126.60, 249.93), FIELD='NGC-5'),
                   dict(stripe(3.29, 315.57, 44.37), FIELD='SGC-1')]
        table = Table({'OBJECT': ['a', 'b', 'c', 'd', 'e'],
                       'IN_IBIS': [1, 0, 1, 0, 0],
                       'RA': [0, 0, 180, 0, 180],
                       'DEC': [6.121, 8.0, 16.0, 6.0, 4.0]})
        _, mask = labels(table, set(), stripes, extend_hsc_north=True)
        np.testing.assert_array_equal(mask, [1, 0, 0, 1, 1])


if __name__ == '__main__':
    unittest.main()
