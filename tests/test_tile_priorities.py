"""Stripe boundary, survey-year, and overlap checks for tile priorities."""

import unittest

import numpy as np
from astropy.table import Table

from update_tile_priorities import priority_values


def stripes(rows):
    return Table(rows=rows, names=['STRIPE', 'DEC_CENTER', 'RA_MIN', 'RA_MAX',
                                  'DEC_MIN', 'DEC_MAX', 'YEAR'])


class PriorityTests(unittest.TestCase):
    def test_years_and_all_ibis_rows(self):
        s = stripes([(f'NGC-{y}', y*10., 120., 250., y*10.-1.5, y*10.+1.5,
                      f'Y{y}') for y in range(1, 6)])
        t = Table({'RA': [180.] * 7, 'DEC': [10., 20., 30., 40., 50., 60., 10.],
                   'IN_IBIS': [1, 1, 1, 1, 1, 1, 0],
                   'IN_HSC': [1, 1, 0, 0, 0, 0, 0],
                   'DONE': [1, 0, 1, 0, 0, 0, 0],
                   'PROGRAM': ['NAOC', '', 'LBNL', 'NAOC', '', '', ''],
                   'PRIORITY': [1.] * 6 + [3.]})
        np.testing.assert_array_equal(priority_values(t, s), [10, 9, 8, 7, 6, 1, 3])

    def test_declination_padding_and_ra_boundaries(self):
        s = stripes([('NGC-7', -3.52, 129.81, 225.58, -5.375, -1.625, 'Y1')])
        t = Table({'RA': [180., 180., 180., 180., 129.81, 225.58, 129.809, 225.581],
                   'DEC': [-5.62, -1.42, -5.621, -1.419, -3.52, -3.52, -3.52, -3.52],
                   'IN_IBIS': [1] * 8})
        np.testing.assert_array_equal(priority_values(t, s), [10, 10, 1, 1, 10, 10, 1, 1])

    def test_highest_overlap_wins_independent_of_row_order(self):
        s = stripes([('NGC-1', 0., 100., 200., -1.5, 1.5, 'Y1'),
                     ('NGC-2', 2.8, 100., 200., 1.3, 4.3, 'Y2'),
                     ('NGC-3', 0., 100., 200., -1.5, 1.5, 'Y5')])
        t = Table({'RA': [150., 150.], 'DEC': [1.4, 3.], 'IN_IBIS': [1, 1]})
        np.testing.assert_array_equal(priority_values(t, s), [10, 9])
        np.testing.assert_array_equal(priority_values(t, s[::-1]), [10, 9])

    def test_sgc_wraparound(self):
        s = stripes([('SGC-7', -13.51, 306.53, 56.64, -15.375, -11.625, 'Y1')])
        t = Table({'RA': [0., 360., 359., 306.53, 56.64, 306.529, 56.641, 180.],
                   'DEC': [-13.51]*8, 'IN_IBIS': [1]*8})
        np.testing.assert_array_equal(priority_values(t, s), [10, 10, 10, 10, 10, 1, 1, 1])

    def test_padding_applies_to_later_years_by_default(self):
        s = stripes([('NGC-2', 0., 100., 200., -1.5, 1.5, 'Y2')])
        t = Table({'RA': [150.], 'DEC': [2.1], 'IN_IBIS': [1]})
        np.testing.assert_array_equal(priority_values(t, s), [9])
        np.testing.assert_array_equal(priority_values(t, s, pad_all=False), [1])

    def test_recomputing_is_idempotent(self):
        s = stripes([('NGC-1', 0., 100., 200., -1.5, 1.5, 'Y1')])
        t = Table({'RA': [150., 150.], 'DEC': [0., 10.], 'IN_IBIS': [1, 1]})
        t['PRIORITY'] = priority_values(t, s)
        np.testing.assert_array_equal(priority_values(t, s), t['PRIORITY'])


if __name__ == '__main__':
    unittest.main()
