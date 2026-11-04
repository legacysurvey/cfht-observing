"""Program membership from the imaging-plan footprint, and PROGRAM labels."""

import unittest

import numpy as np
from astropy.table import Table

from footprint import DECAM_INSIDE, DECAM_OUTER, MEGACAM_INSIDE, MEGACAM_OUTER, Footprint
from update_tile_columns import hsc_flags, labels


class FootprintTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # One POLU stripe (RA 125-250, DEC 0-4) with a Niji rectangle filling RA 130-180,
        # DEC 1-3, and a dropped sliver at RA 125-130, so POLU sky starts at RA 130.
        cls.fp = Footprint(polu_rects=[(125.0, 250.0, 0.0, 4.0)],
                           excluded_rects=[(130.0, 180.0, 1.0, 3.0), (125.0, 130.0, 0.0, 4.0)],
                           niji_rects=[(130.0, 180.0, 1.0, 3.0)])

    def program(self, ra, dec, outer=MEGACAM_OUTER, inside=MEGACAM_INSIDE):
        return self.fp.program_tiles(ra, dec, outer, inside).astype(int).tolist()

    def test_outer_margin_at_stripe_edges(self):
        # The field (half-size 0.5) must overlap the stripe (DEC 0-4): within half a
        # degree of the edge counts; beyond it does not.
        self.assertEqual(self.program([200] * 4, [-0.4, -0.7, 4.4, 4.7]), [1, 0, 1, 0])
        # At RA 130 the POLU sky at DEC 2 is excluded (Niji), so test the RA edge at DEC 0.5.
        self.assertEqual(self.program([129.6, 129.3, 250.4, 250.7], [0.5, 0.5, 2, 2]), [1, 0, 1, 0])

    def test_niji_interior_excluded_but_edges_kept(self):
        # A 1x1 deg field entirely inside the Niji rectangle is excluded; one whose
        # edge reaches POLU sky (within 0.5 deg of the rectangle edge) is kept.
        self.assertEqual(self.program([155] * 4, [2.0, 1.4, 1.6, 2.6]), [0, 1, 0, 1])
        # At the east edge of the rectangle (RA 180) POLU sky continues at the same DEC.
        self.assertEqual(self.program([179.6, 179.4], [2, 2]), [1, 0])
        # At the west edge only the dropped sliver lies beyond, so nothing is kept.
        self.assertEqual(self.program([130.4, 130.6], [2, 2]), [0, 0])

    def test_dropped_sliver_is_not_program_sky(self):
        # Deep in the dropped sliver nothing is near; within half a field of RA 130 the field overlaps.
        self.assertEqual(self.program([126.5, 129.6], [0.5, 0.5]), [0, 1])

    def test_polu_ibis_zero_stripe_contributes_nothing(self):
        fp = Footprint(polu_rects=[], excluded_rects=[], niji_rects=[(0.0, 10.0, 5.0, 8.0)])
        self.assertEqual(fp.program_tiles([5, 5], [6.5, 4.0], MEGACAM_OUTER, MEGACAM_INSIDE).astype(int).tolist(), [0, 0])

    def test_ra_wraparound(self):
        fp = Footprint(polu_rects=[(350.0, 10.0, 0.0, 4.0)], excluded_rects=[], niji_rects=[])
        self.assertEqual(fp.program_tiles([0, 359, 10.7, 349.3, 10.3], [2] * 5, MEGACAM_OUTER, MEGACAM_INSIDE).astype(int).tolist(), [1, 1, 0, 0, 1])

    def test_decam_box_uses_its_own_size(self):
        # DECam's 2.1 x 1.9 deg field: a tile 0.9 deg inside the Niji edge still reaches POLU
        # sky, while MegaCam's 1 x 1 deg field at the same place is Niji-only.
        fp = Footprint(polu_rects=[(125.0, 250.0, -1.0, 5.0)], excluded_rects=[(130.0, 180.0, 0.5, 3.5)],
                       niji_rects=[(130.0, 180.0, 0.5, 3.5)])
        self.assertEqual(fp.program_tiles([155, 155], [1.4, 2.0], DECAM_OUTER, DECAM_INSIDE).astype(int).tolist(), [1, 0])
        self.assertEqual(fp.program_tiles([155, 155], [1.4, 2.0], MEGACAM_OUTER, MEGACAM_INSIDE).astype(int).tolist(), [0, 0])

    def test_area(self):
        fp = Footprint(polu_rects=[(0.0, 10.0, 0.0, 10.0)], excluded_rects=[], niji_rects=[])
        self.assertAlmostEqual(fp.area(fp.polu), 10 * np.rad2deg(np.sin(np.deg2rad(10))), delta=1.0)


class LabelTests(unittest.TestCase):
    def test_program_plan_precedence_and_sticky_labels(self):
        # a: planned program tile; b: program tile; c: planned but outside IBIS (no label);
        # d: outside; e: keeps its LBNL label; f: left the program unobserved (cleared);
        # g: left the program but observed (keeps NAOC).
        table = Table({'OBJECT': list('abcdefg'),
                       'IN_IBIS': [1, 1, 0, 0, 1, 1, 1],
                       'IN_HSC': [0, 0, 0, 0, 0, 1, 1],
                       'DONE': [0, 0, 0, 0, 0, 0, 1],
                       'PROGRAM': ['', '', '', '', 'LBNL', 'NAOC', 'NAOC'],
                       'RA': [180] * 7, 'DEC': [0] * 7})
        np.testing.assert_array_equal(labels(table, {'a', 'c', 'e', 'f'}),
                                      ['NAOC', 'LBNL', '', '', 'LBNL', '', 'NAOC'])
        table.remove_column('PROGRAM')
        np.testing.assert_array_equal(labels(table, {'a', 'c'}),
                                      ['NAOC', 'LBNL', '', '', 'LBNL', '', ''])

    def test_hsc_flags_only_for_in_ibis_rows(self):
        fp = Footprint(polu_rects=[(100.0, 200.0, 0.0, 4.0)], excluded_rects=[(120.0, 180.0, 1.0, 3.0)],
                       niji_rects=[(120.0, 180.0, 1.0, 3.0)])
        table = Table({'OBJECT': list('abcd'), 'IN_IBIS': [1, 1, 0, 1],
                       'RA': [150, 110, 150, 150], 'DEC': [2.0, 2.0, 2.0, 8.0]})
        np.testing.assert_array_equal(hsc_flags(table, fp), [1, 0, 0, 1])


if __name__ == '__main__':
    unittest.main()
