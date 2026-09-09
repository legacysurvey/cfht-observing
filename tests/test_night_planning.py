"""Night dates, twilight roots, timing budgets, and selection priorities."""

import unittest
import warnings
from unittest.mock import patch

import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.utils import iers

from make_cfht_targets import TileTarget
from night_planning import (CFHT, HST, Night, night_bounds, offline_iers,
                            pick_target, schedule_night, slot_count, solar_altitude)


class NightTests(unittest.TestCase):
    def test_current_night_offline_without_warnings(self):
        # In particular, exercise this with the user's Astropy 4.3.1.
        with patch('astropy.utils.data.download_file', side_effect=AssertionError('network')), \
                patch('astropy.utils.iers.iers.download_file', side_effect=AssertionError('network')), \
                warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            night = night_bounds('2026-09-09')
            with offline_iers():
                self.assertFalse(iers.conf.auto_download)
                self.assertEqual(iers.conf.auto_max_age, 30)
            self.assertEqual(slot_count(night.seconds), 198)
        self.assertEqual(caught, [])

    def test_out_of_coverage_rejected_before_solar_calculation(self):
        with patch('night_planning.solar_altitude', side_effect=AssertionError('too early')):
            with self.assertRaisesRegex(ValueError, 'project IERS coverage'):
                night_bounds('2099-09-09')

    def test_nine_hour_capacity(self):
        self.assertEqual(slot_count(9 * 3600), 186)
        self.assertEqual(9 * 3600 - 186 * (130 + 44), 36)
        self.assertEqual(slot_count(173), 0)
        self.assertEqual(slot_count(348), 2)
        self.assertEqual(slot_count(348, 174, 0), 2)
        for args in [(100, 0, 44), (100, 130, -1), (np.nan, 130, 44)]:
            with self.assertRaises(ValueError):
                slot_count(*args)

    def test_hawaii_date_twilight_and_sidereal_wrap(self):
        night = night_bounds('2026-09-09')
        self.assertEqual(night.start.to_datetime(timezone=HST).date().isoformat(), '2026-09-09')
        self.assertEqual(night.end.to_datetime(timezone=HST).date().isoformat(), '2026-09-10')
        self.assertEqual(night.start.utc.isot[:10], '2026-09-10')
        self.assertTrue(9 * 3600 < night.seconds < 10 * 3600)
        self.assertGreater(night.lmst_start, night.lmst_end)
        delta = (night.lmst_end - night.lmst_start) % 360
        self.assertAlmostEqual(delta, night.seconds / 3600 * 15 * 1.002737909, places=4)
        with offline_iers():
            for time in (night.start, night.end):
                self.assertAlmostEqual(solar_altitude(time), -15, places=4)
            self.assertGreater(solar_altitude(night.start - 60 * u.s), -15)
            self.assertLess(solar_altitude(night.start + 60 * u.s), -15)
            self.assertLess(solar_altitude(night.end - 60 * u.s), -15)
            self.assertGreater(solar_altitude(night.end + 60 * u.s), -15)
        deeper = night_bounds('2026-09-09', 18)
        self.assertGreater(deeper.start, night.start)
        self.assertLess(deeper.end, night.end)

    def test_southern_priority_within_window(self):
        targets = [TileTarget('north', 0, 10, 0),
                   TileTarget('south', 0, -10, 359.5),
                   TileTarget('distant', 0, -20, 180)]
        self.assertEqual(pick_target(targets, 0, {0, 1, 2}, 2), [1, 0])
        self.assertEqual(pick_target(targets, 0, {0, 1, 2}, 0.1), [0])

    def test_survey_priority_precedes_declination_and_lmst_distance(self):
        targets = [TileTarget('south', 0, -10, 0, priority=9),
                   TileTarget('north', 0, 10, 355, priority=10),
                   TileTarget('outside', 0, -20, 354.999, priority=20)]
        self.assertEqual(pick_target(targets, 0, {0, 1, 2}, 5), [1, 0])
        self.assertEqual(pick_target(targets, 0, {0, 1, 2}, 2), [0])

    def test_scheduler_uses_five_degree_window_and_records_priority(self):
        night = self.short_night()
        lmst = night.lmst_start
        targets = [TileTarget('south', lmst, -10, lmst, priority=9),
                   TileTarget('north', lmst, 10, (lmst+4) % 360, priority=10)]
        selected, audit = schedule_night(targets, night)
        self.assertEqual([t.object_name for t in selected], ['north', 'south'])
        self.assertEqual(list(audit['PRIORITY']), [10.0, 9.0])
        self.assertEqual(audit.meta['LMST_WINDOW'], 5.0)

    def short_night(self, seconds=348):
        start = Time('2026-09-10T08:00:00')
        with offline_iers():
            lmst = start.sidereal_time('mean', longitude=CFHT.lon).deg
        return Night('2026-09-09', 15, start, start + seconds * u.s, lmst,
                     (lmst + seconds / 3600 * 15.041) % 360)

    def test_slots_order_no_repeats_visibility_and_shortfall(self):
        night = self.short_night()
        lmst = night.lmst_start
        targets = [TileTarget('north', lmst, 10, lmst),
                   TileTarget('south', lmst, -10, lmst),
                   TileTarget('below', (lmst + 180) % 360, -20, lmst)]
        selected, audit = schedule_night(targets, night)
        self.assertEqual([t.object_name for t in selected], ['south', 'north'])
        self.assertEqual(len(audit), 2)
        self.assertTrue(np.all(audit['MIN_ALT'] > 0))
        self.assertLessEqual(Time(audit['SLOT_END_UTC'][-1]), night.end + .001 * u.s)
        short, gaps = schedule_night(targets[:1], night)
        self.assertEqual(len(short), 1)
        self.assertEqual(list(gaps['STATUS']), ['SCHEDULED', 'NO_NEARBY_TARGET'])
        self.assertAlmostEqual(gaps.meta['UNUSED_SECONDS'], 174, places=5)

    def test_nonoverlap_applied_after_southern_priority(self):
        night = self.short_night()
        lmst = night.lmst_start
        targets = [TileTarget('north', lmst, 0.1, lmst),
                   TileTarget('south', lmst, 0, lmst)]
        selected, audit = schedule_night(targets, night, non_overlapping=True)
        self.assertEqual([t.object_name for t in selected], ['south'])
        self.assertEqual(audit.meta['CAPACITY'], 2)

    def test_zero_capacity_and_no_candidates(self):
        selected, audit = schedule_night([], self.short_night(100))
        self.assertEqual(len(selected), 0)
        self.assertEqual(len(audit), 0)
        selected, audit = schedule_night([], self.short_night())
        self.assertEqual(len(selected), 0)
        self.assertEqual(len(audit), 2)


if __name__ == '__main__':
    unittest.main()
