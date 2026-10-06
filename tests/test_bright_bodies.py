import unittest
import warnings

import numpy as np
from astropy import units as u
from astropy.coordinates import get_body, get_sun

from make_cfht_targets import TileTarget
from night_planning import (AVOIDED_PLANETS, BRIGHT_BODIES, CFHT, bright_body_positions,
                            exclude_near_planets, moon_illumination, night_bounds, offline_iers)


class BrightBodyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            cls.night = night_bounds('2026-10-07')
            cls.positions = bright_body_positions(cls.night, steps=5)

    def test_all_bodies_sampled_through_night(self):
        self.assertEqual(set(self.positions), set(BRIGHT_BODIES))
        for ra, dec in self.positions.values():
            self.assertEqual(ra.shape, (5,))
            self.assertTrue(np.all((ra >= 0) & (ra < 360)))
            self.assertTrue(np.all(np.abs(dec) <= 90))

    def test_moon_is_geocentric_direction_not_barycentric(self):
        # The Moon moves roughly half a degree per hour; a barycentric ICRS
        # direction would barely move and would track the Earth instead.
        ra, dec = self.positions['moon']
        motion = np.hypot((ra[-1] - ra[0] + 180) % 360 - 180, dec[-1] - dec[0])
        self.assertGreater(motion, 3.0)
        self.assertLess(motion, 9.0)
        with offline_iers():
            sun = get_sun(self.night.start)
            moon = get_body('moon', self.night.start, location=CFHT)
        # Waning crescent that night: Moon within ~40 degrees of the Sun.
        self.assertLess(moon.separation(sun).deg, 40.0)
        self.assertAlmostEqual(ra[0], moon.ra.deg, places=6)

    def test_exclusion_uses_planet_track(self):
        saturn_ra, saturn_dec = self.positions['saturn']
        near = TileTarget('NEAR', float(saturn_ra[2]) + 0.5, float(saturn_dec[2]), 10.0, 10.0)
        far = TileTarget('FAR', float(saturn_ra[2]) + 5.0, float(saturn_dec[2]), 10.0, 10.0)
        kept, excluded = exclude_near_planets([near, far], self.positions, 1.0)
        self.assertEqual([t.object_name for t in kept], ['FAR'])
        self.assertEqual(excluded, {'NEAR': 'saturn'})
        kept, excluded = exclude_near_planets([near, far], self.positions, 0.0)
        self.assertEqual(len(kept), 2)
        self.assertEqual(excluded, {})
        with self.assertRaises(ValueError):
            exclude_near_planets([near], self.positions, -1.0)

    def test_moon_illumination(self):
        from astropy.time import Time
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            # Known full Moon 2026-01-03 10:03 UT and new Moon 2026-01-18 19:52 UT.
            self.assertGreater(moon_illumination(Time('2026-01-03T10:03:00')), 0.99)
            self.assertLess(moon_illumination(Time('2026-01-18T19:52:00')), 0.01)
            middle = self.night.start + (self.night.end - self.night.start) / 2
            fraction = moon_illumination(middle)
        # Waning crescent on 2026-10-07/08, about 2 days before new Moon (about 5%).
        self.assertGreater(fraction, 0.02)
        self.assertLess(fraction, 0.25)

    def test_avoided_planets_are_the_bright_four(self):
        self.assertEqual(set(AVOIDED_PLANETS), {'venus', 'mars', 'saturn', 'uranus'})


if __name__ == '__main__':
    unittest.main()
