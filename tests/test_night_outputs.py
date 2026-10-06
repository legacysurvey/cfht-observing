import tempfile
import unittest
from pathlib import Path

import numpy as np
from astropy.table import Table

from make_cfht_targets import TileTarget
from night_outputs import lmst_hms, plot_ra, read_tile_positions, write_schedule_text, plot_night


def sample_audit():
    audit = Table(rows=[(1, '2026-10-08T05:04:54.296', '2026-10-08T05:07:04.296',
                         '2026-10-08T05:07:48.296', 'DESI_M4376_1', 297.5, 300.0, 2.5, 50.0, 10.0,
                         'SCHEDULED'),
                        (2, '2026-10-08T05:07:48.296', '2026-10-08T05:09:58.296',
                         '2026-10-08T05:10:42.296', '', 298.2, np.nan, np.nan, np.nan, np.nan,
                         'NO_NEARBY_TARGET')],
                  names=('SLOT', 'START_UTC', 'EXPOSURE_END_UTC', 'SLOT_END_UTC', 'OBJECT', 'LMST',
                         'LMST_DESIGN', 'LMST_OFFSET', 'MIN_ALT', 'PRIORITY', 'STATUS'))
    audit['MAG_AB'] = [24.5, np.nan]
    audit.meta.update(NIGHT='2026-10-07', TWILIGHT=15.0, START_UTC='2026-10-08T05:04:54.296',
                      END_UTC='2026-10-08T15:14:12.896', LMST_START=297.605, LMST_END=90.35,
                      EXPTIME=130.0, OVERHEAD=44.0, CAPACITY=2, SCHEDULED=1,
                      MOON_ILLUMINATION=0.123)
    return audit


class NightOutputTests(unittest.TestCase):
    def test_lmst_hms(self):
        self.assertEqual(lmst_hms(0.0), '00:00:00')
        self.assertEqual(lmst_hms(297.605261), '19:50:25')
        self.assertEqual(lmst_hms(359.9999), '00:00:00')

    def test_plot_ra_runs_from_300_down_to_minus_60(self):
        np.testing.assert_allclose(plot_ra([300.0, 270.0, 0.0, 300.1, 359.0]),
                                   [300.0, 270.0, 0.0, -59.9, -1.0])

    def test_schedule_text_lists_ut_lmst_and_object(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'plan.txt'
            write_schedule_text(path, sample_audit())
            lines = path.read_text().splitlines()
        data = [line for line in lines if not line.startswith('#')]
        self.assertEqual(len(data), 2)
        self.assertEqual(data[0].split(), ['2026-10-08', '05:04:54', '19:50:00', '297.500', 'DESI_M4376_1'])
        self.assertIn('(NO_NEARBY_TARGET)', data[1])
        self.assertTrue(any('LMST at twilight: 19:50:25 to 06:01:24' in line for line in lines))
        self.assertTrue(any('Moon: 12% illuminated' in line for line in lines))

    def test_read_positions_and_plot(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ecsv = root / 'tiles.ecsv'
            ecsv.write_text('# %ECSV 1.0\n# ---\n'
                            '# datatype:\n'
                            '# - {name: OBJECT, datatype: string}\n'
                            '# - {name: RA, datatype: float64}\n'
                            '# - {name: DEC, datatype: float64}\n'
                            '# - {name: FILTER, datatype: string}\n'
                            '# - {name: IN_IBIS, datatype: int16}\n'
                            '# - {name: IN_HSC, datatype: int16}\n'
                            '# - {name: DONE, datatype: int16}\n'
                            '# schema: astropy-2.0\n'
                            'OBJECT RA DEC FILTER IN_IBIS IN_HSC DONE\n'
                            'DESI_M4376_1 10.0 -15.0 M4376 1 0 0\n'
                            'DESI_M4376_2 20.0 -15.0 M4376 1 0 1\n'
                            'DESI_M4112_3 20.0 -15.0 M4112 1 1 0\n')
            tiles = read_tile_positions(ecsv)
            self.assertEqual(len(tiles['RA']), 3)
            self.assertEqual(list(tiles['DONE']), [0, 1, 0])
            self.assertEqual(list(tiles['FILTER']), ['M4376', 'M4376', 'M4112'])
            pdf = root / 'plan.pdf'
            bodies = {'moon': (np.array([160.0, 165.0, 170.0]), np.array([4.0, 3.0, 2.0])),
                      'saturn': (np.array([10.8, 10.8, 10.8]), np.array([1.7, 1.7, 1.7])),
                      'mars': (np.array([128.0, 128.1, 128.2]), np.array([25.0, 25.0, 25.0]))}
            plot_night(pdf, sample_audit(), tiles,
                       [TileTarget('DESI_M4376_1', 10.0, -15.0, 300.0, 10.0)], 'M4376', bodies)
            self.assertGreater(pdf.stat().st_size, 1000)
            self.assertTrue(pdf.read_bytes().startswith(b'%PDF'))


if __name__ == '__main__':
    unittest.main()
