"""Extinction corrections in ordinary and nightly queue products."""

from contextlib import redirect_stdout, redirect_stderr
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from astropy.table import Table

from make_cfht_targets import (DEFAULT_MAG_AB, TileTarget, main, parse_args,
                               target_magnitude)
from night_planning import night_bounds


class MagnitudeTests(unittest.TestCase):
    def test_filter_defaults_and_override(self):
        for band, coefficient in [('M4112', 3.774), ('M4376', 3.576)]:
            tile = TileTarget('target', 0, 0, filter_name=band, ebv_med=0.1)
            self.assertAlmostEqual(target_magnitude(tile), 24.25 + coefficient * 0.1)
            self.assertAlmostEqual(target_magnitude(tile, 25.0), 25.0 + coefficient * 0.1)
            with patch.dict(DEFAULT_MAG_AB, {band: 24.0}):
                self.assertAlmostEqual(target_magnitude(tile), 24.0 + coefficient * 0.1)
            self.assertEqual(target_magnitude(TileTarget('clear', 0, 0, filter_name=band)), 24.25)

    def test_invalid_configuration(self):
        with self.assertRaisesRegex(ValueError, 'No extinction coefficient'):
            target_magnitude(TileTarget('unknown', 0, 0, filter_name='unknown'))
        for options in (['--filter', 'unknown'], ['--mag-ab', 'nan'], ['--mag-ab', 'inf']):
            with patch('sys.argv', ['make_cfht_targets.py'] + options), redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    parse_args()

    def test_written_magnitudes_in_both_filters_and_plan_modes(self):
        night = night_bounds('2026-09-09')
        ra = night.lmst_start
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            for band, coefficient in [('M4112', 3.774), ('M4376', 3.576)]:
                table = Table({'OBJECT': ['clear', 'dusty'], 'RA': [ra, ra],
                               'DEC': [0., 1.], 'FILTER': [band, band],
                               'IN_IBIS': [1, 1], 'IN_HSC': [0, 0], 'DONE': [0, 0],
                               'PRIORITY': [10., 10.], 'PROGRAM': ['LBNL', 'LBNL'],
                               'EBV_MED': [0., 0.1], 'LMST_DESIGN': [ra, ra]})
                source = root / (band + '.ecsv')
                table.write(source, format='ascii.ecsv')
                expected = {'clear': 24.25, 'dusty': 24.25 + coefficient * 0.1}
                for mode in ('ordinary', 'night'):
                    with self.subTest(band=band, mode=mode):
                        prefix = band + '-' + mode
                        argv = ['make_cfht_targets.py', '--input', str(source),
                                '--filter', band.lower(), '--outdir', str(root), '--prefix', prefix]
                        if mode == 'night':
                            argv += ['--night', '2026-09-09']
                        with patch('sys.argv', argv), redirect_stdout(io.StringIO()):
                            main()
                        fits = Table.read(root / (prefix + '.fits'))
                        self.assertEqual(len(fits), 2)
                        names = np.char.strip(np.asarray(fits['NAME'], dtype=str))
                        np.testing.assert_allclose(fits['MAG_AB'], [expected[n] for n in names],
                                                   rtol=0, atol=2e-6)
                        xml = (root / (prefix + '.xml')).read_text()
                        for name, magnitude in expected.items():
                            row = next(line for line in xml.splitlines() if line.startswith(name + ' '))
                            self.assertEqual(row.split('|')[3].strip(), f'{magnitude:.2f}')
                        if mode == 'night':
                            audit = Table.read((root / (prefix + '.schedule.ecsv')).read_text().splitlines(),
                                               format='ascii.ecsv')
                            scheduled = audit['STATUS'] == 'SCHEDULED'
                            self.assertEqual(int(np.sum(scheduled)), 2)
                            np.testing.assert_allclose(audit['MAG_AB'][scheduled],
                                                       [expected[n] for n in audit['OBJECT'][scheduled]])
                            self.assertTrue(np.all(np.isnan(audit['MAG_AB'][~scheduled])))


if __name__ == '__main__':
    unittest.main()
