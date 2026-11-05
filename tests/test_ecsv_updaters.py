"""Tile updaters operate without a FITS catalog and preserve other columns."""

from contextlib import redirect_stdout
import io
from pathlib import Path
import tempfile
import unittest

import numpy as np
from astropy.table import Table

import update_lmst_design
import update_tile_columns
import update_tile_priorities


class EcsvUpdaterTests(unittest.TestCase):
    def test_updates_with_only_ecsv_catalog(self):
        cases = [(update_lmst_design.update, {}, {'LMST_DESIGN'}),
                 (update_tile_priorities.update, {}, {'PRIORITY'}),
                 (update_tile_columns.update, {}, {'PROGRAM', 'IN_HSC', 'IN_IBIS'})]
        for update, kwargs, changed_columns in cases:
            with self.subTest(updater=update.__module__), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                (root / 'obstatus').mkdir()
                (root / 'plans').mkdir()
                path = root / 'obstatus/cfht-tiles.ecsv'
                original = Table({'OBJECT': ['planned', 'done', 'remaining', 'outside'],
                                  'RA': [140., 150., 160., 170.], 'DEC': [2., 6., -3., 10.],
                                  'FILTER': ['M4376'] * 4, 'IN_IBIS': [1, 1, 1, 0],
                                  'IN_HSC': [0, 0, 0, 0], 'DONE': [0, 1, 0, 0],
                                  'PROGRAM': ['LBNL', 'NAOC', 'LBNL', ''],
                                  'PRIORITY': [7.] * 4, 'LMST_DESIGN': [np.nan] * 4,
                                  'EBV_MED': [0.01, 0.02, 0.03, 0.04]},
                                 meta={'NOTE': 'preserve catalog metadata'})
                original.write(path, format='ascii.ecsv')
                stripes = Table({'STRIPE': ['NGC-5'], 'DEC_CENTER': [2.],
                                 'RA_MIN': [130.], 'RA_MAX': [180.],
                                 'DEC_MIN': [-4.], 'DEC_MAX': [7.], 'YEAR': ['Y1']})
                stripes.write(root / 'obstatus/desi2-stripes.ecsv', format='ascii.ecsv')
                # The whole of RA 130-180, DEC 1-3 is an HSC-Niji+VST rectangle:
                # 'planned' (DEC 2) is Niji-only; 'done' (DEC 6) and 'remaining' (DEC -3)
                # are within a degree of the remaining POLU+IBIS sky.
                Table({'NAME': ['NGC-5'], 'KIND': ['niji'], 'RA_MIN': [130.], 'RA_MAX': [180.],
                       'DEC_MIN': [1.], 'DEC_MAX': [3.]}).write(
                    root / 'obstatus/footprint-regions.ecsv', format='ascii.ecsv')
                (root / 'plans/existing.xml').write_text(
                    '<ASTRO><CSV headlines="1" colsep="|">NAME|\nplanned|\n</CSV></ASTRO>')
                with redirect_stdout(io.StringIO()):
                    update(root, **kwargs)
                saved = Table.read(path.read_text().splitlines(), format='ascii.ecsv')
                self.assertEqual(saved.colnames, original.colnames)
                self.assertEqual(saved.meta, original.meta)
                for column in original.colnames:
                    if column not in changed_columns:
                        actual = saved[column]
                        if column == 'PROGRAM' and hasattr(actual, 'filled'):
                            actual = actual.filled('')
                        np.testing.assert_array_equal(actual, original[column])
                self.assertEqual(list(root.rglob('*.fits')), [])
                if 'PRIORITY' in changed_columns:
                    np.testing.assert_array_equal(saved['PRIORITY'], [10., 1., 1., 7.])
                elif 'PROGRAM' in changed_columns:
                    # 'planned' is now Niji-only and unobserved, so its label is cleared;
                    # the other labels are kept.
                    np.testing.assert_array_equal(saved['PROGRAM'].filled(''),
                                                  ['', 'NAOC', 'LBNL', ''])
                    np.testing.assert_array_equal(saved['IN_HSC'], [1, 0, 0, 0])
                    # 'outside' at DEC 10 is beyond the stripe (DEC -4..7): IN_IBIS=0.
                    np.testing.assert_array_equal(saved['IN_IBIS'], [1, 1, 1, 0])
                else:
                    np.testing.assert_array_equal(np.isfinite(saved['LMST_DESIGN']),
                                                  [True, False, True, False])


if __name__ == '__main__':
    unittest.main()
