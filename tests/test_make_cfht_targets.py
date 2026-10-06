"""Selection checks for the CFHT plan generator (standard library only)."""

import csv
from contextlib import redirect_stdout
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from make_cfht_targets import (iter_matching_targets, main, parse_args, read_done_positions,
                               trim_non_overlapping)


class SelectionTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.path = Path(temp.name) / 'tiles.ecsv'
        self.columns = ['OBJECT', 'RA', 'DEC', 'FILTER', 'IN_IBIS',
                        'IN_HSC', 'DONE', 'PROGRAM']
        self.rows = [
            ['hsc', 359, 0, 'M4376', 1, 1, 0, 'LBNL'],
            ['lbl', 359.1, 0, 'M4376', 1, 0, 0, 'LBNL'],
            ['naoc', 0, 0, 'M4376', 1, 0, 0, 'NAOC'],
            ['blank', 1, 0, 'M4376', 1, 0, 0, ''],
            ['future', 2, 0, 'M4376', 1, 0, 0, 'Future Survey'],
            ['outside_ibis', 3, 0, 'M4376', 0, 0, 0, 'LBNL'],
            ['done', 4, 0, 'M4376', 1, 0, 1, 'LBNL'],
            ['wrong_filter', 5, 0, 'M4112', 1, 0, 0, 'LBNL'],
            ['wrong_ra', 180, 0, 'M4376', 1, 0, 0, 'LBNL'],
            ['wrong_dec', 6, 20, 'M4376', 1, 0, 0, 'LBNL'],
        ]
        self.columns.append('PRIORITY')
        for row in self.rows:
            row.append(10.0 if row[0] == 'lbl' else 9.0)
        self.columns.append('EBV_MED')
        for i, row in enumerate(self.rows):
            row.append(i * 0.01)
        self.write_table()

    def write_table(self, omit=None):
        indices = [i for i, c in enumerate(self.columns) if c != omit]
        with self.path.open('w', newline='') as handle:
            handle.write('# %ECSV 1.0\n# ---\n# datatype:\n')
            for i in indices:
                dtype = ('string' if self.columns[i] in ('OBJECT', 'FILTER', 'PROGRAM')
                         else 'float64' if self.columns[i] in ('RA', 'DEC', 'PRIORITY', 'LMST_DESIGN', 'EBV_MED') else 'int16')
                handle.write(f'# - {{name: {self.columns[i]}, datatype: {dtype}}}\n')
            writer = csv.writer(handle, delimiter=' ', quoting=csv.QUOTE_NONNUMERIC)
            writer.writerow([self.columns[i] for i in indices])
            writer.writerows([[row[i] for i in indices] for row in self.rows])

    def targets(self, programs=None):
        return list(iter_matching_targets(self.path, (350, 10), (-2, 2),
                                          'M4376', programs))

    def names(self, programs=None):
        return [target.object_name for target in self.targets(programs)]

    def test_default_includes_blank_and_future_programs_with_all_other_cuts(self):
        self.assertEqual(self.names(), ['lbl', 'naoc', 'blank', 'future'])

    def test_single_multiple_and_unknown_programs(self):
        self.assertEqual(self.names(['lbnl']), ['lbl'])
        self.assertEqual(self.names(['LBNL', 'naoc']), ['lbl', 'naoc'])
        self.assertEqual(self.names(['NotYetAssigned']), [])

    def test_naoc_case_and_duplicates(self):
        self.assertEqual(self.names(['NAOC']), ['naoc'])
        self.assertEqual(self.names(['NAOC', 'naoc']), ['naoc'])

    def test_blank_and_quoted_programs(self):
        self.assertEqual(self.names(['']), ['blank'])
        self.assertEqual(self.names(['', 'future survey']), ['blank', 'future'])

    def test_hsc_is_excluded_before_non_overlap_trim(self):
        # The earlier HSC target is close enough to suppress lbl if it is
        # allowed into the greedy trimming step.
        trimmed = trim_non_overlapping(self.targets(['LBNL']), 1.0)
        self.assertEqual([target.object_name for target in trimmed], ['lbl'])

    def test_avoid_done_blocks_overlap_with_observed_tiles(self):
        # 'done' at RA=4 is observed in M4376; 'blank' (RA=1) and 'future' (RA=2)
        # are not near it. Add an unobserved tile overlapping 'done'.
        self.rows.append(['near_done', 4.5, 0.2, 'M4376', 1, 0, 0, 'LBNL', 9.0, 0.0])
        self.rows.append(['done_other_filter', 7, 0, 'M4112', 1, 0, 1, 'LBNL', 9.0, 0.0])
        self.rows.append(['near_other_filter_done', 7.5, 0, 'M4376', 1, 0, 0, 'LBNL', 9.0, 0.0])
        self.write_table()
        blocked = read_done_positions(self.path, 'M4376')
        self.assertEqual(blocked, [(4.0, 0.0)])
        candidates = self.targets(['LBNL'])
        self.assertIn('near_done', [t.object_name for t in candidates])
        without = trim_non_overlapping(candidates, 1.0)
        with_block = trim_non_overlapping(candidates, 1.0, blocked)
        self.assertIn('near_done', [t.object_name for t in without])
        self.assertNotIn('near_done', [t.object_name for t in with_block])
        # DONE tiles in another filter do not block, and blocked positions are never returned.
        self.assertIn('near_other_filter_done', [t.object_name for t in with_block])
        self.assertTrue(all(t.object_name for t in with_block))

    def test_avoid_done_implies_non_overlapping(self):
        with patch('sys.argv', ['make_cfht_targets.py', '--avoid-done']):
            args = parse_args()
        self.assertTrue(args.avoid_done)
        self.assertTrue(args.non_overlapping)

    def test_missing_hsc_column_fails(self):
        self.write_table(omit='IN_HSC')
        with self.assertRaisesRegex(ValueError, 'missing required column.*IN_HSC'):
            self.targets()

    def test_program_column_required_only_for_program_selection(self):
        self.write_table(omit='PROGRAM')
        self.assertEqual(self.names(), ['lbl', 'naoc', 'blank', 'future'])
        with self.assertRaisesRegex(ValueError, 'missing required column.*PROGRAM'):
            self.targets(['LBNL'])

    def test_cli_program_lists_and_default(self):
        with patch('sys.argv', ['make_cfht_targets.py']):
            self.assertEqual(parse_args().programs, ['LBNL'])
            self.assertEqual(parse_args().lmst_window, 5.0)
        with patch('sys.argv', ['make_cfht_targets.py', '--all-programs']):
            self.assertIsNone(parse_args().programs)
        with patch('sys.argv', ['make_cfht_targets.py', '--program', 'LBNL', 'NAOC',
                                '--program', '', 'Future Survey']):
            self.assertEqual(parse_args().programs, ['LBNL', 'NAOC', '', 'Future Survey'])

    def test_iterator_defaults_to_lbnl(self):
        targets = iter_matching_targets(self.path, (350, 10), (-2, 2), 'M4376')
        self.assertEqual([t.object_name for t in targets], ['lbl'])

    def test_night_requires_design_column(self):
        with self.assertRaisesRegex(ValueError, 'missing required column.*LMST_DESIGN'):
            list(iter_matching_targets(self.path, (350, 10), (-2, 2), 'M4376',
                                       require_lmst=True))

    def test_priority_read_from_catalog(self):
        self.assertEqual([t.priority for t in self.targets()], [10.0, 9.0, 9.0, 9.0])
        self.write_table(omit='PRIORITY')
        with self.assertRaisesRegex(ValueError, 'missing required column.*PRIORITY'):
            self.targets()

    def test_nonfinite_priority_rejected(self):
        self.rows[1][self.columns.index('PRIORITY')] = float('nan')
        self.write_table()
        with self.assertRaisesRegex(ValueError, 'finite PRIORITY'):
            self.targets()

    def test_ordinary_nonoverlap_keeps_higher_priority(self):
        self.rows[0][self.columns.index('IN_HSC')] = 0
        self.write_table()
        argv = ['make_cfht_targets.py', '--input', str(self.path),
                '--outdir', str(self.path.parent), '--non-overlapping',
                '--ra', '350', '10', '--dec', '-2', '2']
        with patch('sys.argv', argv), patch('make_cfht_targets.write_xml'), \
                patch('make_cfht_targets.write_fits') as write, redirect_stdout(io.StringIO()):
            main()
        rows = write.call_args[0][1]
        self.assertEqual([row.name for row in rows], ['lbl'])

    def test_extinction_read_and_validated(self):
        self.assertEqual(self.targets()[0].ebv_med, 0.01)
        self.assertEqual(self.targets()[0].filter_name, 'M4376')
        for bad in ('not-a-number', float('nan'), float('inf')):
            self.rows[1][self.columns.index('EBV_MED')] = bad
            self.write_table()
            with self.assertRaisesRegex(ValueError, 'EBV_MED'):
                self.targets()
        self.write_table(omit='EBV_MED')
        with self.assertRaisesRegex(ValueError, 'missing required column.*EBV_MED'):
            self.targets()


if __name__ == '__main__':
    unittest.main()
