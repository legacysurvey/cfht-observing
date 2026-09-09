"""Selection checks for the CFHT plan generator (standard library only)."""

import csv
from contextlib import redirect_stdout
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from make_cfht_targets import iter_matching_targets, main, parse_args, trim_non_overlapping


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
        self.write_table()

    def write_table(self, omit=None):
        indices = [i for i, c in enumerate(self.columns) if c != omit]
        with self.path.open('w', newline='') as handle:
            handle.write('# %ECSV 1.0\n# ---\n# datatype:\n')
            for i in indices:
                dtype = ('string' if self.columns[i] in ('OBJECT', 'FILTER', 'PROGRAM')
                         else 'float64' if self.columns[i] in ('RA', 'DEC', 'PRIORITY', 'LMST_DESIGN') else 'int16')
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
        self.rows[1][-1] = float('nan')
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


if __name__ == '__main__':
    unittest.main()
