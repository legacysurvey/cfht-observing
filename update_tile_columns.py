#!/usr/bin/env python3
"""Update PROGRAM and IN_HSC in the CFHT ECSV tile catalog from the imaging plan.

Requires numpy and astropy. IN_HSC=1 marks IN_IBIS=1 tiles that CFHT should
not observe because they are not MegaCam program tiles under footprint.py:
either their center is more than one degree from POLU+IBIS sky (for
example the whole of SGC-1, which has HSC-Niji+VST imaging only), or their
1 x 1 degree field lies entirely inside an HSC-Niji+VST rectangle.
PROGRAM is NAOC for exact OBJECT matches to target NAME values in plans/,
otherwise LBNL for program tiles (IN_IBIS=1, IN_HSC=0); labels already
present are kept, except that tiles outside the program with no exposure
(IN_HSC=1 or IN_IBIS=0, DONE=0) carry no label: they are not to be observed.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import tempfile
import xml.etree.ElementTree as ET

import numpy as np
from astropy.table import Table

from footprint import MEGACAM_INSIDE, MEGACAM_OUTER, Footprint

ROOT = Path(__file__).resolve().parent


def planned_objects(plans: Path) -> set[str]:
    """Collect exact target names from both FITS and XML plans."""
    names = set()
    for path in sorted(plans.glob('*.fits')):
        names.update(str(name).strip() for name in Table.read(path)['NAME'])
    for path in sorted(plans.glob('*.xml')):
        for csv in ET.parse(path).getroot().iter('CSV'):
            lines = [line.strip() for line in (csv.text or '').splitlines()
                     if line.strip()]
            for line in lines[int(csv.attrib['headlines']):]:
                name = line.split(csv.attrib['colsep'])[0].strip()
                if name:
                    names.add(name)
    if not names:
        raise ValueError(f'No planned targets found in {plans}')
    return names


def hsc_flags(table: Table, footprint: Footprint) -> np.ndarray:
    """IN_HSC: IN_IBIS=1 tiles that are not MegaCam program tiles."""
    in_ibis = np.asarray(table['IN_IBIS']) == 1
    program = footprint.program_tiles(table['RA'], table['DEC'], MEGACAM_OUTER, MEGACAM_INSIDE)
    return (in_ibis & ~program).astype(np.int16)


def labels(table: Table, names: set[str], in_hsc=None) -> np.ndarray:
    """PROGRAM: keep existing labels; NAOC from plans; LBNL for other program tiles.

    ``in_hsc`` is the IN_HSC column to use (the freshly computed one); rows
    outside the program (IN_IBIS=0 or IN_HSC=1) that have no exposure are
    cleared to '' so that only observed tiles keep a label there.
    """
    objects = np.char.strip(np.asarray(table['OBJECT'], dtype=str))
    missing = names - set(objects)
    if missing:
        raise ValueError(f'{len(missing)} planned targets absent from tiles')
    if 'PROGRAM' in table.colnames:
        existing = table['PROGRAM']
        existing = existing.filled('') if hasattr(existing, 'filled') else existing
        program = np.char.strip(np.asarray(existing, dtype=str)).astype('U4')
    else:
        program = np.full(len(table), '', dtype='U4')
    in_ibis = np.asarray(table['IN_IBIS']) == 1
    if in_hsc is None:
        in_hsc = np.asarray(table['IN_HSC']) == 1 if 'IN_HSC' in table.colnames else np.zeros(len(table), bool)
    in_program = in_ibis & ~np.asarray(in_hsc).astype(bool)
    done = np.asarray(table['DONE']) != 0 if 'DONE' in table.colnames else np.zeros(len(table), bool)
    program[~in_program & ~done] = ''
    unlabeled = program == ''
    program[unlabeled & in_program] = 'LBNL'
    program[unlabeled & in_program & np.isin(objects, sorted(names))] = 'NAOC'
    return program


def update(root: Path) -> None:
    names = planned_objects(root / 'plans')
    footprint = Footprint.from_plan(root)
    ecsv_path = root / 'obstatus/cfht-tiles.ecsv'
    ecsv = Table.read(ecsv_path)
    original_columns = [c for c in ecsv.colnames if c not in ('PROGRAM', 'IN_HSC')]
    in_hsc = hsc_flags(ecsv, footprint)
    program = labels(ecsv, names, in_hsc)
    before = np.asarray(ecsv['IN_HSC']) if 'IN_HSC' in ecsv.colnames else None
    ecsv['PROGRAM'] = program
    ecsv['IN_HSC'] = in_hsc

    # Stage and validate the catalog before replacing the original.
    with tempfile.TemporaryDirectory(prefix='cfht-columns-', dir=root) as temp:
        staged_ecsv = Path(temp) / ecsv_path.name
        ecsv.write(staged_ecsv, format='ascii.ecsv')
        saved_ecsv = Table.read(staged_ecsv)
        for column in original_columns:
            np.testing.assert_array_equal(saved_ecsv[column], ecsv[column])
        np.testing.assert_array_equal(saved_ecsv['PROGRAM'].filled('')
                                      if hasattr(saved_ecsv['PROGRAM'], 'filled')
                                      else saved_ecsv['PROGRAM'], program)
        np.testing.assert_array_equal(saved_ecsv['IN_HSC'], in_hsc)
        staged_ecsv.replace(ecsv_path)
    for value in ('NAOC', 'LBNL', ''):
        print(f'PROGRAM={value or "(blank)"}: {np.count_nonzero(program == value)}')
    print(f'IN_HSC=1: {np.count_nonzero(in_hsc)}')
    if before is not None:
        print(f'IN_HSC changed 0->1: {np.count_nonzero((before == 0) & (in_hsc == 1))}, '
              f'1->0: {np.count_nonzero((before == 1) & (in_hsc == 0))}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    args = parser.parse_args()
    update(args.root)
