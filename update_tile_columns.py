#!/usr/bin/env python3
"""Add PROGRAM and IN_HSC to the CFHT ECSV and FITS tile tables.

Requires numpy and astropy. Survey assignments must be supplied explicitly;
stripe coordinates come from obstatus/desi2-stripes.ecsv.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import tempfile
import xml.etree.ElementTree as ET

import numpy as np
from astropy.io import fits
from astropy.table import Table


ROOT = Path(__file__).resolve().parent
HALF_HEIGHT = 1.4
MARGIN = 1.0


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


def hsc_mask(ra, dec, stripes, northern_limits=None) -> np.ndarray:
    """Inset the union of adjoining stripes by one coordinate degree.

    RA margins use coordinate degrees, including at RA=0. Join declination
    intervals only where the tile's full RA margin fits in every stripe in
    the run. This retains shared interiors while protecting exposed steps
    caused by unequal RA limits. Boundaries are inclusive.
    Optional northern_limits map FIELD names to final, already-inset
    northern edges; only the named stripes are extended.
    """
    ra, dec = np.broadcast_arrays(np.asarray(ra), np.asarray(dec))
    result = np.zeros(ra.shape, dtype=bool)
    stripes = sorted(stripes, key=lambda stripe: float(stripe['DEC']))
    low = np.full(ra.shape, np.inf)
    high = np.full(ra.shape, -np.inf)
    for stripe in stripes:
        width = (float(stripe['RA_MAX']) - float(stripe['RA_MIN'])) % 360
        offset = (ra - float(stripe['RA_MIN'])) % 360
        eligible = (offset >= MARGIN) & (offset <= width - MARGIN)
        next_low = float(stripe['DEC']) - HALF_HEIGHT
        next_high = float(stripe['DEC']) + HALF_HEIGHT
        if northern_limits and stripe['FIELD'] in northern_limits:
            next_high = max(next_high, northern_limits[stripe['FIELD']] + MARGIN)
        separate = eligible & (next_low > high + 1e-10)
        result |= separate & (dec >= low + MARGIN - 1e-10) & (
            dec <= high - MARGIN + 1e-10)
        low = np.where(separate, next_low, low)
        high = np.where(eligible, np.maximum(high, next_high), high)
    result |= (dec >= low + MARGIN - 1e-10) & (dec <= high - MARGIN + 1e-10)
    return result.astype(np.int16)


def labels(table: Table, names: set[str], stripes, extend_hsc_north=False):
    objects = np.char.strip(np.asarray(table['OBJECT'], dtype=str))
    missing = names - set(objects)
    if missing:
        raise ValueError(f'{len(missing)} planned targets absent from tiles')
    program = np.full(len(table), '', dtype='U4')
    program[np.asarray(table['IN_IBIS']) == 1] = 'LBNL'
    program[np.isin(objects, sorted(names))] = 'NAOC'
    northern_limits = None
    if extend_hsc_north:
        if not {'NGC-5', 'SGC-1'} <= {str(s['FIELD']) for s in stripes}:
            raise ValueError('The northern extension requires NGC-5 and SGC-1')
        ra = np.asarray(table['RA']) % 360
        sgc_ibis = ((ra >= 270) | (ra < 90)) & (table['IN_IBIS'] == 1)
        if not np.any(sgc_ibis):
            raise ValueError('No SGC IN_IBIS=1 tiles to define the northern edge')
        northern_limits = {'NGC-5': 4.0,
                           'SGC-1': float(np.max(table['DEC'][sgc_ibis]))}
    return program, hsc_mask(table['RA'], table['DEC'], stripes, northern_limits)


def update(root: Path, fields: list[str], extend_hsc_north=False) -> None:
    stripes = Table.read(root / 'obstatus/desi2-stripes.ecsv')
    stripes.rename_column('STRIPE', 'FIELD')
    stripes.rename_column('DEC_CENTER', 'DEC')
    unknown = set(fields) - set(stripes['FIELD'])
    if unknown:
        raise ValueError(f'Unknown stripe names: {sorted(unknown)}')
    stripes = stripes[np.isin(stripes['FIELD'], fields)]
    names = planned_objects(root / 'plans')
    ecsv_path = root / 'obstatus/cfht-tiles.ecsv'
    fits_path = root / 'obstatus/cfht-tiles.fits'
    ecsv = Table.read(ecsv_path)
    original_columns = [c for c in ecsv.colnames if c not in ('PROGRAM', 'IN_HSC')]
    program, in_hsc = labels(ecsv, names, stripes, extend_hsc_north)
    ecsv['PROGRAM'] = program
    ecsv['IN_HSC'] = in_hsc

    # Stage and validate both products before replacing either original.
    with tempfile.TemporaryDirectory(prefix='cfht-columns-', dir=root) as temp:
        staged_ecsv = Path(temp) / ecsv_path.name
        staged_fits = Path(temp) / fits_path.name
        ecsv.write(staged_ecsv, format='ascii.ecsv')
        original_fits = Table.read(fits_path)
        with fits.open(fits_path) as hdus:
            fits_program, fits_hsc = labels(original_fits, names, stripes,
                                           extend_hsc_north)
            np.testing.assert_array_equal(ecsv['OBJECT'], original_fits['OBJECT'])
            np.testing.assert_array_equal(program, fits_program)
            np.testing.assert_array_equal(in_hsc, fits_hsc)
            retained = [c for c in hdus[1].columns
                        if c.name not in ('PROGRAM', 'IN_HSC')]
            columns = fits.ColDefs(retained) + fits.ColDefs([
                fits.Column(name='PROGRAM', format='4A', array=fits_program),
                fits.Column(name='IN_HSC', format='I', array=fits_hsc),
            ])
            hdus[1] = fits.BinTableHDU.from_columns(columns, header=hdus[1].header)
            hdus.writeto(staged_fits)

        saved_ecsv = Table.read(staged_ecsv)
        saved_fits = Table.read(staged_fits)
        for column in original_columns:
            np.testing.assert_array_equal(saved_ecsv[column], ecsv[column])
            np.testing.assert_array_equal(saved_fits[column], original_fits[column])
        for saved in (saved_ecsv, saved_fits):
            np.testing.assert_array_equal(saved['PROGRAM'].filled('')
                                          if hasattr(saved['PROGRAM'], 'filled')
                                          else saved['PROGRAM'], program)
            np.testing.assert_array_equal(saved['IN_HSC'], in_hsc)
        staged_ecsv.replace(ecsv_path)
        staged_fits.replace(fits_path)
    for value in ('NAOC', 'LBNL', ''):
        print(f'PROGRAM={value or "(blank)"}: {np.count_nonzero(program == value)}')
    print(f'IN_HSC=1: {np.count_nonzero(in_hsc)}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--hsc-stripes', nargs='+', required=True,
                        help='HSC+VST STRIPE names from desi2-stripes.ecsv')
    parser.add_argument('--extend-hsc-north', action='store_true',
                        help='Temporarily extend NGC-5 to DEC=4 and SGC-1 to '
                             'the northernmost SGC IN_IBIS=1 tile; omit to '
                             'restore the original stripe footprint')
    args = parser.parse_args()
    update(args.root, args.hsc_stripes, args.extend_hsc_north)
