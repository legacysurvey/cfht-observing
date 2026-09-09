#!/usr/bin/env python3
"""Populate tile PRIORITY from the YEAR values in obstatus/desi2-stripes.ecsv.

Every IN_IBIS=1 tile participates, regardless of DONE, IN_HSC, filter, or
program. Use the highest priority of all overlapping stripes. Other rows
retain their existing priority, or 1.0 if the column is being created.
"""

import argparse
from pathlib import Path
import re
import tempfile

import numpy as np
from astropy.io import fits
from astropy.table import Column, Table


ROOT = Path(__file__).resolve().parent
DEC_HALF_WIDTH = 1.6 + 0.5
FALLBACK_PRIORITY = 1.0


def priority_values(tiles, stripes, pad_all=True):
    """Return max(11-YEAR) within RA bounds and padded declination bounds."""
    if 'PRIORITY' in tiles.colnames:
        result = np.array(tiles['PRIORITY'], dtype=np.float64, copy=True)
    else:
        result = np.full(len(tiles), FALLBACK_PRIORITY, dtype=np.float64)
    eligible = np.asarray(tiles['IN_IBIS']) == 1
    result[eligible] = -np.inf
    ra = np.asarray(tiles['RA']) % 360
    dec = np.asarray(tiles['DEC'])
    for stripe in stripes:
        year_text = str(stripe['YEAR']).strip()
        if not re.fullmatch(r'Y[1-9][0-9]*', year_text):
            raise ValueError(f'Invalid YEAR {year_text!r} in stripe {stripe["STRIPE"]}')
        year = int(year_text[1:])
        priority = 11.0 - year
        ra_min, ra_max = float(stripe['RA_MIN']), float(stripe['RA_MAX'])
        if ra_max - ra_min >= 360:
            inside_ra = np.ones(len(tiles), dtype=bool)
        else:
            width = (ra_max - ra_min) % 360
            offset = (ra - ra_min) % 360
            inside_ra = offset <= width + 1e-10
        if pad_all or year == 1:
            low = float(stripe['DEC_CENTER']) - DEC_HALF_WIDTH
            high = float(stripe['DEC_CENTER']) + DEC_HALF_WIDTH
        else:
            low, high = float(stripe['DEC_MIN']), float(stripe['DEC_MAX'])
        inside = eligible & inside_ra & (dec >= low - 1e-10) & (dec <= high + 1e-10)
        result[inside] = np.maximum(result[inside], priority)
    result[eligible & np.isneginf(result)] = FALLBACK_PRIORITY
    return result


def update(root=ROOT, pad_all=True):
    stripes = Table.read(root / 'obstatus/desi2-stripes.ecsv')
    ecsv_path = root / 'obstatus/cfht-tiles.ecsv'
    fits_path = root / 'obstatus/cfht-tiles.fits'
    ecsv = Table.read(ecsv_path)
    original_fits = Table.read(fits_path)
    for column in ('RA', 'DEC', 'IN_IBIS'):
        np.testing.assert_array_equal(ecsv[column], original_fits[column])
    np.testing.assert_array_equal(np.char.strip(np.asarray(ecsv['OBJECT'], dtype=str)),
                                  np.char.strip(np.asarray(original_fits['OBJECT'], dtype=str)))
    values = priority_values(ecsv, stripes, pad_all)
    fits_values = priority_values(original_fits, stripes, pad_all)
    np.testing.assert_array_equal(values, fits_values)
    ecsv['PRIORITY'] = Column(values, description=(
        'IN_IBIS=1: highest stripe priority (Y1=10, Y2=9, Y3=8, Y4=7, Y5=6); '
        'fallback 1; RA bounds from desi2-stripes.ecsv; '
        + ('all stripes' if pad_all else 'Y1 stripes') + ' use DEC_CENTER +/-2.1 deg'))
    with tempfile.TemporaryDirectory(prefix='cfht-priority-', dir=root) as tmp:
        new_ecsv, new_fits = Path(tmp) / ecsv_path.name, Path(tmp) / fits_path.name
        ecsv.write(new_ecsv, format='ascii.ecsv')
        with fits.open(fits_path) as hdus:
            if 'PRIORITY' in hdus[1].columns.names:
                hdus[1].data['PRIORITY'] = fits_values
            else:
                columns = hdus[1].columns + fits.ColDefs([
                    fits.Column(name='PRIORITY', format='D', array=fits_values)])
                hdus[1] = fits.BinTableHDU.from_columns(columns, header=hdus[1].header)
            hdus[1].header['PRIMETH'] = ('max(11-YEAR)', 'IN_IBIS=1, max of overlapping stripes')
            hdus[1].header['PRIPAD'] = ('ALL' if pad_all else 'Y1', 'DEC_CENTER +/-2.1 deg')
            hdus.writeto(new_fits)
        for path, original, expected in ((new_ecsv, ecsv, values),
                                          (new_fits, original_fits, fits_values)):
            saved = Table.read(path)
            for column in original.colnames:
                if column != 'PRIORITY':
                    np.testing.assert_array_equal(saved[column], original[column])
            np.testing.assert_array_equal(saved['PRIORITY'], expected)
        new_ecsv.replace(ecsv_path)
        new_fits.replace(fits_path)
    unique, counts = np.unique(values[ecsv['IN_IBIS'] == 1], return_counts=True)
    for value, count in zip(unique[::-1], counts[::-1]):
        print(f'IN_IBIS=1, PRIORITY={value:.1f}: {count} rows')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--y1-only-padding', action='store_true',
                        help='Use +/-2.1 deg for Y1 only, and tabulated DEC_MIN/MAX for later years')
    args = parser.parse_args()
    update(args.root, pad_all=not args.y1_only_padding)
