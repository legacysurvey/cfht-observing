#!/usr/bin/env python3
"""Rebalance LMST_DESIGN in both tile tables using remaining ECSV targets.

Requires numpy and astropy. Angles are degrees. Re-run as DONE flags or
footprints change; PROGRAM labels do not affect the design distribution.
"""

import argparse
from pathlib import Path
import tempfile

import numpy as np
from astropy.io import fits
from astropy.table import Column, Table


ROOT = Path(__file__).resolve().parent


def uniform_lmst(ra):
    """Circular RA ranks on a uniform grid, rotated to mean signed HA=0.

    Cut at the largest empty RA gap before unwrapping. Equal RA values keep
    input order. Rotation preserves the uniform spacing of 360/N degrees.
    """
    ra = np.asarray(ra, dtype=float)
    if not np.all(np.isfinite(ra)):
        raise ValueError('RA values must be finite')
    ra = ra % 360
    if len(ra) == 0:
        return ra.copy()
    sorted_ra = np.sort(ra)
    gaps = np.diff(np.r_[sorted_ra, sorted_ra[0] + 360])
    gap = int(np.argmax(gaps))
    cut = (sorted_ra[gap] + gaps[gap] / 2) % 360
    unwrapped = (ra - cut) % 360
    order = np.argsort(unwrapped, kind='stable')
    grid = (np.arange(len(ra)) + 0.5) * 360 / len(ra)
    grid += np.mean(unwrapped[order] - grid)
    result = np.empty(len(ra))
    result[order] = (grid + cut) % 360
    # With highly concentrated late-survey samples, choose the local phase
    # with zero mean wrapped HA rather than allowing a branch-cut bias.
    for _ in range(100):
        mean_ha = np.mean((result - ra + 180) % 360 - 180)
        if abs(mean_ha) < 1e-10:
            return result
        result = (result - mean_ha) % 360
    raise ValueError('Could not find a zero-mean hour-angle phase')


def design_values(table):
    """Assign each filter independently; excluded/completed rows get NaN."""
    result = np.full(len(table), np.nan)
    eligible = ((table['IN_IBIS'] == 1) & (table['IN_HSC'] == 0)
                & (table['DONE'] == 0))
    for band in np.unique(table['FILTER']):
        indices = np.flatnonzero(eligible & (table['FILTER'] == band))
        result[indices] = uniform_lmst(table['RA'][indices])
    return result


def update(root=ROOT):
    ecsv_path = root / 'obstatus/cfht-tiles.ecsv'
    fits_path = root / 'obstatus/cfht-tiles.fits'
    ecsv = Table.read(ecsv_path)
    original_fits = Table.read(fits_path)
    for name in ('OBJECT', 'RA', 'DEC', 'FILTER', 'IN_IBIS', 'IN_HSC'):
        np.testing.assert_array_equal(ecsv[name], original_fits[name])
    values = design_values(ecsv)
    original_columns = [c for c in ecsv.colnames if c != 'LMST_DESIGN']
    ecsv['LMST_DESIGN'] = Column(values, unit='deg', description=(
        'Uniform circular RA ranks per filter, mean signed HA=0; '
        'remaining IN_IBIS=1, IN_HSC=0, DONE=0 rows in ECSV; otherwise NaN'))

    with tempfile.TemporaryDirectory(prefix='cfht-lmst-', dir=root) as tmp:
        new_ecsv, new_fits = Path(tmp) / ecsv_path.name, Path(tmp) / fits_path.name
        ecsv.write(new_ecsv, format='ascii.ecsv')
        with fits.open(fits_path) as hdus:
            columns = fits.ColDefs([c for c in hdus[1].columns if c.name != 'LMST_DESIGN'])
            columns += fits.ColDefs([fits.Column(name='LMST_DESIGN', format='D',
                                                  unit='deg', array=values)])
            hdus[1] = fits.BinTableHDU.from_columns(columns, header=hdus[1].header)
            hdus[1].header['LMDMETH'] = ('RA rank, mean HA=0', 'Independent per filter')
            hdus[1].header['LMDREF'] = ('ECSV', 'Eligibility and DONE source')
            hdus.writeto(new_fits)
        for path, original in ((new_ecsv, ecsv), (new_fits, original_fits)):
            saved = Table.read(path)
            for name in original_columns:
                np.testing.assert_array_equal(saved[name], original[name])
            np.testing.assert_array_equal(saved['LMST_DESIGN'], values)
        new_ecsv.replace(ecsv_path)
        new_fits.replace(fits_path)
    for band in np.unique(ecsv['FILTER']):
        selected = (ecsv['FILTER'] == band) & np.isfinite(values)
        ha = (values[selected] - ecsv['RA'][selected] + 180) % 360 - 180
        if len(ha):
            print(f'{band}: {len(ha)} designs; mean HA={np.mean(ha):.9f} deg; '
                  f'HA range {np.min(ha):.3f} to {np.max(ha):.3f} deg')
    print(f'Undefined/excluded rows: {np.count_nonzero(~np.isfinite(values))}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    update(parser.parse_args().root)
