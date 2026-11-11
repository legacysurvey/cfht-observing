#!/usr/bin/env python3
"""Make a small dithered set of targets on one field, outside the survey tiling.

Writes one plan per filter (FITS and XML, the formats upload-plan.py and the
Kealahou target import use) into plans/special/, which the tile updaters do
not read: these targets are not survey tiles and are not in cfht-tiles.ecsv.
Targets are named PREFIX_FILTER_N. The default pattern is four pointings
stepped diagonally by 25 arcsec in RA and 100 arcsec in DEC, which exceeds the
MegaCam gaps (up to 13 arcsec between CCD columns and between the two middle
rows, 87 arcsec between the outer and middle rows), so no sky in the central
field falls in a gap in more than two of the four exposures.

MAG_AB is the survey depth for the field's extinction:
DEFAULT_MAG_AB[filter] + coefficient * E(B-V), or --mag-ab + coefficient * E(B-V).

Example:
    python make-dither-targets.py --prefix COSMOS --ra 150.119167 --dec 2.205833 --ebv 0.019
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from make_cfht_targets import (DEFAULT_MAG_AB, DEFAULT_TEMPLATE, EXTINCTION_COEFFICIENT, OutputTarget,
                               format_dec, format_ra, write_fits, write_xml)

DEFAULT_OFFSETS = [(-37.5, -150.0), (-12.5, -50.0), (12.5, 50.0), (37.5, 150.0)]  # arcsec (east, north)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--prefix', required=True, help='Target name prefix, e.g. COSMOS')
    parser.add_argument('--ra', type=float, required=True, help='Field center RA (deg, J2000)')
    parser.add_argument('--dec', type=float, required=True, help='Field center DEC (deg, J2000)')
    parser.add_argument('--ebv', type=float, required=True, help='E(B-V) of the field (SFD)')
    parser.add_argument('--filters', nargs='+', default=['M4112', 'M4376'])
    parser.add_argument('--offsets', type=float, nargs='+', metavar='ARCSEC',
                        help='Dither offsets as east north pairs in arcsec (default: 4-point diagonal 25 x 100)')
    parser.add_argument('--mag-ab', type=float, help='Base AB magnitude before extinction (default per filter)')
    parser.add_argument('--outdir', type=Path, default=Path('plans/special'))
    parser.add_argument('--template', type=Path, default=DEFAULT_TEMPLATE)
    args = parser.parse_args()
    if args.prefix.upper().startswith('DESI'):
        parser.error('--prefix must not look like a survey tile name')
    if args.offsets:
        if len(args.offsets) % 2:
            parser.error('--offsets needs east north pairs')
        offsets = list(zip(args.offsets[0::2], args.offsets[1::2]))
    else:
        offsets = DEFAULT_OFFSETS
    args.outdir.mkdir(parents=True, exist_ok=True)
    cosd = math.cos(math.radians(args.dec))
    for filter_name in (f.upper() for f in args.filters):
        if filter_name not in EXTINCTION_COEFFICIENT:
            parser.error(f'no extinction coefficient for filter {filter_name}')
        base = DEFAULT_MAG_AB[filter_name] if args.mag_ab is None else args.mag_ab
        mag = base + EXTINCTION_COEFFICIENT[filter_name] * args.ebv
        rows = []
        for i, (east, north) in enumerate(offsets, start=1):
            ra = (args.ra + east / 3600.0 / cosd) % 360.0
            dec = args.dec + north / 3600.0
            rows.append(OutputTarget(name=f'{args.prefix}_{filter_name}_{i}', ra_j2000=format_ra(ra),
                                     dec_j2000=format_dec(dec), ra_deg=ra, dec_deg=dec, mag_ab=mag))
        stem = args.outdir / f'targets-{args.prefix}-{filter_name}'
        write_xml(args.template, stem.with_suffix('.xml'), rows)
        write_fits(stem.with_suffix('.fits'), rows)
        print(f'{filter_name}: MAG_AB {mag:.2f}')
        for row in rows:
            print(f'   {row.name}  {row.ra_j2000} {row.dec_j2000}  ({row.ra_deg:.5f}, {row.dec_deg:.5f})')
        print(f'Wrote: {stem}.xml and {stem}.fits')


if __name__ == '__main__':
    main()
