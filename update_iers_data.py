#!/usr/bin/env python3
"""Refresh project Earth-orientation and leap-second data from IERS providers.

Requires numpy and astropy. Downloads and validates both files before
replacing either local copy. No Python packages are installed or upgraded.
"""

from pathlib import Path
import tempfile
from urllib.request import urlopen

from astropy.utils import iers


DATA_DIRECTORY = Path(__file__).resolve().parent / 'data' / 'iers'
SOURCES = {
    'finals2000A.all': 'https://maia.usno.navy.mil/ser7/finals2000A.all',
    'Leap_Second.dat': 'https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat',
}


def update(directory=DATA_DIRECTORY):
    directory.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='iers-download-', dir=directory) as tmp:
        paths = {}
        for name, url in SOURCES.items():
            path = Path(tmp) / name
            with urlopen(url, timeout=60) as response:
                path.write_bytes(response.read())
            paths[name] = path
        table = iers.IERS_A.open(str(paths['finals2000A.all']))
        leaps = iers.LeapSeconds.open(str(paths['Leap_Second.dat']))
        if len(table) < 2 or len(leaps) < 1:
            raise ValueError('Downloaded IERS tables are empty')
        for name, path in paths.items():
            path.replace(directory / name)
    print(f'IERS coverage: MJD {table["MJD"][0]} to {table["MJD"][-1]}')
    print(f'Leap-second table expires: {leaps.expires.isot}')


if __name__ == '__main__':
    update()
