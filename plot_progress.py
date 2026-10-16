#!/usr/bin/env python3
"""Plot the footprint observed to date in one filter and report the completed fraction.

Requires numpy and matplotlib. Reads obstatus/cfht-tiles.ecsv (no astropy needed).
The completed fraction is the fraction of tiles: DONE!=0 tiles divided by the
tiles the program is responsible for in that filter (IN_IBIS=1 and IN_HSC=0,
plus any DONE tiles that fall inside the IN_HSC region so that every
observation counts). Counts for the full IN_IBIS footprint are also printed.
Sky areas are reported for reference only, treating each tile as a square
field of --field-size degrees (default 1.0) rasterized on a fine grid.

Example:
    python plot_progress.py --filter M4376
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from night_outputs import PLOT_DEC_RANGE, PLOT_RA_LEFT, PLOT_RA_RIGHT, plot_ra

DEFAULT_INPUT = Path('obstatus/cfht-tiles.ecsv')
DEFAULT_OUTDIR = Path('obstatus')
COLUMNS = ('OBJECT', 'RA', 'DEC', 'FILTER', 'IN_IBIS', 'IN_HSC', 'DONE', 'PROGRAM', 'DATE')
PROGRAM_COLORS = {'LBNL': 'tab:blue', 'NAOC': 'tab:orange'}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--input', type=Path, default=DEFAULT_INPUT)
    parser.add_argument('--filter', dest='filter_name', default='M4376')
    parser.add_argument('--outdir', type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument('--field-size', type=float, default=1.0, metavar='DEG',
                        help='Side of the square field used for area coverage (default 1.0).')
    parser.add_argument('--cell', type=float, default=0.05, metavar='DEG',
                        help='Raster cell size for the area computation (default 0.05).')
    args = parser.parse_args()
    args.filter_name = args.filter_name.upper()
    if not 0 < args.field_size <= 5 or not 0 < args.cell <= args.field_size / 4:
        parser.error('--field-size must be in (0, 5] and --cell at most a quarter of it')
    return args


def read_tiles(path: Path, filter_name: str) -> dict[str, np.ndarray]:
    values = {name: [] for name in COLUMNS}
    with path.open(newline='', encoding='utf-8') as handle:
        rows = csv.reader((line for line in handle if not line.startswith('#')),
                          delimiter=' ', skipinitialspace=True)
        header = next(rows)
        index = {name: header.index(name) for name in COLUMNS}
        for parts in rows:
            if parts and parts[index['FILTER']].upper() == filter_name:
                for name in COLUMNS:
                    values[name].append(parts[index[name]])
    out = {name: np.asarray(values[name]) for name in ('OBJECT', 'FILTER', 'PROGRAM', 'DATE')}
    out['RA'] = np.asarray(values['RA'], dtype=float)
    out['DEC'] = np.asarray(values['DEC'], dtype=float)
    for name in ('IN_IBIS', 'IN_HSC', 'DONE'):
        out[name] = np.asarray(values[name], dtype=float).astype(int)
    out['PROGRAM'] = np.char.upper(out['PROGRAM'])
    return out


class Raster:
    """Boolean coverage map on an RA/DEC grid; RA is handled modulo 360."""

    def __init__(self, dec_range: tuple[float, float], cell: float):
        self.cell = cell
        self.n_ra = int(round(360.0 / cell))
        self.dec0 = dec_range[0]
        self.n_dec = int(np.ceil((dec_range[1] - dec_range[0]) / cell))
        dec_centers = self.dec0 + (np.arange(self.n_dec) + 0.5) * cell
        # Solid angle of each cell in square degrees.
        self.cell_area = np.cos(np.deg2rad(dec_centers)) * cell * cell

    def union(self, ra: np.ndarray, dec: np.ndarray, size: float) -> np.ndarray:
        covered = np.zeros((self.n_dec, self.n_ra), dtype=bool)
        half = size / 2.0
        for ra_c, dec_c in zip(ra, dec):
            d0 = int(np.floor((dec_c - half - self.dec0) / self.cell))
            d1 = int(np.ceil((dec_c + half - self.dec0) / self.cell))
            d0, d1 = max(d0, 0), min(d1, self.n_dec)
            if d1 <= d0:
                continue
            half_ra = half / max(np.cos(np.deg2rad(dec_c)), 1e-6)
            r0 = int(np.floor((ra_c - half_ra) / self.cell))
            r1 = int(np.ceil((ra_c + half_ra) / self.cell))
            cols = np.arange(r0, r1) % self.n_ra
            covered[d0:d1, cols] = True
        return covered

    def area(self, covered: np.ndarray) -> float:
        return float((covered.sum(axis=1) * self.cell_area).sum())


def main() -> None:
    args = parse_args()
    tiles = read_tiles(args.input, args.filter_name)
    footprint = tiles['IN_IBIS'] == 1
    done = footprint & (tiles['DONE'] != 0)
    hsc = footprint & (tiles['IN_HSC'] == 1)
    if not footprint.any():
        raise SystemExit(f'No IN_IBIS=1 tiles in filter {args.filter_name}')

    dec_lo = float(tiles['DEC'][footprint].min() - args.field_size)
    dec_hi = float(tiles['DEC'][footprint].max() + args.field_size)
    raster = Raster((dec_lo, dec_hi), args.cell)
    cover_all = raster.union(tiles['RA'][footprint], tiles['DEC'][footprint], args.field_size)
    cover_hsc = raster.union(tiles['RA'][hsc], tiles['DEC'][hsc], args.field_size)
    cover_done = raster.union(tiles['RA'][done], tiles['DEC'][done], args.field_size)
    cover_cfht = cover_all & ~cover_hsc

    area_all = raster.area(cover_all)
    area_cfht = raster.area(cover_cfht)
    area_done = raster.area(cover_done & cover_all)
    area_done_cfht = raster.area(cover_done & cover_cfht)
    dates = sorted(set(tiles['DATE'][done]))

    responsible = footprint & (~hsc | done)
    n_total, n_done = int(responsible.sum()), int(done.sum())
    fraction = n_done / n_total
    print(f'Filter {args.filter_name}: {n_done} of {n_total} tiles observed ({100 * fraction:.1f}%)'
          f' [IN_IBIS=1 and IN_HSC=0, plus {int((done & hsc).sum())} DONE tiles inside IN_HSC]')
    print(f'  Full IN_IBIS footprint: {n_done} of {int(footprint.sum())} tiles '
          f'({100 * n_done / footprint.sum():.1f}%); {int(hsc.sum())} tiles are IN_HSC')
    for program in sorted(set(tiles['PROGRAM'][footprint])):
        sel = footprint & (tiles['PROGRAM'] == program)
        print(f'  PROGRAM {program or "(blank)"}: {sel.sum()} tiles ({(sel & ~hsc).sum()} outside IN_HSC), '
              f'{(sel & done).sum()} DONE')
    print(f'Sky area for reference (1 tile = {args.field_size:g} deg square):')
    print(f'  Full IN_IBIS footprint:  {area_all:8.1f} sq deg; observed {area_done:7.1f} sq deg '
          f'({100 * area_done / area_all:.1f}%)')
    print(f'  Excluding IN_HSC region: {area_cfht:8.1f} sq deg; observed {area_done_cfht:7.1f} sq deg '
          f'({100 * area_done_cfht / area_cfht:.1f}%)')
    print(f'Observations on {len(dates)} nights: {dates[0]} to {dates[-1]}' if dates else 'No observations')

    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    fig, ax = plt.subplots(figsize=(14, 3.8))
    x_all = plot_ra(tiles['RA'])
    area_hsc = raster.area(cover_hsc & cover_all)
    for mask, color, label in ((footprint & ~hsc, '0.85', f'Footprint: {int((footprint & ~hsc).sum())} tiles'),
                               (hsc, '#c9d8ea', f'IN_HSC: {int(hsc.sum())} tiles')):
        ax.scatter(x_all[mask], tiles['DEC'][mask], s=4, c=color, marker='s', linewidths=0,
                   rasterized=True)
        ax.scatter([], [], marker='s', s=40, c=color, label=label)
    for program in sorted(set(tiles['PROGRAM'][done])):
        sel = done & (tiles['PROGRAM'] == program)
        color = PROGRAM_COLORS.get(program, 'tab:green')
        for ra_c, dec_c in zip(x_all[sel], tiles['DEC'][sel]):
            half_ra = args.field_size / 2 / np.cos(np.deg2rad(dec_c))
            ax.add_patch(Rectangle((ra_c - half_ra, dec_c - args.field_size / 2), 2 * half_ra,
                                   args.field_size, facecolor=color, edgecolor='none', alpha=0.6))
        ax.scatter([], [], marker='s', s=40, c=color, label=f'DONE {program or "(blank)"}: {sel.sum()} tiles')

    ax.set_xlim(PLOT_RA_LEFT, PLOT_RA_RIGHT)
    ax.set_ylim(*PLOT_DEC_RANGE)
    ticks = np.arange(PLOT_RA_LEFT, PLOT_RA_RIGHT - 1, -30.0)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{int(t % 360)}' for t in ticks])
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('DEC (deg)')
    ax.set_aspect('equal')
    ax.grid(True, linewidth=0.3, alpha=0.5)
    ax.set_title(f'{args.filter_name} observed to date ({dates[-1] if dates else "none"}): '
                 f'{n_done} of {n_total} tiles ({100 * fraction:.1f}%); '
                 f'about {area_done_cfht:.0f} of {area_cfht:.0f} sq deg outside IN_HSC',
                 fontsize=10)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.38), ncol=4, fontsize=7, frameon=False)
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / f'progress-{args.filter_name}.pdf'
    fig.savefig(out, format='pdf', dpi=200, bbox_inches='tight')
    print(f'Wrote: {out}')


if __name__ == '__main__':
    main()
