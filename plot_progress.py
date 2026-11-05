#!/usr/bin/env python3
"""Plot the tiles observed to date, one panel per filter, with the fraction of tiles done.

Requires numpy and matplotlib. Reads obstatus/cfht-tiles.ecsv (no astropy needed).
For each filter the completed fraction is DONE!=0 tiles divided by the tiles
the program is responsible for: IN_IBIS=1 and IN_HSC=0, plus any DONE tiles
that fall inside the IN_HSC region so that every observation counts.

Example:
    python plot_progress.py
    python plot_progress.py --filters M4376 M4112 --outdir obstatus
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from night_outputs import PLOT_DEC_RANGE, PLOT_RA_LEFT, PLOT_RA_RIGHT, plot_ra

DEFAULT_INPUT = Path('obstatus/cfht-tiles.ecsv')
DEFAULT_OUTDIR = Path('obstatus')
DEFAULT_FILTERS = ('M4376', 'M4112')
COLUMNS = ('OBJECT', 'RA', 'DEC', 'FILTER', 'IN_IBIS', 'IN_HSC', 'DONE', 'PROGRAM', 'DATE')
PROGRAM_COLORS = {'LBNL': 'tab:blue', 'NAOC': 'tab:orange'}
QUEUE_COLORS = {'observed': '#1f4e9c', 'pending': '#ff7f0e'}  # dark blue vs orange
FIELD_SIZE_DEG = 1.0  # drawn size of an observed tile


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--input', type=Path, default=DEFAULT_INPUT)
    parser.add_argument('--filters', nargs='+', default=list(DEFAULT_FILTERS), metavar='FILTER',
                        help='Filters to plot, one panel each, top to bottom (default: M4376 M4112).')
    parser.add_argument('--outdir', type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument('--queue', metavar='PROGRAM',
                        help='Instead of all DONE tiles, show every tile of this PROGRAM (e.g. NAOC): '
                             'observed in dark blue, not yet observed in orange. Writes progress-PROGRAM.pdf.')
    args = parser.parse_args()
    args.filters = [name.upper() for name in args.filters]
    if args.queue:
        args.queue = args.queue.upper()
    return args


def draw_tiles(ax, x, dec, facecolor, edgecolor, alpha=0.6, linewidth=0.0):
    from matplotlib.patches import Rectangle
    for ra_c, dec_c in zip(x, dec):
        half_ra = FIELD_SIZE_DEG / 2 / np.cos(np.deg2rad(dec_c))
        ax.add_patch(Rectangle((ra_c - half_ra, dec_c - FIELD_SIZE_DEG / 2), 2 * half_ra, FIELD_SIZE_DEG,
                               facecolor=facecolor, edgecolor=edgecolor, alpha=alpha,
                               linewidth=linewidth))


def read_tiles(path: Path) -> dict[str, np.ndarray]:
    values = {name: [] for name in COLUMNS}
    with path.open(newline='', encoding='utf-8') as handle:
        rows = csv.reader((line for line in handle if not line.startswith('#')),
                          delimiter=' ', skipinitialspace=True)
        header = next(rows)
        index = {name: header.index(name) for name in COLUMNS}
        for parts in rows:
            if parts:
                for name in COLUMNS:
                    values[name].append(parts[index[name]])
    out = {name: np.asarray(values[name]) for name in ('OBJECT', 'DATE')}
    out['FILTER'] = np.char.upper(np.asarray(values['FILTER']))
    out['PROGRAM'] = np.char.upper(np.asarray(values['PROGRAM']))
    out['RA'] = np.asarray(values['RA'], dtype=float)
    out['DEC'] = np.asarray(values['DEC'], dtype=float)
    for name in ('IN_IBIS', 'IN_HSC', 'DONE'):
        out[name] = np.asarray(values[name], dtype=float).astype(int)
    return out


def filter_summary(tiles: dict[str, np.ndarray], filter_name: str) -> dict:
    in_filter = tiles['FILTER'] == filter_name
    footprint = in_filter & (tiles['IN_IBIS'] == 1)
    hsc = footprint & (tiles['IN_HSC'] == 1)
    done = in_filter & (tiles['DONE'] != 0)          # every exposure counts, wherever the tile sits now
    responsible = (footprint & ~hsc) | done
    dates = sorted(set(tiles['DATE'][done]))
    return dict(footprint=footprint, hsc=hsc, done=done, responsible=responsible,
                n_total=int(responsible.sum()), n_done=int(done.sum()), dates=dates)


def main() -> None:
    args = parse_args()
    tiles = read_tiles(args.input)
    summaries = {name: filter_summary(tiles, name) for name in args.filters}
    for name, info in summaries.items():
        if not info['footprint'].any():
            raise SystemExit(f'No IN_IBIS=1 tiles in filter {name}')
        fraction = info['n_done'] / info['n_total']
        outside = int((info["done"] & ~(info["footprint"] & ~info["hsc"])).sum())
        print(f'Filter {name}: {info["n_done"]} of {info["n_total"]} tiles observed ({100 * fraction:.1f}%)'
              f' [program tiles IN_IBIS=1 and IN_HSC=0, plus {outside} observed tiles now outside the program;'
              f' {int(info["footprint"].sum())} IN_IBIS tiles in all]')
        for program in sorted(set(tiles['PROGRAM'][info['footprint']])):
            sel = info['footprint'] & (tiles['PROGRAM'] == program)
            print(f'  PROGRAM {program or "(blank)"}: {sel.sum()} tiles '
                  f'({(sel & ~info["hsc"]).sum()} outside IN_HSC), {(sel & info["done"]).sum()} DONE')
        if info['dates']:
            print(f'  {len(info["dates"])} nights: {info["dates"][0]} to {info["dates"][-1]}')
        else:
            print('  No observations')

    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt

    # Each panel spans 360 x 50 degrees at equal aspect: ~1.9 in tall at 14 in wide.
    fig, axes = plt.subplots(len(args.filters), 1, figsize=(14, 2.7 * len(args.filters)),
                             sharex=True, squeeze=False, constrained_layout=True)
    x_all = plot_ra(tiles['RA'])
    for ax, (name, info) in zip(axes[:, 0], summaries.items()):
        footprint, hsc, done = info['footprint'], info['hsc'], info['done']
        for mask, color, label in ((footprint & ~hsc, '0.85', f'Footprint: {int((footprint & ~hsc).sum())} tiles'),
                                   (hsc, '#c9d8ea', f'IN_HSC: {int(hsc.sum())} tiles')):
            ax.scatter(x_all[mask], tiles['DEC'][mask], s=4, c=color, marker='s', linewidths=0,
                       rasterized=True)
            ax.scatter([], [], marker='s', s=40, c=color, label=label)
        if args.queue:
            queue = footprint & (tiles['PROGRAM'] == args.queue)
            pending, observed = queue & ~done, queue & done
            draw_tiles(ax, x_all[pending], tiles['DEC'][pending], QUEUE_COLORS['pending'], 'none', alpha=0.8)
            draw_tiles(ax, x_all[observed], tiles['DEC'][observed], QUEUE_COLORS['observed'], 'none', alpha=0.9)
            ax.scatter([], [], marker='s', s=40, c=QUEUE_COLORS['pending'],
                       label=f'{args.queue} not yet observed: {int(pending.sum())} tiles')
            ax.scatter([], [], marker='s', s=40, c=QUEUE_COLORS['observed'],
                       label=f'{args.queue} observed: {int(observed.sum())} tiles')
            queue_dates = sorted(set(tiles['DATE'][observed]))
            latest = queue_dates[-1] if queue_dates else 'no observations'
            fraction = observed.sum() / queue.sum() if queue.sum() else float('nan')
            ax.set_title(f'{name}: {args.queue} queue, {int(observed.sum())} of {int(queue.sum())} tiles observed '
                         f'({100 * fraction:.1f}%) through {latest}', fontsize=10)
        else:
            for program in sorted(set(tiles['PROGRAM'][done])):
                sel = done & (tiles['PROGRAM'] == program)
                color = PROGRAM_COLORS.get(program, 'tab:green')
                draw_tiles(ax, x_all[sel], tiles['DEC'][sel], color, 'none')
                ax.scatter([], [], marker='s', s=40, c=color, label=f'DONE {program or "(blank)"}: {sel.sum()} tiles')
            fraction = info['n_done'] / info['n_total']
            latest = info['dates'][-1] if info['dates'] else 'no observations'
            ax.set_title(f'{name}: {info["n_done"]} of {info["n_total"]} tiles observed '
                         f'({100 * fraction:.1f}%) through {latest}', fontsize=10)
        ax.set_xlim(PLOT_RA_LEFT, PLOT_RA_RIGHT)
        ax.set_ylim(*PLOT_DEC_RANGE)
        ax.set_ylabel('DEC (deg)')
        ax.set_aspect('equal')
        ax.grid(True, linewidth=0.3, alpha=0.5)
        # The strip above DEC +20 is empty in every filter, so the legend fits inside the axes.
        ax.legend(loc='upper center', ncol=4, fontsize=7, framealpha=0.9, borderpad=0.3)
    ticks = np.arange(PLOT_RA_LEFT, PLOT_RA_RIGHT - 1, -30.0)
    axes[-1, 0].set_xticks(ticks)
    axes[-1, 0].set_xticklabels([f'{int(t % 360)}' for t in ticks])
    axes[-1, 0].set_xlabel('RA (deg)')
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / (f'progress-{args.queue}.pdf' if args.queue else 'progress.pdf')
    fig.savefig(out, format='pdf', dpi=200, bbox_inches='tight')
    print(f'Wrote: {out}')


if __name__ == '__main__':
    main()
