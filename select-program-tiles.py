#!/usr/bin/env python3
"""Pick non-overlapping program tiles inside given stripes and relabel them for a program.

Candidates are unobserved CFHT program tiles (IN_IBIS=1, IN_HSC=0, DONE=0) of
one filter, currently labelled --from (default LBNL), whose centers lie inside
the tabulated RA and DEC bounds of the named stripes (obstatus/desi2-stripes.ecsv).
They are trimmed greedily (highest PRIORITY first, then southern DEC, then RA)
so that no two selected fields overlap and none overlaps an existing tile of the
target program in that filter (observed or not). The chosen OBJECT names are
written to --out, and with --apply their PROGRAM is set in the catalog.

With --layer the selection is a further layer over sky that is already partly
observed: only the target program's unobserved tiles block overlaps, and
candidates are taken in order of how little of their field has been observed
so far in that filter (by any program), so the new layer fills the least
covered sky first.

Example:
    python select-program-tiles.py --stripes SGC-4 SGC-5 --filter M4376 --to NAOC \\
        --out plans/targets-2026-11-05-SGC45.names --apply
    python select-program-tiles.py --stripes SGC-6 SGC-7 --filter M4376 --to NAOC --layer \
        --out plans/targets-2026-11-11-SGC67.names --apply
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from footprint import read_ecsv_rows
from make_cfht_targets import DEFAULT_MIN_SEPARATION_DEG, TileTarget, trim_non_overlapping

ROOT = Path(__file__).resolve().parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', type=Path, default=ROOT / 'obstatus/cfht-tiles.ecsv')
    parser.add_argument('--stripes', nargs='+', required=True, metavar='STRIPE')
    parser.add_argument('--filter', dest='filter_name', default='M4376')
    parser.add_argument('--from', dest='from_program', default='LBNL')
    parser.add_argument('--to', dest='to_program', default='NAOC')
    parser.add_argument('--min-separation', type=float, default=DEFAULT_MIN_SEPARATION_DEG, metavar='DEG')
    parser.add_argument('--layer', action='store_true',
                        help='Add a layer over observed sky: overlap observed tiles freely, avoid only the '
                             "target program's unobserved tiles, and fill the least-covered sky first.")
    parser.add_argument('--field-size', type=float, default=0.96, metavar='DEG',
                        help='Field size used to measure existing coverage with --layer (default 0.96).')
    parser.add_argument('--out', type=Path, required=True, help='File receiving the selected OBJECT names')
    parser.add_argument('--apply', action='store_true', help='Relabel the selected tiles in the catalog')
    args = parser.parse_args()
    filter_name = args.filter_name.upper()

    stripes = {r['STRIPE']: r for r in read_ecsv_rows(ROOT / 'obstatus/desi2-stripes.ecsv')}
    unknown = set(args.stripes) - set(stripes)
    if unknown:
        raise SystemExit(f'Unknown stripes: {sorted(unknown)}')

    def inside(ra, dec):
        for name in args.stripes:
            s = stripes[name]
            ra0, ra1 = float(s['RA_MIN']), float(s['RA_MAX'])
            if ((ra - ra0) % 360.0) <= ((ra1 - ra0) % 360.0) and float(s['DEC_MIN']) <= dec <= float(s['DEC_MAX']):
                return True
        return False

    candidates, blocked, observed, header = [], [], [], None
    with args.input.open(newline='', encoding='utf-8') as handle:
        rows = csv.reader((line for line in handle if not line.startswith('#')), delimiter=' ', skipinitialspace=True)
        header = next(rows)
        col = {name: i for i, name in enumerate(header)}
        for p in rows:
            if not p or p[col['FILTER']].upper() != filter_name:
                continue
            ra, dec = float(p[col['RA']]), float(p[col['DEC']])
            program = p[col['PROGRAM']]
            is_done = float(p[col['DONE']]) != 0
            if is_done:
                observed.append((ra, dec))
            if program == args.to_program:
                if not (args.layer and is_done):
                    blocked.append((ra, dec))
                continue
            if (program == args.from_program and int(float(p[col['IN_IBIS']])) == 1
                    and int(float(p[col['IN_HSC']])) == 0 and float(p[col['DONE']]) == 0 and inside(ra, dec)):
                candidates.append(TileTarget(p[col['OBJECT']], ra, dec, priority=float(p[col['PRIORITY']])))
    print(f'{len(candidates)} candidate {args.from_program} tiles in {", ".join(args.stripes)} ({filter_name}); '
          f'{len(blocked)} {"unobserved " if args.layer else "existing "}{args.to_program} tiles block overlaps')
    if args.layer:
        covered = existing_coverage(candidates, observed, args.field_size)
        ordered = sorted(candidates, key=lambda t: (round(covered[t.object_name], 2), -t.priority,
                                                    t.dec_deg, t.ra_deg, t.object_name))
    else:
        ordered = sorted(candidates, key=lambda t: (-t.priority, t.dec_deg, t.ra_deg, t.object_name))
    chosen = trim_non_overlapping(ordered, args.min_separation, blocked)
    chosen.sort(key=lambda t: t.ra_deg)
    print(f'selected {len(chosen)} non-overlapping tiles')
    if args.layer and chosen:
        values = np.array([covered[t.object_name] for t in chosen])
        print(f'  mean number of existing exposures over the selected fields: {values.mean():.2f} '
              f'(min {values.min():.2f}, max {values.max():.2f})')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('w', encoding='utf-8') as handle:
        handle.write(f'# {len(chosen)} {filter_name} tiles from {", ".join(args.stripes)} relabelled '
                     f'{args.from_program} -> {args.to_program}; non-overlapping with each other and with '
                     f'existing {args.to_program} tiles (min separation {args.min_separation:g} deg)\n')
        for t in chosen:
            handle.write(f'{t.object_name} {t.ra_deg:.3f} {t.dec_deg:.3f}\n')
    print(f'Wrote: {args.out}')
    if not args.apply:
        return
    names = {t.object_name for t in chosen}
    lines, changed = [], 0
    with args.input.open(encoding='utf-8') as handle:
        for line in handle:
            if not line.startswith('#'):
                parts = line.rstrip('\n').split(' ')
                if parts[col['OBJECT']] in names:
                    assert parts[col['PROGRAM']] == args.from_program, line
                    parts[col['PROGRAM']] = args.to_program
                    line = ' '.join(parts) + '\n'
                    changed += 1
            lines.append(line)
    args.input.write_text(''.join(lines), encoding='utf-8')
    print(f'Relabelled {changed} tiles {args.from_program} -> {args.to_program} in {args.input}')


def existing_coverage(candidates, observed, field_size, cell=0.05):
    """Mean number of observed exposures over each candidate's field (flat-sky boxes)."""
    half = field_size / 2.0
    obs = np.asarray(observed, dtype=float).reshape(-1, 2)
    offsets = (np.arange(int(round(field_size / cell))) + 0.5) * cell - half
    result = {}
    for t in candidates:
        if len(obs) == 0:
            result[t.object_name] = 0.0
            continue
        cosd = max(np.cos(np.deg2rad(t.dec_deg)), 1e-6)
        dra = ((obs[:, 0] - t.ra_deg + 180.0) % 360.0 - 180.0) * cosd
        ddec = obs[:, 1] - t.dec_deg
        near = (np.abs(dra) < field_size) & (np.abs(ddec) < field_size)
        if not near.any():
            result[t.object_name] = 0.0
            continue
        gx, gy = np.meshgrid(offsets, offsets)
        count = np.zeros(gx.shape)
        for x, y in zip(dra[near], ddec[near]):
            count += (np.abs(gx - x) < half) & (np.abs(gy - y) < half)
        result[t.object_name] = float(count.mean())
    return result


if __name__ == '__main__':
    main()
