#!/usr/bin/env python
"""Update the requested magnitude (and optionally the template) of targets already in Kealahou.

For each named survey tile that is not yet observed (DONE=0 in the tile file),
set the Kealahou target's AB magnitude to the survey depth for its extinction,
DEFAULT_MAG_AB[filter] + coefficient * EBV_MED, as make_cfht_targets.py writes
it. Observed tiles are left alone. upload-plan.py does not do this: it skips
targets that already exist. With --ot-name, the OGs of these targets are also
switched to that observing template; OG priorities are never changed.

Example:
    python update-targets.py --run-id NAOC --names-file plans/special/naoc-august-unobserved.names --dry-run
    python update-targets.py --run-id NAOC --names-file plans/special/naoc-august-unobserved.names
"""
import argparse
import sys

import numpy as np
from astropy.table import Table

from kealahou_api import KealahouProgram
from make_cfht_targets import DEFAULT_MAG_AB, EXTINCTION_COEFFICIENT


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--token', help='Kealahou API token; default is $KEALAHOU_TOKEN')
    parser.add_argument('--run-id', default='NAOC', help='NAOC (26BZ01), LBNL (26BZ50), or an ID')
    parser.add_argument('--tile-file', default='obstatus/cfht-tiles.ecsv')
    parser.add_argument('--names-file', help='File of target names, one per line (# comments allowed)')
    parser.add_argument('--ot-name', help='Also switch the OGs of these targets to this observing template')
    parser.add_argument('--dry-run', action='store_true', help='Only report what would change')
    parser.add_argument('names', nargs='*')
    args = parser.parse_args()
    names = list(args.names)
    if args.names_file:
        with open(args.names_file) as f:
            names += [line.split('#', 1)[0].split()[0] for line in f if line.split('#', 1)[0].strip()]
    if not names:
        parser.error('no target names given')

    tiles = Table.read(args.tile_file)
    index = {str(n).strip(): i for i, n in enumerate(tiles['OBJECT'])}
    try:
        kealahou = KealahouProgram(args.run_id, token=args.token)
    except RuntimeError:
        return -1
    targets = {t['name']: t for t in kealahou.fetch_targets()}

    changed = skipped = failed = 0
    selected = []
    for name in names:
        if name not in index:
            print(name, ': not a survey tile; skipped'); skipped += 1; continue
        if name not in targets:
            print(name, ': not in this Kealahou program; skipped'); skipped += 1; continue
        i = index[name]
        if tiles['DONE'][i] != 0:
            print(name, ': already observed; magnitude left alone'); skipped += 1; continue
        selected.append(name)
        filt = str(tiles['FILTER'][i]).upper()
        new = float(np.float32(DEFAULT_MAG_AB[filt] + EXTINCTION_COEFFICIENT[filt] * float(tiles['EBV_MED'][i])))
        target = targets[name]
        old = ((target.get('magnitude') or {}).get('AB') or {}).get('value')
        if old is not None and abs(old - new) < 0.005:
            print('%s : magnitude already %.2f' % (name, old)); continue
        print('%s : magnitude %s -> %.2f' % (name, 'none' if old is None else '%.2f' % old, new))
        if args.dry_run:
            changed += 1; continue
        target['magnitude'] = dict(AB=dict(value=new))
        try:
            kealahou.create_or_update_target(target)
            changed += 1
        except RuntimeError as exc:
            print('   failed:', str(exc)[-300:]); failed += 1

    if args.ot_name:
        ot_map = kealahou.fetch_ot_map()
        if args.ot_name not in ot_map:
            print('Error: unknown OT name "%s"; known: %s' % (args.ot_name, sorted(ot_map)))
            return -1
        want = ot_map[args.ot_name]
        token_to_name = {t['token']: n for n, t in targets.items()}
        for og in kealahou.fetch_ogs():
            comps = og['single_observing_group']['observing_block']['observing_component']
            hits = [c for c in comps if token_to_name.get(c['target_token']) in selected and c['observing_template_token'] != want]
            if not hits:
                continue
            print('%s : OG %s template -> %s' % (token_to_name[hits[0]['target_token']], og.get('label', '?'), args.ot_name))
            if args.dry_run:
                continue
            for c in hits:
                c['observing_template_token'] = want
            try:
                kealahou.create_or_update_og(og)
            except RuntimeError as exc:
                print('   failed:', str(exc)[-300:]); failed += 1
    print('%i magnitudes %s, %i skipped, %i failed' % (changed, 'would change' if args.dry_run else 'changed', skipped, failed))
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
