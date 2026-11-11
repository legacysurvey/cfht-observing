#!/usr/bin/env python
"""Cross-check the Kealahou queues of both programs against the tile catalog (read-only).

Fetches targets, observing groups (OGs) and exposures for NAOC (26BZ01) and
LBNL (26BZ50) and reports, without changing anything:
  * active OGs whose tile is labelled for the other program, is outside the
    CFHT program (IN_IBIS=0 or IN_HSC=1), or is already DONE;
  * tiles active in both programs' queues;
  * unobserved tiles labelled for a program with no active OG there;
  * exposures in Kealahou whose tile is not DONE in the catalog (missed by
    end-of-night.py), and DONE tiles with no exposure in either program;
  * active OGs of targets that are not survey tiles (e.g. COSMOS).
With --write DIR the name lists behind each finding are written as text files.

Example:
    python check-queues.py --write queue-check
"""
import argparse
import os
import sys
from collections import Counter, defaultdict

import numpy as np
from astropy.table import Table

from cfht_common import alias_run_ids
from kealahou_api import KealahouProgram

SNAP_COMMENT = 'MEGACAM SNR SNAP'


def load_program(program, token=None):
    k = KealahouProgram(program, token=token)
    targets = k.fetch_targets()
    token_to_name = {t['token']: t['name'] for t in targets}
    og_state = defaultdict(list)          # target name -> [(og token, priority)]
    for og in k.fetch_ogs():
        comps = (((og.get('single_observing_group') or {}).get('observing_block') or {})
                 .get('observing_component') or [])
        for comp in comps:
            name = token_to_name.get(comp.get('target_token'))
            if name is not None:
                og_state[name].append((og.get('token'), og.get('og_priority')))
    exposures = defaultdict(list)         # target name -> [obsid]
    for e in k.fetch_exposures():
        name = (e.get('target_data') or {}).get('name')
        status = e.get('exposure_status') or {}
        if name is None or status.get('obs_comment') == SNAP_COMMENT:
            continue
        exposures[name].append(int(e.get('obsid') or 0))
    return dict(targets=set(token_to_name.values()), ogs=og_state, exposures=exposures)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--token', help='Kealahou API token; default is $KEALAHOU_TOKEN')
    parser.add_argument('--tile-file', default='obstatus/cfht-tiles.ecsv')
    parser.add_argument('--write', metavar='DIR', help='Write the name lists of each finding into DIR')
    args = parser.parse_args()

    tiles = Table.read(args.tile_file)
    label = np.array(tiles['PROGRAM'].filled('') if hasattr(tiles['PROGRAM'], 'filled') else tiles['PROGRAM']).astype(str)
    names = np.char.strip(np.array(tiles['OBJECT'], dtype=str))
    info = {n: (l, bool(d), bool(i == 1 and h == 0)) for n, l, d, i, h in
            zip(names, label, np.array(tiles['DONE']) != 0, tiles['IN_IBIS'], tiles['IN_HSC'])}

    programs = list(alias_run_ids)        # NAOC, LBNL
    try:
        queue = {p: load_program(p, args.token) for p in programs}
    except RuntimeError:
        return -1
    findings = {}

    def report(key, text, items):
        items = sorted(items)
        findings[key] = items
        flag = 'ok ' if not items else '!! '
        print(f'{flag}{text}: {len(items)}' + (f'   e.g. {", ".join(items[:4])}' if items else ''))

    active = {p: {n for n, ogs in queue[p]['ogs'].items() if any(pr != 'INACTIVE' for _, pr in ogs)} for p in programs}
    for p in programs:
        q = queue[p]
        other = [x for x in programs if x != p][0]
        prio = Counter(pr for ogs in q['ogs'].values() for _, pr in ogs)
        print(f'\n=== {p} ({alias_run_ids[p]}): {len(q["targets"])} targets, {sum(prio.values())} OGs {dict(prio)}, '
              f'{sum(len(v) for v in q["exposures"].values())} science exposures on {len(q["exposures"])} targets ===')
        survey = {n for n in active[p] if n in info}
        report(f'{p}-active-nonsurvey', 'active OGs of targets that are not survey tiles (expected for COSMOS)', active[p] - survey)
        report(f'{p}-active-labelled-{other}', f'active OGs whose tile is labelled {other}', {n for n in survey if info[n][0] == other})
        report(f'{p}-active-unlabelled', 'active OGs whose tile has no program label', {n for n in survey if info[n][0] == ''})
        report(f'{p}-active-outside-footprint', 'active OGs whose tile is outside the CFHT program', {n for n in survey if not info[n][2]})
        report(f'{p}-active-done', 'active OGs whose tile is already DONE', {n for n in survey if info[n][1]})
        mine = {n for n, (l, d, inprog) in info.items() if l == p and not d and inprog}
        report(f'{p}-unobserved-not-active', f'unobserved {p} tiles with no active OG in {p}'
               + (' (expected: most LBNL tiles are not loaded)' if p == 'LBNL' else ''), mine - active[p])
        report(f'{p}-exposed-not-done', 'tiles with an exposure here but DONE=0 in the catalog (missed by end-of-night)',
               {n for n in q['exposures'] if n in info and not info[n][1]})
        report(f'{p}-exposed-labelled-{other}', f'tiles exposed here but labelled {other}',
               {n for n in q['exposures'] if n in info and info[n][0] == other})
        report(f'{p}-multiple-exposures', 'tiles with more than one science exposure here',
               {n for n, v in q['exposures'].items() if n in info and len(v) > 1})
        report(f'{p}-multiple-ogs', 'targets with more than one OG', {n for n, v in q['ogs'].items() if len(v) > 1})
    print('\n=== both programs ===')
    a, b = programs
    report('active-in-both', 'tiles with an active OG in both programs', active[a] & active[b])
    report('target-in-both', 'tiles loaded as targets in both programs', (queue[a]['targets'] & queue[b]['targets']) & set(info))
    exposed = set(queue[a]['exposures']) | set(queue[b]['exposures'])
    report('done-without-exposure', 'DONE tiles with no science exposure in either program', {n for n, (l, d, i) in info.items() if d and n not in exposed})
    report('exposed-in-both', 'tiles with exposures in both programs', set(queue[a]['exposures']) & set(queue[b]['exposures']) & set(info))
    if args.write:
        os.makedirs(args.write, exist_ok=True)
        for key, items in findings.items():
            if items:
                with open(os.path.join(args.write, key + '.names'), 'w') as f:
                    f.write('\n'.join(items) + '\n')
        print('\nWrote name lists to', args.write)
    return 0


if __name__ == '__main__':
    sys.exit(main())
