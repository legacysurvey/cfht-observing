#!/usr/bin/env python
"""Remove targets (and their observing groups) from a Kealahou program.

By default the targets are the tiles in the catalog that carry the program's
label but are no longer CFHT program tiles (IN_HSC=1) and have no exposure
(DONE=0); explicit OBJECT names can be given instead. For each one the OG(s)
referencing the target are deleted first, then the target itself. Targets
with exposures cannot be deleted by Kealahou (database integrity), which is
the intended safeguard. Use --dry-run to list what would be removed.
With --inactive nothing is deleted: the OGs are set to priority INACTIVE
instead, for targets that Kealahou refuses to delete.

Example:
    python remove-targets.py --run-id NAOC --dry-run
    python remove-targets.py --run-id NAOC
    python remove-targets.py --run-id NAOC DESI_M4376_259028 DESI_M4376_259188
    python remove-targets.py --run-id NAOC --inactive DESI_M4376_253234
"""
import argparse
import sys

from astropy.table import Table

from cfht_common import alias_run_ids
from kealahou_api import KealahouProgram


def catalog_targets(tile_file, program):
    tiles = Table.read(tile_file)
    label = tiles['PROGRAM'].filled('') if hasattr(tiles['PROGRAM'], 'filled') else tiles['PROGRAM']
    sel = (label == program) & (tiles['IN_HSC'] == 1) & (tiles['DONE'] == 0)
    return [str(name).strip() for name in tiles['OBJECT'][sel]]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--token', help='Kealahou API token; default is the $KEALAHOU_TOKEN environment variable')
    parser.add_argument('--run-id', default='NAOC', help='NAOC (26BZ01), LBNL (26BZ50), or an ID; default NAOC')
    parser.add_argument('--tile-file', default='obstatus/cfht-tiles.ecsv')
    parser.add_argument('--dry-run', action='store_true', help='Only report what would be deleted')
    parser.add_argument('--inactive', action='store_true',
                        help='Set the OGs to priority INACTIVE instead of deleting OGs and targets')
    parser.add_argument('--names-file', help='File of target names, one per line (# comments allowed)')
    parser.add_argument('names', nargs='*', help='Target names; default: catalog tiles with this PROGRAM, IN_HSC=1, DONE=0')
    args = parser.parse_args()

    program = args.run_id if args.run_id in alias_run_ids else next(
        (k for k, v in alias_run_ids.items() if v == args.run_id), args.run_id)
    names = list(args.names)
    if args.names_file:
        with open(args.names_file) as f:
            names += [line.split('#', 1)[0].split()[0] for line in f if line.split('#', 1)[0].strip()]
    names = names or catalog_targets(args.tile_file, program)
    print(f'{len(names)} target(s) to remove from {program}:')
    for name in names:
        print('  ', name)
    if not names:
        return 0

    try:
        kealahou = KealahouProgram(args.run_id, token=args.token)
    except Exception:
        return -1  # missing token

    print('Fetching targets and OGs...')
    targets = kealahou.fetch_targets()
    name_to_token = kealahou.get_target_name_to_token_map(targets)
    ogs = kealahou.fetch_ogs()
    ogs_by_target = {}
    for og in ogs:
        for comp in og['single_observing_group']['observing_block']['observing_component']:
            ogs_by_target.setdefault(comp['target_token'], []).append(og)

    missing, removed, failed = [], [], []
    for name in names:
        token = name_to_token.get(name)
        if token is None:
            missing.append(name)
            print(f'{name}: not in Kealahou program {program}')
            continue
        for og in ogs_by_target.get(token, []):
            print(f'{name}: OG {og.get("label", "?")} ({og["token"]}, priority {og["og_priority"]})')
            if args.inactive:
                if og['og_priority'] == 'INACTIVE':
                    print('   already INACTIVE')
                elif not args.dry_run:
                    og['og_priority'] = 'INACTIVE'
                    try:
                        kealahou.create_or_update_og(og)
                        print('   set INACTIVE')
                    except RuntimeError as exc:
                        print('   could not set INACTIVE:', str(exc)[:200])
                        failed.append((name, 'OG', og['token']))
                continue
            if not args.dry_run:
                j = kealahou.delete_og(og['token'])
                print('   delete OG ->', j)
                if 'error' in j:
                    failed.append((name, 'OG', og['token']))
        if args.inactive:
            if not args.dry_run and not any(f[0] == name for f in failed):
                removed.append(name)
            continue
        print(f'{name}: target token {token}')
        if not args.dry_run:
            j = kealahou.delete_target(token)
            print('   delete target ->', j)
            if 'error' in j:
                failed.append((name, 'target', token))
            else:
                removed.append(name)
    print(f'{len(removed)} {"set INACTIVE" if args.inactive else "removed"}, {len(missing)} not present, {len(failed)} failed'
          + (' (dry run: nothing deleted)' if args.dry_run else ''))
    for item in failed:
        print('  failed:', item)
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
