#!/usr/bin/env python
import argparse
import requests
import json
import os
import sys
from astropy.table import Table

from kealahou_api import KealahouProgram

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--token', help='Kealahou API token; default is the $KEALAHOU_TOKEN environment variable')
    parser.add_argument('--run-id', help='Run ID: either NAOC (alias for 26BZ01) or LBNL (alias for 26BZ50), or an ID', default='LBNL')
    parser.add_argument('--ot-name', default='M4376_0.6', help='Observing Template (OT) name')
    parser.add_argument('planfile', help='Plan filename (FITS)')

    args = parser.parse_args()

    try:
        kealahou = KealahouProgram(args.run_id, token=args.token)
    except:
        # missing token
        return -1

    plan = Table.read(args.planfile)
    print('Read', len(plan), 'targets from plan file')

    # Fetch all targets for this Run ID.
    print('Fetching targets...')
    targets = kealahou.fetch_targets()
    target_name_to_token = kealahou.get_target_name_to_token_map(targets)
    target_token_to_name = kealahou.get_target_token_to_name_map(targets)
    print(len(target_name_to_token), 'unique target names and',
          len(target_token_to_name), 'unique target tokens')
    if len(target_name_to_token) != len(target_token_to_name):
        print('Warning: target names are not unique?')

    # Add any targets that don't exist yet
    for target_name,ra,dec,target_mag in zip(plan['NAME'], plan['RA'], plan['DEC'],
                                             plan['MAG_AB']):
        target_name = target_name.strip()
        # FIXME?  We could just unconditionally create/update...
        if target_name in target_name_to_token:
            print('Target already exists:', target_name)
            continue
        token = kealahou.token_for_target(target_name)
        target_data = dict(
            token=token,
            name=target_name,
            fixed_target = dict(
                coordinate=dict(ra=float(ra), dec=float(dec)), proper_motion=dict()),
            magnitude = dict(AB=dict(value=float(target_mag))),
            # HACK -- ignoring the pointing_offset in the target file!
            # HACK - hard-coded pointing offset token.
            pointing_offset=dict(token='00AZ00-PO+MEGACAM+1', name='1', offset={}))
        j = kealahou.create_or_update_target(target_data)
        ent = j['entity']
        label = ent['label']
        print('Added target', target_name)
        target_name_to_token[target_name] = token
        target_token_to_name[token] = target_name

    # Fetch all OTs for this Run ID
    print('Fetching OTs...')
    ot_map = kealahou.fetch_ot_map()
    ot_name = args.ot_name
    if not ot_name in ot_map:
        print('Error: unknown OT name "%s"' % ot_name)
        return -1
    new_ot_token = ot_map[ot_name]

    # Fetch all OGs for this Run ID.
    print('Fetching OGs...')
    ogs = kealahou.fetch_ogs()
    print('Found', len(ogs), 'OGs')

    # Parse OGs
    targetname_to_og = {}
    og_priority = {}
    og_ot = {}
    for og in ogs:
        og_token = og['token']
        og_priority[og_token] = og['og_priority']
        components = og['single_observing_group']['observing_block']['observing_component']
        assert(len(components) == 1)
        comp = components[0]
        og_ot[og_token] = comp['observing_template_token']
        target_token = comp['target_token']
        target_name = target_token_to_name[target_token]
        targetname_to_og[target_name] = og_token

    # Find targets without OGs.
    new_targets = []
    inactive_ogs = []
    for target_name in plan['NAME']:
        target_name = target_name.strip()
        if not target_name in targetname_to_og:
            new_targets.append((target_name, None))
        else:
            og_token = targetname_to_og[target_name]
            og_prio = og_priority[og_token]
            if og_prio == 'INACTIVE':
                inactive_ogs.append((target_name, og_token))
    print('Need to create', len(new_targets), 'new OGs')
    print('Need to set', len(inactive_ogs), 'from INACTIVE back to normal')

    # Create an OG for each target in the file that does not already have an OG.
    # It turns out that the same code can be used to set the priority back to MEDIUM.
    for i,(target_name,og_token) in enumerate(new_targets + inactive_ogs):
        target_token = target_name_to_token[target_name]
        if og_token is None:
            og_token = kealahou.token_for_og(ot_name, target_name)
            # We're making a new OG
            ot_token = new_ot_token
        else:
            ot_token = og_ot[og_token]
        og_data = dict(
            token = og_token,
            og_priority = 'MEDIUM',
            target_type = 'OBJECT',
            single_observing_group = dict(
                observing_block = dict(
                    observing_component = [dict(
                        target_token = target_token,
                        observing_template_token = ot_token,
                    )],),),)
        try:
            j = kealahou.create_or_update_og(og_data)
        except:
            import traceback
            traceback.print_exc()
            if i < len(new_targets):
                print('OG creation request failed for target', target_name)
            else:
                print('OG priority update failed for target', target_name)
            return -1
        ent = j['entity']
        label = ent['label']
        if i < len(new_targets):
            print('Added OG', (i+1), 'of', len(new_targets), 'for target', target_name)
            #'labelled', label, 'in Kealahou,
        else:
            print('Updated OG', (i+1-len(new_targets)), 'of', len(inactive_ogs))
    print('Created/updated all OGs')

if __name__ == '__main__':
    sys.exit(main())
