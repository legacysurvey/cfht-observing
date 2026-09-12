#!/usr/bin/env python
import argparse
import requests
import json
import os
import sys
from astropy.table import Table

from cfht_common import alias_run_ids, baseurl

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--token', help='Kealahou API token; default is the $KEALAHOU_TOKEN environment variable')
    parser.add_argument('--run-id', help='Run ID: either NAOC (alias for 26BZ01) or LBNL (alias for 26BZ50), or an ID', default='LBNL')
    parser.add_argument('--ot-name', default='M4376_0.6', help='Observing Template (OT) name')
    parser.add_argument('planfile', help='Plan filename (FITS)')

    args = parser.parse_args()
    token = os.environ.get('KEALAHOU_TOKEN', None)
    if args.token:
        token = args.token
    if token is None:
        print('Kealahou API token is required, either in $KEALAHOU_TOKEN environment variable or --token argument.  Create a token at https://kealahou.cfht.hawaii.edu/account/token-manager')
        return -1

    plan = Table.read(args.planfile)
    print('Read', len(plan), 'targets from plan file')

    headers = {
        'Authorization': f'Bearer {token}',
        'Content-Type': 'application/json',
    }

    run_id = args.run_id
    # lookup alias
    run_id = alias_run_ids.get(run_id, run_id)

    # Fetch all targets for this Run ID.
    print('Fetching targets...')
    r = requests.get(baseurl + 'programs/' + run_id + '/targets', headers=headers)

    target_name_to_token = {}
    target_token_to_name = {}
    for e in r.json()['entity']:
        k = e['name']
        v = e['token']
        target_name_to_token[k] = v
        target_token_to_name[v] = k
    print(len(target_name_to_token), 'unique target names and',
          len(target_token_to_name), 'unique tokens')

    # Add any targets that don't exist yet
    for target_name,ra,dec,target_mag in zip(plan['NAME'], plan['RA'], plan['DEC'],
                                             plan['MAG_AB']):
        target_name = target_name.strip()
        token = run_id + '-FT-' + target_name
        entity = dict(
            token=token,
            name=target_name,
            fixed_target = dict(
                coordinate=dict(ra=float(ra), dec=float(dec)), proper_motion=dict()),
            magnitude = dict(AB=dict(value=float(target_mag))),
            # HACK -- ignoring the pointing_offset in the target file!
            # HACK - hard-coded pointing offset token.
            pointing_offset=dict(token='00AZ00-PO+MEGACAM+1', name='1', offset={}))
        req_data = dict(entity=entity)
        url = baseurl + 'programs/' + run_id + '/targets/' + token
        r = requests.put(url, headers=headers, data=json.dumps(req_data).encode('utf-8'))
        j = r.json()
        print('got', j)
        if not j['success']:
            print('Target creation request failed for target', target_name)
            return -1
        ent = j['entity']
        label = ent['label']
        print('Added target', target_name)

    # Fetch all OTs for this Run ID
    print('Fetching OTs...')
    r = requests.get(baseurl + 'programs/' + run_id + '/observing-templates', headers=headers)
    ot_map = dict()
    for d in r.json()['entity']:
        v = d['token']
        k = d['name']
        ot_map[k] = v
        print('  OT:', k)

    ot_name = args.ot_name
    if not ot_name in ot_map:
        print('Error: unknown OT name "%s"' % ot_name)
        return -1
    ot_token = ot_map[ot_name]

    # Fetch all OGs for this Run ID.
    print('Fetching OGs...')
    r = requests.get(baseurl + 'programs/' + run_id + '/observing-groups', headers=headers)
    ogs = r.json()['entity']
    print('Found', len(ogs), 'OGs')

    # Parse OGs
    og_tokens = []
    og_priority = {}
    for og in ogs:
        print('og', og)
        og_token = og['token']
        og_priority[og_token] = og['og_priority']
        components = og['single_observing_group']['observing_block']['observing_component']
        assert(len(components) == 1)
        comp = components[0]
        ot_token = comp['observing_template_token']
        target_token = comp['target_token']
        og_tokens.append((og_token, ot_token, target_token))
    #print(len(og_tokens), 'OGs parsed')
    targetname_to_og = {}
    for og_token, ot_token, target_token in og_tokens:
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
            # OG token names seem to be arbitrary in Kealahou,
            # but we may as well make them structured.
            og_token = run_id + '-' + ot_name + '-' + target_name
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
        req_data = dict(entity=og_data)
        url = baseurl + 'programs/' + run_id + '/observing-groups/' + og_token
        r = requests.put(url, headers=headers, data=json.dumps(req_data).encode('utf-8'))
        j = r.json()
        if not j['success']:
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
