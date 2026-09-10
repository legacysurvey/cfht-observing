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

    '''
    For information about the API, see:

    https://github.com/CFHT/kealahou-workshop
    including examples like
    https://github.com/CFHT/kealahou-workshop/blob/main/example.py

    https://swagger.cfht.hawaii.edu/


    Overall comments:

    * Fixed targets have RA,Dec,mag -- we can bake in the SFD dust
    correction to the mags.  No need for API control of these - can just use
    the bulk upload functionality.

    * OTs have filter, SNR goal, exposure time, max IQ (in r), max airmass.  We
      expect to use a small number of these, so we will create them by hand.

    * OGs have priority (eg, INACTIVE), target and OT.  This is the main object
      that we need API control over.
    '''

    plan = Table.read(args.planfile)
    print('Read', len(plan), 'targets from plan file')

    headers = {
        'Authorization': f'Bearer {token}',
        'Content-Type': 'application/json',
    }

    # List programs for this user
    # r = requests.get(baseurl + 'programs', headers=headers)

    run_id = args.run_id
    # lookup alias
    run_id = alias_run_ids.get(run_id, run_id)

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

    # Fetch all OGs for this Run ID.
    print('Fetching OGs...')
    r = requests.get(baseurl + 'programs/' + run_id + '/observing-groups', headers=headers)
    ogs = r.json()['entity']
    print('Found', len(ogs), 'OGs')

    # Parse OGs
    og_tokens = []
    for og in ogs:
        og_token = og['token']
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
    for target_name in plan['name']:
        target_name = target_name.strip()
        if not target_name in targetname_to_og:
            new_targets.append(target_name)
    print('Need to create', len(new_targets), 'new OGs')

    # Create an OG for each target in the file that does not already have an OG.
    for i,target_name in enumerate(new_targets):
        target_token = target_name_to_token[target_name]
        # OG token names seem to be arbitrary, but we may as well make them structured.
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
            print('OG creation request failed for target', target_name)
            return -1
        ent = j['entity']
        label = ent['label']
        print('Added OG', (i+1), 'of', len(new_targets),
              'labelled', label, 'in Kealahou, for target', target_name)
    print('Created all new OGs')

if __name__ == '__main__':
    sys.exit(main())
