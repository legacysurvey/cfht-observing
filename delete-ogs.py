#!/usr/bin/env python
import argparse
import requests
import json
import os
import sys
from astropy.table import Table

from kealahou_api import KealahouProgram

'''

Tries to delete all OGs in the given program.

Those that have exposures associated with them will fail (database integrity);
that's good.

'''

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--token', help='Kealahou API token; default is the $KEALAHOU_TOKEN environment variable')
    parser.add_argument('--run-id', help='Run ID: either NAOC (alias for 26BZ01) or LBNL (alias for 26BZ50), or an ID', default='LBNL')

    args = parser.parse_args()

    try:
        kealahou = KealahouProgram(args.run_id, token=args.token)
    except:
        # missing token
        return -1

    # Fetch all OGs for this Run ID.
    print('Fetching OGs...')
    ogs = kealahou.fetch_ogs()
    print('Found', len(ogs), 'OGs')

    # Parse & Delete OGs
    nfail = 0
    nsuccess = 0
    for og in ogs:
        og_token = og['token']
        print('Deleting OG', og_token)
        j = kealahou.delete_og(og_token)
        print('Got', j)
        if 'error' in j:
            nfail += 1
        if 'success' in j:
            nsuccess += 1

    print('%i OG deletions succeeded' % nsuccess)
    print('%i OG deletions failed' % nfail)

if __name__ == '__main__':
    sys.exit(main())
