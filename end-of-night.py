#!/usr/bin/env python
import argparse
import requests
import json
import os
import sys
from astropy.table import Table
from datetime import UTC
from datetime import datetime
import numpy as np

from cfht_common import alias_run_ids, baseurl, datetomjd, mjdtodate

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--token', help='Kealahou API token; default is the $KEALAHOU_TOKEN environment variable')
    parser.add_argument('--run-id', help='Run ID: either NAOC (alias for 26BZ01) or LBNL (alias for 26BZ50), or an ID', default='LBNL')
    parser.add_argument('--tile-file', default='obstatus/cfht-tiles.ecsv', help='Tile filename, default %(default)s')
    parser.add_argument('--night', help='Night to look at, YYYY-MM-DD at sunset, default last night')

    args = parser.parse_args()
    token = os.environ.get('KEALAHOU_TOKEN', None)
    if args.token:
        token = args.token
    if token is None:
        print('Kealahou API token is required, either in $KEALAHOU_TOKEN environment variable or --token argument.  Create a token at https://kealahou.cfht.hawaii.edu/account/token-manager')
        return -1

    tiles = Table.read(args.tile_file)
    print('Read', len(tiles), 'tiles')

    headers = {
        'Authorization': f'Bearer {token}',
        'Content-Type': 'application/json',
    }

    run_id = args.run_id
    # lookup alias
    run_id = alias_run_ids.get(run_id, run_id)

    if args.night:
        yy,mm,dd = args.night.split('-')
        yy = int(yy)
        mm = int(mm, 10)
        dd = int(dd, 10)
        d = datetime(yy, mm, dd, tzinfo=UTC)
        mjd_min = datetomjd(d) + 1.0
    else:
        mjd_now = datetomjd(datetime.now(UTC))
        print('MJD now: %.2f' % mjd_now)
        mjd_min = np.floor(mjd_now - 0.5)
    print('MJD min:', mjd_min)
    mjd_max = mjd_min + 1.0

    r = requests.get(baseurl + 'programs/' + run_id + '/exposures', headers=headers)
    exposures = r.json()['exposure']
    print('Found', len(exposures), 'exposures total')

    keep_exposures = []
    for exp in exposures:
        st = exp['exposure_status']
        mjd = st['exp_date_mjd']
        if mjd < mjd_min or mjd > mjd_max:
            continue
        expnum = int(exp['obsid'])
        print('Exposure number', expnum, ': MJD %.2f,' % mjd, 'Hawaii time:', st['read_obs_date_hst'])
        if st['obs_comment'] == 'MEGACAM SNR SNAP':
            print('  SNR SNAP')
            continue
        exptime = st['exposure_time']['exposure_time_ms'] * 0.001
        date = mjdtodate(mjd - 1.0)
        dd = date.date()
        yymmdd = dd.strftime('%Y-%m-%d')
        target_data = exp['target_data']
        target_name = target_data['name']
        target_token = target_data['token']
        og_data = exp['observing_group_context']
        og_token = og_data['observing_group_token']
        og_label = og_data['observing_group_label']
        ot_data = exp['observing_template']
        ot_token = ot_data['token']
        print('  Target "%s"' % target_name, 'OG', og_label, '= token', og_token)
        keep_exposures.append((target_name, target_token, og_token, ot_token, expnum, yymmdd, exptime))
    print(len(keep_exposures), 'from the night')

    # Mark OGs as priority=INACTIVE
    for target_name, target_token, og_token, ot_token, expnum, yymmdd, exptime in keep_exposures:
        og_data = dict(
            token = og_token,
            og_priority = 'INACTIVE',
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
            print('OG update failed for target', target_name)
            return -1
        ent = j['entity']
        label = ent['label']
        print('OG update successful for OG labelled', label, 'target', target_name)

    # Update tile file - mark these target names as DONE=1
    obj_to_tileindex = dict([(n,i) for i,n in enumerate(tiles['OBJECT'])])
    for target_name, target_token, og_token, ot_token, expnum, yymmdd, exptime in keep_exposures:
        ind = obj_to_tileindex[target_name]
        old_done = tiles['DONE'][ind]
        tiles['DONE'][ind] = 1
        old_expnum = tiles['EXPNUM'][ind]
        tiles['EXPNUM'][ind] = expnum
        old_date = tiles['DATE'][ind]
        tiles['DATE'][ind] = yymmdd
        old_et = tiles['EXPTIME'][ind]
        exptime = np.round(exptime, decimals=3)
        tiles['EXPTIME'][ind] = exptime
        print('Marked tile', target_name, 'DONE from', old_done, 'to', 1,
              'EXPNUM from', old_expnum, 'to', expnum, 'DATE from', old_date, 'to', yymmdd,
              'EXPTIME from', old_et, 'to', exptime)
    tiles.write(args.tile_file, overwrite=True)
    print('Wrote updated tile file')

if __name__ == '__main__':
    sys.exit(main())
