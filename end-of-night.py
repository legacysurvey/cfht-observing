#!/usr/bin/env python
import argparse
import requests
import json
import os
import sys
from astropy.table import Table
from datetime import datetime, timezone
import numpy as np

from cfht_common import SNR_GOAL, alias_run_ids, baseurl, datetomjd, exposure_quality, mjdtodate
from kealahou_api import KealahouProgram

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--token', help='Kealahou API token; default is the $KEALAHOU_TOKEN environment variable')
    parser.add_argument('--run-id', help='Run ID: either NAOC (alias for 26BZ01) or LBNL (alias for 26BZ50), or an ID; default is LBNL', default='LBNL')
    parser.add_argument('--tile-file', default='obstatus/cfht-tiles.ecsv', help='Tile filename, default %(default)s')
    parser.add_argument('--snr-goal', type=float, default=SNR_GOAL,
                        help='SNR goal; an exposure is good when its SNR exceeds this minus 0.4 (default %(default)s)')
    parser.add_argument('--night', help='Night to look at, YYYY-MM-DD at sunset (Hawaii time), default last night')

    args = parser.parse_args()

    try:
        kealahou = KealahouProgram(args.run_id, token=args.token)
    except:
        # missing token
        return -1

    tiles = Table.read(args.tile_file)
    print('Read', len(tiles), 'tiles from', args.tile_file)
    if 'SNR' not in tiles.colnames:
        tiles['SNR'] = np.zeros(len(tiles), dtype=np.float32)

    if args.night:
        yy,mm,dd = args.night.split('-')
        yy = int(yy)
        mm = int(mm, 10)
        dd = int(dd, 10)
        d = datetime(yy, mm, dd, tzinfo=timezone.utc)
        mjd_min = datetomjd(d) + 1.0
    else:
        mjd_now = datetomjd(datetime.now(timezone.utc))
        print('MJD now: %.2f' % mjd_now)
        mjd_min = np.floor(mjd_now - 0.5)
    print('MJD min:', mjd_min)
    mjd_max = mjd_min + 1.0

    print('Fetching exposures from Kealahou...')
    exposures = kealahou.fetch_exposures()
    print('Found', len(exposures), 'exposures total')

    # This magic 1.0 is because the MJD is UTC, and we want the sunset date in Hawaii time.
    MJD_DATE_OFFSET = 1.0

    obj_to_tileindex = dict([(str(n).strip(),i) for i,n in enumerate(tiles['OBJECT'])])
    keep_exposures = []
    not_validated = []
    n_catchup = 0
    exposures_mjd_max = 0.
    for exp in exposures:
        st = exp['exposure_status']
        mjd = st['exp_date_mjd']
        exposures_mjd_max = max(exposures_mjd_max, mjd)
        expnum = int(exp['obsid'])
        if mjd < mjd_min or mjd > mjd_max:
            # Catch-up: a validated exposure from another night that the tile file does
            # not have yet (validated late, or taken after that night's update was run).
            name = (exp.get('target_data') or {}).get('name')
            quality, snr = exposure_quality(st, args.snr_goal)
            if quality is None or name not in obj_to_tileindex:
                continue
            ind = obj_to_tileindex[name]
            if tiles['DONE'][ind] != 0 and (tiles['EXPNUM'][ind] == expnum or (snr or 0) <= tiles['SNR'][ind] + 1e-6):
                continue
            n_catchup += 1
            print('Catch-up: exposure', expnum, 'of', name, 'from', st['read_obs_date_hst'],
                  'counts (%s, SNR %s) but is not in the tile file' % (quality, snr))
        else:
            print('Exposure number', expnum, ': MJD %.2f,' % mjd, 'Hawaii time:', st['read_obs_date_hst'])
        if st['obs_comment'] == 'MEGACAM SNR SNAP':
            print('  SNR SNAP')
            continue
        # An exposure counts when it is good (graded, not marked bad, SNR > goal - 0.4)
        # or when CFHT validated it at lower SNR (counted for now, to be topped up
        # later). Anything else leaves the tile unobserved and its OG active, so that
        # the queue repeats it.
        quality, snr = exposure_quality(st, args.snr_goal)
        if quality is None:
            print('  Not counted: status %s, grade %s, SNR %s' % (st.get('exp_status'), st.get('grade'), snr))
            not_validated.append((expnum, st.get('exp_status'), st.get('grade')))
            continue
        if quality == 'low-snr':
            print('  Validated by CFHT with low SNR %s: counted for now' % snr)
        exptime = st['exposure_time']['exposure_time_ms'] * 0.001
        date = mjdtodate(mjd - MJD_DATE_OFFSET)
        dd = date.date()
        yymmdd = dd.strftime('%Y-%m-%d')
        target_data = exp.get('target_data')
        if not target_data:
            print('  No target attached; skipping')
            continue
        target_name = target_data['name']
        target_token = target_data['token']
        og_data = exp['observing_group_context']
        og_token = og_data['observing_group_token']
        og_label = og_data['observing_group_label']
        og_active = og_data.get('observing_group_priority') != 'INACTIVE'
        ot_data = exp['observing_template']
        ot_token = ot_data['token']
        print('  Target "%s"' % target_name, 'OG', og_label, '= token', og_token)
        keep_exposures.append((target_name, target_token, og_token, ot_token, expnum, yymmdd, exptime, og_active, snr or 0.))
    print('Keeping', len(keep_exposures) - n_catchup, 'exposures from the night and',
          n_catchup, 'caught up from other nights;', len(not_validated), 'from the night not counted')
    if any(grade is None for _, _, grade in not_validated):
        print('Some exposures are not graded yet; a later run (for any night) will pick them up once graded.')
    if len(keep_exposures) == 0:
        print('Latest exposure seen: MJD %.2f, date %s' % (exposures_mjd_max, mjdtodate(exposures_mjd_max - MJD_DATE_OFFSET).strftime('%Y-%m-%d')))
        return -1

    # Mark OGs as priority=INACTIVE
    for target_name, target_token, og_token, ot_token, expnum, yymmdd, exptime, og_active, snr in keep_exposures:
        if not og_active:
            continue
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
        j = kealahou.create_or_update_og(og_data)
        ent = j['entity']
        label = ent['label']
        print('OG update successful for OG labelled', label, 'target', target_name)

    # Update tile file - mark these target names as DONE=1
    for target_name, target_token, og_token, ot_token, expnum, yymmdd, exptime, og_active, snr in keep_exposures:
        if target_name not in obj_to_tileindex:
            # Not a survey tile (e.g. the COSMOS dither targets of make-dither-targets.py).
            record_non_survey_exposure(target_name, expnum, yymmdd, exptime)
            continue
        ind = obj_to_tileindex[target_name]
        old_done = tiles['DONE'][ind]
        tiles['DONE'][ind] = 1
        old_expnum = tiles['EXPNUM'][ind]
        old_date = tiles['DATE'][ind]
        old_et = tiles['EXPTIME'][ind]
        # The tile file records the highest-SNR exposure that counts.
        if old_done == 0 or (expnum != old_expnum and snr > tiles['SNR'][ind] + 1e-6):
            # update
            tiles['EXPNUM'][ind] = expnum
            tiles['DATE'][ind] = yymmdd
            exptime = np.round(exptime, decimals=3)
            tiles['EXPTIME'][ind] = exptime
            tiles['SNR'][ind] = snr
            print('Marked tile', target_name, 'DONE from', old_done, 'to', 1,
                  'EXPNUM from', old_expnum, 'to', expnum, 'DATE from', old_date, 'to', yymmdd,
                  'EXPTIME from', old_et, 'to', exptime)
        else:
            print('Skipped tile', target_name, 'that was already populated:'
                  'EXPNUM', old_expnum, 'vs current', expnum, 'DATE', old_date, 'vs current', yymmdd)
    tiles.write(args.tile_file, overwrite=True)
    print('Wrote updated tile file')

def record_non_survey_exposure(target_name, expnum, yymmdd, exptime,
                               filename='obstatus/non-survey-exposures.txt'):
    '''Log an exposure of a target that is not in the tile file.'''
    header = '# TARGET EXPNUM DATE EXPTIME\n'
    lines = open(filename).readlines() if os.path.exists(filename) else [header]
    if any(line.split()[1:2] == [str(expnum)] for line in lines if not line.startswith('#')):
        print('Non-survey target', target_name, 'exposure', expnum, 'already recorded')
        return
    lines.append('%s %i %s %.3f\n' % (target_name, expnum, yymmdd, exptime))
    with open(filename, 'w') as f:
        f.writelines(lines)
    print('Recorded non-survey target', target_name, 'exposure', expnum, 'in', filename)


if __name__ == '__main__':
    sys.exit(main())
