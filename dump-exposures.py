#!/usr/bin/env python
"""Save everything Kealahou holds for a program as JSON (read-only).

Writes queue-check/exposures-PROGRAM.json, targets-PROGRAM.json and
ogs-PROGRAM.json (git-ignored), for offline analysis and for
`check-queues.py --from-dump queue-check`. By default both programs are saved.

Example:
    python dump-exposures.py
    python dump-exposures.py --run-id NAOC
"""
import argparse
import json
import os
import sys

from kealahou_api import KealahouProgram


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--token', help='Kealahou API token; default is $KEALAHOU_TOKEN')
    parser.add_argument('--run-id', nargs='+', default=['NAOC', 'LBNL'], help='NAOC (26BZ01), LBNL (26BZ50), or IDs; default both')
    parser.add_argument('--outdir', default='queue-check')
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    for run_id in args.run_id:
        try:
            kealahou = KealahouProgram(run_id, token=args.token)
        except RuntimeError:
            return -1
        for kind, fetch in (('exposures', kealahou.fetch_exposures), ('targets', kealahou.fetch_targets),
                            ('ogs', kealahou.fetch_ogs)):
            records = fetch()
            path = os.path.join(args.outdir, '%s-%s.json' % (kind, run_id))
            with open(path, 'w') as f:
                json.dump(records, f, indent=1, sort_keys=True)
            print('Wrote', len(records), kind, 'records to', path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
