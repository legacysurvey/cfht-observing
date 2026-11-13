#!/usr/bin/env python
"""Save every exposure record Kealahou holds for a program as JSON (read-only).

Writes queue-check/exposures-PROGRAM.json (git-ignored), for offline analysis
of exposure times, repeats and validation status.

Example:
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
    parser.add_argument('--run-id', default='NAOC', help='NAOC (26BZ01), LBNL (26BZ50), or an ID')
    parser.add_argument('--outdir', default='queue-check')
    args = parser.parse_args()
    try:
        kealahou = KealahouProgram(args.run_id, token=args.token)
    except RuntimeError:
        return -1
    exposures = kealahou.fetch_exposures()
    os.makedirs(args.outdir, exist_ok=True)
    path = os.path.join(args.outdir, 'exposures-%s.json' % args.run_id)
    with open(path, 'w') as f:
        json.dump(exposures, f, indent=1, sort_keys=True)
    print('Wrote', len(exposures), 'exposure records to', path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
