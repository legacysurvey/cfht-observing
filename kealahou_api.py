import requests
import json
import os

from cfht_common import alias_run_ids, baseurl

class KealahouProgram(object):
    def __init__(self, run_id, token=None):
        self.baseurl = baseurl
        self.run_id = alias_run_ids.get(run_id, run_id)

        mytoken = os.environ.get('KEALAHOU_TOKEN', None)
        if token:
            mytoken = token
        if mytoken is None:
            print('Kealahou API token is required, either in $KEALAHOU_TOKEN environment variable or --token argument.  Create a token at https://kealahou.cfht.hawaii.edu/account/token-manager')
            raise RuntimeError('Kealahou token required')

        self.headers = {
            'Authorization': 'Bearer %s' % mytoken,
            'Content-Type': 'application/json',
        }

    def fetch_exposures(self):
        r = requests.get(self.baseurl + 'programs/' + self.run_id + '/exposures',
                         headers=self.headers)
        exposures = r.json()['exposure']
        return exposures

    def create_or_update_og(self, og_data):
        req_data = dict(entity=og_data)
        og_token = og_data['token']
        url = self.baseurl + 'programs/' + self.run_id + '/observing-groups/' + og_token
        r = requests.put(url, headers=self.headers,
                         data=json.dumps(req_data).encode('utf-8'))
        j = r.json()
        if not j['success']:
            raise RuntimeError('Kealahou: failed to update/create OG:', og_data)
        return j
