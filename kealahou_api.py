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

    def fetch_targets(self):
        r = requests.get(self.baseurl + 'programs/' + self.run_id + '/targets',
                         headers=self.headers)
        return r.json()['entity']

    def fetch_ots(self):
        r = requests.get(self.baseurl + 'programs/' + self.run_id + '/observing-templates',
                         headers=self.headers)
        return r.json()['entity']

    def fetch_ot_map(self):
        ot_map = dict()
        for d in self.fetch_ots():
            v = d['token']
            k = d['name']
            ot_map[k] = v
        return ot_map

    def fetch_ogs(self):
        r = requests.get(self.baseurl + 'programs/' + self.run_id + '/observing-groups',
                         headers=self.headers)
        return r.json()['entity']

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

    def create_or_update_target(self, target_data):
        req_data = dict(entity=target_data)
        token = target_data['token']
        url = self.baseurl + 'programs/' + self.run_id + '/targets/' + token
        r = requests.put(url, headers=self.headers,
                         data=json.dumps(req_data).encode('utf-8'))
        j = r.json()
        if not j['success']:
            raise RuntimeError('Kealahou: failed to update/create Target:', target_data)
        return j

    def delete_target(self, target_token):
        url = self.baseurl + 'programs/' + self.run_id + '/targets/' + target_token
        r = requests.delete(url, headers=self.headers)
        j = r.json()
        return j

    def delete_og(self, og_token):
        url = self.baseurl + 'programs/' + self.run_id + '/observing-groups/' + og_token
        r = requests.delete(url, headers=self.headers)
        j = r.json()
        return j

    def token_for_target(self, target_name):
        # Return the token name to be used for the given target_name.
        token = self.run_id + '-FT-' + target_name
        return token

    def token_for_og(self, ot_name, target_name):
        # OG token names seem to be arbitrary in Kealahou,
        # but we may as well make them structured.
        token = self.run_id + '-' + ot_name + '-' + target_name
        return token

    def get_target_name_to_token_map(self, targets):
        target_name_to_token = {}
        for target in targets:
            k = target['name']
            v = target['token']
            target_name_to_token[k] = v
        return target_name_to_token

    def get_target_token_to_name_map(self, targets):
        target_token_to_name = {}
        for target in targets:
            k = target['name']
            v = target['token']
            target_token_to_name[v] = k
        return target_token_to_name
