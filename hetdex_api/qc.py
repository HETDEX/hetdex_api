import json
import os

try:
    # Python 3
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError
except ImportError:
    # Python 2
    from urllib2 import urlopen, HTTPError


class AmplifierQC():

    def __init__(self, datestr):
        '''
        Initialize the AmplifierQC object for a given date

        Parameters
        ----------
        datestr : str
            The date string to check, must be of format YYYYMMDD,
            otherwise the server will return an error.
        '''

        url = 'https://luna.mpe.mpg.de/qc/' + datestr

        try:
            resp = urlopen(url)
        except HTTPError as e:
            raise Exception(' Failed to retrieve qc data, server '
                            'responded with %d %s' % (e.getcode(), e.reason))

        self._datestr = datestr
        self._qcdata = json.loads(resp.read())
        self._pristine = True

        # Attempt to read the authorization key
        try:
            self._apikey = os.environ['HETDEX_APIKEY']
        except KeyError:
            self._apikey = ''

    def get_amp(self, ifuslot, amp):
        '''
        Retrieve the quality control information for a given IFUslot and
        amplifier.

        Parameters
        ----------
        ifuslot : str
            The IFUslot to check e.g. '094' or '113'
        amp : str
            The amplifier to check one of. 'LL', 'RL', 'LU'. 'RU'

        Returns
        -------
        The stored error code for this ifuslot / amplifier combination, 0
        otherwise.

        '''

        try:
            return self._qcdata[ifuslot][amp]
        except KeyError:
            return 0

    def set_amp(self, ifuslot, amp, value):
        '''
        Update the locally stored quality control data for the given
        IFUslot and amplifier.

        Parameters
        ----------
        ifuslot : str
            The IFUslot to check e.g. '094' or '113'
        amp : str
            The amplifier to check one of. 'LL', 'RL', 'LU'. 'RU' and 'ALL'
            to set all amplifiers at once
        value : int
            The error value for this amplifier.
        '''

        if amp == 'ALL':
            for a in ['LL', 'RL', 'RU', 'LU']:
                self.set_amp(ifuslot, a, value)
        else:
            if ifuslot not in self._qcdata:
                self._qcdata[ifuslot] = {}

            self._qcdata[ifuslot][amp] = value
            self._pristine = False

    def save(self):
        '''
        Save the updated quality control data to the database.
        '''

        if self._pristine:
            return

        if self._apikey == '':
            raise Exception('Cannot save to database without api key. Please'
                            ' set \'HETDEX_APIKEY\' in the shell, or call'
                            'set_apikey.')
        url = 'https://luna.mpe.mpg.de/qc/' + self._datestr

        try:
            req = Request(url=url, data=json.dumps(self._qcdata).encode(),
                          headers={'Authorization':
                                   'FPS_API apikey='+self._apikey},
                          method='PUT')
            urlopen(req)
        except HTTPError as e:
            raise Exception(' Failed to retrieve qc data, server '
                            'responded with %d %s' % (e.getcode(), e.reason))

    def set_apikey(self, key):
        '''
        Update the apikey, if it wasn't set in the HETDEX_APIKEY shell
        environment variable.

        Parameters
        ----------
        key : str
            API key to use for saving the quality control data.

        '''

        self._apikey = key
