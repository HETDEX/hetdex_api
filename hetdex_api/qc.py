import json
import os
import copy
import six

try:
    # Python 3
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError
except ImportError:
    # Python 2
    from urllib2 import urlopen, HTTPError
    import urllib2

    class Request(urllib2.Request):
        def __init__(self, *args, **kwargs):
            if 'method' in kwargs:
                self._method = kwargs['method']
                del kwargs['method']
            else:
                self._method = None
            return urllib2.Request.__init__(self, *args, **kwargs)

        def get_method(self, *args, **kwargs):
            if self._method is not None:
                return self._method
            return urllib2.Request.get_method(self, *args, **kwargs)

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
        self._original = copy.deepcopy(self._qcdata)

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

    def save(self, force=False):
        '''
        Save the updated quality control data to the database.
        '''

        if self._original == self._qcdata:
            print('No changes found, not saving to database.')
            return

        if self._apikey == '':
            raise Exception('Cannot save to database without api key. Please'
                            ' set \'HETDEX_APIKEY\' in the shell, or call'
                            'set_apikey.')
        url = 'https://luna.mpe.mpg.de/qc/' + self._datestr

        try:
            datadict = {'new': self._qcdata,
                        'old': self._original}
            if force:
                datadict['forceupdate'] = True
            req = Request(url=url,
                          data=json.dumps(datadict).encode(),
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


class AllAmplifierQC():

    def __init__(self, filename=None):
        '''
        Initialize the AmplifierQC object for a given date

        Parameters
        ----------
        filename : str, optional
            An optional JSON file to read the data from. Using a file
            disables saving of updates to the database.
        '''

        if filename != None:
            print('Attempting read from file')
            with open(filename, 'r') as f:
                self._qcdata = json.load(f)
            self._original = None
        else:
            url = 'https://luna.mpe.mpg.de/qc/all'

            try:
                resp = urlopen(url)
            except HTTPError as e:
                raise Exception(' Failed to retrieve qc data, server '
                                'responded with %d %s' % (e.getcode(),
                                                          e.reason))

            qcdata = json.loads(resp.read())
            self._qcdata = {}
            for k in qcdata:
                self._qcdata[int(k)] = qcdata[k]
            self._original = copy.deepcopy(self._qcdata)

        # Attempt to read the authorization key
        try:
            self._apikey = os.environ['HETDEX_APIKEY']
        except KeyError:
            self._apikey = ''

    def get_amp(self, date, ifuslot, amp):
        '''
        Retrieve the quality control information for a given IFUslot and
        amplifier.

        Parameters
        ----------
        date : int or str
            Date for the qc data
        ifuslot : str
            The IFUslot to check e.g. '094' or '113'
        amp : str
            The amplifier to check one of. 'LL', 'RL', 'LU'. 'RU'

        Returns
        -------
        The stored error code for this ifuslot / amplifier combination, 0
        otherwise.

        '''
        if isinstance(date, six.string_types):
            date = int(date)

        if date not in self._qcdata:
            print('WARNING! No data for date %d' % date)
            return 0
        try:
            return self._qcdata[date][ifuslot][amp]
        except KeyError:
            return 0

    def set_amp(self, date, ifuslot, amp, value):
        '''
        Update the locally stored quality control data for the given
        IFUslot and amplifier.

        Parameters
        ----------
        date : int or str
            Date for the qc data
        ifuslot : str
            The IFUslot to check e.g. '094' or '113'
        amp : str
            The amplifier to check one of. 'LL', 'RL', 'LU'. 'RU' and 'ALL'
            to set all amplifiers at once
        value : int
            The error value for this amplifier.
        '''
        if isinstance(date, six.string_types):
            date = int(date)

        if date not in self._qcdata:
            self._qcdata[date] = {}
        if amp == 'ALL':
            for a in ['LL', 'RL', 'RU', 'LU']:
                self.set_amp(date, ifuslot, a, value)
        else:
            if ifuslot not in self._qcdata[date]:
                self._qcdata[date][ifuslot] = {}

            self._qcdata[date][ifuslot][amp] = value

    def save(self, force=False):
        '''
        Save the updated quality control data to the database.
        '''

        if self._original is None:
            raise Exception('Cannot save, data was loaded from json file!')
        if self._original == self._qcdata:
            print('No changes found, not saving to database.')
            return

        if self._apikey == '':
            raise Exception('Cannot save to database without api key. Please'
                            ' set \'HETDEX_APIKEY\' in the shell, or call'
                            'set_apikey.')
        url = 'https://luna.mpe.mpg.de/qc/'

        try:
            for d in self._qcdata.keys():
                if d not in self._original:
                    self._original[d] = {}

                if self._qcdata[d] == self._original[d]:
                    continue

                datadict = {'new': self._qcdata[d],
                            'old': self._original[d]}
                if force:
                    datadict['forceupdate'] = True
                req = Request(url=url + str(d),
                              data=json.dumps(datadict).encode(),
                              headers={'Authorization':
                                       'FPS_API apikey='+self._apikey},
                              method='PUT')
                urlopen(req)
        except HTTPError as e:
            raise Exception(' Failed to save qc data, server '
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
