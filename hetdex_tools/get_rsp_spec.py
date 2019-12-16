# -*- coding: utf-8 -*-
"""
Created on March 29 2019
@author: Erin Mentuch Cooper

Extracts fiber spectral info to run rsp3 and rsp3mc
This will be slow to do one aperture at a time. 
A bandaid for the DR1 switch to HDF5 file 
format

In the future we want to optimize to open
a shot and do a bunch of extractions at once

Either do:

python get_rsp_spec.py -ra 150.025513 -dec 2.087768 -s 20180124v010 

This assumes rad=3.0

To work directly in an rsp call 
python get_rsp_spec.py -ra $1 -dec $2 -rad $3 -s $7

Note: id isn't actually needed, everything is done in 
the current working directory
"""

from __future__ import print_function

import glob
import re

import tables as tb
import argparse as ap
import os.path as op
import numpy as np

from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.table import Table

from input_utils import setup_logging
from hetdex_api.shot import *
from hetdex_api.survey import *

def save_rsp_spectrum(self, idx, file='spec.dat'):
    spectab = Table()
    spectab['wave_rect'] = self.wave_rect
    spectab['calfib'] = self.table[idx]['calfib']
    spectab['calfibe'] = self.table[idx]['calfibe']
    spectab['Amp2Amp'] = self.table[idx]['Amp2Amp']
    spectab['Throughput'] = self.table[idx]['Throughput']
    spectab.write(file, format='ascii', overwrite=True)


def main(argv=None):
    ''' Main Function '''
    # Call initial parser from init_utils
    parser = ap.ArgumentParser(description="""Create rsp tmpXXX.datfiles.""",
                               add_help=True)

    parser.add_argument("-s", "--datevobs",
                        help='''ShotID, e.g., 20170321v009, YYYYMMDDvOBS''',
                        type=str, default=None)
    

    parser.add_argument("-ra", "--ra",
                        help='''ra, e.g., right ascension in degrees''',
                        type=float, default=None)

    parser.add_argument("-dec", "--dec",
                        help='''ra, e.g., right ascension in degrees''',
                        type=float, default=None)

    parser.add_argument("-rad", "--rad",
                        help='''radius, e.g., aperture radius in arcsec''',
                        type=float, default=3.0)

    parser.add_argument("-w", "--wave",
                        help='''wavelength in AA''',
                        type=float, default=None)

    parser.add_argument("-dw", "--dwave",
                        help='''delta wavelength in AA''',
                        type=float, default=50.0)

    
    args = parser.parse_args(argv)
    args.log = setup_logging()

    # initiate Fibers and Survey Classes
    print(args)
    fibers = Fibers(args.datevobs)
    survey = Survey('hdr1')
    
    fwhm = survey.fwhm_moffat[survey.datevobs == args.datevobs]
    structaz = survey.structaz[survey.datevobs == args.datevobs]
    ascii.write([fwhm, structaz], 'shot.info', 
                names=['fwhm_moffat', 'structaz'], overwrite=True)

    obj_coords = SkyCoord(args.ra * u.deg, args.dec * u.deg, frame='icrs')
    
    idx = fibers.query_region_idx(obj_coords, radius=(args.rad/3600.))

    output = Table()
    output['ra'] = fibers.coords.ra[idx]*u.deg
    output['dec'] = fibers.coords.dec[idx]*u.deg
    filenames = []

    fileidx = 101
    for i in idx:
        filename = 'tmp' + str(fileidx) + '.dat'
        filenames.append(filename)
        save_rsp_spectrum(fibers, i, file=filename, )
        fileidx += 1

    output['filename'] = np.array(filenames)
    ascii.write(output, 'fib_coords.dat', overwrite=True)
    
if __name__ == '__main__':
    main()
