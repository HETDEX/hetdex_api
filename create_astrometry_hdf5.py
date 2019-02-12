# -*- coding: utf-8 -*-                                                                                                                                                                 
""" 
@author: Erin Mentuch Cooper

This script can add Astromtry group info and data to
an exising HDF5 file or create a new HDF5 file

Example command line use:

python create_astrometry_hdf5.py -d 20181111 -o 15 -of astrometry_20181111v015.h5

or 

python create_astrometry_hdf5.py -d 20181111 -o 15 -of main_20181111v015.h5 --append

"""

import glob
import re

import tables as tb
import argparse as ap
import os.path as op
import numpy as np

from astropy.io import fits
from astropy.io import ascii
from input_utils import setup_logging

class QualityAssessment(tb.IsDescription):
    expn = tb.StringCol((5), pos=0)
    xoffset = tb.Float32Col(pos=1)
    yoffset = tb.Float32Col(pos=2)
    xrms = tb.Float32Col(pos=3)
    yrms = tb.Float32Col(pos=4)
    nstars = tb.Float32Col(pos=5)


def main(argv=None):
    ''' Main Function '''
    # Call initial parser from init_utils
    parser = ap.ArgumentParser(description="""Create HDF5 Astrometry file.""",
                               add_help=True)

    parser.add_argument("-d", "--date",
                        help='''Date, e.g., 20170321, YYYYMMDD''',
                        type=str, default=None)

    parser.add_argument("-o", "--observation",
                        help='''Observation number, "00000007" or "7"''',
                        type=str, default=None)

    parser.add_argument("-r", "--rootdir",
                        help='''Root Directory for Shifts''',
                        type=str, default='/work/00115/gebhardt/maverick/vdrp/shifts/')

    parser.add_argument('-of', '--outfilename', type=str,
                        help='''Relative or absolute path for output HDF5
                        file.''', default=None)

    parser.add_argument('-a', '--append',
                        help='''Appending to existing file.''',
                        action="count", default=0)


    args = parser.parse_args(argv)
    args.log = setup_logging()
    
    # Creates a new file if the "--append" option is not set or the file                                              
    # does not already exist.
    does_exist = False
    if op.exists(args.outfilename) and args.append:
        fileh = tb.open_file(args.outfilename, 'a')
        does_exist = True
    else:
        fileh = tb.open_file(args.outfilename, 'w')

    groupAstrometry = fileh.create_group(fileh.root, 'Astrometry', 'Astrometry Info')
    groupCoadd = fileh.create_group(groupAstrometry, 'CoaddImages', 'Coadd Images')
    groupDithall = fileh.create_group(groupAstrometry, 'Dithall', 'Fiber Astrometry Info')
    groupOffsets = fileh.create_group(groupAstrometry, 'PositionOffsets', 
                                      'Offset is star matches')

    tableQA = fileh.create_table(groupAstrometry, 'QA', QualityAssessment, 
                             'Qulity Assessment')
    filefplane = op.join(args.rootdir, str(args.date) + "v" + str(args.observation).zfill(3),
                         'fplane.txt')
    
    f = ascii.read(filefplane, names=['ifuslot', 'fpx', 'fpy', 'specid',
                                      'specslot', 'ifuid', 'ifurot', 'platesc'])
    fplanetable = fileh.create_table(groupAstrometry, 'fplane', f.as_array())

    file_stars = op.join(args.rootdir, str(args.date) + 'v' + str(args.observation).zfill(3), 
                         'shout.ifustars')
    
    f_stars = ascii.read(file_stars, names=['ignore', 'star_ID', 'ra_cat', 'dec_cat',
                                            'u', 'g', 'r', 'i', 'z'])
    starstable = fileh.create_table(groupAstrometry, 'StarCatalog',f_stars.as_array()) 

    
    for expn in ['exp01','exp02','exp03']:
        fitsfile = op.join(args.rootdir, str(args.date) + 'v' + str(args.observation).zfill(3),
                           str(args.date) + 'v' + str(args.observation).zfill(3)
                           + 'fp_' + expn + '.fits')
        
        if op.exists(fitsfile):
            
            F = fits.open(fitsfile)
            fitsim = fileh.create_array(groupCoadd, expn, F[0].data)
            fitsim.attrs['CLASS'] = 'IMAGE'
            fitsim.attrs['IMAGE_MINMAXRANGE'] = (-1.5, 100)
            fitsim.attrs['HEADER'] = F[0].header
            F.close()


            # populate offset info for catalog matches
            file_getoff = op.join(args.rootdir, str(args.date) + 'v' + str(args.observation).zfill(3),
                          'getoff_' + expn + '.out')
            f_getoff = ascii.read(file_getoff, names=['xoffset', 'yoffset', 'ra_dex', 'dec_dex',
                                                      'ra_cat','dec_cat','ifuslot'])
            getoffinfo = fileh.create_table(groupOffsets, expn, f_getoff.as_array())
            

            # populate fiber astrometry data
            file_dith = op.join(args.rootdir, str(args.date) + 'v' + str(args.observation).zfill(3),
                                'dith_' + expn + '.all')    
            f_dith = ascii.read(file_dith)
            dithinfo = fileh.create_table(groupDithall, expn, f_dith.as_array())

            #populate median and rms in offsets for quality assessment purposes
            file_getoff2 = op.join(args.rootdir, str(args.date) + 'v' + str(args.observation).zfill(3),
                                   'getoff2_' + expn + '.out')
            f_getoff2 = ascii.read(file_getoff2)
            row = tableQA.row
            row['expn'] = expn
            row['xoffset'] = f_getoff2['col1']
            row['yoffset'] = f_getoff2['col2']
            row['xrms'] = f_getoff2['col3']
            row['yrms'] = f_getoff2['col4']
            row['nstars'] = f_getoff2['col5']
            row.append()
            
            tableQA.flush()
    fileh.close()


if __name__ == '__main__':
    main()
