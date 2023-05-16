# -*- coding: utf-8 -*-                                                                                                                                                                 
""" 
@author: Erin Mentuch Cooper

This script adds FullSkyModel  group info to
an exising HDF5 file or create a new HDF5 file

python3 create_fullskymodel_hdf5.py -d 20181111 -o 15 -of 20181111v015.h5 --append

"""

import glob
import re
import os

import tables as tb
import argparse as ap
import os.path as op
import numpy as np

from astropy.table import Table

from hetdex_api.input_utils import setup_logging
from hetdex_api.config import HDRconfig


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
                        help='''Root Directory for SkySub''',
                        type=str,
                        default='/scratch/00115/gebhardt/alldet/output/')

    parser.add_argument('-of', '--outfilename', type=str,
                        help='''Relative or absolute path for output HDF5
                        file.''', default=None)

    parser.add_argument('-a', '--append',
                        help='''Appending to existing file.''',
                        action="count", default=0)

    parser.add_argument("-survey", "--survey",
                        help="""{hdr1, hdr2, hdr2.1, 'hdr3', 'hdr4'}""",
                        type=str, default="hdr4")
    
    
    args = parser.parse_args(argv)
    args.log = setup_logging()

    if op.exists(args.outfilename) and args.append:
        args.log.info('Appending fullskymodel to %s' % args.outfilename)
        fileh = tb.open_file(args.outfilename, 'a')
    
        try:
            fileh.remove_node(fileh.root.FullSkyModel, recursive=True)
        except:
            pass
    else:
        args.log.info('Creating new file for FullSkyModel %s' % args.outfilename)
        fileh = tb.open_file(args.outfilename, 'w')
        
    groupFullSkyModel = fileh.create_group(fileh.root, 'FullSkyModel', 'FullSkyModel')

    datevshot = str(args.date) + 'v' + str(args.observation).zfill(3)
    shotid = int(str(args.date) + str(args.observation).zfill(3))

    #check if shotid is in badlist

    config = HDRconfig(args.survey)
    
    badshots = np.loadtxt(config.badshot, dtype=int)

    badshotflag = False

    if shotid in badshots:
        badshotflag = True
    
    # store shuffle.cfg and DATEvOBS.log files

    for expn in ['exp01', 'exp02', 'exp03']: 

        skyfile = op.join(args.rootdir, 'd' + str(args.date) + 's' +
                          str(args.observation).zfill(3) + expn + 'sky.dat')
        try:
            sky_array = np.loadtxt(skyfile)
            if np.size(sky_array) > 0:
                fileh.create_array(groupFullSkyModel, expn, sky_array)
            else:
                if badshotflag:
                    args.log.warning('File empty %s' % skyfile)
                else:
                    args.log.error('File empty %s' % skyfile)
        except:
            if badshotflag:
                args.log.warning('Could not include %s' % skyfile)
            else:
                args.log.error('Could not include %s' % skyfile)

    fileh.close()

if __name__ == '__main__':
    main()
