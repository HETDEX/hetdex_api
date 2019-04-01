# -*- coding: utf-8 -*-
"""
Created: 2019/03/28

@authors: Erin Mentuch Cooper, Dustin Davis

This file contains limited information related to the ELIXER output
for HETDEX  line detections
It is indexed to match the indexing of the input detect HDF5 file to 
allow for easy querying. Check that detectid's are consistent
when joining the tables

To match the indexing of the Detections database, it requires
a detections database H5 file.

python create_elixer_hdf5.py 

"""

import sys
import os
import os.path as op
import argparse as ap
import re
import glob
import subprocess
import numpy as np
import tables as tb

import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.io import fits
import astropy.units as u
from input_utils import setup_logging

from hetdex_api import config


class Classifications(tb.IsDescription):
    detectid = tb.Int64Col(pos=0)
    plae_poii_hetdex = tb.Float32Col(pos=1)
    plae_poii_aperture = tb.Float32Col(pos=3)
    aperture_filter = tb.StringCol((15), pos=4)
    aperture_mag = tb.Float32Col(pos=2)
    plae_poii_cat = tb.Float32Col(pos=7)
    mag_match = tb.Float32Col(pos=5)
    cat_filter = tb.StringCol((15), pos=6)
    ra = tb.Float32Col()
    dec = tb.Float32Col()
    dist_match = tb.Float32Col()
    z_prelim = tb.Float32Col()
    ra_match = tb.Float32Col()
    dec_match = tb.Float32Col()

def get_elixer_image(detectid, elix_path):
    file_jpg = op.join(elix_path, 'jpgs', str(detectid) + '.jpg')
    file_pdf = op.join(elix_path, 'pdfs', str(detectid) + '.pdf')
    
    if op.exists(file_jpg):
        elixim = plt.imread(file_jpg)
        print np.size(elixim)
    elif op.exists(file_pdf):
        file_png = str(detectid) + '.png'
        os.system('pdftoppm ' + file_pdf + ' ' + str(detectid) + ' -png -singlefile')
        elixim = plt.imread(file_png)
    return elixim

def main(argv=None):
    ''' Main Function '''
    # Call initial parser from init_utils
    parser = ap.ArgumentParser(description="""Create HDF5 file.""",
                               add_help=True)

    parser.add_argument("-df", "--detectfile",
                        help='''Provide HDF5 of detections''',
                        type=str, default=config.detecth5)
    
    parser.add_argument("-dets", "--dets",
                        help='''List of detections in form DATEvSHOT_inputID''',
                        type=str, default=None)

    parser.add_argument('-of', '--outfilename', type=str,
                        help='''Relative or absolute path for output HDF5
                        file.''', default=None)
    
    parser.add_argument('-a', '--append',
                        help='''Appending to existing detection HDF5 file.''',
                        action="count", default=0)
    
    parser.add_argument("-ep", "--elixer_path",
                        help='''Path to elixer output''',
                        type=str, default='/scratch/03261/polonius/')

    parser.add_argument("-cat", "--elixer_cat",
                        help='''Path to elixer catalog''',
                        type=str, default='/work/03261/polonius/stampede2/erin/simple_cat.txt')


    args = parser.parse_args(argv)
    args.log = setup_logging()

    # open elixer catalog file to be ingested                                                

    if op.exists(args.elixer_cat):
        colnames = ['detectid', 'ra', 'dec', 'z_prelim', 'ew_obs',
                    'ew_rest', 'plae_poii_hetdex',
                    'plae_poii_aperture', 'aperture_mag',
                    'aperture_filter', 'plae_poii_cat',
                    'cat_filter', 'dist_match', 'mag_match',
                    'ra_match', 'dec_match']
        elixer_table = ascii.read(args.elixer_cat, names=colnames, comment="#")
                                  #delimiter="\t", guess=False, format='basic')
    else:
        print('Could not open %s' % args.elixer_cat)

    filedet = tb.open_file(config.detecth5, 'r')
    detect_list = filedet.root.Detections.cols.detectid[:]
    inputid = filedet.root.Detections.cols.inputid[:]
    filedet.close()

    # open elixer HDF5 file, if append option is given this will be a group added
    # to an existing exlier HDF5 file
    if args.append:
        try:
            fileh = tb.open_file(args.outfilename, 'a')
        except:
            args.log.warning('Could not open %s to append.', args.outfilename)

    else:
        try:
            fileh = tb.open_file(args.outfilename, 'w')
            groupElix = fileh.create_group(fileh.root, 'Elixer', "ELiXer Summaries")
            fileh.create_table(fileh.root, 'Classifications', Classifications)
        except:
            args.log.warning('Could not open %s.', args.outfilename)

    tableMain = fileh.root.Classifications

    # set array of columns to store

    colkeep = ['ra', 'dec', 'z_prelim', 'plae_poii_hetdex',
                    'plae_poii_aperture', 'aperture_mag',
                    'aperture_filter', 'plae_poii_cat',
                    'cat_filter', 'dist_match', 'mag_match',
                    'ra_match', 'dec_match']

    for index, detect_i in enumerate(detect_list):
        row = tableMain.row
        row['detectid'] = detect_i

        idx = np.where(elixer_table['detectid'] == detect_i)

        if np.size(idx) > 0:
            if np.size(idx > 1):
                print(detect_i, inputid[index])

            for colname_i in colkeep:
                try:
                    if colname_i == 'aperture_filter' or colname_i == 'cat_filter':
                        row[colname_i] = str(elixer_table[colname_i][idx])[-elixer_table[colname_i][idx].size:]
                    else:
                        row[colname_i] = elixer_table[colname_i][idx]
                except Exception as e:
                    args.log.warning('Could not ingest col(%s) detectid(%d)' % (colname_i,detect_i), exc_info=True )
                    #args.log.warning('Could not ingest %s' % detect_i, exc_info=True)
        else:
            print('Could not ingest %s' % detect_i)

        row.append()

    tableMain.flush()
    fileh.close()

if __name__ == '__main__':
    main()
