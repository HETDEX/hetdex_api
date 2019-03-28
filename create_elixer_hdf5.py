# -*- coding: utf-8 -*-
"""
Created: 2019/03/28

@author: Erin Mentuch Cooper

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

from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.io import fits
import astropy.units as u
from input_utils import setup_logging

def read_elixer_cat(catfn):

    """
    read  _cat.txt catalogs
    create Detections for each (combining across and within catalogs as needed)

    :return: list of Detections
    """
    
    detections = []

    with open(catfn,"r") as f:
        for line in f:
            if line[0] == '#':
                continue
            if len(line) < 100:
                continue

            toks = line.split() #white space split

            if len(toks) < 29: #must be at least 29 for _fib.txt
                continue

            pdfname = toks[0]
            detectid = np.int64(toks[1])
            entryid = toks[2]

            ra = float(toks[3])
            dec = float(toks[4])

            num_cat_matches = int(toks[5])

            w = float(toks[6])

            match_ra = float(toks[24])
            match_dec = float(toks[25])
            dist = float(toks[26])
            mag = toks[28]
            filter = toks[29]
            if (toks[32] is not None) and (toks[32].lower() != 'none'): #can be None if could not get an aperture
                plae_poii = float(toks[32])
            else:
                plae_poii = -1 #this does me no good here, so skip this line
                print("Bad (None) PLAE/POIU. Skipping %s %s %s" %(pdfname,str(detectid),entryid))
                continue

            #find your mates:
            #there may be several ... the 1st is the aperture
            #each additional one is a catalog match
            for m in detections:
                if (m.entryid == entryid) and (m.detectid == detectid) and (m.ra == ra) and (m.dec == dec) and (m.w == w):
                    #match
                    if (match_ra == 666) and (match_dec == 666) and (dist == 0.0): #this is the aperture entry
                        m.plae_poii_aperture = plae_poii
                        m.aperture_filter = filter
                        m.aperture_mag = mag
                    else: #this is a matched catalog object
                        n = Neighbor()
                        n.ra = match_ra
                        n.dec = match_dec
                        n.mag = mag
                        n.filter = filter  # filter for the magnitude
                        n.dist = dist  # distance in arcsecs
                        n.plae_poii = plae_poii  # P(LAE)/P(OII)

                        m.num_neighbors += 1
                        m.neighbors.append(n)

    #Now we have consumed the catalog files
    return detections

class Classifications(tb.IsDescription):
    detectid = tb.Int64Col(pos=0)
    plae_poii_hetdex = tb.Float32Col()
    plae_poii_image = tb.Float32Col()
    image_src = tb.StringCol((20))
    plae_poii_cat = tb.Float32Col()
    cat_src = tb.StringCol((20))
    
def main(argv=None):
    ''' Main Function '''
    # Call initial parser from init_utils
    parser = ap.ArgumentParser(description="""Create HDF5 file.""",
                               add_help=True)

    parser.add_argument("-df", "--detectfile",
                        help='''Provide HDF5 of detections''',
                        type=str, default='/work/05350/ecooper/hdr1/detect/detect_hdr1.h5')
    
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
                        type=str, default='/work/05350/ecooper/stampede2/elixer/')

    parser.add_argument("-name", "--elixer_name",
                        help='''Path to elixer output''',
                        type=str, default=None)


    args = parser.parse_args(argv)
    args.log = setup_logging()

    # open elixer _cat.txt file to be ingested                                                
    catfn = op.join(args.elixer_path, str(args.elixer_name) + '_cat.txt')

    if op.exists(catfn):
        elixer_info = read_elixer_cat(catfn)
    else:
        print('Could not open %s' % catfn)

    # open elixer HDF5 file, if append option is given this will be a group added
    # to an existing exlier HDF5 file
    try:
        filedet = tb.open_file(args.detectfile, 'r')
        detect_list = filedet.root.Detections.cols.detectid[:]
        filedet.close()
    except:
        args.log.warning('Could not open %s' % args.detectfile)

    if args.append:
        try:
            fileh = tb.open_file(args.outfilename, 'a')
        except:
            args.log.warning('Could not open %s to append.', args.outfilename)

    else:
        try:
            fileh = tb.open_file(args.outfilename, 'w')
            fileh.create_group(fileh.root, 'Elixer', "Line Emission Classifications")
            fileh.create_table(fileh.root, 'Classifications', Classifications)
        except:
            args.log.warning('Could not open %s.', args.detectfile)

    table = fileh.root.Classifications
  
    for detect in detect_list:
        print detect

    table.flush()
    fileh.close()

if __name__ == '__main__':
    main()
