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

class Detection:
    #represents a single detection (may or may not be an acceptable LAE
    def __init__(self):

        self.detectid = None
        self.entryid = None #unique for a run?
        self.pdfname = None
        self.detectname = None #i.e. 20180123v009_137
        self.ra = None
        self.dec = None
        self.w = None #observed wavelength

        self.z = None #assumes LyA

        self.ew_obs = None
        self.ew_rest = None #assuming LyA z
        self.line_flux = None

        self.sigma = None
        self.chi2 = None
        self.continuum = None

        self.plae_poii_hetdex = None #hetdex only data
        self.plae_poii_aperture = None
        self.aperture_mag = None
        self.aperture_filter = None

        self.neighbors = [] #empty list (of Neighbors

    @property
    def num_neighbors(self):
        try:
            return len(self.neighbors)
        except:
            return -1


class Neighbor:
    #a catalog neighbor
    def __init__(self):
        self.ra = None
        self.dec = None
        self.mag = None
        self.filter = None #filter for the magnitude
        self.dist = None #distance in arcsecs
        self.plae_poii_cat = None #P(LAE)/P(OII)


def read_elixer_catalogs(fibfn, catfn):
    """
    read the _fib.txt and _cat.txt catalogs
    create Detections for each (combining across and within catalogs as needed)

    :return: list of Detections

    Author: Dustin Davis
    """

    detections = []
    fib_version = None #can switch on version number if necessary
    cat_version = None

    #should be one line per detection in fib.txt
    #at least one line but may be more in cat.txt

    #read fib.txt first and build out list of detections
    with open(fibfn,"r") as f:
        for line in f:
            if line[0] == '#': #except get # version 1.5.0a16
                if (fib_version is None) and ('version' in line): #1st line should be the version
                    toks = line.split()
                    if (toks is not None) and (len(toks)== 3):
                        fib_version = toks[2]
                continue

            if len(line) < 100:
                continue

            toks = line.split() #white space split

            if len(toks) < 29: #must be at least 29 for _fib.txt
                continue

            d = Detection()

            d.pdfname = toks[0]
            d.detectid = np.int64(toks[1])
            d.entryid = toks[2]
            d.ra = float(toks[3])
            d.dec = float(toks[4])
            d.w = float(toks[5])
            d.sigma = float(toks[6])
            d.chi2 = float(toks[7])
            d.line_flux = float(toks[10])
            d.continuum = float(toks[12])

            d.plae_poii_hetdex = float(toks[14])

            detections.append(d)


    #then read cat.txt and match up to existing entry in detections list
    with open(catfn,"r") as f:
        for line in f:
            if line[0] == '#':
                if (cat_version is None) and ('version' in line): #1st line should be the version
                    toks = line.split()
                    if (toks is not None) and (len(toks)== 3):
                        cat_version = toks[2]
                continue

            if len(line) < 100:
                continue

            toks = line.split() #white space split

            if len(toks) < 29: #must be at least 29 for _fib.txt
                continue

            
            pdfname = toks[0]
            detectid = np.int64(toks[1])
                
            try:
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
                    plae_poii = -1 
                    
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
            except:
                print('Could not ingest cat info for %s' % detectid) 

    #Now we have consumed the catalog files
    return detections

class Classifications(tb.IsDescription):
    detectid = tb.Int64Col(pos=0)
    plae_poii_hetdex = tb.Float32Col()
    plae_poii_image = tb.Float32Col()
    image_src = tb.StringCol((20))
    plae_poii_cat = tb.Float32Col()
    cat_src = tb.StringCol((20))
    

def get_elixer_image(detectid, elix_path):
    file_jpg = op.join(elix_path, 'jpgs', str(detectid) + '.jpg')
    file_pdf = op.join(elix_path, 'pdfs', str(detectid) + '.pdf')
    if op.exists(file_jpg):
        elixim = plt.imread(file_jpg)
    elif op.exists(file_pdf):
        file_png = detectid + '.png'
        os.system('pdftoppm ' + file_pdf + ' ' + file_png + ' -png -singlefile')
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

    parser.add_argument("-name", "--elixer_name",
                        help='''Path to elixer output''',
                        type=str, default=None)


    args = parser.parse_args(argv)
    args.log = setup_logging()

    # open elixer _cat.txt file to be ingested                                                
    catfn = op.join(args.elixer_path, 'catalogs', '2017xx_cat.txt')
    fibfn = op.join(args.elixer_path, 'catalogs', '2017xx_fib.txt')

    if op.exists(catfn):
        elixer_info = read_elixer_catalogs(fibfn, catfn)
    else:
        print('Could not open %s' % catfn)

    filedet = tb.open_file(config.detecth5, 'r')
    detect_list = filedet.root.Detections.cols.detectid[:]
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

    table = fileh.root.Classifications
  
    for detect_i in detect_list[0:50]:
        
        args.log.info('Ingesting detectid = %s' % detect_i)
        try:
            elixim = get_elixer_image(detect_i, args.elixer_path)
            elixarray = fileh.create_array(groupElix, 'DEX' + str(detect_i), elixim)
            elixarray.attrs['CLASS'] = 'IMAGE'
        except:
            args.log.warning('Could not open elixer image for %s' % detect_i)
            
    table.flush()
    fileh.close()

if __name__ == '__main__':
    main()
