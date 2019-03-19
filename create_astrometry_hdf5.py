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

import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from input_utils import setup_logging

class QualityAssessment(tb.IsDescription):
    expnum = tb.Int32Col(pos=0)
    xoffset = tb.Float32Col(pos=1)
    yoffset = tb.Float32Col(pos=2)
    xrms = tb.Float32Col(pos=3)
    yrms = tb.Float32Col(pos=4)
    nstars = tb.Float32Col(pos=5)

class NominalVals(tb.IsDescription):
    expnum = tb.Int32Col(pos=0)
    ra = tb.Float32Col(pos=1)
    dec = tb.Float32Col(pos=2)
    parangle = tb.Float32Col(pos=3)
    x_dither_pos = tb.Float32Col(pos=4)
    y_dither_pos = tb.Float32Col(pos=5)
    norm = tb.Float32Col(pos=6)
    als_filename = tb.StringCol((23))


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
                                      'Offset in star matches')
    groupMathes = fileh.create_group(groupAstrometry, 'CatalogMatches', 'Match Catalog Info')

    tableQA = fileh.create_table(groupAstrometry, 'QA', QualityAssessment, 
                             'Qulity Assessment')
    tableNV = fileh.create_table(groupAstrometry,'NominalVals', NominalVals,
                                 'Nominal Values')

    # store shuffle.cfg and DATEvOBS.log files as simple text arrays

    fileshuffle = op.join(args.rootdir, str(args.date) + 'v' + str(args.observation).zfill(3),'shuffle.cfg')
    try:
        f = open(fileshuffle, 'r')
        shuffle = fileh.create_array(groupAstrometry, 'ShuffleCfg', f.read())
        shuffle.set_attr('filename', fileshuffle)
        f.close()
    except:
        args.log.warning('Could not include %s' % fileshuffle)

    logfile = op.join(args.rootdir, str(args.date) + 'v' + str(args.observation).zfill(3),
                      str(args.date) + 'v' + str(args.observation).zfill(3) + '.log')
    try: 
        f = open(logfile)
        fileh.create_array(groupAstrometry, 'LogInfo', f.read())
        f.close()
    except: 
        args.log.warning('Could not include %s' % logfile)

        
    # store fplane table

    filefplane = op.join(args.rootdir, str(args.date) + "v" + str(args.observation).zfill(3),
                         'fplane.txt')    
    try:
        f = ascii.read(filefplane, names=['ifuslot', 'fpx', 'fpy', 'specid',
                                          'specslot', 'ifuid', 'ifurot', 'platesc'])
        fplanetable = fileh.create_table(groupAstrometry, 'fplane', f.as_array())
        fplanetable.set_attr('filename', filefplane)
    except:
        args.log.warning('Could not include %s' % filefplane)

    # store catalogs of stars used in astrometric fit

    file_stars = op.join(args.rootdir, str(args.date) + 'v' + str(args.observation).zfill(3), 
                         'shout.ifustars')    
    try:
        f_stars = ascii.read(file_stars, names=['ignore', 'star_ID', 'ra_cat', 'dec_cat',
                                                'u', 'g', 'r', 'i', 'z'])
        starstable = fileh.create_table(groupAstrometry, 'StarCatalog', f_stars.as_array())
        starstable.set_attr('filename',file_stars)
        if any(f_stars['z'] > 0):
            starstable.set_attr('source','SDSS')
        else:
            starstable.set_attr('source','GAIA')
    except:
        args.log.warning('Could not include %s' % file_stars)
    
    pngfiles = glob.glob(op.join(args.rootdir, str(args.date) + 'v' + str(args.observation).zfill(3), '*.png'))
    pngnum = 1

    for pngfile in pngfiles:
        plt_image = plt.imread(pngfile)
        pngim = fileh.create_array(groupCoadd, 'png_exp' + str(pngnum).zfill(2), plt_image)
        pngim.attrs['CLASS'] = 'IMAGE'
        pngim.attrs['filename'] = pngfile
        pngnum += 1

    fileallmch = op.join(args.rootdir, str(args.date) + 'v' + str(args.observation).zfill(3),'all.mch')
    try:
        allmch = ascii.read(fileallmch)
    except:
        args.log.warning('Could not include %s' % fileallmch)

        
    filenorm = op.join(args.rootdir, str(args.date) + 'v' + str(args.observation).zfill(3),'norm.dat')
    try:
        norm = ascii.read(filenorm)
    except:
        args.log.warning('Could not include %s' % filenorm)

    # index over dithers to gather diher specific info    
    for idx, expn in enumerate(['exp01', 'exp02', 'exp03']):

        radecfile = op.join(args.rootdir, str(args.date) + 'v' + str(args.observation).zfill(3), 
                             'radec2_' + expn + 'dat')
        try:
            rowNV = tableNV.row
            radec = ascii.read(radecfile)
            rowNV['expnum'] = int(expn[3:5])
            rowNV['ra'] = radec['col1']
            rowNV['dec'] = radec['col2']
            rowNV['parangle'] = radec['col3']
            rowNV['norm'] = filenorm[idx]
            rowNV['x_dither_pos'] = allmch['col3'][idx]
            rowNV['y_dither_pos'] = allmch['col4'][idx]
            rowNV['als_filename'] = allmch['col1'][idx]
            rowNV.append()
        except:
            args.log.warning('Could not include NV table info for dither '+ expn)


        fitsfile = op.join(args.rootdir, str(args.date) + 'v' + str(args.observation).zfill(3),
                           str(args.date) + 'v' + str(args.observation).zfill(3)
                           + 'fp_' + expn + '.fits')

        if op.exists(fitsfile):
            
            F = fits.open(fitsfile)
            fitsim = fileh.create_array(groupCoadd, expn, F[0].data)
            fitsim.attrs['CLASS'] = 'IMAGE'
            fitsim.attrs['IMAGE_MINMAXRANGE'] = (-1.5, 100)
            fitsim.attrs['HEADER'] = F[0].header
            fitsim.attrs['filename'] = fitsfile
            F.close()


        # populate offset info for catalog matches
        file_getoff = op.join(args.rootdir, str(args.date) + 'v' + str(args.observation).zfill(3),
                              'getoff_' + expn + '.out')
            
        try:
            f_getoff = ascii.read(file_getoff, names=['xoffset', 'yoffset', 'ra_dex',
                                                      'dec_dex', 'ra_cat', 'dec_cat',
                                                      'ifuslot'])
            getoffinfo = fileh.create_table(groupOffsets, expn, f_getoff.as_array())
            getoffinfo.set_attr('filename', file_getoff)
        except:
            args.log.warning('Could not include %s' % file_getoff)
            

        # populate fiber astrometry data
        file_dith = op.join(args.rootdir, str(args.date) + 'v' + str(args.observation).zfill(3),
                            'dith_' + expn + '.all')    
        try:
            f_dith = ascii.read(file_dith)
            dithinfo = fileh.create_table(groupDithall, expn, f_dith.as_array())
            dithinfo.set_attr('filename', file_dith)
        except:
            args.log.warning('Could not include %s' % file_dith)
            
        # populate median and rms in offsets for quality assessment purposes

        file_getoff2 = op.join(args.rootdir, str(args.date) + 'v'
                               + str(args.observation).zfill(3), 'getoff2_' + expn + '.out')
        try:
            f_getoff2 = ascii.read(file_getoff2)
            row = tableQA.row
            row['expnum'] = int(expn[3:5])
            row['xoffset'] = f_getoff2['col1']
            row['yoffset'] = f_getoff2['col2']
            row['xrms'] = f_getoff2['col3']
            row['yrms'] = f_getoff2['col4']
            row['nstars'] = f_getoff2['col5']
            row.append()
        except:
            args.log.warning('Could not include %s' % file_getoff2)
        
        
        file_xy = op.join(args.rootdir, str(args.date) + 'v'
                          + str(args.observation).zfill(3), 'xy_' + expn + '.dat')
        try:
            xy_table = ascii.read(file_xy)
            fileh.create_table(groupMatches, expn, xy.as_array())
        except:
            args.log.warning('Could not include %s' % file_xy)

    tableQA.set_attr('filename', 'getoff2_exp??.out')
    tableNV.set_attr('dither_file', 'all.mch')
    tableNV.set_attr('norm_file', 'norm.dat')
    tableNV.set_attr('radec_file', 'radec2_exp??.dat')
    tableQA.flush()
    tableNV.flush()
    fileh.close()


if __name__ == '__main__':
    main()
