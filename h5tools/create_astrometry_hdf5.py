# -*- coding: utf-8 -*-                                                                                                                                                                
""" 
@author: Erin Mentuch Cooper

This script can add Astromtry group info and data to
an exising HDF5 file or create a new HDF5 file

Example command line use:

python create_astrometry_hdf5.py -d 20181111 -o 15 -of astrometry_20181111v015.h5

or 

python create_astrometry_hdf5.py -d 20181111 -o 15 -of 20181111v015.h5 --append

"""

import glob
import re
import os
import traceback

import tables as tb
import argparse as ap
import os.path as op
import numpy as np

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from hetdex_api.input_utils import setup_logging
from hetdex_api.config import HDRconfig

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

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
    relflux_virus = tb.Float32Col(pos=6)
    als_filename = tb.StringCol((23))


class Dithall(tb.IsDescription):
    ra = tb.Float32Col(pos=0)
    dec = tb.Float32Col(pos=1)
    ifuslot = tb.StringCol((3))
    XS = tb.Float32Col()
    YS = tb.Float32Col()
    xfplane = tb.Float32Col()
    yfplane = tb.Float32Col()
    multifits = tb.StringCol((28))
    timestamp = tb.StringCol((17))
    exposure = tb.StringCol((5))


class StarCatalog(tb.IsDescription):
    ignore = tb.Int64Col(pos=0)
    star_ID = tb.StringCol((28), pos=1)
    ra_cat = tb.Float64Col(pos=2)
    dec_cat = tb.Float64Col(pos=3)
    u = tb.Float64Col(pos=4)
    g = tb.Float64Col(pos=5)
    r = tb.Float64Col(pos=6)
    i = tb.Float64Col(pos=7)
    z = tb.Float64Col(pos=8)
    
    
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
                        type=str, default='/scratch/00115/gebhardt/vdrp/shifts')

    parser.add_argument('-of', '--outfilename', type=str,
                        help='''Relative or absolute path for output HDF5
                        file.''', default=None)

    parser.add_argument('-a', '--append',
                        help='''Appending to existing file.''',
                        action="count", default=0)

    parser.add_argument("-tp", "--tpdir",
                        help='''Directory for Throughput Info''',
                        type=str,
                       default='/scratch/00115/gebhardt/detect')

    parser.add_argument("-detdir", "--detectdir",
                        help='''Directory for Detect Info''',
                        type=str,
                        default='/scratch/projects/hetdex/detect')
    
    parser.add_argument("-survey", "--survey",
                        help="""{hdr1, hdr2, hdr2.1, hdr3, hdr4, hdr5}""",
                        type=str, default="hdr4")
    
    
    args = parser.parse_args(argv)
    args.log = setup_logging()

    # Creates a new file if the "--append" option is not set or the file                                              
    # does not already exist.
    does_exist = False
    if op.exists(args.outfilename) and args.append:
        args.log.info('Appending astrometry to %s' % args.outfilename)
        fileh = tb.open_file(args.outfilename, 'a')
        does_exist = True
        try:
            fileh.remove_node(fileh.root.Astrometry, recursive=True)
        except:
            pass
    else:
        args.log.info('Creating new file for astrometry %s' % args.outfilename)
        fileh = tb.open_file(args.outfilename, 'w')
        
    groupAstrometry = fileh.create_group(fileh.root, 'Astrometry', 'Astrometry Info')
    groupCoadd = fileh.create_group(groupAstrometry, 'CoaddImages', 'Coadd Images')
    groupDithall = fileh.create_group(groupAstrometry, 'Dithall', 'Fiber Astrometry Info')
    groupOffsets = fileh.create_group(groupAstrometry, 'PositionOffsets', 
                                      'Offset in star matches')
    groupMatches = fileh.create_group(groupAstrometry, 'CatalogMatches', 'Match Catalog Info')

    tableQA = fileh.create_table(groupAstrometry, 'QA', QualityAssessment, 
                             'Quality Assessment')
    tableNV = fileh.create_table(groupAstrometry,'NominalVals', NominalVals,
                                 'Nominal Values')

    datevshot = str(args.date) + 'v' + str(args.observation).zfill(3)
    shotid = int(str(args.date) + str(args.observation).zfill(3))

    #check if shotid is in badlist

    config = HDRconfig(args.survey)
    badshots = np.loadtxt(config.badshot, dtype=int)

    badshotflag = False

    if shotid in badshots:
        badshotflag = True
    
    # store shuffle.cfg and DATEvOBS.log files

    shiftsdir = op.join(args.rootdir, str(args.date) + 'v' + str(args.observation).zfill(3))
    fileshuffle = op.join(shiftsdir, 'shuffle.cfg')
    try:
        f = open(fileshuffle, 'r')
        shuffle = fileh.create_array(groupAstrometry, 'ShuffleCfg', f.read().encode())
        shuffle.set_attr('filename','shuffle.cfg')
        f.close()
    except: #alternate location
        shiftsdir = op.join(args.rootdir, "vdrp")
        fileshuffle = op.join(shiftsdir, 'shuffle.cfg')
        try:
            f = open(fileshuffle, 'r')
            shuffle = fileh.create_array(groupAstrometry, 'ShuffleCfg', f.read().encode())
            shuffle.set_attr('filename', 'shuffle.cfg')
            f.close()
        except:
            args.log.warning('Could not include %s' % fileshuffle)

    # store fplane table

    filefplane = op.join(shiftsdir, 'fplane.txt')
    try:
        f = ascii.read(filefplane, names=['ifuslot', 'fpx', 'fpy', 'specid',
                                          'specslot', 'ifuid', 'ifurot', 'platesc'])
        fplanetable = fileh.create_table(groupAstrometry, 'fplane', f.as_array())
        fplanetable.set_attr('filename', 'fplane.txt')
    except:
        args.log.warning('Could not include %s' % filefplane)

    # store catalogs of stars used in astrometric fit

    file_stars = op.join(shiftsdir, 'shout.ifustars')

    try:
        f_stars = Table.read(file_stars,
                             names=['ignore', 'star_ID', 'ra_cat', 'dec_cat',
                                    'u', 'g', 'r', 'i', 'z'], format='ascii')

        starstable = fileh.create_table(groupAstrometry, 'StarCatalog', StarCatalog, 'StarCatalog')

        for row in f_stars:
            if row['g'] > 22:
                continue
            starrow = starstable.row
            
            starrow['star_ID'] = row['star_ID']
            for col in ['ignore','ra_cat', 'dec_cat', 'u', 'g', 'r', 'i', 'z']:
                starrow[col] = row[col]
            starrow.append()
            
        
        starstable.set_attr('filename', 'DATEvOBS.ifu')
        if any(f_stars['z'] > 0):
            starstable.set_attr('catalog', 'SDSS')
        else:
            starstable.set_attr('catalog', 'GAIA')
    except:
        args.log.warning('Could not include %s' % file_stars)
    
    pngfiles = glob.glob(op.join(shiftsdir, '2*.png'))
    pngnum = 1

    for pngfile in pngfiles:
        plt_image = plt.imread(pngfile)
        pngim = fileh.create_array(groupCoadd, 'png_exp' + str(pngnum).zfill(2), plt_image)
        pngim.attrs['CLASS'] = 'IMAGE'
        pngim.attrs['filename'] = pngfile
        pngnum += 1

    fileallmch = op.join(shiftsdir, 'all.mch')
    try:
        allmch = ascii.read(fileallmch)
    except:
        args.log.warning('Could not include %s' % fileallmch)

    filenorm = op.join(shiftsdir, 'norm.dat')
    if not op.exists(filenorm):
        #try alternate location
        filenorm = op.join(args.detectdir, 'norm.dat')

    try:
        norm = Table.read(filenorm, format='ascii.no_header')
    except:
        args.log.warning('Could not include (table) %s' % filenorm)
        args.log.error(traceback.format_exc())

    fns = glob.glob(f"match_pngs/match_*.png")
    exps = [f"exp{str(x).zfill(2)}" for x in np.arange(1,len(fns)+1,1)]


    # index over dithers to gather diher specific info
    for idx, expn in enumerate(exps):
    #for idx, expn in enumerate(['exp01', 'exp02', 'exp03']):

        radecfile = op.join(shiftsdir, 'radec2_' + expn + '.dat')
        rowNV = tableNV.row
        try:
            radec = ascii.read(radecfile)
            rowNV['expnum'] = int(expn[3:5])
            rowNV['ra'] = radec['col1'][0]
            rowNV['dec'] = radec['col2'][0]
            rowNV['parangle'] = radec['col3'][0]
        except:
            args.log.warning('Could not include %s' % radecfile)

        try:
#            sel_datevobs = norm['col1'] == str(args.date) + 'v' + str(args.observation).zfill(3)
            if idx == 0:
                rowNV['relflux_virus'] = norm['col1'][0]
            elif idx == 1:
                rowNV['relflux_virus'] = norm['col2'][0]
            elif idx == 2:
                rowNV['relflux_virus'] = norm['col3'][0]
        except Exception:
            args.log.error(traceback.format_exc())
            args.log.warning('Could not include norm.dat')
        
        try:
            rowNV['x_dither_pos'] = allmch['col3'][idx]
            rowNV['y_dither_pos'] = allmch['col4'][idx]
            rowNV['als_filename'] = allmch['col1'][idx]
        except:
            args.log.warning('Could not include %s' % fileallmch)
        
        rowNV.append()

        fitsfile = op.join(shiftsdir,
                           str(args.date) + 'v' + str(args.observation).zfill(3)
                           + 'fp_' + expn + '.fits')

        if op.exists(fitsfile):
            
            F = fits.open(fitsfile)
            fitsim = fileh.create_array(groupCoadd, expn, F[0].data)
            fitsim.attrs['CLASS'] = 'IMAGE'
            fitsim.attrs['IMAGE_MINMAXRANGE'] = (-1.5, 100)
            fitsim.attrs['HEADER'] = F[0].header
            fitsim.attrs['filename'] = 'DATEvOBSfp_exp??.fits'
            F.close()


        matchpdf = op.join(shiftsdir,
                           'match_' + expn + '.pdf')

        matchpng_orig = op.join( shiftsdir, 'match_{}.png'.format(expn) )
        
        matchpng = 'match_pngs/match_'+ str(args.date) + 'v' + str(args.observation).zfill(3) + '_' + expn + '.png'
        
        if op.exists(matchpdf):
            os.system('convert ' + matchpdf + ' ' + matchpng)  
            plt_matchim = plt.imread(matchpng)
            matchim = fileh.create_array(groupCoadd, 'match_' + expn, plt_matchim)
            matchim.attrs['CLASS'] = 'IMAGE'
            matchim.attrs['filename'] = matchpdf
        elif op.exists( matchpng_orig):
            os.system('cp ' + matchpng_orig + ' ' + matchpng)
            plt_matchim = plt.imread(matchpng)
            matchim = fileh.create_array(groupCoadd, 'match_' + expn, plt_matchim)
            matchim.attrs['CLASS'] = 'IMAGE'
            matchim.attrs['filename'] = matchpng_orig
        else:
            args.log.warning('Could not include %s' % matchpng)

        # populate offset info for catalog matches
        file_getoff = op.join(shiftsdir,
                              'getoff_' + expn + '.out')

        try:
            f_getoff = ascii.read(file_getoff, names=['xoffset', 'yoffset', 'ra_dex',
                                                      'dec_dex', 'ra_cat', 'dec_cat',
                                                      'ifuslot'])
            getoffinfo = fileh.create_table(groupOffsets, expn, f_getoff.as_array())
            getoffinfo.set_attr('filename', 'getoff_exp??.out')
        except:
            args.log.warning('Could not include %s' % file_getoff)
            
        # populate fiber astrometry data
        file_dith = op.join(shiftsdir,
                            'dith_' + expn + '.all')    
        try:
            f_dith = ascii.read(file_dith)
            dithtab = fileh.create_table(groupDithall, expn, Dithall)
            
            for f_dith_row in f_dith:
                dithrow = dithtab.row
                
                dithrow['ra'] = f_dith_row['ra']
                dithrow['dec'] = f_dith_row['dec']
                dithrow['ifuslot'] = f_dith_row['ifuslot']
                dithrow['XS'] = f_dith_row['XS']
                dithrow['YS'] = f_dith_row['YS']
                dithrow['xfplane'] = f_dith_row['xfplane']
                dithrow['yfplane'] = f_dith_row['yfplane']
                dithrow['multifits'] = f_dith_row['multifits']
                dithrow['timestamp'] = f_dith_row['timestamp']
                dithrow['exposure'] = f_dith_row['exposure']
                dithrow.append()
           
            dithtab.set_attr('filename', file_dith)
            dithtab.flush()
        except:
            args.log.warning('Could not include %s' % file_dith)
            
        # populate median and rms in offsets for quality assessment purposes

        file_getoff2 = op.join(shiftsdir, 'getoff2_' + expn + '.out')
        try:
            f_getoff2 = ascii.read(file_getoff2)
            row = tableQA.row
            row['expnum'] = int(expn[3:5])
            row['xoffset'] = f_getoff2['col1'][0]
            row['yoffset'] = f_getoff2['col2'][0]
            row['xrms'] = f_getoff2['col3'][0]
            row['yrms'] = f_getoff2['col4'][0]
            row['nstars'] = f_getoff2['col5'][0]
            row.append()

        except:
            args.log.warning('Could not include %s' % file_getoff2)
        
        
        file_xy = op.join(shiftsdir, 'xy_' + expn + '.dat')
        try:
            xy_table = ascii.read(file_xy)
            tableXY = fileh.create_table(groupMatches, expn, xy_table.as_array())
            tableXY.set_attr('filename', 'xy_exp??.dat')
        except:
            args.log.warning('Could not include %s' % file_xy)
            
    tableQA.set_attr('filename', 'getoff2_exp??.out')
    tableNV.set_attr('dither_file', 'all.mch')
    tableNV.set_attr('norm_file', 'norm.dat')
    tableNV.set_attr('radec_file', 'radec2_exp??.dat')
    tableQA.flush()
    tableNV.flush()

    tableQA = fileh.root.Astrometry.QA
    tableNV = fileh.root.Astrometry.NominalVals

    try:
        radecfinalfile = op.join(shiftsdir, 'radec2_final.dat')
        radectab = ascii.read(radecfinalfile, names=['ra','dec','pa'])

    except:

        if badshotflag:
            args.log.warning('Could not open %s' % radecfinalfile)
        else:
            args.log.error('Could not open %s' % radecfinalfile)
        
    shottable = fileh.root.Shot

    for shot in shottable:
        if op.exists(radecfinalfile):
            shot['ra'] = radectab['ra'][0]
            shot['dec'] = radectab['dec'][0]
            shot['pa'] = radectab['pa'][0]
        else:
            if badshotflag:
                args.log.warning('Could not open %s' % radecfinalfile)
            else:
                args.log.error('Could not open %s' % radecfinalfile)

        try:
            if len(tableQA.cols.xoffset[:]) == 3: #might now be a single exposure, in which case there are no offsets
                shot['xoffset'] = tableQA.cols.xoffset[:]
                shot['yoffset'] = tableQA.cols.yoffset[:]
                shot['xrms'] = tableQA.cols.xrms[:]
                shot['yrms'] = tableQA.cols.yrms[:]
                shot['nstars_fit'] = tableQA.cols.nstars[:]
            else:
                args.log.info(f'Could not include astrometry QA info (cnd 1a) for {datevshot}. Shape mismatch. '
                                 f'len = {len(tableQA.cols.xoffset[:])}. May be single exposure?')
        except:
            if badshotflag:
                args.log.warning('Could not include astrometry QA info (cnd 1) for %s' % datevshot)
            else:
                args.log.error('Could not include astrometry shot info (cnd 2) for %s' % datevshot)

            args.log.warning(traceback.format_exc())
        try:
            if len(tableNV.cols.x_dither_pos[:]) == 3:
                shot['xditherpos'] = tableNV.cols.x_dither_pos[:]
                shot['yditherpos'] = tableNV.cols.y_dither_pos[:]
            else:
                args.log.info(f'Could not include astrometry Norm info, x & y_dither_pos (cnd 3a) for {datevshot}. Shape mismatch. '
                                 f'len = {len(tableNV.cols.x_dither_pos[:])}. May be single exposure?')
        except:
            args.log.warning('Could not include astrometry shot info (cnd 3) for %s' % datevshot)
            args.log.warning(traceback.format_exc())

        try:
            if len(tableNV.cols.relflux_virus[:]) == 3:
                shot['relflux_virus'] = tableNV.cols.relflux_virus[:]
            else:
                args.log.info(f'Could not include astrometry Norm info, relflux_virus, (cnd 3a) for {datevshot}. Shape mismatch. '
                                 f'len = {len(tableNV.cols.relflux_virus[:])}. May be single exposure?')
        except:
            if badshotflag:
                args.log.warning('Could not include relflux_virus info (cnd 1) for %s' % datevshot)
            else:
                args.log.error('Could not include relflux_virus info (cnd 2) for %s' % datevshot)
            args.log.warning(traceback.format_exc())
                
        shot.update()

    fileh.close()

if __name__ == '__main__':
    main()
