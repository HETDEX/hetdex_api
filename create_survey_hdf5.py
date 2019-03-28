# -*- coding: utf-8 -*-
"""
Created: 2019/01/23

@author: Erin Mentuch Cooper

Script to gather information for the Survey Table. This provides
a master look up table for representaive values for each shot
in the HETDEX survey

The script uses the files to generate the Table, but any
list of DATE OBS would work.

/work/03946/hetdex/hdr1/reduction/hdr1.scilist 
/work/03946/hetdex/hdr1/reduction/hdr1.callist 

Requires multi*fits files to gather header info, a shifts/
directory to gather some QA astrometry data, and additional
paths contained in the HETDEX_API/config.py files


To run:

pyton create_survey_hdf5.py -of survey_hdr1.h5

"""

import re
import sys
import glob
import os
import os.path as op
import subprocess
import numpy as np
import tables as tb
import argparse as ap

from astropy.io import ascii
from astropy.io import fits
from astropy.table import vstack, Table
from input_utils import setup_logging

import config


def get_files(path_to_multifits, date, obsid):
    files = glob.glob(op.join(path_to_multifits, str(date), 'virus',
                              'virus%07d' % int(obsid),
                              'exp*', 'virus', 'multi_*.fits'))
    return files

def get_files_exp(path_to_multifits, date, obsid, expn):
    files = glob.glob(op.join(path_to_multifits, str(date), 'virus',
                              'virus%07d' % int(obsid),
                              expn, 'virus', 'multi_*.fits'))
    return files


def define_field(objname):
    if re.match('par', str(objname)):
        field = 'parallel'
    elif re.match('COS|cos|DEXcos', str(objname)):
        field = 'cosmos'
    elif re.match('EGS', str(objname)):
        field = 'egs'
    elif re.match('GN', str(objname)):
        field = 'goods-n'
    elif re.match('DEX0|DEXfl', str(objname)):
        field = 'dex-fall'
    elif re.match('HS|DEXsp', str(objname)):
        field = 'dex-spring'
    else:
        field = 'other'

    return field


class Survey(tb.IsDescription):
    shotid = tb.Int64Col()
    date = tb.Int32Col(pos=1)
    obsid = tb.Int32Col(pos=2)
    objid = tb.StringCol((18), pos=3)
    field = tb.StringCol((12), pos=0)
    ra = tb.Float64Col(pos=4)
    dec = tb.Float64Col(pos=5)
    pa = tb.Float64Col(pos=6)
    n_ifu = tb.Int32Col()
    trajcra = tb.Float32Col()
    trajcdec = tb.Float32Col()
    trajcpa = tb.Float32Col()
    fwhm_flag = tb.Int32Col(pos=9)
    fwhm_gaussian = tb.Float32Col(pos=10)
    fwhm_moffat = tb.Float32Col(pos=11)
    moffat_beta = tb.Float32Col(pos=12)
    RelFlux1 = tb.Float32Col(pos=13)
    RelFlux2 = tb.Float32Col(pos=14)
    RelFlux3 = tb.Float32Col(pos=15)
    response_4540 = tb.Float32Col(pos=8)  # normalized for 360s
    xditherpos = tb.Float32Col((3))
    yditherpos = tb.Float32Col((3))
    xoffset = tb.Float32Col((3))
    yoffset = tb.Float32Col((3))
    xrms = tb.Float32Col((3))
    yrms = tb.Float32Col((3))
    nstars_fit = tb.Int32Col((3))
    datevobs = tb.StringCol((12))
    expnum = tb.Int32Col((3))
    mjd = tb.Float32Col((3))
    exptime = tb.Float32Col((3))
    darktime = tb.Float32Col((3))
    ra_flag = tb.StringCol((1),3)


def main(argv=None):
    ''' Main Function '''
    # Call initial parser from init_utils
    parser = ap.ArgumentParser(description="""Create HDF5 Astrometry file.""",
                               add_help=True)

    parser.add_argument("-r", "--rootdir",
                        help='''Root Directory for Reductions''',
                        type=str, default='/work/03946/hetdex/maverick/red1/reductions/')

    parser.add_argument("-ad", "--astrometry_dir",
                        help='''Directory for Shifts''',
                        type=str, default='/work/00115/gebhardt/maverick/vdrp/shifts/')

    parser.add_argument('-of', '--outfilename', type=str,
                        help='''Relative or absolute path for output HDF5
                        file.''', default=None)

    args = parser.parse_args(argv)
    args.log = setup_logging()

    fileh = tb.open_file(args.outfilename, mode="w", title="HDR1 Survey file ")
    
    tableMain = fileh.create_table(fileh.root, 'Survey', Survey,
                                   'Main Survey Info')

    master_fwhm = ascii.read('/work/03030/grnagara/maverick/getgp2/DR1FWHM.txt')
#    master_fwhm = ascii.read(config.path_gpinfo)

    master_astrometry = ascii.read(config.path_radec)
    ra_flag_table = ascii.read(config.path_acc_flags, names=['date','obsid','exp01','exp02','exp03'])

    survey_table = Table(np.loadtxt(config.survey_list, dtype=np.int, usecols=[0,1]), names=['date','obs'])
    cal_table = Table(np.loadtxt(config.cal_list, dtype=np.int, usecols=[0,1]), names=['date','obs'])
    shotlist = vstack([survey_table, cal_table])

    for idx in np.arange(np.size(shotlist)):
        
        date = shotlist['date'][idx]
        obsid = shotlist['obs'][idx]
        datevshot = str(date) + 'v' + str(obsid).zfill(3)

        # grabbing multi fits files to get header info 
        files = get_files(args.rootdir, date, obsid)

        shotid = int(str(date)+str(obsid).zfill(3))
        
        row = tableMain.row
        row['shotid'] = shotid
        row['date'] = str(row['shotid'])[0:8]
        row['obsid'] = str(row['shotid'])[8:11]
        row['datevobs'] = str(row['shotid'])[0:8] + 'v' + str(row['shotid'])[8:11]
        
        F = fits.open(files[0])
        row['objid'] = F[0].header['OBJECT']
        row['field'] = define_field(row['objid'])

        # call calibration shots 'cal'
        if np.any( (cal_table['date'] == date) & (cal_table['obs'] == obsid)):
            row['field'] = 'cal'

        row['trajcra'] = F[0].header['TRAJCRA']
        row['trajcdec'] = F[0].header['TRAJCDEC']
        row['trajcpa'] = F[0].header['PARANGLE']

        sel2 = np.where((master_astrometry['col1'] == np.int(row['date']))
                        * (master_astrometry['col2'] == np.int(row['obsid'])))

        if np.size(sel2) == 1:
            row['ra'] = master_astrometry['col3'][sel2]
            row['dec'] = master_astrometry['col4'][sel2]
            row['pa'] = master_astrometry['col5'][sel2]
        elif np.size(sel2) > 1:
            args.log.warning('Weird: More than one RA/DEC/PA value for %s' % datevshot)
            row['dec'] = master_astrometry['col4'][sel2][0]
            row['pa'] = master_astrometry['col5'][sel2][0]
        else:
            args.log.warning('Missing RA/DEC/PA for %s' % datevshot)
            row['ra'] = np.nan
            row['dec'] = np.nan
            row['pa'] = np.nan

        sel1 = np.where(master_fwhm['Datevshot'] == datevshot)
        
        if np.size(sel1) == 1:
            row['fwhm_flag'] = np.int(master_fwhm['Good(1)/Approximated(0)'][sel1])
            row['fwhm_moffat'] = np.float(master_fwhm['Moffat_FWHM'][sel1])
            row['fwhm_gaussian'] = np.float(master_fwhm['Gaussian_FWHM'][sel1])
            row['moffat_beta'] = np.float(master_fwhm['Moffat_Beta'][sel1])
            row['RelFlux1'] = np.float(master_fwhm['RelFlux1'][sel1])
            row['RelFlux2'] = np.float(master_fwhm['RelFlux2'][sel1])
            row['RelFlux3'] = np.float(master_fwhm['RelFlux3'][sel1])

        elif np.size(sel1) > 1:
            args.log.warning('Weird: More than one FWHM value for %s' & datevshot)
        else:
            if row['field'] != 'cal':
                args.log.warning('Could not find FWHM for %s' % datevshot)

        tpfile = op.join(config.tp_dir, str(row['date']) + 'v'
                         + str(row['obsid']).zfill(3)+'sedtp_f.dat')
        
        if op.exists(tpfile):
            tp_tab = ascii.read(tpfile)
            tp_4540 = tp_tab['col2'][np.where(tp_tab['col1'] == 4540.)][0]
            row['response_4540'] = tp_4540
        else:
            row['response_4540'] = 0.0
            if row['field'] != 'cal':
                args.log.warning('Could not get response_4540 for %s' % datevshot)
    
        # add RA acceptance flag data
        sel_acc = np.where((ra_flag_table['date'] == row['date']) & (ra_flag_table['obsid'] == row['obsid']))
        if (np.size(sel_acc) > 0):
            row['ra_flag'] = [ra_flag_table['exp01'][sel_acc][0], ra_flag_table['exp02'][sel_acc][0], 
                              ra_flag_table['exp03'][sel_acc][0]]
        else:
            row['ra_flag'] =['','','']
            if row['field'] != 'cal':
                args.log.warning('Could not get RA flags from %s' % datevshot)
        
        allmchfile = op.join(args.astrometry_dir, str(row['date']) + 'v'
                             + str(row['obsid']).zfill(3), 'all.mch')
        if op.exists(allmchfile):
            try: 
                allmch = ascii.read(allmchfile)
                row['xditherpos'] = allmch['col3']
                row['yditherpos'] = allmch['col3']
            except:
                args.log.warning('Could not ingest %s' % allmchfile)
        else:
            args.log.warning('Could not open %s' % allmchfile)
            
        # add in additional dither specific info as an array                                                                          
        expnum_arr = np.zeros(3, dtype=int)
        xoffset_arr = np.full(3, np.nan)
        yoffset_arr = np.full(3, np.nan)
        xrms_arr = np.full(3, np.nan)
        yrms_arr = np.full(3, np.nan)
        nstars_fit_arr = np.full(3, np.nan)
        mjd_arr = np.zeros(3)
        exptime_arr = np.zeros(3)
        darktime_arr = np.zeros(3)

        for idx2, expn in enumerate(['exp01','exp02', 'exp03']):

            files_exp = get_files_exp(args.rootdir, date, obsid, expn)
            
            if expn == 'exp01':
                try:
                    row['n_ifu'] = np.size(files_exp)/4
                except:
                    args.log.warning('Could not get N_ifu for %s' % datevshot)

            if np.size(files_exp) > 0:
                F_exp = fits.open(files_exp[0])
                mjd_arr[idx2] = F_exp[0].header['MJD']
                exptime_arr[idx2] = F_exp[0].header['EXPTIME']
                darktime_arr[idx2] = F_exp[0].header['DARKTIME']
            else:
                args.log.info('Missing some dither info for  %s' % datevshot)
        
            file_getoff2 = op.join(args.astrometry_dir, str(row['date']) + 'v'
                                       + str(row['obsid']).zfill(3), 'getoff2_' + expn + '.out')
            if op.exists(file_getoff2):
                f_getoff2 = ascii.read(file_getoff2)
                if f_getoff2['col5'][0] > 0:
                    expnum_arr[idx2]    = int(expn[3:5])
                    xoffset_arr[idx2] = f_getoff2['col1']
                    yoffset_arr[idx2] = f_getoff2['col2']
                    xrms_arr[idx2] = f_getoff2['col3']
                    yrms_arr[idx2] = f_getoff2['col4']
                    nstars_fit_arr[idx2] = f_getoff2['col5'] 
                else:
                    args.log.warning('Could not get data from %s' % file_getoff2)
                    expnum_arr[idx2] = int(expn[3:5])
                    xoffset_arr[idx2] = np.nan
                    yoffset_arr[idx2] = np.nan
                    xrms_arr[idx2] = np.nan
                    yrms_arr[idx2] = np.nan
                    nstars_fit_arr[idx2] = 0
            else: 
                args.log.warning('Could not open %s' % file_getoff2)

        row['expnum'] = expnum_arr
        row['xoffset'] = xoffset_arr
        row['yoffset'] = yoffset_arr
        row['xrms'] = xrms_arr
        row['yrms'] = yrms_arr
        row['nstars_fit'] = nstars_fit_arr
        row['mjd'] = mjd_arr
        row['exptime'] = exptime_arr
        row['darktime'] = darktime_arr

        row.append()

    tableMain.flush()
    fileh.close()


if __name__ == '__main__':
    main()
