# -*- coding: utf-8 -*-
"""
Created: 2019/01/23

@author: Erin Mentuch Cooper

Note: requires makemaster.sh to collect info from Karl's
/work/00115/gebhardt/maverick/gettar database

"""

import re
import sys
import os
import os.path as op
import subprocess
import numpy as np
import tables as tb

from astropy.io import ascii
from astropy.io import fits

survey_list = '/work/00115/gebhardt/maverick/gettar/hdr1.objlist' 
path_astrometry = "/work/00115/gebhardt/maverick/vdrp/shifts/"
path_gpinfo = '/work/03946/hetdex/hdr1/calib/fwhm_and_fluxes_better.txt'
mastersci_file = "/work/05350/ecooper/hdr1/survey/mastersci_DR1"
path_radec = "/work/03946/hetdex/hdr1/calib/radec.all"
path_throughput = "/work/00115/gebhardt/maverick/detect/tp/"
path_dithers = '/work/00115/gebhardt/maverick/vdrp/shifts/dithoff/dithall.dat'


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
    ra = tb.Float32Col(pos=4)
    dec = tb.Float32Col(pos=5)
    pa = tb.Float32Col(pos=6)
    trajcra = tb.Float32Col()
    trajcdec = tb.Float32Col()
    trajcpa = tb.Float32Col()
    fwhm = tb.Float32Col(pos=7)
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
    timestamp = tb.StringCol((18),3)
    mjd = tb.Float32Col((3))
    exptime = tb.Float32Col((3))
    darktime = tb.Float32Col((3))


fileh = tb.open_file("survey_test.h5", mode="w", title="Survey test file ")

tableMain = fileh.create_table(fileh.root, 'Survey', Survey,
                               'Main Survey Info')

master_sci = ascii.read(mastersci_file)

master_fwhm = ascii.read(path_gpinfo, data_start=2,
                         names=('datevshot', 'survey', 'fwhm', 'fwhm_mof',
                                'beta', 'n', 'rflux1gc1', 'rflux2gc1',
                                'rflux3gc1', 'rflux1g2', 'rflux2gc2',
                                'rflux3gc2'))

master_astrometry = ascii.read(path_radec)
master_dither = ascii.read(path_dithers, names=['x1', 'y1', 'x2', 'y2',
                                                'x3', 'y3', 'date', 'obsid'])

shotlist = ascii.read(survey_list)

for idx in np.arange(np.size(shotlist)):
    date = shotlist['col1'][idx]
    obsid = shotlist['col2'][idx]
    shotid = int(str(date)+str(obsid).zfill(3))

    sel = np.where( (master_sci['DATE'] == date) & (master_sci['EXP'] == obsid) )
    
    row = tableMain.row

    row['shotid'] = shotid
    row['date'] = str(row['shotid'])[0:8]
    row['obsid'] = str(row['shotid'])[8:11]
    row['datevobs'] = str(row['shotid'])[0:8] + 'v' + str(row['shotid'])[8:11]
   

    row['objid'] = master_sci['OBJECT'][sel][0]
    row['field'] = define_field(row['objid'])

    row['trajcra'] = master_sci['TRAJCRA'][sel][0]
    row['trajcdec'] = master_sci['TRAJCDEC'][sel][0]
    row['trajcpa'] = master_sci['PARANGLE'][sel][0]

    sel2 = np.where((master_astrometry['col1'] == np.int(row['date']))
                   * (master_astrometry['col2'] == np.int(row['obsid'])))

    if np.size(sel2) == 1:
        row['ra'] = master_astrometry['col3'][sel2]
        row['dec'] = master_astrometry['col4'][sel2]
        row['pa'] = master_astrometry['col5'][sel2]
    elif np.size(sel2) > 1:
        print "Weird: More than one RA/DEC/PA value for", row['shotid']
        print "Saving the first RA/DEC/PA value only. Make sure this is ok."
        row['ra'] = master_astrometry['ra'][sel2][0]
        row['dec'] = master_astrometry['dec'][sel2][0]
        row['pa'] = master_astrometry['pa'][sel2][0]
    else:
        row['ra'] = np.nan
        row['dec'] = np.nan
        row['pa'] = np.nan

    datevshot = str(row['date'])+'v'+str(row['obsid']).zfill(3)
    sel1 = np.where(master_fwhm['datevshot'] == datevshot)

    if np.size(sel1) == 1:
        row['fwhm'] = np.float(master_fwhm['fwhm_mof'][sel1])
    elif np.size(sel1) > 1:
        print "Weird: More than one FWHM value for", row['shotid']
        print "Saving the first FWHM value only. Make sure this is ok."
        row['fwhm'] = np.float(master_fwhm['fwhm_mof'][sel1][0])
    else:
        row['fwhm'] = np.nan

    if row['fwhm'] < 0:
        row['fwhm'] = np.nan

    tpfile = op.join(path_throughput, str(row['date']) + 'v'
                     + str(row['obsid']).zfill(3)+'sedtp_f.dat')

    if os.path.isfile(tpfile):
        for line in open(tpfile):
            if '4540' in line:
                tp_4540 = (line.split()[1])
                row['response_4540'] = tp_4540
    else:
        row['response_4540'] = np.nan

    if np.size(sel) == 3:
        
        # add in dither info

        sel3 = np.where((master_dither['date'] == row['date'])
                        * (master_dither['obsid'] == row['obsid']))
        dith_arr = ['x1', 'y1', 'x2', 'y2', 'x3', 'y3']
        x1 = master_dither['x1'][sel3]
        y1 = master_dither['y1'][sel3]
        x2 = master_dither['x2'][sel3]
        y2 = master_dither['y2'][sel3]
        x3 = master_dither['x3'][sel3]
        y3 = master_dither['y3'][sel3]
        
        if np.size(sel) == 1:
            row['xditherpos'] = (x1, x2, x3)
            row['yditherpos'] = (y1, y2, y3)
        
        # add in dither specific info as an array
        expnum_arr = np.full(3, np.nan)
        xoffset_arr = np.full(3, np.nan)
        yoffset_arr = np.full(3, np.nan)
        xrms_arr = np.full(3, np.nan)
        yrms_arr = np.full(3, np.nan)
        nstars_fit_arr = np.full(3, np.nan)
        timestamp_arr = np.chararray(3, itemsize=18)
        mjd_arr = np.full(3, np.nan)
        exptime_arr = np.full(3, np.nan)
        darktime_arr = np.full(3, np.nan)

        for idx2, expn in enumerate(['exp01', 'exp02', 'exp03']):
        
            # add in astrometry offset values
            file_getoff2 = op.join(path_astrometry, str(row['date']) + 'v'
                                + str(row['obsid']).zfill(3), 'getoff2_' + expn + '.out')
            if op.exists(file_getoff2):
                try:
                    f_getoff2 = ascii.read(file_getoff2)
                    expnum_arr[idx2] = int(expn[3:5])
                    xoffset_arr[idx2] = f_getoff2['col1']
                    yoffset_arr[idx2] = f_getoff2['col2']
                    xrms_arr[idx2] = f_getoff2['col3']
                    yrms_arr[idx2] = f_getoff2['col4']
                    nstars_fit_arr[idx2] = f_getoff2['col5'] 
                except:
                    print "Couldn't get data from:",file_getoff2
                    expnum_arr = int(expn[3:5])
                    xoffset_arr = np.full(3, np.nan)
                    yoffset_arr = np.full(3, np.nan)
                    xrms_arr = np.full(3, np.nan)
                    yrms_arr = np.full(3, np.nan)
                    nstars_fit_arr = np.full(3, np.nan)
                    timestamp_arr = np.chararray(3, itemsize=18)
                    mjd_arr = np.full(3, np.nan)
                    exptime_arr = np.full(3, np.nan)
                    darktime_arr = np.full(3, np.nan)

            selexp = np.where( (master_sci['DITHER'] == expn) & 
                               (master_sci['SHOT'] == shotid))

            if np.size(selexp) == 1:
                timestamp_arr[idx2] = master_sci['TIMESTAMP'][selexp][0]
                mjd_arr[idx2] = master_sci['MJD'][selexp][0]
                exptime_arr[idx2] = master_sci['EXPTIME'][selexp][0]
                darktime_arr[idx2] = master_sci['DARKTIME'][selexp][0]
               
        row['expnum'] = expnum_arr
        row['xoffset'] = xoffset_arr
        row['yoffset'] = yoffset_arr
        row['xrms'] = xrms_arr
        row['yrms'] = yrms_arr
        row['nstars_fit'] = nstars_fit_arr
        
        row['timestamp'] = timestamp_arr
        row['mjd'] = mjd_arr
        row['exptime'] = exptime_arr
        row['darktime'] = darktime_arr
                
    row.append()

tableMain.flush()
fileh.close()
