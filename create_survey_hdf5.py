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
from astropy.table import vstack, Table

mastersci_file = "/work/05350/ecooper/hdr1/survey/mastersci_DR1"

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
    ra_flag = tb.StringCol((1),3)

fileh = tb.open_file("survey_test.h5", mode="w", title="Survey test file ")

tableMain = fileh.create_table(fileh.root, 'Survey', Survey,
                               'Main Survey Info')

master_sci = ascii.read(mastersci_file)
master_fwhm = ascii.read(path_gpinfo)
master_astrometry = ascii.read(path_radec)
ra_flag_table = ascii.read(path_accflage, names=['date','obsid','exp01','exp02','exp03'])


survey_table = Table(np.loadtxt(survey_list, dtype=np.int, usecols=[0,1]), names=['date','obs'])
cal_table = Table(np.loadtxt(cal_list, dtype=np.int, usecols=[0,1]), names=['date','obs'])
shotlist = vstack([survey_table, cal_table])

for idx in np.arange(np.size(shotlist)):
    
    date = shotlist['date'][idx]
    obsid = shotlist['obs'][idx]

    shotid = int(str(date)+str(obsid).zfill(3))

    sel = np.where( (master_sci['DATE'] == date) & (master_sci['EXP'] ==obsid))
    
    row = tableMain.row

    row['shotid'] = shotid
    row['date'] = str(row['shotid'])[0:8]
    row['obsid'] = str(row['shotid'])[8:11]
    row['datevobs'] = str(row['shotid'])[0:8] + 'v' + str(row['shotid'])[8:11]
   

    row['objid'] = master_sci['OBJECT'][sel][0]
    row['field'] = define_field(row['objid'])
    # call calibration shots 'cal'
    if np.any( (cal_table['date'] == date) & (cal_table['obs'] == obsid)):
        row['field'] = 'cal'

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
        row['ra'] = master_astrometry['col3'][sel2][0]
        row['dec'] = master_astrometry['col4'][sel2][0]
        row['pa'] = master_astrometry['col5'][sel2][0]
    else:
        print "No link to file in master_sci for ", datevshot
        row['ra'] = np.nan
        row['dec'] = np.nan
        row['pa'] = np.nan

    datevshot = str(row['date'])+'v'+str(row['obsid']).zfill(3)
    sel1 = np.where(master_fwhm['Datevshot'] == datevshot)

    if np.size(sel1) == 1:
        row['fwhm'] = np.float(master_fwhm['FWHM_Moffat'][sel1])
    elif np.size(sel1) > 1:
        print "Weird: More than one FWHM value for", row['shotid']
        print "Saving the first FWHM value only. Make sure this is ok."
        row['fwhm'] = np.float(master_fwhm['fwhm_mof'][sel1][0])
    else:
        row['fwhm'] = np.nan
        print "Couldn't find FWHM for", datevshot

    tpfile = op.join(tp_dir, str(row['date']) + 'v'
                     + str(row['obsid']).zfill(3)+'sedtp_f.dat')

    if os.path.isfile(tpfile):
        for line in open(tpfile):
            if '4540' in line:
                tp_4540 = (line.split()[1])
                row['response_4540'] = tp_4540
    else:
        row['response_4540'] = np.nan

    
    # add RA acceptance flag data
    sel_acc = np.where((ra_flag_table['date'] == row['date']) & (ra_flag_table['obsid'] == row['obsid']))
    if (np.size(sel_acc) > 0):
        row['ra_flag'] = [ra_flag_table['exp01'][sel_acc][0], ra_flag_table['exp02'][sel_acc][0], 
                          ra_flag_table['exp03'][sel_acc][0]]
    else:
        row['ra_flag'] =['','','']
        if row['field'] != 'cal':
            print "Couldn't get RA flags from:", row['datevobs']

        
    if np.size(sel) == 3: 
        
        # add in dither info
        allmchfile = op.join(astrometry_dir, str(row['date']) + 'v'
                             + str(row['obsid']).zfill(3), 'all.mch')
        if op.exists(allmchfile):
            allmch = ascii.read(allmchfile)
            row['xditherpos'] = allmch['col3']
            row['yditherpos'] = allmch['col3']
        else: 
            print "Could not open: ", allmchfile 
            
        # add in dither specific info as an array
        expnum_arr = np.zeros(3, dtype=int)
        xoffset_arr = np.full(3, np.nan)
        yoffset_arr = np.full(3, np.nan)
        xrms_arr = np.full(3, np.nan)
        yrms_arr = np.full(3, np.nan)
        nstars_fit_arr = np.full(3, np.nan)
        timestamp_arr = np.chararray(3, itemsize=18)
        mjd_arr = np.full(3, np.nan)
        exptime_arr = np.full(3, np.nan)
        darktime_arr = np.full(3, np.nan)

        for idx2, expn in enumerate(['exp01','exp02', 'exp03']):
        
            # add in astrometry offset values
            file_getoff2 = op.join(astrometry_dir, str(row['date']) + 'v'
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
                    print "Couldn't get data from:",file_getoff2
                    expnum_arr[idx2] = int(expn[3:5])
                    xoffset_arr[idx2] = np.nan
                    yoffset_arr[idx2] = np.nan
                    xrms_arr[idx2] = np.nan
                    yrms_arr[idx2] = np.nan
                    nstars_fit_arr[idx2] = 0
            else: 
                print "could not open getoff2_exp??.out file"

            selexp = np.where( (master_sci['DITHER'] == expn) & 
                               (master_sci['SHOT'] == shotid))

            if np.size(selexp) == 1:
                timestamp_arr[idx2] = master_sci['TIMESTAMP'][selexp][0]
                mjd_arr[idx2] = master_sci['MJD'][selexp][0]
                exptime_arr[idx2] = master_sci['EXPTIME'][selexp][0]
                darktime_arr[idx2] = master_sci['DARKTIME'][selexp][0]
            else:
                print "Couldn't get info from master_sci for," [date, obsid, expn]

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
