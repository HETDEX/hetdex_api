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

path_astrometry = "/work/00115/gebhardt/maverick/vdrp/shifts/"
path_gpinfo = '/work/03946/hetdex/hdr1/calib/fwhm_and_fluxes_better.txt'
mastersci_file = "/work/05350/ecooper/maverick/gettar/mastersci_DR1"
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
    x1 = tb.Float32Col()
    y1 = tb.Float32Col()
    x2 = tb.Float32Col()
    y2 = tb.Float32Col()
    x3 = tb.Float32Col() 
    y3 = tb.Float32Col()
    datevobs = tb.StringCol((12))

class Exp(tb.IsDescription):
    shotid = tb.Int64Col(pos=0)
    date = tb.Int32Col(pos=1)
    obsid = tb.Int32Col(pos=2)
    expnum = tb.Int32Col(pos=3)
    timestamp = tb.StringCol((18), pos=4)
    mjd = tb.Float32Col(pos=5)
    exptime = tb.Float32Col((7), pos=6)
    darktime = tb.Float32Col(pos=7)


fileh = tb.open_file("survey_test.h5", mode="w", title="Survey test file ")

tableMain = fileh.create_table(fileh.root, 'Survey', Survey,
                               'Main Survey Info')
tableExp = fileh.create_table(fileh.root, 'Exp', Exp,
                              'Exposure Specific Info')

master_sci = ascii.read(mastersci_file)

master_fwhm = ascii.read(path_gpinfo, data_start=2,
                         names=('datevshot', 'survey', 'fwhm', 'fwhm_mof',
                                'beta', 'n', 'rflux1gc1', 'rflux2gc1',
                                'rflux3gc1', 'rflux1g2', 'rflux2gc2',
                                'rflux3gc2'))

master_astrometry = ascii.read(path_radec)
master_dither = ascii.read(path_dithers, names=['x1', 'y1', 'x2', 'y2',
                                                'x3', 'y3', 'date', 'obsid'])

shotarray = np.unique(master_sci['SHOT'])

for idx in np.arange(np.size(shotarray)):

    row = tableMain.row

    row['shotid'] = np.int(shotarray[idx])
    row['date'] = str(row['shotid'])[0:8]
    row['obsid'] = str(row['shotid'])[8:11]
    row['datevshot'] = str(row['shotid'])[0:8] + 'v' + str(row['shotid'])[8:11]
    
    sel = np.where(master_sci['SHOT'] == shotarray[idx])
    row['objid'] = master_sci['OBJECT'][sel][0]
    row['field'] = define_field(row['objid'])

    row['trajcra'] = master_sci['TRAJCRA'][idx]
    row['trajcdec'] = master_sci['TRAJCDEC'][idx]
    row['trajcpa'] = master_sci['PARANGLE'][idx]

    sel = np.where((master_astrometry['col1'] == np.int(row['date']))
                   * (master_astrometry['col2'] == np.int(row['obsid'])))

    if np.size(sel) == 1:
        row['ra'] = master_astrometry['col3'][sel]
        row['dec'] = master_astrometry['col4'][sel]
        row['pa'] = master_astrometry['col5'][sel]
    elif np.size(sel) > 1:
        print "Weird: More than one RA/DEC/PA value for", row['shotid']
        print "Saving the first RA/DEC/PA value only. Make sure this is ok."
        row['ra'] = master_astrometry['ra'][sel][0]
        row['dec'] = master_astrometry['dec'][sel][0]
        row['pa'] = master_astrometry['pa'][sel][0]
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

    # add in dither info

    sel = np.where((master_dither['date'] == row['date'])
                   * (master_dither['obsid'] == row['obsid']))
    row['ditherpos'] = master_dither[sel]
    dith_arr = ['x1', 'y1', 'x2', 'y2', 'x3', 'y3']
    if np.size(sel) == 1:
        for dith in dith_arr:
            row[dith] = master_dither[dith][sel]
    elif np.size(sel) == 0:
        for dith in dith_arr:
            row[dith] = np.nan

    row.append()

tableMain.flush()


for idx in np.arange(np.size(master_sci['SHOT'])):

    row = tableExp.row

    row['shotid'] = np.int(master_sci['SHOT'][idx])
    row['date'] = int(str(row['shotid'])[0:8])
    row['obsid'] = int(str(row['shotid'])[8:11])
    row['expnum'] = int(str(master_sci['DITHER'][idx])[3:5])
    row['timestamp'] = master_sci['TIMESTAMP'][idx]
    row['mjd'] = master_sci['MJD'][idx]
    row['exptime'] = master_sci['EXPTIME'][idx]
    row['darktime'] = master_sci['DARKTIME'][idx]

    row.append()

tableExp.flush()

fileh.close()
