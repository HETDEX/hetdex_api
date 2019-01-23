# -*- coding: utf-8 -*-
"""
Created: 2019/01/23

@author: Erin Mentuch Cooper

Note: requires makemaster.sh to collect info from Karl's
/work/00115/gebhardt/maverick/gettar
database and from other lookup tables

"""

import sys
import os
import os.path as op
import subprocess
import numpy as np
import tables as tb

from astropy.io import ascii
from astropy.io import fits

#
# eventually get this from /work/00115/gebhardt/maverick/api/paths.cfg
#

path_astrometry = "/work/00115/gebhardt/maverick/vdrp/shifts/"
path_gpinfo = '/work/05350/ecooper/maverick/gettar/master_fwhm_DR1'
mastersci_file = "/work/05350/ecooper/maverick/gettar/mastersci_DR1"
mastertp_file = '/work/05350/ecooper/maverick/gettar/master_tp'
path_radec = "/work/00115/gebhardt/maverick/getfib/radec.all"
path_fplane = "/work/00115/gebhardt/maverick/vdrp/shifts/"
path_throughput = "/work/00115/gebhardt/maverick/detect/tp/"
path_amp2amp = "/work/00115/gebhardt/maverick/getampnorm/all/"


class Survey(tb.IsDescription):
    shotid = tb.Int64Col()
    date = tb.Time32Col()
    obsid = tb.Int32Col()
    objid = tb.StringCol(18)
    ra = tb.Float32Col()
    dec = tb.Float32Col()
    pa = tb.Float32Col()
    trajcra = tb.Float32Col()
    trajcdec = tb.Float32Col()
    trajcpa = tb.Float32Col()
    fwhm = tb.Float32Col()
    response_4540 = tb.Float32Col()  # normalized for 360s
    twi_shotid = tb.Int64Col()


class Exp(tb.IsDescription):
    shotid = tb.Int64Col()
    date = tb.Time32Col()
    obsid = tb.Int32Col()
    expn = tb.StringCol(5)
    timestamp = tb.StringCol(18)
    mjd = tb.Float32Col()
    exptime = tb.Float32Col(7)
    darktime = tb.Float32Col()
    cofits = tb.Float32Col((1300, 1300))


fileh = tb.open_file("survey_test_new.h5", mode="w", title="Survey test file ")

group = fileh.create_group(fileh.root, 'Info', 'HETDEX Survey Info')
tableMain = fileh.create_table(group, 'Survey', Survey, 'Main Survey Info')
tableExp = fileh.create_table(group, 'Exp', Exp, 'Exposure Specific Info')

master_sci = ascii.read(mastersci_file)
master_fwhm = ascii.read(path_gpinfo)
master_astrometry = ascii.read(path_radec)
master_tp = ascii.read(mastertp_file)

shotarray = np.unique(master_sci['SHOT'])

for idx in np.arange(np.size(shotarray)):

    row = tableMain.row

    row['shotid'] = np.int(shotarray[idx])
    row['date'] = str(row['shotid'])[0:8]
    row['obsid'] = str(row['shotid'])[8:11]

    sel = np.where(master_sci['SHOT'] == shotarray[idx])
    row['objid'] = master_sci['OBJECT'][sel][0]

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

    sel1 = np.where(master_fwhm['SHOT'] == row['shotid'])

    if np.size(sel1) == 1:
        row['fwhm'] = np.float(master_fwhm['FWHM'][sel1])
    elif np.size(sel1) > 1:
        print "Weird: More than one FWHM value for", row['shotid']
        print "Saving the first FWHM value only. Make sure this is ok."
        row['fwhm'] = np.float(master_fwhm['FWHM'][sel1][0])
    else:
        row['fwhm'] = np.nan

    if row['fwhm'] < 0:
        row['fwhm'] = np.nan

    sel2 = np.where((master_tp['DATE'] == row['date'])
                    * (master_tp['OBSID'] == row['obsid']))

    if np.size(sel2) == 1:
        row['response_4540'] = master_tp['TP'][sel2]
    elif np.size(sel2) > 1:
        print "Weird: More than two TPs found for", row['shotid']
        row['response_4540'] = master_tp['TP'][sel2][0]
    else:
        row['response_4540'] = np.nan

    row.append()

tableMain.flush()


for idx in np.arange(np.size(master_sci['SHOT'])):

    row = tableExp.row

    row['shotid'] = np.int(master_sci['SHOT'][idx])
    row['date'] = str(row['shotid'])[0:8]
    row['obsid'] = str(row['shotid'])[8:11]
    row['expn'] = master_sci['DITHER'][idx]
    row['timestamp'] = master_sci['TIMESTAMP'][idx]
    row['mjd'] = master_sci['MJD'][idx]
    row['exptime'] = master_sci['EXPTIME'][idx]
    row['darktime'] = master_sci['DARKTIME'][idx]

    fitsfile = op.join(path_astrometry,
                       str(row['date']) + 'v' + str(row['obsid']).zfill(3),
                       str(row['date']) + 'v' + str(row['obsid']).zfill(3)
                       + 'fp_' + row['expn'] + '.fits')

    if os.path.exists(fitsfile):
        F = fits.open(fitsfile)
        row['cofits'] = F[0].data * 1.

    row.append()

tableExp.flush()

fileh.close()
