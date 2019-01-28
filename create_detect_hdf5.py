# -*- coding: utf-8 -*-
"""
Created: 2019/01/25

@author: Erin Mentuch Cooper

This file takes a list of shots and creates a HDF5 of all line detections

"""

import sys
import os
import os.path as op
import re
import subprocess
import numpy as np
import tables as tb

from astropy.io import ascii
from astropy.io import fits

#
# eventually get this from /work/00115/gebhardt/maverick/api/paths.cfg
#
path_detect = '/work/05178/cxliu/maverick/detect/'
path_elixer = '/work/05350/ecooper/maverick/elixer/'


class Detections(tb.IsDescription):
    shotid = tb.Int64Col()
    date = tb.Time32Col()
    obsid = tb.Int32Col()
    detectid = tb.Int64Col()
    ra = tb.Float32Col()
    dec = tb.Float32Col()
    wave = tb.Float32Col()
    wave_err = tb.Float32Col()
    flux = tb.Float32Col()
    flux_err = tb.Float32Col()
    linewidth = tb.Float32Col()
    linewidth_err = tb.Float32Col()
    continuum = tb.Float32Col()
    continuum_err = tb.Float32Col()
    sn = tb.Float32Col()
    sn_err = tb.Float32Col()
    chi2 = tb.Float32Col()
    chi2_err = tb.Float32Col()
    x_raw = tb.Int32Col()
    y_raw = tb.Int32Col()
    fiber_num = tb.Int32Col()
    multiframe = tb.StringCol(20)


fileh = tb.open_file("detect_test.h5", mode="w", title="Detections test file ")
group = fileh.create_group(fileh.root, 'Info', 'HETDEX Detect Catalog')
tableMain = fileh.create_table(group, 'Detections', Detections,
                               'Line Detection Catalog')
detectfile = '/work/00115/gebhardt/maverick/detect/dexall/xypos/res2018all'

os.system("grep -v '\*\*\*' /work/00115/gebhardt/maverick/detect/dexall/xypos/res2018all > tmp")

detectcat = ascii.read('tmp')


for idx in np.arange(np.size(detectcat)):

    row = tableMain.row

    datevobs = detectcat['col16'][idx]
    p = re.compile('v')
    row['shotid'] = p.sub('', datevobs)
    row['date'] = str(row['shotid'])[0:8]
    row['obsid'] = str(row['shotid'])[8:11]
    row['detectid'] = detectcat['col15'][idx]
    row['ra'] = detectcat['col1'][idx]
    row['dec'] = detectcat['col2'][idx]
    row['wave'] = detectcat['col3'][idx]
    row['wave_err'] = detectcat['col4'][idx]
    row['flux'] = detectcat['col5'][idx]
    row['flux_err'] = detectcat['col6'][idx]
    row['linewidth'] = detectcat['col7'][idx]
    row['linewidth_err'] = detectcat['col8'][idx]
    row['continuum'] = detectcat['col9'][idx]
    row['continuum_err'] = detectcat['col10'][idx]
    row['sn'] = detectcat['col11'][idx]
    row['sn_err'] = detectcat['col12'][idx]
    row['chi2'] = detectcat['col13'][idx]
    row['chi2_err'] = detectcat['col14'][idx]
    row['x_raw'] = detectcat['col17'][idx]
    row['y_raw'] = detectcat['col18'][idx]
    row['fiber_num'] = detectcat['col19'][idx]
    filemulti = detectcat['col20'][idx]
    idx2 = filemulti.find('multi')
    row['multiframe'] = filemulti[idx2:idx2+20]
    row.append()

tableMain.flush()

fileh.close()
