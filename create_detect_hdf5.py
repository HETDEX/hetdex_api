# -*- coding: utf-8 -*-
"""
Created: 2019/01/25

@author: Erin Mentuch Cooper

This file contains all information related to the HETDEX line detections
catalog

"""

import sys
import os
import os.path as op
import re
import glob
import subprocess
import numpy as np
import tables as tb

from astropy.io import ascii
from astropy.io import fits

path_detect = '/tmp/HETDEX'


def build_spec_path(path_detects, date, obsid, detectID):
    specf_str = str(date) + 'v' + str(obsid).zfill(3) + '_' + str(detectID) + 'specf.dat'
    return op.join(path_detects, specf_str)


def build_mcres_path(path_detects, date, obsid, detectID):
    mcres_str = str(date) + 'v' + str(obsid).zfill(3) + '_' + str(detectID) + 'mc.res'
    return op.join(path_detects, mcres_str)


def build_fiberinfo_path(path_detects, date, obsid, detectID):
    info_str = str(date) + 'v' + str(obsid).zfill(3) + '_' + str(detectID) + '.info'
    return op.join(path_detects, info_str)


class Detections(tb.IsDescription):
    shotid = tb.Int64Col(pos=2)
    date = tb.Int32Col(pos=3)
    obsid = tb.Int32Col(pos=4)
    detectid = tb.Int64Col(pos=0) 
    detectname = tb.StringCol((20), pos=1)  
    ra = tb.Float32Col(pos=5)
    dec = tb.Float32Col(pos=6)
    wave = tb.Float32Col(pos=7)
    wave_err = tb.Float32Col(pos=8)
    flux = tb.Float32Col(pos=9)
    flux_err = tb.Float32Col(pos=10)
    linewidth = tb.Float32Col(pos=11)
    linewidth_err = tb.Float32Col(pos=12)
    continuum = tb.Float32Col(pos=13)
    continuum_err = tb.Float32Col(pos=14)
    sn = tb.Float32Col(pos=15)
    sn_err = tb.Float32Col(pos=16)
    chi2 = tb.Float32Col(pos=17)
    chi2_err = tb.Float32Col(pos=18)
    x_raw = tb.Int32Col(pos=21)
    y_raw = tb.Int32Col(pos=22)
    fiber_num = tb.Int32Col(pos=20)
    multiframe = tb.StringCol((20), pos=19)
    specid = tb.StrinCol((3))
    ifuslot = tb.StringCol((3))
    ifuid = tb.StringCol((3))
    amp = tb.StringCol((2))
    inputid = tb.StringCol((32)) 


class Spectra(tb.IsDescription):
    detectid = tb.Int64Col(pos=0)
    wave1d = tb.Float32Col(1036, pos=1)
    spec1d = tb.Float32Col(1036, pos=2)
    spec1d_err = tb.Float32Col(1036, pos=3)
    counts1d = tb.Float32Col(1036, pos=4)
    counts1d_err = tb.Float32Col(1036, pos=5)
    apsum_counts = tb.Float32Col(1036, pos=6)
    apsum_counts_err = tb.Float32Col(1036, pos=7)


class Fibers(tb.IsDescription):
    detectid = tb.Int64Col(pos=0)
    ra = tb.Float32Col(pos=1)
    dec = tb.Float32Col(pos=2)
    x_ifu = tb.Float32Col(pos=5)
    y_ifu = tb.Float32Col(pos=6)
    multiframe = tb.StringCol((20), pos=3)
    fiber_num = tb.Int32Col(pos=4)
    expn = tb.StringCol((5), pos=9)
    distance = tb.Float32Col(pos=10)
    wavein = tb.Float32Col(pos=12)
    timestamp = tb.StringCol((17), pos=11)
    date = tb.Int32Col(pos=7)
    obsid = tb.Int32Col(pos=8)
    flag = tb.Int32Col(pos=13)
    weight = tb.Float32Col(pos=14)
    ADC = tb.Float32Col((5), pos=15)
    specid = tb.StrinCol((3))
    ifuslot = tb.StringCol((3))
    ifuid = tb.StringCol((3))
    amp = tb.StringCol((2))


fileh = tb.open_file("detect_big.h5", mode="w", title="Detections test file ")

tableMain = fileh.create_table(fileh.root, 'Detections', Detections,
                               'HETDEX Line Detection Catalog')

tableFibers = fileh.create_table(fileh.root, 'Fibers', Fibers,
                                 'Fiber info for each detection')

tableSpectra = fileh.create_table(fileh.root, 'Spectra', Spectra,
                                  '1D Spectra for each Line Detection')

detectidx = 1000000000

#resfile = ascii.read('/work/00115/gebhardt/maverick/detect/dexall/select/res_use')
# for this we are using a list in DATEvOBS_DET form as we would get this from rsp3mc 
# calls and can be grabbed from a directory directly 

f = open('detsbig.list',"r")

for line in f:
    datevobs_det = line.rstrip()
    datevobs = datevobs_det[0:12]
    date = datevobs_det[0:8]
    obs = datevobs_det[10:12]
    det = datevobs_det[13:]    

    print("Ingesting input ID:" + str(datevobs_det) + "   Date="
          + date + "  OBSID="+ obs + "  ID=" + det+"\n")

    detectfile = build_mcres_path(path_detect, date, obs, det) 
    
    if op.exists(detectfile):
        
        row = tableMain.row
        
        detectcat = ascii.read(detectfile, delimiter=' ')
        row['detectid'] = detectidx
        p = re.compile('v')
        row['shotid'] = p.sub('', datevobs)
        row['date'] = date
        row['obsid'] = obs
        row['inputid'] = datevobs_det
        row['ra'] = detectcat['col1']
        row['dec'] = detectcat['col2']
        row['wave'] = detectcat['col3']
        row['wave_err'] = detectcat['col4']
        row['flux'] = detectcat['col5']
        row['flux_err'] = detectcat['col6']
        row['linewidth'] = detectcat['col7']
        row['linewidth_err'] = detectcat['col8']
        row['continuum'] = detectcat['col9']
        row['continuum_err'] = detectcat['col10']
        row['sn'] = detectcat['col11']
        row['sn_err'] = detectcat['col12']
        row['chi2'] = detectcat['col13']
        row['chi2_err'] = detectcat['col14']
        row['x_raw'] = detectcat['col17']
        row['y_raw'] = detectcat['col18']
        row['fiber_num'] = detectcat['col19']
        filemulti = detectcat['col20'][0]
        idx = filemulti.find('multi')
        multiframe = filemulti[idx:idx+20]
        row['multiframe'] = multiframe
        row['specid'] = multiframe[6:9]
        row['ifuslot'] = multiframe[10:13]
        row['ifuid'] = multiframe[14:17]
        row['amp'] = multiframe[18:20]


        # now populate table with 1D spectra, queried by detectid
               
        filespec = build_spec_path(path_detect, date, obs, det)
        if op.exists(filespec):
            rowspectra = tableSpectra.row
            dataspec = ascii.read(filespec)
            rowspectra['detectid'] = detectidx
            rowspectra['wave1d'] = dataspec['col1']
            rowspectra['spec1d'] = dataspec['col2']
            rowspectra['spec1d_err'] = dataspec['col3']
            rowspectra['counts1d'] = dataspec['col4']
            rowspectra['counts1d_err'] = dataspec['col5']
#            rowspectra['apsum_counts'] = dataspec['col6']
#            rowspectra['apsum_counts_err'] = dataspec['col7']
            rowspectra.append()
            
        # now populate fiber table with additional info
        filefiberinfo = build_fiberinfo_path(path_detect, date, obs, det) 
        
        if op.exists(filefiberinfo):
            datafiber = ascii.read(filefiberinfo, format='no_header', delimiter=' ', data_end = 3)
                
            for ifiber in np.arange(np.size(datafiber)):
                rowfiber = tableFibers.row
                
                rowfiber['detectid'] = detectidx
                rowfiber['ra'] = datafiber['col1'][ifiber]
                rowfiber['dec'] = datafiber['col2'][ifiber]
                rowfiber['x_ifu'] = datafiber['col3'][ifiber]
                rowfiber['y_ifu'] = datafiber['col4'][ifiber]
                multiname = datafiber['col5'][ifiber]
                multiframe = multiname[0:20]
                rowfiber['multiframe'] = multiframe
                rowfiber['multiframe'] = multiframe
                rowfiber['specid'] = multiframe[6:9]
                rowfiber['ifuslot'] = multiframe[10:13]
                rowfiber['ifuid'] = multiframe[14:17]
                rowfiber['amp'] = multiframe[18:20]
                rowfiber['fiber_num'] = multiname[21:22]
                rowfiber['expn'] = datafiber['col6'][ifiber]
                rowfiber['distance'] = datafiber['col7'][ifiber]
                rowfiber['wavein'] = datafiber['col8'][ifiber]
                rowfiber['timestamp'] = datafiber['col9'][ifiber]
                rowfiber['date'] = datafiber['col10'][ifiber]
                rowfiber['obsid'] = str(datafiber['col11'][ifiber])[0:3]
                rowfiber['flag'] = datafiber['col12'][ifiber]
                rowfiber['weight'] = datafiber['col13'][ifiber]
                rowfiber['ADC'] = [datafiber['col14'][ifiber],
                                   datafiber['col15'][ifiber],
                                   datafiber['col16'][ifiber],
                                   datafiber['col17'][ifiber], 
                                   datafiber['col18'][ifiber]]
        
                rowfiber.append()
        detectidx += 1
        row.append()

tableMain.flush()
tableFibers.flush()
tableSpectra.flush()

#create completely sorted index on the detectid to make queries against that column much faster
tableFibers.cols.detectid.create_csindex()
tableSpectra.cols.detectid.create_csindex()
tableFibers.flush() #just to be safe
tableSpectra.flush()

fileh.close()
