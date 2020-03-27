# -*- coding: utf-8 -*-
"""
Created: 2019/01/25

@author: Erin Mentuch Cooper

This file contains all information related to the HETDEX line detections
catalog

Examples
--------

To create for a month:

>>> python3 create_detect_hdf5.py -m 201901 -of detect_201901.h5

Then once all months are done, merge into one file:

>>> python3 create_detect_hdf5.py --merge -of detect_hdr2.h5

To run continuum sources:

>>> python3 create_detect_hdf5.py -dp /data/00115/gebhardt/cs/spec /
-cs /data/00115/gebhardt/cs/rext1 -of continuum_sources.h5

"""
from __future__ import print_function

import sys
import os
import os.path as op
import argparse as ap
import re
import glob
import subprocess
import numpy as np
import tables as tb

from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
from hetdex_api.input_utils import setup_logging


def build_spec_path(path_detects, date, obsid, detectID):
    specf_str = str(date) + 'v' + str(obsid).zfill(3) + '_' + str(detectID) + 'specf.dat'
    return op.join(path_detects, specf_str)


def build_mcres_path(path_detects, date, obsid, detectID):
    mcres_str = str(date) + 'v' + str(obsid).zfill(3) + '_' + str(detectID) + 'mc.res'
    return op.join(path_detects, mcres_str)


def build_fiberinfo_path(path_detects, date, obsid, detectID):
    info_str = str(date) + 'v' + str(obsid).zfill(3) + '_' + str(detectID) + '.info'
    return op.join(path_detects, info_str)

def emission_line_to_IAU_string(coord, lambda_in_A):
        """
        convert ra,dec coordinates in the form of a astropy.skycoord object 
        and a float wavelength in Angstroms to a IAU-style object name.
        """
        return 'HETDEX J{0}{1}w{2}'.format(coord.ra.to_string(unit=u.hourangle, sep='', precision=2, pad=True), coord.dec.to_string(sep='', precision=2, alwayssign=True, pad=True), '%+4.2f'%lambda_in_A)

class Detections(tb.IsDescription):
    shotid = tb.Int64Col(pos=2)
    date = tb.Int32Col(pos=3)
    obsid = tb.Int32Col(pos=4)
    detectid = tb.Int64Col(pos=0) 
    fiber_id = tb.StringCol((38), pos=1)
    #detectname = tb.StringCol((40))  
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
    inputid = tb.StringCol((40), pos=3)
#    x_raw = tb.Int32Col(pos=21)
#    y_raw = tb.Int32Col(pos=22)
#    x_ifu = tb.Float32Col()
#    y_ifu = tb.Float32Col()
#    weight = tb.Float32Col()
#    fibnum = tb.Int32Col(pos=20)
#    multiframe = tb.StringCol((20), pos=19)
#    specid = tb.StringCol((3))
#    ifuslot = tb.StringCol((3))
#    ifuid = tb.StringCol((3))
#    amp = tb.StringCol((2))
#    expnum = tb.Int32Col()

class Spectra(tb.IsDescription):
    detectid = tb.Int64Col(pos=0)
    wave1d = tb.Float32Col(1036, pos=1)
    spec1d = tb.Float32Col(1036, pos=2)
    spec1d_err = tb.Float32Col(1036, pos=3)
    counts1d = tb.Float32Col(1036, pos=4)
    counts1d_err = tb.Float32Col(1036, pos=5)
    apsum_counts = tb.Float32Col(1036, pos=6)
    apsum_counts_err = tb.Float32Col(1036, pos=7)
    spec1d_nc = tb.Float32Col(1036)
    spec1d_nc_err = tb.Float32Col(1036)
    apcor = tb.Float32Col(1036)

class Fibers(tb.IsDescription):
    detectid = tb.Int64Col(pos=0)
    ra = tb.Float32Col(pos=1)
    dec = tb.Float32Col(pos=2)
    x_ifu = tb.Float32Col(pos=5)
    y_ifu = tb.Float32Col(pos=6)
    multiframe = tb.StringCol((20), pos=3)
    fibnum = tb.Int32Col()
    fiber_id = tb.StringCol((38), pos=4)        
    expnum = tb.Int32Col(pos=9)
    distance = tb.Float32Col(pos=10)
    wavein = tb.Float32Col(pos=12)
    timestamp = tb.StringCol((17), pos=11)
    date = tb.Int32Col(pos=7)
    obsid = tb.Int32Col(pos=8)
    flag = tb.Int32Col(pos=13)
    weight = tb.Float32Col(pos=14)
    ADC = tb.Float32Col((5), pos=15)
    specid = tb.StringCol((3))
    ifuslot = tb.StringCol((3))
    ifuid = tb.StringCol((3))
    amp = tb.StringCol((2))
    x_raw = tb.Int32Col()
    y_raw = tb.Int32Col()


def get_detect_cat(detectidx, catfile, args):
    '''
    This is for Chenxu's latest .cat files
    '''

    detectid = []
    multiframe = []
    fiber_id = []
    ifuslot = []
    specid = []
    ifuid = []
    fibnum = []
    shotid = []
    amp = []
    date = []
    obsid = []
    expnum = []
    
    detectid_i = detectidx

    detectcat = Table.read(catfile, format='ascii')
    
    for row in detectcat:

        detectid.append(detectid_i)
        detectid_i += 1

        p = re.compile('v')
        shotid_i = int(p.sub('', row['datevshot']))

        date.append( int(str(shotid_i)[0:8]))
        obsid.append( int(str(shotid_i)[8:11]))
        shotid.append( shotid_i)

        if True: #args.cxcat:
            continue
        else:
        
            fiber_name = row['fiber_name']
        
            multiframe_i = fiber_name[0:20]
            multiframe.append(multiframe_i)
        
            fibnum_i = int(fiber_name[21:24])
            fibnum.append(fibnum_i)
        
            fiber_id_i = str(shotid_i) + '_' + str(int(row['expnum'][-2:])) + '_' + multiframe_i + '_' + str(fibnum_i).zfill(3)

            fiber_id.append( fiber_id_i)
            specid.append( multiframe_i[6:9])
            ifuslot.append( multiframe_i[10:13])
            ifuid.append( multiframe_i[14:17])
            amp.append( multiframe_i[18:20])
            expnum.append(int(row['expnum'][-2:]))

    detectcat['date'] = date
    detectcat['obsid'] = obsid
    detectcat['detectid'] = detectid
    detectcat['shotid'] = shotid
    detectcat['inputid'] = detectcat['hdr2_id'].astype(bytes)
    
    if True:
        return detectcat 
    else:
        detectcat['expnum'] = np.array(expnum).astype(int)
        detectcat['multiframe'] = np.array(multiframe).astype(bytes)
        detectcat['fiber_id'] = np.array(fiber_id).astype(bytes)
        detectcat['ifuslot'] = np.array(ifuslot).astype(bytes)
        detectcat['specid'] = np.array(specid).astype(bytes)
        detectcat['ifuid'] = np.array(ifuid).astype(bytes)
        detectcat['fibnum'] = fibnum
        detectcat['amp'] = np.array(amp).astype(bytes)
        
        detectcat['X_amp'].name = 'x_raw'
        detectcat['Y_amp'].name = 'y_raw'
        detectcat['X_FP'].name = 'x_ifu'
        detectcat['Y_FP'].name = 'y_ifu'
        
        detectcat.remove_columns(['datevshot', 'original_name', 'fiber_name', 'hdr2_id',
                                  'linewidth_fix', 'linewidth_fix_err','chi2_fix'])
    
        return detectcat
                                                                        
                                                                        
def append_detection(detectidx, date, obs, det, detect_path, tableMain,
                     tableFibers, tableSpectra):
    '''
    This was for HDR1
    '''
    
    print("Ingesting Date=" + date + "  OBSID="+ obs + "  ID=" + det)
    detectfile = build_mcres_path(detect_path, date, obs, det)
    datevobs_det = str(date) + 'v' + str(obs).zfill(3) + '_' + str(det)
    datevobs = str(date) + 'v' + str(obs).zfill(3)
    
    try:
        detectcat = ascii.read(detectfile, delimiter=' ')
        row = tableMain.row
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
        coord = SkyCoord(row['ra'] * u.degree,row['dec'] * u.degree, frame='icrs')
        row['detectname'] = emission_line_to_IAU_string(coord, row['wave'])
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
        row['fibnum'] = detectcat['col19']
        filemulti = detectcat['col20'][0]
        idx = filemulti.find('multi')
        multiframe = filemulti[idx:idx+20]
        row['multiframe'] = multiframe
        row['specid'] = multiframe[6:9]
        row['ifuslot'] = multiframe[10:13]
        row['ifuid'] = multiframe[14:17]
        row['amp'] = multiframe[18:20]
        row.append()
        
        # now populate table with 1D spectra, queried by detectid
        filespec = build_spec_path(detect_path, date, obs, det)
        try:
            dataspec = ascii.read(filespec)
            rowspectra = tableSpectra.row
            dataspec = ascii.read(filespec)
            rowspectra['detectid'] = detectidx
            rowspectra['wave1d'] = dataspec['col1']
            rowspectra['spec1d'] = dataspec['col2']
            rowspectra['spec1d_err'] = dataspec['col3']
            rowspectra['counts1d'] = dataspec['col4']
            rowspectra['counts1d_err'] = dataspec['col5']
            rowspectra['apsum_counts'] = dataspec['col6']
            rowspectra['apsum_counts_err'] = dataspec['col7']
            rowspectra.append()
        except:
            print('Could not ingest %s' % filespec)

        # now populate fiber table with additional info
        filefiberinfo = build_fiberinfo_path(detect_path, date, obs, det)
        try:
            datafiber = ascii.read(filefiberinfo, format='no_header', delimiter=' ')
            
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
                rowfiber['specid'] = multiframe[6:9]
                rowfiber['ifuslot'] = multiframe[10:13]
                rowfiber['ifuid'] = multiframe[14:17]
                rowfiber['amp'] = multiframe[18:20]
                rowfiber['fibnum'] = int(multiname[21:24])
                rowfiber['expnum'] = str(datafiber['col6'][ifiber])[3:5]
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
        except:
            print('Could not ingest %s' % filefiberinfo)
            
        detectidx += 1
    except:
            print('Could not ingest %s' % detectfile)
    
    return detectidx

def main(argv=None):
    ''' Main Function '''
    # Call initial parser from init_utils
    parser = ap.ArgumentParser(description="""Create HDF5 file.""",
                               add_help=True)

    parser.add_argument("-m", "--month",
                        help='''Month to run: 201901''',
                        type=str, default=None)
                        
    parser.add_argument("-d", "--date",
                        help='''Date, e.g., 20170321, YYYYMMDD''',
                        type=str, default=None)

    parser.add_argument("-o", "--observation",
                        help='''Observation number, "00000007" or "7"''',
                        type=str, default=None)

    parser.add_argument("-id", "--inputid", 
                        help='''Detection input ID, "string" or "532"''',
                        type=str, default=None)
    
    parser.add_argument("-dets", "--dets",
                        help='''List of detections in form DATEvSHOT_inputID''',
                        type=str, default=None)

    parser.add_argument('-of', '--outfilename', type=str,
                        help='''Relative or absolute path for output HDF5
                        file.''', default=None)
    
    parser.add_argument('-a', '--append',
                        help='''Appending to existing detection HDF5 file.''',
                        action="count", default=0)
    
    parser.add_argument("-dp", "--detect_path",
                        help='''Path to detections''',
                        type=str, default='/data/05178/cxliu/detect/')

    parser.add_argument("-cs","--contsource",
                        help='''Path to Karls rext1 catalog''',
                        type=str, default=None)

    parser.add_argument("-cxcat", "--cxcat",
                        help='''Path to chenxus catalog''',
                        type=str, default=None)
                         
    parser.add_argument("-md", "--mergedir",
                        help='''Merge all HDF5 files in the defined merge 
                        directory. Can append to existing file using --append option''',
                        type=str, default=os.getcwd())

    parser.add_argument("--merge",
                        "-merge",
                        help="""Boolean trigger to merge all 2*.fits files in cwd""",
                        default=False,
                        required=False,
                        action="store_true")

    args = parser.parse_args(argv)
    args.log = setup_logging()

    # Creates a new file if the "--append" option is not set or the file
    # does not already exist.
    does_exist = False
    if op.exists(args.outfilename) and args.append:
        fileh = tb.open_file(args.outfilename, 'a', 'HDR2 Detections Database')
        does_exist = True
        # initiate new unique detect index
        detectidx = np.max(fileh.root.Detections.cols.detectid) + 1
    else:
        fileh = tb.open_file(args.outfilename, 'w', 'HDR2 Detections Database')
        index_buff = 2000000000
        detectidx = index_buff
        

    if args.merge:
        
        tableMain = fileh.create_table(fileh.root, 'Detections', Detections,
                                       'HETDEX Line Detection Catalog',
                                       expectedrows=1000000)
        tableFibers = fileh.create_table(fileh.root, 'Fibers', Fibers,
                                         'Fiber info for each detection',
                                         expectedrows=1000000)
        tableSpectra = fileh.create_table(fileh.root, 'Spectra', Spectra,
                                          '1D Spectra for each Line Detection',
                                          expectedrows=15000000)
        
        files = sorted(glob.glob(op.join(args.mergedir,'detect_2*.h5')))
        
        detectid_max = 1

        for file in files:

            print('Appending detect H5 file: %s' % file)
            
            fileh_i = tb.open_file(file, 'r')
            
            tableMain_i = fileh_i.root.Detections.read()                
            tableFibers_i = fileh_i.root.Fibers.read()
            tableSpectra_i = fileh_i.root.Spectra.read()

            tableMain_i['detectid'] += detectid_max
            tableFibers_i['detectid'] += detectid_max
            tableSpectra_i['detectid'] += detectid_max

            tableMain.append(tableMain_i)
            tableFibers.append(tableFibers_i)
            tableSpectra.append(tableSpectra_i)
            
            detectid_max = np.max(tableMain.cols.detectid[:]) - index_buff + 1
            
            fileh_i.close()
            tableFibers.flush() #just to be safe
            tableSpectra.flush()
            tableMain.flush()

    else:

        if args.month:
            catfile = op.join(args.detect_path, args.month, args.month + '.cat')

            detectcat = get_detect_cat(detectidx, catfile, args)

        elif args.contsource:
            detectcat = Table.read(args.contsource, format='ascii')
            detectcat.remove_columns(['col1','col4','col5','col6','col9','col10',
                                      'col11','col12','col13','col14'])
            detectcat['col2'].name = 'ra'
            detectcat['col3'].name = 'dec'
            detectcat['col7'].name = 'obnum'
            detectcat['col8'].name = 'datevshot'

            shotid = []
            date = []
            obsid = []
            inputid = []
            detectid = []
            
            detectid_i = detectidx

            for row in detectcat:
                p = re.compile('v')
                shotid_i = int(p.sub('', row['datevshot']))
                inputid_i = str(row['datevshot']) + '_' + str(row['obnum'])
            
                detectid.append(detectid_i)
                inputid.append(inputid_i)
                date.append( int(str(shotid_i)[0:8]))
                obsid.append( int(str(shotid_i)[8:11]))
                shotid.append( shotid_i)
                detectid_i += 1

            detectcat['detectid'] = detectid
            detectcat['inputid'] = inputid
            detectcat['date'] = date
            detectcat['obsid'] = obsid
            detectcat['detectid'] = detectid
            detectcat['shotid'] = shotid

        elif args.cxcat:
            detectcat = get_detect_cat(detectidx, args.cxcat, args)
            
        tableMain = fileh.create_table(fileh.root, 'Detections', Detections,
                                       'HETDEX Line Detection Catalog',
                                       expectedrows=np.size(detectcat))
        tableFibers = fileh.create_table(fileh.root, 'Fibers', Fibers,
                                         'Fiber info for each detection',
                                         expectedrows=np.size(detectcat))
        tableSpectra = fileh.create_table(fileh.root, 'Spectra', Spectra,
                                          '1D Spectra for each Line Detection',
                                          expectedrows=15*np.size(detectcat))
        
        det_cols = fileh.root.Detections.colnames

        for row in detectcat:
            rowMain = tableMain.row
            
            for col in det_cols:
                try:
                    rowMain[col] = row[col]
                except:
                    rowMain[col] = 0.0
            rowMain.append()

        
        # add spectra for each detectid in the detections table
        tableMain.flush()

        for row in tableMain:
            try:
                inputid_i = row['inputid'].decode()

                if args.month:
                    specfile = op.join(args.detect_path, args.month, 'rf', inputid_i + '.spec')
                else:
                    specfile = op.join(args.detect_path, inputid_i + '.spec')

                dataspec = Table.read(specfile, format='ascii.no_header')
                rowspectra = tableSpectra.row
                rowspectra['detectid'] = row['detectid']
                rowspectra['spec1d'] = dataspec['col2']/dataspec['col9']
                rowspectra['spec1d_err'] = dataspec['col3']/dataspec['col9']
                rowspectra['wave1d'] = dataspec['col1']
                rowspectra['spec1d_nc'] = dataspec['col2']
                rowspectra['spec1d_nc_err'] = dataspec['col3']
                rowspectra['counts1d'] = dataspec['col4']
                rowspectra['counts1d_err'] = dataspec['col5']
                rowspectra['apsum_counts'] = dataspec['col6']
                rowspectra['apsum_counts_err'] = dataspec['col7']
                rowspectra['apcor'] = dataspec['col9']
                rowspectra.append()
            except:
                args.log.error('Could not ingest %s' % specfile)
            
        tableSpectra.flush()
        
        # add fiber info for each detection
    
        for row in tableMain:
            inputid_i = row['inputid'].decode()

            if args.month:
                filefiberinfo = op.join(args.detect_path, args.month, 'rf', inputid_i + '.list')
            else:
                filefiberinfo = op.join(args.detect_path, inputid_i + '.list')
                
            try:
                datafiber = Table.read(filefiberinfo, format='ascii.no_header')
            
                for ifiber in np.arange(np.size(datafiber)):
                    rowfiber = tableFibers.row
                    rowfiber['detectid'] = row['detectid']
                    rowfiber['ra'] = datafiber['col1'][ifiber]
                    rowfiber['dec'] = datafiber['col2'][ifiber]
                    rowfiber['x_ifu'] = datafiber['col3'][ifiber]
                    rowfiber['y_ifu'] = datafiber['col4'][ifiber]
                    rowfiber['expnum'] = str(datafiber['col6'][ifiber])[3:5]
                    multiname = datafiber['col5'][ifiber]
                    multiframe = multiname[0:20]
                    fiber_id_i = str(row['shotid']) + '_' + str(int(rowfiber['expnum'])) + '_' + multiframe + '_' + str(int(multiname[21:24])).zfill(3)
                    rowfiber['fiber_id'] = fiber_id_i
                    rowfiber['multiframe'] = multiframe
                    rowfiber['specid'] = multiframe[6:9]
                    rowfiber['ifuslot'] = multiframe[10:13]
                    rowfiber['ifuid'] = multiframe[14:17]
                    rowfiber['amp'] = multiframe[18:20]
                    rowfiber['fibnum'] = int(multiname[21:24])
                    rowfiber['distance'] = datafiber['col7'][ifiber]
                    rowfiber['wavein'] = datafiber['col8'][ifiber]
                    rowfiber['timestamp'] = datafiber['col9'][ifiber]
                    rowfiber['date'] = datafiber['col10'][ifiber]
                    rowfiber['obsid'] = str(datafiber['col11'][ifiber])[0:3]
                    rowfiber['x_raw'] = datafiber['col12'][ifiber]
                    rowfiber['y_raw'] = datafiber['col13'][ifiber]
                    rowfiber['flag'] = datafiber['col15'][ifiber]
                    rowfiber['weight'] = datafiber['col14'][ifiber]
                    rowfiber.append()
            except:
                args.log.error('Could not ingest %s' % filefiberinfo)

        tableFibers.flush()
        
    # create completely sorted index on the detectid 
    # to make queries against that column much faster
    if (args.append): 
        print ("Reindexing the detectid column")
        tableMain.cols.detectid.reindex()
        tableFibers.cols.detectid.reindex()
        tableSpectra.cols.detectid.reindex()
        tableFibers.flush() #just to be safe
        tableSpectra.flush()
        tableMain.flush()
    else:
        tableMain.cols.detectid.create_csindex()
        tableFibers.cols.detectid.create_csindex()
        tableSpectra.cols.detectid.create_csindex()
        tableFibers.flush() #just to be safe                                                            
        tableSpectra.flush()
        tableMain.flush()
    args.log.info('File finished: %s' % args.outfilename)             
    fileh.close()

if __name__ == '__main__':
    main()
