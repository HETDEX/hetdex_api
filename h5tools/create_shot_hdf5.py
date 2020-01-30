# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 09:52:35 2017

@author: gregz
"""
import glob
import re
import os
import shutil
import tarfile
import sys

import tables as tb
import argparse as ap
import os.path as op
import numpy as np

from astropy.io import fits
from input_utils import setup_logging
from astropy.table import Table

# hard coded variable to initialize 'rms', 'chi' arrays
# and remove 'twi_spectrum' for all realeases past hdr1
global hdr_survey
hdr_survey = 'hdr2'

def build_path(reduction_folder, instr, date, obsid, expn):
    folder = op.join(date, instr, "{:s}{:07d}".format(instr, int(obsid)),
                     "exp{:02d}".format(int(expn)), instr)
    return op.join(reduction_folder, folder)


def get_files(args):
    if args.tar == False:
        files = glob.glob(op.join(args.rootdir, args.date, 'virus',
                                  'virus%07d' % int(args.observation),
                                  'exp*', 'virus', 'multi_*.fits'))
    else:

        datestr = 'd%ss%03d' % (args.date, int(args.observation))

        tmppath = 'tmp'
        # remove any old temporary multifits
        try:
            os.mkdir(tmppath)
        except FileExistsError:
            pass
            
        datepath = op.join(tmppath, datestr)
        if op.isdir(datepath):
            shutil.rmtree(datepath, ignore_errors=True)

        os.mkdir(datepath)

        tarfiles = glob.glob(op.join(args.rootdir,
                                     'sci' + str(args.date)[0:6],
                                     datestr + 'exp0?',
                                     datestr + 'exp0?_mu.tar'))

        if np.size(tarfiles) == 3:
            for exp in ['exp01', 'exp02', 'exp03']:
                expfile = op.join(args.rootdir,
                                  'sci' + str(args.date)[0:6],
                                  datestr + exp,
                                  datestr + exp + '_mu.tar')
                tar = tarfile.open(name=expfile, mode='r')
                dir_exp = op.join(datepath, exp)

                os.mkdir(dir_exp)
                tar.extractall(path=dir_exp)

            files = glob.glob(op.join(datepath, 'exp*/multi*.fits'))

        else:
            args.log.error('Could not locate tar files for three dithers.')
            sys.exit()

    return files

class VIRUSFiberIndex(tb.IsDescription):
    multiframe = tb.StringCol((20), pos=0)
    fiber_id = tb.StringCol((38), pos=4)
    fibidx = tb.Int32Col()
    ifux = tb.Float32Col()
    ifuy = tb.Float32Col()
    fpx = tb.Float32Col()
    fpy = tb.Float32Col()
    ra = tb.Float32Col(pos=1)
    dec = tb.Float32Col(pos=2)

class VIRUSFiber(tb.IsDescription):
    obsind = tb.Int32Col()
    multiframe = tb.StringCol((20), pos=0)
    fiber_id = tb.StringCol((38), pos=4)
    fibidx = tb.Int32Col()
    ifux = tb.Float32Col()
    ifuy = tb.Float32Col()
    fpx = tb.Float32Col()
    fpy = tb.Float32Col()
    ra = tb.Float32Col(pos=1)
    dec = tb.Float32Col(pos=2)
    spectrum = tb.Float32Col((1032,))
    wavelength = tb.Float32Col((1032,))
    fiber_to_fiber = tb.Float32Col((1032,))
    global hdr_survey
    if hdr_survey == 'hdr1':
        twi_spectrum = tb.Float32Col((1032,))
    else:
        chi2 = tb.Float32Col((1032,))
        rms = tb.Float32Col((1032,))
    trace = tb.Float32Col((1032,))
    sky_subtracted = tb.Float32Col((1032,))
    sky_spectrum = tb.Float32Col((1032,))
    error1Dfib = tb.Float32Col((1032,))
    calfib = tb.Float32Col((1036,))
    calfibe = tb.Float32Col((1036,))
    Amp2Amp = tb.Float32Col((1036,))
    Throughput = tb.Float32Col((1036,))
    ifuslot = tb.StringCol(3)
    ifuid = tb.StringCol(3)
    specid = tb.StringCol(3)
    contid = tb.StringCol(8)
    amp = tb.StringCol(2)
    expnum = tb.Int32Col()
    

class VIRUSImage(tb.IsDescription):
    obsind = tb.Int32Col()
    multiframe = tb.StringCol((20), pos=0)
    image = tb.Float32Col((1032, 1032))
    error = tb.Float32Col((1032, 1032))
    clean_image = tb.Float32Col((1032, 1032))
    ifuslot = tb.StringCol(3)
    ifuid = tb.StringCol(3)
    specid = tb.StringCol(3)
    contid = tb.StringCol(8)
    amp = tb.StringCol(2)
    expnum = tb.Int32Col()


class VIRUSShot(tb.IsDescription):
    obsind = tb.Int32Col()
    objid = tb.StringCol((18), pos=4)
    date = tb.Int32Col(pos=0)
    mjd = tb.Float32Col(pos=6)
    obsid = tb.Int32Col(pos=1)
    ra = tb.Float32Col(pos=2)
    dec = tb.Float32Col(pos=3)
    pa = tb.Float32Col(pos=5)
    expn = tb.Int32Col()
    time = tb.StringCol(7)
    ambtemp = tb.Float32Col()
    humidity = tb.Float32Col()
    dewpoint = tb.Float32Col()
    pressure = tb.Float32Col()
    exptime = tb.Float32Col()


def append_shot_to_table(shot, shottable, fn, cnt):
    F = fits.open(fn)
    shot['obsind'] = cnt
    shot['date'] = int(''.join(F[0].header['DATE-OBS'].split('-')))
    shot['objid'] = F[0].header['OBJECT']
    shot['time'] = ''.join(re.split('[:,.]', F[0].header['UT']))[:7]
    shot['mjd'] = F[0].header['MJD']
    shot['obsid'] = int(F[0].header['OBSID'])
    shot['ra'] = F[0].header['TRAJCRA'] * 15.
    shot['dec'] = F[0].header['TRAJCDEC']
    shot['pa'] = F[0].header['PARANGLE']
    shot['ambtemp'] = F[0].header['AMBTEMP']
    shot['humidity'] = F[0].header['HUMIDITY']
    shot['dewpoint'] = F[0].header['DEWPOINT']
    shot['pressure'] = F[0].header['BAROMPRE']
    shot['exptime'] = F[0].header['EXPTIME']
    try:
        shot['expn'] = int(op.basename(op.dirname(op.dirname(F.filename())))[-2:])
    except:
        shot['expn'] = int(op.dirname(F.filename())[-2:])
    shottable.attrs['HEADER'] = F[0].header
    shot.append()


def append_fibers_to_table(fibindex, fib, im, fn, cnt, T, args):
    F = fits.open(fn)
    shotid = int(args.date) * 1000 + int(args.observation) 
    idx = fn.find('multi')
    multiframe = fn[idx:idx+20]
    im['multiframe'] = multiframe
    n = F['spectrum'].data.shape[0]
    d = F['spectrum'].data.shape[1]
    if args.survey == 'hdr1':
        attr = ['spectrum', 'wavelength', 'fiber_to_fiber', 'twi_spectrum',
                'sky_spectrum','sky_subtracted', 'trace', 'error1Dfib',
                'calfib', 'calfibe','Amp2Amp', 'Throughput']
        imattr = ['image', 'error', 'clean_image']
    else:
        attr = ['spectrum', 'wavelength', 'fiber_to_fiber',
                'sky_spectrum','sky_subtracted', 'trace', 'error1Dfib',
                'calfib', 'calfibe','Amp2Amp', 'Throughput','chi2','rms']
        imattr = ['image', 'error', 'clean_image']

    for att in imattr:
        if att == 'image':
            im[att] = F['PRIMARY'].data * 1.
        else:
            im[att] = F[att].data * 1.
    mname = op.basename(fn)[:-5]

    if args.tar:
        expn = op.basename((op.dirname(fn)))
    else:
        expn = op.basename(op.dirname(op.dirname(fn)))

    if T is not None:
        sel = T['col8'] == (mname + '_001.ixy')
        sel1 = T['col10'] == expn
        loc = np.where(sel * sel1)[0]
    for i in np.arange(n):
        fib['obsind'] = cnt
        fib['fibidx'] = i
        fib['multiframe'] = multiframe
        fibindex['multiframe'] = fib['multiframe']
        fib['fiber_id'] = str(shotid) + '_' + str(int(expn[-2:])) + '_' + multiframe + '_' + str(i+1).zfill(3)
        fibindex['fiber_id'] = fib['fiber_id']

        if T is not None:
            if len(loc):
                loci = loc[0] + i
                if isinstance(T['col1'][loci], float):
                    fib['ra'] = T['col1'][loci]
                    fib['dec'] = T['col2'][loci]
                    fib['fpx'] = T['col6'][loci]
                    fib['fpy'] = T['col7'][loci]
                else:
                    fib['ra'] = np.nan
                    fib['dec'] = np.nan
                    fib['fpx'] = np.nan
                    fib['fpy'] = np.nan
            else:
                fib['ra'] = np.nan
                fib['dec'] = np.nan
                fib['fpx'] = np.nan
                fib['fpy'] = np.nan
        else:
            fib['ra'] = np.nan
            fib['dec'] = np.nan
            fib['fpx'] = np.nan
            fib['fpy'] = np.nan
        
        fibindex['ra'] = fib['ra']
        fibindex['dec'] = fib['dec']
        fibindex['fpx'] = fib['fpx']
        fibindex['fpy'] = fib['fpy']
        
        fib['ifux'] = F['ifupos'].data[i, 0]
        fib['ifuy'] = F['ifupos'].data[i, 1]
        fibindex['ifux'] = fib['ifux']
        fibindex['ifuy'] = fib['ifuy']

        for att in attr:
            if att in F:
                fib[att] = F[att].data[i, :]
        fib['ifuslot'] = '%03d' % int(F[0].header['IFUSLOT'])
        fib['ifuid'] = '%03d' % int(F[0].header['IFUID'])
        fib['specid'] = '%03d' % int(F[0].header['SPECID'])
        fib['contid'] = F[0].header['CONTID']
        try:
            fib['amp'] = '%s' % F[0].header['amp'][:2]
        except:
            fib['amp'] = '%s' % F.filename()[-7:-5]
        
        fib['expnum'] = int(expn[-2:])
        fib.append()
        fibindex.append()

    im['obsind'] = cnt
    im['ifuslot'] = '%03d' % int(F[0].header['IFUSLOT'])
    im['ifuid'] = '%03d' % int(F[0].header['IFUID'])
    im['specid'] = '%03d' % int(F[0].header['SPECID'])
    im['contid'] = F[0].header['CONTID']
    try:
        im['amp'] = '%s' % F[0].header['amp'][:2]
    except:
        im['amp'] = '%s' % F.filename()[-7:-5]
    im['expnum'] = int(expn[-2:])
    im.append()

    #close the fits file
    F.close()
    return True


def main(argv=None):
    ''' Main Function '''
    # Call initial parser from init_utils
    parser = ap.ArgumentParser(description="""Create HDF5 file.""",
                               add_help=True)

    parser.add_argument("-d", "--date",
                        help='''Date, e.g., 20170321, YYYYMMDD''',
                        type=str, default=None)

    parser.add_argument("-o", "--observation",
                        help='''Observation number, "00000007" or "7"''',
                        type=str, default=None)

    parser.add_argument("-r", "--rootdir",
                        help='''Root Directory for Reductions''',
                        type=str, default='/work/03946/hetdex/maverick/red1/reductions/')

    parser.add_argument('-of', '--outfilename', type=str,
                        help='''Relative or absolute path for output HDF5
                        file.''', default=None)
    
    parser.add_argument('-a', '--append',
                        help='''Appending to existing file.''',
                        action="count", default=0)
    
    parser.add_argument("-dp", "--detect_path",
                        help='''Path to detections''',
                        type=str, default='/work/00115/gebhardt/maverick/detect')

    parser.add_argument("-survey", "--survey", help='''{hdr1, hdr2, hdr3}''',
                        type=str, default='hdr2')
    
    parser.add_argument("-tar", "--tar", help='''Flag to open tarred multifits''',
                        action='store_true')

    args = parser.parse_args(argv)
    args.log = setup_logging()

    global hdr_survey
    hdr_survey = args.survey
    if args.survey != hdr_survey:
        args.log.warning('Hard coded hdr_survey does not match input survey.')
        sys.exit()
        
    # Get the daterange over which reduced files will be collected
    files = get_files(args)
    datestr = '%sv%03d' % (args.date, int(args.observation))
    filepath = '%s/%s/dithall.use' % (args.detect_path, datestr)
    try:
        T = Table.read(filepath, format='ascii')
    except:
        T = None

    # Creates a new file if the "--append" option is not set or the file
    # does not already exist.
    does_exist = False
    if op.exists(args.outfilename) and args.append:
        fileh = tb.open_file(args.outfilename, 'a')
        does_exist = True
    else:
        fileh = tb.open_file(args.outfilename, 'w')
        group = fileh.create_group(fileh.root, 'Data',
                                   'VIRUS Fiber Data and Metadata')
        fileh.create_table(group, 'Fibers', VIRUSFiber, 'Fiber Info')
        fileh.create_table(fileh.root, 'Shot', VIRUSShot, 'Shot Info')
        fileh.create_table(group, 'Images', VIRUSImage, 'Image Info')
        fileh.create_table(group, 'FiberIndex', VIRUSFiberIndex, 'Fiber Coord Info')

    # Grab the fiber table and amplifier table for writing
    fibtable = fileh.root.Data.Fibers
    shottable = fileh.root.Shot
    imagetable = fileh.root.Data.Images
    fibindextable = fileh.root.Data.FiberIndex

    if does_exist:
        cnt = shottable[-1]['obsind']
    else:
        cnt = 1

    shot = shottable.row
    success = append_shot_to_table(shot, shottable, files[0], cnt)
    if success:
        shottable.flush()
    for fn in files:
        args.log.info('Working on %s' % fn)
        fib = fibtable.row
        im = imagetable.row
        fibindex = fibindextable.row
 
        success = append_fibers_to_table(fibindex, fib, im, fn, cnt, T, args)
        if success:
            fibtable.flush()
            imagetable.flush()

    # create completely sorted index on the specid to make queries against that column much faster
    # specid chosen as the old multi*fits naming started with specid and it is fixed vs ifuslot and ifuid
    # for any given shot
    fibtable.cols.multiframe.create_csindex()
    imagetable.cols.multiframe.create_csindex()
    fibtable.flush()
    imagetable.flush()

    fileh.close()

    # remove all temporary multifits
    if args.tar:
        datestr = 'd%ss%03d' % (args.date, int(args.observation))
        tmppath = 'tmp'
        datepath = op.join(tmppath, datestr)
        shutil.rmtree(datepath, ignore_errors=True)

if __name__ == '__main__':
    main()
