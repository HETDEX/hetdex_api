# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 09:52:35 2017

@author: gregz
"""
import glob
import re

import tables as tb
import argparse as ap
import os.path as op
import numpy as np

from astropy.io import fits
from input_utils import setup_logging
from astropy.table import Table


def build_path(reduction_folder, instr, date, obsid, expn):
    folder = op.join(date, instr, "{:s}{:07d}".format(instr, int(obsid)),
                     "exp{:02d}".format(int(expn)), instr)
    return op.join(reduction_folder, folder)


def get_files(args):
    files = glob.glob(op.join(args.rootdir, args.date, 'virus',
                              'virus%07d' % int(args.observation),
                              'exp*', 'virus', 'multi_*.fits'))
    return files


class VIRUSFiber(tb.IsDescription):
    obsind = tb.Int32Col()
    multiframe = tb.StringCol((20), pos=0)
    fibnum = tb.Int32Col()
    ifux = tb.Float32Col()
    ifuy = tb.Float32Col()
    fpx = tb.Float32Col()
    fpy = tb.Float32Col()
    ra = tb.Float32Col(pos=1)
    dec = tb.Float32Col(pos=2)
    spectrum = tb.Float32Col((1032,))
    wavelength = tb.Float32Col((1032,))
    fiber_to_fiber = tb.Float32Col((1032,))
    twi_spectrum = tb.Float32Col((1032,))
    trace = tb.Float32Col((1032,))
    sky_subtracted = tb.Float32Col((1032,))
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


def append_shot_to_table(shot, fn, cnt):
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
    shot['expn'] = int(op.basename(op.dirname(op.dirname(F.filename())))[-2:])
    shot.append()


def append_fibers_to_table(fib, im, fn, cnt, T):
    F = fits.open(fn)
    idx = fn.find('multi')
    multiframe = fn[idx:idx+20]
    im['multiframe'] = multiframe
    n = F['spectrum'].data.shape[0]
    d = F['spectrum'].data.shape[1]
    attr = ['spectrum', 'wavelength', 'fiber_to_fiber', 'twi_spectrum',
            'sky_subtracted', 'trace', 'error1Dfib', 'calfib', 'calfibe',
            'Amp2Amp', 'Throughput']
    imattr = ['image', 'error', 'clean_image']
    for att in imattr:
        if att == 'image':
            im[att] = F['PRIMARY'].data * 1.
        else:
            im[att] = F[att].data * 1.
    mname = op.basename(fn)[:-5]
    expn = op.basename(op.dirname(op.dirname(fn)))
    if T is not None:
        sel = T['col8'] == (mname + '_001.ixy')
        sel1 = T['col10'] == expn
        loc = np.where(sel * sel1)[0]
    for i in np.arange(n):
        fib['obsind'] = cnt
        fib['fibnum'] = i
        if T is not None:
            fib['multiframe'] = multiframe
            loci = loc + i
            if len(loc):
                if isinstance(T['col1'][loci], float):
                    fib['ra'] = T['col1'][loci]
                    fib['dec'] = T['col2'][loci]
                    fib['fpx'] = T['col6'][loci]
                    fib['fpy'] = T['col7'][loci]
            else:
                fib['ra'] = -999.0
                fib['dec'] = -999.0
                fib['fpx'] = -999.0
                fib['fpy'] = -999.0
        fib['ifux'] = F['ifupos'].data[i, 0]
        fib['ifuy'] = F['ifupos'].data[i, 1]
        for att in attr:
            if att in F:
                fib[att] = F[att].data[i, :]
            else:
                fib[att] = np.zeros((d,))
        fib['ifuslot'] = '%03d' % int(F[0].header['IFUSLOT'])
        fib['ifuid'] = '%03d' % int(F[0].header['IFUID'])
        fib['specid'] = '%03d' % int(F[0].header['SPECID'])
        fib['contid'] = F[0].header['CONTID']
        fib['amp'] = '%s' % F[0].header['amp'][:2]
        fib['expnum'] = int(expn[-2:])
        fib.append()
    im['obsind'] = cnt
    im['ifuslot'] = '%03d' % int(F[0].header['IFUSLOT'])
    im['ifuid'] = '%03d' % int(F[0].header['IFUID'])
    im['specid'] = '%03d' % int(F[0].header['SPECID'])
    im['contid'] = F[0].header['CONTID']
    im['amp'] = '%s' % F[0].header['amp'][:2]
    im['expnum'] = int(expn[-2:])
    im.append()
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

    args = parser.parse_args(argv)
    args.log = setup_logging()

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

    # Grab the fiber table and amplifier table for writing
    fibtable = fileh.root.Data.Fibers
    shottable = fileh.root.Shot
    imagetable = fileh.root.Data.Images

    if does_exist:
        cnt = shottable[-1]['obsind']
    else:
        cnt = 1

    shot = shottable.row
    success = append_shot_to_table(shot, files[0], cnt)
    if success:
        shottable.flush()
    for fn in files:
        args.log.info('Working on %s' % fn)
        fib = fibtable.row
        im = imagetable.row
        success = append_fibers_to_table(fib, im, fn, cnt, T)
        if success:
            fibtable.flush()
            imagetable.flush()

    # create completely sorted index on the specid to make queries against that column much faster
    # specid chosen as the old multi*fits naming started with specid and it is fixed vs ifuslot and ifuid
    # for any given shot
    fibtable.cols.specid.create_csindex()
    imagetable.cols.specid.create_csindex()
    fibtable.flush()
    imagetable.flush()

    fileh.close()


if __name__ == '__main__':
    main()
