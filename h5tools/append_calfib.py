# -*- coding: utf-8 -*-                                                                       
"""                                                                                                                 
Created on January 30, 2020

Script to append calfib fits files to Fibers class
                                                                                                                    
@author: Erin Mentuch Cooper                                                                  
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
from hetdex_api.input_utils import setup_logging
from astropy.table import Table, join
from hetdex_api.config import HDRconfig


def get_cal_files(args):
    datestr = '%sv%03d' % (args.date, int(args.observation))

    datemonth = datestr[0:6]
    files = glob.glob(op.join(args.rootdir, datemonth, datestr + '*cal.fits')) 

    return files


def get_cal_table(calfile):
    
    cal_f = fits.open(calfile)

    cal_table = Table([cal_f[0].data, cal_f[1].data, cal_f[2].data,
                       cal_f[3].data, cal_f[4].data],
                      names=['calfib','calfibe', 'calfib_counts',
                             'calfibe_counts','spec_fullsky_sub'])
    multi  = calfile[49:60]

    amp_col = []
    exp_col = []
    fib_col = []
    multi_col = []

    for exp_i in [1,2,3]:
        for amp_i in ['LL','LU','RL','RU']:
            for fib_i in np.arange(0,112):
                amp_col.append(amp_i)
                exp_col.append(exp_i)
                fib_col.append(fib_i)
                multi_col.append('multi_'+ multi + '_' + amp_i)

    cal_table['amp'] = np.array(amp_col)
    cal_table['expnum'] = np.array(exp_col)
    cal_table['fibidx'] = np.array(fib_col)
    cal_table['multiframe'] = np.array(multi_col)
    
    cal_f.close()

    return cal_table

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
                        type=str, default='/data/00115/gebhardt/calfits/')

    parser.add_argument('-of', '--outfilename', type=str,
                        help='''Relative or absolute path for output HDF5                          
                        file.''', default=None)

    parser.add_argument("-survey", "--survey", help='''{hdr1, hdr2, hdr2.1}''',
                        type=str, default='hdr2.1')

    args = parser.parse_args(argv)
    args.log = setup_logging()

    calfiles = get_cal_files(args)

    datestr = '%sv%03d' % (args.date, int(args.observation))

    shotid = int(str(args.date) + str(args.observation).zfill(3))

    #check if shotid is in badlist
    config = HDRconfig(args.survey)
    badshots = np.loadtxt(config.badshot, dtype=int)
    
    badshotflag = False
    
    if shotid in badshots:
        badshotflag = True
    
    if len(calfiles) == 0:
        if badshotflag:
            args.log.warning("No calfits file to append for %s" % datestr)
        else:
            args.log.error("No calfits file to append for %s" % datestr)

        sys.exit('Exiting cal append script for %s' % datestr)

    if op.exists(args.outfilename):
        fileh = tb.open_file(args.outfilename, 'a')
    else:
        args.log.error('Problem opening : ' + args.outfilename)
        sys.exit('Exiting Script')
    
    args.log.info('Appending calibrated fiber arrays to ' + args.outfilename)

    fibtable = fileh.root.Data.Fibers
    
    for calfile in calfiles:

        multi  = calfile[49:60]
        try:
            cal_table = get_cal_table(calfile)
        except:
            continue
            args.log.error('Could not ingest calfile: %s' % calfile)
            
        args.log.info('Working on IFU ' + multi )
        for amp_i in ['LL','LU','RL','RU']:
            
            multiframe_i = 'multi_'+ multi + '_' + amp_i

            for fibrow in fibtable.where('multiframe == multiframe_i'):
                
                idx = (cal_table['expnum'] == fibrow['expnum']) * (cal_table['multiframe'] == fibrow['multiframe'].decode()) * (cal_table['fibidx'] == fibrow['fibidx'])

                if np.sum(idx) >= 1:
                    fibrow['calfib']  = cal_table['calfib'][idx]
                    fibrow['calfibe'] = cal_table['calfibe'][idx]
                    fibrow['calfib_counts'] = cal_table['calfib_counts'][idx]
                    fibrow['calfibe_counts'] = cal_table['calfibe_counts'][idx]
                    fibrow['spec_fullsky_sub'] = cal_table['spec_fullsky_sub'][idx]
                    fibrow.update()
                #else:
                   # args.log.warning("No fiber match for %s" % fibrow['fiber_id'])
                    
    args.log.info('Flushing and closing H5 file')
    fibtable.flush()
    fileh.close()

if __name__ == '__main__':
    main()
