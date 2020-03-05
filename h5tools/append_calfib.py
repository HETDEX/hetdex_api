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


def get_cal_files(args):
    datestr = '%sv%03d' % (args.date, int(args.observation))

    files = glob.glob(op.join(args.rootdir, datestr + '*cal.fits')) 

    return files


def get_cal_table(calfile):
    
    cal_f = fits.open(calfile)

    cal_table = Table([cal_f[0].data, cal_f[1].data, cal_f[2].data, cal_f[3].data], 
                      names=['calfib','calfibe', 'calfib_counts','calfibe_counts'])

    multi = calfile[42:53]

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

    parser.add_argument("-survey", "--survey", help='''{hdr1, hdr2, hdr3}''',
                        type=str, default='hdr2')

    args = parser.parse_args(argv)
    args.log = setup_logging()

    calfiles = get_cal_files(args)

    datestr = '%sv%03d' % (args.date, int(args.observation))

    if op.exists(args.outfilename):
        fileh = tb.open_file(args.outfilename, 'a')
    else:
        args.log.warning('Problem opening : ' + args.outfilename)
        sys.exit('Exiting Script')
    
    args.log.info('Appending calibrated fiber arrays to ' + args.outfilename)

    fibtable = fileh.root.Data.Fibers

    
    for calfile in calfiles:
    
        multi = calfile[42:53]
        cal_table = get_cal_table(calfile)
        args.log.info('Working on IFU ' + multi )
        for amp_i in ['LL','LU','RL','RU']:
            
            multiframe_i = 'multi_'+ multi + '_' + amp_i
        
            for fibrow in fibtable.where('multiframe == multiframe_i'):
                idx = (cal_table['expnum'] == fibrow['expnum']) * (cal_table['multiframe'] == fibrow['multiframe'].decode()) * (cal_table['fibidx'] == fibrow['fibidx'])
                fibrow['calfib']  = cal_table['calfib'][idx]
                fibrow['calfibe'] = cal_table['calfibe'][idx]
                fibrow['calfib_counts'] = cal_table['calfib_counts'][idx]
                fibrow['calfibe_counts'] = cal_table['calfibe_counts'][idx]
                fibrow.update()
    args.log.info('Flushing and closing H5 file')
    fibtable.flush()
    fileh.close()

if __name__ == '__main__':
    main()
