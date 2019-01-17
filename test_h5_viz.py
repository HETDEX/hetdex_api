# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 16:11:17 2019

@author: gregz
"""

import argparse as ap
from input_utils import setup_logging
from tables import open_file
import numpy as np
import pyds9


h5file = open_file('/work/03730/gregz/maverick/test.h5', mode='r')
table = h5file.root.Info.Fibers
spec = np.array(table.cols.spectrum[:])
ds9 = pyds9.DS9()
ds9.set_np2arr(spec)
raw_input('Done Inspecting?')
h5file.close()

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
                        type=str, default='/work/03946/hetdex/maverick')

    parser.add_argument('-of', '--outfilename', type=str,
                        help='''Relative or absolute path for output HDF5
                        file.''', default=None)
    parser.add_argument('-a', '--append',
                        help='''Appending to existing file.''',
                        action="count", default=0)

    args = parser.parse_args(argv)
    args.log = setup_logging()