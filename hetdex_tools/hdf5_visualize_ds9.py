# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 16:11:17 2019

@author: gregz
"""

import argparse as ap
import numpy as np
import os.path as op
import pyds9

from input_utils import setup_logging
from tables import open_file


def main(argv=None):
    ''' Main Function '''
    # Call initial parser from init_utils
    parser = ap.ArgumentParser(description="""Create HDF5 file.""",
                               add_help=True)

    parser.add_argument("hdf5_file", type=str,
                        help='''Name of h5/hdf5 file''')

    parser.add_argument("-e", "--extension",
                        help='''Extension or extensions''',
                        type=str, default=None)

    parser.add_argument("-q", "--query",
                        help='''Query to be applied''',
                        type=str, default=None)

    parser.add_argument('-s', '--show',
                        help='''Show tables/extensions within file''',
                        action="count", default=0)

    args = parser.parse_args(argv)
    args.log = setup_logging()
    try:
        h5file = open_file(args.hdf5_file, mode='r')
    except:
        if not op.exists(args.hdf5_file):
            args.log.error('%s does not exist' % args.hdf5_file)
            return None
        else:
            args.log.error('File exists but could not open %s'
                           % args.hdf5_file)
            return None
    if args.show:
        table_names = ['Shot', 'Fibers', 'Images']
        for kind in table_names:
            print('%s column names:' % kind)
            b = getattr(h5file.root.Info, kind)
            for name in b.colnames:
                base = getattr(b.cols, name)
                shape = str(base.shape)
                print('\t%s: %s %s' % (name, base.type, shape))
        return None

    if args.extension is None:
        args.log.error('No extension provided to display in ds9')
        return None

    table = h5file.root.Info.Fibers
    spec = np.array(table.cols.spectrum[:])
    ds9 = pyds9.DS9()
    ds9.set_np2arr(spec)
    h5file.close()

if __name__ == '__main__':
    main()
