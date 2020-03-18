# -*- coding: utf-8 -*-
"""
Created: 2019/01/23

@author: Erin Mentuch Cooper

Script to gather information for the Survey Table. This provides
a master look up table for representaive values for each shot
in the HETDEX survey. Use the API provided in 
hetdex_api/survey.py to easily access data and query 
astropy coordinates.


The script uses the files to generate the table, but any
list of DATE OBS would work.

/work/03946/hetdex/hdr1/reduction/hdr1.scilist 
/work/03946/hetdex/hdr1/reduction/hdr1.callist 

Requires multi*fits files to gather header info, a shifts/
directory to gather some QA astrometry data, and additional
paths contained in the HETDEX_API/config.py files

Updates for HDR2:
-now ingesting as much as possible from the H5 shot files
-removed TRAJCRA, TRAJCDEC, TRAJPA, STRUCTAZ

To run:

pyton create_survey_hdf5.py -of survey_hdr1.h5

"""

import re
import os
import os.path as op
import numpy as np
import tables as tb
import argparse as ap

from astropy.table import vstack, Table
from hetdex_api.input_utils import setup_logging

from hetdex_api.shot import open_shot_file
from hetdex_api import config


def define_field(objname):
    if re.match('par', str(objname)):
        field = 'parallel'
    elif re.match('COS|cos|DEXcos', str(objname)):
        field = 'cosmos'
    elif re.match('EGS', str(objname)):
        field = 'egs'
    elif re.match('GN', str(objname)):
        field = 'goods-n'
    elif re.match('DEX0|DEXfl', str(objname)):
        field = 'dex-fall'
    elif re.match('HS|DEXsp', str(objname)):
        field = 'dex-spring'
    else:
        field = 'other'

    return field


def main(argv=None):
    ''' Main Function '''
    # Call initial parser from init_utils
    parser = ap.ArgumentParser(description="""Create HDF5 Astrometry file.""",
                               add_help=True)

    parser.add_argument("-r", "--rootdir",
                        help='''Root Directory for Reductions''',
                        type=str, default='/work/03946/hetdex/maverick/red1/reductions/')

    parser.add_argument("-sdir", "--shotdir",
                        help='''Directory for shot H5 files to ingest''',
                        type=str, default='/data/05350/ecooper/hdr2/reduction/data')

    parser.add_argument("-sl", "--shotlist", 
                        help='''Text file of DATE OBS list''',
                        type=str, default='hdr2_shotlist')

    parser.add_argument("-ad", "--astrometry_dir",
                        help='''Directory for Shifts''',
                        type=str, default='/data/00115/gebhardt/vdrp/shifts/')

    parser.add_argument('-of', '--outfilename', type=str,
                        help='''Relative or absolute path for output HDF5
                        file.''', default=None)

    parser.add_argument("-flim", "--flim_dir",
                        help='''Path to flim look up table''',
                        type=str, default='/work/04120/dfarrow/wrangler/flims/hdr1/average_flims_4500_4600.txt')

    args = parser.parse_args(argv)
    args.log = setup_logging()

    fileh = tb.open_file(args.outfilename, mode="w", title= args.survey.upper() + "Survey file ")
    
    shotlist = Table.read(args.shotlist, format='ascii.no_header', names=['date','obs'])

    for shotrow in shotlist:
        try:
            datevshot = str(shotrow['date']) + 'v' + str(shotrow['obs']).zfill(3)
            file_obs = tb.open_file(op.join(args.shotdir, datevshot + '.h5' ), 'r')
            shottable = Table(file_obs.root.Shot.read())
            survey = vstack([survey, shottable])
            file_obs.close()
        except:
            args.log.Error('could not ingest %s' % datevshot)

    tableMain = fileh.create_table(fileh.root, 'Survey', obj=survey.as_array())
    fileh.close()


if __name__ == '__main__':
    main()
