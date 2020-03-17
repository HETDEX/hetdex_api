# -*- coding: utf-8 -*-
"""
Created: 2019/01/23

@author: Erin Mentuch Cooper

This script can add calibration group info and data to
an exising HDF5 file or create a new HDF5 file

Example command line use:

python create_cal_hdf5.py -d 20181111 -o 21 -of cal_20181111v021.h5

or

python create_cal_hdf5.py -d 20181111 -o 21 -of 20181111v021.h5 --append


"""

import glob
import sys
import os
import os.path as op
import argparse as ap
import numpy as np
import tables as tb

import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.io import fits
from hetdex_api.input_utils import setup_logging


#class AmpToAmp(tb.IsDescription):
#    ampid = tb.StringCol(20)
#    ifuslot = tb.StringCol(3)
#    ifuid = tb.StringCol(3)
#    specid = tb.StringCol(3)
#    amp = tb.StringCol(2)
#    normdata = tb.Float32Col((2, 100))


def main(argv=None):
    ''' Main Function '''
    # Call initial parser from init_utils
    parser = ap.ArgumentParser(description="""Create HDF5 Astrometry file.""",
                               add_help=True)

    parser.add_argument("-d", "--date",
                        help='''Date, e.g., 20170321, YYYYMMDD''',
                        type=str, default=None)

    parser.add_argument("-o", "--observation",
                        help='''Observation number, "00000007" or "7"''',
                        type=str, default=None)

    parser.add_argument("-tp", "--tpdir",
                        help='''Directory for Throughput Info''',
                        type=str,
                        default='/data/00115/gebhardt/detect')
                        #default='/work/03946/hetdex/hdr1/reduction/throughput')

    parser.add_argument("-am", "--ampdir",
                        help='''Directory for Amp to Amp''',
                        type=str,
                        default='/work/00115/gebhardt/maverick/getampnorm/all/')

    parser.add_argument('-of', '--outfilename', type=str,
                        help='''Relative or absolute path for output HDF5
                        file.''', default=None)

    parser.add_argument('-a', '--append',
                        help='''Appending to existing shot HDF5 file.''',
                        action="count", default=0)

    args = parser.parse_args(argv)
    args.log = setup_logging()

    # Creates a new file if the "--append" option is not set or the file
    # does not already exist.

    does_exist = False
    if op.exists(args.outfilename) and args.append:
        fileh = tb.open_file(args.outfilename, 'a')
        args.log.info('Appending calibration info to %s'% args.outfilename)
        does_exist = True
    else:
        fileh = tb.open_file(args.outfilename, 'w')
        args.log.info('Writingcalibration info to %s'% args.outfilename)

    try: 
        fileh.remove_node(fileh.root.Calibration, recursive=True)
    except:
        args.log.info('Creating new Calibration group')

    group = fileh.create_group(fileh.root, 'Calibration',
                               'HETDEX Calibration Info')
#    groupAmpToAmp = fileh.create_group(group, 'AmpToAmp',
#                                       'Amp to amp Fiber Normalization')
    groupThroughput = fileh.create_group(group, 'Throughput',
                                         'Throughput Curve')

    # populate the AmpToAmp group with normalization curves and metadata

 #   amptoampfiles = glob.glob(op.join(args.ampdir, 'multi*.norm'))

 #   for filen in amptoampfiles:
 #       idx = filen.find('multi_')
 #       ampid = filen[idx:idx+20]
 #       try:
 #           data = ascii.read(filen, names=['wave', 'norm'])
 #           fileh.create_table(groupAmpToAmp, str(ampid), data.as_array())
 #       except:
 #           args.log.warning('Could not include %s' % filen)
    # populate Throughput group with throughput curves

    datevshot = str(args.date) + 'v' + str(args.observation.zfill(3))
    
    tpfile = op.join(args.tpdir,'tp', datevshot + 'sedtp_f.dat')
    
    try:
        tp_data = ascii.read(tpfile,names=['wavelength','throughput','tp_low', 'tp_high',
                                           'rat_poly', 'tp_gband'])
        tp_4540 = tp_data['throughput'][np.where(tp_data['wavelength'] == 4540.)][0]
        
        tp_array = fileh.create_table(groupThroughput, 'throughput', tp_data.as_array())
        tp_array.set_attr('filename', tpfile)
    except:
        args.log.warning('Could not include %s' % tpfile)

    tppngfile = op.join(args.tpdir,
                        str(args.date) + 'v' + str(args.observation.zfill(3)),
                        'res',
                        str(args.date) + 'v' +
                        str(args.observation.zfill(3)) + 'sedtpa.png')

    try:
        pngimarr = plt.imread(tppngfile)
        pngim = fileh.create_array(groupThroughput, 'tp_png', pngimarr)
        pngim.attrs['CLASS'] = 'IMAGE'
    except:
        args.log.warning('Could not include %s' % tppngfile)


    # add virus FWHM and response_4540 to the Shot table
    
    shottable = fileh.root.Shot
    fwhm_file =  op.join(args.tpdir,
                         str(args.date) + 'v' +
                         str(args.observation.zfill(3)),
                         'fwhm.out')
    try:
        fwhm_arr = np.loadtxt(fwhm_file)

        for row in shottable:
            row['fwhm_virus'] = fwhm_arr[0]
            row['fwhm_virus_err'] = fwhm_arr[1]
            row['nstars_fit_fwhm'] = int(fwhm_arr[2])
            row['response_4540'] = tp_4540
            row.update()

    except:
        args.log.error('Could not include cal info in shot table for %s' % datevshot) 


    fileh.close()

if __name__ == '__main__':
    main()
