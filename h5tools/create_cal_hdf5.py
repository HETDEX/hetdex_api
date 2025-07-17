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

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table

from hetdex_api.input_utils import setup_logging
from hetdex_api.config import HDRconfig

import traceback

def main(argv=None):
    """ Main Function """
    # Call initial parser from init_utils
    parser = ap.ArgumentParser(
        description="""Append HDF5 Calibration info table.""", add_help=True
    )

    parser.add_argument(
        "-d",
        "--date",
        help="""Date, e.g., 20170321, YYYYMMDD""",
        type=str,
        default=None,
    )

    parser.add_argument(
        "-o",
        "--observation",
        help='''Observation number, "00000007" or "7"''',
        type=str,
        default=None,
    )

    parser.add_argument(
        "-tp",
        "--tpdir",
        help="""Directory for Throughput Info""",
        type=str,
        default="/scratch/00115/gebhardt/detect",
    )

    parser.add_argument(
        "-detdir",
        "--detectdir",
        help="""Directory for Detect Info""",
        type=str,
        default="/scratch/projects/hetdex/detect",
    )
        
    parser.add_argument(
        "-of",
        "--outfilename",
        type=str,
        help="""Relative or absolute path for output HDF5file.""",
        default=None,
    )

    parser.add_argument(
        "-a",
        "--append",
        help="""Appending to existing shot HDF5 file.""",
        action="count",
        default=0,
    )

    parser.add_argument("-survey", "--survey",
                        help="""{hdr1, hdr2, hdr2.1, hdr3}""",
                        type=str, default="hdr3")
    
                    
    args = parser.parse_args(argv)
    args.log = setup_logging()

    # Creates a new file if the "--append" option is not set or the file
    # does not already exist.

    does_exist = False
    if op.exists(args.outfilename) and args.append:
        fileh = tb.open_file(args.outfilename, "a")
        args.log.info("Appending calibration info to %s" % args.outfilename)
        does_exist = True
    else:
        fileh = tb.open_file(args.outfilename, "w")
        args.log.info("Writingcalibration info to %s" % args.outfilename)

    shotid = int(str(args.date) + str(args.observation).zfill(3))
    
    #check if shotid is in badlist
    
    config = HDRconfig(args.survey)
    badshots = np.loadtxt(config.badshot, dtype=int)
    
    badshotflag = False
    
    if shotid in badshots:
        badshotflag = True
        
    try:
        fileh.remove_node(fileh.root.Calibration, recursive=True)
    except:
        args.log.info("Creating new Calibration group")

    group = fileh.create_group(fileh.root, "Calibration", "HETDEX Calibration Info")

    groupThroughput = fileh.create_group(group, "Throughput", "Throughput Curve")

    datevshot = str(args.date) + "v" + str(args.observation.zfill(3))

    tpfile = op.join(args.detectdir, "tp", datevshot + "sedtp_f.dat")
    if not op.exists(tpfile):
        #try alternate location
        tpfile = op.join(args.tpdir, "tp", datevshot + "sedtp_f.dat")
        #if this is bad, will fail below

    try:
        tp_data = ascii.read(
            tpfile,
            names=[
                "wavelength",
                "throughput",
                "tp_low",
                "tp_high",
                "rat_poly",
                "tp_gband",
            ],
        )
        tp_4540 = tp_data["throughput"][np.where(tp_data["wavelength"] == 4540.0)][0]

        tp_array = fileh.create_table(groupThroughput, "throughput", tp_data.as_array())
        tp_array.set_attr("filename", tpfile)
    except:
        args.log.warning("Could not include %s" % tpfile)

    tppngfile = op.join(
        args.tpdir,
        str(args.date) + "v" + str(args.observation.zfill(3)),
        "res",
        str(args.date) + "v" + str(args.observation.zfill(3)) + "sedtpa.png",
    )

    try:
        pngimarr = plt.imread(tppngfile)
        pngim = fileh.create_array(groupThroughput, "tp_png", pngimarr)
        pngim.attrs["CLASS"] = "IMAGE"
    except:
        args.log.warning("Could not include %s" % tppngfile)

    # add virus FWHM and response_4540 to the Shot table

    shottable = fileh.root.Shot
    fwhm_file = op.join( args.detectdir, "fwhm.all")

    try:
        
        fwhm_tab = Table.read( fwhm_file, format='ascii.no_header')
        sel_datevobs = fwhm_tab['col1'] == str(args.date) + "v" + str(args.observation.zfill(3))
        
        for row in shottable:
            row["fwhm_virus"] = float(fwhm_tab['col2'][sel_datevobs][0])
            row["fwhm_virus_err"] = float(fwhm_tab['col3'][sel_datevobs][0])
            row["nstars_fit_fwhm"] = int(fwhm_tab['col4'][sel_datevobs][0])
            row["response_4540"] = tp_4540
            row.update()

    except:
        #temp extra logging
        args.log.error(traceback.format_exc())

        if badshotflag:
            args.log.warning("Could not include cal info in shot table for %s" % datevshot)
        else:
            args.log.error("Could not include cal info in shot table for %s" % datevshot)

    fileh.close()


if __name__ == "__main__":
    main()
