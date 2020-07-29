# -*- coding: utf-8 -*-
"""
Created: 2019/01/23

@author: Erin Mentuch Cooper

Script to gather information for the Survey Table. This provides
a master look up table for representaive values for each shot
in the HETDEX survey. Use the API provided in 
hetdex_api/survey.py to easily access data and query 
astropy coordinates.

Updates for HDR2:
-now ingesting as much as possible from the H5 shot files so
this just a table stack of all the HDR2 Shot tables

To run:

python3 create_survey_hdf5.py -of survey_hdr2.h5 -sl dex.hdr2.shotlist 

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

    if re.match("par", str(objname)):
        field = "parallel"
    elif re.match("COS|cos|DEXcos", str(objname)):
        field = "cosmos"
    elif re.match("EGS|DEXeg", str(objname)):
        field = "egs"
    elif re.match("GN|DEXgn", str(objname)):
        field = "goods-n"
    elif re.match("DEX0|DEXfl|HF", str(objname)):
        field = "dex-fall"
    elif re.match("HS|DEXsp", str(objname)):
        field = "dex-spring"
    elif re.match("NEP", str(objname)):
        field = "nep"
    elif re.match("SSA22", str(objname)):
        field = "ssa22"
    else:
        field = "other"
        
    return field


def main(argv=None):
    """ Main Function """
    # Call initial parser from init_utils
    parser = ap.ArgumentParser(
        description="""Create HDF5 Astrometry file.""", add_help=True
    )

    parser.add_argument(
        "-r",
        "--rootdir",
        help="""Root Directory for Reductions""",
        type=str,
        default="/work/03946/hetdex/maverick/red1/reductions/",
    )

    parser.add_argument(
        "-sdir",
        "--shotdir",
        help="""Directory for shot H5 files to ingest""",
        type=str,
        default="/data/05350/ecooper/hdr2.1/reduction/data",
    )

    parser.add_argument(
        "-sl",
        "--shotlist",
        help="""Text file of DATE OBS list""",
        type=str,
        default="dex.hdr2.shotlist",
    )

    parser.add_argument(
        "-ad",
        "--astrometry_dir",
        help="""Directory for Shifts""",
        type=str,
        default="/data/00115/gebhardt/vdrp/shifts/",
    )

    parser.add_argument(
        "-of",
        "--outfilename",
        type=str,
        help="""Relative or absolute path for output HDF5 file.""",
        default=None,
    )

    parser.add_argument(
        "-flim",
        "--flim",
        help="""Path to flim look up table""",
        type=str,
        default="/data/05350/ecooper/hdr2.1/survey/average_one_sigma.txt",
    )

    parser.add_argument("-survey", "--survey", type=str, default="hdr2.1")
    
    args = parser.parse_args(argv)

    print(args)
    
    args.log = setup_logging()

    fileh = tb.open_file(
        args.outfilename, mode="w", title=args.survey.upper() + "Survey file "
    )

    shotlist = Table.read(
        args.shotlist, format="ascii.no_header", names=["date", "obs"]
    )

    survey = Table()

    for shotrow in shotlist:
        datevshot = str(shotrow["date"]) + "v" + str(shotrow["obs"]).zfill(3)

        try:
            args.log.info('Ingesting ' + datevshot)
            file_obs = tb.open_file(op.join(args.shotdir, datevshot + ".h5"), "r")

            shottable = Table(file_obs.root.Shot.read())

            # updating field in survey file 
            shottable['field'] = define_field(str(shottable['objid'][0]))

            survey = vstack([survey, shottable])
            file_obs.close()
        except:
            args.log.error("Could not ingest %s" % datevshot)
    
    tableMain = fileh.create_table(fileh.root, "Survey", obj=survey.as_array())

    tableMain.flush()
    fileh.close()
    
if __name__ == "__main__":
    main()
