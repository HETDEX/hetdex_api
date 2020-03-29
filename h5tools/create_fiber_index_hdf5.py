# -*- coding: utf-8 -*-
"""
Created: 2020/03/29

@author: Erin Mentuch Cooper

Script to gather information for the FiberIndex Table. This provides
a master look up table for of all fibers contained in a 
HETDEX data release. Use the FiberIndex Class API provided in 
hetdex_api/survey.py to easily access data and query 
astropy coordinates.

To run:

python3 create_fiber_index_hdf5.py -of fiber_index_hdr2.h5 

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


def main(argv=None):
    """ Main Function """
    # Call initial parser from init_utils
    parser = ap.ArgumentParser(
        description="""Create HDF5 Astrometry file.""", add_help=True
    )

    parser.add_argument(
        "-sdir",
        "--shotdir",
        help="""Directory for shot H5 files to ingest""",
        type=str,
        default="/data/05350/ecooper/hdr2/reduction/data",
    )

    parser.add_argument(
        "-sl",
        "--shotlist",
        help="""Text file of DATE OBS list""",
        type=str,
        default="dex.hdr2.shotlist",
    )

    parser.add_argument(
        "-of",
        "--outfilename",
        type=str,
        help="""Relative or absolute path for output HDF5
                        file.""",
        default=None,
    )

    parser.add_argument("-survey", "--survey", type=str, default="hdr2")

    args = parser.parse_args(argv)
    args.log = setup_logging()

    fileh = tb.open_file(
        args.outfilename, mode="w", title=args.survey.upper() + " Fiber Index file "
    )

    shotlist = Table.read(
        args.shotlist, format="ascii.no_header", names=["date", "obs"]
    )

    master_fiber = Table()

    for shotrow in shotlist:
        try:
            datevshot = str(shotrow["date"]) + "v" + str(shotrow["obs"]).zfill(3)
            file_obs = tb.open_file(op.join(args.shotdir, datevshot + ".h5"), "r")
            fibertable = Table(file_obs.root.Data.FiberIndex.read())
            master_fiber = vstack([fibertable, master_fiber])
            file_obs.close()
            args.log.info("Ingesting %s" % datevshot)
        except:
            args.log.error("could not ingest %s" % datevshot)

    tableMain = fileh.create_table(fileh.root, "FiberIndex", obj=master_fiber.as_array())
    fileh.close()


if __name__ == "__main__":
    main()
