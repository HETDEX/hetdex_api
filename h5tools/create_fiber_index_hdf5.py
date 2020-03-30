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
import healpy as hp

from astropy.table import vstack, Table
from hetdex_api.input_utils import setup_logging

from hetdex_api.shot import open_shot_file
from hetdex_api.config import HDRconfig


class VIRUSFiberIndex(tb.IsDescription):
    multiframe = tb.StringCol((20), pos=0)
    fiber_id = tb.StringCol((38), pos=4)
    shotid = tb.Int64Col()
    healpix = tb.Int64Col(pos=5)
    ifuslot = tb.StringCol(3)
    ifuid = tb.StringCol(3)
    specid = tb.StringCol(3)
    amp = tb.StringCol(2)
    fibidx = tb.Int32Col()
    fibnum = tb.Int32Col()
    ifux = tb.Float32Col()
    ifuy = tb.Float32Col()
    fpx = tb.Float32Col()
    fpy = tb.Float32Col()
    ra = tb.Float32Col(pos=1)
    dec = tb.Float32Col(pos=2)
    expnum = tb.Int32Col()


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

    tableFibers = fileh.create_table(
        fileh.root,
        "FiberIndex",
        VIRUSFiberIndex,
        "Survey Fiber Coord Info",
        expectedrows=140369264,
    )

    # set up HEALPIX options
    Nside = 2 ** 17
    hp.max_pixrad(Nside, degrees=True) * 3600  # in unit of arcsec

    config = HDRconfig(survey=args.survey)

    badshot = np.loadtxt(config.badshot, dtype=int)

    for shotrow in shotlist:
        datevshot = str(shotrow["date"]) + "v" + str(shotrow["obs"]).zfill(3)
        shotid = int(str(shotrow["date"]) + str(shotrow["obs"]).zfill(3))
        
        try:
            args.log.info("Ingesting %s" % datevshot)
            file_obs = tb.open_file(op.join(args.shotdir, datevshot + ".h5"), "r")
            tableFibers_i = file_obs.root.Data.FiberIndex

            for row_i in tableFibers_i:

                row_main = tableFibers.row

                for col in tableFibers_i.colnames:
                    row_main[col] = row_i[col]

                fiberid = row_i["fiber_id"]

                try:
                    row_main["healpix"] = hp.ang2pix(
                        Nside, row_i["ra"], row_i["dec"], lonlat=True)
                except:
                    row_main["healpix"] = 0

                row_main["shotid"] = int(fiberid[0:10])
                row_main["specid"] = fiberid[20:23]
                row_main["ifuslot"] = fiberid[24:27]
                row_main["ifuid"] = fiberid[28:31]
                row_main["amp"] = fiberid[32:34]
                row_main.append()

            file_obs.close()

        except:            
            if shotid in badshot:
                pass
            else:
                args.log.error("could not ingest %s" % datevshot)

    tableFibers.cols.healpix.create_csindex()
    tableFibers.cols.ra.create_csindex()
    tableFibers.cols.fiber_id.create_csindex()
    tableFibers.flush()
    fileh.close()


if __name__ == "__main__":
    main()
