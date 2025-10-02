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
import sys
import glob
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
    datevobs = tb.StringCol((12))
    date = tb.Int64Col()
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

#same as above but directly adds in the 8 flags
#odering probably not necessary, but force here to match the ordering of the survey level mask file
class VIRUSFiberIndexWithFlags(tb.IsDescription):
    multiframe = tb.StringCol((20), pos=0)
    ra = tb.Float32Col(pos=1)
    dec = tb.Float32Col(pos=2)
    fiber_id = tb.StringCol((38), pos=3)
    healpix = tb.Int64Col(pos=4)
    amp = tb.StringCol(2,pos=5)
    date = tb.Int64Col(pos=6)
    datevobs = tb.StringCol((12),pos=7)
    expnum = tb.Int32Col(pos=8)
    fibidx = tb.Int32Col(pos=9)
    fibnum = tb.Int32Col(pos=10)
    fpx = tb.Float32Col(pos=11)
    fpy = tb.Float32Col(pos=12)
    ifuslot = tb.StringCol(3,pos=13)
    ifuid = tb.StringCol(3,pos=14)
    ifux = tb.Float32Col(pos=15)
    ifuy = tb.Float32Col(pos=16)
    shotid = tb.Int64Col(pos=17)
    specid = tb.StringCol(3,pos=18)
    flag = tb.Int32Col(dflt=1,pos=19) #1 is "good"
    flag_badamp = tb.Int8Col(dflt=1,pos=20) #1 is "good"
    flag_badfib = tb.Int8Col(dflt=1,pos=21) #1 is "good"
    flag_meteor = tb.Int8Col(dflt=1,pos=22) #1 is "good"
    flag_satellite = tb.Int8Col(dflt=1,pos=23) #1 is "good"
    flag_largegal = tb.Int8Col(dflt=1,pos=24) #1 is "good"
    flag_shot = tb.Int8Col(dflt=1,pos=25) #1 is "good"
    flag_throughput = tb.Int8Col(dflt=1,pos=26) #1 is "good"



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
        default="/scratch/projects/hetdex/hdr4/reduction/data",
    )

    parser.add_argument(
        "-sl",
        "--shotlist",
        help="""Text file of DATE OBS list""",
        type=str,
        default="hdr4.shotlist",
    )

    parser.add_argument(
        "-of",
        "--outfilename",
        type=str,
        help="""Relative or absolute path for output HDF5
                        file.""",
        default=None,
    )

    parser.add_argument("-survey", "--survey", type=str, default="hdr4")

    parser.add_argument(
        "-m",
        "--month",
        type=int,
        default=None,
        help="""Create FiberIndex for a single month""",
    )

    parser.add_argument(
        "--merge",
        "-merge",
        help="""Boolean trigger to merge all 2*.fits files in cwd""",
        default=False,
        required=False,
        action="store_true",
    )

    parser.add_argument(
        "-shot_h5",
        "--shot_h5",
        help="""A specific shot HDF5 file to update. If provided, only the shot HDF5 file is updated. There
        is no change to the survey and no separate fiber index file is created.""",
        type=str,
        default=None
    )

    args = parser.parse_args(argv)
    args.log = setup_logging()

    if args.shot_h5 is None: #normal, HETDEX historical mode
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
            expectedrows=300000000,
        )

    else: #single shot h5
        args.merge = False
        args.month = None
        fileh = tb.open_file(args.shot_h5,mode="a")
        shotlist = [fileh.root.Shot.read(field="shotid")[0]]

        if fileh.__contains__("/FiberIndex"):
            tableFibers = fileh.root.FiberIndex
        else:
            tableFibers = fileh.create_table(
                fileh.root,
                "FiberIndex",
                VIRUSFiberIndexWithFlags,
                "Shot Fiber Coord Info",
                expectedrows=105000,
            )


    if args.merge:
        files = glob.glob("mfi*h5")
        files.sort()
        
        for file in files:
            args.log.info("Appending detect H5 file: %s" % file)
            fileh_i = tb.open_file(file, "r")
            tableFibers_i = fileh_i.root.FiberIndex.read()
            tableFibers.append(tableFibers_i)

        tableFibers.flush()
        tableFibers.cols.healpix.create_csindex()
        tableFibers.cols.ra.create_csindex()
        tableFibers.cols.shotid.create_csindex()
        tableFibers.flush()
        fileh.close()
        args.log.info("Completed {}".format(args.outfilename))
        sys.exit()

    # set up HEALPIX options
    Nside = 2 ** 15
    hp.max_pixrad(Nside, degrees=True) * 3600  # in unit of arcsec

    config = HDRconfig(survey=args.survey)

    badshot = np.loadtxt(config.badshot, dtype=int)

    if args.month is not None:
        args.log.info("Working on month {}".format(args.month))
        # if working on a single month downselect
        shotlist["month"] = np.array(shotlist["date"] / 100, dtype=int)
        sel_month = shotlist["month"] == args.month
        shotlist = shotlist[sel_month]

    for shotrow in shotlist:
        try:
            datevshot = str(shotrow["date"]) + "v" + str(shotrow["obs"]).zfill(3)
            shotid = int(str(shotrow["date"]) + str(shotrow["obs"]).zfill(3))

            date = shotrow["date"]
        except: #this is the single shot version with shotids
            datevshot = str(shotrow)[0:8] + "v" + str(shotrow)[-3:]
            shotid = int(shotrow)
            date = str(shotrow)[0:8]

        try:
            args.log.info("Ingesting %s" % datevshot)
            if args.shot_h5 is None:
                file_obs = tb.open_file(op.join(args.shotdir, datevshot + ".h5"), "r")
            else:
                file_obs = fileh #tb.open_file(args.shot_h5, mode="r")

            tableFibers_i = file_obs.root.Data.FiberIndex

            for row_i in tableFibers_i:

                row_main = tableFibers.row
                
                for col in tableFibers_i.colnames:
                    row_main[col] = row_i[col]

                fiberid = row_i["fiber_id"]

                try:
                    row_main["healpix"] = hp.ang2pix(
                        Nside, row_i["ra"], row_i["dec"], lonlat=True
                    )
                except:
                    row_main["healpix"] = 0

                row_main["shotid"] = shotid
                row_main["date"] = date
                row_main["datevobs"] = datevshot

                row_main["specid"] = fiberid[20:23]
                row_main["ifuslot"] = fiberid[24:27]
                row_main["ifuid"] = fiberid[28:31]
                row_main["amp"] = fiberid[32:34]
                row_main.append()

            if file_obs != fileh:
                file_obs.close()

        except:
            if shotid in badshot:
                pass
            else:
                args.log.error("could not ingest %s" % datevshot)

    tableFibers.flush()

    try: #this could be an update, in which case, remove the old index before re-creating
        tableFibers.cols.healpix.remove_index()
        tableFibers.cols.ra.remove_index()
        tableFibers.cols.shotid.remove_index()
        tableFibers.flush()
    except:
        pass

    tableFibers.cols.healpix.create_csindex()
    tableFibers.cols.ra.create_csindex()
    tableFibers.cols.shotid.create_csindex()
    tableFibers.flush()
    fileh.close()
    if args.shot_h5 is None:
        args.log.info("Completed {}".format(args.outfilename))
    else:
        args.log.info("Completed {}".format(args.shot_h5))


if __name__ == "__main__":
    main()
