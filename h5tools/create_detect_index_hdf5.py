# -*- coding: utf-8 -*-
"""
Created: 2020/03/29

@author: Erin Mentuch Cooper

Combine Script to gather information for the FiberIndex Table. This provides
a master look up table for of all fibers contained in a 
HETDEX data release. Use the FiberIndex Class API provcd ..ided in 
hetdex_api/survey.py to easily access data and query 
astropy coordinates.

To run:

python3 create_detect_index_hdf5.py -of detect_index_hdr5.h5

"""

import numpy as np
from numpy.lib.recfunctions import append_fields
import tables as tb
import argparse as ap

import healpy as hp

from astropy.table import Table

from hetdex_api.input_utils import setup_logging


class Detections(tb.IsDescription):
    detectid = tb.Int64Col(pos=0)
    shotid = tb.Int64Col(pos=2)
    ra = tb.Float32Col(pos=3)
    dec = tb.Float32Col(pos=4)
    wave = tb.Float32Col(pos=7)
    sn = tb.Float32Col(pos=15)
    healpix = tb.Int64Col()
    det_type = tb.StringCol((4))
    survey = tb.StringCol((4))
    fiber_id = tb.StringCol((38))


def main(argv=None):
    """ Main Function """
    # Call initial parser from init_utils
    parser = ap.ArgumentParser(
        description="""Create HDF5 Astrometry file.""", add_help=True
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

    args = parser.parse_args(argv)
    args.log = setup_logging()

    fileh = tb.open_file(
        args.outfilename, mode="w", title=args.survey.upper() + " Detect Index file "
    )

    # tableDetects = fileh.create_table(
    #    fileh.root,
    #    "DetectIndex",
    #    Detections,
    #    expectedrows=40000000
    # )

    files = [
        "/scratch/projects/hetdex/hdr3/detect/detect_hdr3.h5",
        "/scratch/projects/hetdex/hdr3/detect/continuum_sources.h5",
        "/scratch/projects/hetdex/hdr4/detect/detect_hdr4.h5",
        "/scratch/projects/hetdex/hdr4/detect/continuum_sources.h5",
        "/scratch/projects/hetdex/hdr5/detect/detect_hdr5.h5",
        "/scratch/projects/hetdex/hdr5/detect/continuum_sources.h5",
    ]

    # set up HEALPIX options
    Nside = 2 ** 15
    firstfile = True

    for file in files:
        args.log.info("Appending detect H5 file: %s" % file)
        fileh_i = tb.open_file(file, "r")
        tableDetects_i = Table(fileh_i.root.Detections.read())
        tableDetects_i["healpix"] = hp.ang2pix(
            Nside, tableDetects_i["ra"], tableDetects_i["dec"], lonlat=True
        )

        if "continuum" in file:
            tableDetects_i["det_type"] = "cont".encode()
        else:
            tableDetects_i["det_type"] = "line".encode()

        if "hdr3" in file:
            tableDetects_i["survey"] = "hdr3".encode()
        elif "hdr4" in file:
            tableDetects_i["survey"] = "hdr4".encode()
        elif "hdr5" in file:
            tableDetects_i["survey"] = "hdr5".encode()
            
        if firstfile:

            tableDetects = fileh.create_table(
                fileh.root,
                "DetectIndex",
                tableDetects_i[list(Detections.columns.keys())].as_array(),
            )
            firstfile = False
        else:
            tableDetects.append(
                tableDetects_i[list(Detections.columns.keys())].as_array()
            )

        fileh_i.close()

    tableDetects.flush()
    tableDetects.cols.detectid.create_csindex()
    tableDetects.cols.healpix.create_csindex()
    tableDetects.cols.ra.create_csindex()
    tableDetects.cols.shotid.create_csindex()
    tableDetects.flush()
    fileh.close()
    args.log.info("Completed {}".format(args.outfilename))


if __name__ == "__main__":
    main()
