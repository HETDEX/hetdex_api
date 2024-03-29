# -*- coding: utf-8 -*-
"""
Created: 2020/07/15

@author: Erin Mentuch Cooper

This file contains all information related to the HETDEX continuum source 
catalog

Examples
--------

To run continuum sources:

>>> python3 create_cont_hdf5.py -dp /data/00115/gebhardt/cs/spec /
-cs /data/00115/gebhardt/cs/rext2 -of continuum_sources.h5

"""
from __future__ import print_function

import sys
import os
import os.path as op
import argparse as ap
import re
import glob
import subprocess
import numpy as np
import tables as tb

from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table, vstack
import astropy.units as u
from hetdex_api.input_utils import setup_logging

import warnings

warnings.filterwarnings("ignore")


def get_detectname(ra, dec):
    """
    convert ra,dec coordinates to a IAU-style object name.
    """
    coord = SkyCoord(ra * u.deg, dec * u.deg)

    return "HETDEX J{0}{1}".format(
        coord.ra.to_string(unit=u.hourangle, sep="", precision=2, pad=True),
        coord.dec.to_string(sep="", precision=1, alwayssign=True, pad=True),
    )


class Detections(tb.IsDescription):
    shotid = tb.Int64Col(pos=2)
    date = tb.Int32Col(pos=5)
    obsid = tb.Int32Col(pos=6)
    detectid = tb.Int64Col(pos=0)
    fiber_id = tb.StringCol((38))
    detectname = tb.StringCol((40))
    ra = tb.Float32Col(pos=3)
    dec = tb.Float32Col(pos=4)
    wave = tb.Float32Col(pos=7)
    wave_err = tb.Float32Col(pos=8)
    flux = tb.Float32Col(pos=9)
    flux_err = tb.Float32Col(pos=10)
    linewidth = tb.Float32Col(pos=11)
    linewidth_err = tb.Float32Col(pos=12)
    continuum = tb.Float32Col(pos=13)
    continuum_err = tb.Float32Col(pos=14)
    sn = tb.Float32Col(pos=15)
    sn_err = tb.Float32Col(pos=16)
    chi2 = tb.Float32Col(pos=17)
    chi2_err = tb.Float32Col(pos=18)
    noise_ratio = tb.Float32Col()
    inputid = tb.StringCol((40))
    x_raw = tb.Int32Col(pos=21)
    y_raw = tb.Int32Col(pos=22)
    x_ifu = tb.Float32Col()
    y_ifu = tb.Float32Col()
    weight = tb.Float32Col()
    fibnum = tb.Int32Col(pos=20)
    multiframe = tb.StringCol((20), pos=19)
    specid = tb.StringCol((3))
    ifuslot = tb.StringCol((3))
    ifuid = tb.StringCol((3))
    amp = tb.StringCol((2))
    expnum = tb.Int32Col()
    chi2fib = tb.Float32Col()


class Spectra(tb.IsDescription):
    detectid = tb.Int64Col(pos=0)
    wave1d = tb.Float32Col(1036, pos=1)
    spec1d = tb.Float32Col(1036, pos=2)
    spec1d_err = tb.Float32Col(1036, pos=3)
    counts1d = tb.Float32Col(1036, pos=4)
    counts1d_err = tb.Float32Col(1036, pos=5)
    apsum_counts = tb.Float32Col(1036, pos=6)
    apsum_counts_err = tb.Float32Col(1036, pos=7)
    spec1d_ffsky = tb.Float32Col(1036)
    spec1d_nc = tb.Float32Col(1036)
    spec1d_nc_err = tb.Float32Col(1036)
    apcor = tb.Float32Col(1036)
    flag_pix = tb.Float32Col(1036)


class Fibers(tb.IsDescription):
    detectid = tb.Int64Col(pos=0)
    ra = tb.Float32Col(pos=1)
    dec = tb.Float32Col(pos=2)
    x_ifu = tb.Float32Col(pos=5)
    y_ifu = tb.Float32Col(pos=6)
    multiframe = tb.StringCol((20), pos=3)
    fibnum = tb.Int32Col()
    fiber_id = tb.StringCol((38), pos=4)
    expnum = tb.Int32Col(pos=9)
    distance = tb.Float32Col(pos=10)
    wavein = tb.Float32Col(pos=12)
    timestamp = tb.StringCol((17), pos=11)
    date = tb.Int32Col(pos=7)
    obsid = tb.Int32Col(pos=8)
    flag = tb.Int32Col(pos=13)
    weight = tb.Float32Col(pos=14)
    ADC = tb.Float32Col((5), pos=15)
    specid = tb.StringCol((3))
    ifuslot = tb.StringCol((3))
    ifuid = tb.StringCol((3))
    amp = tb.StringCol((2))
    x_raw = tb.Int32Col()
    y_raw = tb.Int32Col()


def main(argv=None):
    """ Main Function """
    # Call initial parser from init_utils
    parser = ap.ArgumentParser(description="""Create HDF5 file.""", add_help=True)

    parser.add_argument(
        "-cs",
        "--contsource",
        help="""Path to Karls rext catalog""",
        type=str,
        default=None,
    )

    parser.add_argument(
        "-dp",
        "--detect_path",
        help="""Path to detections""",
        type=str,
        default="/data/00115/gebhardt/alldet/output",
    )

    parser.add_argument(
        "-of",
        "--outfilename",
        type=str,
        help="""Relative or absolute path for output HDF5
        file.""",
        default=None,
    )

    parser.add_argument(
        "-a",
        "--append",
        help="""Appending to existing detection HDF5 file.""",
        action="count",
        default=0,
    )

    parser.add_argument(
        "-sl",
        "--shotlist",
        help="""Text file of DATE OBS list""",
        type=str,
        default="/scratch/03946/hetdex/hdr3/survey/hdr3.shotlist",
    )
   
    args = parser.parse_args(argv)
    args.log = setup_logging()

    outfilename = args.outfilename

    fileh = tb.open_file(outfilename, "w", "HDR3 Continuum Source Database")
    index_buff = 3090000000
    detectidx = index_buff

    detectcat = Table.read(args.contsource, format="ascii.no_header")
    detectcat.remove_columns(
        [
            "col1",
            "col4",
            "col5",
            "col6",
            "col9",
            "col10",
            "col11",
            "col12",
            "col13",
            "col14",
        ]
    )
    detectcat["col2"].name = "ra"
    detectcat["col3"].name = "dec"
    detectcat["col7"].name = "obnum"
    detectcat["col8"].name = "datevshot"

    tableMain = fileh.create_table(
        fileh.root,
        "Detections",
        Detections,
        "HETDEX Continuum Source Catalog",
        expectedrows=np.size(detectcat),
    )
    tableFibers = fileh.create_table(
        fileh.root,
        "Fibers",
        Fibers,
        "Fiber info for each source",
        expectedrows=np.size(detectcat),
    )
    tableSpectra = fileh.create_table(
        fileh.root,
        "Spectra",
        Spectra,
        "1D Spectra for each Line Detection",
        expectedrows=15 * np.size(detectcat),
    )

    shotid = []
    date = []
    obsid = []
    inputid = []
    detectid = []

    detectid_i = detectidx

    for row in detectcat:
        p = re.compile("v")
        shotid_i = int(p.sub("", row["datevshot"]))
        inputid_i = str(row["datevshot"]) + "_" + str(row["obnum"])

        detectid.append(detectid_i)
        inputid.append(inputid_i)
        date.append(int(str(shotid_i)[0:8]))
        obsid.append(int(str(shotid_i)[8:11]))
        shotid.append(shotid_i)
        detectid_i += 1

    detectcat["detectid"] = detectid
    detectcat["inputid"] = inputid
    detectcat["date"] = date
    detectcat["obsid"] = obsid
    detectcat["detectid"] = detectid
    detectcat["shotid"] = shotid

    det_cols = fileh.root.Detections.colnames

    shottab = Table.read(args.shotlist, format='ascii.no_header')
    shotlist = []
    for row in shottab:
        shotlist.append( int(str(row['col1']) + str(row['col2']).zfill(3)))

    for row in detectcat:

        if row['shotid'] not in shotlist:
            continue

        rowMain = tableMain.row

        for col in det_cols:
            try:
                rowMain[col] = row[col]
            except:
                rowMain[col] = 0.0

        rowMain.append()

    tableMain.flush()

    # add spectra for each detectid in the detections table

    for row in tableMain:
        try:

            inputid_i = row["inputid"].decode()
            specfile = op.join(args.detect_path, inputid_i + ".spec")
            dataspec = Table(
                np.loadtxt(specfile),
                names=[
                    "wave1d",
                    "spec1d_nc",
                    "spec1d_nc_err",
                    "counts1d",
                    "counts1d_err",
                    "apsum_counts",
                    "apsum_counts_err",
                    "dummy",
                    "apcor",
                    "flag_pix",
                    "obnum",
                    "spec1d_nc_ffsky",
                ],
            )

            rowspectra = tableSpectra.row
            rowspectra["detectid"] = row["detectid"]
            rowspectra["spec1d"] = dataspec["spec1d_nc"] / dataspec["apcor"]
            rowspectra["spec1d_err"] = dataspec["spec1d_nc_err"] / dataspec["apcor"]
            rowspectra["spec1d_ffsky"] = dataspec["spec1d_nc_ffsky"] / dataspec["apcor"]
            rowspectra["wave1d"] = dataspec["wave1d"]
            rowspectra["spec1d_nc"] = dataspec["spec1d_nc"]
            rowspectra["spec1d_nc_err"] = dataspec["spec1d_nc_err"]
            rowspectra["counts1d"] = dataspec["counts1d"]
            rowspectra["counts1d_err"] = dataspec["counts1d_err"]
            rowspectra["apsum_counts"] = dataspec["apsum_counts"]
            rowspectra["apsum_counts_err"] = dataspec["apsum_counts_err"]
            rowspectra["apcor"] = dataspec["apcor"]
            rowspectra["flag_pix"] = dataspec["flag_pix"]
            rowspectra.append()
        except Exception:
            args.log.error("Could not ingest %s" % specfile)

    tableSpectra.flush()

    # add fiber info for each detection

    for row in tableMain:
        inputid_i = row["inputid"].decode()

        filefiberinfo = op.join(args.detect_path, inputid_i + ".list")

        try:
            datafiber = Table.read(filefiberinfo, format="ascii.no_header")

            for ifiber in np.arange(np.size(datafiber)):
                rowfiber = tableFibers.row
                rowfiber["detectid"] = row["detectid"]
                rowfiber["ra"] = datafiber["col1"][ifiber]
                rowfiber["dec"] = datafiber["col2"][ifiber]
                rowfiber["x_ifu"] = datafiber["col3"][ifiber]
                rowfiber["y_ifu"] = datafiber["col4"][ifiber]
                rowfiber["expnum"] = str(datafiber["col6"][ifiber])[3:5]
                multiname = datafiber["col5"][ifiber]
                multiframe = multiname[0:20]
                fiber_id_i = (
                    str(row["shotid"])
                    + "_"
                    + str(int(rowfiber["expnum"]))
                    + "_"
                    + multiframe
                    + "_"
                    + str(int(multiname[21:24])).zfill(3)
                )
                rowfiber["fiber_id"] = fiber_id_i
                rowfiber["multiframe"] = multiframe
                rowfiber["specid"] = multiframe[6:9]
                rowfiber["ifuslot"] = multiframe[10:13]
                rowfiber["ifuid"] = multiframe[14:17]
                rowfiber["amp"] = multiframe[18:20]
                rowfiber["fibnum"] = int(multiname[21:24])
                rowfiber["distance"] = datafiber["col7"][ifiber]
                rowfiber["wavein"] = datafiber["col8"][ifiber]
                rowfiber["timestamp"] = datafiber["col9"][ifiber]
                rowfiber["date"] = datafiber["col10"][ifiber]
                rowfiber["obsid"] = str(datafiber["col11"][ifiber])[0:3]
                rowfiber["x_raw"] = datafiber["col12"][ifiber]
                rowfiber["y_raw"] = datafiber["col13"][ifiber]
                rowfiber["flag"] = datafiber["col15"][ifiber]
                rowfiber["weight"] = datafiber["col14"][ifiber]
                rowfiber.append()

            # Now append brightest fiber info to Detections table:
            ifiber = 0  # np.argmax(datafiber["col14"])
            multiname = datafiber["col5"][ifiber]
            multiframe = multiname[0:20]
            row["expnum"] = int(str(datafiber["col6"][ifiber])[3:5])
            fiber_id_i = (
                str(row["shotid"])
                + "_"
                + str(row["expnum"])
                + "_"
                + multiframe
                + "_"
                + str(int(multiname[21:24])).zfill(3)
            )
            row["fiber_id"] = fiber_id_i
            row["multiframe"] = multiframe
            row["specid"] = multiframe[6:9]
            row["ifuslot"] = multiframe[10:13]
            row["ifuid"] = multiframe[14:17]
            row["amp"] = multiframe[18:20]
            row["fibnum"] = int(multiname[21:24])
            row["x_raw"] = datafiber["col12"][ifiber]
            row["y_raw"] = datafiber["col13"][ifiber]
            row["x_ifu"] = datafiber["col3"][ifiber]
            row["y_ifu"] = datafiber["col4"][ifiber]
            row["expnum"] = str(datafiber["col6"][ifiber])[3:5]
            row["weight"] = datafiber["col14"][ifiber]
            row.update()

        except Exception:
            args.log.error("Could not ingest %s" % filefiberinfo)

        tableFibers.flush()

    tableMain.cols.detectid.create_csindex()
    tableFibers.cols.detectid.create_csindex()
    tableSpectra.cols.detectid.create_csindex()
    tableFibers.flush()  # just to be safe
    tableSpectra.flush()
    tableMain.flush()
    args.log.info("File finished: %s" % args.outfilename)
    fileh.close()


if __name__ == "__main__":
    main()
