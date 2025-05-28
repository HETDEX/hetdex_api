# -*- coding: utf-8 -*-
"""
Created: 2020/07/15

@author: Erin Mentuch Cooper

This file contains all information related to the HETDEX continuum source 
catalog

Examples
--------

To run continuum sources one month at a time:

#>>> python3 create_cont_hdf5.py -d 20200523 -o 012

To merge into one month

>> python3 create_cont_hdf5.py --merge -m 201901 

To merge into one file:

>> python3 create_cont_hdf5py --merge -of continuum_sources.h5  
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
import tarfile

from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table, vstack
import astropy.units as u
from hetdex_api.input_utils import setup_logging

import warnings
import traceback

warnings.filterwarnings("ignore")


DETECTID_BASE = int(5.09e9) #note the x.09 for continuum sources
DETECTID_STEP = 1
DETECTID_OFFSET = 0        #last detectid (w/o the BASE) from previous release -OR- 0 if starting a new release
                           #start numbering DETECTID_STEP past this value (usually +1)
                           #297876 for HDR4 extension to HDR3


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
        "-survey", "--survey", help="""{hdr1, hdr2, hdr2.1, hdr3, hdr4, hdr5, pdr1}""",
        type=str,
        default="hdr4"
    )

    parser.add_argument(
        "-m", "--month", help="""Month to run: 201901""", type=str, default=None
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
        "-dp",
        "--detect_path",
        help="""Path to detections""",
        type=str,
        default="/scratch/00115/gebhardt/cs/spec",
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

    #parser.add_argument(
    #    "-sl",
    #    "--shotlist",
    #    help="""Text file of DATE OBS list""",
    #    type=str,
    #    default="/scratch/projects/hetdex/hdr4/survey/hdr4.shotlist",
    #)

    parser.add_argument(
        "--merge",
        "-merge",
        help="""Boolean trigger to merge all cont_2*.fits files in cwd""",
        default=False,
        required=False,
        action="store_true",
    )

    parser.add_argument(
        "-md",
        "--mergedir",
        help="""Merge all HDF5 files in the defined merge
        directory. Can append to existing file using --append option""",
        type=str,
        default=os.getcwd(),
    )

    parser.add_argument(
        "--reindex",
        "-reindex",
        help="""Boolean to reindex an append""",
        default=False,
        required=False,
        action="store_true",
    )
    
    args = parser.parse_args(argv)
    args.log = setup_logging()

    index_buff = DETECTID_BASE + DETECTID_OFFSET + DETECTID_STEP # starting count post HDR3 detection indices
    detectidx = index_buff
    bad_ingest_files = []
    
    if args.merge:
        n_size = 300000
        
        fileh = tb.open_file(args.outfilename, "w", "{} Continuum Source Database".format(args.survey.upper()))

        tableMain = fileh.create_table(
            fileh.root,
            "Detections",
            Detections,
            "HETDEX Continuum Source Catalog",
            expectedrows= n_size,
        )
        tableFibers = fileh.create_table(
            fileh.root,
            "Fibers",
            Fibers,
            "Fiber info for each source",
            expectedrows= n_size,
        )
        tableSpectra = fileh.create_table(
            fileh.root,
            "Spectra",
            Spectra,
            "1D Spectra for each Line Detection",
            expectedrows=15 * n_size,
        )
        files = sorted(glob.glob(op.join(args.mergedir, "cont_2*.h5")))

        detectid_max = 0



        for file in files:
            
            args.log.info("Appending detect H5 file: %s" % file)

            try:
                fileh_i = tb.open_file(file, "r")
                tableMain_i = fileh_i.root.Detections.read()
    
                if np.size(tableMain_i) == 0:
                    args.log.error('No detections for %s' % file)
                    continue
                
                tableFibers_i = fileh_i.root.Fibers.read()
                tableSpectra_i = fileh_i.root.Spectra.read()

            except Exception as e:
                args.log.error('Could not append {}'.format(file))
                args.log.error(f"Exception: {e}\n\n{traceback.format_exc()}")
                bad_ingest_files.append(file)
                continue
                
            tableMain_i["detectid"] += detectid_max
            tableFibers_i["detectid"] += detectid_max
            tableSpectra_i["detectid"] += detectid_max
            
            # after first table be sure to add one to the index
            
            detectid_max = 1
            
            tableMain.append(tableMain_i)
            tableFibers.append(tableFibers_i)
            tableSpectra.append(tableSpectra_i)
            
            detectid_max = np.max(tableMain.cols.detectid[:]) - index_buff + DETECTID_STEP

            fileh_i.close()
            tableFibers.flush()  # just to be safe
            tableSpectra.flush()
            tableMain.flush()

        tableMain.cols.shotid.create_csindex()
        tableMain.cols.detectid.create_csindex()
        tableFibers.cols.detectid.create_csindex()
        tableSpectra.cols.detectid.create_csindex()
        tableFibers.flush()  # just to be safe
        tableSpectra.flush()
        tableMain.flush()
        args.log.info("File finished: %s" % args.outfilename)

        try:
            if len(bad_ingest_files) > 0:
                args.log.error("WARNING!!!! Some files failed to merge. See bad_merge_cont.rerun")
                with open("bad_merge_cont.rerun", "w") as f:
                    for bif in bad_ingest_files:
                        bn = os.path.basename(bif).split("_")[1].split(".")[0]
                        date = bn.split("v")[0]
                        obs = bn.split("v")[1]
                        f.write(f"ingest_cont {date} {obs}\n")
            else:
                args.log.info("All files merged.")

        except Exception as e:
            args.log.error(f"Exception: {e}\n\n{traceback.format_exc()}")

        sys.exit()
    # open up datevobs tarball with ingestion data

    datevobs = str(args.date) + 'v' + str(args.observation).zfill(3)
    
    spectarfile = op.join(args.detect_path, "{}cs.tar".format(datevobs))
    
    if not op.exists(spectarfile):
        args.log.error("Could not find {}".format(spectarfile))
        sys.exit()
    
    if args.outfilename:
        outfilename = args.outfilename
    elif args.month and args.merge:
        outfilename = 'cont_month_' + str(args.month) + '.h5'
    else:
        outfilename = 'cont_'+ str(args.date) + 'v' + str(args.observation).zfill(3) + '.h5'

#    contfiles = np.loadtxt('complete_cont_h5list', dtype=str)
#    if outfilename in contfiles:
#        sys.exit('{} already complete'.format(outfilename))

    # open up datevobs tarball with ingestion data      
    spectar = tarfile.open(spectarfile)

    if not os.path.exists('./spec'):
        os.makedirs('./spec')
       
    spectar.extractall('./spec')
    spectar.close()

    try:
        detectcat = Table.read('./spec/{}.rcs'.format(datevobs), format="ascii.no_header")
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
    except:
        args.log.error("Could not read {}.rcs".format(datevobs))
        sys.exit()
        
    detectcat["col2"].name = "ra"
    detectcat["col3"].name = "dec"
    detectcat["col7"].name = "obnum"
    detectcat["col8"].name = "datevshot"

    if args.append:
        fileh = tb.open_file(outfilename, "a", "{} Continuum Source Database".format(args.survey.upper()))
        tableMain = fileh.root.Detections
        tableSpectra = fileh.root.Spectra
        tableFibers = fileh.root.Fibers

    else:
        fileh = tb.open_file(outfilename, "w", "{} Continuum Source Database".format(args.survey.upper()))
        
        tableMain = fileh.create_table(
            fileh.root,
            "Detections",
            Detections,
            "HETDEX Continuum Source Catalog",
            expectedrows= np.size(detectcat),
        )
        tableFibers = fileh.create_table(
            fileh.root,
            "Fibers",
            Fibers,
            "Fiber info for each source",
            expectedrows= np.size(detectcat),
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

    if args.append:
        detectid_i = np.max(tableMain.cols.detectid[:]) + DETECTID_STEP
    else:
        detectid_i = detectidx
    print(detectid_i)
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

    #Not needed as we are generating for an input shot now. Removed 2023-09-22 EMC
    #shottab = Table.read(args.shotlist, format='ascii.no_header')
    #shotlist = []
    #for row in shottab:
    #    shotlist.append( int(str(row['col1']) + str(row['col2']).zfill(3)))

    for row in detectcat:

        #if row['shotid'] not in shotlist:
        #    continue

        rowMain = tableMain.row

        for col in det_cols:
            try:
                rowMain[col] = row[col]
            except:
                rowMain[col] = 0.0

        try:
            inputid_i = row["inputid"]
            specfile = "./spec/{}.spec".format(inputid_i)

            # 20240416 dd - add some flexibility on whether .spec file has 11 or 12 columns
            # check column count
            expected_names = [
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
            ]

            try:
                ncols = -1
                with open(specfile) as f:
                    ncols = len(f.readline().split())

                if ncols == 12:
                    expected_names.append("spec1d_nc_ffsky")
                elif ncols == 11:
                    pass  # all good as is
                else:
                    args.log.warning(f"Unexpected number of columns ({ncols}) in {specfile}")

            except Exception as e:
                args.log.warning('Possible problem with ' + specfile)
                args.log.warning(f"Exception: {e}\n\n{traceback.format_exc()}")


            dataspec = Table(
                np.loadtxt(specfile),
                names=expected_names
            )

            rowspectra = tableSpectra.row
            rowspectra["detectid"] = row["detectid"]
            rowspectra["spec1d"] = dataspec["spec1d_nc"] / dataspec["apcor"]
            rowspectra["spec1d_err"] = dataspec["spec1d_nc_err"] / dataspec["apcor"]
            if "spec1d_nc_ffsky" in dataspec.colnames:
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
        except Exception as e:
            args.log.error("Could not ingest %s" % specfile)
            args.log.error(f"Exception: {e}\n\n{traceback.format_exc()}")


        # add fiber data
        inputid_i = row["inputid"]
        filefiberinfo = "./spec/{}.list".format(inputid_i)

        if True:#try:
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
            ifiber = np.argmax(datafiber["col14"])
            multiname = datafiber["col5"][ifiber]
            multiframe = multiname[0:20]
            rowMain["expnum"] = int(str(datafiber["col6"][ifiber])[3:5])
            fiber_id_i = (
                str(rowMain["shotid"])
                + "_"
                + str(rowMain["expnum"])
                + "_"
                + multiframe
                + "_"
                + str(int(multiname[21:24])).zfill(3)
            )
            rowMain["fiber_id"] = fiber_id_i
            rowMain["multiframe"] = multiframe
            rowMain["specid"] = multiframe[6:9]
            rowMain["ifuslot"] = multiframe[10:13]
            rowMain["ifuid"] = multiframe[14:17]
            rowMain["amp"] = multiframe[18:20]
            rowMain["fibnum"] = int(multiname[21:24])
            rowMain["x_raw"] = datafiber["col12"][ifiber]
            rowMain["y_raw"] = datafiber["col13"][ifiber]
            rowMain["x_ifu"] = datafiber["col3"][ifiber]
            rowMain["y_ifu"] = datafiber["col4"][ifiber]
            rowMain["expnum"] = str(datafiber["col6"][ifiber])[3:5]
            rowMain["weight"] = datafiber["col14"][ifiber]
            rowMain.append()
        else:#except Exception:
            args.log.error("Could not ingest %s" % filefiberinfo)

        tableMain.flush()
        tableSpectra.flush()
        tableFibers.flush()

    if args.append:
        if args.reindex:
            args.log.info("Reindexing the detectid column")
            tableMain.cols.detectid.reindex()
            tableFibers.cols.detectid.reindex()
            tableSpectra.cols.detectid.reindex()
        tableFibers.flush()  # just to be safe
        tableSpectra.flush()
        tableMain.flush()
    else:
        args.log.info("Indexing the detectid column")
        tableMain.cols.detectid.create_csindex()
        tableFibers.cols.detectid.create_csindex()
        tableSpectra.cols.detectid.create_csindex()
        tableFibers.flush()  # just to be safe
        tableSpectra.flush()
        tableMain.flush()
        
    args.log.info("File finished: {}".format(outfilename))
    fileh.close()

    # remove untarred files
    # this was very slow so removed for HDR4. Just remove manually after

    #tarfiles = glob.glob('./spec/{}*'.format(datevobs))
    #
    #for f in tarfiles:
    #    try:
    #        os.remove(f)
    #    except OSError as e:
    #        print("Error: %s : %s" % (f, e.strerror))
                         
if __name__ == "__main__":
    main()
