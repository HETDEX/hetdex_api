# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 09:52:35 2017

@author: gregz
"""
import glob
import re
import os
import shutil
import tarfile
import sys

import tables as tb
import argparse as ap
import os.path as op
import numpy as np

from astropy.io import fits
from hetdex_api.input_utils import setup_logging
from astropy.table import Table

# hard coded variable to initialize 'rms', 'chi' arrays
# and remove 'twi_spectrum' for all realeases past hdr1
global hdr_survey
hdr_survey = "hdr2"


def build_path(reduction_folder, instr, date, obsid, expn):
    folder = op.join(
        date,
        instr,
        "{:s}{:07d}".format(instr, int(obsid)),
        "exp{:02d}".format(int(expn)),
        instr,
    )
    return op.join(reduction_folder, folder)


def get_files(args):
    if args.tar == False:
        files = glob.glob(
            op.join(
                args.rootdir,
                args.date,
                "virus",
                "virus%07d" % int(args.observation),
                "exp*",
                "virus",
                "multi_*.fits",
            )
        )
    else:
        datestr = "d%ss%03d" % (args.date, int(args.observation))

        files = glob.glob(
            op.join(
                args.rootdir,
                "sci" + str(args.date)[0:6],
                datestr + "exp0?",
                datestr + "exp0?_mu.tar",
            )
        )

    return files


def define_field(objname):
    if re.match("par", str(objname)):
        field = "parallel"
    elif re.match("COS|cos|DEXcos", str(objname)):
        field = "cosmos"
    elif re.match("EGS", str(objname)):
        field = "egs"
    elif re.match("GN", str(objname)):
        field = "goods-n"
    elif re.match("DEX0|DEXfl|HF", str(objname)):
        field = "dex-fall"
    elif re.match("HS|DEXsp", str(objname)):
        field = "dex-spring"
    else:
        field = "other"

    return field


class VIRUSFiberIndex(tb.IsDescription):
    multiframe = tb.StringCol((20), pos=0)
    fiber_id = tb.StringCol((38), pos=4)
    fibidx = tb.Int32Col()
    fibnum = tb.Int32Col()
    ifux = tb.Float32Col()
    ifuy = tb.Float32Col()
    fpx = tb.Float32Col()
    fpy = tb.Float32Col()
    ra = tb.Float32Col(pos=1)
    dec = tb.Float32Col(pos=2)
    expnum = tb.Int32Col()


class VIRUSFiber(tb.IsDescription):
    obsind = tb.Int32Col()
    multiframe = tb.StringCol((20), pos=0)
    fiber_id = tb.StringCol((38), pos=4)
    fibidx = tb.Int32Col()
    fibnum = tb.Int32Col()
    ifux = tb.Float32Col()
    ifuy = tb.Float32Col()
    fpx = tb.Float32Col()
    fpy = tb.Float32Col()
    ra = tb.Float32Col(pos=1)
    dec = tb.Float32Col(pos=2)
    spectrum = tb.Float32Col((1032,))
    wavelength = tb.Float32Col((1032,))
    fiber_to_fiber = tb.Float32Col((1032,))
    global hdr_survey
    if hdr_survey == "hdr1":
        twi_spectrum = tb.Float32Col((1032,))
    else:
        chi2 = tb.Float32Col((1032,))
        rms = tb.Float32Col((1032,))
        calfib_counts = tb.Float32Col((1036,))
        calfibe_counts = tb.Float32Col((1036,))
    trace = tb.Float32Col((1032,))
    sky_subtracted = tb.Float32Col((1032,))
    sky_spectrum = tb.Float32Col((1032,))
    error1Dfib = tb.Float32Col((1032,))
    calfib = tb.Float32Col((1036,))
    calfibe = tb.Float32Col((1036,))
    Amp2Amp = tb.Float32Col((1036,))
    Throughput = tb.Float32Col((1036,))
    ifuslot = tb.StringCol(3)
    ifuid = tb.StringCol(3)
    specid = tb.StringCol(3)
    contid = tb.StringCol(8)
    amp = tb.StringCol(2)
    expnum = tb.Int32Col()


class VIRUSImage(tb.IsDescription):
    obsind = tb.Int32Col()
    multiframe = tb.StringCol((20), pos=0)
    image = tb.Float32Col((1032, 1032))
    error = tb.Float32Col((1032, 1032))
    clean_image = tb.Float32Col((1032, 1032))
    ifuslot = tb.StringCol(3)
    ifuid = tb.StringCol(3)
    specid = tb.StringCol(3)
    contid = tb.StringCol(8)
    amp = tb.StringCol(2)
    expnum = tb.Int32Col()


class VIRUSShot(tb.IsDescription):
    shotid = tb.Int64Col(pos=0)
    date = tb.Int32Col(pos=3)
    obsid = tb.Int32Col(pos=4)
    objid = tb.StringCol((18), pos=2)
    field = tb.StringCol((12), pos=1)
    ra = tb.Float64Col(pos=5)
    dec = tb.Float64Col(pos=6)
    pa = tb.Float64Col(pos=7)
    n_ifu = tb.Int32Col(pos=8)
    datevobs = tb.StringCol((12))
    trajcra = tb.Float32Col()
    trajcdec = tb.Float32Col()
    trajcpa = tb.Float32Col()
    structaz = tb.Float32Col()
    time = tb.StringCol(7)
    ambtemp = tb.Float32Col()
    humidity = tb.Float32Col()
    dewpoint = tb.Float32Col()
    pressure = tb.Float32Col()
    expnum = tb.Int32Col((3))
    exptime = tb.Float32Col((3))
    darktime = tb.Float32Col((3))
    mjd = tb.Float32Col((3))
    fwhm_virus = tb.Float32Col(pos=9)
    fwhm_virus_err = tb.Float32Col(pos=10)
    nstars_fit_fwhm = tb.Int32Col()
    # relflux_guider = tb.Float32Col((3),pos=13)
    relflux_virus = tb.Float32Col((3), pos=14)
    response_4540 = tb.Float32Col(pos=11)  # normalized for 360s
    xditherpos = tb.Float32Col((3))
    yditherpos = tb.Float32Col((3))
    xoffset = tb.Float32Col((3))
    yoffset = tb.Float32Col((3))
    xrms = tb.Float32Col((3))
    yrms = tb.Float32Col((3))
    nstars_fit = tb.Int32Col((3))

    obsind = tb.Int32Col()


def append_shot_to_table(shot, shottable, fn, cnt):
    F = fits.open(fn)

    try:
        filename = fn.name
    except:
        filename = fn

    idx = filename.find("exp")
    expn = filename[idx : idx + 5]

    if expn == "exp01":
        shot["obsind"] = cnt
        date = int("".join(F[0].header["DATE-OBS"].split("-")))
        obsid = int(F[0].header["OBSID"])
        shot["date"] = date
        shot["obsid"] = obsid

        shot["datevobs"] = str(date) + "v" + str(obsid).zfill(3)
        shot["shotid"] = int(str(date) + str(obsid).zfill(3))

        shot["objid"] = F[0].header["OBJECT"]
        shot["field"] = define_field(F[0].header["OBJECT"])

        shot["time"] = "".join(re.split("[:,.]", F[0].header["UT"]))[:7]

        shot["trajcra"] = F[0].header["TRAJCRA"] * 15.0
        shot["trajcdec"] = F[0].header["TRAJCDEC"]
        shot["trajcpa"] = F[0].header["PARANGLE"]
        shot["structaz"] = F[0].header["STRUCTAZ"]
        shot["ambtemp"] = F[0].header["AMBTEMP"]
        shot["humidity"] = F[0].header["HUMIDITY"]
        shot["dewpoint"] = F[0].header["DEWPOINT"]
        shot["pressure"] = F[0].header["BAROMPRE"]
        shot["expnum"] = shot["expnum"] + [int(expn[3:5]), 0, 0]
        shot["darktime"] = shot["darktime"] + [F[0].header["DARKTIME"], 0, 0]
        shot["exptime"] = shot["exptime"] + [F[0].header["EXPTIME"], 0, 0]
        shot["mjd"] = shot["mjd"] + [F[0].header["MJD"], 0, 0]

    elif expn == "exp02":
        shot["expnum"] = shot["expnum"] + [0, int(expn[3:5]), 0]
        shot["darktime"] = shot["darktime"] + [0, F[0].header["DARKTIME"], 0]
        shot["exptime"] = shot["exptime"] + [0, F[0].header["EXPTIME"], 0]
        shot["mjd"] = shot["mjd"] + [0, F[0].header["MJD"], 0]

    elif expn == "exp03":
        shot["expnum"] = shot["expnum"] + [0, 0, int(expn[3:5])]
        shot["darktime"] = shot["darktime"] + [0, 0, F[0].header["DARKTIME"]]
        shot["exptime"] = shot["exptime"] + [0, 0, F[0].header["EXPTIME"]]
        shot["mjd"] = shot["mjd"] + [0, 0, F[0].header["MJD"]]

    header = "header_" + expn
    shottable.attrs[header] = F[0].header

    return True


def append_fibers_to_table(fibindex, fib, im, fn, cnt, T, args):
    F = fits.open(fn)
    shotid = int(args.date) * 1000 + int(args.observation)

    if args.tar == True:
        ifuslot = "%03d" % int(F[0].header["IFUSLOT"])
        ifuid = "%03d" % int(F[0].header["IFUID"])
        specid = "%03d" % int(F[0].header["SPECID"])
        amp = F[0].header["NAME0"][21:23]
        multiframe = "multi_" + specid + "_" + ifuslot + "_" + ifuid + "_" + amp
    else:
        idx = fn.find("multi")
        multiframe = fn[idx : idx + 20]

    im["multiframe"] = multiframe
    n = F["spectrum"].data.shape[0]
    d = F["spectrum"].data.shape[1]
    if args.survey == "hdr1":
        attr = [
            "spectrum",
            "wavelength",
            "fiber_to_fiber",
            "twi_spectrum",
            "sky_spectrum",
            "sky_subtracted",
            "trace",
            "error1Dfib",
            "calfib",
            "calfibe",
            "Amp2Amp",
            "Throughput",
        ]
        imattr = ["image", "error", "clean_image"]
    else:
        attr = [
            "spectrum",
            "wavelength",
            "fiber_to_fiber",
            "sky_spectrum",
            "sky_subtracted",
            "trace",
            "error1Dfib",
            "calfib",
            "calfibe",
            "Amp2Amp",
            "Throughput",
            "chi2",
            "rms",
        ]
        imattr = ["image", "error", "clean_image"]

    for att in imattr:
        if att == "image":
            im[att] = F["PRIMARY"].data * 1.0
        else:
            im[att] = F[att].data * 1.0

    if args.tar:
        expn = fn.name[-12:-7]
    else:
        expn = op.basename(op.dirname(op.dirname(fn)))

    if T is not None:
        sel = T["col8"] == (multiframe + "_001.ixy")
        sel1 = T["col10"] == expn
        loc = np.where(sel * sel1)[0]
    for i in np.arange(n):
        fib["obsind"] = cnt
        fib["fibidx"] = i
        fib["fibnum"] = i + 1
        fibindex["fibidx"] = i
        fibindex["fibnum"] = i + 1
        fib["multiframe"] = multiframe
        fibindex["multiframe"] = multiframe
        fiberid = (
            str(shotid)
            + "_"
            + str(int(expn[-2:]))
            + "_"
            + multiframe
            + "_"
            + str(i + 1).zfill(3)
        )
        fib["fiber_id"] = fiberid
        fibindex["fiber_id"] = fiberid

        if T is not None:
            if len(loc):
                loci = loc[0] + i
                if isinstance(T["col1"][loci], float):
                    fib["ra"] = T["col1"][loci]
                    fib["dec"] = T["col2"][loci]
                    fib["fpx"] = T["col6"][loci]
                    fib["fpy"] = T["col7"][loci]
                else:
                    fib["ra"] = np.nan
                    fib["dec"] = np.nan
                    fib["fpx"] = np.nan
                    fib["fpy"] = np.nan
            else:
                fib["ra"] = np.nan
                fib["dec"] = np.nan
                fib["fpx"] = np.nan
                fib["fpy"] = np.nan
        else:
            fib["ra"] = np.nan
            fib["dec"] = np.nan
            fib["fpx"] = np.nan
            fib["fpy"] = np.nan

        fibindex["ra"] = fib["ra"]
        fibindex["dec"] = fib["dec"]
        fibindex["fpx"] = fib["fpx"]
        fibindex["fpy"] = fib["fpy"]

        fib["ifux"] = F["ifupos"].data[i, 0]
        fib["ifuy"] = F["ifupos"].data[i, 1]
        fibindex["ifux"] = fib["ifux"]
        fibindex["ifuy"] = fib["ifuy"]

        for att in attr:
            if att in F:
                fib[att] = F[att].data[i, :]
        fib["ifuslot"] = "%03d" % int(F[0].header["IFUSLOT"])
        fib["ifuid"] = "%03d" % int(F[0].header["IFUID"])
        fib["specid"] = "%03d" % int(F[0].header["SPECID"])
        fib["contid"] = F[0].header["CONTID"]
        try:
            fib["amp"] = "%s" % F[0].header["amp"][:2]
        except:
            fib["amp"] = "%s" % F[0].header["NAME0"][21:23]

        fib["expnum"] = int(expn[-2:])
        fibindex["expnum"] = int(expn[-2:])
        fib.append()
        fibindex.append()

    im["obsind"] = cnt
    im["ifuslot"] = "%03d" % int(F[0].header["IFUSLOT"])
    im["ifuid"] = "%03d" % int(F[0].header["IFUID"])
    im["specid"] = "%03d" % int(F[0].header["SPECID"])
    im["contid"] = F[0].header["CONTID"]
    try:
        im["amp"] = "%s" % F[0].header["amp"][:2]
    except:
        im["amp"] = "%s" % F.filename()[-7:-5]
    im["expnum"] = int(expn[-2:])
    im.append()

    # close the fits file
    F.close()
    return True


def main(argv=None):
    """ Main Function """
    # Call initial parser from init_utils
    parser = ap.ArgumentParser(description="""Create HDF5 file.""", add_help=True)

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
        "-r",
        "--rootdir",
        help="""Root Directory for Reductions""",
        type=str,
        default="/data/05350/ecooper/",
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
        help="""Appending to existing file.""",
        action="count",
        default=0,
    )

    parser.add_argument(
        "-dp",
        "--detect_path",
        help="""Path to detections""",
        type=str,
#        default="/work/00115/gebhardt/maverick/detect",
        default="/data/00115/gebhardt/detect")

    parser.add_argument(
        "-survey", "--survey", help="""{hdr1, hdr2, hdr3}""", type=str, default="hdr2"
    )

    parser.add_argument(
        "-tar", "--tar", help="""Flag to open tarred multifits""", action="store_true"
    )

    parser.add_argument("-tp", "--tmppath", type=str, default=os.getcwd())

    args = parser.parse_args(argv)
    args.log = setup_logging()

    global hdr_survey
    hdr_survey = args.survey
    if args.survey != hdr_survey:
        args.log.warning("Hard coded hdr_survey does not match input survey.")
        sys.exit()

    # Get the daterange over which reduced files will be collected
    files = get_files(args)
    datestr = "%sv%03d" % (args.date, int(args.observation))
    filepath = "%s/%s/dithall.use" % (args.detect_path, datestr)
    try:
        T = Table.read(filepath, format="ascii")
    except:
        T = None
        args.log.error("Could not open the dithall file from %s" % filepath)

    # Creates a new file if the "--append" option is not set or the file
    # does not already exist.
    does_exist = False
    if op.exists(args.outfilename) and args.append:
        fileh = tb.open_file(args.outfilename, "a")
        does_exist = True
    else:
        outfile = op.join(args.tmppath, args.outfilename)
        fileh = tb.open_file(outfile, "w")
        group = fileh.create_group(fileh.root, "Data", "VIRUS Fiber Data and Metadata")
        fileh.create_table(group, "Fibers", VIRUSFiber, "Fiber Info")
        fileh.create_table(fileh.root, "Shot", VIRUSShot, "Shot Info")
        fileh.create_table(group, "Images", VIRUSImage, "Image Info")
        fileh.create_table(group, "FiberIndex", VIRUSFiberIndex, "Fiber Coord Info")

    # Grab the fiber table and amplifier table for writing
    fibtable = fileh.root.Data.Fibers
    shottable = fileh.root.Shot
    imagetable = fileh.root.Data.Images
    fibindextable = fileh.root.Data.FiberIndex

    if does_exist:
        cnt = shottable[-1]["obsind"]
    else:
        cnt = 1

    if args.tar == True:
        
        shot = shottable.row
        n_ifu = {}
        for file_i in files:
            tar = tarfile.open(name=file_i, mode="r")

            members = tar.getmembers()
            fn = tar.extractfile(members[0])

            filename = fn.name
            idx = filename.find("exp")
            expn = filename[idx : idx + 5]

            n_ifu[expn] = int(len(members) / 4)

            success = append_shot_to_table(shot, shottable, fn, cnt)

            for member in members:
                fn = tar.extractfile(member)
                args.log.info("Working on %s" % member.name)
                fib = fibtable.row
                im = imagetable.row
                fibindex = fibindextable.row

                success = append_fibers_to_table(fibindex, fib, im, fn, cnt, T, args)
                if success:
                    fibtable.flush()
                    imagetable.flush()
                    fibindextable.flush()
                    
        shot["n_ifu"] = n_ifu["exp01"]
        shot.append()

    else:

        shot = shottable.row
        success = append_shot_to_table(shot, shottable, files[0], cnt)

        for fn in files:
            args.log.info("Working on %s" % fn)
            fib = fibtable.row
            im = imagetable.row
            fibindex = fibindextable.row

            success = append_fibers_to_table(fibindex, fib, im, fn, cnt, T, args)
            if success:
                fibtable.flush()
                imagetable.flush()
                fibindextable.flush()

    # create completely sorted index on the specid to make queries against that column much faster
    # specid chosen as the old multi*fits naming started with specid and it is fixed vs ifuslot and ifuid
    # for any given shot
    fibtable.cols.ra.create_csindex()
    fibindextable.cols.ra.create_csindex()

    imagetable.cols.multiframe.create_csindex()
    fibindextable.cols.multiframe.create_csindex()
    fibtable.cols.multiframe.create_csindex()

    shottable.flush()
    fibtable.flush()
    fibindextable.flush()
    imagetable.flush()

    fileh.close()

    # remove all temporary multifits
    if args.tar:
        datestr = "d%ss%03d" % (args.date, int(args.observation))
        datepath = op.join(args.tmppath, datestr)
        shutil.rmtree(datepath, ignore_errors=True)
        outfile = op.join(args.tmppath, args.outfilename)
        try:
            shutil.move(outfile, args.outfilename)
        except:
            os.remove(args.outfilename)
            shututil.move(outfile, args.outfilename)


if __name__ == "__main__":
    main()
