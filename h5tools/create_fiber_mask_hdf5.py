# -*- coding: utf-8 -*-
"""
Created: 2024/02/01

@author: Erin Mentuch Cooper

This script applies the HETDEX masking model to each
calibrated fiber in the data release

Example command line use:

python create_fiber_mask_hdf5.py -s 20170222007

"""

import glob
import sys
import os
import os.path as op
import argparse as ap
import numpy as np
import tables
import tables as tb

from scipy.interpolate import interp1d
from scipy.ndimage import binary_dilation

from astropy.table import Table, hstack
from astropy.nddata import bitmask
from astropy.nddata.bitmask import BitFlagNameMap

from hetdex_api.input_utils import setup_logging
from hetdex_api.config import HDRconfig
from hetdex_api.shot import get_fibers_table
from hetdex_api.survey import FiberIndex


# Set up bitflagnamemap for calfib data quality bitflag dictionary
class CALFIB_DQ(BitFlagNameMap):
    MAIN = 1
    FTF = 2
    CHI2FIB = 4
    BADPIX = 8
    BADAMP = 16
    LARGEGAL = 32
    METEOR = 64
    BADSHOT = 128
    THROUGHPUT = 256
    BADFIB = 512
    SATELLITE = 1024
    BADCAL = 2048
    PIXMASK = 4096
    
def main(argv=None):
    """Main Function"""
    # Call initial parser from init_utils
    parser = ap.ArgumentParser(
        description="""Append HDF5 Calibration info table.""", add_help=True
    )

    parser.add_argument(
        "-s",
        "--shotid",
        help="""shotid, e.g., 20230719012, YYYYMMDDOOO""",
        type=int,
        default=None,
    )

    parser.add_argument(
        "-v",
        "--version",
        help="""version: 1,2,3""",
        type=int,
        default=1,
    )

    parser.add_argument(
        "-survey",
        "--survey",
        help="""{hdr1, hdr2, hdr2.1, hdr3, hdr4, hdr5, pdr1}""",
        type=str,
        default="hdr5",
    )

    parser.add_argument(
        "-shot_h5",
        "--shot_h5",
        help="""A specific shot HDF5 file to update. If provided, only the shot HDF5 file is updated. There
        is no change to the survey and no separate fiber mask file is created.""",
        type=str,
        default=None
    )

    args = parser.parse_args(argv)
    args.log = setup_logging()

    if args.shot_h5 is not None:
        try:
            h5 = tables.open_file(args.shot_h5,mode="a")
        except:
            args.log(f"Fatal. Could not open {args.shot_h5} for appending.")
            sys.exit()


        if args.shotid is None: #shot_h5 is normally like "20250828v011.h5" ... if this fails, though, we want the exception
            args.shotid = int(os.path.basename(args.shot_h5)[:-3].replace("v",""))

        shotid_obj = args.shotid

        if args.survey is None:
            args.survey = "hdr5"

        config = HDRconfig(args.survey)

        # open bad pixel and bad fiber tables
        badpix = Table.read(
            config.badpix, format="ascii", names=["multiframe", "x1", "x2", "y1", "y2"]
        )

        #badshot = np.loadtxt(config.badshot, dtype=int)

        # Intiate the FiberIndex class from hetdex_api.survey:
        #F = FiberIndex(args.survey)

        # initiate rectified wavelength array
        wave_rect = np.linspace(3470, 5540, 1036)

        try:
            cols_to_keep = ['fiber_id', 'multiframe', 'fibnum', 'calfib', 'calfibe', 'trace', 'chi2', 'fiber_to_fiber',
                            'wavelength', 'spectrum']
            spec_tab = Table(h5.root.Data.Fibers.read())[cols_to_keep] #this is the long read

            if len(spec_tab) == 0 or spec_tab is None:
                args.log(f"Fatal. Could not read the Data.Fiber group (or is empty) in {args.shot_h5}")
                h5.close()
                sys.exit()
        except:
            #this is fatal
            args.log(f"Fatal. Could not read the Data.Fiber group in {args.shot_h5}")
            h5.close()
            sys.exit()

        #make a dummy (for now) mask_tab ... fiber level (not wavelength level)
        num_fibs = len(spec_tab)
        true_array = np.full(num_fibs,True)

        #add columns to spec_tab
        spec_tab['flag'] = true_array
        spec_tab['flag_badamp'] = true_array
        spec_tab['flag_badfib'] = true_array
        spec_tab['flag_meteor'] = true_array
        spec_tab['flag_satellite'] = true_array
        spec_tab['flag_largegal'] = true_array
        spec_tab['flag_shot'] = true_array
        spec_tab['flag_throughput'] = true_array

        del true_array
        #
        # #remember, True here means it is good
        # mask_tab = Table([true_array,true_array,true_array,true_array,true_array,true_array,true_array,true_array],
        #                  names=['flag', 'flag_badamp', 'flag_badfib', 'flag_meteor', 'flag_satellite',
        #                         'flag_largegal', 'flag_shot', 'flag_thoughput'])
        #
        # #spec_tab is same as alternate fib_tab, but without healpix
        # #the mask_tab is uniformly True, so don't need to match up to specific fibers
        # fib_tab = hstack([spec_tab, mask_tab])
        # del mask_tab

        # still todo is to set the flags we can in the now fib_tab
        # we rejoin the common processing farther below (after the else statement)

    else: #standard (original HETDEX)
        shotid_obj = args.shotid
        args.log.info("Working on {}".format(shotid_obj))

        date = str(shotid_obj)[0:8]
        obs = str(shotid_obj)[-3:]
        outfilename = "m{}v{}.h5".format(date, obs)

        fileh = tb.open_file(
            outfilename,
            "w",
            title="{} calibrated mask model. Version={}".format(args.survey, args.version),
        )

        config = HDRconfig(args.survey)

        # open bad pixel and bad fiber tables
        badpix = Table.read(
            config.badpix, format="ascii", names=["multiframe", "x1", "x2", "y1", "y2"]
        )

        badshot = np.loadtxt(config.badshot, dtype=int)

        # Intiate the FiberIndex class from hetdex_api.survey:
        F = FiberIndex(args.survey)

        # initiate rectified wavelength array
        wave_rect = np.linspace(3470, 5540, 1036)

        # acccess fiber data with selected columns to save memory
        try:
            cols_to_keep = ['fiber_id','multiframe','fibnum','calfib','calfibe', 'trace','chi2','fiber_to_fiber', 'wavelength', 'spectrum']
            spec_tab = get_fibers_table(shotid_obj, survey=args.survey)[cols_to_keep]
        except:
            if shotid_obj not in badshot:
                args.log.error("get_fibers_table failed for {}. Exiting".format(shotid_obj))
            fileh.close()
            F.close()
            args.log.info('Finished {}'.format(shotid_obj))
            sys.exit()

        if len(spec_tab) == 0 or spec_tab is None:
            # check if shotid is in badshot list
            if shotid_obj not in badshot:
                args.log.error("No fibers in {}. Exiting".format(shotid_obj))
            fileh.close()
            F.close()
            args.log.info('Finished {}'.format(shotid_obj))
            sys.exit()

        # access fiber index masks to propagate to fiber spectral masks
        idx = F.hdfile.root.FiberIndex.get_where_list("shotid==shotid_obj")
        fib_tab = Table(F.hdfile.root.FiberIndex.read_coordinates(idx))
        mask_tab = Table(F.fibermaskh5.root.Flags.read_coordinates(idx))
        fib_tab = hstack([fib_tab, mask_tab])

        del mask_tab
        # check to see that the hstacked table has matching fiber_ids
        for row in fib_tab:
            if row["fiber_id_1"] != row["fiber_id_2"]:
                args.log.error(
                    "Something wrong. Fiber id did not match up for {}:{}".format(
                        row["fiber_id_1"], row["fiber_id_2"]
                    )
                )

        fib_tab["fiber_id"] = fib_tab["fiber_id_1"].astype(str)
        fib_tab.remove_column("fiber_id_1")
        fib_tab.remove_column("fiber_id_2")

        # combined fiber table and associated masking with spectra table from get_fibers_table
        for i in np.arange(0, len(fib_tab)):
            if fib_tab["fiber_id"][i] != spec_tab["fiber_id"][i]:
                args.log.error(
                    "Something wrong with matching fiber_table and spec_table for {}".format(
                        shotid_obj
                    )
                )

        tab = hstack(
            [
                spec_tab,
                fib_tab[
                    "fiber_id",
                    "flag",
                    "flag_badamp",
                    "flag_badfib",
                    "flag_meteor",
                    "flag_satellite",
                    "flag_largegal",
                    "flag_shot",
                    "flag_throughput",
                ],
            ]
        )

        tab["fiber_id"] = tab["fiber_id_1"]
        tab.remove_column("fiber_id_1")
        tab.remove_column("fiber_id_2")

        del fib_tab
        spec_tab = tab
        del tab

#break here and recombine logic
    nfibs = np.shape(spec_tab)[0]
    wavelength = np.array(spec_tab["wavelength"])
    calfib = np.array(spec_tab["calfib"])
    calfibe = np.array(spec_tab["calfibe"])
    ftf = np.array(spec_tab["fiber_to_fiber"])

    chi2fib = np.zeros_like(calfib)
    trace_data = np.zeros_like(calfib)
    spectrum = np.zeros_like(calfib)

    # interpolate chi2, trace, spectrum into rectified wavelength
    for i in np.arange(0, nfibs):
        chi2_interp = interp1d(
            np.array(spec_tab["wavelength"][i]),
            np.array(spec_tab["chi2"][i]),
            fill_value="extrapolate",
            kind="nearest",
        )
        chi2fib[i] = chi2_interp(wave_rect)

        trace_interp = interp1d(
            np.arange(1, 1033),
            np.array(spec_tab["trace"][i]),
            fill_value="extrapolate",
            kind="nearest",
        )
        trace_data[i] = trace_interp(np.arange(1, 1037))

        spec_interp = interp1d(
            np.array(spec_tab["wavelength"][i]),
            np.array(spec_tab["spectrum"][i]),
            fill_value="extrapolate",
            kind="nearest",
        )

        spectrum[i] = spec_interp(wave_rect)

    # Get native pixel mask when spectrum==0, expand +/- 1 in wave dim
    pixmask = spectrum != 0
    structure = np.array([[1, 1, 1]])
    # Apply binary dilation along axis=1
    mask_pixmask = binary_dilation(pixmask, structure=structure)

    # Get bad pix mask
    mask_badpix = np.ones_like(calfib, dtype=bool)
    x_rect = [np.linspace(1, 1032, 1036) for row in calfib]

    mf_array = np.dstack([spec_tab["multiframe"].astype(str)] * 1036)[0]

    for row in badpix:
        sel_mf = mf_array == row["multiframe"]
        if np.sum(sel_mf) > 0:
            sel_reg_x = (x_rect >= row["x1"]) & (x_rect <= row["x2"])
            sel_reg_y = (trace_data >= row["y1"]) & (trace_data <= row["y2"])
            mask_badpix[sel_mf * sel_reg_x * sel_reg_y] = False

    mask_main = (calfibe > 1e-8) | (calfib != 0.0)
    mask_ftf_per_fib = np.median(ftf, axis=1) > 0.5  # only one per fiber
    mask_ftf = np.dstack([mask_ftf_per_fib] * 1036)[0]

    mask_chi2fib = chi2fib < 5 #Much Stricter. Updated 2025-03-28 EMC (from 150)

    # add badcal mask
    mask_badcal =  np.ones_like(calfib, dtype=bool)

    sel_wave =  (wave_rect >= 3534) * (wave_rect <= 3556)
    mask_badcal[:, sel_wave] = 0


    cal5200_tab = Table.read(
        config.cal5200,
        format="ascii",
        names=["shotid", "multiframe", "expnum"],
    )
    cal5460_tab = Table.read(
        config.cal5460,
        format="ascii",
        names=["shotid", "multiframe", "expnum"],
    )

    mfs = np.array( spec_tab['multiframe']).astype(str)

    for	row in cal5200_tab[ cal5200_tab['shotid'] == shotid_obj]:
        sel_mf = mfs == row["multiframe"]
        if np.sum(sel_mf)>0:
            mask_badcal[sel_mf, 862:868] = False

    for row in cal5460_tab[ cal5460_tab['shotid'] == shotid_obj]:
        sel_mf = mfs == row["multiframe"]
        mask_badcal[sel_mf, 993: 998] = 0

    # Add full fiber flags
    mask_amp = np.dstack([spec_tab["flag_badamp"]] * 1036)[0]
    mask_gal = np.dstack([spec_tab["flag_largegal"]] * 1036)[0]
    mask_meteor = np.dstack([spec_tab["flag_meteor"]] * 1036)[0]
    mask_shot = np.dstack([spec_tab["flag_shot"]] * 1036)[0]
    mask_throughput = np.dstack([spec_tab["flag_throughput"]] * 1036)[0]
    mask_satellite = np.dstack([spec_tab['flag_satellite']] * 1036)[0]
    mask_badfib = np.dstack([spec_tab['flag_badfib']] * 1036)[0]

    CALFIB_NET = (
        CALFIB_DQ.MAIN * np.invert(mask_main)
        + CALFIB_DQ.FTF * np.invert(mask_ftf)
        + CALFIB_DQ.CHI2FIB * np.invert(mask_chi2fib)
        + CALFIB_DQ.BADPIX * np.invert(mask_badpix)
        + CALFIB_DQ.BADAMP * np.invert(mask_amp)
        + CALFIB_DQ.LARGEGAL * np.invert(mask_gal)
        + CALFIB_DQ.METEOR * np.invert(mask_meteor)
        + CALFIB_DQ.BADSHOT * np.invert(mask_shot)
        + CALFIB_DQ.THROUGHPUT * np.invert(mask_throughput)
        + CALFIB_DQ.BADFIB * np.invert(mask_badfib)
        + CALFIB_DQ.SATELLITE * np.invert(mask_satellite)
        + CALFIB_DQ.BADCAL * np.invert(mask_badcal)
        + CALFIB_DQ.PIXMASK * np.invert(mask_pixmask)
    )

    if args.shot_h5 is not None:
        #add new table to h5 file
        flags = h5.create_table(
            h5.root, #put in the same spot so other callers can use the shot_h5 with the same table name resolution
            "CalfibDQ",
            Table(
                [spec_tab["fiber_id"], CALFIB_NET.astype(np.int16)],
                names=["fiber_id", "calfib_dq"],
            ).as_array(),
        )
        flags.flush()
        h5.close()

        args.log.info(f"Finished {args.shot_h5}")
    else:
        flags = fileh.create_table(
            fileh.root,
            "CalfibDQ",
            Table(
                [spec_tab["fiber_id"], CALFIB_NET.astype(np.int16)],
                names=["fiber_id", "calfib_dq"],
            ).as_array(),
        )
        flags.flush()
        fileh.close()
        F.close()
        args.log.info('Finished {}'.format(shotid_obj))

if __name__ == "__main__":
    main()
