# -*- coding: utf-8 -*-
"""
Created on October 8 2019
@author: Erin Mentuch Cooper

Extracts 1D spectrum at specified RA/DEC.

Provide either an ID/RA/DEC as input or a file with a list of ID, RA, DEC. Can
be in the format of an astropy table (fits, ascii and no_header are currently
supported) provided columns are appropriately marked 'ID','ra', 'dec' or 'RA',
'DEC'. Coordinates can be in any astropy unit coordinate in the input table.

Option can be to provide a single shotid/datevobs to work on for batch
processing, or can let program figure out which shots to extract spectra on.

Parameters
----------

ra
    RA for single source in an astropy unit, or assumed as deg
dec
    DEC for single source in an astropy unit, or assumed as deg
rad
    radius size for aperture in arcsec or an astropy unit, defaults to 3arcsec
ID
    object ID if using a single source
shotid
    use if you are running on just a single shotid (or datevobs)
input
    path to input catalog of ID/RA/DEC, can use any astropy units or program
    will assume degree
output
    name of output pickle file
multiprocess
    flag to use python multiprocess, don't use in a slurm job
single
    flag to write out several astropy tables for each ID/shot spectra
merge
    will combine all .pkl files in a directory, for use in slurm job cleanup
mergepath
    use if you want to combine pickle files in another directory

Examples
--------

For a single coordinate, you can search and extract all shots in HDR1

python3 get_spec.py --ra 150.02548 --dec 2.087987 --ID cosmos_LAE --outfile cosmos_LAE

To speed up, you can use the --multiprocess option. This works for interactive use
on a node either through a jupyter notebook or an interactive node called by idev

python3 get_spec.py --multiprocess --ra 150.02548 --dec 2.087987 --ID cosmos_LAE --outfile cosmos_LAE  

For working on a single shot only you can use the -s option. 
It is recommended to batch process each 
SHOTID individually and provide the same source list in each call.

python3 get_spec.py -ra 8.86535 -dec 0.59352  -s 20190104008

This can also work on an input file. It should either have the format of
3 columns = ID/RA/DEC or be an astropy table with columns 
labeled 'ID', 'ra' or 'RA', 'dec' or 'DEC'

python3 get_spec.py --multiprocess -i '3dhst_input.cat' -o '3dhst' 


"""

import sys
import argparse as ap
import os
import os.path as op
import glob
import re

import numpy as np
import pickle
import warnings
import logging
from hetdex_api.input_utils import setup_logging

import types

from astropy.coordinates import SkyCoord
from astropy.table import Table, Column, vstack
import astropy.units as u

from hetdex_api.extract import Extract
from hetdex_api.survey import Survey

from copy import deepcopy
from collections.abc import Mapping

from multiprocessing import Process, Manager
import time

if not sys.warnoptions:
    warnings.simplefilter("ignore")


def merge(dict1, dict2):
    """ Return a new dictionary by merging two dictionaries recursively. """
    result = deepcopy(dict1)

    for key, value in dict2.items():
        if isinstance(value, Mapping):
            result[key] = merge(result.get(key, {}), value)
        else:
            result[key] = deepcopy(dict2[key])

    return result


def get_source_spectra(shotid, args):
    E = Extract()
    source_dict = {}

    if args.survey == "hdr1":
        source_num_switch = 20
    else:
        source_num_switch = 0

    if len(args.matched_sources[shotid]) > 0:
        args.log.info("Working on shot: %s" % shotid)
        if args.survey == "hdr1":
            fwhm = args.survey_class.fwhm_moffat[args.survey_class.shotid == shotid][0]
        else:
            fwhm = args.survey_class.fwhm_virus[args.survey_class.shotid == shotid][0]

        moffat = E.moffat_psf(fwhm, 10.5, 0.25)

        if len(args.matched_sources[shotid]) > source_num_switch:
            E.load_shot(shotid, fibers=True, survey=args.survey)
        else:
            E.load_shot(shotid, fibers=False, survey=args.survey)

        for ind in args.matched_sources[shotid]:
            if np.size(args.coords) > 1:
                info_result = E.get_fiberinfo_for_coord(
                    args.coords[ind], radius=args.rad, ffsky=args.ffsky
                )
            else:
                info_result = E.get_fiberinfo_for_coord(
                    args.coords, radius=args.rad, ffsky=args.ffsky
                )
            if info_result is not None:
                if np.size(args.ID) > 1:
                    args.log.info("Extracting %s" % args.ID[ind])
                else:
                    args.log.info("Extracting %s" % args.ID)

                ifux, ifuy, xc, yc, ra, dec, data, error, mask = info_result

                weights = E.build_weights(xc, yc, ifux, ifuy, moffat)
                result = E.get_spectrum(data, error, mask, weights)
                spectrum_aper, spectrum_aper_error = [res for res in result]

                #add in the total weight of each fiber (as the sum of its weight per wavebin)
                if args.fiberweights:
                    try:
                        fiber_weights = np.array( [x for x in zip(ra, dec, np.sum(weights*mask, axis=1))])
                    except:
                        fiber_weights = []
                else:
                    fiber_weights = []

                if np.size(args.ID) > 1:
                    if args.ID[ind] in source_dict:
                        source_dict[args.ID[ind]][shotid] = [
                            spectrum_aper,
                            spectrum_aper_error,
                            weights.sum(axis=0),
                            fiber_weights,
                        ]
                    else:
                        source_dict[args.ID[ind]] = dict()
                        source_dict[args.ID[ind]][shotid] = [
                            spectrum_aper,
                            spectrum_aper_error,
                            weights.sum(axis=0),
                            fiber_weights,
                        ]
                else:
                    if args.ID in source_dict:
                        source_dict[args.ID][shotid] = [
                            spectrum_aper,
                            spectrum_aper_error,
                            weights.sum(axis=0),
                            fiber_weights,
                        ]
                    else:
                        source_dict[args.ID] = dict()
                        source_dict[args.ID][shotid] = [
                            spectrum_aper,
                            spectrum_aper_error,
                            weights.sum(axis=0),
                            fiber_weights,
                        ]

        E.shoth5.close()
        
    return source_dict
        

def get_source_spectra_mp(source_dict, shotid, manager, args):
    E = Extract()

    if args.survey == "hdr1":
        source_num_switch = 20
    else:
        source_num_switch = 0

    if len(args.matched_sources[shotid]) > 0:
        args.log.info("Working on shot: %s" % shotid)
        if args.survey == "hdr1":
            fwhm = args.survey_class.fwhm_moffat[args.survey_class.shotid == shotid][0]
        else:
            fwhm = args.survey_class.fwhm_virus[args.survey_class.shotid == shotid][0]

        moffat = E.moffat_psf(fwhm, 10.5, 0.25)

        if len(args.matched_sources[shotid]) > source_num_switch:
            E.load_shot(shotid, fibers=True, survey=args.survey)
        else:
            E.load_shot(shotid, fibers=False, survey=args.survey)

        for ind in args.matched_sources[shotid]:
            if np.size(args.coords) > 1:
                info_result = E.get_fiberinfo_for_coord(
                    args.coords[ind], radius=args.rad, ffsky=args.ffsky
                )
            else:
                info_result = E.get_fiberinfo_for_coord(
                    args.coords, radius=args.rad, ffsky=args.ffsky
                )
            if info_result is not None:
                if np.size(args.ID) > 1:
                    args.log.info("Extracting %s" % args.ID[ind])
                else:
                    args.log.info("Extracting %s" % args.ID)

                ifux, ifuy, xc, yc, ra, dec, data, error, mask = info_result
                weights = E.build_weights(xc, yc, ifux, ifuy, moffat)

                result = E.get_spectrum(data, error, mask, weights)

                spectrum_aper, spectrum_aper_error = [res for res in result]
                sel = np.isfinite(spectrum_aper)

                #add in the total weight of each fiber (as the sum of its weight per wavebin)
                if args.fiberweights:
                    
                    try:
                        fiber_weights = np.array( [x for x in zip(ra,dec,np.sum(weights*mask,axis=1))])
                    except:
                        fiber_weights = []
                else:
                    fiber_weights = []

                if np.size(args.ID) > 1:
                    if args.ID[ind] in source_dict:
                        source_dict[args.ID[ind]][shotid] = [
                            spectrum_aper,
                            spectrum_aper_error,
                            weights.sum(axis=0),
                            fiber_weights,
                        ]
                    else:
                        source_dict[args.ID[ind]] = manager.dict()
                        source_dict[args.ID[ind]][shotid] = [
                            spectrum_aper,
                            spectrum_aper_error,
                            weights.sum(axis=0),
                            fiber_weights,
                        ]
                else:
                    if args.ID in source_dict:
                        source_dict[args.ID][shotid] = [
                            spectrum_aper,
                            spectrum_aper_error,
                            weights.sum(axis=0),
                            fiber_weights,
                        ]
                    else:
                        source_dict[args.ID] = manager.dict()
                        source_dict[args.ID][shotid] = [
                            spectrum_aper,
                            spectrum_aper_error,
                            weights.sum(axis=0),
                            fiber_weights,
                        ]

        E.shoth5.close()
        
    return source_dict
        

def return_astropy_table(Source_dict, fiberweights=False):
    """Returns an astropy table fom a source dictionary"""

    id_arr = []
    shotid_arr = []
    wave_arr = []
    spec_arr = []
    spec_err_arr = []
    weights_arr = []
    fiber_weights_arr = []

    # loop over every ID/observation combo:

    for ID in Source_dict.keys():

        for shotid in Source_dict[ID].keys():
            wave_rect = 2.0 * np.arange(1036) + 3470.0
            spec = Source_dict[ID][shotid][0]
            spec_err = Source_dict[ID][shotid][1]
            weights = Source_dict[ID][shotid][2]
            if fiberweights:
                fiber_weights = Source_dict[ID][shotid][3]
            
            sel = np.isfinite(spec)

            if np.sum(sel) > 0:
                id_arr.append(ID)
                shotid_arr.append(shotid)
                wave_arr.append(wave_rect)
                spec_arr.append(spec)
                spec_err_arr.append(spec_err)
                weights_arr.append(weights)
                if fiberweights:
                    fiber_weights_arr.append(fiber_weights)

    output = Table()
    fluxden_u = 1e-17 * u.erg * u.s ** (-1) * u.cm ** (-2) * u.AA ** (-1)

    output.add_column(Column(id_arr), name="ID")
    output.add_column(Column(shotid_arr), name="shotid")
    output.add_column(Column(wave_arr, unit=u.AA, name="wavelength"))
    output.add_column(Column(spec_arr, unit=fluxden_u, name="spec"))
    output.add_column(Column(spec_err_arr, unit=fluxden_u, name="spec_err"))
    output.add_column(Column(weights_arr), name="weights")
    if fiberweights:
        output.add_column(Column(fiber_weights_arr), name="fiber_weights")

    return output


def get_spectra_dictionary(args):

    args.matched_sources = {}
    shots_of_interest = []

    count = 0

    # this radius applies to the inital shot search and requires a large
    # aperture for the wide FOV of VIRUS
    max_sep = 11.0 * u.arcmin

    args.log.info("Finding shots of interest")

    for i, coord in enumerate(args.survey_class.coords):
        dist = args.coords.separation(coord)
        sep_constraint = dist < max_sep
        shotid = args.survey_class.shotid[i]
        idx = np.where(sep_constraint)[0]

        if np.size(idx) > 0:

            args.matched_sources[shotid] = idx
            count += np.size(idx)
            if len(idx) > 0:
                shots_of_interest.append(shotid)

    args.log.info("Number of shots of interest: %i" % len(shots_of_interest))
    args.log.info("Extracting %i sources" % count)

    if args.multiprocess:

        manager = Manager()
        Sources = manager.dict()

        start = time.time()
        njobs = np.size(shots_of_interest)

        ntasks = 32

        for i in np.arange(0, njobs, ntasks):
            jobs = [
                Process(
                    target=get_source_spectra_mp, args=(Sources, shotid, manager, args)
                )
                for shotid in shots_of_interest[i : np.minimum(njobs, i + ntasks)]
            ]

            for j in jobs:
                j.start()
            for j in jobs:
                j.join()

        end = time.time()
        args.log.info(
            "Extraction of sources completed in %.2f minutes." % ((end - start) / 60.0)
        )

        # convert manager dict to a regular dictionary object to pickle
        Source_dict = {}
        count = 0
        for ID in Sources.keys():
            Source_dict[ID] = {}
            for shotid in Sources[ID].keys():
                count += 1
                Source_dict[ID][shotid] = list(Sources[ID][shotid])

    else:
        Source_dict = {}
        for shotid in shots_of_interest:
            shot_source_dict = get_source_spectra(shotid, args)

            if shot_source_dict:
                Source_dict = merge(Source_dict, shot_source_dict)

    return Source_dict


def get_parser():
    """ function that returns a parser from argparse """

    parser = ap.ArgumentParser(
        description="""Extracts 1D spectrum at specified RA/DEC""", add_help=True
    )

    parser.add_argument(
        "-s",
        "--shotid",
        help="""ShotID, e.g., 20170321v009, YYYYMMDDvOBS""",
        type=str,
        default=None,
    )

    parser.add_argument(
        "-i", "--infile", help="""File with table of ID/RA/DEC""", default=None
    )

    parser.add_argument(
        "-o",
        "--outfile",
        help="""File to store dictionary output""",
        default="output",
        type=str,
    )

    parser.add_argument(
        "-ra",
        "--ra",
        help="""ra, e.g., right ascension in degrees""",
        type=float,
        default=None,
    )

    parser.add_argument(
        "-dec",
        "--dec",
        help="""ra, e.g., right ascension in degrees""",
        type=float,
        default=None,
    )

    parser.add_argument(
        "-rad",
        "--rad",
        help="""radius, e.g., aperture radius in arcsec""",
        type=float,
        default=3.0,
    )

    parser.add_argument("-id", "--ID", help="""source name""", type=str, default=None)

    parser.add_argument(
        "--multiprocess",
        "-mp",
        help="""Multiprocessing Flag""",
        default=False,
        required=False,
        action="store_true",
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
        "--mergepath",
        "-mpath",
        help="""Path to location of files to merge""",
        default=os.getcwd(),
        type=str,
    )

    parser.add_argument(
        "--single",
        "-single",
        help="""Select to create single astropy tables for each ID/SHOT""",
        default=False,
        required=False,
        action="store_true",
    )

    parser.add_argument(
        "--fits",
        "-fit",
        help="""Flag to spectra in an astropy table saved as a fits file.""",
        default=True,
        action="store_true",
    )

    parser.add_argument(
        "--pickle",
        "-pkl",
        help="""Flag to store spectra in a pkl file""",
        default=False,
        action="store_true",
    )

    parser.add_argument(
        "--survey",
        "-survey",
        type=str,
        help="""Data Release you want to access""",
        default="hdr2.1",
    )

    parser.add_argument("-tpmin", "--tpmin", type=float, default=None)

    parser.add_argument(
        "--ffsky",
        "-ffsky",
        help="""Set to True to use the full frame sky sutraction.""",
        default=False,
        required=False,
        action="store_true",
    )

    parser.add_argument(
        "--fiberweights",
        "-fw",
        help="""Set to retrieve fibers weight. Warning this is slow on large sourcelists.""",
        default=False,
        required=False,
        action="store_true",
    )

        
    return parser


def main(argv=None):
    """ Main Function """

    parser = get_parser()
    args = parser.parse_args(argv)
    args.log = setup_logging()

    if args.pickle:
        args.fits = False

    if args.merge:

        if args.fits:
            master_table = Table()
            files = glob.glob(op.join(args.mergepath, "*.fits"))
            args.log.info("Merging all fits files in " + args.mergepath)

            for file in files:
                file_table = Table.read(open(file, "rb"))
                if np.size(file_table) > 0:
                    master_table = vstack([master_table, file_table])
            outfile = args.outfile + ".fits"
            master_table.write(outfile, format="fits", overwrite=True)

        else:
            all_source_dict = {}
            files = glob.glob(op.join(args.mergepath, "*.pkl"))
            args.log.info("Merging all pickle files in " + args.mergepath)
            for file in files:
                file_dict = pickle.load(open(file, "rb"))
                if len(file_dict) > 0:
                    all_source_dict = merge(all_source_dict, file_dict)

            outfile = args.outfile + ".pkl"
            pickle.dump(all_source_dict, open(outfile, "wb"))

        args.log.info("Saved output file to " + outfile)
        sys.exit("Exiting")

    if args.infile:

        args.log.info("Loading External File")

        try:
            try:
                table_in = Table.read(args.infile, format="ascii")
                if table_in.colnames == ["col1", "col2", "col3"]:
                    table_in["col1"].name = "ID"
                    table_in["col2"].name = "ra"
                    table_in["col3"].name = "dec"
                elif np.size(table_in.colnames) != 3:
                    args.log.info("Input file not in three column format")
            except Exception:
                pass
            try:
                table_in = Table.read(args.infile, format="fits")
            except Exception:
                pass
        except Exception:
            if op.exists(args.infile):
                args.log.warning("Could not open input file")
                sys.exit("Exiting")
            else:
                args.log.warning("Input file not found")
                sys.exit("Exiting")
        try:
            args.ID = table_in["ID"]
        except:
            args.ID = table_in["id"]

        try:
            args.ra = table_in["ra"]
            args.dec = table_in["dec"]
        except:
            args.ra = table_in["RA"]
            args.dec = table_in["DEC"]

    else:
        if args.ID == None:
            if np.size(args.ra) > 1:
                args.ID = str(np.arange(1, np.size(table_in) + 1)).zfill(9)
            else:
                args.ID = 1

        args.log.info("Extracting for ID: %s" % args.ID)

    # generate astropy coordinates object for searching

    if re.search(":", str(args.ra)):
        args.coords = SkyCoord(args.ra, args.dec, unit=(u.hourangle, u.deg))
    else:
        args.coords = SkyCoord(args.ra, args.dec, unit=u.deg)

    S = Survey(args.survey)
    
    ind_good_shots = S.remove_shots()

    if args.tpmin:
        ind_tp = S.response_4540 > args.tpmin
        args.survey_class = S[ind_good_shots * ind_tp]
    else:
        args.survey_class = S[ind_good_shots]
        
    # if args.shotidid exists, only select those shots

    if args.shotid:
        try:
            sel_shot = args.survey_class.shotid == int(args.shotid)
        except Exception:
            sel_shot = args.survey_class.datevobs == str(args.shotid)

        args.survey_class = args.survey_class[sel_shot]

    else:
        pass

    # main function to retrieve spectra dictionary
    Source_dict = get_spectra_dictionary(args)

    args.survey_class.close()

    if args.pickle:
        outfile = args.outfile + ".pkl"
        pickle.dump(Source_dict, open(outfile, "wb"))

    if args.single:
        # loop over every ID/observation combo:
        fluxden_u = 1e-17 * u.erg * u.s ** (-1) * u.cm ** (-2) * u.AA ** (-1)
        for ID in Source_dict.keys():
            for shotid in Source_dict[ID].keys():

                wave_rect = 2.0 * np.arange(1036) + 3470.0
                spec = Source_dict[ID][shotid][0]
                spec_err = Source_dict[ID][shotid][1]
                weights = Source_dict[ID][shotid][2]

                sel = np.isfinite(spec)
                if np.sum(sel) > 0:
                    output = Table()

                    output.add_column(Column(wave_rect, name="wavelength", unit=u.AA))
                    output.add_column(Column(spec, name="spec", unit=fluxden_u))
                    output.add_column(Column(spec_err, name="spec_err", unit=fluxden_u))
                    output.add_column(Column(weights, name="weights"))

                    output.write(
                        "spec_" + str(ID) + "_" + str(shotid) + ".tab", format="ascii"
                    )

    if args.fits:
        output = return_astropy_table(Source_dict, fiberweights=args.fiberweights)
        if args.fiberweights:
            # cannot save fiberweights to a fits file
            output.remove_column('fiber_weights')
        output.write(args.outfile + ".fits", format="fits", overwrite=True)


# tables.file._open_files.close_all()

if __name__ == "__main__":
    main()


def get_spectra(
    coords,
    ID=None,
    rad=3.0,
    multiprocess=True,
    shotid=None,
    survey="hdr2.1",
    tpmin=None,
    ffsky=False,
    fiberweights=False
):
    """
    Function to retrieve PSF-weighted, ADR and aperture corrected
    spectral extractions of HETDEX fibers. It will search all shots
    within a specific HETDEX Data Release and return a table of
    spectra for each extraction per shot in which more than 7 fibers
    are found in order to generate an extracted spectrum.

    Parameters
    ----------
    coords
        list astropy coordinates
    ID
        list of ID names (must be same length as coords). Will
        generate a running index if no ID is given
    rad
        radius of circular aperture to be extracted in arcsec.
        Default is 3.0
    multiprocess
        boolean flag to use multiprocessing. This will greatly
        speed up its operation as it will extract on 32 shots at
        time. But only use this when on a compute node. Use
        idev, a jupyter notebook, or submit the job as a single
        python slurm job.
    shotid
        list of integer shotids to do extractions on. By default
        it will search the whole survey except for shots located
        in the bad.shotlist file
    survey
        Survey you want to access. User note that HDR1 extractions
        are much slower compared to HDR2.
    tpmin
        Include only shots above tpmin. Default is None.
    ffsky
        Use the full frame 2D sky subtraction model. Default is
        to use the local sky subtracted, flux calibrated fibers.
    fiberweights
        Boolean flag to include fiber_weights tuple in source
        dictionary. This is used in Elixer, but is slow
        when used on large source lists.
    
    Returns
    -------
    sources
        an astropy table object of source spectra for all input
        coords/ID that have spectra in the survey shots. There
        is one row per source ID/shotid observation.
    """

    args = types.SimpleNamespace()

    args.multiprocess = multiprocess
    args.coords = coords
    args.rad = rad * u.arcsec
    args.survey = survey

    args.ffsky = ffsky
    args.fiberweights = fiberweights

    S = Survey(survey)
    ind_good_shots = S.remove_shots()

    if tpmin:
        ind_tp = S.response_4540 > tpmin
        args.survey_class = S[ind_good_shots * ind_tp]
    else:
        args.survey_class = S[ind_good_shots]
        
    if shotid is not None:
        try:
            if np.size(shotid) == 1:
                sel_shot = args.survey_class.shotid == int(shotid)
                # shut off multiproces flag if its just one shot
                args.multiprocess = False
            else:
                sel_shot = np.zeros(np.size(args.survey_class.shotid), dtype=bool)

                for shot_i in shotid:

                    sel_i = args.survey_class.shotid == int(shot_i)
                    sel_shot = np.logical_or(sel_shot, sel_i)

        except Exception:
            sel_shot = args.survey_class.datevobs == str(shotid)

        args.survey_class = args.survey_class[sel_shot]
    else:
        pass
        # sel_shot = args.survey_class.shotid > 20171200000
        # args.survey_class = args.survey_class[sel_shot]

    args.log = setup_logging()

    args.log.setLevel(logging.INFO)

    nobj = np.size(args.coords)

    if ID is None:
        if nobj > 1:
            args.ID = np.arange(1, nobj + 1)
        else:
            args.ID = 1
    else:
        args.ID = ID

    Source_dict = get_spectra_dictionary(args)

    args.survey_class.close()

    output = return_astropy_table(Source_dict, fiberweights=args.fiberweights)

    args.log.info("Retrieved " + str(np.size(output)) + " spectra.")

    return output
