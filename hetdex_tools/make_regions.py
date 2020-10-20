# -*- coding: utf-8 -*-
"""
Created on September 23 2020
@author: Erin Mentuch Cooper
"""

import sys
import argparse as ap

from astropy.coordinates import SkyCoord
from astropy import units as u


from regions import RectangleSkyRegion, write_ds9

from hetdex_api.survey import Survey
from hetdex_api.flux_limits.hdf5_sensitivity_cubes import (
    SensitivityCubeHDF5Container,
    return_sensitivity_hdf_path,
)
from hetdex_api.input_utils import setup_logging


def get_regions_from_flim(shotid):

    date = str(shotid)[0:8]
    obs = str(shotid)[8:]
    datevobs = str(date) + "v" + str(obs).zfill(3)
    shotid = int(str(date) + str(obs).zfill(3))

    hdf_filename, mask_fn = return_sensitivity_hdf_path(
        datevobs, release="hdr2.1", return_mask_fn=True
    )

    hdfcont = SensitivityCubeHDF5Container(
        filename=hdf_filename,
        aper_corr=1.0,
        flim_model="hdr2pt1",
        mask_filename=mask_fn,
    )

    ifu_name_list = []
    ifu_ra = []
    ifu_dec = []
    ifuregions = []

    for ifu_name, tscube in hdfcont.itercubes():
        shape = tscube.sigmas.shape
        ra, dec, lambda_ = tscube.wcs.all_pix2world(
            shape[2] / 2.0, shape[1] / 2.0, shape[0] / 2.0, 0
        )
        pa = -tscube.header["CROTA2"]
        ifu_name_list.append(ifu_name)
        coord = SkyCoord(ra, dec, unit="deg")
        ifu_ra.append(ra)
        ifu_dec.append(dec)
        ifuregions.append(
            RectangleSkyRegion(
                center=coord,
                width=1.0 * u.arcmin,
                height=1.0 * u.arcmin,
                angle=pa * u.deg,
            )
        )
    hdfcont.close()

    return ifuregions


def main(argv=None):
    """ Main Function """
    # Call initial parser from init_utils
    parser = ap.ArgumentParser(description="""Create Region Files.""", add_help=True)

    parser.add_argument(
        "-s",
        "--shotid",
        help="""Shot identifier, an integer""",
        type=int,
        default=None,
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
        "-f", "--field", help="""Options=""", type=str, default=None,
    )

    args = parser.parse_args(argv)
    args.log = setup_logging()

    if args.field is not None:

        S = Survey()
        survey_table = S.return_astropy_table()
        sel_field = survey_table["field"] == args.field
        ifuregions = []

        for row in survey_table[sel_field]:
            args.log.info("Working on " + str(row["shotid"]))
            ifuregions.extend(get_regions_from_flim(row["shotid"]))

        outname = args.field

    elif args.shotid:
        ifuregions = get_regions_from_flim(args.shotid)
        outname = str(args.shotid)

    elif args.date is not None:
        if args.observation is not None:
            shotid = int(str(args.date) + str(args.observation).zfill(3))
            ifuregions = get_regions_from_flim(shotid)
            outname = str(shotid)

    region_file = outname + '.reg'
    write_ds9(ifuregions, region_file)


if __name__ == "__main__":
    main()
