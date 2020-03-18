# -*- coding: utf-8 -*-
"""

get_shots_of_interest.py
========================

For a catalog of coordinates, find the shots
of interest to search for extrations using
get_spec.py or get_spec2D.py. -i or --infile
should be a file containing eithe 3 columns of ID/RA/DEC
or an astropy table file with columns labeled as
'ID', 'ra', 'dec'.

Examples
--------

python3 get_shots_of_interest.py -i '3dhst_intput.cat'

Returns
-------
This will produce a txt file of shotids 'shotlist.txt'

"""

import argparse as ap
import numpy as np

from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord

from hetdex_api.survey import Survey
from hetdex_api.input_utils import setup_logging


def get_parser():
    parser = ap.ArgumentParser(
        description="""Create a list of shotids for an input source list""",
        add_help=True,
    )

    parser.add_argument(
        "-i",
        "--infile",
        help="""File with table of ID/RA/DEC""",
        default=None
    )

    parser.add_argument(
        "-o",
        "--outfile",
        help="""File to store shotlist info""",
        default="shotlist",
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

    parser.add_argument("-id",
                        "--ID",
                        help="""source name""",
                        type=str,
                        default=None)
    return parser


def main(argv=None):

    parser = get_parser()
    args = parser.parse_args(argv)
    args.log = setup_logging()

    if args.infile:

        args.log.info("Loading External File")

        table_in = Table.read(args.infile, format="ascii")
        args.ID = table_in["ID"]
        args.ra = table_in["ra"]
        args.dec = table_in["dec"]
    else:
        if args.ID is None:
            args.ID = "DEX" + str(args.ra).zfill(4) + "_"
            + str(args.dec).zfill(4)

        args.log.info("Extracting for ID: %s" % args.ID)

    args.coords = SkyCoord(args.ra * u.deg, args.dec * u.deg)

    args.survey = Survey("hdr1")

    args.matched_sources = {}
    shots_of_interest = []

    count = 0

    # this radius applies to the inital shot search and requires a large
    # aperture for the wide FOV of VIRUS

    max_sep = 11.0 * u.arcminute

    args.log.info("Finding shots of interest")

    for i, coord in enumerate(args.survey.coords):
        dist = args.coords.separation(coord)
        sep_constraint = dist < max_sep
        shotid = args.survey.shotid[i]
        idx = np.where(sep_constraint)[0]
        if np.size(idx) > 0:
            args.matched_sources[shotid] = idx
            count += np.size(idx)
            if len(idx) > 0:
                shots_of_interest.append(shotid)

    args.log.info("Number of shots of interest: %i" % len(shots_of_interest))
    args.log.info("Saved shot list to file " + str(args.outfile))
    np.savetxt("shotlist", shots_of_interest, fmt="%i")


if __name__ == "__main__":
    main()
