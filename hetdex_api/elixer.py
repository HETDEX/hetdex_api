# -*- coding: utf-8 -*-
"""

Script to read in exlier catalogs and write
to machine readable format.

Created: 2019/03/28

Authors: Erin Mentuch Cooper, Dustin Davis

To make an astropy Table from elixer catalogs:

elixer_table = read_elixer_catalogs(fibfn, catfn)

Note: This is only used for HDR1
"""
from __future__ import print_function

import numpy as np

from astropy.table import Table


def read_elixer_catalogs(fibfn, catfn):
    """
    read the _fib.txt and _cat.txt catalogs
    create Detections for each (combining across and within catalogs as needed)
    :return: list of Detections
    """

    fib_version = None
    cat_version = None

    # should be one line per detection in fib.txt
    # at least one line but may be more in cat.txt
    elixer_table = Table()
    detectid = []
    detectidcat = []
    ra = []
    dec = []
    ew_obs = []
    plae_poii_hetdex = []
    plae_poii = []
    filt = []
    dist = []
    cat_mag = []
    ra_match = []
    dec_match = []

    with open(fibfn, "r") as f:
        for line in f:
            if line[0] == "#":  # except get # version 1.5.0a16
                if (fib_version is None) and ("version" in line):
                    toks = line.split()
                    if (toks is not None) and (len(toks) == 3):
                        fib_version = toks[2]
                continue

            if len(line) < 100:
                continue

            toks = line.split()

            if len(toks) < 29:  # must be at least 29 for _fib.txt
                continue

            detectid.append(np.int64(toks[1]))
            ra.append(float(toks[3]))
            dec.append(float(toks[4]))
            ew_obs.append(float(toks[13]))
            plae_poii_hetdex.append(float(toks[14]))

    data = [
        np.array(detectid),
        np.array(ra),
        np.array(dec),
        np.array(ew_obs),
        np.array(plae_poii_hetdex),
    ]
    colnames = ["detectid", "ra", "dec", "ew_obs", "plae_poii_hetdex"]
    elixer_table = Table(data, names=colnames)

    with open(catfn, "r") as f:
        for line in f:
            if line[0] == "#":
                if (cat_version is None) and ("version" in line):
                    toks = line.split()
                    if (toks is not None) and (len(toks) == 3):
                        cat_version = toks[2]
                continue

            if len(line) < 100:
                continue

            toks = line.split()

            if len(toks) < 29:
                continue

            detectidcat.append(np.int64(toks[1]))
            ra_match.append(float(toks[24]))
            dec_match.append(float(toks[25]))
            dist.append(float(toks[25]))
            cat_mag.append(toks[28])
            filt.append(toks[29])
            if (toks[32] is not None) and (toks[32].lower() != "none"):
                plae_poii.append(float(toks[32]))
            else:
                plae_poii.append(-1)

    colnames2 = [
        "detectidcat",
        "ra_match",
        "dec_match",
        "dist_match",
        "filt" "cat_mag",
        "plae_poii",
    ]
    data2 = [
        np.array(ra_match),
        np.array(dec_match),
        np.array(dist),
        np.array(filt),
        np.array(cat_mag),
        np.array(plae_poii),
    ]
    cat_table = Table(data2, names=colnames2)

    for detectid_i in enumerate(np.array(detectid)):
        idx = np.where(
            (cat_table["detectidcat"] == detectid_i)
            & (cat_table["ra_match"] == 666)
            & (cat_table["dec_match"] == 666)
        )
        print(idx)

    # Now we have consumed the catalog files
    return elixer_table, cat_table
