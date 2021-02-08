import sys
import os
import os.path as op
import subprocess

import numpy as np
import pandas as pd

from astropy.table import Table, Column, vstack, join, hstack

from hetdex_api.config import HDRconfig

survey = "hdr2.1"

config = HDRconfig(survey)

update = True
create = False


if create:
    # I start with Karl's calibration summary file
    #! cp /data/00115/gebhardt/skysub/analysis/all.hdr2.1 all.hdr2.1
    #! awk '{printf "%20s %14s %8.3f %5i %8.3f %8.3f %8.3f %8.3f %4i \n", $1, $2, $3, $4, $5, $6, $7, $8, $9}' all.hdr2.1 > all.hdr2.1_emc
    #! sed -i s/'s'/''/ all.hdr2.1_emc
    #! sed -i s/'d'/''/ all.hdr2.1_emc
    #! sed -i s/'exp0'/' '/ all.hdr2.1_emc

    amp_stats = Table.read(
        "all.hdr2.1_emc",
        format="ascii.no_header",
        fast_reader=False,
        guess=False,
        names=[
            "shotid",
            "expnum",
            "multiframe",
            "norm",
            "N_cont",
            "background",
            "sky_sub_rms",
            "wave0",
            "wave1",
            "nfib_bad",
        ],
    )

    multiname = []
    for row in amp_stats:
        multiname.append("multi_" + row["multiframe"])

    amp_stats["multiframe"] = multiname

    shotlist = np.unique(amp_stats["shotid"])
    median_sky_rms = []
    stdev_sky_rms = []

    for shoti in shotlist:
        sel = amp_stats["shotid"] == shoti
        median_sky_rms.append(np.nanmedian(amp_stats["sky_sub_rms"][sel]))
        stdev_sky_rms.append(np.nanstd(amp_stats["sky_sub_rms"][sel]))

    date_stats = pd.to_datetime((shotlist / 1000.0).astype(int).astype(str))

    shot_stats = pd.DataFrame(
        np.transpose(np.array([date_stats, shotlist, median_sky_rms, stdev_sky_rms])),
        columns=["Date", "shotid", "median_sky_rms", "stdev_sky_rms"],
    )


if update:

    amp_stats_fn = op.join(
        config.hdr_dir[ survey], "survey", "amp_stats_{}.tab".format(survey)
    )
    amp_stats = Table.read(amp_stats_fn, format="ascii")

    # remove time dependent amps found manually
    badamps = Table.read(config.badamp2, format="ascii")

    amp_stats["date"] = (amp_stats["shotid"] / 1000).astype(int)
    amp_stats["flag_manual"] = np.ones(np.size(amp_stats), dtype=int)

    amps_removed = 0

    for row in badamps:
        selmf = amp_stats["multiframe"] == row["multiframe"]
        seldate = (amp_stats["date"] >= row["date_start"]) * (
            amp_stats["date"] <= row["date_end"]
        )
        amps_removed += np.sum(selmf * seldate)
        amp_stats["flag_manual"][seldate * selmf] = 0

    # remove amps flagged in just a single observation
    badamps_single = Table.read(config.badamp_single, format="ascii")

    amps_removed = 0
    for row in badamps_single:
        selmf = amp_stats["multiframe"] == row["multiframe"]
        selshot = amp_stats["shotid"] == row["shotid"]
        amps_removed += np.sum(selmf * selshot)
        amp_stats["flag_manual"][selshot * selmf] = 0

    # remove all amps in a bad shot
    badshots = Table.read(config.badshot, names=["shotid"], format="ascii.no_header")

    amps_removed = 0
    for row in badshots:
        selshot = amp_stats["shotid"] == row["shotid"]
        amps_removed += np.sum(selshot)
        amp_stats["flag_manual"][selshot] = 0

    sel1 = (amp_stats["background"].astype(float) > -10) * (
        amp_stats["background"].astype(float) < 100
    )
    sel2 = amp_stats["sky_sub_rms_rel"] < 1.5
    sel3 = amp_stats["sky_sub_rms"] > 0.2
    sel4 = (amp_stats["im_median"] > 0.05) | (np.isnan(amp_stats["im_median"]))
    sel5 = (amp_stats["MaskFraction"] < 0.5) | (np.isnan(amp_stats["MaskFraction"]))
    sel6 = amp_stats["N_cont"] < 35
    sel7 = amp_stats["nfib_bad"] <= 1
    sel8 = amp_stats["norm"] > 0.5

    sel = sel1 * sel2 * sel3 * sel4 * sel5 * sel6 * sel7 * sel8

    sel_wiggles = amp_stats["ft_flag"] == 1
    sel_manual = amp_stats["flag_manual"] == 1
    amp_stats["flag"] = (sel * sel_manual * sel_wiggles).astype(int)

    amp_stats_by_exp = amp_stats.group_by(["shotid", "multiframe"])
    amp_stats_tab = amp_stats_by_exp.groups.aggregate(np.min)
    amp_stats_tab["shotid", "multiframe", "flag"].write(
        "amp_flag.tab", format="ascii", overwrite=True
    )
    amp_stats_tab["shotid", "multiframe", "flag"].write("amp_flag.fits", overwrite=True)

    amp_stats.write('amp_stats_{}.tab'.format(survey), format='ascii', overwrite=True)
    amp_stats.write('amp_stats_{}.fits'.format(survey), overwrite=True)
