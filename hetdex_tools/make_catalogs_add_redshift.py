#!/usr/bin/env python
# coding: utf-8

import sys
import numpy as np
import os.path as op
import time
import pandas as pd
import tables as tb
import argparse as ap


from astropy.table import Table, unique, join, Column, vstack
from astropy.cosmology import Planck18 as cosmo
import astropy.units as u

from hetdex_api.config import HDRconfig
import hetdex_tools.fof_kdtree as fof

# import hetdex_api.wave as wv

from multiprocessing import Pool


def add_z_guess(source_id):
    global source_table, agn_tab, config

    #    if agn_tab is None:
    #        agn_tab = Table.read(config.agncat, format="ascii")

    sel_group = source_table["source_id"] == source_id

    group = source_table[sel_group]
    z_guess = -1.0
    s_type = "none"
    agn_flag = -1
    z_conf = -1
    z_src = ""

    # print(group['best_z', 'best_pz','detectid'])
    if True:
        # Check if any member is an AGN first
        for det in group["detectid"]:
            if agn_tab is not None:
                if det in agn_tab["detectid"]:
                    agn_flag = 1
                    agn_det = det

        if agn_flag == 1:
            # get proper z's from Chenxu's catalog
            sel_det = agn_tab["detectid"] == agn_det
            z_guess = agn_tab["z"][sel_det][0]
            # agn_flag = agn_tab['zflag'][sel_det][0]
            z_conf = 1  # agn_tab['zflag'][sel_det][0]
            s_type = "agn"
            z_src = "liu_agn"

        elif np.any((group["cls_diagnose"] == "STAR") & (group["gmag"] < 22)):
            s_type = "star"
            sel_best = (group["cls_diagnose"] == "STAR") & (group["gmag"] < 22)
            s_type = "star"
            closest = np.argmin(group["src_separation"][sel_best])
            z_guess = group["z_diagnose"][sel_best][closest]
            z_conf = 0.9
            z_src = "diagnose"

        # elif np.any((group["cls_diagnose"] == "QSO") & (group["gmag"] < 22)):
        #    sel_best = (group["cls_diagnose"] == "QSO") & (group["gmag"] < 22)
        #    s_type = "agn"
        #    closest = np.argmin(group["src_separation"][sel_best])
        #    z_guess = group["z_diagnose"][sel_best][closest]
        #    z_conf = 0.9
        #    z_src = "diagnose"

        elif np.any((group["cls_diagnose"] == "GALAXY") & (group["gmag"] < 22)):
            sel_best = (group["cls_diagnose"] == "GALAXY") & (group["gmag"] < 22)
            closest = np.argmin(group["src_separation"][sel_best])
            brightest = np.argmin(group["gmag"][sel_best])
            z_guess = group["z_diagnose"][sel_best][brightest]
            z_conf = 0.9
            z_src = "diagnose"
            if np.any(np.abs(group["wave"] - waveoii * (1 + z_guess)) < 10):
                s_type = "oii"
            else:
                s_type = "lzg"

        elif np.any(
            (group["cls_diagnose"] == "GALAXY") & (group["plya_classification"] < 0.5)
        ):
            sel_best = (group["cls_diagnose"] == "GALAXY") & (
                group["plya_classification"] < 0.5
            )
            closest = np.argmin(group["src_separation"][sel_best])
            z_guess = group["z_diagnose"][sel_best][closest]
            z_conf = 0.9
            z_src = "diagnose"
            if np.any(np.abs(group["wave"] - waveoii * (1 + z_guess)) < 10):
                s_type = "oii"
            else:
                s_type = "lzg"

        elif np.any(
            (group["cls_diagnose"] == "GALAXY") * np.isnan(group["plya_classification"])
        ):
            sel_best = group["cls_diagnose"] == "GALAXY"
            closest = np.argmin(group["src_separation"][sel_best])
            z_guess = group["z_diagnose"][closest]
            z_conf = 0.9
            z_src = "diagnose"

            if np.any(np.abs(group["wave"] - waveoii * (1 + z_guess)) < 10):
                s_type = "oii"
            else:
                s_type = "lzg"
        else:
            if np.size(group) == 1:
                if np.isfinite(group["z_elixer"]):
                    z_guess = float(group["z_elixer"])
                    z_conf = float(group["best_pz"])

                    if 0.0 < z_guess < 0.5:
                        s_type = "oii"
                    elif 1.8 < z_guess < 3.7:
                        s_type = "lae"
                    elif z_guess < 0:
                        if float(group["wave"]) < waveoii:
                            # assume this should actually be Lya
                            z_guess = float(group["wave"]) / wavelya - 1
                            s_type = "lae"
                    elif 0.5 < z_guess < 1.9:
                        if float(group["plya_classification"]) >= 0.5:
                            z_guess = float(group["wave"]) / wavelya - 1
                            s_type = "lae"
                        else:
                            z_guess = float(group["wave"]) / waveoii - 1
                            s_type = "oii"
                else:
                    # use plae_classification
                    if group["plya_classification"] >= 0.5:
                        z_guess = float(group["wave"]) / wavelya - 1
                        s_type = "lae"
                    else:
                        z_guess = float(group["wave"]) / waveoii - 1
                        s_type = "oii"
                z_src = "elixer"
            else:
                if np.any(group["wave_group_id"] > 0):
                    sel_wave_group = group["wave_group_id"] > 0
                    sel_most_conf = np.argmin(group["src_separation"][sel_wave_group])
                    z_guess = group["z_elixer"][sel_wave_group][sel_most_conf]
                    z_conf = group["best_pz"][sel_wave_group][sel_most_conf]
                    # assign s_type
                    if 0.0 < z_guess < 0.5:
                        s_type = "oii"
                    elif 1.88 < z_guess < 3.6:
                        s_type = "lae"
                    else:
                        s_type = "gal"

                    z_src = "elixer"

                elif np.std(group["z_elixer"]) < 0.02:
                    # print(np.array(group['detectid']))
                    sel_most_conf = np.argmin(group["src_separation"])
                    z_guess = group["z_elixer"][sel_most_conf]
                    z_conf = group["best_pz"][sel_most_conf]
                    if 0.0 < z_guess < 0.5:
                        s_type = "oii"
                    elif 1.88 < z_guess < 3.6:
                        s_type = "lae"
                    else:
                        s_type = "gal"

                    z_src = "elixer"

                else:
                    # assign redshift of highest best_pz line
                    # or use highest S/N line
                    if np.any(group["best_pz"] > 0.6):
                        sel_sn = np.argmax(group["best_pz"])
                    else:
                        sel_sn = np.argmax(group["sn"])

                    z_guess = group["z_elixer"][sel_sn]
                    z_conf = group["best_pz"][sel_sn]
                    z_src = "elixer"
                    if 0.0 < z_guess < 0.5:
                        s_type = "oii"
                    elif 1.88 < z_guess < 3.6:
                        z_guess = -2
                        s_type = "lae"
                    else:
                        s_type = "gal"

                # check to see if its a Lya/CIV match, assume blue line is Lya
                #                try:
                #                    sel_sn = group['sn'] > 5.5
                #                    wavelya_obs = np.min( group['wave'][sel_sn] )
                #                except:
                #                    wavelya_obs = np.min( group['wave'] )
                #
                #                    zlya_obs = wavelya_obs/wavelya - 1
                sel_line = group["det_type"] == "line"
                zlya_obs = z_guess
                if np.any(
                    np.abs(group["wave"][sel_line] - waveciv * (1 + zlya_obs)) < 20
                ):
                    print("lya, civ match for", list(group["detectid"][sel_line]))
                    z_guess = zlya_obs
                    z_conf = 0.95
                    s_type = "agn"
                    agn_flag = 1
                    z_src = "2em"
                elif np.any(
                    np.abs(group["wave"][sel_line] - waveoii * (1 + z_guess)) < 20
                ):
                    s_type = "oii"
    else:  # except:
        print("could not classify {}".format(source_id))

    if np.isnan(z_guess):
        z_guess = -3
        z_conf = 0
        s_type = "none"
        z_src = ""

    return z_guess, z_conf, s_type, agn_flag, z_src


def zcluster_forshotid(shotid, star=False):
    global source_table

    sel = source_table["shotid"] == shotid

    uniq_table = unique(source_table[sel], keys=["source_id"])

    kdtree, r = fof.mktree(
        np.array(uniq_table["ra_mean"]),
        np.array(uniq_table["dec_mean"]),
        np.array(uniq_table["z_hetdex"]),
        dsky=6,
        dwave=0.01,
    )

    zfriend_lst = fof.frinds_of_friends(kdtree, r, Nmin=2)

    if len(zfriend_lst) == 0:
        return None

    zfriend_table = fof.process_group_list(
        zfriend_lst,
        np.array(uniq_table["source_id"]),
        np.array(uniq_table["ra_mean"]),
        np.array(uniq_table["dec_mean"]),
        np.array(uniq_table["z_hetdex"]),
        np.array(uniq_table["flux_g"]),
    )

    memberlist = []
    friendlist = []

    for row in zfriend_table:
        friendid = row["id"]
        members = np.array(row["members"])
    friendlist.extend(friendid * np.ones_like(members))
    memberlist.extend(members)
    zfriend_table.remove_column("members")

    zdetfriend_tab = Table()
    zdetfriend_tab.add_column(Column(np.array(friendlist), name="id"))
    zdetfriend_tab.add_column(Column(memberlist, name="source_id"))
    zdetfriend_all = join(zdetfriend_tab, zfriend_table, keys="id")
    return zdetfriend_all

def zclusteragn_forshotid(shotid):
    global source_table

    sel = source_table["shotid"] == shotid

    uniq_table = unique(source_table[sel], keys=["source_id"])

    kdtree, r = fof.mktree(
        np.array(uniq_table["ra_mean"]),
        np.array(uniq_table["dec_mean"]),
        np.array(uniq_table["z_hetdex"]),
        dsky=10,
        dwave=0.1,
    )

    zfriend_lst = fof.frinds_of_friends(kdtree, r, Nmin=2)

    if len(zfriend_lst) == 0:
        return None

    zfriend_table = fof.process_group_list(
        zfriend_lst,
        np.array(uniq_table["source_id"]),
        np.array(uniq_table["ra_mean"]),
        np.array(uniq_table["dec_mean"]),
        np.array(uniq_table["z_hetdex"]),
        np.array(uniq_table["flux_g"]),
    )

    memberlist = []
    friendlist = []
    for row in zfriend_table:
        friendid = row["id"]
        members = np.array(row["members"])
        friendlist.extend(friendid * np.ones_like(members))
        memberlist.extend(members)

    zfriend_table.remove_column("members")

    zdetfriend_tab = Table()
    zdetfriend_tab.add_column(Column(np.array(friendlist), name="id"))
    zdetfriend_tab.add_column(Column(memberlist, name="source_id"))
    zdetfriend_all = join(zdetfriend_tab, zfriend_table, keys="id")
    return zdetfriend_all

def zclusterstars_forshotid(shotid):
    global source_table

    sel = source_table["shotid"] == shotid

    uniq_table = unique(source_table[sel], keys=["source_id"])

    kdtree, r = fof.mktree(
        np.array(uniq_table["ra_mean"]),
        np.array(uniq_table["dec_mean"]),
        np.array(uniq_table["z_hetdex"]),
        dsky=30,
        dwave=0.0001,
    )

    zfriend_lst = fof.frinds_of_friends(kdtree, r, Nmin=2)

    if len(zfriend_lst) == 0:
        return None

    zfriend_table = fof.process_group_list(
        zfriend_lst,
        np.array(uniq_table["source_id"]),
        np.array(uniq_table["ra_mean"]),
        np.array(uniq_table["dec_mean"]),
        np.array(uniq_table["z_hetdex"]),
        np.array(uniq_table["flux_g"]),
    )

    memberlist = []
    friendlist = []
    for row in zfriend_table:
        friendid = row["id"]
        members = np.array(row["members"])
        friendlist.extend(friendid * np.ones_like(members))
        memberlist.extend(members)

    zfriend_table.remove_column("members")

    zdetfriend_tab = Table()
    zdetfriend_tab.add_column(Column(np.array(friendlist), name="id"))
    zdetfriend_tab.add_column(Column(memberlist, name="source_id"))
    zdetfriend_all = join(zdetfriend_tab, zfriend_table, keys="id")
    return zdetfriend_all


def update_table(zdetfriend_all):
    global source_table

    for zid in np.unique(zdetfriend_all["id"]):
        sel_zid = zdetfriend_all["id"] == zid
        sids = zdetfriend_all["source_id"][sel_zid]

        sid_main = np.min(sids)

        for sid_i in sids:
            sel_sid = source_table["source_id"] == sid_i
            source_table["source_id"][sel_sid] = sid_main

        sid_ind = np.where(source_table["source_id"] == sid_main)

        res_tab = fof.process_group_list(
            sid_ind,
            source_table["detectid"],
            source_table["ra"],
            source_table["dec"],
            source_table["z_hetdex"],
            source_table["flux_g"],
        )

        res_tab.rename_column("size", "n_members")
        res_tab.rename_column("icx", "ra_mean")
        res_tab.rename_column("icy", "dec_mean")

        for col in res_tab.colnames:
            if col in ["id", "members"]:
                continue
                source_table[col][sid_ind] = res_tab[col]

            if res_tab["izz"] > 0.001:
                print("big redshift error on {}".format(sid_main))
                print(res_tab["icz"], res_tab["izz"])
            else:
                source_table["z_hetdex"][sid_ind] = res_tab["icz"]

        # update source type assignment
        stype_main = source_table["source_type"][source_table["source_id"] == sid_main]
        source_table["source_type"][sid_ind] = stype_main
    return


# Enter the catalog version


def get_parser():
    """function that returns a parser from argparse"""

    parser = ap.ArgumentParser(
        description="""Extracts 1D spectrum at specified RA/DEC""", add_help=True
    )

    parser.add_argument(
        "-v",
        "--version",
        help="""source catalog version you want to create""",
        type=str,
        default=None,
    )

    parser.add_argument(
        "-a",
        "--add_redshifts",
        help="""boolean trigger to add redshifts and do first part of good""",
        default=False,
        required=False,
        action="store_true",
    )
    return parser


def main(argv=None):
    """Main Function"""

    parser = get_parser()
    args = parser.parse_args(argv)
    version = args.version

    global wavelya, waveoii, waveciv, waveheii

    wavelya = 1215.67  # in vacuum
    waveoii = 3727.8  # in air
    waveciv = 1549.5
    waveheii = 1640.4

    global source_table, uniq_table, config, agn_tab

    print("Adding redshifts to catalog: {}".format(version))

    config = HDRconfig("hdr4")
    agn_tab = Table.read(config.agncat, format="ascii")
    # agn_tab = None

    if args.add_redshifts:
        # catfile = op.join(config.detect_dir, 'catalogs', 'source_catalog_' + version + '.fits')

        catfile = "source_catalog_{}.fits".format(version)

        source_table = Table.read(catfile)

        print("Source catalog was found at {}".format(catfile))

        # match Diagnose classification table

        try:
            print("Removing existing redshift and classification info")
            for col in [
                "z_hetdex",
                "z_hetdex_src",
                "z_hetdex_conf",
                "source_type",
                "agn_flag",
                "z_diagnose",
                "cls_diagnose",
                "stellartype",
                "selected_det",
            ]:
                source_table.remove_column(col)
        except:
            print("Could not remove existing redshift and classification info")
            pass

        diagnose_tab_hdr3 = Table.read(
            "/work/05350/ecooper/stampede2/redshift-tests/hdr3.0.3_lt23/diagnose_3.0.3_lt23.fits"
        )
        diagnose_tab_hdr4 = Table.read(
            "/work/05350/ecooper/stampede2/redshift-tests/hdr4.0.0_lt23/diagnose_4.0.0_lt23.fits"
        )
        diagnose_tab_hdr5 = Table.read(
            "/work/05350/ecooper/stampede2/hdr5/diagnose/diagnose_5.0.0_lt23.fits"
        )
        diagnose_tab_sa22 = Table.read(
            "/work/05350/ecooper/stampede2/redshift-tests/sa22/diagnose_sa22_lt23.fits"
            )
        diagnose_tab = unique( vstack([diagnose_tab_sa22, diagnose_tab_hdr3, diagnose_tab_hdr4, diagnose_tab_hdr5]), keys='detectid')

        diagnose_tab.rename_column("z_best", "z_diagnose")
        diagnose_tab.rename_column("classification", "cls_diagnose")

        combined = join(
            source_table,
            diagnose_tab["detectid", "z_diagnose", "cls_diagnose", "stellartype"],
            join_type="left",
            keys=["detectid"],
        )
        source_table = combined.copy()

        del combined, diagnose_tab

        # assign redshifts to sources with single detections
        sel = source_table["n_members"] == 1

        # assign AGN redshifts in AGN catalog

        if agn_tab is not None:
            sel = (source_table["n_members"] == 1) & (source_table["z_agn"] > 0)

            agn_assign = source_table["source_id", "z_agn"][sel].copy()
            agn_assign.rename_column("z_agn", "z_hetdex")
            agn_assign["z_hetdex_conf"] = 0.9
            agn_assign["z_hetdex_src"] = "liu_agn"
            agn_assign["source_type"] = "agn"
            agn_assign["agn_flag"] = 1

            sel = (source_table["n_members"] == 1) & np.invert(
                source_table["z_agn"] > 0
            )

        # Take Diagnose for gmag < 22 and continuum for remaining single detection sources
        sel_gmag = (
            (source_table["gmag"] < 22) | (source_table["det_type"] == "cont")
        ) & (source_table["cls_diagnose"] != "UNKNOWN")

        diagnose_assign = source_table["source_id", "z_diagnose", "cls_diagnose"][
            sel & sel_gmag
        ].copy()
        diagnose_assign.rename_column("z_diagnose", "z_hetdex")
        diagnose_assign["source_type"] = "          "
        diagnose_assign["source_type"][
            diagnose_assign["cls_diagnose"] == "STAR"
        ] = "star"
        diagnose_assign["source_type"][diagnose_assign["cls_diagnose"] == "QSO"] = "agn"
        diagnose_assign["source_type"][
            diagnose_assign["cls_diagnose"] == "GALAXY"
        ] = "lzg"
        diagnose_assign["z_hetdex_src"] = "diagnose"
        diagnose_assign["z_hetdex_conf"] = 0.8
        diagnose_assign["agn_flag"] = -1
        diagnose_assign["agn_flag"][diagnose_assign["cls_diagnose"] == "QSO"] = 1
        diagnose_assign.remove_column("cls_diagnose")

        # Take Elixer best_z for rest

        elixer_assign = source_table["source_id", "z_elixer", "best_pz"][
            sel & np.invert(sel_gmag) & (source_table['det_type'] == 'line')
        ]
        elixer_assign.rename_column("z_elixer", "z_hetdex")
        elixer_assign.rename_column("best_pz", "z_hetdex_conf")
        elixer_assign["z_hetdex_src"] = "elixer"
        elixer_assign["agn_flag"] = -1
        sel_lae = elixer_assign["z_hetdex"] > 1.88
        sel_oii = elixer_assign["z_hetdex"] < 0.5
        elixer_assign["source_type"] = "gal"
        elixer_assign["source_type"][sel_lae] = "lae"
        elixer_assign["source_type"][sel_oii] = "oii"

        sel = source_table["n_members"] == 1

        uniq_table = unique(source_table[np.invert(sel)], keys="source_id")

        src_list = uniq_table["source_id"]

        ntask = 20
        print("Adding z_hetdex using {} cores".format(ntask))
        t0 = time.time()
        p = Pool(ntask)
        res = p.map(add_z_guess, uniq_table["source_id"])
        p.close()
        t1 = time.time()

        print("Adding z_hetdex complete in {:5.3} m".format((t1 - t0) / 60))

        z_hetdex = []
        z_conf = []
        s_type = []
        agn_flag = []
        z_src = []
        for r in res:
            z_hetdex.append(float(r[0]))
            z_conf.append(float(r[1]))
            s_type.append((r[2]))
            agn_flag.append(float(r[3]))
            z_src.append(r[4])

        z_table = Table(
            [uniq_table["source_id"], z_hetdex, z_src, z_conf, s_type, agn_flag],
            names=[
                "source_id",
                "z_hetdex",
                "z_hetdex_src",
                "z_hetdex_conf",
                "source_type",
                "agn_flag",
            ],
        )

        if agn_tab is not None:
            z_stack = vstack([agn_assign, diagnose_assign, elixer_assign, z_table])
        else:
            z_stack = vstack([diagnose_assign, elixer_assign, z_table])

        z_stack.write("z_stack_{}.fits".format(version), overwrite=True)
        source_table.write("source_catalog_{}.y.fits".format(version), overwrite=True)
        sys.exit()
    else:
        catfile = "source_catalog_{}.y.fits".format(version)

        source_table = Table.read(catfile)
        z_stack = Table.read("z_stack_{}.fits".format(version))
        print("Source catalog was found at {}".format(catfile))

    all_table = join(source_table, z_stack, join_type="left")

    source_table = all_table

    # uncluster LAEs whose line pairings are not supported by best_z

    clustered_lae_index = np.where(source_table["z_hetdex"] == -2)[0]

    print(
        "Updating z_hetdex to best_z for {} detectids".format(len(clustered_lae_index))
    )
    
    vs = version.split('.')
    sid_index = int( vs[0]) * 10**13 + int(vs[1])*10**12 + int(vs[2])*10**11 + 3*10**9
    #sid_index = 4010030000000

    for c_ind in clustered_lae_index:
        source_table["source_id"][c_ind] = sid_index
        source_table["z_hetdex"][c_ind] = source_table["z_elixer"][c_ind]
        source_table["z_hetdex_conf"][c_ind] = source_table["best_pz"][c_ind]
        source_table["z_hetdex_src"][c_ind] = "elixer"
        sid_index += 1

    sel_star = source_table["source_type"] == "star"
    sel_oii = source_table["source_type"] == "oii"
    sel_lae = source_table["source_type"] == "lae"
    sel_agn = source_table["source_type"] == "agn"
    sel_lzg = source_table["source_type"] == "lzg"

    # assign Line ID
    spectral_lines = Table.read(
        "/work2/05350/ecooper/stampede2/redshift-tests/spectral_lines.txt",
        format="ascii",
    )

    source_table["line_id"] = np.chararray(len(source_table), 20, unicode=True)
    source_table['line_id'] = 'n/a'
    
    source_table["dwave"] = np.zeros_like(source_table["wave"])

    for row in spectral_lines:
        sel = (
            np.abs(
                source_table["wave"]
                - (1 + source_table["z_hetdex"]) * row["lambda_vac"]
            )
            < 15
        )
        source_table["line_id"][sel] = str(row["species"])

        # allow AGN to have broader fit
        sel = (
            np.abs(
                source_table["wave"]
                - (1 + source_table["z_hetdex"]) * row["lambda_vac"]
            )
            < 30
        )
        source_table["line_id"][sel & sel_agn] = str(row["species"])

    print(
        "There are {} stars, {} OII emitters, {} LAEs and {} Low-z gals".format(
            np.sum(sel_star), np.sum(sel_oii), np.sum(sel_lae), np.sum(sel_lzg)
        )
    )

    # sort table by source wavelength closest to z_hetdex and position
    src_wave = np.zeros_like(source_table["z_hetdex"])
    src_wave[sel_oii] = (1 + source_table["z_hetdex"][sel_oii]) * waveoii
    src_wave[sel_lae] = (1 + source_table["z_hetdex"][sel_lae]) * wavelya

    sel_z = (source_table["z_hetdex"] >= 1.9) * (source_table["z_hetdex"] <= 3.6)
    src_wave[sel_agn * sel_z] = (
        1 + source_table["z_hetdex"][sel_agn * sel_z]
    ) * wavelya

    sel = sel_oii | sel_lae | sel_agn
    source_table["dwave"][sel] = np.abs(src_wave[sel] - source_table["wave"][sel])

    # Cluster in redshift space to group those missed by wave clustering

    print("Clustering in redshift space")
    shotid_list = np.unique(source_table["shotid"])

    t0 = time.time()
    p = Pool(16)
    res = p.map(zcluster_forshotid, shotid_list)
    p.close()

    for r in res:
        if r is not None:
            update_table(r)

    print("Clustering stars in redshift space")

    shotid_list = np.unique(source_table["shotid"])

    p = Pool(16)
    res = p.map(zclusterstars_forshotid, shotid_list)
    p.close()

    for r in res:
        if r is not None:
            update_table(r)
    t1 = time.time()

    print("Done clustering in redshift space in {:5.2f} minutes".format((t1 - t0) / 60))

    source_table.sort(["source_id", "gmag"])

    # change lzgs to oii if an OII line is found
    sel_oii_line = source_table["line_id"] == "OII"
    sel_lzg = source_table["source_type"] == "lzg"

    sids_to_switch = np.unique(source_table["source_id"][sel_oii_line & sel_lzg])

    for sid in sids_to_switch:
        sel = source_table["source_id"] == sid
        source_table["source_type"][sel] = "oii"

    #update OII line list
    sel_oii_line = source_table["line_id"] == "OII"
    
    sel_star = source_table["source_type"] == "star"
    sel_oii = source_table["source_type"] == "oii"
    sel_lae = source_table["source_type"] == "lae"
    sel_agn = source_table["source_type"] == "agn"
    sel_lzg = source_table["source_type"] == "lzg"

    print(
        "There are {} stars, {} OII emitters, {} LAEs, {} AGN and {} Low-z detections".format(
            np.sum(sel_star),
            np.sum(sel_oii),
            np.sum(sel_lae),
            np.sum(sel_agn),
            np.sum(sel_lzg),
        )
    )

    # now assigned selected det flag to LAE sample
    # for agn, stars, lzg and oii galaxies will use detectid closest to FOF center

    sel_lae_line = (source_table["line_id"] == "Lya") & (
        source_table["source_type"] == "lae"
    )

    lae_tab = source_table[sel_lae_line]
    lae_tab.sort("flux", reverse=True)

    uniq_obs_lae = unique(lae_tab, keys="source_id")
    uniq_obs_lae["selected_det"] = True

    oii_tab = source_table[sel_oii_line]
    oii_tab.sort("flux", reverse=True)

    uniq_obs_oii = unique(oii_tab, keys="source_id")
    uniq_obs_oii["selected_det"] = True

    sel_agn_with_lya = (source_table["line_id"] == "Lya") & (
        source_table["source_type"] == "agn"
    )

    agn_with_lya = source_table[sel_agn_with_lya]
    agn_with_lya.sort("flux", reverse=True)

    uniq_agn_with_lya = unique(agn_with_lya, keys="source_id")
    uniq_agn_with_lya["selected_det"] = True

    sel_agn = source_table["source_type"] == "agn"

    # find source_ids not in above list

    agn_without_lya_sid = np.setdiff1d(
        np.unique(source_table["source_id"][sel_agn]),
        uniq_agn_with_lya["source_id"],
        assume_unique=True,
    )

    agn_without_lya = join(
        Table([agn_without_lya_sid], names=["source_id"]),
        source_table,
        join_type="left",
    )
    agn_without_lya.sort("gmag")

    uniq_agn_without_lya = unique(agn_without_lya, keys="source_id")

    uniq_agn_without_lya["selected_det"] = True

    # now assign brightest detectid to star, lzg groups
    sel_rest = (
        (source_table["source_type"] != "lae")
        & (source_table["source_type"] != "oii")
        & (source_table["source_type"] != "agn")
    )

    # sort the rest of catalog using the brightest gmag 
    rest_tab = source_table[sel_rest]
    rest_tab.sort("gmag")

    uniq_rest = unique(rest_tab, keys=["source_id"])
    uniq_rest["selected_det"] = True

    selected_table = vstack(
        [uniq_obs_lae, uniq_obs_oii, uniq_agn_with_lya, uniq_agn_without_lya, uniq_rest]
    )

    source_table2 = join(
        source_table, selected_table["detectid", "selected_det"], join_type="left"
    )

    source_table = source_table2.copy()

    del source_table2, oii_tab, lae_tab, rest_tab
    del uniq_obs_oii, uniq_obs_lae, uniq_rest

    source_table["selected_det"] = source_table["selected_det"].filled(False)

    # add luminosities

    source_table.add_column(Column(name="lum_lya", length=len(source_table)))
    source_table.add_column(Column(name="lum_lya_err", length=len(source_table)))
    source_table.add_column(Column(name="lum_oii", length=len(source_table)))
    source_table.add_column(Column(name="lum_oii_err", length=len(source_table)))

    sel_oii = (source_table["selected_det"] == True) & (
        source_table["source_type"] == "oii"
    )

    lum_dist = cosmo.luminosity_distance(source_table["z_hetdex"][sel_oii]).to(u.cm)
    fac = 10 ** (-17) * 4.0 * np.pi * lum_dist**2

    source_table["lum_oii"][sel_oii] = source_table["flux"][sel_oii] * fac
    source_table["lum_oii_err"][sel_oii] = source_table["flux_err"][sel_oii] * fac

    print("Number of OII assigned detections")
    # print(np.sum(sel_oii), np.sum(sel_oii_aper), np.sum( (source_table['selected_det'] == True) & (source_table['source_type'] == 'oii')))
    print(
        np.sum(sel_oii),
        np.sum(
            (source_table["selected_det"] == True)
            & (source_table["source_type"] == "oii")
        ),
    )
    sel_lae = (source_table["selected_det"] == True) & (
        source_table["source_type"] == "lae"
    )

    lum_dist = cosmo.luminosity_distance(source_table["z_hetdex"][sel_lae]).to(u.cm)
    fac = 10 ** (-17) * 4.0 * np.pi * lum_dist**2

    source_table["lum_lya"][sel_lae] = source_table["flux"][sel_lae] * fac
    source_table["lum_lya_err"][sel_lae] = source_table["flux_err"][sel_lae] * fac

    sel_agn_lya = (
        (source_table["selected_det"] == True)
        & (source_table["source_type"] == "agn")
        & (source_table["line_id"] == "Lya")
    )

    lum_dist = cosmo.luminosity_distance(source_table["z_hetdex"][sel_agn_lya]).to(u.cm)
    fac = 10 ** (-17) * 4.0 * np.pi * lum_dist**2

    source_table["lum_lya"][sel_agn_lya] = source_table["flux"][sel_agn_lya] * fac
    source_table["lum_lya_err"][sel_agn_lya] = (
        source_table["flux_err"][sel_agn_lya] * fac
    )

    source_table.sort("source_id")

    source_table["source_type"] = source_table["source_type"].astype(str)

    # add DEE classifications
    dee = unique(
        Table.read(
            "/scratch/projects/hetdex/hdr4/catalogs/ml/zoo_tsne_rfnumber_20231107.csv"
        ),
        keys="detectid",
    )
    dee["dee_prob"] = dee["rf_number"].filled(-1)
    sel_bad = (dee["tsne_x"] > 5) & (dee["tsne_y"] < 0)
    dee["flag_dee_tsne"] = np.invert(sel_bad).astype(int)
    dee["flag_dee"] = np.invert(
        sel_bad * (dee["dee_prob"] <= 0.3) * (dee["dee_prob"] >= 0.0)
    ).astype(int)

    combined = join(
        source_table,
        dee["detectid", "dee_prob", "tsne_x", "tsne_y", "flag_dee_tsne", "flag_dee"],
        join_type="left",
    )

    print(len(source_table), len(combined))

    # add 2D profile fits

    fit_2D = Table.read(
        "/work/05350/ecooper/stampede2/hdr4/fit_profile_output_4.0.0.fits"
    )
    fit_2D.rename_column('sn_max', 'sn_moffat')
    source_table2 = join(
        combined, fit_2D["detectid", "sn_im", "sn_moffat", "chi2_moffat"], join_type="left"
    )

    # add sn_elixer for p_real

    elixh5 = tb.open_file(config.elixerh5, 'r')
    elix_line_fit = Table( elixh5.root.SpectraLines.read())
    elix_line_fit.rename_column('sn','sn_elix')
    elixh5.close()
    
    src2 = join(source_table, elix_line_fit['detectid','wavelength', 'sn_elix'], keys='detectid')
    sel_wave_match = np.abs( src2['wavelength'] - src2['wave'] ) < 10

    elix_lines = src2['detectid','sn_elix'][sel_wave_match]

    source_table3 = join( source_table2, elix_lines, join_type='left')

    # fill in sn value for any mask sn_elix values
    for i in np.where( source_table3['sn_elix'].mask == True ):
        source_table3['sn_elix'][i] = source_table3['sn'][i]
        
    for i in np.where( np.isnan( source_table3['sn_elix'] )):
        source_table3['sn_elix'][i] = source_table3['sn'][i]
    
    print(len(combined), len(source_table2), len(source_table3))

    source_table = unique(source_table3, keys="detectid")
    # finalize flags

    source_table["dee_prob"] = source_table["dee_prob"].filled(-1)
    source_table["flag_dee"] = source_table["flag_dee"].filled(-1)
    source_table["flag_dee"][source_table["source_type"] == "agn"] = 1
    source_table["flag_dee_tsne"] = source_table["flag_dee_tsne"].filled(-1)
    source_table["flag_dee_tsne"][source_table["source_type"] == "agn"] = 1

    source_table["flag_lowlw"] = (
        (source_table["linewidth"] >= 1.7) | (source_table["det_type"] == "cont")
    ).astype(int)
    
    source_table["flag_apcor"] = np.invert(source_table["apcor"] < 0.45) & (
        source_table["det_type"] == "line"
    ).astype(int)

    # we can keep the AGN vetted detections with low apcor
    source_table["flag_apcor"][source_table["z_agn"] > 0] = 1
    # also keep all continuum sources
    source_table["flag_apcor"][source_table["det_type"] == "cont"] = 1

    source_table['flag_continuum'] = (source_table['continuum'] > -1).astype(int)
    
    source_table["flag_wave"] = (
        (source_table["wave"] >= 3510) * (source_table["wave"] <= 5496)
    ).astype(int)
    source_table["flag_wave"][source_table["det_type"] == "cont"] = 1

    source_table["flag_3540"] = (
        (source_table["wave"] < 3534)
        | (source_table["wave"] > 3556)  # changed 2024-10-28
    ).astype(int)

    source_table["flag_3540"][source_table["z_agn"] > 0] = 1

    source_table["flag_seldet"] = source_table["selected_det"].astype(int)
    source_table["flag_fwhm"] = (source_table["fwhm"] <= 2.66).astype(int)

    # add RAIC continuum data

    raic_cat_cont = Table.read(
        "/scratch/projects/hetdex/hdr4/catalogs/ml/raic_results_cont_labels_20231108_hdr4.0.0.fits"
    )

    # sort by score, will take unique label for highest score
    raic_cat_cont.sort("Score")
    raic_cat_cont.reverse()

    raic_cat_cont_uniq = unique(raic_cat_cont, keys="detectid")
    raic_cat_cont_uniq.rename_column("Score", "RAIC_Score")
    raic_cat_cont_uniq.rename_column("Label", "RAIC_Label")

    # add RAIC line data
    raic_cat = Table.read(
        "/scratch/projects/hetdex/hdr4/catalogs/ml/raic_results_line_labels_20231107_hdr4.0.0.fits"
    )

    raic_cat.sort("Score")
    raic_cat.reverse()

    raic_cat_uniq = unique(raic_cat, keys="detectid")
    raic_cat_uniq.rename_column("Score", "RAIC_Score")
    raic_cat_uniq.rename_column("Label", "RAIC_Label")

    raic_cat_stack = vstack([raic_cat_cont_uniq, raic_cat_uniq])

    source_table2 = join(
        source_table,
        raic_cat_stack["detectid", "RAIC_Score", "RAIC_Label"],
        join_type="left",
    )

    print(
        "Lengths before:after raic join: {}:{}".format(
            len(source_table), len(source_table2)
        )
    )

    source_table = source_table2

    source_table["raic_streak"] = (
        (
            (
                source_table["RAIC_Score"] > 0.35
            )  # note this applies to both continuum and line emission detections
            & (source_table["RAIC_Label"] == "streak")
        )
        .astype(int)
        .filled(-1)
    )

    source_table["raic_meteor"] = (
        ((source_table["RAIC_Score"] > 0.5) & (source_table["RAIC_Label"] == "meteor"))
        .astype(int)
        .filled(-1)
    )

    source_table["raic_satellite"] = (
        (
            (source_table["RAIC_Score"] > 0.5)
            * (source_table["RAIC_Label"] == "satellite")
        )
        .astype(int)
        .filled(-1)
    )
    source_table["raic_cont"] = (
        (
            (source_table["RAIC_Score"] > 0.3)
            * (
                (source_table["RAIC_Label"] == "calibissue")
                | (source_table["RAIC_Label"] == "baddither")
                | (source_table["RAIC_Label"] == "badfiber")  # Score > 0.3
                # | ( (source_table["RAIC_Label"] == "streak")
                # | (source_table["RAIC_Label"] == "lowcounts")
            )
        )
        .astype(int)
        .filled(-1)
    )
    source_table["flag_raic"] = np.invert(
        (source_table["raic_streak"] == 1)
        | (source_table["raic_cont"] == 1)
        | (source_table["raic_satellite"] == 1)
        | (source_table["raic_meteor"] == 1)
    ).astype(int)

    # suggested parameter cuts by Erin
    sel1 = (
        (source_table["chi2"] <= 1.5)
        * (source_table["linewidth"] <= 6)
        * (source_table["sn"] < 5.5)
    )
    sel2 = source_table["sn"] >= 5.5
    sel_2D = np.invert(
        (source_table["sn_moffat"] < 3.9)
        * ((source_table["sn"] < 8) | (source_table["wave"] < 3600))
    ).filled(1)

    flag_erin = (sel1 | sel2) * sel_2D

    source_table.add_column(Column(flag_erin, name="flag_erin_cuts", dtype=int))

    flag_best = (
        source_table["flag_badamp"]
        * source_table["flag_badfib"]
        * source_table["flag_badpix"]
        * source_table["flag_apcor"]
        * source_table["flag_wave"]
        * source_table["flag_3540"]
        * source_table['flag_continuum']
        * source_table["flag_baddet"]
        * source_table["flag_meteor"]
        * source_table["flag_largegal"]
        * source_table["flag_raic"]
        * source_table["flag_chi2fib"]
        * source_table["flag_pixmask"]
        * source_table["flag_lowlw"]
        * source_table['flag_satellite']
        * source_table['flag_cal']
        * np.invert(source_table["flag_dee"] == 0)
    )
    # update all AGN detections inspected by Chenxu to good
    flag_best[source_table["z_agn"] > 0] = 1

    flag_best  = flag_best.filled(1)
    
    source_table.add_column(Column(flag_best, name="flag_best", dtype=int), index=3)


    rres_vals = Table.read('/work/05350/ecooper/stampede2/hdr4/all_rres_4.0.1.txt', 
                       format='ascii', 
                       names=['detectid', 'rres_line', 'rres_cont'])

    source_table2 = join(source_table, rres_vals, join_type='left')

    print('Table sizes before and after rres join: {} {}'.format(len(source_table), len(source_table2)))
    source_table = source_table2
    
    sel_rres = source_table['rres_line']>=1.05
    source_table['sn_rres'] = source_table['sn']
    source_table['sn_rres'][sel_rres] = source_table['sn'][sel_rres]/source_table['rres_line'][sel_rres]

    for col in source_table.columns:

        if source_table[col].dtype == ">f8":
            if col in ["lum_lya", "lum_lya_err", "lum_oii", "lum_oii_err"]:
                pass
            else:
                source_table[col] = source_table[col].astype(np.float32)
            try:
                source_table[col] = source_table[col].filled(-999.0)
            except AttributeError:
                pass
        elif source_table[col].dtype == ">f4":
            try:
                source_table[col] = source_table[col].filled(-999.0)
            except AttributeError:
                pass
        elif (source_table[col].dtype == ">i4") | (source_table[col].dtype == ">i8"):
            try:
                source_table[col] = source_table[col].filled(-999)
            except AttributeError:
                pass
        elif source_table[col].dtype == "bool":
            continue
        elif col in ['rres_line', 'rres_cont']:
            source_table[col] = source_table.filled(1.0)
        else:
            try:
                source_table[col] = source_table[col].filled(" ")
            except:
                pass

    import joblib 

    clf = joblib.load("/work/05350/ecooper/stampede2/cosmos/calibration/rf_clf_2.0_20241106.joblib") 

    X =  np.array( [x[:] for x in source_table[clf.columns]] )

    p_real = clf.predict_proba(X)[:,1]

    source_table['p_real'] = p_real

    # remove nonsense metadata
    source_table.meta = {}

    # can include SA22 as of 2024-11-05
#    sel_sa22 = source_table["field"] == "ssa22"
#    sel_notsa22 = np.invert(sel_sa22)

    # temp solution for issue on 046, 036
    sel_ifuslot = (source_table['ifuslot'] == '046') | (source_table['ifuslot']=='036') 
    sel_date = (source_table['date'] >= 20240501)
    source_table['flag_best'][sel_ifuslot*sel_date] = -1
    
    source_table.write(
        "source_catalog_{}.z.fits".format(version), overwrite=True
    )
    source_table.write(
        "source_catalog_{}.z.tab".format(version), format="ascii", overwrite=True
    )
#    source_table[sel_sa22].write(
#        "source_catalog_{}_sa22.fits".format(version), overwrite=True
#    )

if __name__ == "__main__":
    main()
