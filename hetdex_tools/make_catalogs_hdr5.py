"""
Cloned from make_catalogs_hdr4.py on 2024-11-07. This replaces previous clone from 2024-10-03
as many updates were done

"""

import sys
import numpy as np
import numpy.ma as ma
import time
import argparse as ap
import os.path as op
import tables as tb
import gc

from astropy.table import Table, unique, vstack, join, Column, hstack
from astropy.coordinates import SkyCoord
import astropy.units as u

import matplotlib

matplotlib.use("agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from multiprocessing import Pool

from regions import LineSkyRegion, PixCoord, LinePixelRegion

from hetdex_api.config import HDRconfig
from hetdex_api.survey import Survey
from hetdex_api.detections import Detections
import hetdex_tools.fof_kdtree as fof

from elixer import catalogs

from hetdex_api.extinction import *
import extinction

catlib = catalogs.CatalogLibrary()
config = HDRconfig("hdr5")

wavelya = 1215.67
waveoii = 3727.8

deth5 = None
conth5 = None
add_agn = True


def make_friend_table_for_shot(shotid):
    global detect_table, dsky_2D
    sel_shot = detect_table["shotid"] == shotid
    kdtree, r = fof.mktree(
        detect_table["ra"][sel_shot],
        detect_table["dec"][sel_shot],
        np.zeros_like(detect_table["ra"][sel_shot]),
        dsky=dsky_2D,
    )
    friend_lst = fof.frinds_of_friends(kdtree, r, Nmin=1)

    friend_table = fof.process_group_list(
        friend_lst,
        detect_table["detectid"][sel_shot],
        detect_table["ra"][sel_shot],
        detect_table["dec"][sel_shot],
        0.0 * detect_table["wave"][sel_shot],
        detect_table["flux_g"][sel_shot],
    )

    memberlist = []
    friendlist = []
    for row in friend_table:
        friendid = row["id"]
        members = np.array(row["members"])
        friendlist.extend(friendid * np.ones_like(members))
        memberlist.extend(members)

    friend_table.remove_column("members")

    detfriend_tab = Table()
    detfriend_tab.add_column(Column(np.array(friendlist), name="id"))
    detfriend_tab.add_column(Column(memberlist, name="detectid"))

    detfriend_shot = join(detfriend_tab, friend_table, keys="id")

    return detfriend_shot


def make_wfriend_table_for_shot(shotid):
    global detect_table, dsky_3D, dwave

    sel_shot = (detect_table["shotid"] == shotid) & (detect_table["det_type"] == "line")
    kdtree, r = fof.mktree(
        detect_table["ra"][sel_shot],
        detect_table["dec"][sel_shot],
        detect_table["wave"][sel_shot],
        dsky=dsky_3D,
        dwave=dwave,
    )

    wfriend_lst = fof.frinds_of_friends(kdtree, r, Nmin=2)

    if len(wfriend_lst) > 0:
        wfriend_table = fof.process_group_list(
            wfriend_lst,
            detect_table["detectid"][sel_shot],
            detect_table["ra"][sel_shot],
            detect_table["dec"][sel_shot],
            detect_table["wave"][sel_shot],
            detect_table["flux"][sel_shot],
        )

        memberlist = []
        friendlist = []
        for row in wfriend_table:
            friendid = row["id"]
            members = np.array(row["members"])
            friendlist.extend(friendid * np.ones_like(members))
            memberlist.extend(members)

        wfriend_table.remove_column("members")

        wdetfriend_tab = Table()
        wdetfriend_tab.add_column(Column(np.array(friendlist), name="id"))
        wdetfriend_tab.add_column(Column(memberlist, name="detectid"))

        wdetfriend_shot = join(wdetfriend_tab, wfriend_table, keys="id")
        #        wdetfriend_shot.write('wfriend/detfriend_{}.tab'.format(shotid), format='ascii.no_header')
        return wdetfriend_shot
    else:
        return None


def add_elixer_cat_info(detect_table):
    global config

#    elixer_file1 = op.join(config.hdr_dir['hdr4'],'detect', "elixer_hdr3_hdr4_all_cat.h5")
#    elixer_file2 = op.join(config.hdr_dir['hdr5'],'detect', "elixer_hdr5_all_cat.h5")
    
#    elixer_cat1 = tb.open_file(elixer_file1, "r")
#    elixer_cat2 = tb.open_file(elixer_file2, "r")

#    elix_tab1 = Table(elixer_cat1.root.Detections.read())
#    elix_tab2 = Table(elixer_cat2.root.Detections.read())

    elixer_cat = tb.open_file( config.elixerh5, 'r')
    elix_tab = Table( elixer_cat.root.Detections.read())

    sel_det_col = [
        "detectid",
        "mag_g_wide",
        "plya_classification",
        "z_best_2",
        "z_best_pz_2",
        "combined_plae",
        "combined_plae_lo",
        "combined_plae_hi",
        "flags",
        "classification_labels",
    ]

    catalog = join(
        detect_table,
        unique(elix_tab[sel_det_col], keys="detectid"),
        keys="detectid",
        join_type="left",
    )

    catalog.rename_column("mag_g_wide", "gmag")
    catalog.rename_column("z_best_2", "best_z")
    catalog.rename_column("z_best_pz_2", "best_pz")
    
#    extracted_objects1 = Table(elixer_cat1.root.ExtractedObjects.read())
#    extracted_objects2 = Table(elixer_cat2.root.ExtractedObjects.read())

#    extracted_objects = vstack( [ extracted_objects1, extracted_objects2] )
    extracted_objects = Table( elixer_cat.root.ExtractedObjects.read())
                      
    selected1 = (extracted_objects["selected"] == True) & (
        extracted_objects["filter_name"] == b"r"
    )
    selected2 = (extracted_objects["catalog_name"] == b"HSC-DEX") | (
        extracted_objects["catalog_name"] == b"HSC-SSP"
    )
    selected = selected1 & selected2

    eo_tab = unique(extracted_objects[selected], keys="detectid")

    eo_tab.rename_column("flags", "flags_eo")
    eo_tab.rename_column("ra", "counterpart_ra")
    eo_tab.rename_column("dec", "counterpart_dec")
    eo_tab.rename_column("mag", "counterpart_mag")
    eo_tab.rename_column("mag_err", "counterpart_mag_err")
    eo_tab.rename_column("dist_baryctr", "counterpart_dist")
    eo_tab.rename_column("catalog_name", "counterpart_catalog_name")
    eo_tab.rename_column("filter_name", "counterpart_filter_name")
    eo_tab.rename_column("image_depth_mag", "counterpart_image_depth_mag")
    eo_tab["elixer_cat_src"] = "ExtractedObjects"

#    aper_tab1 = Table(elixer_cat1.root.Aperture.read())
#    aper_tab2 = Table(elixer_cat2.root.Aperture.read())
    aper_tab = Table( elixer_cat.root.Aperture.read())
                      
#   aper_tab = vstack( [ aper_tab1, aper_tab2] )
    # check if there is an aperture tab value
    sel_cat = (aper_tab["catalog_name"] == "HSC-SSP") | (
        aper_tab["catalog_name"] == "HSC-DEX"
    )
    sel_filter = aper_tab["filter_name"] == "r"

    a_tab = unique(aper_tab[sel_cat & sel_filter], keys="detectid")
    a_tab.rename_column("ra", "counterpart_ra")
    a_tab.rename_column("dec", "counterpart_dec")
    a_tab.rename_column("mag", "counterpart_mag")
    a_tab.rename_column("mag_err", "counterpart_mag_err")
    a_tab.rename_column("catalog_name", "counterpart_catalog_name")
    a_tab.rename_column("filter_name", "counterpart_filter_name")
    a_tab.rename_column("image_depth_mag", "counterpart_image_depth_mag")
    a_tab["counterpart_dist"] = 0.0
    a_tab["elixer_cat_src"] = "Aperture"

    cols_to_use = [
        "detectid",
        "counterpart_ra",
        "counterpart_dec",
        "counterpart_mag",
        "counterpart_mag_err",
        "counterpart_dist",
        "counterpart_catalog_name",
        "counterpart_filter_name",
        "counterpart_image_depth_mag",
        "elixer_cat_src",
    ]

    elix_stack_tab = unique(
        vstack([eo_tab[cols_to_use], a_tab[cols_to_use]]), keys="detectid"
    )

    catalog2 = join(catalog, elix_stack_tab, keys="detectid", join_type="left")

#    elixer_cat1.close()
#    elixer_cat2.close()
    elixer_cat.close()
                      
    catalog2.rename_column("best_z", "z_elixer")

    # fill mask values with nans
    for col in catalog2.columns:
        try:
            catalog2[col] = catalog2[col].filled(np.nan)
        except Exception:
            pass

    return catalog2


def merge_wave_groups(wid):
    global expand_table

    try:
        sel_wid = expand_table["wave_group_id"] == wid
        grp = expand_table[sel_wid]
        sid, ns = np.unique(grp["source_id"], return_counts=True)

        sid_main = np.min(sid)

        sid_ind = []
        for sid_i in sid:
            for ind in np.where(expand_table["source_id"] == sid_i)[0]:
                sid_ind.append(ind)

        # now find any other wave groups and their associated source_id info to merge
        other_wids = np.unique(expand_table["wave_group_id"][sid_ind])

        for wid_i in other_wids:
            if wid_i == 0:
                continue
            elif wid_i == wid:
                continue

            sel_wid_i = expand_table["wave_group_id"] == wid_i
            grp = expand_table[sel_wid_i]
            sid, ns = np.unique(grp["source_id"], return_counts=True)

            for sid_i in sid:
                if sid_i == sid_main:
                    continue
                for ind in np.where(expand_table["source_id"] == sid_i)[0]:
                    sid_ind.append(ind)

        if np.size(np.unique(expand_table["source_id"][sid_ind])) > 1:
            return sid_ind
        else:
            return None
    except Exception:
        print("Merge wave group failed for {}".format(wid))
        return None


def create_source_catalog(version="5.0.0", update=False):
    global config

    if update:
        print("Creating curated detection catalog version={}".format(version))

        print("Opening HDR3 Lines Database")
        D_hdr3 = Detections(survey="hdr3", loadtable=True)

        # get masking info for each detection
        mask_badamp = D_hdr3.remove_bad_amps()
        mask_badshots = D_hdr3.remove_shots()
        mask_tp = D_hdr3.throughput >= 0.08
        # downselect badamps and badshots

        sel_date = D_hdr3.date <= 20210901

        D_hdr3 = D_hdr3[mask_badamp & mask_badshots & mask_tp & sel_date]

        mask_baddet_hdr3 = D_hdr3.remove_bad_detects()

        sel_cut1 = (D_hdr3.sn >= 4.8) & (D_hdr3.chi2 <= 2.5)
        sel_cont = D_hdr3.continuum > -3
        sel_chi2fib = D_hdr3.chi2fib < 4.5
        sel_lw = (D_hdr3.linewidth <= 14) & (D_hdr3.linewidth >= 1.6)

        sel_cat = sel_cut1 & sel_cont & sel_chi2fib * sel_lw

        detects_line_table_hdr3 = D_hdr3[sel_cat].return_astropy_table()
        detects_line_table_hdr3.add_column(
            Column(str("line"), name="det_type", dtype=str)
        )

        detects_line_table_hdr3.add_column(
            Column(str("hdr3"), name="survey", dtype=str)
        )

        detects_line_table_hdr3.add_column(
            Column(mask_baddet_hdr3[sel_cat].astype(int), name="flag_baddet", dtype=int)
        )

        print("Opening HDR4 Lines Database")

        D_hdr4 = Detections(survey="hdr4", loadtable=True)

        # get masking info for each detection
        mask_badamp = D_hdr4.remove_bad_amps()
        mask_badshots = D_hdr4.remove_shots()
        mask_tp = D_hdr4.throughput >= 0.08
        # downselect badamps and badshots
        D_hdr4 = D_hdr4[mask_badamp & mask_badshots & mask_tp]

        mask_baddet_hdr4 = D_hdr4.remove_bad_detects()

        sel_cut1 = (D_hdr4.sn >= 4.8) & (D_hdr4.chi2 <= 2.5)
        sel_cont = D_hdr4.continuum > -3
        sel_chi2fib = D_hdr4.chi2fib < 4.5
        sel_lw = (D_hdr4.linewidth <= 14) & (D_hdr4.linewidth >= 1.6)

        sel_cat = sel_cut1 & sel_cont & sel_chi2fib * sel_lw

        detects_line_table_hdr4 = D_hdr4[sel_cat].return_astropy_table()

        detects_line_table_hdr4.add_column(
            Column(str("line"), name="det_type", dtype=str)
        )
        detects_line_table_hdr4.add_column(
            Column(str("hdr4"), name="survey", dtype=str)
        )
        detects_line_table_hdr4.add_column(
            Column(mask_baddet_hdr4[sel_cat].astype(int), name="flag_baddet", dtype=int)
        )

        print("Opening HDR5 Lines Database")

        D_hdr5 = Detections(survey="hdr5", loadtable=True)

        # get masking info for each detection                                                                                 
        mask_badamp = D_hdr5.remove_bad_amps()
        mask_badshots = D_hdr5.remove_shots()
        mask_tp = D_hdr5.throughput >= 0.08
        # downselect badamps and badshots                                                                                     
        D_hdr5 = D_hdr5[mask_badamp & mask_badshots & mask_tp]

        mask_baddet_hdr5 = D_hdr5.remove_bad_detects()

        sel_cut1 = (D_hdr5.sn >= 4.8) & (D_hdr5.chi2 <= 2.5)
        sel_cont = D_hdr5.continuum > -3
        sel_chi2fib = D_hdr5.chi2fib < 4.5
        sel_lw = (D_hdr5.linewidth <= 14) & (D_hdr5.linewidth >= 1.6)

        sel_cat = sel_cut1 & sel_cont & sel_chi2fib * sel_lw

        detects_line_table_hdr5 = D_hdr5[sel_cat].return_astropy_table()

        detects_line_table_hdr5.add_column(
            Column(str("line"), name="det_type", dtype=str)
        )
        detects_line_table_hdr5.add_column(
            Column(str("hdr5"), name="survey", dtype=str)
        )
        detects_line_table_hdr5.add_column(
            Column(mask_baddet_hdr5[sel_cat].astype(int), name="flag_baddet", dtype=int)
        )

        detects_line_table = vstack([detects_line_table_hdr3, detects_line_table_hdr4, detects_line_table_hdr5])

        print("Adding Elixer Info")
        print(len(detects_line_table))

        detects_line_table = add_elixer_cat_info(detects_line_table)
        print(len(detects_line_table))

        detects_line_table.write("detect_hdr{}.fits".format(version), overwrite=True)
        detects_line_table.write(
            "detect_hdr{}.tab".format(version), format="ascii", overwrite=True
        )
    else:
        detects_line_table = Table.read("detect_hdr{}.fits".format(version))


    if update:
        print("Opening HDR3 Continuum Database")
        detects_cont_hdr3 = Detections(
            catalog_type="continuum", survey="hdr3", loadtable=True
        )

        mask_badamp_cont = detects_cont_hdr3.remove_bad_amps()
        mask_badshot_cont = detects_cont_hdr3.remove_shots()
        mask_tp_cont = detects_cont_hdr3.throughput > 0.08
        sel_date = detects_cont_hdr3.date <= 20210901

        detects_cont_hdr3 = detects_cont_hdr3[
            mask_badamp_cont & mask_badshot_cont & mask_tp_cont & sel_date
        ]

        mask_baddet_cont = detects_cont_hdr3.remove_bad_detects()

        detects_cont_table_hdr3 = detects_cont_hdr3.return_astropy_table()

        detects_cont_table_hdr3.add_column(
            Column(str("cont"), name="det_type", dtype=str)
        )
        detects_cont_table_hdr3.add_column(
            Column(str("hdr3"), name="survey", dtype=str)
        )
        detects_cont_table_hdr3.add_column(
            Column(mask_baddet_cont.astype(int), name="flag_baddet", dtype=int)
        )

        print("Opening HDR4 Continuum Database")
        detects_cont_hdr4 = Detections(
            catalog_type="continuum", survey="hdr4", loadtable=True
        )

        mask_badamp_cont = detects_cont_hdr4.remove_bad_amps()
        mask_badshot_cont = detects_cont_hdr4.remove_shots()
        mask_tp_cont = detects_cont_hdr4.throughput > 0.08

        detects_cont_hdr4 = detects_cont_hdr4[
            mask_badamp_cont & mask_badshot_cont & mask_tp_cont
        ]

        mask_baddet_cont = detects_cont_hdr4.remove_bad_detects()

        detects_cont_table_hdr4 = detects_cont_hdr4.return_astropy_table()

        detects_cont_table_hdr4.add_column(
            Column(str("cont"), name="det_type", dtype=str)
        )
        detects_cont_table_hdr4.add_column(
            Column(str("hdr4"), name="survey", dtype=str)
        )
        detects_cont_table_hdr4.add_column(
            Column(mask_baddet_cont.astype(int), name="flag_baddet", dtype=int)
        )

        print("Opening HDR5 Continuum Database")

        detects_cont_hdr5 = Detections(
            catalog_type="continuum", survey="hdr5", loadtable=True
        )

        mask_badamp_cont = detects_cont_hdr5.remove_bad_amps()
        mask_badshot_cont = detects_cont_hdr5.remove_shots()
        mask_tp_cont = detects_cont_hdr5.throughput > 0.08

        detects_cont_hdr5 = detects_cont_hdr5[
            mask_badamp_cont & mask_badshot_cont & mask_tp_cont
        ]

        mask_baddet_cont = detects_cont_hdr5.remove_bad_detects()

        detects_cont_table_hdr5 = detects_cont_hdr5.return_astropy_table()

        detects_cont_table_hdr5.add_column(
            Column(str("cont"), name="det_type", dtype=str)
        )
        detects_cont_table_hdr5.add_column(
            Column(str("hdr5"), name="survey", dtype=str)
        )
        detects_cont_table_hdr5.add_column(
            Column(mask_baddet_cont.astype(int), name="flag_baddet", dtype=int)
        )

        print("Combining Continuum Databases")

        detects_cont_table = vstack([detects_cont_table_hdr3, detects_cont_table_hdr4, detects_cont_table_hdr5])

        # set columns to 0 that are not relevent to continuum catalog
        for col in ["apcor", "flux_noise_1sigma", "sn_3fib", "sn_3fib_cen", "sn_cen"]:
            detects_cont_table[col] = 0.0

        full_cont_table = detects_cont_table.copy()

        print("Adding Elixer info to combined continuum catalog")

        detects_cont_table = add_elixer_cat_info(detects_cont_table)
        print(len(detects_cont_table))
        detects_cont_table.write("continuum_" + version + ".fits", overwrite=True)
        detects_cont_table.write(
            "continuum_" + version + ".tab", overwrite=True, format="ascii"
        )
        detects_cont_hdr3.close()
        detects_cont_hdr4.close()
        detects_cont_hdr5.close()

    else:
        detects_cont_table = Table.read("continuum_" + version + ".fits")

    if update:
        # create an agn table with detection info added
        if add_agn:
            print('Combining AGN catalog with full database info.')
            full_line_table_hdr3 = D_hdr3.return_astropy_table()
            full_line_table_hdr3.add_column(
                Column(str("hdr3"), name="survey", dtype=str)
            )
            full_line_table_hdr3.add_column(
                Column(mask_baddet_hdr3.astype(int), name="flag_baddet", dtype=int)
            )

            full_line_table_hdr4 = D_hdr4.return_astropy_table()
            full_line_table_hdr4.add_column(
                Column(str("hdr4"), name="survey", dtype=str)
            )
            full_line_table_hdr4.add_column(
                Column(mask_baddet_hdr4.astype(int), name="flag_baddet", dtype=int)
            )

            full_line_table_hdr5 = D_hdr5.return_astropy_table()
            full_line_table_hdr5.add_column(
                Column(str("hdr5"), name="survey", dtype=str)
            )
            
            full_line_table_hdr5.add_column(
                Column(mask_baddet_hdr5.astype(int), name="flag_baddet", dtype=int)
            )
            full_line_table = vstack(
                [full_line_table_hdr3,
                 full_line_table_hdr4,
                 full_line_table_hdr5]
            )

            full_line_table.add_column(Column(str("line"), name="det_type", dtype=str))

            agn_tab = Table.read(
                config.agncat,
                format="ascii",
                include_names=["detectid", "vis_class", "z"],
            )
            agn_tab.rename_column("vis_class", "agn_vis_class")
            agn_tab.rename_column("z", "z_agn")

            detects_agn = join(
                agn_tab['detectid', 'agn_vis_class','z_agn'],
                vstack([full_cont_table, full_line_table]),
                join_type="inner",
            )
            detects_agn = add_elixer_cat_info(detects_agn)

            detects_agn.sort("gmag")
            detects_agn = unique(detects_agn, keys="detectid")
            print("Final AGN length: {}".format(len(detects_agn)))
            detects_agn.write("agn_" + version + ".fits", overwrite=True)
        else:
            print("No updated AGN catalog created")
        D_hdr3.close()
        D_hdr4.close()
        D_hdr5.close()

    else:
        if add_agn:
            detects_agn = Table.read("agn_" + version + ".fits")
        else:
            detects_agn = None
    if update:
        return

    global detect_table

    if add_agn:
        print('Adding AGN table')
        
        detect_table = unique(
            vstack([detects_agn, detects_cont_table, detects_line_table]),
            keys="detectid",
        )
        detect_table["z_agn"] = detect_table["z_agn"].filled(-1)
        detect_table["agn_vis_class"] = detect_table["agn_vis_class"].filled(-1)

        del detects_cont_table, detects_line_table, detects_agn
    else:
        detect_table = unique(
            vstack([detects_cont_table, detects_line_table]), keys="detectid"
        )

        del detects_cont_table, detects_line_table

    # save the continuum and line emission dets in separate files for det flag calculations

    for survey in ["hdr3", "hdr4", "hdr5"]:
        np.savetxt(
            "line_{}_{}.dets".format(survey, version),
            detect_table["detectid"][
                (detect_table["det_type"] == "line")
                * (detect_table["survey"] == survey)
            ],
            fmt="%i",
        )
        np.savetxt(
            "cont_{}_{}.dets".format(survey, version),
            detect_table["detectid"][
                (detect_table["det_type"] == "cont")
                * (detect_table["survey"] == survey)
            ],
            fmt="%i",
        )
        np.savetxt(
            "shots_{}_{}.txt".format(survey, version),
            np.unique(detect_table["shotid"][detect_table["survey"] == survey]),
            fmt="%i",
        )

    gc.collect()

    # add det flags table

    detflags_tab = Table.read('/scratch/projects/hetdex/hdr5/catalogs/det_flags_{}.fits'.format(version) )

    print("Adding det flags to full catalog. Size is {}".format(len(detect_table)))
    detect_table2 = join(detect_table, detflags_tab, keys="detectid", join_type="left")
    detect_table = unique(detect_table2, keys="detectid")

    print("Size after combining with detflags_tab: {}".format(len(detect_table)))

    # fill mask values with nans
    for col in detect_table.columns:
        try:
            if col in [
                "flag_pixmask",
                "flag_badamp",
                "flag_badpix",
                "flag_badfib",
                "flag_meteor",
                "flag_largegal",
                "flag_chi2fib",
                "flag_baddetect",
            ]:
                detect_table[col] = detect_table[col].filled(1)
                print("Replacing masked values with 1 for {}".format(col))
            else:
                detect_table[col] = detect_table[col].filled(np.nan)
                print("Replacing masked values with np.nan for {}".format(col))
        except Exception:
            pass
            print("no", col)

    # calculate ebv and av for every detections
    all_coords = SkyCoord(ra=detect_table["ra"], dec=detect_table["dec"], unit="deg")
    sfd = SFDQuery()
    ebv = sfd(all_coords)
    Rv = 3.1
    corr_SF2011 = 2.742  # Landolt V
    ext = []

    Av = corr_SF2011 * ebv

    detect_table["Av"] = Av
    detect_table["ebv"] = ebv

    # get fluxes to derive flux-weighted distribution of group

    detect_table["gmag"][np.isnan(detect_table["gmag"])] = 27
    gmag = detect_table["gmag"] * u.AB
    detect_table["flux_g"] = gmag.to(u.Jy).value

    # first cluster in positional/wavelegnth space in a larger
    # linking length

    shotlist = np.unique(detect_table["shotid"])

    print(
        "Performing FOF in 3D space wiht dlink={} and dwave={}".format(dsky_3D, dwave)
    )
    t0 = time.time()
    P = Pool(20)
    res = P.map(make_wfriend_table_for_shot, shotlist)
    P.close()

    wstart = int(version.replace(".", "", 2)) * 10**6
    wdetfriend_all = Table()

    firstid = True
    for r in res:
        if r is None:
            continue
        else:
            if firstid:
                r["id"] = r["id"] + wstart
                firstid = False
            else:
                r["id"] = r["id"] + np.max(wdetfriend_all["id"]) + 1
        wdetfriend_all = vstack([wdetfriend_all, r])

    wdetfriend_all.rename_column("id", "wave_group_id")
    wdetfriend_all.rename_column("size", "wave_group_size")
    wdetfriend_all.rename_column("a", "wave_group_a")
    wdetfriend_all.rename_column("b", "wave_group_b")
    wdetfriend_all.rename_column("pa", "wave_group_pa")
    wdetfriend_all.rename_column("icx", "wave_group_ra")
    wdetfriend_all.rename_column("icy", "wave_group_dec")
    wdetfriend_all.rename_column("icz", "wave_group_wave")

    w_keep = wdetfriend_all[
        "detectid",
        "wave_group_id",
        "wave_group_a",
        "wave_group_b",
        "wave_group_pa",
        "wave_group_ra",
        "wave_group_dec",
        "wave_group_wave",
    ]

    t1 = time.time()

    print("3D FOF analysis complete in {:3.2f} minutes \n".format((t1 - t0) / 60))

    print("Performing FOF in 2D space with dlink={}".format(dsky_2D))

    t0 = time.time()
    P = Pool(20)
    res = P.map(make_friend_table_for_shot, shotlist)
    P.close()

    sid_start = int(version.replace(".", "", 2)) * 10**10

    firstid = True
    detfriend_all = Table()
    for r in res:
        if firstid:
            r["id"] = r["id"] + sid_start
            firstid = False
        else:
            r["id"] = r["id"] + np.max(detfriend_all["id"]) + 1
        detfriend_all = vstack([detfriend_all, r])

    t1 = time.time()
    print("2D FOF analysis complete in {:3.2f} minutes \n".format((t1 - t0) / 60))

    # match detectids at large linking length if a wave group exists
    joinfriend = join(detfriend_all, w_keep, keys="detectid", join_type="left")

    global expand_table

    expand_table = join(joinfriend, detect_table, keys="detectid")
    expand_table["wave_group_id"] = expand_table["wave_group_id"].filled(0)
    expand_table.rename_column("id", "source_id")

    try:
        expand_table.remove_column("detectname")
    except:
        pass

    # combine common wavegroups to the same source_id
    # update source properties

    del joinfriend

    print("Combining nearby wavegroups and detections")

    t0 = time.time()
    sel = expand_table["wave_group_id"] > 0
    wid_list = np.unique(expand_table["wave_group_id"][sel])
    p = Pool(16)#reduced from 20 20241109
    res = p.map(merge_wave_groups, wid_list)
    p.close()

    ind_list = []
    for r in res:
        if r is None:
            continue
        else:
            r.sort()
            if r in ind_list:  # skip if indices are a duplicate
                continue
            else:
                ind_list.append(r)

    res_tab = fof.process_group_list(
        ind_list,
        expand_table["detectid"],
        expand_table["ra"],
        expand_table["dec"],
        0.0 * expand_table["wave"],
        expand_table["flux_g"],
    )

    for i, ind in enumerate(ind_list):
        sid_main = np.min(expand_table["source_id"][ind])
        expand_table["source_id"][ind] = sid_main

        for col in res_tab.colnames:
            if col in ["id", "members"]:
                continue
            expand_table[col][ind] = res_tab[col][i]

    t1 = time.time()
    print("Done combining wavegroups in {:4.2f} min".format((t1 - t0) / 60))

    del detfriend_all, detect_table, detect_table2

    expand_table.rename_column("size", "n_members")
    expand_table.rename_column("icx", "ra_mean")
    expand_table.rename_column("icy", "dec_mean")

    source_name = []

    for row in expand_table:
        source_name.append(get_source_name(row["ra_mean"], row["dec_mean"]))
    try:
        expand_table.add_column(source_name, name="source_name", index=1)
    except Exception:
        print("messed up source name again")

    expand_table.sort("source_id")

    # fill mask values with nans
    for col in expand_table.columns:
        try:
            expand_table[col] = expand_table[col].filled(np.nan)
            print("yes", col)
        except:
            pass
            print("no", col)
    return expand_table


def plot_source_group(
    source_id=None, source_table=None, k=3.5, vmin=3, vmax=99, label=True, save=False
):
    """
    Plot a unique source group from the HETDEX
    unique source catalog

    Parameters
    ----------
    source_id: int

    """

    if source_table is None:
        print("Please provide source catalog (an astropy table)")
    else:
        sel = source_table["source_id"] == source_id
        group = source_table[sel]

    if source_id is None:
        print("Please provide source_id (an integer)")

    ellipse = False

    # get ellipse parameters if more than 1 source
    if np.size(group) > 1:
        ellipse = True

    if ellipse:
        a, b, pa, a2, b2, pa2 = (
            group["a"][0],
            group["b"][0],
            group["pa"][0],
            group["a2"][0],
            group["b2"][0],
            group["pa2"][0],
        )

        cosd = np.cos(np.deg2rad(group["dec_mean"][0]))
        imsize = np.max(
            [
                np.max(
                    [
                        (np.max(group["ra"]) - np.min(group["ra"])) * cosd,
                        (np.max(group["dec"]) - np.min(group["dec"])),
                    ]
                )
                * 1.7,
                10.0 / 3600.0,
            ]
        )

    else:
        imsize = 10.0 / 3600.0

    coords = SkyCoord(ra=group["ra_mean"][0] * u.deg, dec=group["dec_mean"][0] * u.deg)

    cutout = catlib.get_cutouts(
        position=coords,
        side=imsize,
        filter=["f606W", "r", "g"],
        first=True,
        allow_bad_image=False,
        allow_web=True,
    )[0]

    im = cutout["cutout"].data
    wcs = cutout["cutout"].wcs
    plt.figure(figsize=(8, 8))

    plt.subplot(projection=wcs)
    ax = plt.gca()
    ax.coords[0].set_major_formatter("d.dddd")
    # ax.coords[0].set_ticks(spacing=1. * u.arcsec)
    ax.coords[1].set_major_formatter("d.dddd")
    # ax.coords[0].set_ticks(spacing=1. * u.arcsec)

    impix = im.shape[0]
    pixscale = imsize / impix  # in degrees/pix
    m = np.percentile(im, (vmin, vmax))
    plt.imshow(im, vmin=m[0], vmax=m[1], origin="lower", cmap="gray_r")
    plt.text(
        0.95,
        0.05,
        cutout["instrument"] + cutout["filter"],
        transform=ax.transAxes,
        fontsize=20,
        color="red",
        horizontalalignment="right",
    )

    # plot the group members
    sel_line = group["det_type"] == "line"
    if np.sum(sel_line) >= 1:
        plt.scatter(
            group["ra"][sel_line],
            group["dec"][sel_line],
            transform=ax.get_transform("world"),
            marker="o",
            color="orange",
            linewidth=4,
            s=group["sn"][sel_line],
            zorder=100,
            label="line emission",
        )

    sel_cont = group["det_type"] == "cont"
    if np.sum(sel_cont) >= 1:
        plt.scatter(
            group["ra"][sel_cont],
            group["dec"][sel_cont],
            transform=ax.get_transform("world"),
            marker="o",
            color="green",
            linewidth=4,
            s=10,
            label="continuum",
        )

    sel_agn = group["det_type"] == "agn"
    if np.sum(sel_agn) >= 1:
        plt.scatter(
            group["ra"][sel_agn],
            group["dec"][sel_agn],
            transform=ax.get_transform("world"),
            marker="o",
            color="red",
            linewidth=4,
            s=10,
            label="agn",
        )

    # plot and elliptical kron-like aperture representing the group. Ellipse
    # in world coords does not work, so plot in pixelcoordinates...
    # East is to the right in these plots, so pa needs transform

    if ellipse:
        ellipse = Ellipse(
            xy=(impix // 2, impix // 2),
            width=k * a2 / pixscale,
            height=k * b2 / pixscale,
            angle=180 - pa2,
            edgecolor="r",
            fc="None",
            lw=1,
        )
        ax.add_patch(ellipse)

        ellipse = Ellipse(
            xy=(impix // 2, impix // 2),
            width=a / pixscale,
            height=b / pixscale,
            angle=180 - pa,
            edgecolor="b",
            fc="None",
            lw=1,
        )
        ax.add_patch(ellipse)

    # add 5 arcsec scale bar
    x1 = 0.05 * impix
    y1 = 0.05 * impix
    y2 = y1 + (5.0 / 3600) / pixscale
    start = PixCoord(x=x1, y=y1)
    end = PixCoord(x=x1, y=y2)
    reg = LinePixelRegion(start=start, end=end)
    plt.text(x1, 0.025 * impix, "5 arcsec", color="blue")
    patch = reg.as_artist(facecolor="none", edgecolor="blue", lw=4)
    ax.add_patch(patch)

    if label:
        # plot detecid labels
        for row in group:
            plt.text(
                row["ra"],
                row["dec"],
                str(row["detectid"]),
                transform=ax.get_transform("world"),
                fontsize=9,
                color="red",
            )
    try:
        z_hetdex = group["z_hetdex"][0]
    except:
        z_hetdex = group["z_guess"][0]

    plt.title(
        "source_id:%d n:%d ra:%6.3f dec:%6.3f z:%5.3f"
        % (
            source_id,
            group["n_members"][0],
            group["ra_mean"][0],
            group["dec_mean"][0],
            z_hetdex,
        )
    )

    plt.xlabel("RA")
    plt.ylabel("DEC")
    plt.legend()

    if save:
        plt.savefig("figures/source-%03d.png" % source_id, format="png")
        plt.close()
    else:
        plt.show()


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
        "-dsky_3D",
        "--dsky_3D",
        type=float,
        help="""Spatial linking length in arcsec""",
        default=6.0,
    )

    parser.add_argument(
        "-dsky_2D",
        "--dsky_2D",
        type=float,
        help="""Spatial linking length in arcsec""",
        default=3.0,
    )

    parser.add_argument(
        "-dwave",
        "--dwave",
        type=float,
        help="""Wavelength linking length in A""",
        default=8.0,
    )

    parser.add_argument(
        "--update",
        "-u",
        help="""Trigger to update detect and continuum catalogs""",
        default=False,
        required=False,
        action="store_true",
    )
    return parser


def get_source_name(ra, dec):
    """
    convert ra,dec coordinates to a IAU-style object name.
    """
    coord = SkyCoord(ra * u.deg, dec * u.deg)
    return "HETDEX_J{0}{1}".format(
        coord.ra.to_string(unit=u.hourangle, sep="", precision=2, pad=True),
        coord.dec.to_string(sep="", precision=1, alwayssign=True, pad=True),
    )


def main(argv=None):
    """Main Function"""

    parser = get_parser()
    args = parser.parse_args(argv)

    global source_table, dsky_2D, dsky_3D, dwave
    dsky_2D = args.dsky_2D
    dsky_3D = args.dsky_3D
    dwave = args.dwave

    if args.update:
        source_table = create_source_catalog(version=args.version, update=args.update)
        sys.exit("Done updating det and cont catalogs. Re-run without --update")

    source_table = create_source_catalog(version=args.version)

    # add z_spec from other catlogs if it exists
    # from Steve F.
    goods_z = Table.read(
        op.join(config.imaging_dir, "catalogs", "goods_n_specz_1018_no_header.txt"),
        names=[
            "ra_zspec",
            "dec_zspec",
            "zspec",
            "z_quality",
            "zspec_catalog",
            "Symbol",
        ],
        format="ascii.no_header",
    )
    # remove z_quality = 1,2, low quality zspec
    sel_bad1 = goods_z["z_quality"] == "1"
    sel_bad2 = goods_z["z_quality"] == "2"
    sel_good_goods = np.invert(sel_bad1 | sel_bad2)
    goods_z = goods_z[sel_good_goods]

    # DEIMOS 10k (Hasinger et al. 2018) z_spec up to ~6
    # https://cosmos.astro.caltech.edu/news/65
    # Secure if comprehensive quality flag (Q) = 2
    deimos = Table.read(
        op.join(config.imaging_dir, "catalogs", "deimos_redshifts.tbl"), format="ascii"
    )
    deimos.rename_column("Ra", "ra_zspec")
    deimos.rename_column("Dec", "dec_zspec")
    deimos["zspec_catalog"] = "CosmosDeimos"
    sel_good_deimos = deimos["Q"] == 2
    deimos = deimos[sel_good_deimos]

    # Kriek et al. (2015)
    # http://mosdef.astro.berkeley.edu/for-scientists/data-releases/
    # 5: Redshift is based on single emission feature detected with 2<=S/N<3,
    # and is within 95% confidence interval of photo-z or within delta(z)=0.05
    # of pre-MOSFIRE spec-z (if it exists)
    mosdef = Table.read(
        op.join(config.imaging_dir, "catalogs", "mosdef_zcat.final_slitap.fits")
    )
    mosdef.rename_column("RA", "ra_zspec")
    mosdef.rename_column("DEC", "dec_zspec")
    mosdef.rename_column("Z_MOSFIRE", "zspec")
    mosdef["zspec_catalog"] = "MOSFIRE"
    sel_good_mosdef = mosdef["Z_MOSFIRE_ZQUAL"] >= 5
    mosdef = mosdef[sel_good_mosdef]

    # VUDS (Tasca et al. 2017), z_spec up to ~6
    # http://cesam.lam.fr/vuds/DR1/
    # zCOSMOS Spectroscopic Redshift Survey 2009ApJS..184..218L
    # https://www.eso.org/qi/catalog/show/65

    # Flag 4’s, 3’s, 2.5, 2.4, 1.5, 9.5, 9.4, 9.3 are considered secure.
    zcosbright = Table.read(
        op.join(
            config.imaging_dir,
            "catalogs",
            "cesam_zcosbrightspec20k_dr3_catalog_1616073679.txt",
        ),
        format="ascii",
    )
    zcosbright.rename_column("zpec", "zspec")
    zcosbright.rename_column("ra", "ra_zspec")
    zcosbright.rename_column("dec", "dec_zspec")
    zcosbright["zspec_catalog"] = "zCosmosBright"
    sel_zcos_good1 = zcosbright["cc"] >= 9.3
    sel_zcos_good2 = (zcosbright["cc"] >= 3.0) & (zcosbright["cc"] < 5)
    sel_zcos_good3 = zcosbright["cc"] == 1.5
    sel_zcos_good4 = (zcosbright["cc"] >= 2.4) & (zcosbright["cc"] <= 2.5)
    sel_zcos_good = sel_zcos_good1 | sel_zcos_good2 | sel_zcos_good3 | sel_zcos_good4
    zcosbright = zcosbright[sel_zcos_good]

    deep_specz = Table.read(
        op.join(
            config.imaging_dir,
            "catalogs",
            "DEEP_zcosmos_spectroscopy_one_v2.6_data+header.cat",
        ),
        format="ascii",
        data_start=100,
    )
    deep_specz.rename_column("col1", "zCOSMOS-deepID")
    deep_specz.rename_column("col2", "zspec")
    deep_specz.rename_column("col3", "flag")
    deep_specz.rename_column("col4", "zphot")
    deep_specz.rename_column("col5", "ra_zspec")
    deep_specz.rename_column("col6", "dec_zspec")
    deep_specz["zspec_catalog"] = "DEEP_zcosmos"
    sel_deepcos_good1 = deep_specz["flag"] >= 9.3
    sel_deepcos_good2 = (deep_specz["flag"] >= 3.0) & (deep_specz["flag"] < 5)
    sel_deepcos_good3 = deep_specz["flag"] == 1.5
    sel_deepcos_good4 = (deep_specz["flag"] >= 2.4) & (deep_specz["flag"] <= 2.5)
    sel_deepcos_good = (
        sel_deepcos_good1 | sel_deepcos_good2 | sel_deepcos_good3 | sel_deepcos_good4
    )
    deep_specz = deep_specz[sel_deepcos_good]

    sdssfile = op.join(config.imaging_dir, "sdss", "specObj-dr16.fits")

    sdss_specz = Table.read(sdssfile)
    sel_good_sdss = sdss_specz["ZWARNING"] == 0
    sdss_specz = sdss_specz[sel_good_sdss]

    sdss_specz.rename_column("PLUG_RA", "ra_zspec")
    sdss_specz.rename_column("PLUG_DEC", "dec_zspec")
    sdss_specz.rename_column("CLASS", "sdss_class")
    sdss_specz.rename_column("Z", "zspec")
    sdss_specz.rename_column("Z_ERR", "z_sdss_err")
    sdss_specz["zspec_catalog"] = "sdss-dr16"

    dtab = Table.read(op.join(config.imaging_dir, "catalogs", "desi-hetdex-v2.0.fits"))

    sel_good_hetdex = dtab["wave"] >= 3640
    sel_good_desi = dtab["COADD_FIBERSTATUS"] == 0
    sel_sample = sel_good_desi * sel_good_hetdex
    sel_conf = dtab["VI_quality"] >= 3
    dtab = dtab[sel_conf & sel_sample]
    dtab.rename_column("TARGET_RA", "ra_zspec")
    dtab.rename_column("TARGET_DEC", "dec_zspec")
    dtab["zspec_catalog"] = "DESI"
    dtab.rename_column("VI_z", "zspec")

    specz_catalogs = vstack(
        [
            goods_z["zspec", "ra_zspec", "dec_zspec", "zspec_catalog"],
            deimos["zspec", "ra_zspec", "dec_zspec", "zspec_catalog"],
            mosdef["zspec", "ra_zspec", "dec_zspec", "zspec_catalog"],
            zcosbright["zspec", "ra_zspec", "dec_zspec", "zspec_catalog"],
            deep_specz["zspec", "ra_zspec", "dec_zspec", "zspec_catalog"],
            sdss_specz["zspec", "ra_zspec", "dec_zspec", "zspec_catalog"],
            dtab["zspec", "ra_zspec", "dec_zspec", "zspec_catalog"],
        ]
    )
    sel = specz_catalogs["zspec"] >= 0
    specz_catalogs = specz_catalogs[sel]

    specz_coords = SkyCoord(
        ra=specz_catalogs["ra_zspec"], dec=specz_catalogs["dec_zspec"], unit="deg"
    )
    source_coords = SkyCoord(ra=source_table["ra"], dec=source_table["dec"], unit="deg")

    idx, d2d, d3d = source_coords.match_to_catalog_sky(specz_coords)

    catalog_matches = specz_catalogs[idx]
    catalog_matches["zspec_dist"] = d2d.to_value(u.arcsec)

    matched_catalog = hstack([source_table, catalog_matches])

    sel_close = matched_catalog["zspec_dist"] < 1.5  # u.arcsec

    print(
        "There are {} zspec matches within 5 arcsec".format(
            np.size(np.unique(matched_catalog["source_id"][sel_close]))
        )
    )

    sel_remove = matched_catalog["zspec_dist"] > 1.5  # u.arcsec

    matched_catalog["zspec"][sel_remove] = np.nan
    matched_catalog["ra_zspec"][sel_remove] = np.nan
    matched_catalog["dec_zspec"][sel_remove] = np.nan
    matched_catalog["zspec_dist"][sel_remove] = np.nan
    matched_catalog["zspec_catalog"][sel_remove] = ""

    source_table = matched_catalog

    # Clear up memory

    for name in dir():
        if source_table:
            continue
        elif not name.startswith("_"):
            del name

    import gc

    gc.collect()

    # sort table closest to group mean position
    # so unique will produce the closest match

    src_coord = SkyCoord(
        ra=source_table["ra_mean"], dec=source_table["dec_mean"], unit="deg"
    )
    det_coord = SkyCoord(ra=source_table["ra"], dec=source_table["dec"], unit="deg")

    source_table["src_separation"] = det_coord.separation(src_coord)
    source_table.sort(["src_separation"])

    source_table.meta = {}

    # Apply dust corretion to flux, continuum, flux_1sigma

    source_table["flux_obs"] = source_table["flux"].copy()
    source_table["flux_err_obs"] = source_table["flux_err"].copy()
    source_table["continuum_obs"] = source_table["continuum"].copy()
    source_table["continuum_err_obs"] = source_table["continuum_err"].copy()
    source_table["flux_noise_1sigma_obs"] = source_table["flux_noise_1sigma"].copy()

    Rv = 3.1
    ext = []

    for index in np.arange(np.size(source_table["detectid"])):
        src_wave = np.array([np.double(source_table["wave"][index])])
        ext_i = extinction.fitzpatrick99(src_wave, source_table["Av"][index], Rv)[0]
        ext.append(ext_i)

    deredden = 10 ** (0.4 * np.array(ext))
    deredden[np.isnan(deredden)] = 0
    source_table["flux"] = deredden * source_table["flux_obs"]
    source_table["flux_err"] = deredden * source_table["flux_err_obs"]
    source_table["continuum"] = deredden * source_table["continuum_obs"]
    source_table["continuum_err"] = deredden * source_table["continuum_err_obs"]
    source_table["flux_noise_1sigma"] = deredden * source_table["flux_noise_1sigma_obs"]

    source_table.write("source_catalog_{}.fits".format(args.version), overwrite=True)
    source_table.write(
        "source_catalog_{}.tab".format(args.version), overwrite=True, format="ascii"
    )


if __name__ == "__main__":
    main()
