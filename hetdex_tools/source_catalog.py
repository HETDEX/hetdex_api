import numpy as np
import numpy.ma as ma
import time
import argparse as ap
import os.path as op

import tables as tb

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

from hetdex_api.elixer_widget_cls import ElixerWidget
from hetdex_api.flux_limits.hdf5_sensitivity_cubes import SensitivityCubeHDF5Container

from elixer import catalogs

catlib = catalogs.CatalogLibrary()
config = HDRconfig()

agn_tab = None
cont_gals = None
cont_stars = None

wavelya = 1215.67
waveoii = 3727.8

deth5 = None

def get_flux_noise_1sigma(detid, mask=False, add_apcor_fix=True):

    """ For a detid get the the f50_1sigma value
    
    No need to mask the flux limits if you are using the
    curated catalog which is already masked. This avoid a
    few good sources on the edges of masks that get flagged
    with a high flux limit value.
    """
    
    global config, detect_table, deth5
    
    sncut = 1

    if add_apcor_fix and deth5 is None:
        deth5 = tb.open_file(config.detecth5, 'r')
    
    sel_det = detect_table["detectid"] == detid
    shotid = detect_table["shotid"][sel_det][0]
    ifuslot = detect_table["ifuslot"][sel_det][0]
    
    det_table_here = detect_table[sel_det]
    
    datevobs = str(shotid)[0:8] + "v" + str(shotid)[8:11]
    
    fn = op.join(config.flim_dir, datevobs + "_sensitivity_cube.h5")
    if mask:
        mask_fn = op.join(config.flimmask, datevobs + "_mask.h5")
        sscube = SensitivityCubeHDF5Container(fn, mask_filename=mask_fn, flim_model="hdr1")
    else:
        sscube = SensitivityCubeHDF5Container(fn, flim_model="hdr1")

    if add_apcor_fix:

        if detect_table['det_type'][sel_det][0] == 'cont':
            flim_update = np.nan
            flim = np.nan
            
        else:
            apcor = detect_table['apcor'][sel_det][0]

            if True:
                scube = sscube.extract_ifu_sensitivity_cube("ifuslot_{}".format(ifuslot))
                flim = scube.get_f50(
                    det_table_here["ra"],
                    det_table_here["dec"],
                    det_table_here["wave"],
                    sncut
                )[0]
                
                if apcor >= 0.5:
                    flim_update = flim
                else:
                    fiber_row = deth5.root.Fibers.read_where("detectid == detid")
                    #print(fiber_row)
                    max_idx = np.argmax(fiber_row["weight"])
                    ra = fiber_row['ra'][max_idx]
                    dec = fiber_row['dec'][max_idx]
                    flim_update = scube.get_f50(
                        ra, dec, det_table_here["wave"], sncut
                    )[0]
                
            else:#except Exception:
                flim = np.nan
                flim_update = np.nan

        sscube.close()

        flim_best = np.minimum(flim, flim_update)
        return flim, flim_best

    else:
        try:
            scube = sscube.extract_ifu_sensitivity_cube("ifuslot_{}".format(ifuslot))
            flim = scube.get_f50(
                det_table_here["ra"],
                det_table_here["dec"],
                det_table_here["wave"],
                sncut
            )[0]
            sscube.close()
            return flim
        except Exception:
            return np.nan


def add_elixer_cat_info(det_table, version):

    global config
    
    elixer_file = op.join(config.detect_dir, "catalogs", "elixer_{}_cat.h5".format(version))
    elixer_cat = tb.open_file(elixer_file, "r")

    cls = []
    mlname = []
    mlz = []
    mlprob = []
    best_z = []
    best_pz = []
    flags = []
    
    counterpart_mag = []
    counterpart_mag_err = []
    counterpart_dist = []
    counterpart_catalog_name = []
    counterpart_filter_name = []
    
    fixed_mag = []
    fixed_mag_err = []
    fixed_catalog_name = []
    fixed_filter_name = []
    fixed_radius = []

    for row in det_table:
        detectid_obj = row["detectid"]
        try:
            elix_row = elixer_cat.root.Detections.read_where("detectid == detectid_obj")
            row["plae_classification"] = elix_row["plae_classification"]
            row["combined_plae"] = elix_row["combined_plae"]
            row["combined_plae_err"] = elix_row["combined_plae_err"]
            mlname.append(elix_row["multiline_name"][0].decode())
            cls.append(elix_row["classification_labels"][0].decode())
            mlz.append(elix_row["multiline_z"][0])
            mlprob.append(elix_row["multiline_prob"][0])
            best_z.append(elix_row['best_z'][0])
            best_pz.append(elix_row['best_pz'][0])
            flags.append(elix_row['flags'][0])
            
        except Exception:
            mlname.append("")
            cls.append("")
            mlz.append(False)
            mlprob.append(np.nan)
            best_z.append(np.nan)
            best_pz.append(np.nan)
            flags.append(np.nan)
            
        # append nearest source extracted neighbour match
        try:
        
            elix_row = elixer_cat.root.ExtractedObjects.read_where(
                "(detectid == detectid_obj) & (selected == True)"
            )
            
            if np.size(elix_row) == 0:
                elix_row = elixer_cat.root.ExtractedObjects.read_where(
                    "detectid == detectid_obj"
                )
            if np.size(elix_row) > 1:
                sel_r = elix_row["filter_name"] == b"r"
                if np.sum(sel_r) == 1:
                    counterpart_mag.append(elix_row["mag"][sel_r][0])
                    counterpart_mag_err.append(elix_row["mag_err"][sel_r][0])
                    counterpart_dist.append(elix_row["dist_baryctr"][sel_r][0])
                    counterpart_catalog_name.append(
                        elix_row["catalog_name"][sel_r][0].decode()
                    )
                    counterpart_filter_name.append(
                        elix_row["filter_name"][sel_r][0].decode()
                    )
                else:
                    counterpart_mag.append(elix_row["mag"][0])
                    counterpart_mag_err.append(elix_row["mag_err"][0])
                    counterpart_dist.append(elix_row["dist_baryctr"][0])
                    counterpart_catalog_name.append(
                        elix_row["catalog_name"][0].decode()
                    )
                    counterpart_filter_name.append(elix_row["filter_name"][0].decode())
            elif np.size(elix_row) == 1:
                counterpart_mag.append(elix_row["mag"][0])
                counterpart_mag_err.append(elix_row["mag_err"][0])
                counterpart_dist.append(elix_row["dist_baryctr"][0])
                counterpart_catalog_name.append(elix_row["catalog_name"][0].decode())
                counterpart_filter_name.append(elix_row["filter_name"][0].decode())
            else:
                counterpart_mag.append(np.nan)
                counterpart_mag_err.append(np.nan)
                counterpart_dist.append(np.nan)
                counterpart_catalog_name.append("")
                counterpart_filter_name.append("")
        except:
            counterpart_mag.append(np.nan)
            counterpart_mag_err.append(np.nan)
            counterpart_dist.append(np.nan)
            counterpart_catalog_name.append("")
            counterpart_filter_name.append("")

        # append fixed aperture mag
        try:
            elix_tab = elixer_cat.root.ElixerApertures.read_where(
                ("detectid == detectid_obj")
            )
            sel_r = elix_tab["filter_name"] == b"r"
            sel_g = elix_tab["filter_name"] == b"g"
            
            if np.any(sel_r):
                elix_r = elix_tab[sel_r]
                fixed_mag.append(elix_r["mag"][-1])
                fixed_mag_err.append(elix_r["mag_err"][-1])
                fixed_catalog_name.append(elix_r["catalog_name"][-1].decode())
                fixed_filter_name.append(elix_r["filter_name"][-1].decode())
                fixed_radius.append(elix_r["radius"][-1])
            elif np.any(sel_g):
                elix_g = elix_tab[sel_g]
                fixed_mag.append(elix_g["mag"][-1])
                fixed_mag_err.append(elix_g["mag_err"][-1])
                fixed_catalog_name.append(elix_g["catalog_name"][-1].decode())
                fixed_filter_name.append(elix_g["filter_name"][-1].decode())
                fixed_radius.append(elix_g["radius"][-1])
            else:
                sel = elix_tab["radius"] < 3
                elix_sel = elix_tab[sel]
                fixed_mag.append(elix_sel["mag"][-1])
                fixed_mag_err.append(elix_sel["mag_err"][-1])
                fixed_catalog_name.append(elix_sel["catalog_name"][-1].decode())
                fixed_filter_name.append(elix_sel["filter_name"][-1].decode())
                fixed_radius.append(elix_sel["radius"][-1])
        except:
            fixed_mag.append(np.nan)
            fixed_mag_err.append(np.nan)
            fixed_catalog_name.append("")
            fixed_filter_name.append("")
            fixed_radius.append(np.nan)

    try:
        det_table.add_column(best_z, name='best_z')
        det_table.add_column(best_pz, name='best_pz')
        det_table.add_column(flags, name='flags_elixer')
        det_table.add_column(mlname, name="multiline_name")
        det_table.add_column(cls, name="classification_labels")
        det_table.add_column(counterpart_mag, name="counterpart_mag")
        det_table.add_column(counterpart_mag_err, name="counterpart_mag_err")
        det_table.add_column(counterpart_dist, name="counterpart_dist")
        det_table.add_column(counterpart_catalog_name, name="counterpart_catalog_name")
        det_table.add_column(counterpart_filter_name, name="counterpart_filter_name")
        det_table.add_column(fixed_mag, name="forced_mag")
        det_table.add_column(fixed_mag_err, name="forced_mag_err")
        det_table.add_column(fixed_catalog_name, name="forced_catalog_name")
        det_table.add_column(fixed_filter_name, name="forced_filter_name")
        det_table.add_column(fixed_radius, name="forced_radius")
    except:
        det_table.remove_column('best_z')
        det_table.remove_column('best_pz')
        det_table.remove_column('flags_elixer')
        det_table.remove_column("multiline_name")
        det_table.remove_column("classification_labels")
        det_table.remove_column("counterpart_mag")
        det_table.remove_column("counterpart_mag_err")
        det_table.remove_column("counterpart_dist")
        det_table.remove_column("counterpart_catalog_name")
        det_table.remove_column("counterpart_filter_name")
        det_table.remove_column("forced_mag")
        det_table.remove_column("forced_mag_err")
        det_table.remove_column("forced_catalog_name")
        det_table.remove_column("forced_filter_name")
        det_table.remove_column("forced_radius")
        
        det_table.add_column(best_z, name='best_z')
        det_table.add_column(best_pz, name='best_pz')
        det_table.add_column(flags, name='flags_elixer')
        det_table.add_column(mlname, name="multiline_name")
        det_table.add_column(cls, name="classification_labels")
        det_table.add_column(counterpart_mag, name="counterpart_mag")
        det_table.add_column(counterpart_mag_err, name="counterpart_mag_err")
        det_table.add_column(counterpart_dist, name="counterpart_dist")
        det_table.add_column(counterpart_catalog_name, name="counterpart_catalog_name")
        det_table.add_column(counterpart_filter_name, name="counterpart_filter_name")
        det_table.add_column(fixed_mag, name="forced_mag")
        det_table.add_column(fixed_mag_err, name="forced_mag_err")
        det_table.add_column(fixed_catalog_name, name="forced_catalog_name")
        det_table.add_column(fixed_filter_name, name="forced_filter_name")
        det_table.add_column(fixed_radius, name="forced_radius")

    return det_table


def create_source_catalog(
        version="2.1.3",
        make_continuum=True,
        save=True,
        dsky=4.0):

    global config

    detects_line_table = Table.read('detect_hdr{}.fits'.format(version))
#    detects = Detections(curated_version=version)
#    detects_line_table = detects.return_astropy_table()

#    detects_line_table.write('test.tab', format='ascii')
    detects_line_table.add_column(Column(str("line"), name="det_type", dtype=str))

    detects_cont = Detections(catalog_type="continuum")

    sel1 = detects_cont.remove_bad_amps()
    sel2 = detects_cont.remove_meteors()
    sel3 = detects_cont.remove_shots()
    sel4 = detects_cont.remove_bad_detects()
    sel5 = detects_cont.remove_large_gal()

    sel6 = detects_cont.throughput > 0.08

    detects_cont_table = detects_cont[sel1 * sel2 * sel3 * sel4 * sel5 * sel6].return_astropy_table()
    detects_cont_table.add_column(Column(str("cont"), name="det_type", dtype=str))

    if make_continuum:
        detects_cont_table.write("continuum_" + version + ".fits", overwrite=True)

    dets_all = Detections().refine()
    sel_tp = dets_all.throughput > 0.08
    dets_all_table = dets_all[sel_tp].return_astropy_table()
    dets_all_table.add_column(Column(str("line"), name="det_type", dtype=str))
    agn_tab = Table.read(config.agncat, format="ascii", include_names=['detectid','flux_LyA'])

    # add in continuum sources to match to Chenxu's combined catalog
    #detects_cont_table_orig = detects_cont[sel1 * sel2 * sel3].return_astropy_table()
    
    dets_all_table = vstack([dets_all_table, detects_cont_table])

    detects_broad_table = join(
        agn_tab, dets_all_table, join_type="inner", keys=["detectid"]
    )
    #detects_broad_table.remove_column('det_type')
    #detects_broad_table.add_column(Column(str("agn"), name="det_type", dtype=str))
    dets_all.close()
    del dets_all_table

    global detect_table

    detect_table = unique(
        vstack([detects_broad_table, detects_cont_table, detects_line_table]),
        keys='detectid')

    detect_table.write('test.fits', overwrite=True)

    del detects_cont_table, detects_broad_table

    detect_table = add_elixer_cat_info(detect_table, version)
    
    print('Adding 1sigma noise from flux limits')
    p = Pool(24)
    res = p.map(get_flux_noise_1sigma, detect_table['detectid'])
    p.close()

    flim = []
    flim_update = []
    for r in res:
        flim.append(r[0])
        flim_update.append(r[1])
        
    detect_table['flux_noise_1sigma'] = flim_update
    detect_table['flux_noise_1sigma_orig'] = flim
    
    detect_table.write('test2.fits', overwrite=True)

    print("Performing FOF in 3D space with linking length=10 arcsec")

    # get fluxes to derive flux-weighted distribution of group
    
    detect_table["gmag"][np.isnan(detect_table["gmag"])] = 27
    gmag = detect_table["gmag"] * u.AB
    flux = gmag.to(u.Jy).value

    sel_line = detect_table['wave'] > 0 # remove continuum sources

    # first cluster in positional/wavelegnth space in a larger
    # linking length

    kdtree, r = fof.mktree(
                detect_table["ra"][sel_line],
                detect_table["dec"][sel_line],
                detect_table["wave"][sel_line],
                dsky=6.0, dwave=6.0)

    t0 = time.time()
    print("starting fof ...")
    wfriend_lst = fof.frinds_of_friends(kdtree, r, Nmin=2)
    t1 = time.time()

    wfriend_table = fof.process_group_list(
        wfriend_lst,
        detect_table["detectid"][sel_line],
        detect_table["ra"][sel_line],
        detect_table["dec"][sel_line],
        detect_table["wave"][sel_line],
        detect_table['flux'][sel_line],
    )
    print("Generating combined table \n")
        
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

    wdetfriend_all = join(wdetfriend_tab, wfriend_table, keys="id")

    wdetfriend_all['wave_group_id'] = wdetfriend_all['id'] + 213000000

    wdetfriend_all.rename_column('size', 'wave_group_size')
    wdetfriend_all.rename_column('a', 'wave_group_a')
    wdetfriend_all.rename_column('b', 'wave_group_b')
    wdetfriend_all.rename_column('pa', 'wave_group_pa')
    wdetfriend_all.rename_column('icx', 'wave_group_ra')
    wdetfriend_all.rename_column('icy', 'wave_group_dec')
    wdetfriend_all.rename_column('icz', 'wave_group_wave')

    w_keep = wdetfriend_all['detectid',
                            'wave_group_id',
                            'wave_group_a',
                            'wave_group_b',
                            'wave_group_pa',
                            'wave_group_ra',
                            'wave_group_dec',
                            'wave_group_wave']

    print("3D FOF analysis complete in {:3.2f} minutes \n".format((t1 - t0) / 60))

    print("Performing FOF in 2D space with dlink=3.5 arcsec")

    kdtree, r = fof.mktree(
        detect_table["ra"],
        detect_table["dec"],
        np.zeros_like(detect_table["ra"]),
        dsky=3.5,
    )
    t0 = time.time()
    print("starting fof ...")
    friend_lst = fof.frinds_of_friends(kdtree, r, Nmin=1)
    t1 = time.time()

    print("FOF analysis complete in {:3.2f} minutes \n".format((t1 - t0) / 60))

    friend_table = fof.process_group_list(
                friend_lst,
                detect_table["detectid"],
                detect_table["ra"],
                detect_table["dec"],
                0.0 * detect_table["wave"],
                flux,
            )
 
    print("Generating combined table \n")
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

    detfriend_all = join(detfriend_tab, friend_table, keys="id")

    del detfriend_tab

    # match detectids at large linking length if a wave group exists
    joinfriend = join(detfriend_all, w_keep, keys='detectid', join_type='left')
    grp_by_id = joinfriend.group_by('id')
    sum_grp = grp_by_id.groups.aggregate(np.sum)
    spatial_id = grp_by_id.groups.keys
    spatial_id_to_keep1 = spatial_id[np.isfinite(sum_grp['wave_group_id'])]

    # also match if a group member is brighter than gmag=18 (testing this)
    sel_bright = detect_table['gmag'] < 18
    gmag_join = join( joinfriend, detect_table['detectid','gmag'][sel_bright])
    spatial_id_to_keep2 = gmag_join['id']

    spatial_id_to_keep = vstack([spatial_id_to_keep1, spatial_id_to_keep2])
    detfriend_1 = join(unique(spatial_id_to_keep), joinfriend)

    # link the rest of the detectids with smaller linking length

    keep_row = np.ones(np.size(detect_table), dtype=bool)

    for i, det in enumerate(detect_table['detectid']):
        if det in detfriend_1['detectid']:
            keep_row[i] = 0

    print("Performing FOF in 2D space with dlink=3.0 arcsec")

    kdtree, r = fof.mktree(
        detect_table["ra"][keep_row],
        detect_table["dec"][keep_row],
        np.zeros_like(detect_table["ra"][keep_row]),
        dsky=3.0,
    )
        
    t0 = time.time()
    print("starting fof ...")
    friend_lst = fof.frinds_of_friends(kdtree, r, Nmin=1)
    t1 = time.time()

    print("Final FOF analysis complete in {:3.2f} minutes \n".format((t1 - t0) / 60))

    friend_table = fof.process_group_list(
            friend_lst,
            detect_table["detectid"][keep_row],
            detect_table["ra"][keep_row],
            detect_table["dec"][keep_row],
            0.0 * detect_table["wave"][keep_row],
            flux[keep_row],
        )

    print("Generating combined table \n")
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
     
    detfriend_2 = join(detfriend_tab, friend_table, keys="id")
     
    starting_id_1 = int(version.replace('.', '', 2))*10**10
    starting_id_2 = starting_id_1 + 10**8
    detfriend_1.add_column(
        Column(detfriend_1["id"] + starting_id_1, name="source_id"), index=0
    )
    detfriend_2.add_column(
        Column(detfriend_2["id"] + starting_id_2, name="source_id"), index=0
    )
    detfriend_1.remove_column("id")
    detfriend_2.remove_column("id")

    detfriend_all = vstack([detfriend_1, detfriend_2])
    expand_table = join(detfriend_all, detect_table, keys="detectid")
    expand_table['wave_group_id'] = expand_table['wave_group_id'].filled(0)
    expand_table.write('test3.fits', overwrite=True)

    # combine common wavegroups to the same source_id

    sel = expand_table['wave_group_id'] > 0
    s_by_id = expand_table[sel].group_by('wave_group_id')

    for grp in s_by_id.groups:
        sid, ns = np.unique(grp['source_id'], return_counts=True)
        wid = grp['wave_group_id'][0]
        if np.size(sid) > 1:
            sel_wid = expand_table['wave_group_id'] == wid
            expand_table['source_id'][sel_wid] = sid[0]

    del detfriend_all, detect_table, friend_table

    gaia_stars = Table.read(config.gaiacat)

    gaia_coords = SkyCoord(ra=gaia_stars["ra"] * u.deg, dec=gaia_stars["dec"] * u.deg)
    src_coords = SkyCoord(
        ra=expand_table["ra"] * u.deg, dec=expand_table["dec"] * u.deg
    )

    idx, d2d, d3d = src_coords.match_to_catalog_sky(gaia_coords)

    sel = d2d < 5.0 * u.arcsec

    gaia_match_name = np.zeros_like(expand_table["source_id"], dtype=int)
    gaia_match_name[sel] = gaia_stars["source_id"][idx][sel]

    gaia_match_dist = np.zeros_like(expand_table["source_id"], dtype=float)
    gaia_match_dist[sel] = d2d[sel].to_value(u.arcsec)
    
    expand_table["gaia_match_id"] = gaia_match_name
    expand_table["gaia_match_dist"] = gaia_match_dist
    
    expand_table.rename_column("size", "n_members")
    expand_table.rename_column("icx", "ra_mean")
    expand_table.rename_column("icy", "dec_mean")

    expand_table.sort("source_id")

    return expand_table


def z_guess_3727(group, cont=False):
    sel_good_lines = None
    if np.any((group["sn"] > 15) * (group["wave"] > waveoii)):
        sel_good_lines = (group["sn"] > 15) * (group["wave"] > waveoii)
    elif np.any((group["sn"] > 10) * (group["wave"] > waveoii)):
        sel_good_lines = (group["sn"] > 8) * (group["wave"] > waveoii)
    elif np.any((group["sn"] > 8) * (group["wave"] > waveoii)):
        sel_good_lines = (group["sn"] > 8) * (group["wave"] > waveoii)
    elif np.any((group["sn"] > 6) * (group["wave"] > waveoii)):
        if cont:
            pass
        else:
            sel_good_lines = (group["sn"] > 6) * (group["wave"] > waveoii)
    else:
        if cont:
            pass
        else:
            sel_good_lines = group["wave"] > waveoii

    if sel_good_lines is not None:
        try:
            wave_guess = np.min(group["wave"][sel_good_lines])
            z_guess = wave_guess / waveoii - 1
        except Exception:
            z_guess = -3.0
    else:
        z_guess = -2.0

    return z_guess


def guess_source_wavelength(source_id):

    global source_table, agn_tab, cont_gals, cont_stars, config

    if agn_tab is None:
        agn_tab = Table.read(config.agncat, format="ascii")
    if cont_gals is None:
        cont_gals = np.loadtxt(config.galaxylabels, dtype=int)
    if cont_stars is None:
        cont_stars = np.loadtxt(config.starlabels, dtype=int)
        
    sel_group = source_table["source_id"] == source_id
    group = source_table[sel_group]
    z_guess = -1.0
    s_type = "none"
    z_flag = -1

    # Check if any member is an AGN first
    agn_flag = False
    for det in group['detectid']:
        if det in agn_tab['detectid']:
            agn_flag = True
            agn_det = det

    if agn_flag:
        # get proper z's from Chenxu's catalog
        sel_det = agn_tab["detectid"] == agn_det
        z_guess = agn_tab["zguess"][sel_det][0]
        z_flag = agn_tab['zflag'][sel_det][0]
        s_type = "agn"

    elif np.any(group["det_type"] == "cont"):
        sel_cont = group['det_type'] == "cont"
        det = group['detectid'][sel_cont]
        if det in cont_gals:
            if np.any(group['det_type']=="line"):
                z_guess = z_guess_3727(group, cont=True)
                if z_guess > 0:
                    s_type = "oii"
                else:
                    s_type = "unsure-cont-line"
            else:
                s_type = 'lzg'
                z_guess = -4.0

        elif det in cont_stars:
            if np.any(group['det_type']=="line"):
                z_guess = z_guess_3727(group, cont=True)
                if z_guess > 0:
                    s_type = "oii"
                else:
                    s_type = "star"
                    z_guess = 0.0
            else:
                z_guess = 0.0
                s_type = "star"

        elif np.any(group["gaia_match_id"] > 0):
            z_guess = z_guess_3727(group, cont=True)
            if z_guess > 0:
                s_type = "oii"
            else:
                s_type = "star"
                z_guess = 0.0
        else:
            if np.any(group['det_type'] == 'line'):
                z_guess = z_guess_3727(group, cont=True)
                if z_guess > 0:
                    s_type = "oii"
                else:
                    s_type = "unsure-cont-line"
            else:
                s_type = "unsure-cont"
                z_guess = -5.0
                
    elif np.any(group["plae_classification"] < 0.1):
        z_guess = z_guess_3727(group)
        if z_guess > 0:
            s_type = "oii"
        else:
            s_type = "unsure-line"

    elif (np.nanmedian(group["plae_classification"]) < 0.4):
        z_guess = z_guess_3727(group)
        if z_guess > 0:
            s_type = "oii"
        else:
            s_type = "unsure-line"

    elif (np.nanmedian(group["plae_classification"]) >= 0.4):
        if np.any(group["sn"] > 15):
            sel_good_lines = group["sn"] > 15
            wave_guess = np.min(group["wave"][sel_good_lines])
        elif np.any(group["sn"] > 10):
            sel_good_lines = group["sn"] > 10
            wave_guess = np.min(group["wave"][sel_good_lines])
        elif np.any(group["sn"] > 8):
            sel_good_lines = group["sn"] > 8
            wave_guess = np.min(group["wave"][sel_good_lines])
        elif np.any(group["sn"] > 6.5):
            sel_good_lines = group["sn"] > 6.5
            wave_guess = np.min(group["wave"][sel_good_lines])
        elif np.any(group["sn"] > 5.5):
            sel_good_lines = group["sn"] > 5.5
            wave_guess = np.min(group["wave"][sel_good_lines])
        else:
            argmaxsn = np.argmax(group["sn"])
            wave_guess = group["wave"][argmaxsn]
        z_guess = wave_guess / wavelya - 1
        s_type = "lae"
    else:
        argmaxsn = np.argmax(group["sn"])
        wave_guess = group["wave"][argmaxsn]
        z_guess = wave_guess / wavelya - 1
        s_type = "lae"
        
    return z_guess, z_flag, str(s_type),


def add_z_guess(source_table):

    try:
        # remove z_guess column if it exists
        print("Removing existing z_guess column")
        source_table.remove_column("z_guess")
    except Exception:
        pass

    try:
        source_table.remove_column("source_type")
    except Exception:
        pass
    print("Assigning a best guess redshift")

    from multiprocessing import Pool

    src_list = np.unique(source_table["source_id"])

    #z_guess = []
    #s_type = []
    #for src in src_list:
    #    z_g, s_t = guess_source_wavelength(src)
    #    z_guess.append(z_g)
    #    s_type.append(s_t)
        
    t0 = time.time()
    p = Pool(6)
    res = p.map(guess_source_wavelength, src_list)
    t1 = time.time()
    z_guess = []
    s_type = []
    z_flag = []
    for r in res:
        z_guess.append(r[0])
        z_flag.append(r[1])
        s_type.append(r[2])
    p.close()

    print("Finished in {:3.2f} minutes".format((t1 - t0) / 60))

    z_table = Table(
        [src_list, z_guess, z_flag, s_type],
        names=["source_id", "z_guess", "z_agn_flag", "source_type"]
    )

    all_table = join(source_table, z_table, join_type="left")
    del source_table, z_table
    
    return all_table


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
    #    z_guess = guess_source_wavelength(source_id, source_table)
    z_guess = group["z_guess"][0]

    plt.title(
        "source_id:%d n:%d ra:%6.3f dec:%6.3f z:%4.3f"
        % (
            source_id,
            group["n_members"][0],
            group["ra_mean"][0],
            group["dec_mean"][0],
            z_guess,
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
    """ function that returns a parser from argparse """

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
        "-dsky",
        "--dsky",
        type=float,
        help="""Spatial linking length in arcsec""",
        default=4.0,
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
    """ Main Function """

    parser = get_parser()
    args = parser.parse_args(argv)

    print('Combining catalogs')
    global source_table
    print(args.dsky, args.version)

    source_table = create_source_catalog(version=args.version, dsky=args.dsky)
    source_table.write('source_cat_tmp.fits', overwrite=True)
#    source_table = Table.read('source_cat_tmp.fits')
    # add source name
    source_name = []
    for row in source_table:
        source_name.append(get_source_name(row["ra_mean"], row["dec_mean"]))
    try:
        source_table.add_column(source_name,
                                name="source_name",
                                index=1)
    except:
        print('messed up source name again')

    # match to SDSS
    source_table_coords = SkyCoord(source_table['ra_mean'],
                                   source_table['dec_mean'],
                                   unit='deg')
    sdssfile = op.join( config.imaging_dir, 'sdss', 'specObj-dr16-trim.fits')
    sdsstab = Table.read(sdssfile)
    sdsstab = Table.read( sdssfile)
    sdss_coords = SkyCoord(ra=sdsstab['PLUG_RA'], dec=sdsstab['PLUG_DEC'], unit='deg')
    idx, d2d, d3d = source_table_coords.match_to_catalog_sky(sdss_coords)
    
    catalog_matches = sdsstab[idx]
    catalog_matches['sdss_dist'] = d2d.to_value(u.arcsec)

    catalog_matches.rename_column('PLUG_RA', 'ra_sdss')
    catalog_matches.rename_column('PLUG_DEC', 'dec_sdss')
    catalog_matches.rename_column('CLASS', 'sdss_class')
    catalog_matches.rename_column('Z', 'z_sdss')
    catalog_matches.rename_column('Z_ERR', 'z_sdss_err')
    
    matched_catalog = hstack([source_table, catalog_matches])

    sel_remove = matched_catalog['sdss_dist'] > 5

    matched_catalog['ra_sdss'][sel_remove] = np.nan
    matched_catalog['dec_sdss'][sel_remove] = np.nan
    matched_catalog['sdss_dist'][sel_remove] = np.nan
    matched_catalog['sdss_class'][sel_remove] = ''
    matched_catalog['z_sdss'][sel_remove] = np.nan
    matched_catalog['z_sdss_err'][sel_remove] = np.nan

    source_table = matched_catalog
                             
    # match band-merged WISE catalog

    wise_catalog = Table.read('wise-hetdexoverlap.fits')
    source_table_coords = SkyCoord(source_table['ra_mean'],
                                   source_table['dec_mean'],
                                   unit='deg')
    wise_coords = SkyCoord(ra=wise_catalog['ra'], dec=wise_catalog['dec'], unit='deg')
    idx, d2d, d3d = source_table_coords.match_to_catalog_sky(wise_coords)

    catalog_matches = wise_catalog[idx]
    catalog_matches['wise_dist'] = d2d.to_value(u.arcsec)
    
    keep_wise = catalog_matches['ra','dec','primary','unwise_objid','flux', 'wise_dist']
    keep_wise.rename_column('flux','wise_fluxes')
    keep_wise.rename_column('ra','ra_wise')
    keep_wise.rename_column('dec','dec_wise')

    matched_catalog = hstack([source_table, keep_wise])

    w1 = []
    w2 = []
    for row in matched_catalog:
            w1.append(row['wise_fluxes'][0])
            w2.append(row['wise_fluxes'][1])

    matched_catalog['flux_w1'] = w1
    matched_catalog['flux_w2'] = w2
    matched_catalog.remove_column('wise_fluxes')
    sel_close = matched_catalog['wise_dist'] < 5 #arcsec
    
    print('There are {} wise matches'.format(np.size(np.unique(matched_catalog['source_id'][sel_close]))))
    # remove column info for WISE matches more than 5 arcsec away

    sel_remove = matched_catalog['wise_dist'] > 5 #arcsec

    matched_catalog['ra_wise'][sel_remove] = np.nan
    matched_catalog['dec_wise'][sel_remove] = np.nan
    matched_catalog['wise_dist'][sel_remove] = np.nan
    matched_catalog['primary'][sel_remove] = -1
    matched_catalog['unwise_objid'][sel_remove] = np.nan
    matched_catalog['flux_w1'][sel_remove] = np.nan
    matched_catalog['flux_w2'][sel_remove] = np.nan
    
    source_table = matched_catalog

    # add z_spec from other catlogs if it exists:
    goods_z = Table.read('catalogs/goods_n_specz_1018_no_header.txt',
                         names=['ra_zspec','dec_zspec','zspec',
                                'z_quality','zspec_catalog','Symbol'],
                         format='ascii.no_header')

    #DEIMOS 10k (Hasinger et al. 2018) z_spec up to ~6
    #https://cosmos.astro.caltech.edu/news/65
    deimos = Table.read('catalogs/deimos_redshifts.tbl', format='ascii')
    deimos.rename_column('Ra', 'ra_zspec')
    deimos.rename_column('Dec', 'dec_zspec')
    deimos['zspec_catalog'] = 'CosmosDeimos'

    #Kriek et al. (2015)
    #http://mosdef.astro.berkeley.edu/for-scientists/data-releases/
    mosdef = Table.read('catalogs/mosdef_zcat.final_slitap.fits')
    mosdef.rename_column('RA','ra_zspec')
    mosdef.rename_column('DEC','dec_zspec')
    mosdef.rename_column('Z_MOSFIRE', 'zspec')
    mosdef['zspec_catalog'] = 'MOSFIRE'

    #VUDS (Tasca et al. 2017), z_spec up to ~6
    #http://cesam.lam.fr/vuds/DR1/
    zcosbright = Table.read('catalogs/cesam_zcosbrightspec20k_dr3_catalog_1616073679.txt', format='ascii')
    zcosbright.rename_column('zpec','zspec')
    zcosbright.rename_column('ra','ra_zspec')
    zcosbright.rename_column('dec','dec_zspec')
    zcosbright['zspec_catalog'] = 'zCosmosBright'

    deep_specz = Table.read('catalogs/DEEP_zcosmos_spectroscopy_one_v2.6_data.cat',
                            format='ascii', data_start=100)
    deep_specz.rename_column('col1','zCOSMOS-deepID')
    deep_specz.rename_column('col2','zspec')
    deep_specz.rename_column('col3','flag')
    deep_specz.rename_column('col4','zphot')
    deep_specz.rename_column('col5','ra_zspec')
    deep_specz.rename_column('col6','dec_zspec')
    deep_specz['zspec_catalog'] = 'DEEP_zcosmos'

    sdssfile = op.join( config.imaging_dir, 'sdss', 'specObj-dr16-trim.fits')

    sdss_specz = Table.read(sdssfile)
    sdss_specz.rename_column('PLUG_RA', 'ra_zspec')
    sdss_specz.rename_column('PLUG_DEC', 'dec_zspec')
    sdss_specz.rename_column('CLASS', 'sdss_class')
    sdss_specz.rename_column('Z', 'zspec')
    sdss_specz.rename_column('Z_ERR', 'z_sdss_err')
    sdss_specz['zspec_catalog'] = 'sdss-dr16'

    specz_catalogs = vstack([goods_z['zspec','ra_zspec','dec_zspec','zspec_catalog'],
                             deimos['zspec','ra_zspec','dec_zspec','zspec_catalog'],
                             mosdef['zspec','ra_zspec','dec_zspec','zspec_catalog'],
                             zcosbright['zspec','ra_zspec','dec_zspec','zspec_catalog'],
                             deep_specz['zspec','ra_zspec','dec_zspec','zspec_catalog'],
                             sdss_specz['zspec','ra_zspec','dec_zspec','zspec_catalog'],
                             ])
    sel=specz_catalogs['zspec'] >= 0
    specz_catalogs = specz_catalogs[sel]
    
    specz_coords = SkyCoord(ra=specz_catalogs['ra_zspec'],
                            dec=specz_catalogs['dec_zspec'],
                            unit='deg')
    source_coords = SkyCoord(ra=source_table['ra'],
                             dec=source_table['dec'],
                             unit='deg')
            
    idx, d2d, d3d = source_coords.match_to_catalog_sky(specz_coords)
    
    catalog_matches = specz_catalogs[idx]
    catalog_matches['zspec_dist'] = d2d.to_value(u.arcsec)

    matched_catalog = hstack([source_table, catalog_matches])

    sel_close = matched_catalog['zspec_dist'] < 5 #u.arcsec

    print('There are {} zspec matches within 5 arcsec'.format(np.size(np.unique(matched_catalog['source_id'][sel_close]))))

    sel_remove = matched_catalog['zspec_dist'] > 5 #u.arcsec

    matched_catalog['zspec'][sel_remove] = np.nan
    matched_catalog['ra_zspec'][sel_remove] = np.nan
    matched_catalog['dec_zspec'][sel_remove] = np.nan
    matched_catalog['zspec_dist'][sel_remove] = np.nan
    matched_catalog['zspec_catalog'][sel_remove] = ''

    #add desi confirmed redshifts

    dtab = Table.read('catalogs/desi-hetdex.fits')

    sel_good_hetdex = (dtab['artifact_flag'] == 0) * (dtab['wave'] >= 3610 )
    sel_good_desi = (dtab['FIBERSTATUS'] == 0)
    sel_sample = sel_good_desi*sel_good_hetdex
    sel_conf = dtab['VI_quality'] >= 1
   
    hetdex_coords = SkyCoord(ra=dtab['ra'], dec=dtab['dec'], unit='deg')
    desi_coords = SkyCoord(ra=dtab['TARGET_RA'], dec=dtab['TARGET_DEC'], unit='deg')
    dtab['zspec_dist'] = hetdex_coords.separation(desi_coords).arcsec

    zspec = []
    zspec_dist = []

    desi_matches = dtab['detectid'][sel_sample*sel_conf]
    for row in dtab[sel_sample*sel_conf]:
        if row['VI_z'] > 1.9:
            wave_z = (1 + row['VI_z'])*wavelya
            if (np.abs(wave_z - row['wave']) < 10):
                zspec.append(row['VI_z'])
                zspec_dist.append(row['zspec_dist'])
            else:
                zspec.append(np.nan)
                zspec_dist.append(np.nan)

        elif row['VI_z'] < 0.5:
            wave_z = (1 + row['VI_z'])*waveoii
            if (np.abs(wave_z - row['wave']) < 10):
                zspec.append(row['VI_z'])
                zspec_dist.append(row['zspec_dist'])
            else:
                zspec.append(np.nan)
                zspec_dist.append(np.nan)
        else:
            zspec.append(np.nan)
            zspec_dist.append(np.nan)

    for i in np.arange(len(desi_matches)):
        sel_det = matched_catalog['detectid'] == desi_matches[i]
        matched_catalog['zspec'][sel_det] = zspec[i]
        matched_catalog['zspec_dist'][sel_det] = zspec_dist[i]
        matched_catalog['zspec_catalog'][sel_det] = 'DESI'

    source_table = matched_catalog
   
    # Clear up memory

    for name in dir():
        if source_table:
            continue
        elif not name.startswith('_'):
            del name
            
    import gc
    gc.collect()
    
    # add a guess on redshift and source type
    out_table = add_z_guess(source_table)

    #names=['zCOSMOS-deepID','zspec','flag','zphot','ra','dec'],
    
    sel_star = out_table['source_type'] == 'star'
    sel_oii = out_table['source_type'] == 'oii'
    sel_lae = out_table['source_type'] == 'lae'
    sel_agn = out_table['source_type'] == 'agn'

    print('There are {} stars, {} OII emitters and {} LAEs'.format(np.sum(sel_star), np.sum(sel_oii), np.sum(sel_lae)))

    # sort table by source wavelength closest to z_guess and position
    # so unique will produce the closest match    
    src_coord = SkyCoord(ra=out_table['ra_mean'], dec=out_table['dec_mean'], unit='deg')
    det_coord = SkyCoord(ra=out_table['ra'], dec=out_table['dec'], unit='deg')

    src_wave = np.zeros_like(out_table['z_guess'])
    src_wave[sel_oii] = (1 + out_table['z_guess'][sel_oii]) * waveoii
    src_wave[sel_lae] = (1 + out_table['z_guess'][sel_lae]) * wavelya

    sel_z = ( out_table['z_guess'] >= 1.9) * (out_table['z_guess'] <= 3.6)
    src_wave[sel_agn*sel_z] = (1 + out_table['z_guess'][sel_agn*sel_z]) * wavelya

    out_table['src_separation'] = det_coord.separation(src_coord)
    out_table['dwave'] = np.abs(src_wave - source_table['wave'])

    out_table.sort(['dwave', 'src_separation'])

    print('Filling masked values with NaNs')
    
    for col in out_table.columns:
        try:
            out_table[col] = out_table[col].filled(np.nan)
            print('yes', col)
        except:
            print('no', col)

    out_table.write("source_catalog_{}.fits".format(args.version),
                    overwrite=True)
    out_table.write("source_catalog_{}.tab".format(args.version),
                    format='ascii', overwrite=True)


if __name__ == "__main__":
    main()
