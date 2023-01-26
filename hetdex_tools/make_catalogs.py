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
config = HDRconfig('hdr3')

wavelya = 1215.67
waveoii = 3727.8

deth5 = None
conth5 = None
add_agn = True

def make_friend_table_for_shot(shotid):

    global detect_table, dsky_2D
    sel_shot = detect_table['shotid'] == shotid
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
        detect_table['flux_g'][sel_shot],
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
    
    sel_shot = (detect_table['shotid'] == shotid) & (detect_table['det_type'] == 'line')
    kdtree, r = fof.mktree(
        detect_table["ra"][sel_shot],
        detect_table["dec"][sel_shot],
        detect_table["wave"][sel_shot],
        dsky=dsky_3D, dwave=dwave)
    
    wfriend_lst = fof.frinds_of_friends(kdtree, r, Nmin=2)

    if len(wfriend_lst) > 0:
        
        wfriend_table = fof.process_group_list(
            wfriend_lst,
            detect_table["detectid"][sel_shot],
            detect_table["ra"][sel_shot],
            detect_table["dec"][sel_shot],
            detect_table["wave"][sel_shot],
            detect_table['flux'][sel_shot],
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


def return_fiber_ratio(det, det_type):

    global deth5, conth5

    if deth5 is None:
        deth5 = tb.open_file(config.detecth5, 'r')
    if conth5 is None:
        conth5 = tb.open_file(config.contsourceh5, 'r')
    
    if det_type == 'line':
        fiber_tab = deth5.root.Fibers.read_where("detectid == det")
        weights = np.sort(fiber_tab["weight"])
        fiber_ratio = weights[-1] / weights[-2]
    else:
        fiber_tab = conth5.root.Fibers.read_where("detectid == det")
        weights = np.sort(fiber_tab["weight"])
        fiber_ratio = weights[-1] / weights[-2]
    return fiber_ratio

                                                                        
def add_elixer_cat_info(detect_table, det_type='line'):

    global config

    if det_type == 'line':
        elixer_file = op.join(config.detect_dir, "elixer_hdr3_emis_cat.h5")
        elixer_cat = tb.open_file(elixer_file, "r")
    else:
        elixer_file = op.join(config.detect_dir, 'elixer_hdr3_cont_cat.h5')
        elixer_cat = tb.open_file(elixer_file, 'r')

    elix_tab = Table(elixer_cat.root.Detections.read())

    sel_det_col = ['detectid',
                   'mag_g_wide',
                   'plya_classification',
                   'best_z',
                   'best_pz',
                   'combined_plae',
                   'combined_plae_lo',
                   'combined_plae_hi',
                   'flags',
                   'classification_labels',
    ]

    catalog = join(detect_table, unique( elix_tab[sel_det_col], keys='detectid'),  keys='detectid', join_type='left')

    catalog.rename_column('mag_g_wide', 'gmag')

    extracted_objects = Table(elixer_cat.root.ExtractedObjects.read())
    selected1 = (extracted_objects['selected'] == True) & (extracted_objects['filter_name'] == b'r')
    selected2 = (extracted_objects['catalog_name'] == b'HSC-DEX') \
                | (extracted_objects['catalog_name'] == b'HSC-SSP')
    selected = selected1 & selected2

    eo_tab = unique( extracted_objects[selected], keys='detectid')

    eo_tab.rename_column('flags','flags_eo')
    eo_tab.rename_column('ra', 'counterpart_ra')
    eo_tab.rename_column('dec', 'counterpart_dec')
    eo_tab.rename_column('mag', 'counterpart_mag')
    eo_tab.rename_column('mag_err', 'counterpart_mag_err')
    eo_tab.rename_column('dist_baryctr', 'counterpart_dist')
    eo_tab.rename_column('catalog_name','counterpart_catalog_name')
    eo_tab.rename_column('filter_name', 'counterpart_filter_name')
    eo_tab.remove_column('flux_cts')
    eo_tab.remove_column('flux_err')
    
    catalog2 = join(catalog, eo_tab, keys='detectid', join_type='left')

    elixer_cat.close()

    catalog2.rename_column('best_z','z_elixer')
    
    # fill mask values with nans
    for col in catalog2.columns:
        try:
            catalog2[col] = catalog2[col].filled(np.nan)
            print('Replacing masked values with np.nan for {}'.format(col))
        except:
            pass
            #print('no', col)   
    return catalog2


def merge_wave_groups(wid):

    global expand_table
    
    try:
        sel_wid = expand_table['wave_group_id'] == wid
        grp = expand_table[sel_wid]
        sid, ns = np.unique(grp['source_id'], return_counts=True)
        
        sid_main = np.min(sid)
        
        sid_ind = []
        for sid_i in sid:
            for ind in np.where( expand_table['source_id'] == sid_i)[0]:
                sid_ind.append(ind)
        
        # now find any other wave groups and their associated source_id info to merge
        other_wids = np.unique(expand_table['wave_group_id'][sid_ind])
        
        for wid_i in other_wids:
            if wid_i == 0:
                continue
            elif wid_i == wid:
                continue
                
            sel_wid_i = expand_table['wave_group_id'] == wid_i
            grp = expand_table[sel_wid_i]
            sid, ns = np.unique(grp['source_id'], return_counts=True)
            
            for sid_i in sid:
                if sid_i == sid_main:
                    continue
                for ind in np.where( expand_table['source_id'] == sid_i)[0]:
                    sid_ind.append(ind)
                    
        if np.size( np.unique( expand_table['source_id'][sid_ind])) > 1:
            return sid_ind
        else:
            return None
    except Exception:
        print('Merge wave group failed for {}'.format(wid))
        return None
        

def create_source_catalog(
        version="3.0.0",
        update=False):

    global config

    if update:

        print('Creating curated detection catalog version={}'.format(version))
        print('Updating chi2fib values')
        # Note chi2fib calculated in
        #/work/05350/ecooper/stampede2/hdr3/fiber_chi2_check-hdr3.ipynb
        chi2fib_tab = Table.read('/work/05350/ecooper/stampede2/hdr3/chi2fib_3.0.0.txt', format='ascii')
        
        D = Detections(survey='hdr3', loadtable=True)
        sel_cut1 = (D.sn >= 7) & (D.chi2 <= 2.5)
        sel_cut2 = (D.sn >= 4.8) & (D.sn < 7) & (D.chi2 <= 1.2)
        
        sel_cont = D.continuum > -3
        sel_chi2fib = D.chi2fib < 4
        sel_tp = D.throughput >= 0.08
        
        sel = sel_cont & sel_chi2fib & sel_tp & (sel_cut1 | sel_cut2) #&sel_field
    
        sel_wave = (D.wave >= 3550) & (D.wave <= 5460)
        sel_lw = (D.linewidth <= 14) & (D.linewidth > 6) & (D.sn >= 6.5)

        sel1 = sel & sel_wave & sel_lw
       
        sel_wave = (D.wave >= 3510) & (D.wave <= 5490)
        sel_lw = (D.linewidth <= 6) & (D.linewidth>= 1.7)
       
        sel2 = sel & sel_wave & sel_lw
        
        sel_cat = sel1 | sel2

        # downselect catalog criteria
        D = D[sel_cat]
        mask_badamp = D.remove_bad_amps()
        mask_badshots = D.remove_shots()

        # downselect badamps and badshots
        D = D[mask_badamp & mask_badshots]
        
        mask_baddet = D.remove_bad_detects()
        mask_badpix = D.remove_bad_pix()
        mask_meteor = D.remove_meteors()
        galmask = D.remove_large_gal(d25scale=3.0)

        print(len(D))
        detects_line_table = D.return_astropy_table()
        detects_line_table.add_column(Column(str("line"), name="det_type", dtype=str))
        print(len(detects_line_table))
        detects_line_table.add_column(Column(mask_baddet.astype(int), name="flag_baddet", dtype=int))
        detects_line_table.add_column(Column(mask_badpix.astype(int), name="flag_badpix", dtype=int))
        detects_line_table.add_column(Column(mask_meteor.astype(int), name="flag_meteor", dtype=int))
        detects_line_table.add_column(Column(galmask.astype(int), name="flag_gal", dtype=int))

        detects_line_table = add_elixer_cat_info(detects_line_table, det_type='line')
        print(len(detects_line_table))
        #apply additional chi2fib cut
        detects_line_table2 = join(detects_line_table, chi2fib_tab, keys='detectid', join_type='left')
        detects_line_table3 = unique(detects_line_table2, keys='detectid')
        
        print('Size after combining with chi2fib_tab: {}'.format(len(detects_line_table3)))

        sel_chi2fib = (detects_line_table3['chi2fib_1'] < 4) & (detects_line_table3['chi2fib_2'] < 4)
        detects_line_table = detects_line_table3[sel_chi2fib]
        print('Size after chi2fib cut: {}'.format(len(detects_line_table)))
        detects_line_table.write('detect_hdr{}.fits'.format(version), overwrite=True)
        detects_line_table.write('detect_hdr{}.tab'.format(version),
                                 format='ascii',
                                 overwrite=True)
    else:
        detects_line_table = Table.read('detect_hdr{}.fits'.format(version))

    if update:
        detects_cont = Detections(catalog_type="continuum", survey='hdr3', loadtable=True)
        
        sel1 = detects_cont.remove_bad_amps()
        sel2 = detects_cont.remove_meteors()
        sel3 = detects_cont.remove_shots()
        sel4 = detects_cont.remove_bad_detects()
        sel5 = detects_cont.remove_large_gal()
        
        sel6 = detects_cont.throughput > 0.08
        
        detects_cont_table = detects_cont[sel1 & sel2 & sel3 & sel4 & sel5 & sel6].return_astropy_table()
        detects_cont_table.add_column(Column(str("cont"), name="det_type", dtype=str))
        print(len(detects_cont_table))
        detects_cont_table = add_elixer_cat_info(detects_cont_table, det_type='cont')
        print(len(detects_cont_table))
        detects_cont_table.write("continuum_" + version + ".fits", overwrite=True)
        detects_cont_table.write("continuum_" + version + ".tab", overwrite=True, format='ascii')
        detects_cont.close()
    else:
        detects_cont_table = Table.read("continuum_" + version + ".fits")

    if update:
        # create an agn table with detection info added
        if add_agn:
            full_line_table = D.return_astropy_table()
            full_line_table.add_column(Column(str("line"), name="det_type", dtype=str))
            agn_tab = Table.read(config.agncat, format="ascii", include_names=['detectid','vis_class','z'])
            agn_tab.rename_column('vis_class', 'agn_vis_class')
            agn_tab.rename_column('z','z_agn')
            detects_agn = join(agn_tab, vstack([detects_cont_table, full_line_table]), join_type='inner')
            detects_agn.write("agn_" + version + ".fits", overwrite=True)
        else:
            print('No updated AGN catalog created')
        D.close()
    else:
        detects_agn = Table.read("agn_" + version + ".fits")

    if update:
        return
    
    global detect_table

    if add_agn:
        detect_table = unique(
            vstack([detects_agn, detects_cont_table, detects_line_table]),
            keys='detectid')
    else:
        
        detect_table = unique(
            vstack([detects_cont_table, detects_line_table]),
            keys='detectid')

    append_fiber_ratio = False
    if append_fiber_ratio:
        # add fiber_ratio
        fiber_ratio = []
        for row in detect_table:
            det = row['detectid']
            det_type = row['det_type']
            try:
                fiber_ratio.append( return_fiber_ratio(det, det_type))
            except:
                fiber_ratio.append(np.nan)
                print('fiber_ratio failed for {}'.format(det))
        detect_table['fiber_ratio'] = fiber_ratio

    del detects_cont_table, detects_line_table

    gc.collect()
    
    # calculate ebv and av for every detections
    all_coords = SkyCoord(ra=detect_table['ra'], dec=detect_table['dec'], unit='deg')
    sfd = SFDQuery()
    ebv = sfd(all_coords)
    Rv = 3.1
    corr_SF2011 = 2.742  # Landolt V
    ext = []

    Av = corr_SF2011*ebv

    detect_table['Av'] = Av
    detect_table['ebv'] = ebv

    # get fluxes to derive flux-weighted distribution of group
    
    detect_table["gmag"][np.isnan(detect_table["gmag"])] = 27
    gmag = detect_table["gmag"] * u.AB
    detect_table['flux_g'] = gmag.to(u.Jy).value

    # first cluster in positional/wavelegnth space in a larger
    # linking length

    shotlist = np.unique(detect_table['shotid'])

    print("Performing FOF in 3D space wiht dlink={} and dwave={}".format(dsky_3D, dwave))
    t0 = time.time()
    P = Pool(24)
    res = P.map(make_wfriend_table_for_shot, shotlist)
    P.close()
    
    wstart = int(version.replace('.', '', 2))*10**6
    wdetfriend_all = Table()

    firstid = True
    for r in res:

        if r is None:
            continue
        else:
            if firstid:
                r['id'] = r['id'] + wstart
                firstid = False
            else:
                r['id'] = r['id'] + np.max(wdetfriend_all['id']) + 1
        wdetfriend_all = vstack([wdetfriend_all, r])
        
    wdetfriend_all.rename_column('id','wave_group_id')
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

    t1 = time.time()
    
    print("3D FOF analysis complete in {:3.2f} minutes \n".format((t1 - t0) / 60))

    print("Performing FOF in 2D space with dlink={}".format(dsky_2D))

    t0 = time.time()
    P = Pool(24)
    res = P.map(make_friend_table_for_shot, shotlist)
    P.close()

    sid_start = int(version.replace('.', '', 2))*10**10

    firstid = True
    detfriend_all = Table()
    for r in res:

        if firstid:
            r['id'] = r['id'] + sid_start
            firstid = False
        else:
            r['id'] = r['id'] + np.max(detfriend_all['id']) + 1
        detfriend_all = vstack([detfriend_all, r])

    t1 = time.time()
    print("2D FOF analysis complete in {:3.2f} minutes \n".format((t1 - t0) / 60))
    
    # match detectids at large linking length if a wave group exists
    joinfriend = join(detfriend_all, w_keep, keys='detectid', join_type='left')

    if False:# will hopefully remove all this soon
        grp_by_id = joinfriend.group_by('id')
        sum_grp = grp_by_id.groups.aggregate(np.sum)
        spatial_id = grp_by_id.groups.keys
        spatial_id_to_keep1 = spatial_id[np.isfinite(sum_grp['wave_group_id'])]
        
        # also match if a group member is brighter than gmag=22 (testing this)
        sel_bright = detect_table['gmag'] < 22
        gmag_join = join( joinfriend, detect_table['detectid','gmag'][sel_bright])
        spatial_id_to_keep2 = gmag_join['id']
        
        spatial_id_to_keep = vstack([spatial_id_to_keep1, spatial_id_to_keep2])
        detfriend_1 = join(unique(spatial_id_to_keep), joinfriend)
        
        # link the rest of the detectids with smaller linking length
        
        keep_row = np.ones(np.size(detect_table), dtype=bool)
        
        for i, det in enumerate(detect_table['detectid']):
            if det in detfriend_1['detectid']:
                keep_row[i] = 0

        print("Performing FOF in 2D space with dlink=2.0 arcsec")

        kdtree, r = fof.mktree(
            detect_table["ra"][keep_row],
            detect_table["dec"][keep_row],
            np.zeros_like(detect_table["ra"][keep_row]),
            dsky=1.0,
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
            detect_table['flux_g'][keep_row],
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

    global expand_table
    
    expand_table = join(joinfriend, detect_table, keys="detectid")
    expand_table['wave_group_id'] = expand_table['wave_group_id'].filled(0)
    expand_table.rename_column('id','source_id')
    
    try:
        expand_table.remove_column('detectname')
    except:
        pass

    # combine common wavegroups to the same source_id
    # update source properties

    print('Combining nearby wavegroups and detections')
    
    t0 = time.time()
    sel = expand_table['wave_group_id'] > 0
    wid_list = np.unique(expand_table['wave_group_id'][sel])
    p = Pool(32)
    res = p.map(merge_wave_groups, wid_list)
    p.close()

    ind_list = []
    for r in res:
        if r is None:
            continue
        else:
            r.sort()
            if r in ind_list:# skip if indices are a duplicate
                continue
            else:
                ind_list.append(r)

    res_tab = fof.process_group_list(
            ind_list,
            expand_table["detectid"],
            expand_table["ra"],
            expand_table["dec"],
            0.0* expand_table["wave"],
            expand_table['flux_g'])

    for i, ind in enumerate(ind_list):
        sid_main = np.min(expand_table['source_id'][ind])
        expand_table['source_id'][ind] = sid_main

        for col in res_tab.colnames:
            if col in ['id', 'members']:
                continue
            expand_table[col][ind] = res_tab[col][i]

    t1 = time.time()                                               
    print('Done combining wavegroups in {:4.2f} min'.format( (t1-t0)/60))

    del detfriend_all, detect_table

    gaia_stars = Table.read(config.gaiacat)

    gaia_coords = SkyCoord(ra=gaia_stars["ra"] * u.deg, dec=gaia_stars["dec"] * u.deg)
    src_coords = SkyCoord(
        ra=expand_table["ra"] * u.deg, dec=expand_table["dec"] * u.deg
    )

    idx, d2d, d3d = src_coords.match_to_catalog_sky(gaia_coords)

    sel = d2d < 1.5 * u.arcsec

    gaia_match_name = np.zeros_like(expand_table["source_id"], dtype=int)
    gaia_match_name[sel] = gaia_stars["source_id"][idx][sel]

    gaia_match_dist = np.zeros_like(expand_table["source_id"], dtype=float)
    gaia_match_dist[sel] = d2d[sel].to_value(u.arcsec)
    
    expand_table["gaia_match_id"] = gaia_match_name
    expand_table["gaia_match_dist"] = gaia_match_dist
        
    expand_table.rename_column("size", "n_members")
    expand_table.rename_column("icx", "ra_mean")
    expand_table.rename_column("icy", "dec_mean")

    source_name = []

    for row in expand_table:
        source_name.append(get_source_name(row["ra_mean"], row["dec_mean"]))
    try:
        expand_table.add_column(source_name,
                                name="source_name",
                                index=1)
    except Exception:
        print('messed up source name again')

    expand_table.sort("source_id")

    # fill mask values with nans
    for col in expand_table.columns:
        try:
            expand_table[col] = expand_table[col].filled(np.nan)
            print('yes', col)
        except:
            pass
            #print('no', col)
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
        default=2.0,
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
    """ Main Function """

    parser = get_parser()
    args = parser.parse_args(argv)

    global source_table, dsky_2D, dsky_3D, dwave
    dsky_2D = args.dsky_2D
    dsky_3D = args.dsky_3D
    dwave = args.dwave

    if args.update:
        source_table = create_source_catalog(version=args.version,
                                             update=args.update)
        sys.exit('Done updating det and cont catalogs. Re-run without --update')
        
    source_table = create_source_catalog(version=args.version)

    # match band-merged WISE catalog

    wise_catalog = Table.read(op.join(
        config.imaging_dir,
        'catalogs',
        'wise-hetdexoverlap.fits'))
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
    sel_close = matched_catalog['wise_dist'] < 1.5 #arcsec
    
    print('There are {} wise matches'.format(np.size(np.unique(matched_catalog['source_id'][sel_close]))))
    # remove column info for WISE matches more than 5 arcsec away

    sel_remove = matched_catalog['wise_dist'] > 1.5 #arcsec

    matched_catalog['ra_wise'][sel_remove] = np.nan
    matched_catalog['dec_wise'][sel_remove] = np.nan
    matched_catalog['wise_dist'][sel_remove] = np.nan
    matched_catalog['primary'][sel_remove] = -1
    matched_catalog['unwise_objid'][sel_remove] = np.nan
    matched_catalog['flux_w1'][sel_remove] = np.nan
    matched_catalog['flux_w2'][sel_remove] = np.nan
    
    source_table = matched_catalog

    # add z_spec from other catlogs if it exists
    # from Steve F.
    goods_z = Table.read(op.join(
        config.imaging_dir,
        'catalogs',
        'goods_n_specz_1018_no_header.txt'),
                         names=['ra_zspec','dec_zspec','zspec',
                                'z_quality','zspec_catalog','Symbol'],
                         format='ascii.no_header')
    # remove z_quality = 1,2, low quality zspec
    sel_bad1 = (goods_z['z_quality'] == '1')
    sel_bad2 = (goods_z['z_quality'] == '2')
    sel_good_goods = np.invert(sel_bad1 | sel_bad2)
    goods_z = goods_z[sel_good_goods]
    
    #DEIMOS 10k (Hasinger et al. 2018) z_spec up to ~6
    #https://cosmos.astro.caltech.edu/news/65
    #Secure if comprehensive quality flag (Q) = 2
    deimos = Table.read(op.join(
        config.imaging_dir,
        'catalogs',
        'deimos_redshifts.tbl'), format='ascii')
    deimos.rename_column('Ra', 'ra_zspec')
    deimos.rename_column('Dec', 'dec_zspec')
    deimos['zspec_catalog'] = 'CosmosDeimos'
    sel_good_deimos = deimos['Q'] == 2
    deimos = deimos[sel_good_deimos]
                         
    #Kriek et al. (2015)
    #http://mosdef.astro.berkeley.edu/for-scientists/data-releases/
    #5: Redshift is based on single emission feature detected with 2<=S/N<3,
    #and is within 95% confidence interval of photo-z or within delta(z)=0.05
    #of pre-MOSFIRE spec-z (if it exists)
    mosdef = Table.read(op.join(
        config.imaging_dir,
        'catalogs',
        'mosdef_zcat.final_slitap.fits'))
    mosdef.rename_column('RA','ra_zspec')
    mosdef.rename_column('DEC','dec_zspec')
    mosdef.rename_column('Z_MOSFIRE', 'zspec')
    mosdef['zspec_catalog'] = 'MOSFIRE'
    sel_good_mosdef = mosdef['Z_MOSFIRE_ZQUAL'] >= 5
    mosdef = mosdef[sel_good_mosdef]
    
    #VUDS (Tasca et al. 2017), z_spec up to ~6
    #http://cesam.lam.fr/vuds/DR1/
    # zCOSMOS Spectroscopic Redshift Survey 2009ApJS..184..218L
    #https://www.eso.org/qi/catalog/show/65

    #Flag 4’s, 3’s, 2.5, 2.4, 1.5, 9.5, 9.4, 9.3 are considered secure.
    zcosbright = Table.read(op.join(
        config.imaging_dir,
        'catalogs',
        'cesam_zcosbrightspec20k_dr3_catalog_1616073679.txt'), format='ascii')
    zcosbright.rename_column('zpec','zspec')
    zcosbright.rename_column('ra','ra_zspec')
    zcosbright.rename_column('dec','dec_zspec')
    zcosbright['zspec_catalog'] = 'zCosmosBright'
    sel_zcos_good1 = zcosbright['cc']>= 9.3
    sel_zcos_good2 = (zcosbright['cc'] >= 3.0 ) & (zcosbright['cc'] < 5)
    sel_zcos_good3 = zcosbright['cc'] == 1.5
    sel_zcos_good4 = (zcosbright['cc'] >= 2.4) & (zcosbright['cc'] <= 2.5)
    sel_zcos_good = sel_zcos_good1 | sel_zcos_good2 | sel_zcos_good3 | sel_zcos_good4
    zcosbright = zcosbright[sel_zcos_good]

    deep_specz = Table.read(op.join(
        config.imaging_dir,
        'catalogs',
        'DEEP_zcosmos_spectroscopy_one_v2.6_data+header.cat'),
        format='ascii', data_start=100)
    deep_specz.rename_column('col1','zCOSMOS-deepID')
    deep_specz.rename_column('col2','zspec')
    deep_specz.rename_column('col3','flag')
    deep_specz.rename_column('col4','zphot')
    deep_specz.rename_column('col5','ra_zspec')
    deep_specz.rename_column('col6','dec_zspec')
    deep_specz['zspec_catalog'] = 'DEEP_zcosmos'
    sel_deepcos_good1 = deep_specz['flag']>= 9.3
    sel_deepcos_good2 = (deep_specz['flag'] >= 3.0 ) & (deep_specz['flag'] < 5)
    sel_deepcos_good3 = deep_specz['flag'] == 1.5
    sel_deepcos_good4 = (deep_specz['flag'] >= 2.4) & (deep_specz['flag'] <= 2.5)
    sel_deepcos_good = sel_deepcos_good1 | sel_deepcos_good2 | sel_deepcos_good3 | sel_deepcos_good4
    deep_specz = deep_specz[sel_deepcos_good]

    sdssfile = op.join( config.imaging_dir, 'sdss', 'specObj-dr16.fits')

    sdss_specz = Table.read(sdssfile)
    sel_good_sdss = sdss_specz['ZWARNING'] == 0
    sdss_specz = sdss_specz[sel_good_sdss]

    sdss_specz.rename_column('PLUG_RA', 'ra_zspec')
    sdss_specz.rename_column('PLUG_DEC', 'dec_zspec')
    sdss_specz.rename_column('CLASS', 'sdss_class')
    sdss_specz.rename_column('Z', 'zspec')
    sdss_specz.rename_column('Z_ERR', 'z_sdss_err')
    sdss_specz['zspec_catalog'] = 'sdss-dr16'

    dtab = Table.read(op.join(
        config.imaging_dir,
        'catalogs',
        'desi-hetdex-v2.0.fits'))

    sel_good_hetdex = (dtab['wave'] >= 3640 )
    sel_good_desi = (dtab['COADD_FIBERSTATUS'] == 0)
    sel_sample = sel_good_desi*sel_good_hetdex
    sel_conf = dtab['VI_quality'] >= 3
    dtab = dtab[sel_conf&sel_sample]                    
    dtab.rename_column('TARGET_RA','ra_zspec')
    dtab.rename_column('TARGET_DEC','dec_zspec')
    dtab['zspec_catalog'] = 'DESI'
    dtab.rename_column('VI_z', 'zspec')

    specz_catalogs = vstack([goods_z['zspec','ra_zspec','dec_zspec','zspec_catalog'],
                             deimos['zspec','ra_zspec','dec_zspec','zspec_catalog'],
                             mosdef['zspec','ra_zspec','dec_zspec','zspec_catalog'],
                             zcosbright['zspec','ra_zspec','dec_zspec','zspec_catalog'],
                             deep_specz['zspec','ra_zspec','dec_zspec','zspec_catalog'],
                             sdss_specz['zspec','ra_zspec','dec_zspec','zspec_catalog'],
                             dtab['zspec','ra_zspec','dec_zspec','zspec_catalog'],
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

    sel_close = matched_catalog['zspec_dist'] < 1.5 #u.arcsec

    print('There are {} zspec matches within 5 arcsec'.format(np.size(np.unique(matched_catalog['source_id'][sel_close]))))

    sel_remove = matched_catalog['zspec_dist'] > 1.5 #u.arcsec

    matched_catalog['zspec'][sel_remove] = np.nan
    matched_catalog['ra_zspec'][sel_remove] = np.nan
    matched_catalog['dec_zspec'][sel_remove] = np.nan
    matched_catalog['zspec_dist'][sel_remove] = np.nan
    matched_catalog['zspec_catalog'][sel_remove] = ''

    if False:#doing above because we need to cross-match by coordinates for HDR3
        #add desi confirmed redshifts

        dtab = Table.read(op.join(
            config.imaging_dir,
            'catalogs',
            'desi-hetdex-v1.0.fits'))

        sel_good_hetdex = (dtab['wave'] >= 3640 )
        sel_good_desi = (dtab['COADD_FIBERSTATUS'] == 0)
        sel_sample = sel_good_desi*sel_good_hetdex
        sel_conf = dtab['VI_quality'] >= 3

        hetdex_coords = SkyCoord(ra=dtab['RA_HETDEX'], dec=dtab['DEC_HETDEX'], unit='deg')
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
    
    # sort table closest to group mean position
    # so unique will produce the closest match

    src_coord = SkyCoord(ra=source_table['ra_mean'], dec=source_table['dec_mean'], unit='deg')
    det_coord = SkyCoord(ra=source_table['ra'], dec=source_table['dec'], unit='deg')

    source_table['src_separation'] = det_coord.separation(src_coord)
    source_table.sort(['src_separation'])

    print('Filling masked values with NaNs')
    
    for col in source_table.columns:
        try:
            source_table[col] = source_table[col].filled(np.nan)
            print('yes', col)
        except:
            pass
            #print('no', col)
    #remove nonsense metadata
    source_table.meta = {}

    # Apply dust corretion to flux, continuum, flux_1sigma

    source_table['flux_obs'] = source_table['flux'].copy()
    source_table['flux_err_obs'] = source_table['flux_err'].copy()
    source_table['continuum_obs'] = source_table['continuum'].copy()
    source_table['continuum_err_obs'] = source_table['continuum_err'].copy()
    source_table['flux_noise_1sigma_obs'] = source_table['flux_noise_1sigma'].copy()

    Rv = 3.1
    ext = []

    for index in np.arange( np.size(source_table['detectid'])):
        src_wave = np.array([np.double(source_table['wave'][index])])
        ext_i = extinction.fitzpatrick99(src_wave, source_table['Av'][index], Rv)[0]
        ext.append(ext_i)
       
    deredden = 10**(0.4*np.array(ext))
    deredden[np.isnan(deredden)] = 0
    source_table['flux'] = deredden * source_table['flux_obs']
    source_table['flux_err'] = deredden * source_table['flux_err_obs']
    source_table['continuum'] = deredden * source_table['continuum_obs']
    source_table['continuum_err'] = deredden * source_table['continuum_err_obs']
    source_table['flux_noise_1sigma'] = deredden * source_table['flux_noise_1sigma_obs']
    
    source_table.write("source_catalog_{}.fits".format(args.version),
                       overwrite=True)
    source_table.write("source_catalog_{}.tab".format(args.version),
                       overwrite=True, format='ascii')


if __name__ == "__main__":
    main()
