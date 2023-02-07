#!/usr/bin/env python
# coding: utf-8

import sys
import numpy as np
import os.path as op
import time
import pandas as pd
import tables as tb

from astropy.table import Table, unique, join, Column, vstack
from astropy.cosmology import Planck18 as cosmo
import astropy.units as u

from hetdex_api.config import HDRconfig
import hetdex_tools.fof_kdtree as fof
#import hetdex_api.wave as wv

from multiprocessing import Pool


def zcluster_forshotid(shotid):

    global source_table
    
    sel = source_table['shotid'] == shotid
    
    uniq_table = unique(source_table[sel], keys=['source_id'])
    
    kdtree, r = fof.mktree(
        np.array( uniq_table["ra_mean"]),
        np.array( uniq_table["dec_mean"]),
        np.array(uniq_table["z_hetdex"]),
        dsky=6, dwave=0.01)
    
    zfriend_lst = fof.frinds_of_friends(kdtree, r, Nmin=2)
    
    if len(zfriend_lst)  == 0:
        return None
        
    zfriend_table = fof.process_group_list(
        zfriend_lst,
        np.array( uniq_table["source_id"]),
        np.array( uniq_table["ra_mean"]),
        np.array( uniq_table["dec_mean"]),
        np.array( uniq_table["z_hetdex"]),
        np.array( uniq_table['flux_g']),
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
    
    for zid in np.unique( zdetfriend_all['id']):
        sel_zid = zdetfriend_all['id'] == zid
        sids = zdetfriend_all['source_id'][sel_zid]
        
        sid_main = np.min(sids)
        
        for sid_i in sids:
            sel_sid = source_table['source_id'] == sid_i
            source_table['source_id'][sel_sid] = sid_main
            
        sid_ind = np.where( source_table['source_id'] == sid_main)
        
        res_tab = fof.process_group_list(
            sid_ind,
            source_table["detectid"],
            source_table["ra"],
            source_table["dec"],
            source_table["z_hetdex"],
            source_table['flux_g'])
        
        res_tab.rename_column('size', 'n_members')
        res_tab.rename_column("icx", "ra_mean")
        res_tab.rename_column("icy", "dec_mean")
        
        for col in res_tab.colnames:
            if col in ['id', 'members']:
                continue
                source_table[col][sid_ind] = res_tab[col]
                
            if res_tab['izz'] > 0.001:
                print('big redshift error on {}'.format(sid_main))
                print(res_tab['icz'],res_tab['izz'])
            else:
                source_table['z_hetdex'][sid_ind] = res_tab['icz']

        # update source type assignment
        stype_main = source_table['source_type'][source_table['source_id'] == sid_main]
        source_table['source_type'][sid_ind] = stype_main
    return

# Enter the catalog version


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

version = sys.argv[1]

print('Adding redshifts to catalog: {}'.format(version))

config = HDRconfig('hdr3')
#catfile = op.join(config.detect_dir, 'catalogs', 'source_catalog_' + version + '.fits')

catfile = 'source_catalog_{}.fits'.format(version)
source_table = Table.read(catfile)

agn_tab = Table.read(config.agncat, format='ascii')
#agn_tab = None
print('Source catalog was found at {}'.format(catfile))

wavelya = 1215.67 #in vacuum
waveoii = 3727.8 #in air
waveciv = 1549.5
waveheii = 1640.4

# match Diagnose classification table

try:
    print('Removing existing redshift and classification info')
    for col in ["z_hetdex",
                "z_hetdex_src",
                "z_hetdex_conf",
                "source_type",
                "agn_flag",
                "z_diagnose",
                "cls_diagnose",
                "stellartype",
                "selected_det"]:
        source_table.remove_column(col)
except:
    print('Could not remove existing redshift and classification info')
    pass

try:
    print('removing resolved oii values')
    for col in ['ra_aper',
                'dec_aper',
                'mag_aper',
                'mag_aper_err',
                'dist_aper',
                'catalog_name_aper',
                'filter_name_aper',
                'flux_aper',
                'flux_aper_err',
                'bkg_median_aper',
                'bkg_stddev_aper',
                'lum_lya',
                'lum_lya_err',
                'lum_oii',
                'lum_oii_err']:
        source_table.remove_column(col)
except:
    print('could not remove resolved oii values')

diagnose_tab = Table.read( '/work/05350/ecooper/stampede2/redshift-tests/hdr3.0.0_lt23/diagnose_3.0.0_lt23.fits')
diagnose_tab.rename_column('z_best', 'z_diagnose')
diagnose_tab.rename_column('classification', 'cls_diagnose')

combined = join( source_table, diagnose_tab['detectid', 'z_diagnose', 'cls_diagnose', 'stellartype'], join_type='left', keys=['detectid'])
source_table = combined.copy()

del combined, diagnose_tab

# assign redshifts to sources with single detections
sel = source_table['n_members'] == 1

# assign AGN redshifts in AGN catalog

if agn_tab is not None:
    agn_assign = source_table['source_id','z_agn'][ (~source_table['z_agn'].mask) & sel].copy()
    agn_assign.rename_column('z_agn','z_hetdex')
    agn_assign['z_hetdex_conf'] = 0.9
    agn_assign['z_hetdex_src'] = 'liu_agn'
    agn_assign['source_type'] = 'agn'
    agn_assign['agn_flag'] = 1

    sel = (source_table['n_members'] == 1) & (source_table['z_agn'].mask)

# Take Diagnose for gmag < 22 for remaining single detection sources
sel_gmag = (source_table['gmag'] < 22) & (source_table['cls_diagnose'] != 'UNKNOWN' )

diagnose_assign = source_table['source_id', 'z_diagnose', 'cls_diagnose'][sel & sel_gmag].copy()
diagnose_assign.rename_column('z_diagnose', 'z_hetdex')
diagnose_assign['source_type'] = '          '
diagnose_assign['source_type'][diagnose_assign['cls_diagnose'] == 'STAR'] = 'star'
diagnose_assign['source_type'][diagnose_assign['cls_diagnose'] == 'QSO'] = 'agn'
diagnose_assign['source_type'][diagnose_assign['cls_diagnose'] == 'GALAXY'] = 'lzg'
diagnose_assign['z_hetdex_src'] = 'diagnose'
diagnose_assign['z_hetdex_conf'] = 0.8
diagnose_assign['agn_flag'] = -1
diagnose_assign['agn_flag'][diagnose_assign['cls_diagnose'] == 'QSO'] = 1
diagnose_assign.remove_column('cls_diagnose')

# Take Elixer best_z for rest

elixer_assign = source_table['source_id','z_elixer','best_pz'][sel & np.invert(sel_gmag)]
elixer_assign.rename_column('z_elixer', 'z_hetdex')
elixer_assign.rename_column('best_pz', 'z_hetdex_conf')
elixer_assign['z_hetdex_src'] = 'elixer'
elixer_assign['agn_flag'] = -1
sel_lae = elixer_assign['z_hetdex'] > 1.88
sel_oii = elixer_assign['z_hetdex'] < 0.5
elixer_assign['source_type'] = 'gal'
elixer_assign['source_type'][sel_lae] = 'lae'
elixer_assign['source_type'][sel_oii] = 'oii'

sel = source_table['n_members'] == 1

uniq_table = unique(source_table[np.invert(sel)], keys='source_id')

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
    z_src = ''

    #print(group['best_z', 'best_pz','detectid'])
    if True:
        # Check if any member is an AGN first
        for det in group['detectid']:
            if agn_tab is not None:
                if det in agn_tab['detectid']:
                    agn_flag = 1
                    agn_det = det
            
        if agn_flag == 1:
        # get proper z's from Chenxu's catalog
            sel_det = agn_tab["detectid"] == agn_det
            z_guess = agn_tab["z"][sel_det][0]
            #agn_flag = agn_tab['zflag'][sel_det][0]
            z_conf = 1 #agn_tab['zflag'][sel_det][0]
            s_type = "agn"
            z_src = 'liu_agn'
            
        elif np.any( (group["cls_diagnose"] == "STAR" ) & (group['gmag'] < 22) ):
            s_type = "star"
            sel_best = (group["cls_diagnose"] == 'STAR') & (group['gmag'] < 22)
            s_type = "star"
            closest = np.argmin(group['src_separation'][sel_best])
            z_guess = group['z_diagnose'][sel_best][closest]
            z_conf = 0.9
            z_src = 'diagnose'

        elif np.any( (group["cls_diagnose"] == "QSO") & (group['gmag'] < 22) ):
            sel_best = (group["cls_diagnose"] == 'QSO') & (group['gmag'] < 22)
            s_type = "agn"
            closest = np.argmin(group['src_separation'][sel_best])
            z_guess = group['z_diagnose'][sel_best][closest]
            z_conf = 0.9
            z_src = 'diagnose'
                                                                        
        elif np.any( (group["cls_diagnose"] == "GALAXY" ) & (group['gmag'] < 22) ):
            sel_best = (group["cls_diagnose"] == "GALAXY" ) & (group['gmag'] < 22)
            closest = np.argmin(group['src_separation'][sel_best])
            brightest = np.argmin(group['gmag'][sel_best])
            z_guess = group['z_diagnose'][sel_best][brightest]
            z_conf = 0.9
            z_src = 'diagnose'
            if np.any( np.abs( group['wave'] - waveoii*(1+z_guess)) < 10):
                s_type = "oii"
            else:
                s_type = 'lzg'
                
        elif np.any( (group["cls_diagnose"] == "GALAXY" ) & (group['plya_classification'] < 0.5)):
            sel_best = (group["cls_diagnose"] == "GALAXY" ) & (group['plya_classification'] < 0.5)
            closest = np.argmin(group['src_separation'][sel_best])
            z_guess = group['z_diagnose'][sel_best][closest]
            z_conf = 0.9
            z_src = 'diagnose'            
            if np.any( np.abs( group['wave'] - waveoii*(1+z_guess)) < 10):
                s_type = "oii"
            else:
                s_type = 'lzg'
                
        elif np.any( (group["cls_diagnose"] == "GALAXY" ) * np.isnan(group['plya_classification'])):
            sel_best = (group["cls_diagnose"] == "GALAXY")
            closest = np.argmin(group['src_separation'][sel_best])
            z_guess = group['z_diagnose'][closest]
            z_conf = 0.9
            z_src = 'diagnose'

            if np.any( np.abs( group['wave'] - waveoii*(1+z_guess)) < 10):
                s_type = "oii"
            else:
                s_type = 'lzg'
        else:
            if np.size(group) == 1:
                if np.isfinite(group['z_elixer']):
                    z_guess = float( group['z_elixer'])
                    z_conf = float( group['best_pz'])

                    if 0.0 < z_guess < 0.5:
                        s_type = 'oii'
                    elif 1.8 < z_guess < 3.7:
                        s_type = 'lae'
                    elif z_guess < 0:
                        if float(group['wave']) < waveoii:
                            #assume this should actually be Lya
                            z_guess = float(group['wave'])/wavelya - 1
                            s_type = 'lae'
                    elif 0.5 < z_guess < 1.9:
                        if float(group['plya_classification']) >= 0.5:
                            z_guess = float(group['wave'])/wavelya - 1
                            s_type = 'lae'
                        else:
                            z_guess = float(group['wave'])/waveoii - 1
                            s_type = 'oii'
                else:
                    #use plae_classification
                    if group['plya_classification'] >= 0.5:
                        z_guess = float(group['wave'])/wavelya - 1
                        s_type = 'lae'
                    else:
                        z_guess = float(group['wave'])/waveoii - 1
                        s_type = 'oii'
                z_src = 'elixer'
            else:
                if np.any( group['wave_group_id'] > 0):
                    sel_wave_group = group['wave_group_id'] > 0
                    sel_most_conf = np.argmin( group['src_separation'][sel_wave_group] )
                    z_guess = group['z_elixer'][sel_wave_group][sel_most_conf]
                    z_conf = group['best_pz'][sel_wave_group][sel_most_conf]
                    # assign s_type
                    if 0.0 < z_guess < 0.5:
                        s_type='oii'
                    elif 1.88 < z_guess < 3.6:
                        s_type='lae'
                    else:
                        s_type='gal'

                    z_src = 'elixer'
                    
                elif np.std(group['z_elixer']) < 0.02:
                    #print(np.array(group['detectid']))
                    sel_most_conf = np.argmin( group['src_separation'] )
                    z_guess = group['z_elixer'][sel_most_conf]
                    z_conf = group['best_pz'][sel_most_conf]
                    if 0.0 < z_guess < 0.5:
                        s_type='oii'
                    elif 1.88 < z_guess < 3.6:
                        s_type='lae'
                    else:
                        s_type='gal'
                        
                    z_src  = 'elixer'
                    
                else:                    
                    # assign redshift of highest best_pz line
                    # or use highest S/N line
                    if np.any(group['best_pz'] > 0.6):
                        sel_sn = np.argmax(group['best_pz'])
                    else:
                        sel_sn = np.argmax(group['sn'])
                        
                    z_guess = group['z_elixer'][sel_sn]
                    z_conf = group['best_pz'][sel_sn]
                    z_src = 'elixer'
                    if 0.0 < z_guess < 0.5:
                        s_type='oii'
                    elif 1.88 < z_guess < 3.6:
                        z_guess = -2
                        s_type='lae'
                    else:
                        s_type='gal'

                # check to see if its a Lya/CIV match, assume blue line is Lya
#                try:
#                    sel_sn = group['sn'] > 5.5
#                    wavelya_obs = np.min( group['wave'][sel_sn] )
#                except:
#                    wavelya_obs = np.min( group['wave'] )
#                    
#                    zlya_obs = wavelya_obs/wavelya - 1
                sel_line = group['det_type'] == 'line'
                zlya_obs = z_guess
                if np.any( np.abs( group['wave'][sel_line] - waveciv * (1 + zlya_obs )) < 20):
                    print('lya, civ match for', list(group['detectid'][sel_line]))
                    z_guess = zlya_obs
                    z_conf = 0.95
                    s_type = 'agn'
                    agn_flag = 1
                    z_src = '2em'
                elif np.any( np.abs( group['wave'][sel_line] - waveoii * (1 + z_guess)) < 20):
                    s_type = 'oii'
    else:#except:
        print('could not classify {}'.format(source_id))

    if np.isnan(z_guess):
        z_guess = -3
        z_conf = 0
        s_type = 'none'
        z_src = ''
        
    return z_guess, z_conf, s_type, agn_flag, z_src        

src_list = uniq_table['source_id']

ntask = 32
print('Adding z_hetdex using {} cores'.format(ntask))
t0 = time.time()
p = Pool(ntask)
res = p.map(add_z_guess, uniq_table['source_id'])
p.close()
t1 = time.time()

print('Adding z_hetdex complete in {:5.3} m'.format( (t1-t0) / 60))

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
        [uniq_table['source_id'], z_hetdex, z_src, z_conf, s_type, agn_flag],
        names=["source_id", "z_hetdex", "z_hetdex_src", "z_hetdex_conf", "source_type", "agn_flag"]
    )

if agn_tab is not None:
    z_stack = vstack([agn_assign, diagnose_assign, elixer_assign, z_table])
else:
    z_stack = vstack([diagnose_assign, elixer_assign, z_table])


all_table = join(source_table, z_stack, join_type="left")
source_table = all_table.copy()


# uncluster LAEs whose line pairings are not supported by best_z

clustered_lae_index = np.where( source_table['z_hetdex'] == -2)[0]

print('Updating z_hetdex to best_z for {} detectids'.format(len(clustered_lae_index)))

sid_index = 3010030000000

for c_ind in clustered_lae_index:
    source_table['source_id'][c_ind] = sid_index
    source_table['z_hetdex'][c_ind] = source_table['z_elixer'][c_ind]
    source_table['z_hetdex_conf'][c_ind] = source_table['best_pz'][c_ind]
    source_table['z_hetdex_src'][c_ind] = 'elixer'
    sid_index += 1

sel_star = source_table['source_type'] == 'star'
sel_oii = source_table['source_type'] == 'oii'
sel_lae = source_table['source_type'] == 'lae'
sel_agn = source_table['source_type'] == 'agn'
sel_lzg = source_table['source_type'] == 'lzg'

# assign Line ID
spectral_lines = Table.read('/work2/05350/ecooper/stampede2/redshift-tests/spectral_lines.txt', format='ascii')

source_table['line_id'] = np.chararray(len(source_table), 20, unicode=True)

source_table['dwave'] = np.zeros_like(source_table['wave'])

for row in spectral_lines:
    
    sel = np.abs( source_table['wave'] - (1 + source_table['z_hetdex']) * row['lambda_vac']) < 15
    source_table['line_id'][sel] = str(row['species'])

    #allow AGN to have broader fit
    sel = np.abs( source_table['wave'] - (1 + source_table['z_hetdex']) * row['lambda_vac']) < 30
    source_table['line_id'][sel&sel_agn] = str(row['species'])

print('There are {} stars, {} OII emitters, {} LAEs and {} Low-z gals'.format(
    np.sum(sel_star), np.sum(sel_oii), np.sum(sel_lae), np.sum(sel_lzg)))

# sort table by source wavelength closest to z_hetdex and position         
src_wave = np.zeros_like(source_table['z_hetdex'])
src_wave[sel_oii] = (1 + source_table['z_hetdex'][sel_oii]) * waveoii
src_wave[sel_lae] = (1 + source_table['z_hetdex'][sel_lae]) * wavelya

sel_z = ( source_table['z_hetdex'] >= 1.9) * (source_table['z_hetdex'] <= 3.6)
src_wave[sel_agn*sel_z] = (1 + source_table['z_hetdex'][sel_agn*sel_z]) * wavelya

sel = sel_oii | sel_lae | sel_agn
source_table['dwave'][sel] = np.abs(src_wave[sel] - source_table['wave'][sel])

# Cluster in redshift space to group those missed by wave clustering

print('Clustering in redshift space')
shotid_list = np.unique(source_table['shotid'])

t0 = time.time()
p = Pool(24)
res = p.map(zcluster_forshotid, shotid_list)
p.close()

for r in res:
    if r is not None:
        update_table(r)
t1 = time.time()

print('Done clustering in redshift space in {:5.2f} minutes'.format((t1-t0)/60))
      
source_table.sort('source_id', 'gmag')

add_oii = False
if add_oii:
    im_info = unique( Table.read('/work/05350/ecooper/stampede2/oii/im_info.txt', format='ascii'), keys='detectid')

    sel_z = (source_table['z_hetdex'] < 0.5) * ( (source_table['source_type'] == 'oii') | (source_table['source_type'] == 'lzg'))
    oii_table = source_table[sel_z]
    oii_table.sort('gmag')
    uniq_oii_table = unique(oii_table, keys=['source_id','shotid'])
    uniq_oii_table.sort('gmag')
    
    col_keep = ['detectid',
                'ra',
                'dec',
                'catalog_name',
                'filter_name',
                'dist_baryctr',
                'mag',
                'mag_err',
                'major',
                'minor',
                'theta',
                'im_flux',
                'im_flux_err',
                'bkg_median',
                'bkg_stddev',
                'im_apcor']

    im_info = im_info[col_keep]

    im_info.rename_column('ra','ra_aper')
    im_info.rename_column('dec','dec_aper')
    im_info.rename_column('mag','mag_aper')
    im_info.rename_column('mag_err','mag_aper_err')
    im_info.rename_column('dist_baryctr','dist_aper')
    im_info.rename_column('catalog_name','catalog_name_aper')
    im_info.rename_column('filter_name', 'filter_name_aper')
    im_info.rename_column('im_flux', 'flux_aper_obs')
    im_info.rename_column('im_flux_err','flux_aper_obs_err')
    im_info.rename_column('bkg_median', 'bkg_median_aper')
    im_info.rename_column('bkg_stddev', 'bkg_stddev_aper')

    oii_to_keep = join( uniq_oii_table['detectid', 'source_id'], im_info,
                        keys='detectid',
                        join_type='left')

    oii_to_keep.remove_column('source_id')

    oii_to_keep['selected_det'] = True
    
    det_oii_to_keep = oii_to_keep['detectid']
    
    print(len(uniq_oii_table), len(oii_to_keep) )

    source_table_oii = join(source_table, oii_to_keep, keys=['detectid'], join_type='left')
    print(len(source_table_oii))
    source_table = source_table_oii.copy()

    # get intrinsic flux_aper by dereddening
    
    import extinction
    
    print('applying extinction to aperture fluxes')
    
    Rv = 3.1
    ext = []

    sel_aper = source_table['flux_aper_obs'] > 0

    for index in np.arange( np.size(source_table['detectid'])):
        if source_table['flux_aper_obs'][index] > 0:
            oii_wave = np.array([np.double( waveoii * ( 1. + source_table['z_hetdex'][index]))])
            ext_i = extinction.fitzpatrick99(oii_wave, source_table['Av'][index], Rv)[0]
            ext.append(ext_i)
        else:
            ext.append(np.nan)
        
    deredden = 10**(0.4*np.array(ext))

    source_table['flux_aper'] = np.nan
    source_table['flux_aper_err'] = np.nan

    source_table['flux_aper'][sel_aper] = deredden[sel_aper] * source_table['flux_aper_obs'][sel_aper]
    source_table['flux_aper_err'][sel_aper] = deredden[sel_aper] * source_table['flux_aper_obs_err'][sel_aper]
                    
    print(len(source_table))

    sel_to_switch = (source_table['source_type'] == 'lzg') & (source_table['flux_aper'] > 0)
    sids_to_switch = np.unique(source_table['source_id'][sel_to_switch])

    for sid in sids_to_switch:
        sel = source_table['source_id'] == sid
        source_table['source_type'][sel] = 'oii'

#change lzgs to oii if an OII line is found 
sel_oii_line = (source_table['line_id'] == 'OII')
sel_lzg = source_table['source_type'] == 'lzg'

sids_to_switch = np.unique( source_table['source_id'][sel_oii_line & sel_lzg])

for sid in sids_to_switch:
    sel = source_table['source_id'] == sid
    source_table['source_type'][sel] = 'oii'

sel_star = source_table['source_type'] == 'star'
sel_oii = source_table['source_type'] == 'oii'
sel_lae = source_table['source_type'] == 'lae'
sel_agn = source_table['source_type'] == 'agn'
sel_lzg = source_table['source_type'] == 'lzg'

print('There are {} stars, {} OII emitters, {} LAEs, {} AGN and {} Low-z detections'.format(
        np.sum(sel_star), np.sum(sel_oii), np.sum(sel_lae), np.sum(sel_agn), np.sum(sel_lzg)))

# now assigned selected det flag to LAE sample
#for lzg and oii galaxies, this will be the brightest OII line (not doing resolved OII right now)

sel_lae_line = (source_table['line_id'] == 'Lya') & (source_table['source_type'] == 'lae')
#sel_z = source_table['z_hetdex'] > 1.87
lae_tab = source_table[sel_lae_line ]
lae_tab.sort('sn')
lae_tab.reverse()

uniq_obs_lae = unique(lae_tab, keys='source_id' )
uniq_obs_lae['selected_det'] = True

oii_tab = source_table[sel_oii_line]
oii_tab.sort('sn')
oii_tab.reverse()

uniq_obs_oii = unique(oii_tab, keys='source_id')
uniq_obs_oii['selected_det'] = True

# now assign brightest detectid to star, agn, lzg groups
sel_rest = (source_table['source_type'] != 'lae') & (source_table['source_type'] != 'oii')
rest_tab = source_table[sel_rest]
rest_tab.sort('gmag')

uniq_rest = unique(rest_tab, keys=['source_id'])
uniq_rest['selected_det'] = True

selected_table = vstack([uniq_obs_lae, uniq_obs_oii, uniq_rest])

source_table2 = join(source_table, selected_table['detectid','selected_det'], join_type='left')

source_table = source_table2.copy()

del source_table2, oii_tab, lae_tab, rest_tab
del uniq_obs_oii, uniq_obs_lae, uniq_rest

source_table['selected_det'] = source_table['selected_det'].filled(False)    

# add luminosities

source_table.add_column(Column(name='lum_lya', length=len(source_table)))
source_table.add_column(Column(name='lum_lya_err', length=len(source_table)))
source_table.add_column(Column(name='lum_oii', length=len(source_table)))
source_table.add_column(Column(name='lum_oii_err', length=len(source_table)))

# oii lum comes from flux_aper unless flux_aper is negative

#sel_oii_aper = (source_table['selected_det'] == True) & (source_table['source_type'] == 'oii') & (source_table['flux_aper'] > 0) & (source_table['major'] >= 2)

#source_table['flag_aper'] = -1

#lum_dist = cosmo.luminosity_distance(source_table['z_hetdex'][sel_oii_aper]).to(u.cm)
#fac = 10**(-17) * 4.* np.pi * lum_dist**2

#source_table['lum_oii'][sel_oii_aper] = source_table['flux_aper'][sel_oii_aper]*fac
#source_table['lum_oii_err'][sel_oii_aper] = source_table['flux_aper_err'][sel_oii_aper]*fac
#source_table['flag_aper'][sel_oii_aper] = 1

#sel_oii_1 = (source_table['selected_det'] == True) & (source_table['source_type'] == 'oii') & (source_table['flux_aper'] <= 0)
#sel_oii_2 = (source_table['selected_det'] == True) & (source_table['source_type'] == 'oii') & (source_table['major'] <= 2)
#sel_oii = sel_oii_1 | sel_oii_2
#source_table['flag_aper'][sel_oii] = 0

sel_oii = (source_table['selected_det'] == True) & (source_table['source_type'] == 'oii')

lum_dist = cosmo.luminosity_distance(source_table['z_hetdex'][sel_oii]).to(u.cm)
fac = 10**(-17) * 4.* np.pi * lum_dist**2

source_table['lum_oii'][sel_oii] = source_table['flux'][sel_oii]*fac
source_table['lum_oii_err'][sel_oii] = source_table['flux_err'][sel_oii]*fac

print('Number of OII assigned detections')
#print(np.sum(sel_oii), np.sum(sel_oii_aper), np.sum( (source_table['selected_det'] == True) & (source_table['source_type'] == 'oii')))
print(np.sum(sel_oii), np.sum( (source_table['selected_det'] == True) & (source_table['source_type'] == 'oii')))  
sel_lae = (source_table['selected_det'] == True) & (source_table['source_type'] == 'lae')

lum_dist = cosmo.luminosity_distance(source_table['z_hetdex'][sel_lae]).to(u.cm)
fac = 10**(-17) * 4.* np.pi * lum_dist**2

source_table['lum_lya'][sel_lae] = source_table['flux'][sel_lae]*fac
source_table['lum_lya_err'][sel_lae] = source_table['flux_err'][sel_lae]*fac

#fill mask values and force to numpy float 32 
for col in source_table.columns:
    try:
        source_table[col] = source_table[col].filled(np.nan)
        print('yes', col)
    except:
        pass
        
    if source_table[col].dtype == '>f8':
        print(col)
        source_table[col] = source_table[col].astype(np.float32)

source_table.sort('source_id')

source_table['source_type'] = source_table['source_type'].astype(str)

sel_sa22 = source_table['field'] == 'ssa22'
sel_notsa22 = np.invert(sel_sa22)

#add sn_im info:

fileh = tb.open_file(op.join( config.hdr_dir['hdr3'], 'catalogs/line_images_3.0.0.h5'),'r')
sn_im_table = unique( Table( fileh.root.Info.read()), keys='detectid')
fileh.close()

source_table2 = join(source_table, sn_im_table['detectid','sn_im'], join_type='left')

print(len(source_table), len(source_table2))

# add DEE classifications
dee = unique( Table.read('/work/05350/ecooper/stampede2/hdr3/catalogs/Source_tsne_DEE.csv'), keys='detectid') # there are some duplicate rows
dee['detectid'] = dee['detectid'].astype(int)
sel_bad = (dee['tsne_x'] > 3) & (dee['tsne_y'] < -1)
dee['flag_dee_tsne'] = np.invert( sel_bad).astype(int)
dee['dee_prob'] = dee['prob'].filled(-1)

combined = join(source_table2, dee['detectid','dee_prob', 'tsne_x','tsne_y', 'flag_dee_tsne'], join_type='left')

print(len(source_table2), len(combined))

#add Ben Thomas + Ben Ayers ML plae predictions

pred_df = pd.read_pickle('/work/05350/ecooper/wrangler/team_classify/shared/MLwork/hdr3-plya_vs_pae/all_predictions_hdr3.pkl')
pred_tab = Table.from_pandas(pred_df)
pred_tab = unique(pred_tab, keys='detectid')
pred_tab['detectid'] = pred_tab['detectid'].astype(int)

source_table2 = join(combined, pred_tab['detectid','pred_prob_lae'], join_type='left')

print(len(combined), len(source_table2))

source_table = unique( source_table2, keys='detectid')

sel_sa22 = source_table['field'] == 'ssa22'
sel_notsa22 = np.invert(sel_sa22)

# finalize flags

source_table['flag_apcor'] = np.invert( (source_table['apcor'] < 0.5) & (source_table['sn'] < 6)).astype(int)
#we can keep the AGN vetted detections with low apcor
source_table['flag_apcor'][source_table['z_agn'] > 0] = 1

source_table['flag_3540'] = ((source_table['wave'] < 3538) | (source_table['wave'] > 3545)).astype(int)
source_table['flag_3540'][source_table['z_agn'] > 0] = 1

source_table['flag_baddet'] = source_table['flag_baddet'].filled(1)
source_table['flag_badpix'] = source_table['flag_badpix'].filled(1)
source_table['flag_meteor'] = source_table['flag_meteor'].filled(1)
source_table['flag_gal'] = source_table['flag_gal'].filled(1)
source_table['flag_seldet'] = source_table['selected_det'].astype(int)
source_table['flag_fwhm'] = (source_table['fwhm'] <= 2.66).astype(int)
source_table['flag_dee_tsne'] = source_table['flag_dee_tsne'].filled(-1)

source_table['flag_best'] = source_table['flag_badpix'] * source_table['flag_apcor'] * source_table['flag_3540'] * source_table['flag_baddet'] * source_table['flag_meteor'] * source_table['flag_gal']

for col in source_table.columns:
    try:
        source_table[col] = source_table[col].filled(np.nan)
        print("yes", col)
    except:
        print('no', col)
# remove nonsense metadata
source_table.meta = {}
                                                                                        
source_table[sel_notsa22].write('source_catalog_{}.z.fits'.format(version), overwrite=True)
source_table[sel_notsa22].write('source_catalog_{}.z.tab'.format(version), format='ascii', overwrite=True)

source_table[sel_sa22].write('source_catalog_{}_sa22.fits'.format(version), overwrite=True)
