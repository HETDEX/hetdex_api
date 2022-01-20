#!/usr/bin/env python
# coding: utf-8

import sys
import numpy as np
import os.path as op
import time

from astropy.table import Table, unique, join, Column
from astropy.cosmology import Planck18 as cosmo
import astropy.units as u

from hetdex_api.config import HDRconfig
import hetdex_tools.fof_kdtree as fof

from multiprocessing import Pool

# Enter the catalog version

version = '2.1.3'

config = HDRconfig()
#catfile = op.join(config.detect_dir, 'catalogs', 'source_catalog_' + version + '.fits')

catfile = 'source_catalog_2.1.4_orig.fits'
source_table = Table.read(catfile)

#fudging this for now
sel = source_table['source_id'] < 2140000000000

source_table['source_id'][sel] = source_table['source_id'][sel] + 10000000000

agn_tab = Table.read(config.agncat, format='ascii')
print('Source catalog was found at {}'.format(catfile))

wavelya = 1215.67
waveoii = 3727.8
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
                "stellartype"]:
        source_table.remove_column(col)
except:
    pass

if True:

    diagnose_tab = Table.read( '/work/05350/ecooper/stampede2/redshift-tests/diagnose_2.1.3.fits')
    diagnose_tab.rename_column('z_best', 'z_diagnose')
    diagnose_tab.rename_column('classification', 'cls_diagnose')
    
    combined = join( source_table, diagnose_tab['detectid', 'z_diagnose', 'cls_diagnose', 'stellartype'], join_type='left', keys=['detectid'])
    source_table = combined.copy()

uniq_table = unique(source_table, keys='source_id')

def add_z_guess(source_id):
    
    global source_table, agn_tab, config

    if agn_tab is None:
        agn_tab = Table.read(config.agncat, format="ascii")

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
            if det in agn_tab['detectid']:
                agn_flag = 1
                agn_det = det
            
        if agn_flag == 1:
        # get proper z's from Chenxu's catalog
            sel_det = agn_tab["detectid"] == agn_det
            z_guess = agn_tab["z"][sel_det][0]
            agn_flag = agn_tab['zflag'][sel_det][0]
            z_conf = agn_tab['zflag'][sel_det][0]
            s_type = "agn"
            z_src = 'Liu+2022'
            
        elif np.any( (group["cls_diagnose"] == "STAR" ) * (group['gmag'] < 22) ):
            s_type = "star"
            z_guess = 0.0
            z_conf = 0.9
            z_src = 'diagnose'
            
        elif np.any( (group["cls_diagnose"] == "GALAXY" ) * (group['plae_classification'] < 0.5) ):
            sel_best = (group["cls_diagnose"] == "GALAXY" ) * (group['plae_classification'] < 0.5)
            closest = np.argmin(group['src_separation'][sel_best])
            z_guess = group['z_diagnose'][sel_best][closest]
            z_conf = 0.9
            z_src = 'diagnose'            
            if np.any( np.abs( group['wave'] - waveoii*(1+z_guess)) < 10):
                s_type = "oii"
            else:
                s_type = 'lzg'
                
        elif np.any( (group["cls_diagnose"] == "GALAXY" ) * np.isnan(group['plae_classification'])):
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
                if np.isfinite(group['best_z']):
                    z_guess = float( group['best_z'])
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
                        if float(group['plae_classification']) >= 0.5:
                            z_guess = float(group['wave'])/wavelya - 1
                            s_type = 'lae'
                        else:
                            z_guess = float(group['wave'])/waveoii - 1
                            s_type = 'oii'
                else:
                    #use plae_classification
                    if group['plae_classification'] >= 0.5:
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
                    z_guess = group['best_z'][sel_wave_group][sel_most_conf]
                    z_conf = group['best_pz'][sel_wave_group][sel_most_conf]
                    # assign s_type
                    if 0.0 < z_guess < 0.5:
                        s_type='oii'
                    elif 1.8 < z_guess < 3.6:
                        s_type='lae'
                    z_src = 'elixer'
                    
                elif np.any( group['plae_classification'] < 0.4):
                    
                    sel_best = (group['plae_classification'] < 0.4) * np.isfinite( group['best_pz'])
                    sel_most_conf = np.argmin( group['src_separation'][sel_best] )
                    z_guess = group['best_z'][sel_best][sel_most_conf]
                    z_conf = group['best_pz'][sel_best][sel_most_conf]
                    s_type='oii'
                    z_src = 'elixer'
                    
                elif np.any( group['plae_classification'] >= 0.5):
                    
                    # check to see if its a Lya/CIV match, assume blue line is Lya
                    try:
                        sel_sn = group['sn'] > 5.5
                        wavelya_obs = np.min( group['wave'][sel_sn] )
                    except:
                        wavelya_obs = np.min( group['wave'] )
                    zlya_obs = wavelya_obs/wavelya - 1

                    sel_line = group['det_type'] == 'line'
                    if np.any( np.abs( group['wave'][sel_line] - waveciv * (1 + zlya_obs )) < 20):
                        print('lya, civ match for', list(group['detectid'][sel_line]))
                        z_guess = zlya_obs
                        z_conf = 0.95
                    elif np.any( np.abs( group['wave'][sel_line] - waveheii * (1 + zlya_obs )) < 20):
                        print('lya, heii match for', list(group['detectid'][sel_line]))
                        z_guess = zlya_obs
                        z_conf = 0.95
                    elif np.std(group['best_z']) < 0.02:
                        print(np.array(group['detectid']))
                        sel_most_conf = np.argmin( group['src_separation'] )
                        z_guess = group['best_z'][sel_most_conf]
                        z_conf = group['best_pz'][sel_most_conf]
                    else:
                        z_guess = -2
                        z_conf = -2
                    s_type = 'lae'
                    z_src  = 'elixer'
                else:
                    z_guess = -2
                    z_src = 'elixer'
                    s_type = 'lae'
                
    else:#except:
        print('could not classify {}'.format(source_id))

    if np.isnan(z_guess):
        z_guess = -3
        z_conf = 0
        s_type = 'none'
        z_src = ''
        
    return z_guess, z_conf, s_type, agn_flag, z_src        

src_list = uniq_table['source_id']

ntask = 48
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

all_table = join(source_table, z_table, join_type="left")
source_table = all_table.copy()

# uncluster LAEs whose line pairings are not supported by best_z

clustered_lae_index = np.where( source_table['z_hetdex'] == -2)[0]

print('Updating z_hetdex to best_z for {} detectids'.format(len(clustered_lae_index)))

sid_index = 2140300000000

for c_ind in clustered_lae_index:
    source_table['source_id'][c_ind] = sid_index
    source_table['z_hetdex'][c_ind] = source_table['best_z'][c_ind]
    source_table['z_hetdex_conf'][c_ind] = source_table['best_pz'][c_ind]
    source_table['z_hetdex_src'][c_ind] = 'elixer'
    sid_index += 1
    
source_table['z_guess'] = source_table['z_hetdex']

# assign Line ID
spectral_lines = Table.read('spectral_lines.txt', format='ascii')

source_table['line_id'] = np.chararray(len(source_table), 20, unicode=True)

for row in spectral_lines:
    sel = np.abs( source_table['wave'] - (1 + source_table['z_hetdex']) * row['lambda_vac']) < 15
    source_table['line_id'][sel] = row['species']

source_table['dwave'] = np.zeros_like(source_table['wave'])

sel_star = source_table['source_type'] == 'star'
sel_oii = source_table['source_type'] == 'oii'
sel_lae = source_table['source_type'] == 'lae'
sel_agn = source_table['source_type'] == 'agn'
sel_lzg = source_table['source_type'] == 'lzg'

print('There are {} stars, {} OII emitters, {} LAEs and {} Low-z gals'.format(
    np.sum(sel_star), np.sum(sel_oii), np.sum(sel_lae), np.sum(sel_lzg)))

#update clusering 


# sort table by source wavelength closest to z_hetdex and position         
src_wave = np.zeros_like(source_table['z_hetdex'])
src_wave[sel_oii] = (1 + source_table['z_hetdex'][sel_oii]) * waveoii
src_wave[sel_lae] = (1 + source_table['z_hetdex'][sel_lae]) * wavelya

sel_z = ( source_table['z_hetdex'] >= 1.9) * (source_table['z_hetdex'] <= 3.6)
src_wave[sel_agn*sel_z] = (1 + source_table['z_hetdex'][sel_agn*sel_z]) * wavelya

sel = sel_oii | sel_lae | sel_agn
source_table['dwave'][sel] = np.abs(src_wave[sel] - source_table['wave'][sel])

# Cluster in redshift space to group those missed by wave clustering

sel = source_table['z_hetdex'] >= -1
uniq_table = unique(source_table[sel], keys=['source_id'])

kdtree, r = fof.mktree(
                uniq_table["ra_mean"],
                uniq_table["dec_mean"],
                uniq_table["z_hetdex"],
                dsky=6, dwave=0.01)
t0 = time.time()
print("starting fof ...")
zfriend_lst = fof.frinds_of_friends(kdtree, r, Nmin=2)
t1 = time.time()
print('Done fof in {}s'.format(t1-t0))

zfriend_table = fof.process_group_list(
        zfriend_lst,
        uniq_table["source_id"],
        uniq_table["ra_mean"],
        uniq_table["dec_mean"],
        uniq_table["z_hetdex"],
        uniq_table['flux_g'],
    )
print("Generating combined table \n")

memberlist = []
friendlist = []
print(zfriend_table.colnames)

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

# Update group info for new groups

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

# fix for DAR apcor issue
source_table['flux_2.1.3'] = source_table['flux'].copy()
source_table['flux_err_2.1.3'] = source_table['flux_err'].copy()
source_table['continuum_2.1.3'] = source_table['continuum'].copy()
source_table['continuum_err_2.1.3'] = source_table['continuum_err'].copy()

source_table['flux'] = source_table['flux_2.1.3'] * source_table['apcor'] / source_table['apcor_api']
source_table['flux_err'] = source_table['flux_err_2.1.3'] * source_table['apcor'] / source_table['apcor_api']
source_table['continuum'] = source_table['continuum_2.1.3'] * source_table['apcor'] / source_table['apcor_api']
source_table['continuum_err'] = source_table['continuum_err_2.1.3'] * source_table['apcor'] / source_table['apcor_api']

source_table.rename_column('apcor', 'apcor_2.1.3')
source_table.rename_column('apcor_api', 'apcor')

source_table.sort('source_id', 'gmag')

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
im_info.rename_column('im_flux', 'flux_aper')
im_info.rename_column('im_flux_err','flux_aper_err')
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
print(len(source_table))

sel_to_switch = (source_table['source_type'] == 'lzg') & (source_table['flux_aper'] > 0)
sids_to_switch = np.unique(source_table['source_id'][sel_to_switch])

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

sel_lae_line = source_table['line_id'] == 'Ly$\\\\alpha$'
sel_z = source_table['z_hetdex'] > 1.87
lae_tab = source_table[sel_lae_line & sel_z]

lae_tab.sort('sn')

uniq_obs_lae = unique(lae_tab, keys=['source_id', 'shotid'] )

det_lae_selected = uniq_obs_lae['detectid']

for det in det_lae_selected:
    sel_det = source_table['detectid'] == det
    source_table['selected_det'][sel_det] = True


# now assign brightest detectid to star, agn, lzg groups

sel_rest = (source_table['source_type'] != 'lae') * (source_table['source_type'] != 'oii')

rest_tab = source_table[sel_rest]
rest_tab.sort('gmag')
uniq_rest = unique( rest_tab, keys=['source_id','shotid'])

for det in uniq_rest['detectid']:
    sel_det = source_table['detectid'] == det
    source_table['selected_det'][sel_det] = True

source_table['selected_det'] = source_table['selected_det'].filled(False)    

# add luminosities

source_table.add_column(Column(name='lum_lya', length=len(source_table)))
source_table.add_column(Column(name='lum_lya_err', length=len(source_table)))
source_table.add_column(Column(name='lum_oii', length=len(source_table)))
source_table.add_column(Column(name='lum_oii_err', length=len(source_table)))

# oii lum comes from flux_aper unless flux_aper is negative

sel_oii = (source_table['selected_det'] == True) & (source_table['source_type'] == 'oii') & (source_table['flux_aper'] > 0)

lum_dist = cosmo.luminosity_distance(source_table['z_hetdex'][sel_oii]).to(u.cm)
fac = 10**(-17) * 4.* np.pi * lum_dist**2

source_table['lum_oii'][sel_oii] = source_table['flux_aper'][sel_oii]*fac
source_table['lum_oii_err'][sel_oii] = source_table['flux_aper_err'][sel_oii]*fac


sel_oii = (source_table['selected_det'] == True) & (source_table['source_type'] == 'oii') & (source_table['flux_aper'] <= 0)

lum_dist = cosmo.luminosity_distance(source_table['z_hetdex'][sel_oii]).to(u.cm)
fac = 10**(-17) * 4.* np.pi * lum_dist**2

source_table['lum_oii'][sel_oii] = source_table['flux'][sel_oii]*fac
source_table['lum_oii_err'][sel_oii] = source_table['flux_err'][sel_oii]*fac

sel_lae = (source_table['selected_det'] == True) & (source_table['source_type'] == 'lae')

lum_dist = cosmo.luminosity_distance(source_table['z_hetdex'][sel_lae]).to(u.cm)
fac = 10**(-17) * 4.* np.pi * lum_dist**2

source_table['lum_lya'][sel_lae] = source_table['flux'][sel_lae]*fac
source_table['lum_lya_err'][sel_lae] = source_table['flux_err'][sel_lae]*fac

# update flux_noise_1sigma masked values:

sel_bad = source_table['flux_noise_1sigma'] > 999

source_table['flux_noise_1sigma'][sel_bad] = 999
source_table['flux_noise_1sigma_obs'][sel_bad] = 999

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
source_table.write('source_catalog_2.1.4.fits', overwrite=True)
source_table.write('source_catalog_2.1.4.tab', format='ascii', overwrite=True)

#uniq_table = unique(source_table, keys='source_id')
#uniq_table.write('source_unique_2.1.4.tab', format='ascii', overwrite=True)
#uniq_obs = unique(source_table, keys=['source_id', 'shotid'])
#uniq_obs.write('source_unique_obs_2.1.4.tab', format='ascii', overwrite=True)
