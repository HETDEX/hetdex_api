#!/usr/bin/env python
# coding: utf-8

import numpy as np
import os.path as op
import time

from astropy.table import Table, unique, join

from hetdex_tools.source_catalog import plot_source_group
from hetdex_api.elixer_widget_cls import ElixerWidget
from hetdex_api.query_widget import QueryWidget
from hetdex_api.amp_widget import AmpWidget
from hetdex_api.config import HDRconfig

import matplotlib.pyplot as plt
from multiprocessing import Pool

# Enter the catalog version

version = '2.1.3'

config = HDRconfig()
#catfile = op.join(config.detect_dir, 'catalogs', 'source_catalog_' + version + '.fits')

catfile = 'source_catalog_2.1.3.fits'
source_table = Table.read(catfile)
agn_tab = Table.read(config.agncat, format='ascii')
print('Source catalog was found at {}'.format(catfile))

wavelya = 1215.67
waveoii = 3727.8
waveciv = 1549.5
waveheii = 1640.4

# match Diagnose classification table

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
            z_src = 'chenxu'
            
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
            sel_best = (group["cls_diagnose"] == "GALAXY" )
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
                z_guess = float( group['best_z'])
                z_conf = float( group['best_pz'])
                if 0.0 < z_guess < 0.5:
                    s_type='oii'
                elif 1.8 < z_guess < 3.6:
                    s_type='lae'
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
                    
                elif np.any( group['plae_classification'] < 0.5):
                    
                    sel_best = (group['plae_classification'] < 0.5) * np.isfinite( group['best_pz'])
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
                    
                    if np.any( np.abs( group['wave'] - waveciv * (1 + zlya_obs )) < 20):
                        print('lya, civ match for', list(group['detectid']))
                        z_guess = zlya_obs
                        z_conf = 0.95
                    elif np.any( np.abs( group['wave'] - waveheii * (1 + zlya_obs )) < 20):
                        print('lya, heii match for', list(group['detectid']))
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
                    s_type='lae'
                    z_src = 'elixer'
                else:
                    z_guess = -3
                    
                
    else:#except:
        print('could not classify {}'.format(source_id))
            
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
    z_hetdex.append(np.float(r[0]))
    z_conf.append(np.float(r[1]))
    s_type.append((r[2]))
    agn_flag.append(np.float(r[3]))
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

sid_index = 2130200000000
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

# sort table by source wavelength closest to z_hetdex and position         
src_wave = np.zeros_like(source_table['z_hetdex'])
src_wave[sel_oii] = (1 + source_table['z_hetdex'][sel_oii]) * waveoii
src_wave[sel_lae] = (1 + source_table['z_hetdex'][sel_lae]) * wavelya

sel_z = ( source_table['z_hetdex'] >= 1.9) * (source_table['z_hetdex'] <= 3.6)
src_wave[sel_agn*sel_z] = (1 + source_table['z_hetdex'][sel_agn*sel_z]) * wavelya

sel = sel_oii | sel_lae | sel_agn
source_table['dwave'][sel] = np.abs(src_wave[sel] - source_table['wave'][sel])

source_table.sort(['ra','dec','src_separation'])

source_table.write('source_catalog_2.1.3.fits', overwrite=True)
source_table.write('source_catalog_2.1.3.tab', format='ascii', overwrite=True)
uniq_table = unique(source_table, keys='source_id')
uniq_table.write('source_unique_2.1.3.tab', format='ascii', overwrite=True)
uniq_obs = unique(source_table, keys=['source_id', 'shotid'])
uniq_obs.write('source_unique_obs_2.1.3.tab', format='ascii', overwrite=True)
