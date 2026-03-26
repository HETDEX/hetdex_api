# script to update p_conf value and apply other updates
import sys

import numpy as np
import astropy.table
from astropy.table import Table, join, vstack
import tables as tb
import joblib

from astropy.cosmology import Planck18 as cosmo


version = '5.0.2'

BADVAL = -999.0

update_det_flags = True

source_table = Table.read("source_catalog_{}.yy.fits".format(version))

# add in apcor for continuum sources

elixh5 = tb.open_file('/scratch/projects/hetdex/hdr5/detect/elixer_hdr345_cluster_cat.h5', 'r')

elix_cont_tab = Table( elixh5.root.Detections.read())

print('length of source table before apcor fix', len(source_table))

tmp = join(source_table, elix_cont_tab['detectid', 'apcor_4500'],
           keys='detectid', join_type='left')

mask = source_table['apcor'] == 0
source_table['apcor'][mask] = tmp['apcor_4500'][mask]
print('length of source table after apcor fix', len(source_table))

# fix some bad OII aper fluxes

sel_bad_aper = ( (source_table['flux_aper']> 1e6) | ( source_table['flux_aper'] < 0) ) 
source_table['flux_aper'][sel_bad_aper] = BADVAL
source_table['flux_aper_err'][sel_bad_aper] = BADVAL

# then fix the good oiis to use pipeline flux for these


sel_oii_update = sel_bad_aper & (source_table['source_type'] == 'oii') & (source_table['flag_aper'] == 1)

print( 'Updating flux and lum values for {} OIIs'.format(np.sum(sel_oii_update)))

source_table['flux_oii'][sel_oii_update] = source_table['flux'][sel_oii_update]
source_table['flux_oii_err'][sel_oii_update] = source_table['flux_err'][sel_oii_update]
source_table['flag_aper'][sel_oii_update] = 1

lum_dist = cosmo.luminosity_distance(source_table["z_hetdex"][sel_oii_update]).to(u.cm)
fac = 10 ** (-17) * 4.0 * np.pi * lum_dist**2

source_table["lum_oii"][sel_oii_update] = source_table["flux_oii"][sel_oii_update] * fac
source_table["lum_oii_err"][sel_oii_update] = source_table["flux_oii_err"][sel_oii_update] * fac

# force other weird negative flux measurements to same BADVAL

for col in ['flux_aper', 'flux_oii', 'flux_oii_err', 'lum_oii', 'lum_oii_err', 'lum_lya', 'lum_lya_err']:
    sel = (source_table[col] < 0) & (source_table[col] > BADVAL)
    print(col, np.sum(sel))
    source_table[col][sel] = BADVAL

# add in the bad qso dets

baddets = np.loadtxt('/home1/05350/ecooper/hetdex_api/known_issues/hdr3/bad_qso_diagnose_check.dets', dtype=int)

# Create a mask for all rows where detectid is in dets
mask = np.isin(source_table['detectid'], baddets)

# Apply updates using the mask
source_table['flag_baddet'][mask] = 0
source_table['flag_best'][mask] = 0

# if a continuum source is faint, but part of a larger group then keep it

sel_faint_cont_keep = (source_table['flag_faint_cont'] == 0) * (source_table['n_members']>1)
source_table['flag_faint_cont'][sel_faint_cont_keep] = 1
source_table['flag_best'][sel_faint_cont_keep] = 1

if update_det_flags:
    print('Updating det flags')
    if version == '5.0.2':
        detflags_tab = Table.read('/scratch/projects/hetdex/hdr5/catalogs/det_flags_5.0.1.fits')
    else:
        detflags_tab = Table.read('/scratch/projects/hetdex/hdr5/catalogs/det_flags_{}.fits'.format(version))

    flag_cols = detflags_tab.colnames

    # Convert detectid to a NumPy array for fast searching
    source_ids = source_table['detectid'].data
    detflag_ids = detflags_tab['detectid'].data

    # Create a mapping from detectid to row index in source_table
    source_index_map = {id_: i for i, id_ in enumerate(source_ids)}

    # Find the common indices in both tables
    mask = np.isin(detflag_ids, source_ids)
    detflag_common = detflags_tab[mask]  # Filter only existing IDs

    # Get corresponding indices in source_table
    source_indices = np.array([source_index_map[id_] for id_ in detflag_common['detectid']])

    # Perform vectorized update only where values differ
    for col in flag_cols[1:]:  # Skip "detectid"
        source_values = source_table[col][source_indices]
        detflag_values = detflag_common[col]

        # Identify where values have changed
        changed_mask = source_values != detflag_values

        print('updating for ', col, np.sum(changed_mask))
        # Only update where changes exist
        source_table[col][source_indices[changed_mask]] = detflag_values[changed_mask]

        # for now will assume any changed value is updated to a flag so flag_best==0
        source_table['flag_best'][source_indices[changed_mask]] = 0

    print('Done updating det flags')

    #sys.exit()

for c in source_table.colnames:
    
    if isinstance(source_table[c], astropy.table.column.MaskedColumn):
        print(f"Converting masked column: {c}")
        
        if 'flag' in c:
            source_table[c] = np.nan_to_num(np.array(source_table[c]),nan=1)
        else:
            source_table[c] = np.nan_to_num(np.array(source_table[c]),nan=-99.9)
                
#clf = joblib.load("/work/05350/ecooper/stampede2/cosmos/calibration/rf_clf_3.0_20250216_lowsn.joblib")
clf = joblib.load("/work/05350/ecooper/stampede2/cosmos/calibration/rf_clf_20251006.joblib")
X = np.zeros(shape=( len(source_table), len( clf.columns) ))

for i in range(len( clf.columns)):
    X[:, i]= source_table[clf.columns[i]]
    
p_conf = clf.predict_proba(X)[:,1]

source_table['p_conf'] = 1.0
sel_sample = (source_table['gmag']> 22) * (source_table['agn_flag']==-1)
source_table['p_conf'][sel_sample] = p_conf[sel_sample]

#clf = joblib.load("/work/05350/ecooper/stampede2/cosmos/calibration/rf_clf_3.0_20250216_highsn.joblib")
#X = np.zeros(shape=( len(source_table), len( clf.columns) ))

#for i in range(len( clf.columns)):
#    X[:, i]= source_table[clf.columns[i]]
    
#p_conf = clf.predict_proba(X)[:,1]
#sel_sample = (source_table['gmag']> 22) * (source_table['sn'] >= 6.5) * (source_table['agn_flag']==-1)
#source_table['p_conf'][sel_sample] = p_conf[sel_sample]

# add in Shiro's CNN Score Values
cnn_score_lae = Table.read(
    '/scratch/projects/hetdex/hdr5/catalogs/ml/cnn_mukae/cnn_lae_20260316.txt', #updated 2026-03-26
#    '/scratch/projects/hetdex/hdr5/catalogs/ml/cnn_mukae/cnn_5.0.3_lae.txt', # updated 2025-09-26, matches Karls
#    '/scratch/projects/hetdex/hdr5/catalogs/ml/cnn_mukae/cnn_2D_Spectra_hdr5.0.0.txt', #8.0.22-k-fold updated 2025-05-17
#    '/scratch/projects/hetdex/hdr5/catalogs/ml/cnn_mukae/cnn_2D_Spectra_hdr5.0.0_dex.txt',#8.0.22_k-fold updated 2025-05-07
#    '/scratch/projects/hetdex/hdr5/catalogs/ml/cnn_mukae/lae_model_hdr5.0.0_lae_dex_small.txt',
    format='ascii'
)
#cnn_score_oii = Table.read('/scratch/projects/hetdex/hdr5/catalogs/ml/cnn_mukae/cnn_2D_Spectra_hdr5.0.1_oii.txt', format='ascii')
#cnn_score_oii = Table.read('/scratch/projects/hetdex/hdr5/catalogs/ml/cnn_mukae/cnn_5.0.3_oii.txt', format='ascii')

cnn_score_oii = Table.read('/scratch/projects/hetdex/hdr5/catalogs/ml/cnn_mukae/cnn_oii_20260316.txt', format='ascii')

cnn_score = vstack([cnn_score_lae, cnn_score_oii])

source_table2 = join( source_table, cnn_score, join_type='left')

source_table = source_table2.copy()
source_table['CNN_Score_2D_Spectra'] = source_table['CNN_Score_2D_Spectra'].filled(1)
source_table['CNN_Score_2D_Spectra'][source_table['CNN_Score_2D_Spectra'] < 0] = 1.0

# add Yuting Shen's DensFlow value
densflow = Table.read('/scratch/projects/hetdex/hdr5/catalogs/ml/densflow_yuting/predictions_all.txt',
                      format='ascii')
densflow["detectid"] = densflow["detectid"].astype(int)
densflow["densflow_score"] = densflow["score"].astype(np.float32)

source_table2 = join( source_table, densflow['detectid', 'densflow_score'], join_type='left')

source_table = source_table2.copy()
source_table['densflow_score'] = source_table['densflow_score'].filled(1)

 # remove nonsense metadata
source_table.meta = {}

for s in [20231014013,
          20231017015,
          20231104011,
          20231104015,
          20231108013,
          20231204011,
          20231206015]:
    source_table['field'][ source_table['shotid'] == s] = "dex-fall"

# call all EGS frames dex-spring

source_table['field'][ source_table['field'] == 'egs'] = 'dex-spring'

for col in source_table.colnames:
    if col == 'selected_det':
        continue
    if isinstance(source_table[col][0], str):
        source_table[col] = source_table[col].astype(str)
        print(col, ' is string')
    try:
        source_table[col] = source_table[col].filled(BADVAL)
        #print(col)
    except:
        pass
    
    if col in ['sn','flux','flux_err','flux_aper','flux_aper_err','flux_lya','flux_lya_err','lum_lya','lum_lya_err','flux_oii','flux_oii_err','lum_oii','lum_oii_err']:
        try:
            sel_zero = source_table[col] == 0.0
            source_table[col][sel_zero] = BADVAL
            print('zero fill', col, np.sum(sel_zero))
        except:
            pass

#add column to indicate if shot is used for cosmology
#if version=='5.0.2':
#    survey_use = Table.read( '/scratch/projects/hetdex/hdr5/survey/survey_use_{}.txt'.format('5.0.1'), format='ascii')
#else:
# created survey_use file on 2025-11-17 so can use it now
survey_use = Table.read( '/scratch/projects/hetdex/hdr5/survey/survey_use_{}.txt'.format( version), format='ascii')

for s in [20170621008, 20180224008, 20180805010, 20180210012]:
    sel_shot = source_table['shotid'] == s
    print('Removing {} detections from {}'. format(np.sum(sel_shot), s))
    print(len(source_table))
    source_table = source_table[np.invert( sel_shot)]
    print(len(source_table))
          
survey_use['flag_shot_cosmology'] = 1

source_table2 = join( source_table, survey_use['shotid', 'flag_shot_cosmology'], join_type='left')
source_table = source_table2.copy()
source_table['flag_shot_cosmology'] = source_table['flag_shot_cosmology'].filled(0)

# check to see if flag_best==0 and flag_seldet==1 and n_members>0

sel_other = source_table["field"] == "other"
sel_notother = np.invert(sel_other)

source_table[sel_notother].write(
    "source_catalog_{}.z.fits".format(version), overwrite=True
)
source_table[sel_notother].write(
    "source_catalog_{}.z.tab".format(version), format="ascii", overwrite=True
)
source_table[sel_other].write(
    "source_catalog_{}_other.fits".format(version), overwrite=True
)

# create H5 version

fileh = tb.open_file('source_catalog_{}.h5'.format(version), 'w')

tableMain = fileh.create_table(fileh.root, "SourceCatalog", obj=source_table.as_array())

tableMain.cols.detectid.create_csindex()
tableMain.cols.shotid.create_csindex()
tableMain.cols.source_id.create_csindex()
tableMain.cols.wave_group_id.create_csindex()
tableMain.cols.ra.create_csindex()

tableMain.flush()
fileh.close()
