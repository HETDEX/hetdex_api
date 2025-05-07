# script to update p_conf value and apply other updates
import sys

import numpy as np
import astropy.table
from astropy.table import Table, join
import tables as tb
import joblib

version = '5.0.1'

update_det_flags = True

source_table = Table.read("source_catalog_{}.yy.fits".format(version))

if update_det_flags:
    print('Updating det flags')
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
                
clf = joblib.load("/work/05350/ecooper/stampede2/cosmos/calibration/rf_clf_3.0_20250216_lowsn.joblib")
    
X = np.zeros(shape=( len(source_table), len( clf.columns) ))

for i in range(len( clf.columns)):
    X[:, i]= source_table[clf.columns[i]]
    
p_conf = clf.predict_proba(X)[:,1]

source_table['p_conf'] = 1.0
sel_sample = (source_table['gmag']> 22) * (source_table['sn'] < 6.5) * (source_table['agn_flag']==-1)
source_table['p_conf'][sel_sample] = p_conf[sel_sample]

clf = joblib.load("/work/05350/ecooper/stampede2/cosmos/calibration/rf_clf_3.0_20250216_highsn.joblib")
X = np.zeros(shape=( len(source_table), len( clf.columns) ))

for i in range(len( clf.columns)):
    X[:, i]= source_table[clf.columns[i]]
    
p_conf = clf.predict_proba(X)[:,1]
sel_sample = (source_table['gmag']> 22) * (source_table['sn'] >= 6.5) * (source_table['agn_flag']==-1)
source_table['p_conf'][sel_sample] = p_conf[sel_sample]

# add in Shiro's CNN Score Values
cnn_score = Table.read(
    '/scratch/projects/hetdex/hdr5/catalogs/ml/cnn_mukae/cnn_2D_Spectra_hdr5.0.0_dex.txt',#8.0.22_k-fold updated 2025-05-07
#    '/scratch/projects/hetdex/hdr5/catalogs/ml/cnn_mukae/lae_model_hdr5.0.0_lae_dex_small.txt',
    format='ascii'
)

source_table2 = join( source_table, cnn_score, join_type='left')

source_table = source_table2.copy()
source_table['CNN_Score_2D_Spectra'] = source_table['CNN_Score_2D_Spectra'].filled(1)

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

#add column to indicate if shot is used for cosmology
survey_use = Table.read( '/scratch/projects/hetdex/hdr5/survey/survey_use_{}.txt'.format( version), format='ascii')
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
