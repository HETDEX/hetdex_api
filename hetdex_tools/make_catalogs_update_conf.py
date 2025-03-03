# script to update p_conf value

import numpy as np
import astropy.table
from astropy.table import Table

import joblib

version = '5.0.0'

source_table = Table.read("source_catalog_{}.yy.fits".format(version))

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
