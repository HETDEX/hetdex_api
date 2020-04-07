import numpy as np

from hetdex_api.detections import Detections
from hetdex_api.config import HDRconfig
import tables as tb
from astropy.table import Table, Column, join

class Elixer(tb.IsDescription):
    detectid = tb.Int64Col(pos=0)
    mag_sdss_g = tb.Float32Col()
    mag_sdss_g_err = tb.Float32Col()
    plae_sdss_g = tb.Float32Col()
    plae_sdss_g_max = tb.Float32Col()
    plae_sdss_g_min = tb.Float32Col()
    #ELiXer combined (rules, inv variance, weights and Bayes) classification info
    combined_plae = tb.Float32Col()
    combined_plae_err = tb.Float32Col()
    plae_classification = tb.Float32Col()
    combined_continuum = tb.Float32Col()
    combined_continuum_err = tb.Float32Col()
                                        

config = HDRconfig('hdr2')

filedet = tb.open_file('/data/05350/ecooper/hdr2/detect/detect_hdr2_copy.h5', 'a')

fileelix = tb.open_file(config.elixerh5, 'r')

detectid = filedet.root.Detections.cols.detectid[:]

elix_table = Table(fileelix.root.Detections.read())

try:
    filedet.remove_node(filedet.root.Elixer, recursive=True)
except:
    pass
    
tableElixer = filedet.create_table(filedet.root, "Elixer", Elixer, "Elixer Calibration Info")

for detid in detectid:
    row = tableElixer.row

    sel_det = elix_table['detectid'] == detid

    if np.sum(sel_det) == 1:
        for col in Elixer().columns:
            row[col] = elix_table[col][sel_det]
        row.append()
    elif np.sum(sel_det) == 0:
        row.append()
        print('No elixer match for %s' % detid)
    else:
        row.append()
        print('More than 1 elixer match for %s' % detid)

tableElixer.flush()
tableElixer.cols.detectid.create_csindex()

fileelix.close()
filedet.close()
