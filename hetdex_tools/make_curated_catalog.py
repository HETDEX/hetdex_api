#python3 make_curated_catalog.py 2.1.2
#
import sys
import os.path as op
import tables as tb
from astropy.table import Table, join, Column
import numpy as np

from hetdex_api.detections import Detections
from hetdex_api.config import HDRconfig

version = str(sys.argv[1])

config = HDRconfig()

# Note because refine is constantly updated, it isn't possible to
# truly replicate older catalogs. TODO for HDR3

detects = Detections(survey='hdr2.1', catalog_type='lines').refine()

sel_field = (detects.field == 'cosmos') |(detects.field == 'dex-fall') | (detects.field == 'dex-spring') | (detects.field == 'egs') | (detects.field == 'goods-n')

if version == '2.1.1':
    sel_chi2 = detects.chi2 < 1.2
    sel_wave = ( detects.wave >= 3510 ) * (detects.wave <= 5490)
    sel_lw = (detects.linewidth <= 6)
    sel_cont = detects.continuum > -3
    sel_sn = detects.sn >= 4.8
    sel_chi2fib = (detects.chi2fib < 4.5)

    sel_cat = sel_field * sel_chi2 * sel_wave * sel_lw * sel_cont * sel_sn * sel_chi2fib

elif version == '2.1.2':
    sel_cut1 = (detects.sn>=7) * (detects.chi2<=2.5)
    sel_cut2 = (detects.sn>=4.8) *  (detects.sn<7) * (detects.chi2<=1.2)

    sel_cont = detects.continuum > -3
    sel_chi2fib = (detects.chi2fib < 4.5)
    sel_tp = detects.throughput >= 0.08

    sel = sel_field * sel_cont * sel_chi2fib * sel_tp * (sel_cut1 | sel_cut2)

    sel_wave = ( detects.wave >= 3550 ) * (detects.wave <= 5470)
    sel_lw = (detects.linewidth <= 14) *(detects.linewidth > 6)

    sel1 = sel * sel_wave * sel_lw

    sel_wave = ( detects.wave >= 3510 ) * (detects.wave <= 5490)
    sel_lw = (detects.linewidth <= 6)

    sel2 = sel * sel_wave * sel_lw

    sel_cat = sel1 | sel2

    det_table = detects[sel_cat].return_astropy_table()

    elixer_file = op.join(config.detect_dir, 'catalogs', 'elixer.2.1.2.h5')
    elixer_cat = tb.open_file(elixer_file, 'r')
    
    cls = []
    mlname = []
    mlz = []
    mlprob = []

    for row in det_table:
        detectid_obj = row['detectid']
        try:
            elix_row = elixer_cat.root.Detections.read_where('detectid == detectid_obj')
            row['plae_classification'] = elix_row['plae_classification']
            row['combined_plae'] = elix_row['combined_plae']
            row['combined_plae_err'] = elix_row['combined_plae_err']
            mlname.append(elix_row['multiline_name'][0].decode())
            cls.append(elix_row['classification_labels'][0].decode())
            mlz.append(elix_row["multiline_z"][0])
            mlprob.append(elix_row['multiline_prob'][0])
        except Exception:
            mlname.append('')
            cls.append('')
            mlz.append(False)
            mlprob.append(0.0)

    det_table.add_column(mlname, name='multiline_name')
    det_table.add_column(cls, name='classification_labels')

else:
    print("Provide a version : eg. 2.1.2")
    sys.exit()

det_table.write('detect_hdr{}.tab'.format(version),
                         format='ascii',
                         overwrite=True)
det_table.write('detect_hdr{}.fits'.format(version), overwrite=True)
