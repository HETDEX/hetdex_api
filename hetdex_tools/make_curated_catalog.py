import sys
from astropy.table import Table, join, Column
import numpy as np

from hetdex_api.detections import Detections
from hetdex_api.config import HDRconfig

version = sys.argv[0]

config = HDRconfig()

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
    sel_chi2 = (detects.chi2 <= 1.2)
    sel_cont = detects.continuum > -3
    sel_sn = detects.sn >= 4.8
    sel_chi2fib = (detects.chi2fib < 4.5)
    sel_tp = detects.throughput >= 0.08

    sel = sel_field * sel_chi2 * sel_cont * sel_sn * sel_chi2fib * sel_tp

    sel_wave = ( detects.wave >= 3550 ) * (detects.wave <= 5470)
    sel_lw = (detects.linewidth <= 14) *(detects.linewidth > 6)

    sel1 = sel * sel_wave * sel_lw

    sel_wave = ( detects.wave >= 3510 ) * (detects.wave <= 5490)
    sel_lw = (detects.linewidth <= 6)

    sel2 = sel * sel_wave * sel_lw

    sel_cat = sel1 | sel2

else:
    print("Provide a version : eg. '2.1.2'")

det_table = detects.return_astropy_table()

det_table[sel_cat].write('detect_hdr{}.tab'.format(version), format='ascii', overwrite=True)
det_table[sel_cat].write('detect_hdr{}.fits'.format(version, overwrite=True)
