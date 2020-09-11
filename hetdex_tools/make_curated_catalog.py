from astropy.table import Table, join, Column
import numpy as np

from hetdex_api.detections import Detections

detects = Detections(survey='hdr2.1', catalog_type='lines').refine()

sel_field = (detects.field == 'cosmos') |(detects.field == 'dex-fall') | (detects.field == 'dex-spring') | (detects.field == 'egs') | (detects.field == 'goods-n')
sel_chi2 = detects.chi2 < 1.2
sel_wave = ( detects.wave >= 3510 ) * (detects.wave <= 5490)
sel_lw = (detects.linewidth <= 6)
sel_cont = detects.continuum > -3
sel_sn = detects.sn >= 4.8
sel_chi2fib = (detects.chi2fib < 4.5)

sel = sel_field * sel_chi2 * sel_wave * sel_lw * sel_cont * sel_sn * sel_chi2fib

det_table = detects.return_astropy_table()

det_table[sel].write('detect_hdr2.1.1.tab', format='ascii', overwrite=True)
det_table[sel].write('detect_hdr2.1.1.fits', overwrite=True)
