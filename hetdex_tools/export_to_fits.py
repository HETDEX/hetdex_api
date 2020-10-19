"""
Code to export fiber spectra from the
Fibers table in a shot h5 file

Argument
--------
arg1 = h5 file path

Examples
--------
python3 export_to_fits.py

"""
import tables
import sys
from astropy.io import fits
from astropy.table import Table
import os.path as op
import numpy as np

filename = sys.argv[1]

t = tables.open_file(filename, 'r')
ra = t.root.Data.Fibers.cols.ra[:]
dec = t.root.Data.Fibers.cols.dec[:]
fwhm = t.root.Shot.cols.fwhm_virus[0]

multiframe = t.root.Data.Fibers.cols.multiframe[:]
fibnum = t.root.Data.Fibers.cols.fibnum[:]

multiname = []

for mf, fibn in zip(multiframe, fibnum):
    multiname.append(mf.decode('utf-8') + '_' + str(fibn).zfill(3))

spectra = t.root.Data.Fibers.cols.calfib[:]
error = t.root.Data.Fibers.cols.calfibe[:]
spectra[np.isnan(spectra)] = 0
error[np.isnan(error)] = 0

T = Table([ra, dec, multiname], names=['ra', 'dec', 'multiname'])
h = fits.PrimaryHDU()
h.header['fwhm'] = fwhm
L = fits.HDUList([h, fits.ImageHDU(spectra),
                  fits.ImageHDU(error), fits.BinTableHDU(T)])
fitsname = op.basename(filename).split('.h5')[0] + '.fits'
L.writeto(fitsname, overwrite=True)
t.close()
