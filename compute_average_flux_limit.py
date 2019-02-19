"""

Compute biweight flux limits at a wavelength and
optionally plot the histogram of flux limits
at a chosen wavelength. Used to compute 
e.g. fluxlimit_4540 type columns

AUTHOR: Daniel Farrow (MPE)

"""


from __future__ import print_function

import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.io.fits import getdata
from astropy.stats import biweight_location, biweight_midvariance
from astropy.table import Table

# Parse the arguments
parser = argparse.ArgumentParser(description="Compute biweight flux limits and plot histograms")
parser.add_argument("--wl", default=4540, type=int, help="Wavelength to compute median for")
parser.add_argument("--hist", action="store_true", help="Plot histograms of flux limit")
parser.add_argument("--hist_all", action="store_true", help="Plot histograms flux limit for all inputs")
parser.add_argument("--fout", default=None, type=str, help="Ascii file to save results to")
parser.add_argument("--fn-shot-average", default=None, type=str, help="Ascii file to append shot average flim to")
parser.add_argument("files", help="Flux limit cubes to consider", nargs='+')
opts = parser.parse_args()

bins = np.linspace(0.0, 5e-16, 40)
bcens = 0.5*(bins[1:] + bins[:-1])

print("Using wavelength {:f} AA".format(opts.wl))

# Loop over the files producing median and percentiles
biwt_ls = []
biwt_vars = []
dateshot = []
ifu = []
flims_all = []

for fn in opts.files:
     
    try:
        data, header = getdata(fn, header=True)
    except IOError as e:
        print("{:s} failed!".format(fn))    
        continue

    try:
        windex = (opts.wl - header["CRVAL3"])/header["CDELT3"] - 1
    except KeyError as e:
        print("{:s} failed!".format(fn))
        continue

    # Grab a wavelength slice and collapse it
    wavelength_slice = data[windex,:,:] 
    slice_flattened = wavelength_slice.flatten()
    ztrimmed_slice = slice_flattened[slice_flattened > 0.0]

    # XXX Rescale as Karl got wrong SNR in old code version
    flims = (6.0/8.0)*1.0e-17/ztrimmed_slice

    flims_all.extend(flims)    
    biwt_ls.append(biweight_location(flims))
    biwt_vars.append(biweight_midvariance(flims))

    if opts.hist:
        n, b, p = plt.hist(flims, bins=bins, histtype='step', label=fn)
        plt.axvline(biweight_location(flims), linestyle="--", color=p[0].get_edgecolor())
  
    # Find the dateshot and the IFU 
    add_ifu_datashot = True
    if add_ifu_datashot:
        m = re.search('([0-9]*v[0-9]{3})_multi_([0-9]{3})', fn)
        ifu.append(m.group(2))
        dateshot.append(m.group(1))


n, b, p = plt.hist(flims_all, bins=bins, histtype='step')
peak = bcens[np.argmax(n)]

if opts.fn_shot_average:
    with open(opts.fn_shot_average, 'a') as fn:
        fn.write("{:s} {:e} {:e}\n".format(dateshot[0], biweight_location(flims_all), peak))

if opts.hist_all:
    plt.axvline(biweight_location(flims_all), linestyle="--", color=p[0].get_edgecolor()) 
    plt.axvline(peak, linestyle=":", color=p[0].get_edgecolor())

table = Table([dateshot, ifu, biwt_ls, np.sqrt(biwt_vars)], names=["dateshot", "ifu", "biwt_loc", "sqrt_biwt_mdvar"]) 


if opts.fout:
    table.write(opts.fout, format='ascii')
else:
    table.pprint(max_lines=-1)

if opts.hist or opts.hist_all:
    plt.legend(loc="best")
    plt.xlabel("Flux (cgs units)")
    plt.ylabel("Number")
    plt.show()
    
