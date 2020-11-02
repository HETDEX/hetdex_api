# -*- coding: utf-8 -*-
"""
Created on March 29 2019
@author: Erin Mentuch Cooper

Extracts 2D contour images for a specific detection

You can specify the size in pixels..more options to come

python3 get_countour.py -d 1000597882 -size 15 -step 0.5

"""

from __future__ import print_function

import sys
import argparse as ap
import os
import numpy as np
import tables as tb
from scipy.stats import norm

from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, join, Column

import matplotlib
matplotlib.use("agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.visualization import ZScaleInterval

from hetdex_api.config import HDRconfig
from hetdex_api.input_utils import setup_logging
from hetdex_api.detections import Detections
from hetdex_tools.get_spec import get_spectra
from hetdex_tools.line_fitting import *

from astropy.modeling import models
from specutils import Spectrum1D, SpectralRegion
from specutils.analysis import equivalent_width, line_flux
from specutils.fitting import estimate_line_parameters
from specutils.manipulation import extract_region

from elixer import catalogs

plt.style.use('fivethirtyeight')

def main(argv=None):
    """ Main Function """
    # Call initial parser from init_utils
    parser = ap.ArgumentParser(
        description='''Create a contour plot for a given det or ra/dec/wave''',
        add_help=True
    )

    parser.add_argument(
        "--gridstep",
        "--step",
        help="Step size for grid in arcsec",
        default=0.5,
        type=float
    )

    parser.add_argument(
        "--gridsize",
        "--size",
        help="Step size for grid in arcsec",
        default=15,
        type=float
    )

    parser.add_argument(
        "-s",
        "--shotid",
        help="""ShotID, e.g., 20170321v009, YYYYMMDDvOBS""",
        type=str,
        default=None,
    )

    parser.add_argument(
        "-d",
        "--detectid",
        type=int,
        help="""detectid""",
        default=None
    )

    parser.add_argument(
        "-dets",
        "--dets",
        help="""filelist of detectids""",
        default=None
    )

    parser.add_argument(
        "--ID",
        help="""Object Name if not a detection""",
        default=None)

    parser.add_argument(
        "-ra",
        "--ra",
        help="""ra, e.g., right ascension in degrees""",
        type=float,
        default=None,
    )

    parser.add_argument(
        "-dec",
        "--dec",
        help="""ra, e.g., right ascension in degrees""",
        type=float,
        default=None,
    )

    parser.add_argument(
        "-rad",
        "--rad",
        help="""radius, e.g., aperture radius in arcsec""",
        type=float,
        default=3.0,
    )

    parser.add_argument(
        "-w",
        "--wave",
        help="""wavelength in AA""",
        type=float,
        default=None
    )

    parser.add_argument(
        "-dw",
        "--dwave",
        help="""delta wavelength in AA""",
        type=float, default=50.0
    )

    parser.add_argument(
        "-p",
        "--path",
        help="""Path to save output""",
        default=os.getcwd(),
        type=str
    )

    parser.add_argument(
        "-i",
        "--infile",
        help="""File with table of ID/RA/DEC""",
        default=None
    )

    parser.add_argument(
        "--survey",
        "-survey",
        type=str,
        help="""Data Release you want to access""",
        default="hdr2.1",
    )

    parser.add_argument(
        "--ffsky",
        "-ffsky",
        help="""Set to True to use the full frame sky sutraction.""",
        default=False,
        required=False,
        action="store_true",
    )
        
    args = parser.parse_args(argv)
    args.log = setup_logging()

    print(args)

    if args.detectid:
        config = HDRconfig()
        fileh = tb.open_file(config.detecth5)
        det_tab = fileh.root.Detections
        detectid_obj = args.detectid      
        det_row = det_tab.read_where('detectid == detectid_obj')
        shot_obj = det_row['shotid'][0]
        wave_obj = det_row['wave'][0]
        ra_cen = det_row['ra'][0]
        dec_cen = det_row['dec'][0]
        fileh.close()
    else:
        detectid_obj = args.ID
        ra_cen = args.ra
        dec_cen = args.dec
        wave_obj = args.wave
        shot_obj = args.shotid

    w = wcs.WCS(naxis=2)
    ndim = int(2 * args.gridsize / args.gridstep + 1)
    center = int(ndim / 2)
    w.wcs.crval = [ra_cen, dec_cen]
    w.wcs.crpix = [center, center]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cdelt = [-args.gridstep / 3600.0, args.gridstep / 3600.0]
    pixs = np.arange(0, ndim, 1)
    ra, dec = w.wcs_pix2world(pixs, pixs, 0)

    count = 1
    tab = Table(names=("ID", "ra", "dec"), dtype=(int, float, float))

    for i in np.arange(ndim):
        for j in np.arange(ndim):
            ra_i, dec_j = w.wcs_pix2world(i, j, 0)
            tab.add_row([count, ra_i, dec_j])
            count += 1

    coords = SkyCoord(tab["ra"], tab["dec"], unit=u.deg)
    id_list = np.array(tab["ID"])

    sources = get_spectra(coords, shotid=shot_obj, multiprocess=False,
                          survey=args.survey, ffsky=args.ffsky)

    tab.add_column(wave_obj, name='wave')
    
    output = make_line_catalog(tab, sources)

    amp_array = np.zeros([ndim, ndim])
    wave_array = np.zeros([ndim, ndim])
    sigma_array = np.zeros([ndim, ndim])
    cont_array = np.zeros([ndim, ndim])
    lf_array = np.zeros([ndim, ndim])
    
    count = 1
    for i in np.arange(ndim):
        for j in np.arange(ndim):
            sel = np.where(output["ID"] == count)[0]
            if sel:
                amp_array[j, i] = output["amp_fit"][sel]
                wave_array[j, i] = output["wave_fit"][sel]
                sigma_array[j, i] = output["sigma_fit"][sel]
                lf_array[j, i] = output["line_flux_data"][sel]
                cont_array[j,i] = output["cont_fit"][sel]
            else:
                amp_array[j, i] = np.nan
                wave_array[j, i] = np.nan
                sigma_array[j, i] = np.nan
                lf_array[j, i] = np.nan
                cont_array[j,i] = np.nan
            count += 1


    header = w.to_header()
    hdu = fits.PrimaryHDU(lf_array, header=header)
    hdu.writeto( str(detectid_obj) + '.fits' )
    
    catlib = catalogs.CatalogLibrary()

    try:
        cutout = catlib.get_cutouts(
            position=SkyCoord(ra_cen * u.deg, dec_cen * u.deg),
            radius=args.gridsize,
            aperture=None,
            dynamic=False,
            filter=["r", "g", "f606W"],
            first=True,
            allow_bad_image=False,
            allow_web=True,
        )[0]
    except:
        print("Could not get imaging for " + str(detectid_obj) )
        
    zscale = ZScaleInterval(contrast=0.5, krej=2.5)
    vmin, vmax = zscale.get_limits(values=cutout["cutout"].data)
    fig = plt.figure(figsize=(25, 6))

    ax1 = fig.add_subplot(131, projection=cutout["cutout"].wcs)
    plt.imshow(
        cutout["cutout"].data,
        vmin=vmin,
        vmax=vmax,
        origin="lower",
        cmap=plt.get_cmap("gray"),
        interpolation="none",
    )
    
    plt.text( 0.8, 0.9, cutout['instrument'] + cutout['filter'], transform=ax1.transAxes)
    
    inten =  hdu.data[ np.isfinite( hdu.data)].flatten()
    sel = inten < 3.*np.median(inten)

    (mu, sigma) = norm.fit(inten[sel])
    cont_levels = [0,1,2,3,4,5,7,10,15]

    contplot = plt.contour(
        hdu.data/sigma,
        levels=cont_levels,
        cmap="Greens",
        transform=ax1.get_transform(wcs.WCS(hdu.header)),
    )
    plt.clabel(contplot)
    plt.colorbar()
    plt.xlabel("RA")
    plt.ylabel("DEC")
    
    plt.title( 'filter='+ cutout['instrument'] + '-' +  cutout['filter'] )

    ax2 = fig.add_subplot(132, projection=cutout["cutout"].wcs)

    zscale = ZScaleInterval(contrast=0.25, krej=2.5)
    vmin, vmax = zscale.get_limits(values=lf_array)

    plt.imshow(lf_array, vmin=vmin, vmax=vmax)
    plt.colorbar(label='Line Flux (ergs/s/cm^2)')
    plt.xlabel("RA")
    plt.ylabel("DEC")
    plt.title( str(detectid_obj) + '  z = ' + '{:4.2f}   '.format( wave_obj/1216 - 1.))

    ax3 = fig.add_subplot(133, projection=cutout["cutout"].wcs)

    mask = hdu.data/sigma > 2

    plt.imshow(wave_array*mask,vmin=wave_obj-5, vmax=wave_obj+5)
    plt.colorbar(label='wave (AA)')
    plt.title('fitted central wavelength')
    plt.xlabel("RA")
    plt.ylabel("DEC")
    
    plt.savefig("im_" + str(detectid_obj) + ".png")

if __name__ == "__main__":
    main()
