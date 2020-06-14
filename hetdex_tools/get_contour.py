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

from hetdex_api.input_utils import setup_logging
from hetdex_tools.get_spec import get_spectra
from astropy.modeling import models
from specutils import Spectrum1D, SpectralRegion
from specutils.analysis import equivalent_width, line_flux
from specutils.fitting import estimate_line_parameters
from specutils.manipulation import extract_region

from elixer import catalogs 


def line_fit(index, wave_obj, sources):
    spectrum = Spectrum1D(
        flux=sources["spec"][index] * u.erg * u.cm ** -2 / u.s / u.AA,
        spectral_axis=(2.0 * np.arange(1036) + 3470.0) * u.AA,
    )

    sub_region = SpectralRegion((wave_obj - 10) * u.AA, (wave_obj + 10) * u.AA)

    sub_spectrum = extract_region(spectrum, sub_region)

    line_param = estimate_line_parameters(sub_spectrum, models.Gaussian1D())

    lf = line_flux(sub_spectrum).to(u.erg * u.cm**-2 * u.s**-1)
    
    return line_param.mean.value, line_param.amplitude.value, line_param.fwhm.value, lf.value


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

    args = parser.parse_args(argv)
    args.log = setup_logging()

    print(args)

    cat = Table.read(
        "/work/05350/ecooper/hdr1/catalogs/hdr1_sngt6pt5_for.tab",
        format="ascii"
    )

    if args.detectid:
        detectid_obj = args.detectid
        sel_obj = cat["detectid"] == detectid_obj
        shot_obj = cat["shotid"][sel_obj][0]
        wave_obj = cat["wave"][sel_obj][0]
        ra_cen = cat["ra"][sel_obj][0]
        dec_cen = cat["dec"][sel_obj][0]
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

    sources = get_spectra(coords, shotid=shot_obj, multiprocess=False, survey='hdr2')

    wave_fit = []
    amp_fit = []
    fwhm_fit = []
    lf_fit = []

    for i in np.arange(np.size(sources)):
        line_param = line_fit(i, wave_obj, sources)
        lf_fit.append(line_param[3])
        wave_fit.append(line_param[0])
        amp_fit.append(line_param[1])
        fwhm_fit.append(line_param[2] / 2.355)

    sources.add_column(Column(np.array(lf_fit)), name='lf_fit')
    sources.add_column(Column(np.array(wave_fit)), name="wave_fit")
    sources.add_column(Column(np.array(amp_fit)), name="amp_fit")
    sources.add_column(Column(np.array(fwhm_fit)), name="fwhm_fit")

    # to match up with coordinates:
    # out_table = join(tab, sources)
    # outfile = 'output_' + str(detectid_obj) + '.fits'
    # out_table.write(outfile, overwrite=True)

    amp_array = np.zeros([ndim, ndim])
    wave_array = np.zeros([ndim, ndim])
    fwhm_array = np.zeros([ndim, ndim])
    lf_array = np.zeros([ndim, ndim])
    
    count = 1
    for i in np.arange(ndim):
        for j in np.arange(ndim):
            sel = np.where(sources["ID"] == count)[0]
            if sel:
                amp_array[j, i] = sources["amp_fit"][sel]
                wave_array[j, i] = sources["wave_fit"][sel]
                fwhm_array[j, i] = sources["fwhm_fit"][sel]
                lf_array[j, i] = sources["lf_fit"][sel]
            else:
                amp_array[j, i] = np.nan
                wave_array[j, i] = np.nan
                fwhm_array[j, i] = np.nan
                lf_array[j, i] = np.nan
            count += 1

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 6))
    cs1 = ax1.contour(ra, dec, amp_array)
    cs2 = ax2.contour(
        ra, dec, wave_array,
        levels=np.arange(wave_obj - 4, wave_obj + 4, 1, dtype=int)
    )
    cs3 = ax3.contour(ra, dec, fwhm_array)

    ax1.set_xlabel("RA")
    ax1.set_ylabel("Dec")
    ax2.set_xlabel("RA")
    ax2.set_ylabel("Dec")
    ax3.set_xlabel("RA")
    ax3.set_ylabel("Dec")

    ax1.set_title("Line amplitude (10^-17 erg /cm^2/s/AA) )")
    ax2.set_title("wavelength center from fit (AA)")
    ax3.set_title("fwhm (AA)")

    ax1.clabel(cs1, inline=1, fontsize=10, fmt="%i")
    ax2.clabel(cs2, inline=1, fontsize=10, fmt="%i")
    ax3.clabel(cs3, inline=1, fontsize=10, fmt="%i")

    imfile = "contours_" + str(detectid_obj) + ".png"
    plt.savefig(imfile)

    header = w.to_header()
    hdu = fits.PrimaryHDU(lf_array, header=header)
    #    fitsfile = 'imcont_' + str(detectid_obj) + '.fits'
    #    hdu.writeto(fitsfile, overwrite=True)

    # grab image from elixer API and plot with amp contour

    catlib = catalogs.CatalogLibrary()

    cutout = catlib.get_cutouts(
        position=SkyCoord(ra_cen * u.deg, dec_cen * u.deg),
        radius=args.gridsize,
        aperture=None,
        dynamic=False,
        filter=["r", "g", "f606W"],
        first=True,
    )[0]

    if cutout:
        zscale = ZScaleInterval(contrast=0.5, krej=2.5)
        vmin, vmax = zscale.get_limits(values=cutout["cutout"].data)
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection=cutout["cutout"].wcs)
        plt.imshow(
            cutout["cutout"].data,
            vmin=vmin,
            vmax=vmax,
            origin="lower",
            cmap=plt.get_cmap("gray"),
            interpolation="none",
        )
    else:
        print("No image data found")

    plt.contour(
        hdu.data,
        alpha=0.8,
        cmap="Reds",
        transform=ax.get_transform(wcs.WCS(hdu.header)),
    )
    plt.xlabel("RA")
    plt.ylabel("DEC")
    plt.savefig("im_" + str(detectid_obj) + ".png")


if __name__ == "__main__":
    main()
