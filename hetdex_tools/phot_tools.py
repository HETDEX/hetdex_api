#!/usr/bin/env python

import numpy as np
import tables as tb
import os
import os.path as op
import time

import matplotlib
matplotlib.use("agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from photutils.isophote import EllipseGeometry, Ellipse, IsophoteList
from photutils import EllipticalAperture, SkyEllipticalAperture
from photutils import SkyCircularAperture, SkyCircularAnnulus
from photutils import aperture_photometry
from photutils.isophote import build_ellipse_model
from photutils.utils import calc_total_error
from photutils.detection import find_peaks

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import wcs
from astropy.io import fits
from astropy.table import Table, unique, join
from astropy.modeling import models, fitting
from astropy.cosmology import WMAP9 as cosmo
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.convolution import Moffat2DKernel
from astropy.convolution import convolve_models
from astropy.convolution import Gaussian1DKernel
from astropy.visualization import ZScaleInterval
from astropy.stats import sigma_clipped_stats

from hetdex_api.config import HDRconfig
from hetdex_tools.interpolate import make_narrowband_image
from hetdex_api.extract import Extract

from elixer import catalogs

plt.style.use("fivethirtyeight")

LATEST_HDR_NAME = HDRconfig.LATEST_HDR_NAME

config = HDRconfig()
surveyh5 = tb.open_file(config.surveyh5, "r")
deth5 = tb.open_file(config.detecth5, "r")
conth5 = tb.open_file(config.contsourceh5, "r")

pixscale = 0.25
imsize = 40.0


def FitCircularAperture(
    hdu=None,
    coords=None,
    radius=1.5 * u.arcsec,
    annulus=[5, 7] * u.arcsec,
    plot=False,
    plottitle=None,
):
    """
    Fit a circular aperture with either an HDU object or and
    """

    im = hdu[0].data
    error = hdu[1].data
    w = wcs.WCS(hdu[0].header)

    # create mask (set True for where you want to mask)
    im_mask = im == 0

    # define circular aperture and annulus aperture for background subtraction
    aper = SkyCircularAperture(coords, r=radius)
    aper_annulus = SkyCircularAnnulus(coords, annulus[0], annulus[1])

    mask = aper_annulus.to_pixel(w).to_mask(method="center").data
    annulus_data = aper_annulus.to_pixel(w).to_mask(method="center").multiply(im)
    annulus_data_1d = annulus_data[mask > 0]
    annulus_mask = aper_annulus.to_pixel(w).to_mask(method="center").multiply(im_mask)
    annulus_mask_1d = annulus_mask[mask > 0]

    # determine fractional fiber coverage
    apcor_im = aper.to_pixel(w).to_mask(method="center").multiply(
        im_mask
    ) / aper.to_pixel(w).to_mask(method="center").multiply(np.ones_like(im))
    apcor = np.sum(apcor_im == 0) / np.sum(np.isfinite(apcor_im))

    # get median and standard deviation in background
    mean_sigclip, median_sigclip, stddev_sigclip = sigma_clipped_stats(
        annulus_data_1d, mask=annulus_mask_1d
    )
    bkg_median = median_sigclip * aper.to_pixel(w).area * apcor
    bkg_stddev = stddev_sigclip * aper.to_pixel(w).area * apcor

    phottable = aperture_photometry(
        hdu[0].data,
        [aper, aper_annulus],
        error=hdu[1].data,
        mask=im_mask,
        wcs=wcs.WCS(hdu[0].header),
    )
    if np.abs(bkg_median) > 2 * bkg_stddev:
        flux = (phottable["aperture_sum_0"][0] - bkg_median) * u.Unit(
            "10^-17 erg cm-2 s-1"
        )
    else:
        flux = (phottable["aperture_sum_0"][0]) * u.Unit("10^-17 erg cm-2 s-1")

    flux_err = phottable["aperture_sum_err_0"][0] * u.Unit("10^-17 erg cm-2 s-1")

    if plot:
        plt.subplot(111, projection=w)
        plt.imshow(im, vmin= -1 * stddev_sigclip, vmax=4 * stddev_sigclip)
        aper.to_pixel(w).plot(color="white")  # for SkyCircularAperture
        aper_annulus.to_pixel(w).plot(color="red", linestyle="dashed")
        plt.xlabel("RA")
        plt.ylabel("Dec")
        plt.colorbar()
        if plottitle is not None:
            plt.title(plottitle)

    return flux, flux_err, bkg_stddev * u.Unit("10^-17 erg cm-2 s-1"), apcor


def get_flux_for_index(index, plot=True, convolve_image=True):
    global table
    "For a table with RA/DEC/WAVE/LINEWIDTH/SHOTID Info"
    det_info = table[index]
    det_obj = det_info["detectid"]
    coords = SkyCoord(det_info["ra"], det_info["dec"], unit="deg")
    wave_obj = det_info["wave"]
    linewidth = det_info["linewidth"]
    shot_det = det_info["shotid"]
    shotid_obj = det_info["shotid_obs"]
    fwhm = surveyh5.root.Survey.read_where("shotid == shotid_obj")["fwhm_virus"][0]

    try:
        hdu = make_narrowband_image(
            coords=coords,
            shotid=shotid_obj,
            imsize=20 * u.arcsec,
            pixscale=0.25 * u.arcsec,
            convolve_image=convolve_image,
            wave_range=[wave_obj - 2.0 * linewidth, wave_obj + 2.0 * linewidth],
            subcont=True,
            dcont=50,
            include_error=True,
        )
    except:
        return np.nan, np.nan, np.nan, np.nan

    if plot:
        plt.figure()
        plottitle = "d={} s={} s_i={}".format(det_obj, shot_det, shotid_obj)
        flux, flux_err, bkg_stddev, apcor = FitCircularAperture(
            hdu=hdu, coords=coords, plot=True, plottitle=plottitle
        )
        plt.text(
            2,
            2,
            "S/N={:3.2f}".format(flux.value / bkg_stddev.value),
            size=18,
            color="w",
        )
        plt.savefig("im_2sigma/{}_{}.png".format(det_obj, shotid_obj))
    else:
        flux, flux_err, bkg_stddev, apcor = FitCircularAperture(
            hdu=hdu, coords=coords, plot=False
        )
    return flux.value, flux_err.value, bkg_stddev.value, apcor


def get_flux_for_source(
    detectid,
    coords=None,
    radius=1 * u.arcsec,
    annulus=[5, 7] * u.arcsec,
    shotid=None,
    wave=None,
    linewidth=None,
    plot=False,
    convolve_image=False,
):

    global deth5, conth5
    
    if detectid is not None:
        detectid_obj = detectid
        
        if detectid_obj <= 2190000000:
            det_info = deth5.root.Detections.read_where('detectid == detectid_obj')[0]
            linewidth = det_info['linewidth']
            wave_obj = det_info['wave']
            redshift = wave_obj/(1216) - 1
        else:
            det_info = conth5.root.Detections.read_where('detectid == detectid_obj')[0]
            redshift = 0
            wave_obj = 4500
            
        coords_obj = SkyCoord(det_info['ra'], det_info['dec'], unit='deg')
        shotid_obj = det_info['shotid']

    if coords is not None:
        coords_obj = coords

    if shotid is not None:
        shotid_obj = shotid

    #fwhm = surveyh5.root.Survey.read_where("shotid == shotid_obj")["fwhm_virus"][0]

    if wave is not None:
        wave_obj = wave

    if linewidth is not None:
        linewidth_obj = linewidth

    try:
        hdu = make_narrowband_image(
            coords=coords_obj,
            shotid=shotid_obj,
            imsize=20 * u.arcsec,
            pixscale=0.25 * u.arcsec,
            convolve_image=convolve_image,
            dcont=50,
            wave_range=[wave_obj - 2.0 * linewidth_obj, wave_obj + 2.0 * linewidth_obj],
            subcont=True,
            include_error=True,
        )
    except:
        print('Could not make narrowband image for {}'.format(detectid_obj))
        return np.nan, np.nan, np.nan, np.nan
        
    if plot:
        plottitle = "{} {}".format(detectid_obj, shotid_obj)
        flux, flux_err, bkg_stddev, apcor = FitCircularAperture(
            hdu=hdu, coords=coords_obj, plot=True, plottitle=plottitle,
            radius=radius, annulus=annulus
        )
        plt.text(
            2,
            2,
            "S/N={:3.2f}".format(flux.value / bkg_stddev.value),
            size=18,
            color="w",
        )
    else:
        flux, flux_err, bkg_stddev, apcor = FitCircularAperture(
            hdu=hdu, coords=coords_obj,
            raidus=radius, annulus=annulus, plot=False
        )
    return flux.value, flux_err.value, bkg_stddev.value, apcor


def plot_friends(friendid, friend_cat, cutout, k=3.5, ax=None, label=True):
    """
    Plot a friends group from the fof_kdtree clustering output.
    Adapted from Niv Drory's original code
    
    Parameters
    ----------
    friends: int
        unique group identifier
    friend_cat:
        the combined group+detections catalog
    ax
        optional axis object to project coords onto
    
    """
    im = cutout["cutout"].data
    wcs = cutout["cutout"].wcs

    if ax is None:
        ax = fig.add_subplot(111, project=wcs)

    im_zscale = ZScaleInterval(contrast=0.5, krej=2.5)
    im_vmin, im_vmax = im_zscale.get_limits(values=cutout["cutout"].data)

    plt.imshow(
        cutout["cutout"].data,
        vmin=im_vmin,
        vmax=im_vmax,
        origin="lower",
        cmap=plt.get_cmap("gray"),
        interpolation="none",
    )

    sel = friend_cat["friendid"] == friendid
    group = friend_cat[sel]

    coords = SkyCoord(ra=group["icx"][0] * u.deg, dec=group["icy"][0] * u.deg)

    # get ellipse parameters in degrees
    a, b, pa, a2, b2, pa2 = (
        group["a"][0],
        group["b"][0],
        group["pa"][0],
        group["a2"][0],
        group["b2"][0],
        group["pa2"][0],
    )

    # plot detections for group
    plt.scatter(
        group["ra"],
        group["dec"],
        transform=ax.get_transform("fk5"),
        marker="x",
        color="orange",
        linewidth=1,
        s=group["flux"],
    )

    # plot and elliptical kron-like aperture representing the group
    # theta is measured
    aper_group = SkyEllipticalAperture(
        coords, a * u.deg, b * u.deg, theta=(90 - pa) * u.deg
    )
    aper_group_k = SkyEllipticalAperture(
        coords, k * a * u.deg, k * b * u.deg, theta=(90 - pa) * u.deg
    )
    aper_group.to_pixel(wcs).plot(color="blue")
    aper_group_k.to_pixel(wcs).plot(color="red")

    if label:
        # plot detecid labels
        for row in group:
            plt.text(
                row["ra"],
                row["dec"],
                str(row["detectid"]),
                transform=ax.get_transform("world"),
                fontsize=9,
                color="red",
            )


def get_line_image(
    friendid=None,
    detectid=None,
    coords=None,
    shotid=None,
    subcont=True,
    convolve_image=False,
    pixscale=pixscale,
    imsize=imsize,
    wave_range=None,
    return_coords=False,
):

    if detectid is not None:

        global deth5

        detectid_obj = detectid

        if detectid_obj <= 2190000000:
            det_info = deth5.root.Detections.read_where("detectid == detectid_obj")[0]
            linewidth = det_info["linewidth"]
            wave_obj = det_info["wave"]
            redshift = wave_obj / (1216) - 1
        else:
            det_info = conth5.root.Detections.read_where("detectid == detectid_obj")[0]
            redshift = 0
            wave_obj = 4500

        coords_obj = SkyCoord(det_info["ra"], det_info["dec"], unit="deg")

        shotid_obj = det_info["shotid"]
        fwhm = surveyh5.root.Survey.read_where("shotid == shotid_obj")["fwhm_virus"][0]
        amp = det_info["multiframe"]

        if wave_range is not None:
            wave_range_obj = wave_range
        else:
            if detectid_obj <= 2190000000:
                wave_range_obj = [wave_obj - 2 * linewidth, wave_obj + 2 * linewidth]
            else:
                wave_range_obj = [4100, 4200]

        if coords is not None:
            coords_obj = coords

        if shotid is not None:
            shotid_obj = shotid
            fwhm = surveyh5.root.Survey.read_where("shotid == shotid_obj")[
                "fwhm_virus"
            ][0]

        try:
            hdu = make_narrowband_image(
                coords=coords_obj,
                shotid=shotid_obj,
                wave_range=wave_range_obj,
                imsize=imsize * u.arcsec,
                pixscale=pixscale * u.arcsec,
                subcont=subcont,
                convolve_image=convolve_image,
                include_error=True,
            )

        except Exception:
            print("Could not make narrowband image for {}".format(detectid))
            return None

    elif friendid is not None:

        global friend_cat

        sel = friend_cat["friendid"] == friendid
        group = friend_cat[sel]
        coords_obj = SkyCoord(ra=group["icx"][0] * u.deg, dec=group["icy"][0] * u.deg)
        wave_obj = group["icz"][0]
        redshift = wave_obj / (1216) - 1
        linewidth = group["linewidth"][0]
        shotid_obj = group["shotid"][0]
        fwhm = group["fwhm"][0]
        amp = group["multiframe"][0]

        if wave_range is not None:
            wave_range_obj = wave_range
        else:
            wave_range_obj = [wave_obj - 2 * linewidth, wave_obj + 2 * linewidth]

        if shotid is not None:
            shotid_obj = shotid
            fwhm = surveyh5.root.Survey.read_where("shotid == shotid_obj")[
                "fwhm_virus"
            ][0]

        try:
            hdu = make_narrowband_image(
                coords=coords_obj,
                shotid=shotid_obj,
                wave_range=wave_range_obj,
                imsize=imsize * u.arcsec,
                pixscale=pixscale * u.arcsec,
                subcont=subcont,
                convolve_image=convolve_image,
                include_error=True,
            )
        except:
            print("Could not make narrowband image for {}".format(friendid))
            return None

    elif coords is not None:
        coords_obj = coords

        if wave_range is not None:
            wave_range_obj = wave_range
        else:
            print(
                "You need to supply wave_range=[wave_start, wave_end] for collapsed image"
            )

        if shotid is not None:
            shotid_obj = shotid
            fwhm = surveyh5.root.Survey.read_where("shotid == shotid_obj")[
                "fwhm_virus"
            ][0]
        else:
            print("Enter the shotid to use (eg. 20200123003)")

        hdu = make_narrowband_image(
            coords=coords_obj,
            shotid=shotid_obj,
            wave_range=wave_range_obj,
            imsize=imsize * u.arcsec,
            pixscale=pixscale * u.arcsec,
            subcont=subcont,
            convolve_image=convolve_image,
            include_error=True,
        )
    else:
        print("You must provide a detectid, friendid or coords/wave_range/shotid")
        return None
    if return_coords:
        return hdu, coords_obj
    else:
        return hdu


def measure_aper_flux(hdu, aper):

    # create mask (set True for where you want to mask)
    im_mask = im == 0

    # define circular aperture and annulus aperture for background subtraction
    aper = SkyCircularAperture(coords, r=radius)
    aper_annulus = SkyCircularAnnulus(coords, annulus[0], annulus[1])

    mask = aper_annulus.to_pixel(w).to_mask(method="center").data
    annulus_data = aper_annulus.to_pixel(w).to_mask(method="center").multiply(im)
    annulus_data_1d = annulus_data[mask > 0]
    annulus_mask = aper_annulus.to_pixel(w).to_mask(method="center").multiply(im_mask)
    annulus_mask_1d = annulus_mask[mask > 0]

    # determine fractional fiber coverage
    apcor_im = aper.to_pixel(w).to_mask(method="center").multiply(
        im_mask
    ) / aper.to_pixel(w).to_mask(method="center").multiply(np.ones_like(im))
    apcor = np.sum(apcor_im == 0) / np.sum(np.isfinite(apcor_im))

    # get median and standard deviation in background
    mean_sigclip, median_sigclip, stddev_sigclip = sigma_clipped_stats(
        annulus_data_1d, mask=annulus_mask_1d
    )
    bkg_median = median_sigclip * aper.to_pixel(w).area * apcor
    bkg_stddev = stddev_sigclip * aper.to_pixel(w).area * apcor

    phottable = aperture_photometry(
        hdu[0].data,
        [aper, aper_annulus],
        error=hdu[1].data,
        mask=im_mask,
        wcs=wcs.WCS(hdu[0].header),
    )
    if np.abs(bkg_median) > 2 * bkg_stddev:
        flux = (phottable["aperture_sum_0"][0] - bkg_median) * u.Unit(
            "10^-17 erg cm-2 s-1"
        )
    else:
        flux = (phottable["aperture_sum_0"][0]) * u.Unit("10^-17 erg cm-2 s-1")

    flux_err = phottable["aperture_sum_err_0"][0] * u.Unit("10^-17 erg cm-2 s-1")

    if plot:
        plt.subplot(111, projection=w)
        plt.imshow(im, vmin=0 * stddev_sigclip, vmax=3 * stddev_sigclip)
        aper.to_pixel(w).plot(color="white")  # for SkyCircularAperture
        aper_annulus.to_pixel(w).plot(color="red", linestyle="dashed")
        plt.xlabel("RA")
        plt.ylabel("Dec")
        plt.colorbar()
        if plottitle is not None:
            plt.title(plottitle)

    return flux, flux_err, bkg_stddev * u.Unit("10^-17 erg cm-2 s-1"), apcor


def fit_ellipse_for_source(
    friendid=None,
    detectid=None,
    coords=None,
    shotid=None,
    subcont=True,
    convolve_image=False,
    pixscale=pixscale,
    imsize=imsize,
    wave_range=None,
):

    if detectid is not None:

        global deth5

        detectid_obj = detectid

        if detectid_obj <= 2190000000:
            det_info = deth5.root.Detections.read_where("detectid == detectid_obj")[0]
            linewidth = det_info["linewidth"]
            wave_obj = det_info["wave"]
            redshift = wave_obj / (1216) - 1
        else:
            det_info = conth5.root.Detections.read_where("detectid == detectid_obj")[0]
            redshift = 0
            wave_obj = 4500

        coords_obj = SkyCoord(det_info["ra"], det_info["dec"], unit="deg")

        shotid_obj = det_info["shotid"]
        fwhm = surveyh5.root.Survey.read_where("shotid == shotid_obj")["fwhm_virus"][0]
        amp = det_info["multiframe"]

        if wave_range is not None:
            wave_range_obj = wave_range
        else:
            if detectid_obj <= 2190000000:
                wave_range_obj = [wave_obj - 2 * linewidth, wave_obj + 2 * linewidth]
            else:
                wave_range_obj = [4100, 4200]

        if coords is not None:
            coords_obj = coords

        if shotid is not None:
            shotid_obj = shotid
            fwhm = surveyh5.root.Survey.read_where("shotid == shotid_obj")[
                "fwhm_virus"
            ][0]

        try:
            hdu = make_narrowband_image(
                coords=coords_obj,
                shotid=shotid_obj,
                wave_range=wave_range_obj,
                imsize=imsize * u.arcsec,
                pixscale=pixscale * u.arcsec,
                subcont=subcont,
                convolve_image=convolve_image,
                include_error=True,
            )

        except:
            print("Could not make narrowband image for {}".format(detectid))
            return np.nan, np.nan

    elif friendid is not None:

        global friend_cat

        sel = friend_cat["friendid"] == friendid
        group = friend_cat[sel]
        coords_obj = SkyCoord(ra=group["icx"][0] * u.deg, dec=group["icy"][0] * u.deg)
        wave_obj = group["icz"][0]
        redshift = wave_obj / (1216) - 1
        linewidth = group["linewidth"][0]
        shotid_obj = group["shotid"][0]
        fwhm = group["fwhm"][0]
        amp = group["multiframe"][0]

        if wave_range is not None:
            wave_range_obj = wave_range
        else:
            wave_range_obj = [wave_obj - 2 * linewidth, wave_obj + 2 * linewidth]

        if shotid is not None:
            shotid_obj = shotid
            fwhm = surveyh5.root.Survey.read_where("shotid == shotid_obj")[
                "fwhm_virus"
            ][0]

        try:
            hdu = make_narrowband_image(
                coords=coords_obj,
                shotid=shotid_obj,
                wave_range=wave_range_obj,
                imsize=imsize * u.arcsec,
                pixscale=pixscale * u.arcsec,
                subcont=subcont,
                convolve_image=convolve_image,
                include_error=True,
            )
        except:
            print("Could not make narrowband image for {}".format(friendid))
            return None

    elif coords is not None:
        coords_obj = coords

        if wave_range is not None:
            wave_range_obj = wave_range
        else:
            print(
                "You need to supply wave_range=[wave_start, wave_end] for collapsed image"
            )

        if shotid is not None:
            shotid_obj = shotid
            fwhm = surveyh5.root.Survey.read_where("shotid == shotid_obj")[
                "fwhm_virus"
            ][0]
        else:
            print("Enter the shotid to use (eg. 20200123003)")

        hdu = make_narrowband_image(
            coords=coords_obj,
            shotid=shotid_obj,
            wave_range=wave_range_obj,
            imsize=imsize * u.arcsec,
            pixscale=pixscale * u.arcsec,
            subcont=subcont,
            convolve_image=convolve_image,
            include_error=True,
        )
    else:
        print("You must provide a detectid, friendid or coords/wave_range/shotid")
        return np.nan, np.nan

    w = wcs.WCS(hdu[0].header)

    if friendid is not None:

        sel_friend_group = friend_cat["friendid"] == friendid
        group = friend_cat[sel_friend_group]
        eps = 1 - group["a2"][0] / group["b2"][0]
        pa = group["pa"][0] * np.pi / 180.0 - 90
        sma = group["a"][0] * 3600 / pixscale

        coords = SkyCoord(ra=group["icx"][0] * u.deg, dec=group["icy"][0] * u.deg)
        wave_obj = group["icz"][0]
        redshift = wave_obj / (1216) - 1
        linewidth = np.nanmedian(group["linewidth"])
        shotid_obj = group["shotid"][0]
        fwhm = group["fwhm"][0]

        geometry = EllipseGeometry(
            x0=w.wcs.crpix[0], y0=w.wcs.crpix[0], sma=sma, eps=eps, pa=pa
        )
    else:
        geometry = EllipseGeometry(
            x0=w.wcs.crpix[0], y0=w.wcs.crpix[0], sma=20, eps=0.2, pa=20.0
        )

    geometry = EllipseGeometry(
        x0=w.wcs.crpix[0], y0=w.wcs.crpix[0], sma=20, eps=0.2, pa=20.0
    )
    # geometry.find_center(hdu.data)
    # aper = EllipticalAperture((geometry.x0, geometry.y0), geometry.sma,
    #                          geometry.sma*(1 - geometry.eps), geometry.pa)

    # plt.imshow(hdu.data, origin='lower')
    # aper.plot(color='white')

    ellipse = Ellipse(hdu[0].data)
    isolist = ellipse.fit_image()
    iso_tab = isolist.to_table()

    if len(iso_tab) == 0:
        geometry.find_center(hdu[0].data, verbose=False, threshold=0.5)
        ellipse = Ellipse(hdu[0].data, geometry)
        isolist = ellipse.fit_image()
        iso_tab = isolist.to_table()

    if len(iso_tab) == 0:
        return np.nan, np.nan, np.nan

        try:
            # compute iso's manually in steps of 3 pixels
            ellipse = Ellipse(hdu[0].data)  # reset ellipse
            iso_list = []
            for sma in np.arange(1, 60, 2):
                iso = ellipse.fit_isophote(sma)
                if np.isnan(iso.intens):
                    # print('break at {}'.format(sma))
                    break
                else:
                    iso_list.append(iso)
            isolist = IsophoteList(iso_list)
            iso_tab = isolist.to_table()
        except:
            return np.nan, np.nan, np.nan

    try:
        model_image = build_ellipse_model(hdu[0].data.shape, isolist)
        residual = hdu[0].data - model_image
    except:
        return np.nan, np.nan, np.nan

    sma = iso_tab["sma"] * pixscale
    const_arcsec_to_kpc = cosmo.kpc_proper_per_arcmin(redshift).value / 60.0

    def arcsec_to_kpc(sma):
        dist = const_arcsec_to_kpc * sma
        return dist

    def kpc_to_arcsec(dist):
        sma = dist / const_arcsec_to_kpc
        return sma

    dist_kpc = (
        sma * u.arcsec.to(u.arcmin) * u.arcmin * cosmo.kpc_proper_per_arcmin(redshift)
    )
    dist_arcsec = kpc_to_arcsec(dist_kpc)

    # print(shotid_obj, fwhm)
    # s_exp1d = models.Exponential1D(amplitude=0.2, tau=-50)

    alpha = 3.5
    s_moffat = models.Moffat1D(
        amplitude=1,
        gamma=(0.5 * fwhm) / np.sqrt(2 ** (1.0 / alpha) - 1.0),
        x_0=0.0,
        alpha=alpha,
        fixed={"amplitude": False, "x_0": True, "gamma": True, "alpha": True},
    )

    s_init = models.Exponential1D(amplitude=0.2, tau=-50)

    fit = fitting.LevMarLSQFitter()
    s_r = fit(s_init, dist_kpc, iso_tab["intens"])

    # Fitting can be done using the uncertainties as weights.
    # To get the standard weighting of 1/unc^2 for the case of
    # Gaussian errors, the weights to pass to the fitting are 1/unc.
    # fitted_line = fit(line_init, x, y, weights=1.0/yunc)

    # s_r = fit(s_init, dist_kpc, iso_tab['intens'])#, weights=iso_tab['intens']/iso_tab['intens_err'] )

    print(s_r)
    try:
        r_n = -1.0 * s_r.tau  # _0 #* const_arcsec_to_kpc
    except:
        r_n = np.nan  # r_n = -1. * s_r.tau_0
    try:
        sel_iso = np.where(dist_kpc >= 2 * r_n)[0][0]
    except:
        sel_iso = -1

    aper = EllipticalAperture(
        (isolist.x0[sel_iso], isolist.y0[sel_iso]),
        isolist.sma[sel_iso],
        isolist.sma[sel_iso] * (1 - isolist.eps[sel_iso]),
        isolist.pa[sel_iso],
    )

    phottable = aperture_photometry(hdu[0].data, aper, error=hdu[1].data)
    flux = phottable["aperture_sum"][0] * 10 ** -17 * u.erg / (u.cm ** 2 * u.s)
    flux_err = phottable["aperture_sum_err"][0] * 10 ** -17 * u.erg / (u.cm ** 2 * u.s)

    lum_dist = cosmo.luminosity_distance(redshift).to(u.cm)
    lum = flux * 4.0 * np.pi * lum_dist ** 2
    lum_err = flux_err * 4.0 * np.pi * lum_dist ** 2

    if detectid:
        name = detectid
    elif friendid:
        name = friendid

    # Get Image data from Elixer
    catlib = catalogs.CatalogLibrary()
    try:
        cutout = catlib.get_cutouts(
            position=coords_obj,
            side=imsize,
            aperture=None,
            dynamic=False,
            filter=["r", "g", "f606W"],
            first=True,
            allow_bad_image=False,
            allow_web=True,
        )[0]
    except:
        print("Could not get imaging for " + str(name))

    zscale = ZScaleInterval(contrast=0.5, krej=1.5)
    vmin, vmax = zscale.get_limits(values=hdu[0].data)

    fig = plt.figure(figsize=(20, 12))
    fig.suptitle(
        "{}  ra={:3.2f}, dec={:3.2f}, wave={:5.2f}, z={:3.2f}, mf={}".format(
            name, coords_obj.ra.value, coords_obj.dec.value, wave_obj, redshift, amp
        ),
        fontsize=22,
    )

    ax1 = fig.add_subplot(231, projection=w)
    plt.imshow(hdu[0].data, vmin=vmin, vmax=vmax)
    plt.xlabel("RA")
    plt.ylabel("Dec")
    plt.colorbar()
    plt.title("Image summed across 4*linewidth")

    ax2 = fig.add_subplot(232, projection=w)
    plt.imshow(model_image, vmin=vmin, vmax=vmax)
    plt.xlabel("RA")
    plt.ylabel("Dec")
    plt.colorbar()
    plt.title("model")

    ax3 = fig.add_subplot(233, projection=w)
    plt.imshow(residual, vmin=vmin, vmax=vmax)
    plt.xlabel("RA")
    plt.ylabel("Dec")
    plt.colorbar()
    plt.title("residuals (image-model)")
    # fig = plt.figure(figsize=(10,5))

    im_zscale = ZScaleInterval(contrast=0.5, krej=2.5)
    im_vmin, im_vmax = im_zscale.get_limits(values=cutout["cutout"].data)

    ax4 = fig.add_subplot(234, projection=cutout["cutout"].wcs)

    plt.imshow(
        cutout["cutout"].data,
        vmin=im_vmin,
        vmax=im_vmax,
        origin="lower",
        cmap=plt.get_cmap("gray"),
        interpolation="none",
    )

    plt.text(
        0.8,
        0.9,
        cutout["instrument"] + cutout["filter"],
        transform=ax4.transAxes,
        fontsize=20,
        color="w",
    )
    plt.contour(hdu[0].data, transform=ax4.get_transform(w))
    plt.xlabel("RA")
    plt.ylabel("Dec")
    aper.plot(
        color="white", linestyle="dashed", linewidth=2, transform=ax4.get_transform(w)
    )

    ax5 = fig.add_subplot(235)
    plt.errorbar(
        dist_kpc.value,
        iso_tab["intens"],
        yerr=iso_tab["intens_err"] * iso_tab["intens"],
        linestyle="none",
        marker="o",
        label="Lya SB profile",
    )

    plt.plot(dist_kpc, s_r(dist_kpc), color="r", label="Lya exp SB model", linewidth=2)

    plt.xlabel("Semi-major axis (kpc)")
    # plt.xlabel('Semi-major axis (arcsec)')
    plt.ylabel("Flux ({})".format(10 ** -17 * (u.erg / (u.s * u.cm ** 2))))
    plt.text(0.4, 0.7, "r_n={:3.2f}".format(r_n), transform=ax5.transAxes, fontsize=16)
    plt.text(
        0.4, 0.6, "L_lya={:3.2e}".format(lum), transform=ax5.transAxes, fontsize=16
    )
    secax = ax5.secondary_xaxis("top", functions=(kpc_to_arcsec, kpc_to_arcsec))
    secax.set_xlabel("Semi-major axis (arcsec)")
    # secax.set_xlabel('Semi-major axis (kpc)')
    plt.xlim(0, 100)

    # plt.plot(sma, s_r(sma), label='moffat psf')

    # plt.plot(dist_kpc.value, s1(kpc_to_arcsec(dist_kpc.value)),
    #        linestyle='dashed', linewidth=2,
    #         color='green', label='PSF seeing:{:3.2f}'.format(fwhm))

    # These two are the exact same
    # s1 = models.Moffat1D()
    # s1.amplitude = iso_tab['intens'][0]
    # alpha=3.5
    # s1.gamma = 0.5*(fwhm*const_arcsec_to_kpc)/ np.sqrt(2 ** (1.0 / alpha) - 1.0)
    # s1.alpha = alpha
    # plt.plot(r_1d, moffat_1d, color='orange')

    #    plt.plot(dist_kpc.value, (s1(dist_kpc.value)),
    #            linestyle='dashed', linewidth=2,
    #             color='blue', label='PSF seeing:{:3.2f}'.format(fwhm))

    E = Extract()
    E.load_shot(shotid_obj)
    moffat_psf = E.moffat_psf(seeing=fwhm, boxsize=imsize, scale=pixscale)
    moffat_shape = np.shape(moffat_psf)
    xcen = int(moffat_shape[1] / 2)
    ycen = int(moffat_shape[2] / 2)
    moffat_1d = (
        moffat_psf[0, xcen:-1, ycen] / moffat_psf[0, xcen, ycen] * iso_tab["intens"][0]
    )
    r_1d = moffat_psf[1, xcen:-1, ycen]
    E.close()

    plt.plot(
        arcsec_to_kpc(pixscale * np.arange(80)),
        iso_tab["intens"][0] * (moffat_psf[0, 80:-1, 80] / moffat_psf[0, 80, 80]),
        linestyle="dashed",
        color="green",
        label="PSF seeing:{:3.2f}".format(fwhm),
    )
    plt.legend()

    if friendid is not None:
        ax6 = fig.add_subplot(236, projection=cutout["cutout"].wcs)
        plot_friends(friendid, friend_cat, cutout, ax=ax6, label=False)
    plt.savefig("fit2d_{}.png".format(name))

    # filename = 'param_{}.txt'.format(name)
    # np.savetxt(filename, (r_n.value, lum.value))

    return r_n, lum, lum_err


def get_sn_for_aperture_range(
    hdu, find_peak=False, r_list=np.arange(0.5, 3.1, 0.1), dtype=float
):

    im = hdu[0].data
    im_mask = im == 0
    mean, median, std = sigma_clipped_stats(im[im_mask == 0], sigma=3.0)
    threshold = 3 * std

    if find_peak:
        r = np.sqrt(hdu[2].data ** 2 + hdu[3].data ** 2)
        searchregion = r < 3
        tbl = find_peaks(
            im,
            threshold,
            footprint=searchregion,
            mask=im_mask,
            wcs=wcs.WCS(hdu[0].header),
        )
        if tbl is None:
            coords_peak = None
        else:
            coords_peak = tbl["skycoord_peak"][np.argmax(tbl["peak_value"])]
    else:
        coords_peak = None

    if coords_peak is None:
        coords_center = SkyCoord(
            ra=hdu[0].header["CRVAL1"], dec=hdu[0].header["CRVAL2"], unit="deg"
        )
    else:
        coords_center = coords_peak

    flux_list = []
    flux_err_list = []
    bkg_list = []
    apcor_list = []

    for r in r_list:

        flux, flux_err, bkg_stddev, apcor = FitCircularAperture(
            hdu=hdu,
            coords=coords_center,
            plot=False,
            radius=r * u.arcsec,
            annulus=[r * 2, r * 3] * u.arcsec,
        )
        flux_list.append(flux.value)
        flux_err_list.append(flux_err.value)
        bkg_list.append(bkg_stddev.value)
        apcor_list.append(apcor)

    flux_arr = np.array(flux_list)
    flux_err_arr = np.array(flux_err_list)
    bkg_arr = np.array(bkg_list)
    apcor_arr = np.array(apcor_list)
    sn = flux_arr / bkg_arr

    return r_list, sn, coords_center


def fit_growing_aperture(detectid,
                         shotid=None,
                         plot=True,
                         img_dir='line_images'):

    global imsize
    if plot:
        if op.exists(img_dir):
            pass
        else:
            os.makedirs(img_dir)
    try:
        if shotid is None:
            hdu, coords = get_line_image(detectid=detectid,
                                         imsize=imsize/2,
                                         return_coords=True)
        else:
            hdu, coords = get_line_image(detectid=detectid,
                                         shotid=shotid,
                                         imsize=imsize/2,
                                         return_coords=True)
        pixsize=imsize/2/pixscale
        
    except:
        print('Could not get image for {}'.format(detectid))

    r_list, sn, coords_center = get_sn_for_aperture_range(hdu)

    index_max = np.argmax(sn)
    sn_max = sn[index_max]
    sn_upper_limit = sn_max
    r_snmax = r_list[index_max]

    try:
        index_2sigma = np.max(np.where(sn > 2))
        sn_2sigma = sn[index_2sigma]
        r_2sigma = r_list[index_2sigma]
    except:
        sn_2sigma = np.nan
        r_2sigma = np.nan

    sn_upper_limit = sn[-1]

    if sn[-1] > 2:
        # increase image cutout and aperture search to larger radii
        if shotid is None:
            hdu, coords = get_line_image(detectid=detectid,
                                         imsize=imsize,
                                         shotid=shotid,
                                         return_coords=True)
        else:
            hdu, coords = get_line_image(detectid=detectid,
                                         shotid=shotid,
                                         imsize=imsize,
                                         return_coords=True)

        r_list, sn, coords_center = get_sn_for_aperture_range(
            hdu, r_list=np.arange(1.5, 8.2, 0.2, dtype=float)
        )
        
        pixsize=imsize/pixscale
        
        try:
            index_2sigma = np.max(
                np.where((sn > 2) & (sn < sn_upper_limit) & np.isfinite(sn))
            )
            # print(index_2sigma)
            sn_2sigma = sn[index_2sigma]
            r_2sigma = r_list[index_2sigma]

            if np.max(sn) > sn_max:
                index_max = np.argmax(sn)
                sn_max = sn[index_max]
                r_snmax = r_list[index_max]
        except:
            pass

    # stop growing at 8 arcsec for now
    # get flux info at sn_max and sn2sigma

    flux_snmax, flux_err_snmax, bkg_stddev_snmax, apcor_snmax = FitCircularAperture(
                hdu=hdu,
                coords=coords_center,
                plot=False,
                radius=r_snmax * u.arcsec,
                annulus=[r_snmax, r_snmax * 3] * u.arcsec,
            )

    if np.isfinite(r_2sigma):
        if plot:
            plt.figure(figsize=(6,6))
            if shotid is not None:
                plottitle = "{} {} r_eff={:3.2}".format(detectid,
                                                        shotid,
                                                        r_2sigma)
            else:
                plottitle = " {}  r_eff={:3.2} r_snmax={:3.2}".format(detectid,
                                                                      r_2sigma,
                                                                      r_snmax)
            (
                flux_2sigma,
                flux_err_2sigma,
                bkg_stddev_2sigma,
                apcor_2sigma,
            ) = FitCircularAperture(
                hdu=hdu,
                coords=coords_center,
                radius=r_2sigma * u.arcsec,
                annulus=[r_2sigma * 2, r_2sigma * 3] * u.arcsec,
                plot=True,
                plottitle=plottitle,
            )
            plt.text(
                2,
                2,
                "S/Nmax={:3.2f}".format(flux_snmax.value / bkg_stddev_snmax.value),
                size=18,
                color="w",
            )
            plt.text(
                0.55*pixsize,
                2,
                "S/Nsig={:3.2f}".format(flux_2sigma.value / bkg_stddev_2sigma.value),
                size=18,
                color="w",
            )
            if shotid is not None:
                plt.savefig(op.join(img_dir,"{}_{}.png".format(detectid,
                                                               shotid)))
            else:
                plt.savefig(op.join(img_dir,"{}.png".format(detectid)))
        else:
            (
                flux_2sigma,
                flux_err_2sigma,
                bkg_stddev_2sigma,
                apcor_2sigma,
            ) = FitCircularAperture(
                hdu=hdu,
                coords=coords_center,
                radius=r_2sigma * u.arcsec,
                annulus=[r_2sigma * 2, r_2sigma * 3] * u.arcsec,
                plot=False
            )
    else:
        #plot SNmax if r_2sigma is not defined
        plt.figure()
        if shotid is not None:
            plottitle = "{} {} r_eff={:3.2}".format(detectid,
                                                    shotid,
                                                    r_2sigma)
        else:
            plottitle = " {}  r_eff={:3.2} r_snmax={:3.2}".format(detectid,
                                                                  r_2sigma,
                                                                  r_snmax)
        (
            flux_snmax,
            flux_err_snmax,
            bkg_stddev_snmax,
            apcor_snmax,
        ) = FitCircularAperture(
            hdu=hdu,
            coords=coords_center,
            radius=r_snmax * u.arcsec,
            annulus=[r_snmax * 2, r_snmax * 3] * u.arcsec,
            plot=True,
            plottitle=plottitle,
        )
        plt.text(
            2,
            2,
            "S/N={:3.2f}".format(flux_snmax.value / bkg_stddev_snmax.value),
            size=18,
            color="w",
        )

        if shotid is not None:
            plt.savefig(op.join(img_dir,"{}_{}.png".format(detectid, shotid)))
        else:
            plt.savefig(op.join(img_dir,"{}.png".format(detectid)))
                                                                    
        flux_2sigma = np.nan
        flux_err_2sigma = np.nan
        bkg_stddev_2sigma = np.nan
        apcor_2sigma = np.nan

    return (
        r_2sigma,
        sn_2sigma,
        r_snmax,
        sn_max,
        flux_2sigma,
        flux_err_2sigma,
        bkg_stddev_2sigma,
        apcor_2sigma,
        flux_snmax,
        flux_err_snmax,
        bkg_stddev_snmax,
        apcor_snmax,
    )


def make_im_catalog(detlist,
                    shotlist=None,
                    filename="imflux.tab",
                    save=True,
                    plot=True,
                    img_dir='line_images'):

    if plot:
        if op.exists(img_dir):
            pass
        else:
            os.makedirs(img_dir)
                                                            
    t0 = time.time()
    imflux = Table(
        names=[
            "r_2sigma",
            "sn_2sigma",
            "r_snmax",
            "sn_max",
            "flux_2sigma",
            "flux_err_2sigma",
            "bkg_stddev_2sigma",
            "apcor_2sigma",
            "flux_snmax",
            "flux_err_snmax",
            "bkg_stddev_snmax",
            "apcor_snmax",
        ]
    )

    detcheck = []
    for i, det in enumerate( detlist):

        detcheck.append(det)
        
        try:
            if shotlist is None:
                (
                    r_2sigma,
                    sn_2sigma,
                    r_snmax,
                    sn_max,
                    flux_2sigma,
                    flux_err_2sigma,
                    bkg_stddev_2sigma,
                    apcor_2sigma,
                    flux_snmax,
                    flux_err_snmax,
                    bkg_stddev_snmax,
                    apcor_snmax,
                ) = fit_growing_aperture(det,
                                         plot=plot,
                                         img_dir=img_dir)
            else:
                shotid = shotlist[i]
                (
                    r_2sigma,
                    sn_2sigma,
                    r_snmax,
                    sn_max,
                    flux_2sigma,
                    flux_err_2sigma,
                    bkg_stddev_2sigma,
                    apcor_2sigma,
                    flux_snmax,
                    flux_err_snmax,
                    bkg_stddev_snmax,
                    apcor_snmax,
                ) = fit_growing_aperture(det,
                                         plot=plot,
                                         shotid=shotid,
                                         img_dir=img_dir)

            imflux.add_row(
                [
                    np.float32(r_2sigma),
                    np.float32(sn_2sigma),
                    np.float32(r_snmax),
                    np.float32(sn_max),
                    np.float32(flux_2sigma),
                    np.float32(flux_err_2sigma),
                    np.float32(bkg_stddev_2sigma),
                    np.float32(apcor_2sigma),
                    np.float32(flux_snmax),
                    np.float32(flux_err_snmax),
                    np.float32(bkg_stddev_snmax),
                    np.float32(apcor_snmax),
                ]
            )
        except:
            imflux.add_row(
                [
                 np.float32(np.nan),
                 np.float32(np.nan),
                 np.float32(np.nan),
                 np.float32(np.nan),
                 np.float32(np.nan),
                 np.float32(np.nan),
                 np.float32(np.nan),
                 np.float32(np.nan),
                 np.float32(np.nan),
                 np.float32(np.nan),
                 np.float32(np.nan),
                 np.float32(np.nan),
             ]
            )

    print("Done in {:4.2f} s".format(time.time() - t0))
    
    imflux.add_column(detcheck, name='detectid', index=0)

    if shotlist is not None:
        imflux.add_column(shotlist, name="shotid_obs", index=1)

    imflux.write(filename, format="ascii", overwrite=True)

    if np.sum(np.array(detcheck)-np.array(detlist)) >0:
        print('detlist failed detcheck')
        
    return imflux
