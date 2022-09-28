#!/usr/bin/env python
# coding: utf-8
#!pip3 install pyimfit
#!cd /home/jovyan/software/hetdex_api; git stash; git pull

# if running on wave_group_id
# do_pyimfit(objectid)
# do_pyimfit(detectid=objectid, wave_range=[4000,5000], subcont=False, star=True)
# command line
# python3 lya_pyimfit.py --detectid 3090001360 --wave_range 4000 5000 --subcont False --star

import numpy as np
import tables as tb
import os.path as op

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import wcs
from astropy.stats import sigma_clipped_stats

from photutils import EllipticalAnnulus, CircularAnnulus
from photutils import aperture_photometry

from astropy.cosmology import Planck18 as cosmo

from hetdex_api.config import HDRconfig
from hetdex_tools.interpolate import make_narrowband_image
from hetdex_api.extract import Extract
from hetdex_tools.get_spec import get_spectra

from elixer import catalogs

import pyimfit
import argparse as ap

import matplotlib

matplotlib.use("agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.visualization import ZScaleInterval
from matplotlib import gridspec

plot_both_models = False
plot_log = True

LATEST_HDR_NAME = HDRconfig.LATEST_HDR_NAME

config = HDRconfig()
surveyh5 = tb.open_file(config.surveyh5, "r")
deth5 = tb.open_file(config.detecth5, "r")
conth5 = tb.open_file(config.contsourceh5, "r")

version = "3.0.1"
catfile = op.join(
    config.hdr_dir["hdr3"], "catalogs", "source_catalog_" + version + ".fits"
)
source_table = None


def do_pyimfit(
    wave_group_id=None,
    detectid=None,
    coords=None,
    shotid=None,
    subcont=True,
    convolve_image=False,
    pixscale=0.25,
    imsize=30.0,
    wave_range=None,
    star=False,
):

    if detectid is not None:

        global deth5

        detectid_obj = detectid
        name = detectid
        if detectid_obj < 3090000000:
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
            pass  # return np.nan, np.nan

    elif wave_group_id is not None:

        global source_table

        if source_table is None:
            source_table = Table.read(catfile)

        sel = source_table["wave_group_id"] == wave_group_id
        name = wave_group_id

        group = source_table[sel]

        coords_obj = SkyCoord(
            ra=group["wave_group_ra"][0] * u.deg, dec=group["wave_group_dec"][0] * u.deg
        )
        wave_obj = group["wave_group_wave"][0]
        redshift = wave_obj / (1216) - 1
        linewidth = group["linewidth"][0]
        shotid_obj = group["shotid"][0]
        fwhm = group["fwhm"][0]
        amp = group["multiframe"][0]

        if wave_range is not None:
            wave_range_obj = wave_range
        else:
            wave_range_obj = [wave_obj - 2 * linewidth, wave_obj + 2 * linewidth]

        fwhm = surveyh5.root.Survey.read_where("shotid == shotid_obj")["fwhm_virus"][0]

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

    elif coords is not None:
        coords_obj = coords
        name = "ra={:4.3f}_dec={4.3f}".format(coords.ra, coords.dec)
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
        print("You must provide a detectid, wave_group_id or coords/wave_range/shotid")

    w = wcs.WCS(hdu[0].header)

    spec_table = get_spectra(
        coords_obj, shotid=shotid_obj, multiprocess=False, loglevel="WARNING"
    )
    const_arcsec_to_kpc = cosmo.kpc_proper_per_arcmin(redshift).value / 60.0

    def arcsec_to_kpc(sma):
        dist = const_arcsec_to_kpc * sma
        return dist

    def kpc_to_arcsec(dist):
        sma = dist / const_arcsec_to_kpc
        return sma

    pix_to_kpc = (
        pixscale
        * u.arcsec.to(u.arcmin)
        * u.arcmin
        * cosmo.kpc_proper_per_arcmin(redshift)
    ).value

    image_data = hdu[0].data

    E = Extract()
    E.load_shot(shotid_obj)
    moffat_psf = E.moffat_psf(seeing=fwhm, boxsize=imsize, scale=pixscale)
    moffat_shape = np.shape(moffat_psf)
    psf_xcen = int(moffat_shape[1] / 2)
    psf_ycen = int(moffat_shape[2] / 2)
    E.close()

    xpix, ypix = np.shape(image_data)
    xcen = int(xpix / 2)
    ycen = int(ypix / 2)
    dx = 20
    dy = 20

    mask = image_data == 0

    sky = np.median(image_data[~mask])
    mean, median, stddev = sigma_clipped_stats(hdu[0].data[~mask], sigma=2, maxiters=5)
    masked_data = hdu[0].data
    masked_data = np.ma.masked_where((masked_data - mean) < 3 * stddev, masked_data)

    model_desc = pyimfit.SimpleModelDescription()

    model_desc.x0.setValue(xcen, [xcen - dx, xcen + dx])
    model_desc.y0.setValue(ycen, [ycen - dy, ycen + dx])

    expmodel = pyimfit.make_imfit_function("Exponential")
    # set initial values, lower and upper limits for central surface brightness I_0, scale length h;
    # specify that ellipticity is to remain fixed
    expmodel.I_0.setValue(10, [0, 1000])
    expmodel.h.setValue(10, [0, 100])
    expmodel.PA.setValue(40, [0, 180])
    expmodel.ell.setValue(0.5, [0, 1])

    model_desc.addFunction(expmodel)

    imfit_fitter = pyimfit.Imfit(model_desc, psf=moffat_psf[0])

    imfit_fitter.fit(image_data, mask=mask, error=hdu[1].data)
    chi2 = imfit_fitter.reducedFitStatistic

    fit_params = imfit_fitter.getRawParameters()

    exp_x0 = fit_params[0]
    exp_y0 = fit_params[1]
    exp_PA = fit_params[2]
    exp_ell = fit_params[3]
    exp_I_0 = fit_params[4]
    exp_h = fit_params[5]

    moffat_model_desc = pyimfit.SimpleModelDescription()

    moffat_model_desc.x0.setValue(xcen, [xcen - dx, xcen + dx])
    moffat_model_desc.y0.setValue(ycen, [ycen - dy, ycen + dx])

    # create an Exponential image function, then define the parameter initial values and limits
    moffatmodel = pyimfit.make_imfit_function("Moffat")
    # set initial values, lower and upper limits for central surface brightness I_0, scale length h;
    # specify that ellipticity is to remain fixed
    moffatmodel.PA.setValue(0, fixed=True)
    moffatmodel.ell.setValue(0, fixed=True)
    moffatmodel.fwhm.setValue(fwhm / pixscale, fixed=True)
    moffatmodel.beta.setValue(3.5, fixed=True)
    moffatmodel.I_0.setValue(10, [0.001, 10000])
    moffat_model_desc.addFunction(moffatmodel)

    imfit_fitter_moffat = pyimfit.Imfit(moffat_model_desc)
    imfit_fitter_moffat.fit(image_data, mask=mask, error=hdu[1].data)
    chi2_moffat = imfit_fitter_moffat.reducedFitStatistic

    fit_params_moffat = imfit_fitter_moffat.getRawParameters()

    moffat_x0 = fit_params_moffat[0]
    moffat_y0 = fit_params_moffat[1]

    r_n = (
        exp_h
        * pixscale
        * u.arcsec.to(u.arcmin)
        * u.arcmin
        * cosmo.kpc_proper_per_arcmin(redshift)
    ).value

    output = np.append(r_n, chi2)
    output = np.append(output, chi2_moffat)
    output = np.append(output, fit_params)

    model_im = imfit_fitter.getModelImage()
    moffat_im = imfit_fitter_moffat.getModelImage()

    dr = 1
    sb = []
    sb_error = []
    sb_model_im = []
    sb_moffat_im = []

    r_in_array = np.arange(1, xcen, dr)
    r_in_array = np.logspace(-3, np.log10(xcen))
    for r_in in r_in_array:

        if star:
            r_out = r_in + dr
            aper = CircularAnnulus((moffat_x0, moffat_y0), r_in, r_out)
        else:
            a_in = r_in
            a_out = r_in + dr
            b_out = a_out * (1 - exp_ell ** 2)

            aper = EllipticalAnnulus((exp_x0, exp_y0), a_in, a_out, b_out, theta=exp_PA)

        phot_table = aperture_photometry(
            image_data, aper, mask=mask, error=hdu[1].data + stddev
        )

        sb.append(phot_table["aperture_sum"][0] / aper.area)

        sb_error.append(phot_table["aperture_sum_err"][0] / aper.area)

        phot_table = aperture_photometry(model_im, aper)
        sb_model_im.append(phot_table["aperture_sum"][0] / aper.area)

        phot_table = aperture_photometry(moffat_im, aper)
        sb_moffat_im.append(phot_table["aperture_sum"][0] / aper.area)

    sb = 10 ** -17 * np.array(sb) / (pixscale ** 2)
    sb_error = 10 ** -17 * np.array(sb_error) / (pixscale ** 2)
    sb_model_im = 10 ** -17 * np.array(sb_model_im) / (pixscale ** 2)
    sb_moffat_im = 10 ** -17 * np.array(sb_moffat_im) / (pixscale ** 2)

    plt.style.use("default")
    plt.rcParams["axes.linewidth"] = 2
    plt.rcParams.update({"font.size": 12})

    plt.rcParams["lines.linewidth"] = 2

    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.fontset"] = "dejavuserif"
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["xtick.labelsize"] = 12.0
    plt.rcParams["ytick.labelsize"] = 12.0

    w1 = 3540
    w2 = 5450

    fig = plt.figure(figsize=(14, 14))
    gs = gridspec.GridSpec(3, 3)
    fig.suptitle(
        "{}  ra={:3.2f}, dec={:3.2f}, wave={:5.1f}-{:5.1f}, z={:3.2f}, mf={}".format(
            name,
            coords_obj.ra.value,
            coords_obj.dec.value,
            wave_range_obj[0],
            wave_range_obj[1],
            redshift,
            amp,
        ),
        fontsize=16,
    )

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

    im_zscale = ZScaleInterval(contrast=0.5, krej=2.5)
    im_vmin, im_vmax = im_zscale.get_limits(values=cutout["cutout"].data)

    ax4 = fig.add_subplot(gs[0, 0], projection=w)  # projection=cutout['cutout'].wcs)

    plt.imshow(
        cutout["cutout"].data,
        vmin=im_vmin,
        vmax=im_vmax,
        origin="lower",
        cmap=plt.get_cmap("gray"),
        interpolation="none",
        transform=ax4.get_transform(cutout["cutout"].wcs),
    )

    plt.text(
        0.6,
        0.9,
        cutout["instrument"] + cutout["filter"],
        transform=ax4.transAxes,
        fontsize=20,
        color="w",
    )

    mean, median, stddev = sigma_clipped_stats(hdu[0].data[~mask], sigma=2, maxiters=5)
    masked_data = hdu[0].data
    masked_data = np.ma.masked_where((masked_data - mean) < 3 * stddev, masked_data)

    plt.imshow(
        masked_data, alpha=0.6
    )  # , cmap=plt.get_cmap("Blues"))#, transform=ax4.get_transform(w))

    plt.xlabel("RA")
    plt.ylabel("Dec")
    lon = ax4.coords[0]
    lat = ax4.coords[1]
    lon.set_axislabel("RA", minpad=0.2)
    lat.set_axislabel("Dec", minpad=-0.2)
    lon.set_ticklabel(exclude_overlapping=True)

    # aper.plot(color='white', linestyle='dashed', linewidth=4, transform=ax4.get_transform(w))

    ax3 = fig.add_subplot(gs[1, :])
    plt.plot(
        spec_table["wavelength"][0],
        spec_table["spec"][0],
        linewidth=1.2,
        color="tab:blue",
    )  # , yerr=spec_table['spec1d_err'])
    plt.xlabel(r"$\lambda$ ($\AA$)")
    plt.ylabel("f$_\\lambda$ (10$^{-17}$ ergs/s/cm$^2$/$\AA$)")

    plt.xlim(w1, w2)
    # plt.axvline(x=det_info['wave'], color='r', linestyle='dashed')
    selw = (spec_table["wavelength"] > w1) & (spec_table["wavelength"] < w2)
    y2 = np.max(spec_table["spec"][selw])
    y1 = np.min(spec_table["spec"][selw])
    plt.ylim(y1, 1.1 * y2)
    plt.axhline(0, color="tab:grey", linestyle="dashed")
    plt.text(
        0.05,
        0.7,
        "z={:6.4f}".format(redshift),
        transform=ax3.transAxes,
        color="tab:red",
        fontsize=20,
    )

    plt.axvspan(wave_range_obj[0], wave_range_obj[1], color="orange", alpha=0.4)

    ax5 = fig.add_subplot(gs[0, 1:])

    if star:
        plt.errorbar(
            r_in_array * pixscale,
            sb,
            yerr=sb_error,
            label="data",
            linestyle="none",
            color="tab:blue",
            marker="o",
        )
        if plot_both_models:
            plt.plot(
                r_in_array * pixscale,
                sb_model_im,
                color="tab:red",
                label="Exp Model ($\chi2$={:3.2f})".format(chi2),
            )

        plt.plot(
            r_in_array * pixscale,
            sb_moffat_im,
            label="Moffat Model ($\chi2$={:3.2f})".format(chi2_moffat),
            color="tab:orange",
        )
        plt.plot(
            np.arange(psf_xcen) * pixscale,
            np.max(sb)
            * (
                moffat_psf[0, psf_xcen:-1, psf_ycen] / moffat_psf[0, psf_xcen, psf_ycen]
            ),
            linestyle="dashed",
            color="tab:green",
            label="Extract Moffat PSF (FWHM:{:3.2f})".format(fwhm),
        )
        plt.yscale("log")
        plt.xlim(0, 15)
        # plt.ylim(10**-19,10**-14)
        plt.xlabel("Semi-major axis (arcsec)")
    else:
        plt.errorbar(
            r_in_array * pix_to_kpc,
            sb,
            yerr=sb_error,
            label="data",
            linestyle="none",
            color="tab:blue",
            marker="o",
        )

        plt.plot(
            r_in_array * pix_to_kpc,
            sb_model_im,
            label="Exp Model ($\chi2$={:3.2f})".format(chi2),
            color="r",
        )

        if plot_both_models:
            plt.plot(
                r_in_array * pix_to_kpc,
                sb_moffat_im,
                label="Moffat Model ($\chi2$={:3.2f})".format(chi2_moffat),
                color="tab:orange",
            )

        plt.plot(
            pix_to_kpc * np.arange(psf_xcen),
            np.max(sb)
            * (
                moffat_psf[0, psf_xcen:-1, psf_ycen] / moffat_psf[0, psf_xcen, psf_ycen]
            ),
            linestyle="dashed",
            color="tab:green",
            label="Extract Moffat PSF (FWHM:{:3.2f})".format(fwhm),
        )
        plt.xlim(0, 50)
        plt.xlabel("Semi-major axis (kpc)")
        plt.text(
            0.6,
            0.6,
            "r_n={:3.2f} kpc".format(r_n),
            transform=ax5.transAxes,
            fontsize=16,
        )
        secax = ax5.secondary_xaxis("top", functions=(kpc_to_arcsec, kpc_to_arcsec))
        secax.set_xlabel("Semi-major axis (arcsec)")
        plt.axvline(x=r_n, color="orange", linestyle="dashed")

    if plot_log:
        plt.yscale("log")
        plt.ylim(10 ** -19, 10 ** -14)

    plt.text(
        0.05,
        0.85,
        "{}".format(name),
        transform=ax3.transAxes,
        fontsize=14,
        color="black",
    )
    plt.ylabel(r"SB (erg/s/cm$^2$/arcsec$^2$")
    plt.legend(markerscale=1, fontsize=12, labelspacing=0.4)

    zscale = ZScaleInterval(contrast=0.5, krej=1.5)
    vmin, vmax = zscale.get_limits(values=image_data)

    ax6 = plt.subplot(gs[2, 0])
    plt.imshow(image_data, vmin=vmin, vmax=vmax)
    plt.text(
        0.5,
        0.9,
        r"Flux Ly$\alpha$",
        transform=ax6.transAxes,
        fontsize=16,
        ha="center",
        color="white",
    )
    ax7 = plt.subplot(gs[2, 1])

    if star:
        plt.imshow(moffat_im, vmin=vmin, vmax=vmax)
        plt.text(
            0.5,
            0.9,
            r"Moffat Model",
            transform=ax7.transAxes,
            fontsize=16,
            ha="center",
            color="white",
        )
    else:
        plt.imshow(model_im, vmin=vmin, vmax=vmax)
        plt.text(
            0.5,
            0.9,
            r"Exponential Model",
            transform=ax7.transAxes,
            fontsize=16,
            ha="center",
            color="white",
        )

    ax8 = plt.subplot(gs[2, 2])
    plt.text(
        0.5,
        0.9,
        r"Residual",
        transform=ax8.transAxes,
        fontsize=16,
        ha="center",
        color="white",
    )
    plt.imshow(image_data - model_im, vmin=vmin, vmax=vmax)
    plt.subplots_adjust(wspace=0.3,)
    if star:
        plt.savefig(
            "pyimfit_figs/star_{}.png".format(name), dpi=150, bbox_inches="tight"
        )
    else:
        plt.savefig("pyimfit_figs/{}.png".format(name), dpi=150, bbox_inches="tight")

    if star == False:
        np.save("pyimfit_output/params_{}".format(name), output)

    return


def get_parser():
    """ function that returns a parser from argparse """

    parser = ap.ArgumentParser(
        description="""Extracts 1D spectrum at specified RA/DEC""", add_help=True
    )

    parser.add_argument(
        "--wave_group_id", help="""wave_group_id""", type=int, default=None,
    )
    parser.add_argument(
        "--detectid", type=int, help="""detectid""", default=None,
    )
    parser.add_argument(
        "--coords",
        help="""Coordinates in RA DEC in degrees""",
        nargs=2,
        type=float,
        default=None,
        required=False,
    )
    parser.add_argument(
        "--wave_range",
        help="""Wave range to collapse image""",
        nargs=2,
        type=float,
        default=None,
        required=None,
    )
    parser.add_argument(
        "--shotid",
        help="""Shotid to use. This is required for coords/wave_range option.""",
        type=int,
        default=None,
        required=None,
    )
    parser.add_argument(
        "--pixscale", type=float, help="""pixel scale in arcsec""", default=0.25
    )
    parser.add_argument(
        "--imsize",
        type=float,
        help="""Image side (always a square) in arcsec""",
        default=30.0,
    )
    parser.add_argument(
        "--convolve_image",
        help="""Trigger to conolve image to FWHM seeing""",
        default=False,
        required=False,
        action="store_true",
    )

    subcont_parser = parser.add_mutually_exclusive_group(required=False)
    subcont_parser.add_argument('--subcont',
                                dest='subcont',
                                action='store_true')
    subcont_parser.add_argument('--no_subcont',
                                dest='subcont',
                                action='store_false')
    parser.set_defaults(subcont=True)
    
    parser.add_argument(
        "--star",
        help="""Trigger to fit a moffat profile only""",
        default=False,
        required=False,
        action="store_true",
    )
    return parser


def main(argv=None):
    """ Main Function """

    parser = get_parser()
    args = parser.parse_args(argv)

    print(args)
    # convert coords to a SkyCoord object if assigned
    if args.coords is not None:
        skycoord = SkyCoord(ra=args.coords[0] * u.deg, dec=args.coords[1] * u.deg)
    else:
        skycoord = None

    do_pyimfit(
        wave_group_id=args.wave_group_id,
        detectid=args.detectid,
        coords=skycoord,
        shotid=args.shotid,
        subcont=args.subcont,
        convolve_image=args.convolve_image,
        pixscale=args.pixscale,
        imsize=args.imsize,
        wave_range=args.wave_range,
        star=args.star,
    )


if __name__ == "__main__":
    main()
