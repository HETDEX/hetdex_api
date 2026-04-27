import sys
import numpy as np
import os
import shutil

import os.path as op
import tables as tb

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
from astropy.table import Table, unique, vstack
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from photutils.aperture import aperture_photometry, SkyEllipticalAperture, SkyEllipticalAnnulus

from hetdex_api.config import HDRconfig
from hetdex_tools.phot_tools import get_line_image
from hetdex_api.extract import Extract

from elixer import catalogs

def reset_folder(shotid):
    parent = "im_figs"
    folder = os.path.join(parent, str(shotid))

    # Make sure parent exists
    os.makedirs(parent, exist_ok=True)

    # If subfolder exists, remove and recreate
    if os.path.exists(folder):
        shutil.rmtree(folder)
    os.makedirs(folder)

    return folder

waveoii = 3727.8

imsize_arcsec = 20
imsize = imsize_arcsec / 3600.0

catlib = catalogs.CatalogLibrary()
# set up catalog files to save memory
config = HDRconfig('hdr5')

global elixh5

elixh5 = tb.open_file(config.elixerh5, "r")

def measure_app_flux(row, E=None):
    global elixh5
    
    detectid_obj = row['detectid']
    z_hetdex = row["z_hetdex"]
    shotid_obj = row["shotid"]
    flux = row["flux"]
    flux_err = row["flux_err"]

    wave_1 = (waveoii - 15) * (1 + z_hetdex)
    wave_2 = (waveoii + 15) * (1 + z_hetdex)

    obj = Table(elixh5.root.ExtractedObjects.read_where("detectid == detectid_obj"))
    obj.sort('dist_baryctr')

    selected = obj["selected"] == True

    if np.sum(selected) == 0:
        return None
    elif np.sum(selected) == 1:
        pass
    elif np.sum(selected) >= 2:
        filter_order = [b"r", b"f606w", b"g"]

        for filt in filter_order:
            selected = (obj["selected"] == True) & (obj["filter_name"] == filt)

            if np.sum(selected) == 1:
                break
            elif np.sum(selected) >= 2:
                selected = (
                    (obj["selected"] == True)
                    & (obj["filter_name"] == filt)
                    & ((obj["catalog_name"] == b"HSC-DEX") | (obj["catalog_name"] == b"HSC-SSP"))
                )
                if np.sum(selected) == 1:
                    break
    else:
        filter_order = [b"r", b"f606w", b"g"]

        for filt in filter_order:
            selected = obj["filter_name"] == filt

            if np.sum(selected) == 1:
                break
            elif np.sum(selected) >= 2:
                selected = (obj["filter_name"] == filt) & (
                    (obj["catalog_name"] == b"HSC-DEX") | (obj["catalog_name"] == b'HSC-SSP')
                )
                if np.sum(selected) == 1:
                    break
                else:
                    selected = obj["filter_name"] == filt
                    if np.sum(selected) >= 1:
                        break
            else:
                selected = obj["filter_name"] == filt
                if np.sum(selected) >= 1:
                    break
                
    ra = obj["ra"][selected][0]
    dec = obj["dec"][selected][0]

    major = obj["major"][selected][0] * u.arcsec
    minor = obj["minor"][selected][0] * u.arcsec
    theta = (np.pi / 2 + obj["theta"][selected][0]) * u.rad

    coords = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)

    aper = SkyEllipticalAperture(coords, major, minor, theta)

    cutout = catlib.get_cutouts(
        position=coords,
        side=imsize,
        filter=["f606W", "r", "g"],
        first=True,
        allow_bad_image=False,
        allow_web=True,
    )[0]

    im = cutout["cutout"].data
    wcs = cutout["cutout"].wcs
    plt.figure(figsize=(8, 8))

    plt.subplot(projection=wcs)
    ax = plt.gca()
    ax.coords[0].set_major_formatter("d.dddd")
    # ax.coords[0].set_ticks(spacing=1. * u.arcsec)
    ax.coords[1].set_major_formatter("d.dddd")
    # ax.coords[0].set_ticks(spacing=1. * u.arcsec)

    vmin = 3
    vmax = 99
    impix = im.shape[0]
    pixscale = imsize / impix  # in degrees/pix
    m = np.percentile(im, (vmin, vmax))
    plt.imshow(im, vmin=m[0], vmax=m[1], origin="lower", cmap="gray_r")
    plt.text(
        0.95,
        0.05,
        cutout["instrument"] + cutout["filter"],
        transform=ax.transAxes,
        fontsize=20,
        color="red",
        horizontalalignment="right",
    )

    plt.text(
        0.95,
        0.95,
        str(detectid_obj),
        transform=ax.transAxes,
        fontsize=20,
        color="red",
        horizontalalignment="right",
    )

    aper.to_pixel(wcs).plot(color="r", linewidth=2)

    hdu = get_line_image(
        coords=coords,
        imsize=imsize_arcsec,
        wave_range=[wave_1, wave_2],
        shotid=shotid_obj,
        extract_class=E, #pass in Extract class to limit I/O
    )
    w = WCS(hdu['DATA'].header)
    plt.contour(hdu['DATA'].data, transform=ax.get_transform(w))

    plt.xlabel("RA")
    plt.ylabel("Dec")
    # create mask (set True for where you want to mask)
    im_mask = hdu['DATA'].data == 0

    phottable = aperture_photometry(
        hdu['DATA'].data, aper, error=hdu['ERROR'].data, mask=im_mask, wcs=w,
    )
    im_flux = phottable["aperture_sum"][0]  # * u.Unit("10^-17 erg cm-2")
    im_flux_err = phottable["aperture_sum_err"][0]  # * u.Unit("10^-17 erg cm-2")

    apcor_im = aper.to_pixel(w).to_mask(method="center").multiply(
        im_mask
    ) / aper.to_pixel(w).to_mask(method="center").multiply(np.ones_like(im))

    apcor = np.sum(apcor_im == 0) / np.sum(np.isfinite(apcor_im))

    aper_annulus = SkyEllipticalAnnulus(coords,
                                        1.2*major,
                                        1.5*major,
                                        1.5*minor,
                                        b_in=1.2*minor,
                                        theta=theta)

    mask = aper_annulus.to_pixel(w).to_mask(method="center").data
    annulus_data = aper_annulus.to_pixel(w).to_mask(method="center").multiply(hdu['DATA'].data)
    annulus_data_1d = annulus_data[mask > 0]
    annulus_mask = aper_annulus.to_pixel(w).to_mask(method="center").multiply(im_mask)
    annulus_mask_1d = annulus_mask[mask > 0]

    # get median and standard deviation in annulus
    mean_sigclip, median_sigclip, stddev_sigclip = sigma_clipped_stats(
            annulus_data_1d, mask=annulus_mask_1d
    )
    bkg_median = median_sigclip * aper.to_pixel(w).area * apcor
    bkg_stddev = stddev_sigclip * aper.to_pixel(w).area * apcor
    
    plt.title("Aperture Line Flux: {:4.2f} +/- {:4.2f} 10^-17 erg/cm^2 \n".format( im_flux/apcor, im_flux_err/apcor) + "HDR5 Line Flux: {:4.2f} +/- {:4.2f} 10^-17 erg/cm^2".format(flux/apcor, flux_err/apcor))
    aper_annulus.to_pixel(wcs).plot(color="red", linestyle="dashed")

    im_fig_name = 'im_figs/{}/im_flux_{}.png'.format(shotid_obj, detectid_obj)
    plt.savefig(im_fig_name, bbox_inches='tight', dpi=150)

    data = Table( obj[selected][0])
    data["im_flux"] = im_flux/apcor
    data["im_flux_err"] = im_flux_err/apcor
    data["bkg_median"] = bkg_median/apcor
    data["bkg_stddev"] = bkg_stddev/apcor
    data["im_apcor"] = apcor
    
    return data

def main(argv=None):

    """ Main Function """
    shotid = int(sys.argv[1])

    # exit if output folder is already generated

    if op.exists( "output/im_info_{}.txt".format(shotid)):
        print(f'Output file for {shotid} already exists. Exiting.')
        sys.exit()
    # create images directory, clear it if it already exists
    
    folder_path = reset_folder(shotid)    
    print(f"Images Folder Ready at: {folder_path}")

    E = Extract()
    E.load_shot(shotid)
    
    output = Table()

    dets_table = Table.read('detlists/{}.txt'.format(shotid), format='ascii')

    for det_row in dets_table:
        try:
            data = measure_app_flux(det_row, E=E)
            
            if data is not None:
                output = vstack([output, data])
        except Exception:
            print("Could not measure aperture flux for {}".format(det_row['detectid']))

    output.write("output/im_info_{}.txt".format(shotid), format="ascii", overwrite=True)
    E.close()
    elixh5.close()
if __name__ == "__main__":
    main()
        
