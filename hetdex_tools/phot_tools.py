#Erin's phot tools
#Not for general use ATM.. talk to Erin

import numpy as np
import tables as tb
import os.path as op
import time

from photutils.isophote import EllipseGeometry, Ellipse, IsophoteList
from photutils import EllipticalAperture, SkyEllipticalAperture
from photutils import SkyCircularAperture
from photutils import aperture_photometry
from photutils.isophote import build_ellipse_model

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import wcs
from astropy.io import fits
from astropy.table import Table, unique
from astropy.modeling import models, fitting
from astropy.cosmology import WMAP9 as cosmo
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.convolution import Moffat2DKernel
from astropy.convolution import convolve_models
from astropy.convolution import Gaussian1DKernel
from astropy.visualization import ZScaleInterval

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from hetdex_api.config import HDRconfig
from hetdex_api.detections import Detections
from hetdex_tools.interpolate import make_narrowband_image

from elixer import catalogs


plt.style.use('fivethirtyeight') 


LATEST_HDR_NAME = HDRconfig.LATEST_HDR_NAME

config = HDRconfig()
surveyh5 = tb.open_file(config.surveyh5, "r")
deth5 = tb.open_file(config.detecth5, "r")
conth5 = tb.open_file(config.contsourceh5, "r")

def FitCircularAperture(hdu=None, coords=None, radius=3*u.arcsec, plot=False):
    """
    Fit a circular aperture with either an HDU object or and
    """
    
    im = hdu[0].data
    w = wcs.WCS(hdu[0].header)
    aper = SkyCircularAperture(coords, r=radius)
    if plot:
        plt.subplot(111, projection=w)
        plt.imshow(im)
        aper.to_pixel(w).plot(color='white') # for SkyCircularAperture
        plt.xlabel('RA')
        plt.ylabel('Dec')
        
    phottable = aperture_photometry(hdu[0].data, aper,
                                    error=hdu[1].data,
                                    wcs=wcs)
    flux = phottable['aperture_sum'][0] * u.Unit('10^-17 erg cm-2 s-1')
    flux_err = phottable['aperture_sum_err'][0] * u.Unit('10^-17 erg cm-2 s-1')
        
    return flux, flux_err


def get_flux_for_detectid(detectid):
    detectid_obj = detectid
    if detectid_obj <= 2190000000:
        det_info = deth5.root.Detections.read_where('detectid == detectid_obj')[0]
    else:
        det_info = conth5.root.Detections.read_where('detectid == detectid_obj')[0]
        
    coords = SkyCoord(det_info['ra'], det_info['dec'], unit='deg')
    wave_obj = det_info['wave']
    #redshift = wave_obj/(1216) - 1
    linewidth = det_info['linewidth']
    
    shotid_obj = det_info['shotid']
    fwhm = surveyh5.root.Survey.read_where('shotid == shotid_obj')['fwhm_virus'][0]
    amp = det_info['multiframe']

    try:
        hdu = make_narrowband_image(detectid = detectid_obj,
                                    imsize=10*u.arcsec,
                                    pixscale=0.25*u.arcsec,
                                    convolve_image=True,
                                    subcont=True, include_error=True)
    except:
        print('Could not make narrowband image for {}'.format(detectid_obj))
        return np.nan, np.nan
    flux, flux_err = FitCircularAperture(hdu=hdu, coords=coords, plot=False)
    return flux, flux_err



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
        
    plt.imshow(cutout["cutout"].data,
               vmin=im_vmin,
               vmax=im_vmax,
               origin="lower",
               cmap=plt.get_cmap("gray"),
               interpolation="none")
    
    sel = friend_cat['friendid'] == friendid
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
    
    #plot detections for group
    plt.scatter(group['ra'],
                group['dec'],
                transform=ax.get_transform('fk5'),
                marker='x',
                color='orange',
                linewidth=1,
                s=group['flux'],
    )
    
    # plot and elliptical kron-like aperture representing the group
    # theta is measured
    aper_group = SkyEllipticalAperture(coords,
                                       a*u.deg,
                                       b*u.deg,
                                       theta=(90-pa)*u.deg)
    aper_group_k = SkyEllipticalAperture(coords,
                                         k*a*u.deg,
                                         k*b*u.deg,
                                         theta=(90-pa)*u.deg)
    aper_group.to_pixel(wcs).plot(color='blue')
    aper_group_k.to_pixel(wcs).plot(color='red')
    
    if label:
        # plot detecid labels
        for row in group:
            plt.text(
                row['ra'],
                row['dec'],
                str(row['detectid']),
                transform=ax.get_transform("world"),
                fontsize=9,
                color="red",
            )
def fit_ellipse_for_id(friendid=None,
                       detectid=None,
                       subcont=True,
                       convolve_image=True,
                       pixscale=pixscale,
                       imsize=imsize,
):
    global deth5
    
    if detectid is not None:
        
        detectid_obj=detectid
        
        if detectid_obj <= 2190000000:
            det_info = deth5.root.Detections.read_where('detectid == detectid_obj')[0]
        else:
            det_info = conth5.root.Detections.read_where('detectid == detectid_obj')[0]
            
        coords = SkyCoord(det_info['ra'], det_info['dec'], unit='deg')
        wave_obj = det_info['wave']
        redshift = wave_obj/(1216) - 1
        linewidth = det_info['linewidth']
        
        shotid_obj = det_info['shotid']
        fwhm = surveyh5.root.Survey.read_where('shotid == shotid_obj')['fwhm_virus'][0]
        amp = det_info['multiframe']
        
        try:
            hdu = make_narrowband_image(detectid = detectid_obj,
                                        imsize=imsize*u.arcsec,
                                        pixscale=pixscale*u.arcsec,
                                        convolve_image=convolve_image,
                                        subcont=subcont, include_error=True)
        except:
            print('Could not make narrowband image for {}'.format(detectid))
            return np.nan, np.nan
            
    elif friendid is not None:
        global friend_cat
        
        sel = friend_cat['friendid'] == friendid
        group = friend_cat[sel]
        coords = SkyCoord(ra=group["icx"][0] * u.deg, dec=group["icy"][0] * u.deg)
        wave_obj = group['icz'][0]
        redshift = wave_obj/(1216) - 1
        linewidth = group['linewidth'][0]
        shotid_obj = group['shotid'][0]
        fwhm = group['fwhm'][0]
        amp = group['multiframe'][0]
        
        try:
            hdu = make_narrowband_image(coords=coords,
                                        shotid=shotid_obj,
                                        wave_range=[wave_obj-4,
                                                    wave_obj+4],
                                        imsize=imsize*u.arcsec,
                                        pixscale=pixscale*u.arcsec,
                                        subcont=subcont,
                                        convolve_image=convolve_image,
                                        include_error=True)
        except:
            print('Could not make narrowband image for {}'.format(friendid))
            return np.nan, np.nan
        else:
            print('You must provide a detectid or friendid')
            
            
        w= wcs.WCS(hdu[0].header)
        
        if friendid is not None:
            
            sel_friend_group = friend_cat['friendid'] == friendid
            group = friend_cat[sel_friend_group]
            eps = 1-group['a2'][0]/group['b2'][0]
            pa = group['pa'][0]*np.pi/180.-90
            sma = group['a'][0]*3600/pixscale
            
            coords = SkyCoord(ra=group["icx"][0] * u.deg, dec=group["icy"][0] * u.deg)
            wave_obj = group['icz'][0]
            redshift = wave_obj/(1216) - 1
            linewidth = np.nanmedian( group['linewidth'] )
            shotid_obj = group['shotid'][0]
            fwhm = group['fwhm'][0]
            
            geometry = EllipseGeometry(x0=w.wcs.crpix[0], y0=w.wcs.crpix[0],
                                       sma=sma, eps=eps, pa=pa)
        else:
            geometry = EllipseGeometry(x0=w.wcs.crpix[0],
                                       y0=w.wcs.crpix[0],
                                       sma=20, eps=0.2, pa=20.0)
            
        geometry = EllipseGeometry(x0=w.wcs.crpix[0],
                                   y0=w.wcs.crpix[0],
                                   sma=20, eps=0.2, pa=20.0)
        #geometry.find_center(hdu.data)
        #aper = EllipticalAperture((geometry.x0, geometry.y0), geometry.sma,
        #                          geometry.sma*(1 - geometry.eps), geometry.pa)
        
        #plt.imshow(hdu.data, origin='lower')
        #aper.plot(color='white')
        
        ellipse = Ellipse(hdu[0].data)
        isolist = ellipse.fit_image()
        iso_tab = isolist.to_table()
        #print(iso_tab)
        
        if len(iso_tab) == 0:
            geometry.find_center(hdu[0].data, verbose=False, threshold=0.5)
            ellipse = Ellipse(hdu[0].data, geometry)
            isolist = ellipse.fit_image()
            iso_tab = isolist.to_table()
            #print(iso_tab)
            
        if len(iso_tab) == 0:
            #exit program
            return np.nan, np.nan
            
            try:
                #compute iso's manually in steps of 3 pixels
                ellipse = Ellipse(hdu[0].data)#reset ellipse
                iso_list = []
                for sma in np.arange(1,60,2):
                    iso = ellipse.fit_isophote(sma)
                    if np.isnan(iso.intens):
                        #print('break at {}'.format(sma))
                        break
                    else:
                        iso_list.append(iso)
                        isolist = IsophoteList(iso_list)
                        iso_tab = isolist.to_table()
            except:
                return np.nan, np.nan
                
            try:
                model_image = build_ellipse_model(hdu[0].data.shape, isolist)
                residual = hdu[0].data - model_image
            except:
                return np.nan, np.nan
            
        sma = iso_tab['sma']* pixscale
        const_arcsec_to_kpc = cosmo.kpc_proper_per_arcmin(redshift).value/60.
    
        def arcsec_to_kpc(sma):
            dist = const_arcsec_to_kpc*sma
            return dist
            
        def kpc_to_arcsec(dist):
            sma = dist/const_arcsec_to_kpc
            return sma
            
        dist_kpc = sma*u.arcsec.to(u.arcmin)*u.arcmin * cosmo.kpc_proper_per_arcmin(redshift)
        dist_arcsec = kpc_to_arcsec(dist_kpc)
        
