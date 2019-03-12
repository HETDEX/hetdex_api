#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 11:48:55 2019

@author: gregz
"""

import numpy as np

from astropy.convolution import Gaussian2DKernel, convolve
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.modeling.models import Gaussian2D
from astropy.modeling.fitting import LevMarLSQFitter
from astropy import units as u
from scipy.interpolate import griddata
from shot import Fibers
import imp

input_utils = imp.load_source('input_utils', '/work/03946/hetdex/hdr1/software/HETDEX_API/input_utils.py')
#
#import warnings
## astroquery emits warning for importing SDSS, ignore that
#with warnings.catch_warnings():
#    warnings.simplefilter("ignore")
#    from astroquery.sdss import SDSS
#
#def querySDSS(ra, dec, radius):
#    ''' Using astroquery sdss system '''
#    pos = SkyCoord(ra*u.deg, dec*u.deg, frame='fk5')
#    table = SDSS.query_region(pos, radius=radius*u.deg,
#                              photoobj_fields=['ra', 'dec', 'objid', 'type',
#                                               'u', 'g', 'r', 'i', 'z'])
#    return table

def get_wave():
    ''' Return implicit wavelength solution for calfib extension '''
    return np.linspace(3470, 5540, 1036)


def get_ADR(wave, angle=0.):
    ''' Use default ADR with angle = 0 along x-direction '''
    wADR = [3500., 4000., 4500., 5000., 5500.]
    ADR = [-0.74, -0.4, -0.08, 0.08, 0.20]
    ADR = np.polyval(np.polyfit(wADR, ADR, 3), wave)
    ADRx = np.cos(np.deg2rad(angle)) * ADR
    ADRy = np.sin(np.deg2rad(angle)) * ADR
    return ADRx, ADRy


def model_source(data, mask, xloc, yloc, wave, chunks=11):
    G = Gaussian2D()
    fitter = LevMarLSQFitter()
    wchunk = np.array([np.mean(chunk)
                       for chunk in np.array_split(wave, chunks)])
    xc, yc, xs, ys = [i * wchunk for i in [0., 0., 0., 0.]]

    A = np.zeros((chunks, len(xloc)))
    B = np.zeros((chunks, len(xloc)))
    i = 0
    for chunk, maskchunk in zip(np.array_split(data, chunks, axis=1),
                                np.array_split(mask, chunks, axis=1)):
        c = np.ma.array(chunk, mask=maskchunk==0.)
        image = np.ma.median(c, axis=1)
        y = image.data
        ind = np.ma.argmax(image)
        dist = np.sqrt((xloc - xloc[ind])**2 + (yloc - yloc[ind])**2)
        inds = (dist < 3.) * (~image.mask)
        x_centroid = np.sum(y[inds] * xloc[inds]) / np.sum(y[inds])
        y_centroid = np.sum(y[inds] * yloc[inds]) / np.sum(y[inds])
        G.amplitude.value = y[ind]
        G.x_mean.value = x_centroid
        G.y_mean.value = y_centroid
        fit = fitter(G, xloc[inds], yloc[inds], y[inds])
        xc[i] = fit.x_mean.value * 1.
        yc[i] = fit.y_mean.value * 1.
        xs[i] = fit.x_stddev.value * 1.
        ys[i] = fit.y_stddev.value * 1.
        A[i] = image.data / np.ma.sum(image)
        B[i] = image.mask
        i += 1
    return A, B, xc, yc, xs, ys, wchunk

    


def extract_source(xc, yc, xloc, yloc, data, mask, Dx, Dy,
                   scale=0.25, seeing_fac=1.8, fcor=1.,
                   boxsize=4.):
    ''' Extract a single source using rectification technique '''
    seeing = seeing_fac / scale
    a, b = data.shape
    xl = xc - boxsize / 2.
    xh = xc + boxsize / 2. + scale
    yl = yc - boxsize / 2.
    yh = yc + boxsize / 2. + scale
    x = np.arange(xl, xh, scale)
    y = np.arange(yl, yh, scale)
    xgrid, ygrid = np.meshgrid(x, y)
    zgrid = np.zeros((b,)+xgrid.shape)
    area = np.pi * 0.75**2
    G = Gaussian2DKernel(seeing / 2.35)
    S = np.zeros((data.shape[0], 2))
    
    for k in np.arange(b):
        S[:, 0] = xloc - Dx[k]
        S[:, 1] = yloc - Dy[k]
        sel = np.isfinite(data[:, k]) * (mask[:, k] != 0.0)
        if sel.sum()>4:
            grid_z = (griddata(S[sel], data[sel, k], (xgrid, ygrid),
                               method='cubic') * scale**2 / area)
            zgrid[k, :, :] = convolve(grid_z, G, boundary='extend',
                                      nan_treatment='interpolate')
    ml = [np.nanmedian(chunk, axis=0)
          for chunk in np.array_split(zgrid, 11, axis=0)]
    model = np.nansum(ml, axis=0) / np.nansum(ml) * fcor
    model[~np.isfinite(model)] = 0.
    spec = (np.nansum(zgrid * model[np.newaxis, :, :], axis=(1, 2)) /
            np.sum(model**2))
    return spec, zgrid, model, xgrid, ygrid


def get_new_ifux_ifuy(expn, ifux, ifuy, ra, dec, rac, decc):
    ''' Get ifux and ifuy from RA (correct for dither pattern) '''
    s = np.where(expn == 1.)[0]
    if len(s) < 2.:
        return None
    ifu_vect = np.array([ifuy[s[1]] - ifuy[s[0]], ifux[s[1]] - ifux[s[0]]])
    radec_vect = np.array([(ra[s[1]] - ra[s[0]]) * np.cos(np.deg2rad(dec[s[0]])),
                          dec[s[1]] - dec[s[0]]])
    V = np.sqrt(ifu_vect[0]**2 + ifu_vect[1]**2)
    W = np.sqrt(radec_vect[0]**2 + radec_vect[1]**2)
    scale_vect = np.array([3600., 3600.])
    v = ifu_vect * np.array([1., 1.])
    w = radec_vect * scale_vect
    W =  np.sqrt(np.sum(w**2))
    ang1 = np.arctan2(v[1] / V, v[0] / V)
    ang2 = np.arctan2(w[1] / W, w[0] / W)
    ang1 += (ang1 < 0.) * 2. * np.pi
    ang2 += (ang2 < 0.) * 2. * np.pi
    theta = ang1 - ang2
    dra = (ra - ra[s[0]]) * np.cos(np.deg2rad(dec[s[0]])) * 3600.
    ddec = (dec - dec[s[0]]) * 3600.
    dx = np.cos(theta) * dra - np.sin(theta) * ddec
    dy = np.sin(theta) * dra + np.cos(theta) * ddec
    yy = dx + ifuy[s[0]]
    xx = dy + ifux[s[0]]
    dra = (rac - ra[s[0]]) * np.cos(np.deg2rad(dec[s[0]])) * 3600.
    ddec = (decc - dec[s[0]]) * 3600.
    dx = np.cos(theta) * dra - np.sin(theta) * ddec
    dy = np.sin(theta) * dra + np.cos(theta) * ddec
    yc = dx + ifuy[s[0]]
    xc = dy + ifux[s[0]]
    return xx, yy, xc, yc

def do_extraction(coord, fibers, ADRx, ADRy, radius=8.):
    ''' Grab fibers and do extraction '''
    boxsize = radius*np.sqrt(2.) / 2. * 0.9
    idx = fibers.query_region_idx(coord, radius=radius/3600.)
    fiber_lower_limit = 5
    if len(idx) < fiber_lower_limit:
        return None
    ifux = fibers.table.read_coordinates(idx, 'ifux')
    ifuy = fibers.table.read_coordinates(idx, 'ifuy')
    ra = fibers.table.read_coordinates(idx, 'ra')
    dec = fibers.table.read_coordinates(idx, 'dec')
    spec = fibers.table.read_coordinates(idx, 'calfib')
    spece = fibers.table.read_coordinates(idx, 'calfibe')
    ftf = fibers.table.read_coordinates(idx, 'fiber_to_fiber')
    mask = fibers.table.read_coordinates(idx, 'Amp2Amp')
    mask = (mask > 1e-8) * (np.median(ftf, axis=1) > 0.5)[:, np.newaxis]
    expn = fibers.table.read_coordinates(idx, 'expnum')
    ifux, ifuy, xc, yc = get_new_ifux_ifuy(expn, ifux, ifuy, ra, dec,
                                           coord.ra.deg, coord.dec.deg)
    #exspec, zgrid, model, xg, yg = extract_source(xc, yc, ifux, ifuy, spec,
    #                                              mask, ADRx, ADRy,
    #                                              boxsize=boxsize)
    #return exspec, zgrid, model, xg, yg
    result = model_source(spec, mask, ifux, ifuy, wave, chunks=11)
    return result

def write_cube(wave, xgrid, ygrid, zgrid, outname):
    hdu = fits.PrimaryHDU(np.array(zgrid, dtype='float32'))
    hdu.header['CRVAL1'] = xgrid[0, 0]
    hdu.header['CRVAL2'] = ygrid[0, 0]
    hdu.header['CRVAL3'] = wave[0]
    hdu.header['CRPIX1'] = 1
    hdu.header['CRPIX2'] = 1
    hdu.header['CRPIX3'] = 1
    hdu.header['CTYPE1'] = 'pixel'
    hdu.header['CTYPE2'] = 'pixel'
    hdu.header['CTYPE3'] = 'pixel'
    hdu.header['CDELT1'] = xgrid[0, 1] - xgrid[0, 0]
    hdu.header['CDELT2'] = ygrid[1, 0] - ygrid[0, 0]
    hdu.header['CDELT3'] = wave[1] - wave[0]
    hdu.writeto(outname, overwrite=True)

wave = get_wave()
ADRx, ADRy = get_ADR(wave)
shotv = '20190208v024'
log = input_utils.setup_logging('extraction')

log.info('Getting HDF5 file')
fibers = Fibers(shotv)

log.info('Getting stars in astrometry catalog')
ras = fibers.hdfile.root.Astrometry.StarCatalog.cols.ra_cat[:]
decs = fibers.hdfile.root.Astrometry.StarCatalog.cols.dec_cat[:]
gmag = fibers.hdfile.root.Astrometry.StarCatalog.cols.g[:]
starid = fibers.hdfile.root.Astrometry.StarCatalog.cols.star_ID[:]


coords = SkyCoord(ras*u.deg, decs*u.deg, frame='fk5')
log.info('Number of stars to extract: %i' % len(coords))
L, M = ([], [])
for i, coord in enumerate(coords):
    log.info("Star %i, g' magnitude: %0.2f" % (starid[i], gmag[i]))
    log.info('Extracting coordinate #%i' % (i+1))
    result = do_extraction(coord, fibers, ADRx, ADRy)
    if result is not None:
        P = []
        for j in np.arange(len(result)):
            if j == 0:
                f = fits.PrimaryHDU
            else:
                f = fits.ImageHDU
            P.append(f(result[j]))
        fits.HDUList(P).writeto('test_model_%i.fits' %(i+1))
#    if result is not None:
#        spectrum, cube, weights, xg, yg = result
#        log.info('Making cube for coordinate #%i' % (i+1))
#        write_cube(wave, xg, yg, cube, 'test_cube_%i.fits' % (i+1))
#        L.append(spectrum)
#        M.append(weights)
#fits.PrimaryHDU(np.vstack(L)).writeto('allspec.fits', overwrite=True)
#fits.PrimaryHDU(np.hstack(M)).writeto('allweights.fits', overwrite=True)

