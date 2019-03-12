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
from astropy import units as u
from scipy.interpolate import griddata
from shot import Fibers
import imp

setup_logging = imp.load_source('setup_logging', '/work/03946/hetdex/hdr1/software/HETDEX_API/input_utils.py')
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


def extract_source(xc, yc, xloc, yloc, data, mask, Dx, Dy,
                   scale=0.25, seeing_fac=1.8, fcor=1.,
                   boxsize=4.):
    ''' Extract a single source using rectification technique '''
    seeing = seeing_fac / scale
    a, b = data.shape
    xl = xc - boxsize / 2.
    xh = xc - boxsize / 2. + scale
    yl = yc - boxsize / 2.
    yh = yc - boxsize / 2. + scale
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
        if np.any(sel):
            grid_z = (griddata(S[sel], data[sel, k], (xgrid, ygrid),
                               method='cubic') * scale**2 / area)
            zgrid[k, :, :] = convolve(grid_z, G, boundary='extend')
    ml = [np.median(chunk, axis=0)
          for chunk in np.array_split(zgrid, 11, axis=0)]
    model = np.sum(ml, axis=0) / np.sum(ml) * fcor
    spec = (np.sum(zgrid * model[np.newaxis, :, :], axis=(1, 2)) /
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

def do_extraction(coord, fibers, ADRx, ADRy, radius=6.):
    ''' Grab fibers and do extraction '''
    boxsize = radius*np.sqrt(2.)/2.
    idx = fibers.query_region_idx(coord, radius=radius/3600.)
    if len(idx) < fiber_lower_limit:
        return None
    ifux = fibers.table.read_coordinates(idx, 'ifux')
    ifuy = fibers.table.read_coordinates(idx, 'ifuy')
    ra = fibers.table.read_coordinates(idx, 'ra')
    dec = fibers.table.read_coordinates(idx, 'dec')
    spec = fibers.table.read_coordinates(idx, 'calfib')
    ftf = fibers.table.read_coordinates(idx, 'fiber_to_fiber')
    mask = fibers.table.read_coordinates(idx, 'Amp2Amp')
    mask = (mask > 1e-8) * (ftf > 0.5)
    expn = fibers.table.read_coordinates(idx, 'expnum')
    ifux, ifuy, xc, yc = get_new_ifux_ifuy(expn, ifux, ifuy, ra, dec,
                                           coord.ra.deg, coord.dec.deg)
    exspec, zgrid, model, xg, yg = extract_source(xc, yc, ifux, ifuy, spec,
                                                  mask, ADRx, ADRy,
                                                  boxsize=boxsize)
    return exspec, zgrid, model, xg, yg

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
log = setup_logging('extraction')

log.info('Getting HDF5 file')
fibers = Fibers(shotv)

fiber_lower_limit = 3
log.info('Getting stars in astrometry catalog')
ras = fibers.hdfile.root.Astrometry.StarCatalog.cols.ra_cat[:]
decs = fibers.hdfile.root.Astrometry.StarCatalog.cols.dec_cat[:]
coords = SkyCoord(ras*u.deg, decs*u.deg, frame='icrs')
log.info('Number of stars to extract: %i' % len(coords))

for i, coord in enumerate(coords):
    log.info('Extracting %ith coordinate' % (i+1))
    spectrum, cube, weights, xg, yg = do_extraction(coord, fibers, ADRx, ADRy)
    log.info('Making cube for %ith coordinate' % (i+1))
    write_cube(wave, xg, yg, cube, 'test_cube_%i.fits' % (i+1))
