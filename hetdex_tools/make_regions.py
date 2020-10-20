# -*- coding: utf-8 -*-
"""
Created on September 23 2020
@author: Erin Mentuch Cooper
"""
import sys

from astropy.coordinates import SkyCoord
from astropy import units as u

from regions import RectangleSkyRegion, write_ds9

from hetdex_api.flux_limits.hdf5_sensitivity_cubes \
    import (SensitivityCubeHDF5Container, return_sensitivity_hdf_path)


date = sys.argv[1]
obs = sys.argv[2]

datevobs = str(date) + 'v' + str(obs).zfill(3)
shotid = int(str(date) + str(obs).zfill(3))

hdf_filename, mask_fn = return_sensitivity_hdf_path(datevobs, release="hdr2.1",
                                                    return_mask_fn=True)

hdfcont = SensitivityCubeHDF5Container(filename=hdf_filename,
                                       aper_corr=1.0,
                                       flim_model="hdr2pt1",
                                       mask_filename=mask_fn)

ifu_name_list = []
ifu_ra = []
ifu_dec = []
sncut = 6
ifuregions = []

for ifu_name, tscube in hdfcont.itercubes():
    slice_ = tscube.f50_from_noise(tscube.sigmas[200, :, :], sncut)
    shape = tscube.sigmas.shape
    ra, dec, lambda_ = tscube.wcs.all_pix2world(shape[2]/2.,
                                                shape[1]/2.,
                                                shape[0]/2.,
                                                0)
    pa = -tscube.header['CROTA2']
    ifu_name_list.append(ifu_name)
    coord = SkyCoord(ra, dec, unit='deg')
    ifu_ra.append(ra)
    ifu_dec.append(dec)
    ifuregions.append( RectangleSkyRegion(center=coord, width=1.*u.arcmin, height=1.*u.arcmin, angle=pa*u.deg))

region_file = datevobs + '.reg'

write_ds9(ifuregions, region_file)
