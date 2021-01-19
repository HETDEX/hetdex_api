# -*- coding: utf-8 -*-
"""
Created on September 23 2020
@author: Erin Mentuch Cooper

"""
import sys
import numpy as np

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

from astropy.wcs import WCS
from astropy.visualization import ZScaleInterval
from astropy.visualization.wcsaxes import WCSAxes
from astropy.coordinates import SkyCoord
from astropy.io import fits

from reproject import reproject_interp, reproject_exact
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject.mosaicking import reproject_and_coadd

from hetdex_api.survey import Survey
from hetdex_api.config import HDRconfig
from hetdex_api.mask import *
from hetdex_api.flux_limits.hdf5_sensitivity_cubes import (SensitivityCubeHDF5Container, return_sensitivity_hdf_path)
from hetdex_api.detections import Detections

date = sys.argv[1]
obs = sys.argv[2]

try:
    wave_slice = float(sys.argv[3])
except:
    wave_slice = 4540.
    
datevobs = str(date) + 'v' + str(obs).zfill(3)
shotid = int(str(date) + str(obs).zfill(3))

survey = Survey()

shot_coords = survey.coords[survey.shotid == shotid][0]
hdf_filename_hdr2pt1, mask_fn = return_sensitivity_hdf_path(datevobs, release="hdr2.1", return_mask_fn=True)

hdfcont_hdr2 = SensitivityCubeHDF5Container(filename=hdf_filename_hdr2pt1, aper_corr=1.0,
                                            flim_model="hdr2pt1",  mask_filename=mask_fn)

#overplot detections
detects = Detections(curated_version='2.1.2')
sel_shot = detects.shotid == shotid
detects_shot = detects[sel_shot]

# sel slice to plot.. can update this later as a variable
wave = 3470 + np.arange(0 , 1036)*2.
sel_slice = np.where(wave >= wave_slice)[0][0]

hdus = []
hdus_mask = []

ifu_name_list = []
ifu_ra = []
ifu_dec = []

sncut=6

for ifu_name, tscube in hdfcont_hdr2.itercubes():
    slice_ = tscube.f50_from_noise(tscube.sigmas[sel_slice, :, :], sncut)
    hdus.append( fits.PrimaryHDU( slice_*1e17, header=tscube.wcs.celestial.to_header()))
    hdus_mask.append( fits.PrimaryHDU( slice_.mask.astype(int), header=tscube.wcs.celestial.to_header()))

    shape = tscube.sigmas.shape
    ra, dec, lambda_ = tscube.wcs.all_pix2world(shape[2]/2., shape[1]/2., shape[0]/2., 0)
    ifu_name_list.append(ifu_name)
    ifu_ra.append(ra)
    ifu_dec.append(dec)

wcs_out, shape_out = find_optimal_celestial_wcs(hdus, reference=shot_coords)
wcs_mask_out, shape_mask_out = find_optimal_celestial_wcs(hdus_mask, reference=shot_coords)

array, footprint = reproject_and_coadd(hdus,
                                       wcs_out,
                                       shape_out=shape_out,
                                       reproject_function=reproject_exact)

#mask_array, footprint = reproject_and_coadd(hdus_mask,
#                                            wcs_mask_out,
#                                            shape_out=shape_mask_out,
#                                            reproject_function=reproject_exact)

config = HDRconfig()
galaxy_cat = Table.read(config.rc3cat, format='ascii')
gal_coords = SkyCoord(galaxy_cat['Coords'], frame='icrs')
sel_reg = np.where(shot_coords.separation(gal_coords) < 1.*u.deg)[0]

gal_regions = []
for idx in sel_reg:
    gal_regions.append( create_gal_ellipse(galaxy_cat, row_index=idx, d25scale=3))

plt.figure(figsize=(15,12))
plt.rcParams.update({'font.size': 22})
ax = plt.subplot(111, projection=wcs_out)
#plt.imshow(mask_array, cmap='Oranges')
plt.imshow(array, cmap='Greys')#, cmap='BuGn')
plt.colorbar( label="50% Detection Flux $10^{-17}$ erg/s/cm$^2$")
plt.clim(2.0, 20)
plt.xlabel('RA')
plt.ylabel('DEC')

for gal_region in gal_regions:
    gal_pix = gal_region.to_pixel(wcs_out)
    gal_pix.plot(ax=ax, color="blue", linewidth=2)
    
plt.scatter(detects_shot.ra*u.deg, detects_shot.dec*u.deg, transform=ax.get_transform('fk5'), s=5,
                      edgecolor='red', facecolor='none')
#plot up IFUnumber:
for i_ifu in np.arange( 0, np.size(ifu_name_list)):
    plt.text( ifu_ra[i_ifu]-0.012, ifu_dec[i_ifu]+0.008, ifu_name_list[i_ifu][9:],
              transform=ax.get_transform('fk5'), fontsize=12)

plt.title( "Flim" + str(datevobs) + ' at ' + str(int(wave[sel_slice])))
plt.savefig("FlimSlice" + str(datevobs)  + 'w' + str(int(wave[sel_slice])) + '.png')
