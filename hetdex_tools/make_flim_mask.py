# -*- coding: utf-8 -*-
"""
Created: 2020/09/09

@author: Erin Mentuch Cooper

Quick script to create h5 files of flim masks
"""

import sys
import numpy as np
import os.path as op

import tables as tb
from astropy.wcs import WCS
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord

from hetdex_api.survey import Survey, FiberIndex
from hetdex_api.config import HDRconfig
from hetdex_api.mask import *
from hetdex_api.flux_limits.hdf5_sensitivity_cubes import *

date = sys.argv[1]
obs = sys.argv[2]

datevshot = str(date) + "v" + str(obs).zfill(3)
shotid = int(str(date) + str(obs).zfill(3))

# open the config class and FiberIndex class so they don't get repeatedly initiated
config = HDRconfig()

LATEST_HDR_NAME = HDRconfig.LATEST_HDR_NAME

# check if file is in badshot list
badshot = np.loadtxt(config.badshot, dtype=int)

if shotid in badshot:
    print("Shot is in badshot list. Making mask zero everywhere")
    badshot = True
else:
    badshot = False

# also check if shot is in bad throughput list
# and set mask to zero everywhere
badtpshots = np.loadtxt(config.lowtpshots, dtype=int)
if shotid in badtpshots:
    badshot = True
    print("Shot has bad throughput. Setting flux limit mask to 0")

bad_amps = Table.read(config.badamp)

FibIndex = FiberIndex()

try:
    hdf_filename = return_sensitivity_hdf_path(datevshot,
                                               release=LATEST_HDR_NAME)
except NoFluxLimsAvailable:
    sys.exit("No flux limit file found for " + datevshot)

flimhdf = SensitivityCubeHDF5Container(
    filename=hdf_filename, aper_corr=1.0, flim_model="hdr2pt1"
)

hdf_outfilename = datevshot + "_mask.h5"

fileh = tb.open_file(hdf_outfilename, "w")

groupMask = fileh.create_group(fileh.root, "Mask", "Flux limit masks")

# check if there are any meteors in the shot:
met_tab = Table.read(config.meteor, format="ascii")

if shotid in met_tab["shotid"]:
    check_meteor = True
else:
    check_meteor = False

# down-select bad_amps_table
sel_shot = bad_amps["shotid"] == shotid
bad_amps_table = bad_amps[sel_shot]

# get ifu list
ifu_list = []

for row in bad_amps_table:
    ifu_list.append(int(row["multiframe"][10:13]))

bad_amps_table.add_column(Column(np.array(ifu_list), name="ifuslot"))

# grab galaxy region files if any near shot:
check_gal = False
galaxy_cat = Table.read(config.rc3cat, format='ascii')
gal_coords = SkyCoord(galaxy_cat['Coords'], frame='icrs')

survey = Survey()
shot_coords = survey.coords[survey.shotid == shotid][0]
sel_reg = np.where(shot_coords.separation(gal_coords) < 1.*u.deg)[0]

if len(sel_reg) > 0:
    gal_regions = []
    for idx in sel_reg:
        gal_regions.append( create_gal_ellipse(galaxy_cat, row_index=idx, d25scale=3.))
    check_gal = True

for ifu_name, tscube in flimhdf.itercubes():
    print("Working on " + datevshot + "_" + ifu_name)
    ifuslot = int(ifu_name[8:11])

    # check to see if there are any bad amps
    sel_ifu = bad_amps_table["ifuslot"] == ifuslot

    if np.any(bad_amps_table["flag"][sel_ifu] == 0):
        check_amp = True
    elif np.all(bad_amps_table["flag"][sel_ifu] == 0):
        mask = np.zeros_like(slice_, dtype=int)
        check_amp = False
        check_meteor = False
    else:
        check_amp = False

    wcs = tscube.wcs.celestial
    wslice = 400
    sncut = 6
    slice_ = tscube.f50_from_noise(tscube.sigmas[wslice, :, :], sncut)

    if badshot:
        mask = np.zeros_like(slice_, dtype=int)
        check_amp = False
        check_meteor = False
    else:
        mask = np.ones_like(slice_, dtype=int)

    if check_amp or check_meteor:
        for i in np.arange(0, 31):
            for j in np.arange(0, 31):

                ra, dec = wcs.wcs_pix2world(i, j, 0)
                coords = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
                
                if check_amp:
                    flag_amp = amp_flag_from_closest_fiber(
                        coords,
                        FibIndex,
                        bad_amps_table,
                        maxdistance=11.0 * u.arcsec,
                        shotid=shotid,
                    )
                else:
                    flag_amp = True

                if check_meteor:
                    flag_meteor = meteor_flag_from_coords(coords, shotid)
                    if flag_amp is not None:
                        mask[j, i] = flag_amp * flag_meteor
                    else:
                        mask[j, i] = flag_meteor
                else:
                    if flag_amp is not None:
                        mask[j, i] = flag_amp

    if check_gal:
        # make galaxy mask
        NAXIS1, NAXIS2 = np.shape(slice_)
        
        x = np.arange(NAXIS1)
        y = np.arange(NAXIS2)
        X, Y = np.meshgrid(x, y)
        gal_mask = np.ones_like(X)
        ra, dec = wcs.wcs_pix2world(X, Y, 0)
        
        for gal_region in gal_regions:
            dummy_wcs = create_dummy_wcs(gal_region.center,
                                         imsize=2*gal_region.height)
            gal_mask = gal_mask * np.invert(gal_region.contains(
                SkyCoord(ra, dec, unit='deg'),
                dummy_wcs))

        mask = mask * gal_mask
            
    # append array to h5 file
    fileh.create_array(groupMask, ifu_name, mask)

flimhdf.close()
FibIndex.hdfile.close()
fileh.close()
