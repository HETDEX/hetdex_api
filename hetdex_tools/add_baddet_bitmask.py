#!/usr/bin/env python
# coding: utf-8

# This code was added to place a 5x5x5 bitmask at the location of
# a bad detection that is not identified by the previous fiber mask
# options


from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table, unique
import numpy as np
import os
import os.path as op
import tables as tb
from hetdex_api.config import HDRconfig
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
import sys
import time

MANUAL_BIT = 13
MANUAL_FLAG = np.uint16(1 << MANUAL_BIT)  # 8192


def _paint_3x3x3(mask, z, y, x, flag_val):
    nz, ny, nx = mask.shape
    z0, z1 = max(z - 1, 0), min(z + 2, nz)
    y0, y1 = max(y - 1, 0), min(y + 2, ny)
    x0, x1 = max(x - 1, 0), min(x + 2, nx)
    mask[z0:z1, y0:y1, x0:x1] = np.bitwise_or(mask[z0:z1, y0:y1, x0:x1], flag_val)


def _paint_5x5x5(mask, z, y, x, flag_val):
    nz, ny, nx = mask.shape

    z0, z1 = max(z - 2, 0), min(z + 3, nz)
    y0, y1 = max(y - 2, 0), min(y + 3, ny)
    x0, x1 = max(x - 2, 0), min(x + 3, nx)

    mask[z0:z1, y0:y1, x0:x1] = np.bitwise_or(mask[z0:z1, y0:y1, x0:x1], flag_val)


def _as_arrays(table_or_arrays, ra_col="ra", dec_col="dec", wave_col="wave"):
    """
    Returns (ra, dec, wave) as 1D numpy arrays.
    """
    if isinstance(table_or_arrays, tuple) or isinstance(table_or_arrays, list):
        if len(table_or_arrays) != 3:
            raise ValueError("Pass (ra, dec, wave) as a 3-tuple/list.")
        ra, dec, wave = table_or_arrays
        return np.asarray(ra), np.asarray(dec), np.asarray(wave)

    # Table / DataFrame / dict-like
    ra = np.asarray(table_or_arrays[ra_col])
    dec = np.asarray(table_or_arrays[dec_col])
    wave = np.asarray(table_or_arrays[wave_col])
    return ra, dec, wave


def add_manual_flags_to_cube(
    shotid,
    ifuslot,
    table_or_arrays,
    output_fits=None,
    ra_col="ra",
    dec_col="dec",
    wave_col="wave",
    tile_shape=(16, 32, 32),
    compression_type="RICE_1",
    chunk_size=20000,
):
    """
    Read existing cube, paint MANUAL bit (bit 13) in MASK for each (ra,dec,wave),
    then write a new cube.

    chunk_size keeps memory stable if you have a lot of points.
    """
    input_fits = f"{pdr_dir}/datacubes/{shot_i}/dex_cube_{shot_i}_{ifuslot}.fits"
    output_fits = f"datacubes/{shot_i}/dex_cube_{shot_i}_{ifuslot}.fits"

    # Make sure output directory exists
    outdir = op.dirname(output_fits)
    if outdir and not op.exists(outdir):
        os.makedirs(op.dirname(output_fits), exist_ok=True)

    ra, dec, wave = _as_arrays(table_or_arrays, ra_col, dec_col, wave_col)
    if not (ra.shape == dec.shape == wave.shape):
        raise ValueError("ra/dec/wave must have the same shape/length.")

    with fits.open(input_fits, memmap=True) as hdul:
        data_hdu = hdul["DATA"]
        err_hdu = hdul["ERROR"]
        mask_hdu = hdul["MASK"]

        data = np.array(data_hdu.data, copy=False).astype(np.float32, copy=False)
        error = np.array(err_hdu.data, copy=False).astype(np.float32, copy=False)
        bitmask = np.array(mask_hdu.data, copy=True).astype(np.uint16, copy=False)

        data_hdr = data_hdu.header.copy()
        err_hdr = err_hdu.header.copy()
        mask_hdr = mask_hdu.header.copy()
        mask_hdr.pop("BUNIT", None)

    # WCS from existing cube
    w = WCS(mask_hdr)

    # Apply in chunks
    n = len(ra)
    for i0 in range(0, n, chunk_size):
        i1 = min(i0 + chunk_size, n)

        # returns x,y,z in pixel coordinates (0-based)
        # x, y, z = w.all_world2pix(ra[i0:i1], dec[i0:i1], wave[i0:i1], 0)
        z, y, x = w.world_to_array_index_values(ra, dec, wave)

        # paint expects integer array indices in (z,y,x)
        z = np.rint(z).astype(int)
        y = np.rint(y).astype(int)
        x = np.rint(x).astype(int)

        for zz, yy, xx in zip(z.astype(int), y.astype(int), x.astype(int)):
            # _paint_3x3x3(bitmask, zz, yy, xx, MANUAL_FLAG)
            _paint_5x5x5(bitmask, zz, yy, xx, MANUAL_FLAG)

    out_hdul = fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.CompImageHDU(
                data,
                header=data_hdr,
                name="DATA",
                compression_type=compression_type,
                tile_shape=tile_shape,
            ),
            fits.CompImageHDU(
                error,
                header=err_hdr,
                name="ERROR",
                compression_type=compression_type,
                tile_shape=tile_shape,
            ),
            fits.CompImageHDU(
                bitmask,
                header=mask_hdr,
                name="MASK",
                compression_type=compression_type,
                tile_shape=tile_shape,
            ),
        ]
    )
    out_hdul.writeto(output_fits, checksum=True, overwrite=True)
    return output_fits


t0 = time.time()

config = HDRconfig()
version = "5.0.2"
hdrv = "hdr{}".format(version[0])

catfile = op.join(config.hdr_dir[hdrv], "catalogs", "source_catalog_" + version + ".h5")
source_tableh5 = tb.open_file(catfile, "r")

pdr_dir = "/corral-repl/utexas/Hobby-Eberly-Telesco/public/HETDEX/internal/pdr1"
# pdr_dir = '/home/jovyan/Hobby-Eberly-Public/HETDEX/internal/pdr1'
ifu_data = Table.read(op.join(pdr_dir, "ifu-index.fits"))

shotlist = np.unique(ifu_data["shotid"])

shot_i = int(sys.argv[1])
print("Working on ", shot_i)
# read source table for given shotid (its not indexed on ifuslot atm)
source_table = Table(source_tableh5.root.SourceCatalog.read_where("shotid == shot_i"))

# No need to flag sources flagged for other reasons already
flag_best_current_model = (
    source_table["flag_badamp"]
    * source_table["flag_badfib"]
    * source_table["flag_badpix"]
    * source_table["flag_wave"]
    * source_table["flag_3540"]
    * source_table["flag_continuum"]
    * source_table["flag_baddet"]
    * source_table["flag_meteor"]
    * source_table["flag_largegal"]
    * source_table["flag_chi2fib"]
    * source_table["flag_pixmask"]
    * source_table["flag_satellite"]
    * source_table["flag_cal"]
)

flag_best_baddet = (
    source_table["flag_raic"]
    * source_table["flag_lowlw"]
    * np.invert(source_table["flag_dee"] == 0)
    * source_table["flag_faint_cont"]
)

ifu_list = np.unique(source_table["ifuslot"])

for ifuslot in ifu_list:
    ifu_baddets = (
        (flag_best_current_model == 1)
        & (source_table["det_type"] == "line")
        & (source_table["ifuslot"] == ifuslot)
        & (source_table["flag_best"] == 0)
    )

#    print(
#        "Number of baddets not yet flagged for ifuslot={} is {}.".format(
#            ifuslot, np.sum(ifu_baddets)
#        )
#    )
    ra_i = source_table["ra"][ifu_baddets]
    dec_i = source_table["dec"][ifu_baddets]
    wave_i = source_table["wave"][ifu_baddets]
    sn_i = source_table["sn"][ifu_baddets]

    # convert wave to meters to match wcs
    wave_m = wave_i * 1e-10
    out = add_manual_flags_to_cube(shot_i, ifuslot, (ra_i, dec_i, wave_m))

t1 = time.time()
print("Done shotid={} in {:3.2f}m".format(shot_i, (t1 - t0) / 60.0))
source_tableh5.close()
