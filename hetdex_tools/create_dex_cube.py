import sys
import time
import argparse
import os
import os.path as op

import numpy as np
from astropy.io import fits
from astropy.io.fits import CompImageHDU
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
from astropy.table import Table

from hetdex_api.config import HDRconfig
from hetdex_api.shot import Fibers
from hetdex_api.extract import Extract
from hetdex_tools.interpolate import make_data_cube

import tables as tb

MANUAL_BIT = 13
MANUAL_FLAG = np.uint16(1 << MANUAL_BIT)  # 8192


def _paint_5x5x5(mask, z, y, x, flag_val):
    nz, ny, nx = mask.shape
    z0, z1 = max(z - 2, 0), min(z + 3, nz)
    y0, y1 = max(y - 2, 0), min(y + 3, ny)
    x0, x1 = max(x - 2, 0), min(x + 3, nx)
    mask[z0:z1, y0:y1, x0:x1] = np.bitwise_or(mask[z0:z1, y0:y1, x0:x1], flag_val)


def paint_manual_baddet_flags(bitmask, mask_hdr, ra, dec, wave):
    """
    Paint MANUAL_FLAG into a 5x5x5 region around each (ra, dec, wave).

    Parameters
    ----------
    bitmask : ndarray
        Cube mask array with shape (nz, ny, nx), dtype uint16
    mask_hdr : fits.Header
        Header for the MASK extension
    ra, dec, wave : array-like
        World coordinates for detections.
        wave must be in the same units expected by the cube WCS.
    """
    if len(ra) == 0:
        return bitmask

    w = WCS(mask_hdr)

    # world_to_array_index_values returns array indices in (z, y, x)
    z, y, x = w.world_to_array_index_values(ra, dec, wave)

    z = np.asarray(z, dtype=int)
    y = np.asarray(y, dtype=int)
    x = np.asarray(x, dtype=int)

    for zz, yy, xx in zip(z, y, x):
        _paint_5x5x5(bitmask, zz, yy, xx, MANUAL_FLAG)

    return bitmask


def save_data_cube(
    shotid=None,
    ifuslot=None,
    filepath="datacubes",
    ffsky=False,
    F=None,
    E=None,
    bad_det_source=None,
):
    """
    Create and save a datacube for a given IFU slot.

    Parameters
    ----------
    shotid : int
        The shot ID to process.
    ifuslot : str
        The IFU slot identifier. Defaults to None.
    filepath : str, optional
        Directory to save the output FITS files. Defaults to 'datacubes'.
    ffsky : bool, optional
        Whether to apply full sky fiber flat sky subtraction. Defaults to False (local sky subtraction).
    F : Fibers
        A `Fibers` class object containing fiber position data for the shot.
    E : Extract
        An `Extract` object with configuration for cube extraction.
    ADD_BADDET: bool
        This flag will open source catalog 5.0.2 and add baddets. Created 2026-04-02 by EMC. 

    Notes
    -----
    This function creates a directory structure of the form `filepath/shotid/`
    and saves three FITS files for each IFU slot: the datacube, error cube, and mask cube.
    """

    # Compute center coordinates for the IFU slot
    ra_cen = np.mean(F.ra[F.ifuslot == ifuslot])
    dec_cen = np.mean(F.dec[F.ifuslot == ifuslot])

    hdu = make_data_cube(
        shotid=shotid,
        coords=SkyCoord(ra=ra_cen, dec=dec_cen, unit="deg"),
        survey="hdr5",
        pixscale=0.5 * u.arcsec,
        imsize=52.0 * u.arcsec,
        dwave=2.0,  # AA
        fill_value=0.0,
        ffsky=ffsky,
        apply_mask=True,
        convolve_image=False,
        include_error=True,
        include_bitmask=True,
        extract_class=E,
    )

    data = hdu[1].data.astype(np.float32, copy=False)
    error = hdu[2].data.astype(np.float32, copy=False)
    bitmask = hdu[3].data.astype(np.uint16, copy=False)
    bitmask = bitmask.astype(np.uint16, copy=False)

    data_hdr = hdu[1].header.copy()
    error_hdr = hdu[2].header.copy()
    mask_hdr = hdu[3].header.copy()
    mask_hdr.pop("BUNIT", None)

    # add manual BADDET flag

    if bad_det_source is not None:
        sel_ifuslot = bad_det_source['ifuslot'] == ifuslot

        if np.sum(sel_ifuslot) > 0:
            ra_i = bad_det_source['ra'][sel_ifuslot]
            dec_i = bad_det_source['dec'][sel_ifuslot]
            wave_i = bad_det_source['wave'][sel_ifuslot]
            
            # convert wave to meters to match wcs
            wave_m = wave_i * 1e-10
        
            bitmask = paint_manual_baddet_flags(bitmask, mask_hdr, ra_i, dec_i, wave_m)
        
    # Create compressed HDUs with TILE compression
    TILE = (16, 32, 32)

    data_hdu = fits.CompImageHDU(
        data,
        header=data_hdr,
        name="DATA",
        compression_type="RICE_1",
        tile_shape=TILE,
    )
    error_hdu = fits.CompImageHDU(
        error,
        header=error_hdr,
        name="ERROR",
        compression_type="RICE_1",
        tile_shape=TILE,
    )
    mask_hdu = fits.CompImageHDU(
        bitmask,
        header=mask_hdr,
        name="MASK",
        compression_type="RICE_1",
        tile_shape=TILE,
    )

    # Create HDU list
    hdu_list = fits.HDUList([fits.PrimaryHDU(), data_hdu, error_hdu, mask_hdu])

    # Write compressed FITS file
    output_filename = op.join(filepath, f"dex_cube_{shotid}_{ifuslot}.fits")
    hdu_list.writeto(output_filename, overwrite=True, checksum=True)

    return hdu_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create and save HETDEX datacubes.")
    parser.add_argument("shotid", type=int, help="Shot ID to process")
    parser.add_argument(
        "--filepath",
        type=str,
        default="datacubes",
        help="Directory to save output files",
    )
    parser.add_argument(
        "--ffsky", action="store_true", help="Apply fiber flat sky subtraction"
    )

    args = parser.parse_args()

    start = time.time()
    print("Working on", args.shotid)

    # Load fiber and extraction data
    E = Extract()
    E.load_shot(args.shotid, survey="hdr5")

    # pull out IFUslot info
    E.fibers.ifuslot = np.array([mf.split("_")[2] for mf in E.fibers.multiframe])
    F = E.fibers

    # Ensure shotid-specific directory exists
    output_dir = os.path.join(args.filepath, str(args.shotid))
    os.makedirs(output_dir, exist_ok=True)

    # open 5.0.2 source catalog h5 file                                                                          
    config = HDRconfig()
    version = "5.0.2"
    hdrv = "hdr{}".format(version[0])

    catfile = op.join(config.hdr_dir[hdrv], "catalogs", "source_catalog_" + version + ".h5")

    source_tableh5 = tb.open_file(catfile, "r")

    # read source table for given shotid (its not indexed on ifuslot atm)
    shot_i = args.shotid
    
    source_table = Table(source_tableh5.root.SourceCatalog.read_where("shotid == shot_i"))

    flag_best_current_model = (
        source_table["flag_badamp"]
        & source_table["flag_badfib"]
        & source_table["flag_badpix"]
        & source_table["flag_wave"]
        & source_table["flag_3540"]
        & source_table["flag_continuum"]
        & source_table["flag_baddet"]
        & source_table["flag_meteor"]
        & source_table["flag_largegal"]
        & source_table["flag_chi2fib"]
        & source_table["flag_pixmask"]
        & source_table["flag_satellite"]
        & source_table["flag_cal"])

    sel_baddets = (flag_best_current_model == 1) & (source_table["det_type"] == "line") & (source_table["flag_best"] == 0)
    source_baddets = source_table['detectid','ifuslot','ra','dec','wave'][ sel_baddets]

    print(f'Manually masking {np.sum(sel_baddets)} bad detectids for shotid={args.shotid}')
    # Loop through all IFU slots for this shot
    for ifuslot in np.unique(F.ifuslot):
        save_data_cube(
            shotid=args.shotid,
            ifuslot=ifuslot,
            filepath=output_dir,
            ffsky=args.ffsky,
            F=F,
            E=E,
            bad_det_source = source_baddets,
        )
    print(f"Done! Total time: {time.time() - start:.2f} seconds")
