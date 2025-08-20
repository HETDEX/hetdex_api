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

from hetdex_api.shot import Fibers
from hetdex_api.extract import Extract
from hetdex_tools.interpolate import make_data_cube


def save_data_cube(
    shotid=None, ifuslot=None, filepath="datacubes", ffsky=False, F=None, E=None
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
        dwave=2,  # AA
        fill_value=0.0,
        ffsky=ffsky,
        apply_mask=True,
        convolve_image=False,
        include_error=True,
        include_bitmask=True,
        extract_class=E,
    )
    data    = hdu[1].data.astype(np.float32, copy=False)
    error   = hdu[2].data.astype(np.float32, copy=False)
    bitmask = hdu[3].data.astype(np.uint16,  copy=False)
    bitmask = bitmask.astype(np.uint16, copy=False)

    data_hdr  = hdu[1].header.copy()
    error_hdr = hdu[2].header.copy()
    mask_hdr = hdu[3].header.copy(); mask_hdr.pop('BUNIT', None)

    # Create compressed HDUs with TILE compression
    TILE = (16, 32, 32)

    data_hdu = fits.CompImageHDU(
        data, header=data_hdr, name="DATA", compression_type="RICE_1", tile_shape=TILE,
    )
    error_hdu = fits.CompImageHDU(
        error, header=error_hdr, name="ERROR", compression_type="RICE_1", tile_shape=TILE,
    )
    mask_hdu = fits.CompImageHDU(
        bitmask, header=mask_hdr, name="MASK", compression_type="RICE_1", tile_shape=TILE
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

    # Loop through all IFU slots for this shot
    for ifuslot in np.unique(F.ifuslot):
        save_data_cube(
            shotid=args.shotid,
            ifuslot=ifuslot,
            filepath=output_dir,
            ffsky=args.ffsky,
            F=F,
            E=E,
        )
    print(f"Done! Total time: {time.time() - start:.2f} seconds")
