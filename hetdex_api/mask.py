# -*- coding: utf-8 -*-
"""

API to masking HETDEX data products

Created on 2020/09/08
Update on 2020/11/10: Updating to create galaxy regions from rc3 catalog

@author: Erin Mentuch Cooper

"""

from __future__ import print_function

import numpy as np

from astropy.table import Table, join
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy import wcs

from regions import EllipseSkyRegion, EllipsePixelRegion

from hetdex_api.config import HDRconfig
from hetdex_api.survey import FiberIndex

config = HDRconfig()


def amp_flag_from_coords(coords, FibIndex, bad_amps_table, radius=None, shotid=None):
    """
    Returns a boolean flag whether the amp has been flagged usable

    Parameters
    ----------
    coords
        an astropy.coordinates SkyCoord object
    FibIndex
        a hetdex_api.survey FiberIndex class object
    bad_amps_table
        astropy table containing the bad amp flag values. This can
        be retrieved from config.badamp
    radius
        radius to search for fibers
    shotid
        shotid to search. If none it will search all shots at once. If
        any are flagged bad then it will return False for all.

    Returns
    -------
    None if no matching fiber is found in the aperture

    True if all fibers in aperture are on good amps

    False if any fiber in the defined region is flagged in a bad amp

    Examples
    --------
    from hetdex_api.config import HDRconfig
    from hetdex_api.mask import *

    config = HDRconfig()
    bad_amps_table = Table.read(config.badamp)
    FibIndex = FiberIndex()

    coords = SkyCoord(11.628530 * u.deg, 0.081790 * u.deg, frame='icrs')

    flag = amp_flag_from_coords(coords, FibIndex,  bad_amps_table,
                                radius=1.5*u.arcsec, shotid=20181003009)

    """

    if radius is None:
        radius = 3.0 * u.arcsec

    fiber_table = FibIndex.query_region(coords, radius=radius, shotid=shotid)
    if np.size(fiber_table) > 0:
        mf_list = np.unique(fiber_table["multiframe"]).astype(str)

        flags = []
        for mf in mf_list:
            sel = bad_amps_table["multiframe"] == mf

            if shotid is not None:
                sel_shot = bad_amps_table["shotid"] == shotid
                sel = sel * sel_shot

            flags.append(bad_amps_table["flag"][sel][0])

        amp_flag = np.all(flags)
    else:
        amp_flag = None

    return amp_flag


def amp_flag_from_fiberid(fiberid, bad_amps_table):
    """
    Returns a boolean flag whether the amp has been flagged usable

    Parameters
    ----------
    fiberid
        string with the unique fiber_id info

    bad_amps_table
        astropy table containing the bad amp flag values. This can
        be retrieved from config.badamp

    Returns
    -------
    None if no matching fiber is found to the coords

    True if amp is usable

    False if any fiber in the defined region is flagged in a bad amp

    """

    shotid = int(fiberid[0:11])
    mf = fiberid[14:34]

    sel = (bad_amps_table["shotid"] == shotid) * (bad_amps_table["multiframe"] == mf)

    if np.sum(sel) > 0:
        return bad_amps_table["flag"][sel][0]
    else:
        print('No match for amplifier in amp_flag.fits', shotid, mf)
        return True

def amp_flag_from_closest_fiber(
    coords, FibIndex, bad_amps_table, shotid=None, maxdistance=None
):
    """
    Function to retrieve the amp flag for the closest fiberid in a shot

    Parameters
    ----------
    self
        the FiberIndex class for a specific survey
    coords
        coordinate you want to search for the closest fiber.
        This should be an astropy SkyCoord object
    FibIndex
        a hetdex_api.survey FiberIndex class object
    bad_amps_table
        astropy table containing the bad amp flag values. This can
        be retrieved from config.badamp
    shotid
        Specific shotid (dtype=int) you want
    maxdistance
        The max distance you want to search for a nearby fiber.
        Default is 8.*u.arcsec

    Returns
    -------
    bool
        None if no matching fiber is found to the coords
        True if amp is usable
        False if any fiber in the defined region is flagged in a bad amp

    """
    if maxdistance is None:
        maxdistance = 8.0 * u.arcsec

    fiberid = FibIndex.get_closest_fiberid(
        coords, shotid=shotid, maxdistance=maxdistance
    )
    if fiberid is not None:
        try:
            flag = amp_flag_from_fiberid(fiberid, bad_amps_table)
        except:
            flag = False
    else:
        flag = None

    return flag


def cal_flag5200_for_amp(amp, shotid, expnum=None, cal5200_tab=None):
    """
    Function to indicate whether to flag the amp at 5200
    """

    if cal5200_tab is None:
        # open table if not yet opened
        cal5200_tab = Table.read(
            config.cal5200, format="ascii", names=["shotid", "multiframe", "expnum"]
        )

    flag = True

    if expnum is not None:
        sel_row = (
            (cal5200_tab["shotid"] == shotid)
            * (cal5200_tab["multiframe"] == amp)
            * (cal5200_tab["expnum"] == expnum)
        )
    else:
        sel_row = (cal5200_tab["shotid"] == shotid) * (cal5200_tab["multiframe"] == amp)

    if np.sum(sel_row) > 0:
        flag = False

    return flag


def cal_flag5460_for_amp(amp, shotid, expnum=None, cal5460_tab=None):
    """
    Function to indicate whether to flag the amp at 5460
    """
    global config

    if cal5460_tab is None:
        # open table if not yet opened
        cal5460_tab = Table.read(
            config.cal5460, format="ascii", names=["shotid", "multiframe", "expnum"]
        )

    flag = True

    if expnum is not None:
        sel_row = (
            (cal5460_tab["shotid"] == shotid)
            * (cal5460_tab["multiframe"] == amp)
            * (cal5460_tab["expnum"] == expnum)
        )
    else:
        sel_row = (cal5460_tab["shotid"] == shotid) * (cal5460_tab["multiframe"] == amp)

    if np.sum(sel_row) > 0:
        flag = False

    return flag


def cal_flag_for_amp_wave(
    wave, amp, shotid, expnum=None, cal5460_tab=None, cal5200_tab=None
):
    """
    For an input amp/shotid/expnum and wavelength flag whether the detection wavelength should be flagged
    """
    flag = True
    flag5200 = cal_flag5200_for_amp(amp, shotid, expnum, cal5200_tab=cal5200_tab)

    if flag5200 == False:
        if ((wave >= 5194) * (wave <= 5197)) | ((wave >= 5200) * (wave <= 5205)):
            flag = False

    flag5460 = cal_flag5460_for_amp(amp, shotid, expnum, cal5460_tab=cal5460_tab)

    if flag5460 == False:
        if (wave >= 5456) * (wave <= 5466):
            flag = False

    # flag for 3540 sky line
    if (wave >= 3534) * (wave <= 3556):
        flag = False
    
    return flag


def meteor_flag_from_coords(coords, shotid=None, streaksize=None):
    """
    Returns a boolean flag value to mask out meteors

    Parameters
    ----------
    coords
        an astropy.coordinates SkyCoord object
    shotid
        shotid to search. If none it will search all shots at once. If
        any are flagged bad then it will return False for all.
    streaksize
        an astropy quantity object defining how far off the
        perpendicular line of the meteor streak to mask out. Default
        is 12*u.arcsec

    Returns
    -------
    bool
        True if no meteors fall in the aperture
        False if a meteor falls in the aperture

    Example
    -------

    """

    global config

    if streaksize is None:
        streaksize = 12.0 * u.arcsec
    # meteors are found with +/- X arcsec of the line DEC=a+RA*b in this file

    met_tab = Table.read(config.meteor, format="ascii")
    sel_shot = met_tab["shotid"] == shotid

    if np.sum(sel_shot) > 0:
        a = met_tab["a"][sel_shot]# this is the intercept
        b = met_tab["b"][sel_shot]# this is the slope

        if np.abs(b) < 500:
            ra_met = coords.ra + np.arange(-180, 180, 0.1) * u.arcsec
        else:# for very large slopes we needed higher resolution
            ra_met = coords.ra + np.arange(-30, 30, 0.001) * u.arcsec
        dec_met = (a + ra_met.deg * b) * u.deg

        # added 20241115 to handle vertical streak on 
        sel_met = (dec_met > -90*u.deg) * (dec_met< 90*u.deg)

        met_coords = SkyCoord(ra=ra_met[sel_met], dec=dec_met[sel_met])

        meteor_match = met_coords.separation(coords) < streaksize

        if np.any(meteor_match):
            flag = False
        else:
            flag = True
    else:
        flag = True

    return flag


def satellite_flag_from_coords(coords, shotid=None, streaksize=None):
    """
    Returns a boolean flag value to mask out satellites

    Parameters
    ----------
    coords
        an astropy.coordinates SkyCoord object
    shotid
        shotid to search. If none it will search all shots at once. If
        any are flagged bad then it will return False for all.
    streaksize
        an astropy quantity object defining how far off the
        perpendicular line of the meteor streak to mask out. Default
        is 6*u.arcsec

    Returns
    -------
    bool
        True if no satellite track falls in the aperture
        False if a satellite track falls in the aperture

    Example
    -------

    """

    global config

    flag = True

    if streaksize is None:
        streaksize = 6.0 * u.arcsec
    # satellites are found with +/- X arcsec of the line DEC=a+RA*b in this file

    sat_tab = Table.read(
        config.satellite,
        format="ascii",
        names=["shotid", "expnum", "slope", "intercept"],
    )
    sel_shot = sat_tab["shotid"] == shotid

    if np.sum(sel_shot) > 0:
        for row in sat_tab[sel_shot]:
            slope = row["slope"]
            intercept = row["intercept"]

            ra_sat = coords.ra + np.arange(-180, 180, 0.1) * u.arcsec
            dec_sat = (intercept + ra_sat.deg * slope) * u.deg

            sat_coords = SkyCoord(ra=ra_sat, dec=dec_sat)

            sat_match = sat_coords.separation(coords) < streaksize

            if np.any(sat_match):
                flag = False
                return flag

    return flag


def create_gal_ellipse(galaxy_cat, row_index=None, pgcname=None, d25scale=1.5):
    """
    Similar to galmask.py/ellreg but can take a galaxy name as input.

    Create galaxy ellipse region using an input galaxy catalog_table (likely need
    to change this API as it heavily follows the RC3 catalog format)

    Parameters
    ----------
    galaxy_catalog_table: an astropy table
        table of nearby galaxy positions and sizes. Must have a central
        coordinate and the SemiMajorAxis, SemiMinorAxis, Position Angle
        info as table column names
    row_index: int
        a row index in the galaxy catalog to create the ellipse for
    pgcname: str
        the PGCNAME in the RC3 cat. This is a string. eg. "PGC 43255" for NGC 4707
    d25scale: float
        how many times D25 should to scale the region. Default is 1.5 based on
        visual exploration
    """
    if row_index is not None:
        index = row_index
    elif pgcname is not None:
        index = np.where(galaxy_cat["PGC"] == pgcname)[0][0]

    coords = SkyCoord(galaxy_cat["Coords"][index], frame="icrs")

    # The ellipse region uses the major and minor axes, so we have to multiply by
    # two first, before applying any user scaling.

    major = 2.0 * (galaxy_cat["SemiMajorAxis"][index]) * d25scale * u.arcmin
    minor = 2.0 * (galaxy_cat["SemiMinorAxis"][index]) * d25scale * u.arcmin
    pa = (galaxy_cat["PositionAngle"][index]) * u.deg
    ellipse_reg = EllipseSkyRegion(center=coords, height=major, width=minor, angle=pa)

    return ellipse_reg


def create_dummy_wcs(coords, pixscale=None, imsize=None):
    """
    Create a simple fake WCS in order to use the regions subroutine.
    Adapted from John Feldmeiers galmask.py

    Parameters
    ----------
    coords: a SkyCoord object
        center coordinates of WCS
    pixscale: astropy quantity
        pixel scale of WCS in astropy angle quantity units
    imsize: astropy quantity
        size of WCS in astropy angle quanity units
    """
    if imsize is None:
        imsize = 60.0 * u.arcmin
    if pixscale is None:
        pixscale = 0.5 * u.arcsec

    gridsize = imsize.to_value("arcsec")
    gridstep = pixscale.to_value("arcsec")

    # Create coordinate center
    ra_cen = coords.ra.deg
    dec_cen = coords.dec.deg

    ndim = int(2 * gridsize / gridstep + 1)
    center = ndim / 2
    w = wcs.WCS(naxis=2)
    w.wcs.crval = [ra_cen, dec_cen]
    w.wcs.crpix = [center, center]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cunit = ["deg", "deg"]
    w.wcs.cdelt = [-gridstep / gridsize, gridstep / gridsize]
    w.array_shape = [ndim, ndim]

    return w


def gal_flag_from_coords(coords, galaxy_cat, d25scale=1.5):
    """
    Returns a boolean flag value to mask sources near large galaxies

    Adapted from John Feldmeier's hetdex_tools/galmask.py

    Parameters
    ----------
    coords
        an astropy.coordinates SkyCoord object
    galaxy_cat
        an astropy table containing the large galaxy parameters. This
        is catered for the RC3 catalog stored in config.rc3cat
    d25scale
        The scaling of ellipses.  1.0 means use the ellipse for D25.
        Experimentation shows value of 1.5 works well

    Returns
    -------
    flag - boolean
        True if the source is not in the galaxy mask
        False if the source is within the scaling of the galaxy mask

    """

    gal_coords = SkyCoord(galaxy_cat["Coords"])

    # calculate angular distances to all of the sources, and pick out the n closest ones
    d2d = coords.separation(gal_coords)

    sel_rows = np.where(d2d < 1.0 * u.deg)[0]

    # create fake WCS for regions use
    mywcs = create_dummy_wcs(coords)

    flag = True

    if len(sel_rows) > 0:
        for row_index in sel_rows:
            ellipse = create_gal_ellipse(
                galaxy_cat, row_index=row_index, d25scale=d25scale
            )
            if ellipse.contains(coords, mywcs):
                flag = False

    return flag
