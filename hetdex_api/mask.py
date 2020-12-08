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

def amp_flag_from_coords(coords, FibIndex, bad_amps_table, radius=3.*u.arcsec, shotid=None):
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

    
    fiber_table = FibIndex.query_region(coords,
                                        radius=radius,
                                        shotid=shotid)
    if np.size(fiber_table) > 0:
        mf_list = np.unique(fiber_table['multiframe']).astype(str)
        
        flags = []
        for mf in mf_list:

            sel = bad_amps_table['multiframe'] == mf
            
            if shotid is not None:
                sel_shot = bad_amps_table['shotid'] == shotid
                sel = sel*sel_shot

            flags.append(bad_amps_table['flag'][sel][0])

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
    
    sel = (bad_amps_table['shotid'] == shotid) * (bad_amps_table['multiframe'] == mf)

    return bad_amps_table['flag'][sel][0]


def amp_flag_from_closest_fiber(coords, FibIndex, bad_amps_table,
                                shotid=None,
                                maxdistance=8.*u.arcsec):
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

    fiberid = FibIndex.get_closest_fiberid(coords, shotid=shotid,
                                           maxdistance=maxdistance)
    if fiberid is not None:
        try:
            flag = amp_flag_from_fiberid(fiberid, bad_amps_table)
        except:
            flag = False
    else:
        flag = None
        
    return flag

    
def meteor_flag_from_coords(coords, shotid=None, streaksize=12.*u.arcsec):
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

    # meteors are found with +/- X arcsec of the line DEC=a+RA*b in this file

    met_tab = Table.read(config.meteor, format='ascii')
    sel_shot = met_tab['shotid'] == shotid

    if np.sum(sel_shot) > 0:
        a = met_tab['a'][sel_shot]
        b = met_tab['b'][sel_shot]

        ra_met = coords.ra + np.arange(-180, 180, 0.1)*u.arcsec
        dec_met = (a + ra_met.deg*b ) * u.deg

        met_coords = SkyCoord(ra=ra_met, dec=dec_met)

        meteor_match = met_coords.separation(coords) < streaksize

        
        if np.any(meteor_match):
            flag = False
        else:
            flag = True
    else:
        flag = True

    return flag


def create_gal_ellipse(galaxy_cat, row_index=None, pgcname=None, d25scale=3.):
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
        how many times D25 should to scale the region
    """
    if row_index is not None:
        index = row_index
    elif pgcname is not None:
        index = np.where(galaxy_cat['PGC'] == pgcname)[0][0]
        
    coords = SkyCoord(galaxy_cat['Coords'][index],frame='icrs')
    
    # The ellipse region uses the major and minor axes, so we have to multiply by
    # two first, before applying any user scaling.
    
    major = (galaxy_cat['SemiMajorAxis'][index]) * d25scale * u.arcmin
    minor = (galaxy_cat['SemiMinorAxis'][index]) * d25scale * u.arcmin
    pa    = (galaxy_cat['PositionAngle'][index]) * u.deg
    ellipse_reg = EllipseSkyRegion(center=coords, height=major, width=minor, angle=pa)

    return ellipse_reg


def create_dummy_wcs(coords, pixscale=0.5*u.arcsec, imsize=60.*u.arcmin):
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

    gridsize = imsize.to_value('arcsec')
    gridstep = pixscale.to_value('arcsec')
    
    # Create coordinate center
    ra_cen = coords.ra.deg
    dec_cen = coords.dec.deg
    
    ndim = np.int(2 * gridsize / gridstep + 1)
    center = ndim / 2
    w = wcs.WCS(naxis=2)
    w.wcs.crval = [ra_cen, dec_cen]
    w.wcs.crpix = [center, center]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cunit = ['deg','deg']
    w.wcs.cdelt = [-gridstep / gridsize, gridstep / gridsize]
    w.array_shape = [ndim, ndim]
    
    return w


def gal_flag_from_coords(coords, galaxy_cat, d25scale=3., nmatches=1):
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
        Experimentation showed a value of 1.75 might be more appropriate
    nmatches
        the closest nmatches are searched for.  nmatches = 1 means
        search the closest coordinate only.  nmatches = 3 is recommended

    Returns
    -------
    flag - boolean
        True if the source lies within the ellipse defined
        False if the source is not within the scaling

    """
   
    gal_coords = SkyCoord(galaxy_cat["Coords"])
   
    # calculate angular distances to all of the sources, and pick out the n closest ones
    d2d = coords.separation(gal_coords)
    ind1 = np.argsort(d2d)  # sort from closest to farthest away
    id_close = ind1[0:nmatches]

    # create fake WCS for regions use
    mywcs = create_dummy_wcs(coords)

    flag = False

    for idnum in id_close:
        ellipse = create_gal_ellipse(galaxy_cat, row_index=idnum, d25scale=d25scale)
        if ellipse.contains(coords, mywcs):
            flag = True

    return flag
