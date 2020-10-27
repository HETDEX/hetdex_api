# -*- coding: utf-8 -*-
"""

API to masking HETDEX data products

Created on 2020/09/08

@author: Erin Mentuch Cooper

"""

from __future__ import print_function

import numpy as np

from astropy.table import Table, join
import astropy.units as u
from astropy.coordinates import SkyCoord

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

    
def meteor_flag_from_coords(coords, shotid=None, streaksize=9*u.arcsec):
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
        is 9*u.arcsec
    
    Returns
    -------
    bool
        True if no meteors fall in the aperture
        False if a meteor falls in the aperture

    Example
    -------
    
    """

    global config

    # meteors are found with +/- 8 arcsec of the line DEC=a+RA*b in this file
    met_tab = Table.read('/work/05350/ecooper/wrangler/mask/meteor.list', format='ascii')
    sel_shot = met_tab['shotid'] == shotid

    if np.sum(sel_shot) > 0:
        a = met_tab['a'][sel_shot]
        b = met_tab['b'][sel_shot]

        ra_met = coords.ra + np.arange(-90,90)*u.arcsec
        dec_met = (a + ra_met.deg*b ) * u.deg
        
        met_coords = SkyCoord(ra=ra_met, dec=dec_met)

        meteor_match = coords.separation(met_coords) < streaksize

        if np.any(meteor_match):
            flag = False
        else:
            flag = True
    else:
        flag = True

    return flag
    
