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
    radius
        radius to search for fibers
    shotid
        shotid to search. If none it will search all shots at once. If
        any are flagged bad then it will return False for all.
    
    """

    fiber_region = FibIndex.query_region(coords,
                                         radius=radius,
                                         shotid=shotid)
    if np.size(fiber_region) > 0:
        mf_list = np.unique(fiber_table['multiframe']).astype(str)
        
        flags = []
        for mf in mf_list:

            sel = bad_amps_table['multiframe'] == mf
            
            if shotid is not None:
                sel_shot = bad_amps_table['shotid'] == shotid
                sel = sel*sel_shot

            flags.append(bad_amps_table['flag'][sel][0])

        amp_flag = np.all(flag)
    else:
        amp_flag = None

    return amp_flag

def amp_flag_from_fiberid(fiberid, bad_amps_table):
    shotid = int(fiberid[0:11])
    mf = fiberid[14:34]
    
    sel = (bad_amps_table['shotid'] == shotid) * (bad_amps_table['multiframe'] == mf)

    return bad_amps_table['flag'][sel][0]
