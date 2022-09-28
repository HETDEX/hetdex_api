# tools to convert HETDEX air wavelength to vacuum to calculate LyA redshifts
# Author: EMC
# Created: 2022-09-20

import numpy as np

from hetdex_api.survey import Survey
from astropy import time, coordinates as coord, units as u
from astropy import constants as const

S = Survey('hdr3')
mcdonald = coord.EarthLocation.of_site('mcdonald')

def get_bary_corr(shotid, units='km/s'):
    """
    For a given HETDEX observation, return the barycentric radial velocity correction.
    
    Parameters
    ----------
    shotid
    interger observation ID
    units
    string indicating units to be used. Must be readable to astropy units.
                Defaults to km/s

    Return
    ------
    vcor
    radial velocity correction in units
    """

    global S, mcdonald
    sel_shot = S.shotid == shotid
    coords = S.coords[sel_shot]
    mjds = S.mjd[sel_shot]
    t = time.Time(np.average( mjds), format='mjd')
    
    vcor = coords.radial_velocity_correction(kind='barycentric', obstime=t, location=mcdonald).to(units)
    
    return vcor.value
    
def convert_wave_vacuum_to_air(wave):
    """
    Converts a vacuum wavelength to the air wavelength.
    
    http://www.sdss3.org/dr9/spectro/spectro_basics.php
    
    Parameters
    ----------
    
    wave: vacuum wavelength in AA
    
    Return
    ------
    wave_air: air wavelength in AA
    
    """
    wave_air = wave / (1.0 + 2.735182*10**(-4) + 131.4182 / wave**2 + 2.76249*10**8 / wave**4)
    
    return wave_air


def get_wave_corr(wave_air, shotid=None, units='km/s'):
    """
    Converts an observed HETDEX wavelength that is measured in air from the
    the McDonald Observatory to a vacuum wavelength with a barycentric motion
    corretion applied if the shotid is given.

    Parameters
    ----------
    wave_air
    wavelength as measured in air at the observatory in AA.
    Can be a float or numpy float array
    shotid
    integer observation ID. If provided this will apply the correction
    for barycentric motion
    units
    string indicating units to be used. Must be readable to astropy units.
    Defaults to km/s
    """
    
    c = const.c.to(units).value
    
    if shotid is None:
        vcor=0
    else:
        vcor=get_bary_corr(shotid, units=units)
        
    wave_corr = wave_air+1.+(wave_air-3500)/4000.+vcor/(c / wave_air)
                                                                                    
    return wave_corr
