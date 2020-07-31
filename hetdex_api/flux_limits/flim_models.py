"""

This module stores different models 
to convert between the values in the 
sensitivity cubes and the flux at
50% detection completeness

.. moduleauthor:: Daniel Farrow <dfarrow@mpe.mpg.de>

"""

from numpy import polyval

def hdr2pt1_f50_from_noise(noise, sncut):
    """
    Return the 50% completeness
    flux given noise and 
    S/N cut. This uses an empirical
    model calibrated on simulated
    LAE detections inserted into
    ~200 different shots (see
    Figures in data release 
    document)

    Parameters
    ----------
    noise : float
        the noise from the
        sensitivity cubes
    sncut : float
        the signal to noise
        cut to assume

    Returns
    -------
    f50s : array
       the fluxes at 50%
       completeness
    """

    if sncut < 4.5 or sncut > 7.5:
        print("WARNING: model not calibrated for this S/N range")

    snslope=1.0
    intercept_poly=[0.11268546,  -1.75103671,
                    9.00946424, -15.71204789]
    intercept = 1e-17*polyval(intercept_poly, sncut)

    f50s = snslope*sncut*noise + intercept
        
    return f50s

def hdr1_f50_from_noise(noise, sncut):
    """
    Return the 50% completeness
    flux given noise and 
    S/N cut. This just assumes
    50% completeness is at
    sncut*noise 

    Parameters
    ----------
    noise : float
        the noise from the
        sensitivity cubes
    sncut : float
        the signal to noise
        cut to assume

    Returns
    -------
    f50s : array
       the fluxes at 50%
       completeness
    """

    f50s = sncut*noise 

    return f50s

