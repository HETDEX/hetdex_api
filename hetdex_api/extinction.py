"""
Scripts to deal with extinction.

This includes HDR2.1 flat extinction fix

"""

import numpy as np
from scipy import interpolate
import extinction
from hetdex_api.config import HDRconfig

LATEST_HDR_NAME = HDRconfig.LATEST_HDR_NAME
config = HDRconfig(LATEST_HDR_NAME)

from dustmaps.config import config as dustmaps_config
if dustmaps_config['data_dir'] is None:
    dustmaps_config['data_dir'] = config.dustmaps        
from dustmaps.sfd import SFDQuery


def dustmaps_setup():
    # this is how I intially set up dustmaps
    # Need to put maps in config.dustmaps

    from dustmaps.config import config as dustmaps_config
    import dustmaps
    config = HDRconfig()
    dustmaps_config['data_dir'] = config.dustmaps
    dustmaps.sfd.fetch() #don't need to do this


def get_2pt1_extinction_fix(pad=True):
    """
    This is to fix the E(B-V)=0.02
    flat extinction applied to 2.1

    Curve from Karl
    /work/00115/gebhardt/fluxcor/extinction

    Paramaters
    ----------
    pad  bool
    This is to pad curve to 3470 and 5400 to match
    wave_rect
    """
    config = HDRconfig()

    karl_data = np.loadtxt( config.extinction_fix)

    if pad:
        karl_wave = np.concatenate(([3450.],
                                    karl_data[:, 0],
                                    [5550.]))
        karl_curve = np.concatenate(([karl_data[0, 1]],
                                     karl_data[:, 1],
                                     [karl_data[-1, 1]]))
    else:
        karl_wave = karl_data[:, 0]
        karl_curve = karl_data[:, 1]

    correction = interpolate.interp1d( karl_wave, karl_curve)

    return correction


def deredden_spectra(wave, coords):
    """
    Apply S&F 2011 Extinction from SFD Map
    https://iopscience.iop.org/article/10.1088/0004-637X/737/2/103#apj398709t6

    Paramters
    ---------
    wave      array
        wavelength to apply correction
    coords    SkyCoord object
        sky coordinates
    """
    Rv = 3.1
    corr_SF2011 = 2.742  # Landolt V
    
    sfd = SFDQuery()
    ebv = sfd(coords)
    Av = corr_SF2011 * ebv
    ext = extinction.fitzpatrick99(np.array(wave, dtype=np.double), Av, Rv)

    deredden = 10**(0.4*np.array(ext))

    return deredden
