"""
Sensitivity Cube Reader

Read in Karl's sensitivity cubes and 
produce expected detection fraction
from the Fleming (Fleming+ 1995) parameterisation

References:

Fleming et al. 1995 
http://adsabs.harvard.edu/abs/1995AJ....109.1044F


See also:

Donghui Jeong's explanation 
https://luna.mpe.mpg.de/wikihetdex/index.php/Flim_files_and_Fleming_curve

.. moduleauthor:: Daniel Farrow <dfarrow@mpe.mpg.de>

"""
from __future__ import (absolute_import, print_function)
from numpy import (rint, array, around, multiply, 
                   sqrt, divide, linspace, ones, log10)
from scipy.interpolate import interp1d
import astropy.io.fits as fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

class WavelengthException(Exception):
    """ Exception to raise when suspicious 
        wavelength is passed"""
    pass

def fleming_function(flux, flim, alpha):
    """
    Implementation of the Fleming+ 1995
    completeness function

    Parameters
    ----------
    flux : float
        flux of source
    flim : float
        flux at 50% completeness
    alpha : float
        parameter controlling
        width of fall off

    Returns
    -------
    val : float
        Fleming+ 1995 function at this
        flux


    """

    fdiff = multiply(alpha, (-2.5*log10(flux) + 2.5*log10(flim)))

    return 0.5*(1.0 + divide(fdiff, sqrt(1.0 + fdiff*fdiff)))


def read_cube(fn, datascale=1e-17):
   """
   Read a Karl sensitivity cube and
   return flux limits and WCS for
   ra, dec conversions

   Parameters
   ----------
   fn : str 
       filename of cube to read

   datascale : float
       the scaling to apply to the
       inverse of the cube values 
       (Optional, default 1e-17)

   Return
   ------
   f50s : array
       the datacube of the flux limits
   WCS : astropy.wcs:WCS 
       the WCS to do x,y,z to ra,dec,wl
       conversions
   """

   hdus = fits.open(fn)
   f50s = datascale/(hdus[0].data)
   header = hdus[0].header

   # Fix issue with header
   header["CD3_3"] = header["CDELT3"]
   header["CD3_1"] = 0.0
   header["CD3_2"] = 0.0
   header["CD2_3"] = 0.0
   header["CD1_3"] = 0.0

   wcs = WCS(hdus[0].header)

   return f50s, wcs 


class SensitivityCube(object):
    """
    Read in a sensitivity cube
    and Fleming parameters
    for an IFU

    Parameters
    ----------
    fn_sensitivity_cube : str
        the file name of a cube
        containing the limiting
        magnitude

    wavelengths, alphas : array
        arrays of the wavelength in
        Angstrom and the alpha parameter
        of the Fleming+ 1995 function

    Attributes
    ----------
    f50vals : array
        50% 'limiting flux' to
        be used in the Fleming function,
        stored as lambda_, dec_, ra_

    wcs : astropy.wcs.WCS
        world coordinate system for
        cube

    alpha_func : callable
        returns the Fleming alpha
        for an input wavelength

    """
    def __init__(self, fn_sensitivity_cube, wavelengths, alphas, 
                 datascale=1e-17):

        self.f50vals, self.wcs = read_cube(fn_sensitivity_cube, 
                                        datascale=datascale)

        self.alpha_func = interp1d(wavelengths, alphas, 
                                   fill_value="extrapolate")


    def radecwltoxyz(self, ra, dec, lambda_):
        """
        Convert ra, dec position to
        x,y coordinate of cube

        Parameter
        ---------
        ra, dec : arrays
            right ascension &
            declination of source
 
        lambda_ : array
            wavelength in Angstrom

       
        Returns
        -------
        ix,iy,iz : arrays of int
            indices of arrays for 
            datacube

        """

        lambda_ = array(lambda_)
        ix,iy,iz = self.wcs.wcs_world2pix(ra, dec, lambda_, 0)

        return array(around(ix), dtype=int), array(around(iy), dtype=int), \
               array(around(iz), dtype=int)

    def get_f50(self, ra, dec, lambda_):
        """
        Get flim from cube at
        ra, dec, lambda

        Parameters
        ----------
        ra, dec : array
            right ascension and dec in degrees
        lambda_ : array
            wavelength in Angstrom

        Raises
        ------
        WavelengthException :
            Annoys user if they pass
            wavelength outside of
            VIRUS range

        """

        ix, iy, iz = self.radecwltoxyz(ra, dec, lambda_)
        f50s = self.f50vals[iz, iy, ix]

        return f50s

    def return_completeness(self, flux, ra, dec, lambda_):
        """
        Return completeness at a 3D position

        Parameters
        ----------
        flux : array
            fluxes of objects
        ra, dec : array
            right ascension and dec in degrees
        lambda_ : array
            wavelength in Angstrom

        Return
        ------
        fracdet : array
            fraction detected 

        """
 
        if lambda_[0] < 3000.0 or lambda_[0] > 6000.0:

            raise WavelengthException("""Odd wavelength value. Are you
                                         sure it's in Angstrom?""")

        f50s = self.get_f50(ra, dec, lambda_)
        alphas = self.alpha_func(lambda_)

        return fleming_function(flux, f50s, alphas)


def plot_completeness(args=None):
    """
    Plot the completeness curve at
    a particular place
    """
 
    import matplotlib.pyplot as plt
    import argparse
    parser = argparse.ArgumentParser(description="""Plot the Fleming fit to
                                                    completeness""")
    parser.add_argument("filename", type=str)
    parser.add_argument("ra", type=str, help="RA of location to plot (HHhMMmSSs)")
    parser.add_argument("dec", type=str, help="DEC of location to plot (DDdMMmSSs)")
    parser.add_argument("lambda_", type=float, help="Wavelength to plot (A)")
    parser.add_argument("alpha", type=float, help="Alpha for Fleming")
    parser.add_argument("--fout", type=str, help="Filename to output to", 
                        default=None)
    opts = parser.parse_args(args=args)

    coord = SkyCoord(opts.ra, opts.dec)

    scube = SensitivityCube(opts.filename, [3500.0, 5500.0], [opts.alpha, opts.alpha])
    f50 = scube.get_f50([coord.ra.deg], [coord.dec.deg], [opts.lambda_])
    print(f50)

    fluxes = linspace(0.01*f50, 4*f50, 100)
 
    compl = scube.return_completeness(array(fluxes), coord.ra.deg*ones(len(fluxes)), 
                                      coord.dec.deg*ones(len(fluxes)), 
                                      opts.lambda_*ones(len(fluxes)))
 
    plt.plot(fluxes/1e-17, compl, "r-")
    plt.xlabel("Flux $10^{-17}\,$(erg/s/cm$^2$/A)", fontsize=14.0)
    plt.ylabel("Completeness", fontsize=14.0)

    if opts.fout:
        plt.savefig(opts.fout)
    else:
        plt.show()


def plot_completeness_versus_wl(args=None):
    """
    Plot the completeness curve versus wavelength
    at a particular ra, dec
    """
 
    import matplotlib.pyplot as plt
    import argparse
    parser = argparse.ArgumentParser(description="""Plot the Fleming fit to
                                                    completeness""")
    parser.add_argument("filename", type=str)
    parser.add_argument("ra", type=str, help="RA of location to plot (HHhMMmSSs)")
    parser.add_argument("dec", type=str, help="DEC of location to plot (DDdMMmSSs)")
    #parser.add_argument("alphas", type=float, help="Alpha for Fleming")
    parser.add_argument("--fout", type=str, help="Filename to output to", 
                        default=None)
    opts = parser.parse_args(args=args)


    coord = SkyCoord(opts.ra, opts.dec)
    scube = SensitivityCube(opts.filename, [3500.0, 5500.0], [-3.1, -3.1])

    wls = linspace(3500, 5490.0, 1000)
    f50 = scube.get_f50(coord.ra.deg*ones(len(wls)), coord.dec.deg*ones(len(wls)), wls)

   # fluxes = linspace(0.01*f50, 4*f50, 100)
   # compl = scube.return_completeness(array(fluxes), coord.ra.deg*ones(len(fluxes)), 
   #                                   coord.dec.deg*ones(len(fluxes)), 
   #                                   opts.lambda_*ones(len(fluxes)))
 
    plt.plot(wls, f50/1e-16, "k-", label="Flux at 50% completeness")
    plt.ylabel("Flux $10^{-16}\,$(erg/s/cm$^2$/A)", fontsize=14.0)
    plt.xlabel("Wavelength (A)", fontsize=14.0)
    plt.legend(loc="upper right")
 
    if opts.fout:
        plt.savefig(opts.fout)
    else:
        plt.show()
 
