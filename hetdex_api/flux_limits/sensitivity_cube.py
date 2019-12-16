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
from numpy import any as nany
from scipy.interpolate import interp1d
import astropy.io.fits as fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord


class WavelengthException(Exception):
    """ Exception to raise when suspicious wavelength is passed"""
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
        the scaling to apply to the inverse of the cube values (Optional,
        default 1e-17)
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

    return f50s, header


class SensitivityCube(object):
    """
    Deals with flux limit cubes

    Parameters
    ----------
    f50vals : array
        3D datacube of datascale/flux_limit
    wcs : astropy.wcs:WCS
        world coordinate system to convert between ra, dec, lambda
        and pixel
    header : dict
        a dictionary of the headervalues to be stored in a
        FITS file
    wavelengths, alphas : array
        arrays of the wavelength in
        Angstrom and the alpha parameter
        of the Fleming+ 1995 function
    aper_corr : float
        Aperture correction to multiply
        the cubes with. If None, read
        from header. If not in header
        and aper_corr=None do nothing

    Attributes
    ----------
    f50vals : array
        50% 'limiting flux' to be used in the Fleming function, stored as
        lambda, dec, ra

    alpha_func : callable
        returns the Fleming alpha
        for an input wavelength

    """

    def __init__(self, f50vals, header, wavelengths, alphas, aper_corr=None):

        self.f50vals = f50vals

        # Fix issue with header
        if not "CD3_3" in header:
            header["CD3_3"] = header["CDELT3"]
            header["CD3_1"] = 0.0
            header["CD3_2"] = 0.0
            header["CD2_3"] = 0.0
            header["CD1_3"] = 0.0

        self.wcs = WCS(header)
        self.header = header

        # Deal with aperture corrections
        if aper_corr:
            self.aper_corr = aper_corr
        elif "APCOR" in self.header:
            self.aper_corr = self.header["APCOR"]
        else:
            self.aper_corr = 1.0        

        self.f50vals = self.f50vals*self.aper_corr

        self.alphas = alphas
        self.wavelengths = wavelengths
        self.alpha_func = interp1d(wavelengths, alphas, 
                                   fill_value="extrapolate")


    @classmethod
    def from_file(cls, fn_sensitivity_cube, wavelengths, alphas, 
                 datascale=1e-17, **kwargs):

        """
        Read in a sensitivity cube
        from a file
        
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
        datascale : float (optional)
            the values stored are 
            this_value/flim
        **kwargs :
            these are passed to the SensitivityCube init
        """

        f50vals, header = read_cube(fn_sensitivity_cube, datascale=datascale)

        return SensitivityCube(f50vals, header, wavelengths, alphas, **kwargs)

    def radecwltoxyz(self, ra, dec, lambda_):
        """
        Convert ra, dec position to
        x,y coordinate of cube

        Parameters
        ----------
        ra, dec : arrays
            right ascension &
            declination of source
        lambda_ : array
            wavelength in Angstrom

        Returns
        -------
        ix,iy,iz : arrays of int
            indices of arrays for datacube
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
            wavelength in Angstroms

        Returns
        -------
        f50s : array
            flux limits. If outside
            of cube return 999

        Raises
        ------
        WavelengthException :
            Annoys user if they pass
            wavelength outside of
            VIRUS range

        """

        ix, iy, iz = self.radecwltoxyz(ra, dec, lambda_)

        # Check for stuff outside of cube
        bad_vals = (ix >= self.f50vals.shape[2]) | (ix < 0) 
        bad_vals = bad_vals | (iy >= self.f50vals.shape[1]) | (iy < 0) 
        bad_vals = bad_vals | (iz >= self.f50vals.shape[0]) | (iz < 0) 
 
        ix[(ix >= self.f50vals.shape[2]) | (ix < 0)] = 0
        iy[(iy >= self.f50vals.shape[1]) | (iy < 0)] = 0
        iz[(iz >= self.f50vals.shape[0]) | (iz < 0)] = 0

        f50s = self.f50vals[iz, iy, ix]

        # Support arrays and floats
        try:
            f50s[bad_vals] = 999.0
        except TypeError:
            if nany(bad_vals):
                f50s = 999.0

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

        try:
            if lambda_[0] < 3000.0 or lambda_[0] > 6000.0:

                raise WavelengthException("""Odd wavelength value. Are you
                                             sure it's in Angstrom?""")
        except TypeError as e:
             if lambda_ < 3000.0 or lambda_ > 6000.0:

                raise WavelengthException("""Odd wavelength value. Are you
                                             sure it's in Angstrom?""")

        f50s = self.get_f50(ra, dec, lambda_)
        alphas = self.alpha_func(lambda_)

        return fleming_function(flux, f50s, alphas)

    def write(self, filename, datascale=1e-17, **kwargs):
        """
        Write the sensitivity cube to a FITS file. If any 
        aperture correction was applied, this is removed
        such that the saved data file should be identical 
        to the input (within numerical accuracy).

        Parameters
        ----------
        filename : str 
            Filename to write to
        datascale : float
           the scaling to apply to the
           inverse of the cube values 
           (Optional, default 1e-17)
        **kwargs :
            passed to the astropy.io.fits:writeto
            function
        """
        fits.writeto(filename, self.aper_corr*datascale/self.f50vals,
                     header=self.header, **kwargs)


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
    parser.add_argument("ra",
                        type=str,
                        help="RA of location to plot (HHhMMmSSs)")
    parser.add_argument("dec",
                        type=str,
                        help="DEC of location to plot (DDdMMmSSs)")
    parser.add_argument("lambda_", type=float, help="Wavelength to plot (A)")
    parser.add_argument("alpha", type=float, help="Alpha for Fleming")
    parser.add_argument("--fout", type=str, help="Filename to output to",
                        default=None)
    opts = parser.parse_args(args=args)

    coord = SkyCoord(opts.ra, opts.dec)

    scube = SensitivityCube.from_file(opts.filename, [3500.0, 5500.0], [opts.alpha, opts.alpha])
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
    scube = SensitivityCube.from_file(opts.filename, [3500.0, 5500.0], [-3.1, -3.1])

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
