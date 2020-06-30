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
import matplotlib.pyplot as plt
from numpy import (rint, array, around, multiply, isnan, meshgrid, mean,
                   median, sqrt, divide, linspace, ones, log10, loadtxt, polyval)
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
    sigmas : array
        3D datacube of datascale/noise
        where noise is the noise on
        a point source detection
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
    nsigma : float
        If the cubes don't contain
        1 sigma noise (e.g. in the HDR1
        cubes it's 6 sigma) specify it here
    conversion_poly : array
        this polynomial is used to 
        convert between the cut in
        detection significance and
        50% completeness values. If
        you specify it (as an array to
        feed into numpy.ployval) this
        value is used. If None try to
        read this info from the header
        if that fails default to assuming
        detection significance = Signal/Noise.

    Attributes
    ----------
    sigmas : array
        an array of the noise values

    alpha_func : callable
        returns the Fleming alpha
        for an input wavelength

    """
    # XXX might need a better name than nsigma as this is already used elsewhere
    def __init__(self, sigmas, header, wavelengths, alphas, aper_corr=None, 
                 nsigma=1.0, conversion_poly=None):

        self.sigmas = sigmas/nsigma
        self.nsigma = nsigma

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
        elif "APCOR0" in self.header: 
            self.aper_corr = self.header["APCOR0"]
            # XXX HACK HACK HACK
            # self.aper_corr = 1.0
        else:
            self.aper_corr = 1.0        
 
        if conversion_poly:
            self.conversion_poly = conversion_poly
        elif "CONVPOLY" in header:
            raise Exception("Using the header not implemented")
        else:
            # Assume sn = Signal/Noise
            self.conversion_poly = [1.0, 0.0]

        self.sigmas = self.sigmas*self.aper_corr

        self.alphas = array(alphas)
        self.wavelengths = wavelengths

        # Depends if alphas depend on wavelength or
        # is specified per cube cell
        if len(self.alphas.shape) == 3:
            self.alpha_is_cube = True 
        else: 
            self.alpha_is_cube = False 
            self.alpha_func = interp1d(wavelengths, alphas, 
                                       fill_value="extrapolate")

    def get_alpha(self, ra, dec, lambda_):

        # If alpha is just an array versus wavelength
        # return the value here
        if not self.alpha_is_cube:
            return self.alpha_func(lambda_)

        # Alpha stored in a cube 
        ix, iy, iz = self.radecwltoxyz(ra, dec, lambda_)

        # Check for stuff outside of cube
        bad_vals = (ix >= self.alphas.shape[2]) | (ix < 0) 
        bad_vals = bad_vals | (iy >= self.alphas.shape[1]) | (iy < 0) 
        bad_vals = bad_vals | (iz >= self.alphas.shape[0]) | (iz < 0) 
 
        ix[(ix >= self.alphas.shape[2]) | (ix < 0)] = 0
        iy[(iy >= self.alphas.shape[1]) | (iy < 0)] = 0
        iz[(iz >= self.alphas.shape[0]) | (iz < 0)] = 0

        alphas_here = self.alphas[iz, iy, ix]

        # Support arrays and floats
        try:
            alphas_here[bad_vals] = 999.0
        except TypeError:
            if isnan(bad_vals):
                aphas_here = 999.0

        return alphas_here


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

        sigmas, header = read_cube(fn_sensitivity_cube, datascale=datascale)

        return SensitivityCube(sigmas, header, wavelengths, alphas, **kwargs)


    def apply_flux_recalibration(self, rescale, flux_calib_correction_file=None):
        """
        Apply a recalibration of the fluxes to the 
        cube

        Parameters
        ----------
        rescale : float 
           value to multiply the flux limit cubes
           to rescale

        flux_calib_correction_file : str (optional)
           filename containing a polynomial
           fit (HETDEX - TRUTH)/HETDEX versus
           wavelength to correct for 
           problems with the flux
           calibration. Should be a polynomial
           centered on 4600, i.e. input to
           polyval(pvals, wl - 4600.0)
        """

        if flux_calib_correction_file:
            pvals = loadtxt(flux_calib_correction_file)

        for iz in range(self.sigmas.shape[0]):
            ra, dec, wl = self.wcs.wcs_pix2world(0, 0, iz, 0)

            if wl < 3850.0:
                wl = 3850.0

            if flux_calib_correction_file:
                self.sigmas[iz, :, :] = rescale*self.sigmas[iz, :, :]*(1.0 - polyval(pvals, wl - 4600.0)) 
            else: 
                self.sigmas[iz, :, :] = rescale*self.sigmas[iz, :, :]  

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

    def get_f50(self, ra, dec, lambda_, sncut):
        """
        Get 50% completeness flux from the cube at
        ra, dec, lambda

        Parameters
        ----------
        ra, dec : array
            right ascension and dec in degrees
        lambda_ : array
            wavelength in Angstroms
        sncut : float
            cut in detection significance 
            that defines this catalogue

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
        bad_vals = (ix >= self.sigmas.shape[2]) | (ix < 0) 
        bad_vals = bad_vals | (iy >= self.sigmas.shape[1]) | (iy < 0) 
        bad_vals = bad_vals | (iz >= self.sigmas.shape[0]) | (iz < 0) 
 
        ix[(ix >= self.sigmas.shape[2]) | (ix < 0)] = 0
        iy[(iy >= self.sigmas.shape[1]) | (iy < 0)] = 0
        iz[(iz >= self.sigmas.shape[0]) | (iz < 0)] = 0

        f50s = polyval(self.conversion_poly, sncut*self.sigmas[iz, iy, ix])

        # Support arrays and floats
        try:
            f50s[bad_vals] = 999.0
        except TypeError:
            if isnan(bad_vals):
                f50s = 999.0

        return f50s

    def return_completeness_hdr1(self, flux, ra, dec, lambda_):
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

        # hard coded sncut for HDR1
        f50s = self.get_f50(ra, dec, lambda_, 6.0)
        alphas = self.get_alpha(ra, dec, lambda_)

        return fleming_function(flux, f50s, alphas)

    def compute_snr(self, flux, ra, dec, lambda_):
        """
        Compute the flux divided by the noise for 
        a given source. This differs from the
        detection significance in that there
        is no PSF model weighting

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
        snr : array
            signal divided by noise


        """
        ix, iy, iz = self.radecwltoxyz(ra, dec, lambda_)

        # Check for stuff outside of cube
        bad_vals = (ix >= self.sigmas.shape[2]) | (ix < 0) 
        bad_vals = bad_vals | (iy >= self.sigmas.shape[1]) | (iy < 0) 
        bad_vals = bad_vals | (iz >= self.sigmas.shape[0]) | (iz < 0) 

        ix[(ix >= self.sigmas.shape[2]) | (ix < 0)] = 0
        iy[(iy >= self.sigmas.shape[1]) | (iy < 0)] = 0
        iz[(iz >= self.sigmas.shape[0]) | (iz < 0)] = 0

        noise = self.sigmas[iz, iy, ix]
        snr = flux/noise

        # Support arrays and floats
        try:
            snr[bad_vals] = 0.0
        except TypeError:
            if isnan(snr):
                snr = 0.0

        return snr


    def return_completeness(self, flux, ra, dec, lambda_, sncut):
        """
        Return completeness at a 3D position. This computes stuff
        in terms of SNR in the Fleming (1995) function but
        it's actually equivalent to computing the 50% flux
        and computing everthing in terms of flux

        Parameters
        ----------
        flux : array
            fluxes of objects
        ra, dec : array
            right ascension and dec in degrees
        lambda_ : array
            wavelength in Angstrom
        sncut : float
            the detection significance (S/N) cut
            applied to the data

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
        

        snr = self.compute_snr(flux, ra, dec, lambda_)
        snr50 = polyval(self.conversion_poly, sncut)
        alphas = self.get_alpha(ra, dec, lambda_)

        return fleming_function(snr, snr50, alphas)


    def return_wlslice_completeness(self, flux, lambda_low, lambda_high, 
                                    sncut, noise_cut=1e-12):
        """
        Return completeness of a wavelength slice. This computes 
        stuff in terms of SNR in the Fleming (1995) function but
        it's actually equivalent to computing the 50% flux
        and computing everthing in terms of flux

        Parameters
        ----------
        flux : array
            fluxes of objects
        lambda_low, lambda_high : float
            wavelength slice in Angstrom
            (includes these slices)
        sncut : float
            the detection significance (S/N) cut
            applied to the data
        noise_cut : float
            remove areas with more noise
            than this
 
        Return
        ------
        fracdet : array
            fraction detected in this slice

        """

        if lambda_low < 3000.0 or lambda_low > 6000.0:
            raise WavelengthException("""Odd wavelength value. Are you
                                         sure it's in Angstrom?""")

        snr50 = polyval(self.conversion_poly, sncut)
        ix, iy, izlo = self.radecwltoxyz(self.wcs.wcs.crval[0], self.wcs.wcs.crval[1], lambda_low)
        ix, iy, izhigh = self.radecwltoxyz(self.wcs.wcs.crval[0], self.wcs.wcs.crval[1], lambda_high)
        noise = self.sigmas[izlo:(izhigh + 1), :, :]

        if len(self.alphas.shape) > 1:
            alphas = self.alphas[izlo:(izhigh + 1), :, :] 
        else:
            # rough approximation to lambda varying across window
            alphas = self.alpha_func(0.5*(lambda_low + lambda_high))

        compls = []
        for f in flux:
            snr = f/noise
            compl = fleming_function(snr, snr50, alphas)       

            # works so long as pixels equal area
            compls.append(mean(compl[noise < noise_cut]))

        return array(compls)

    def return_wlslice_f50(self, lambda_low, lambda_high, 
                           sncut, noise_cut=1e-12):
        """
        Return 50% completeness of a wavelength slice.  

        Parameters
        ----------
        lambda_low, lambda_high : float
            wavelength slice in Angstrom
            (includes these slices)
        sncut : float
            the detection significance (S/N) cut
            applied to the data
        noise_cut : float
            remove areas with more noise
            than this
 
        Return
        ------
        fracdet : array
            fraction detected in this slice

        """

        try:
            if lambda_low < 3000.0 or lambda_low > 6000.0:
                raise WavelengthException("""Odd wavelength value. Are you
                                             sure it's in Angstrom?""")
        except ValueError:
            if any(lambda_low < 3000.0) or any(lambda_low > 6000.0):
                raise WavelengthException("""Odd wavelength value. Are you
                                             sure it's in Angstrom?""")
 
        snr50 = polyval(self.conversion_poly, sncut)
        ix, iy, izlo = self.radecwltoxyz(self.wcs.wcs.crval[0], self.wcs.wcs.crval[1], lambda_low)
        ix, iy, izhigh = self.radecwltoxyz(self.wcs.wcs.crval[0], self.wcs.wcs.crval[1], lambda_high)
        noise = self.sigmas[izlo:(izhigh + 1), :, :]
        noise = noise[noise < noise_cut]

        return sncut*median(noise) 

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

        fits.writeto(filename, self.aper_corr*datascale/self.sigmas,
                     header=self.header, **kwargs)

def plot_completeness(args=None):
    """
    Plot the completeness curve at
    a particular place
    """
    import matplotlib as mpl
    mpl.use("agg")
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


def plot_slice_of_cube(axes, scube, ramin, ramax, decmin, decmax, 
                       wl, sncut, cmap="gist_rainbow", n=100, cax=None):
    """
    Plot a slice of a sensitivity cube

    Parameters
    ----------
    axes : matplotlib.pyplot.axes
        axes to put plot on
    scube : SensitivityCube
        cube to plot
    wl : float
        wavelength slice to plot
    n : int
        resolution
    """
    rar = linspace(ramin, ramax, n)
    decr = linspace(decmin, decmax, n)

    ras, decs = meshgrid(rar, decr)
    wls = wl*ones(len(ras))

    flims = scube.get_f50(ras, decs, wls, sncut)
    flims = 1e17*flims
    img = flims.reshape(n,n)

    im = axes.imshow(img, extent=[ramin, ramax, decmin, decmax], origin="lower left", 
                     aspect="auto", cmap=cmap, vmin=0.0, vmax=30)
    plt.colorbar(im, cax=cax, label="flux limit [10$^{-17}$ erg/s/cm$^{2}$]", pad=0.0)
    axes.set_xlabel("RA (deg)")
    axes.set_ylabel("Dec (deg)")


def plot_completeness_versus_wl(args=None):
    """
    Plot the completeness curve versus wavelength
    at a particular ra, dec
    """
    import matplotlib as mpl
    mpl.use("agg") 
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
