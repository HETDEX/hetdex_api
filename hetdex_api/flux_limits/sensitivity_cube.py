"""
Sensitivity Cube Reader

Read in Karl's sensitivity cubes and 
produce expected detection fractions
from the Fleming (Fleming et al 1995) 
parameterisation

References
----------

- Fleming et al. 1995 

  Fleming D.~E.~B., Harris W.~E., Pritchet C.~J., Hanes D.~A., 1995, 
  AJ, 109, 1044. doi:10.1086/117340

  http://adsabs.harvard.edu/abs/1995AJ....109.1044F

- Donghui Jeong's explanation (internal to HETDEX)
  https://luna.mpe.mpg.de/wikihetdex/index.php/Flim_files_and_Fleming_curve

.. moduleauthor:: Daniel Farrow <dfarrow@mpe.mpg.de>
"""

from __future__ import (absolute_import, print_function)
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from numpy import (rint, array, around, multiply, isnan, meshgrid, mean, isfinite,
                   median, sqrt, divide, linspace, ones, log10, loadtxt, polyval, inf,
                   repeat, newaxis, logical_not, arange, tile)
from numpy.ma import filled
from numpy.ma import array as maskedarray
from numpy.random import normal
from numpy import any as nany
from scipy.interpolate import interp1d
import astropy.io.fits as fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from hetdex_api.flux_limits import flim_models, flim_model_cache
from hetdex_api.config import HDRconfig

class WavelengthException(Exception):
    """ Exception to raise when suspicious wavelength is passed"""
    pass


def fleming_function(flux, flim, alpha):
    """
    Implementation of the Fleming et al 1995
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
    header : dict
        the header of the cube's
        FITS file
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
    header : dict
        a dictionary of the headervalues to be stored in a
        FITS file
    wavelengths, alphas : array
        arrays of the wavelength in
        Angstrom and the alpha parameter
        of the Fleming+ 1995 function
    aper_corr : float (Optional)
        Aperture correction to multiply
        the cubes with. If None, read
        from header. If not in header
        and aper_corr=None do nothing.
        Default is 1.0.
    nsigma : float
        If the cubes don't contain
        1 sigma noise (e.g. in the HDR1
        cubes it's 6 sigma) specify it here
    flim_model : str
        the name of the flux limit model 
        to use. This can either be 'hdr1',
        'hdr2pt1', hdr2pt1pt1'
    mask : array (optional)
        a spatial ra, dec mask with the same
        WCS and dimensions as the data (default:
        None)

    cache_sim_interp : bool (optional)
        cache the SimulationInterpolator,
        so if you use another SensitivityCube
        it will use the same model from before
        (hdr2pt1pt1 or later only, only works
         if you don't change flim_model)



    Attributes
    ----------
    sigmas : array
        an array of the noise values
    alpha_func : callable
        returns the Fleming alpha
        for an input wavelength
    wcs : astropy.wcs:WCS
        world coordinate system to convert between ra, dec, lambda
        and pixel
    f50_from_noise : callable
        function that converts the values
        in `sigmas` to flux values at 
        50% completeness

    """
    def __init__(self, sigmas, header, wavelengths, alphas, aper_corr=1.0, 
                 nsigma=1.0, flim_model="hdr2pt1pt3", mask=None, 
                 cache_sim_interp = True): 

        # Note: flux limit model is also passed here by the HDF5 container class
        #print("Flux limit model: ", flim_model) 

        if type(mask) != type(None):
            mask = logical_not(mask)
            mask3d = repeat(mask[newaxis, :, :], sigmas.shape[0], axis=0)
            self.sigmas = maskedarray(sigmas/nsigma, mask=mask3d, fill_value=999.0)
        else:
            self.sigmas = maskedarray(sigmas/nsigma, fill_value=999.0)

        # collapse the data to create a continuum mask
        self.collapsed_data = filled(self.sigmas, 0).sum(axis=0)

        self.nsigma = nsigma

        # Decide between Fleming function or interpolation
        if (flim_model == "hdr1") or (flim_model == "hdr2pt1"):
            self.sinterp = None	
        elif flim_model_cache.cached_model == flim_model and cache_sim_interp:
            self.sinterp = flim_model_cache.cached_sim_interp
        else:
            conf = HDRconfig()
            fdir = conf.flim_sim_completeness 
            if "snbased" in flim_model:
                print("Using S/N based model!") 
                self.sinterp = flim_models.SimulationInterpolator(fdir, snmode=True)
            else:
                self.sinterp = flim_models.SimulationInterpolator(fdir, snmode=False)

        # cache the model to save reinitialising the object
        if cache_sim_interp: 
            flim_model_cache.cached_sim_interp = self.sinterp
            flim_model_cache.cached_model = flim_model 

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
        else:
            self.aper_corr = 1.0

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

        # Select flux limit model
        try:
            self.f50_from_noise = getattr(flim_models, "{:s}_f50_from_noise".format(flim_model)) 
        except AttributeError as e:
            print("Chosen flux limit model not found!")
            raise e


    def get_alpha(self, ra, dec, lambda_):
        """
        Return the parameter controlling
        the slope of the Fleming+ (1995) function

        """

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
        Convert ra, dec, wavelength position to
        x,y, z coordinate of cube

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

    def get_average_f50(self, ra, dec, lambda_, sncut, npix=1):
        """
        Get the maximum 50% completeness flux from the cube in
        an npix box around and ra, dec, lambda

        Parameters
        ----------
        ra, dec : array
            right ascension and dec in degrees
        lambda_ : array
            wavelength in Angstroms
        sncut : float
            cut in detection significance 
            that defines this catalogue
        npix : int
            the box will be 2*npix + 1 on
            a side, i.e. number of pixels
            around the position to 
            consider.
 
        Returns
        -------
        f50s : array
            max flux limits in cubes. If outside
            of cube return 999
        """

        ixc, iyc, izc = self.radecwltoxyz(ra, dec, lambda_)
  
        na = int(2*npix + 1)
 
        # [x1-1, x1, x1+1, x2-1, x2, x2+1, .....]
        offsets = arange(-1.0*npix, npix + 1, 1, dtype=int)
        ix = ixc.repeat(na) + tile(offsets, len(ixc))

        # same x for all x, y in loop
        ix = ix.repeat(na*na)

        iy = iyc.repeat(na) + tile(offsets, len(iyc))

        # same y for all z values in loop
        iy = iy.repeat(na)

        # tile full y-loop for each x-value
        iy = tile(iy.reshape(len(iyc), na*na), na)
        iy = iy.flatten()

        # [z1-1, z1, z1+1, z2-1, z2, z2+1, .....]
        iz = izc.repeat(len(offsets)) + tile(offsets, len(izc))

        # z axis fastest repeating, tile z loop for every x and y value
        iz = tile(iz.reshape(len(izc), na), na*na)
        iz = iz.flatten()

        # Check for stuff outside of cube
        bad_vals = (ix >= self.sigmas.shape[2]) | (ix < 0) 
        bad_vals = bad_vals | (iy >= self.sigmas.shape[1]) | (iy < 0) 
        bad_vals = bad_vals | (iz >= self.sigmas.shape[0]) | (iz < 0) 
 
        ix[(ix >= self.sigmas.shape[2]) | (ix < 0)] = 0
        iy[(iy >= self.sigmas.shape[1]) | (iy < 0)] = 0
        iz[(iz >= self.sigmas.shape[0]) | (iz < 0)] = 0

        f50s = self.f50_from_noise(self.sigmas.filled()[iz, iy, ix], sncut)

        # Support arrays and floats
        f50s[bad_vals] = 999.0

        #print(ix)
        #print(iy)
        #print(iz)

        # return the average flim in the area
        f50s = f50s*f50s
        return sqrt(f50s.reshape(len(ra), na*na*na).mean(axis=1))

    def get_collapsed_value(self, ra, dec):

        ix, iy, iz = self.radecwltoxyz(ra, dec, 4500.)

        # Check for stuff outside of cube
        bad_vals = (ix >= self.sigmas.shape[2]) | (ix < 0) 
        bad_vals = bad_vals | (iy >= self.sigmas.shape[1]) | (iy < 0) 
 
        ix[(ix >= self.sigmas.shape[2]) | (ix < 0)] = 0
        iy[(iy >= self.sigmas.shape[1]) | (iy < 0)] = 0

        noise = self.collapsed_data[iy, ix]
        noise[bad_vals] = 999.0

        return noise
   

    def get_local_max_f50(self, ra, dec, lambda_, sncut, npix=1):
        """
        Get the maximum 50% completeness flux from the cube in
        an npix box around and ra, dec, lambda

        Parameters
        ----------
        ra, dec : array
            right ascension and dec in degrees
        lambda_ : array
            wavelength in Angstroms
        sncut : float
            cut in detection significance 
            that defines this catalogue
        npix : int
            the box will be 2*npix + 1 on
            a side, i.e. number of pixels
            around the position to 
            consider.
 
        Returns
        -------
        f50s : array
            max flux limits in cubes. If outside
            of cube return 999
        """

        ixc, iyc, izc = self.radecwltoxyz(ra, dec, lambda_)
  
        na = int(2*npix + 1)
 
        # [x1-1, x1, x1+1, x2-1, x2, x2+1, .....]
        offsets = arange(-1.0*npix, npix + 1, 1, dtype=int)
        ix = ixc.repeat(na) + tile(offsets, len(ixc))

        # same x for all x, y in loop
        ix = ix.repeat(na*na)

        iy = iyc.repeat(na) + tile(offsets, len(iyc))

        # same y for all z values in loop
        iy = iy.repeat(na)

        # tile full y-loop for each x-value
        iy = tile(iy.reshape(len(iyc), na*na), na)
        iy = iy.flatten()

        # [z1-1, z1, z1+1, z2-1, z2, z2+1, .....]
        iz = izc.repeat(len(offsets)) + tile(offsets, len(izc))

        # z axis fastest repeating, tile z loop for every x and y value
        iz = tile(iz.reshape(len(izc), na), na*na)
        iz = iz.flatten()

        # Check for stuff outside of cube
        bad_vals = (ix >= self.sigmas.shape[2]) | (ix < 0) 
        bad_vals = bad_vals | (iy >= self.sigmas.shape[1]) | (iy < 0) 
        bad_vals = bad_vals | (iz >= self.sigmas.shape[0]) | (iz < 0) 
 
        ix[(ix >= self.sigmas.shape[2]) | (ix < 0)] = 0
        iy[(iy >= self.sigmas.shape[1]) | (iy < 0)] = 0
        iz[(iz >= self.sigmas.shape[0]) | (iz < 0)] = 0

        f50s = self.f50_from_noise(self.sigmas.filled()[iz, iy, ix], sncut)

        # Support arrays and floats
        f50s[bad_vals] = 999.0

        #print(ix)
        #print(iy)
        #print(iz)

        # return the max value in area
        return f50s.reshape(len(ra), na*na*na).max(axis=1) 

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

        """

        ix, iy, iz = self.radecwltoxyz(ra, dec, lambda_)

        # Check for stuff outside of cube
        bad_vals = (ix >= self.sigmas.shape[2]) | (ix < 0) 
        bad_vals = bad_vals | (iy >= self.sigmas.shape[1]) | (iy < 0) 
        bad_vals = bad_vals | (iz >= self.sigmas.shape[0]) | (iz < 0) 
 
        ix[(ix >= self.sigmas.shape[2]) | (ix < 0)] = 0
        iy[(iy >= self.sigmas.shape[1]) | (iy < 0)] = 0
        iz[(iz >= self.sigmas.shape[0]) | (iz < 0)] = 0

        f50s = self.f50_from_noise(self.sigmas.filled()[iz, iy, ix], sncut)

        # Support arrays and floats
        try:
            f50s[bad_vals] = 999.0
        except TypeError:
            if bad_vals:
                f50s = 999.0

        return f50s

    def compute_snr(self, flux, ra, dec, lambda_):
        """
        Compute the flux divided by the noise for 
        a given source. 

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

        noise = self.sigmas.filled()[iz, iy, ix]
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
        Return completeness at a 3D position as an array. 
        If for whatever reason the completeness is NaN, it's
        replaced by 0.0. 

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

        Raises
        ------
        WavelengthException :
            Annoys user if they pass
            wavelength outside of
            VIRUS range
        """

        try: 
            if lambda_[0] < 3000.0 or lambda_[0] > 6000.0:

                raise WavelengthException("""Odd wavelength value. Are you
                                             sure it's in Angstrom?""")
        except TypeError as e:
             if lambda_ < 3000.0 or lambda_ > 6000.0:

                raise WavelengthException("""Odd wavelength value. Are you
                                             sure it's in Angstrom?""")
        

        f50s = self.get_f50(ra, dec, lambda_, sncut)

        if self.sinterp:
            # interpolate over the simulation
            fracdet = self.sinterp(flux, f50s, lambda_, sncut)
        else:
            alphas = self.get_alpha(ra, dec, lambda_)
            fracdet = fleming_function(flux, f50s, alphas)

        try:
            fracdet[isnan(fracdet)] = 0.0
        except TypeError:
            if isnan(fracdet):
                fracdet = 0.0

        return fracdet


    def return_wlslice_completeness(self, flux, lambda_low, lambda_high, 
                                    sncut, noise_cut=1e-15, pixlo=9, 
                                    pixhi=22, return_vals = False):
        """
        Return completeness of a wavelength slice. NaN completeness
        values are replaced with zeroes, noise values greater than
        noise cut or NaN noise values are simply excluded from the 
        mean 

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
            than this. Default: 1e-16 erg/s/cm2
        return_vals : bool (optional)
            if true alse return an array
            of the noise values
 
        Return
        ------
        fracdet : array
            fraction detected in this slice

        """

        if lambda_low < 3000.0 or lambda_low > 6000.0:
            raise WavelengthException("""Odd wavelength value. Are you
                                         sure it's in Angstrom?""")

        ix, iy, izlo = self.radecwltoxyz(self.wcs.wcs.crval[0], self.wcs.wcs.crval[1], lambda_low)
        ix, iy, izhigh = self.radecwltoxyz(self.wcs.wcs.crval[0], self.wcs.wcs.crval[1], lambda_high)

        if izlo < 0:
            print("Warning! Lower wavelength below range")
            izlo = 0        

        if izhigh > self.sigmas.shape[0] - 1:
            print("Warning! Upper wavelength above range")
            izhigh = self.sigmas.shape[0] - 1

        izlo = int(izlo)
        izhigh = int(izhigh)

        # remove pixel border and select wavelength slice
        noise = self.sigmas.filled()[izlo:(izhigh + 1), pixlo:pixhi, pixlo:pixhi]
 
        # Test what happens with fixed noise
        #noise = noise*0 + normal(loc=1e-17, scale=2e-18, 
        #                         size=noise.shape[0]*noise.shape[1]).reshape(noise.shape[0], noise.shape[1])

        # create a cube of the wavelengths
        r, d, wl_1d = self.wcs.wcs_pix2world(ones(1 + izhigh - izlo), ones(1 + izhigh - izlo), 
                                            range(izlo, izhigh + 1), 0)
        waves = wl_1d.repeat(noise.shape[1]*noise.shape[2])
 
        try:
            waves = waves.reshape(noise.shape) 
        except ValueError as e:
            print(noise.shape)
            print(len(wl_1d))
            print(izlo, izhigh)

        # remove masked data and bad data
        sel = (noise < noise_cut) & isfinite(noise) 
        noise = noise[sel] 
        waves = waves[sel]

        # Test for fixed lambda
        # waves = waves*0 + lambda_low

        if len(noise) == 0:
            if return_vals:
                return [], []
            else:
                return []

        f50s = self.f50_from_noise(noise, sncut)

        if type(self.sinterp) == type(None):
            if len(self.alphas.shape) > 1:
                alphas = self.alphas[izlo:(izhigh + 1), :, :] 
            else:
                # rough approximation to lambda varying across window
                alphas = self.alpha_func(0.5*(lambda_low + lambda_high))

        compls = []
        for f in flux:
            if self.sinterp:
                compl = self.sinterp(f, f50s.flatten(), waves.flatten(), sncut)
            else:
                compl = fleming_function(f, f50s, alphas)       

            compl[isnan(compl)] = 0.0

            # works so long as pixels equal area
            if len(compl) > 0:
                compls.append(mean(compl))
            else:
                compls.append(0.0)

        if return_vals:
            return array(compls), noise.flatten()
        else:
            return array(compls)

    def return_wlslice_f50(self, lambda_low, lambda_high, 
                           sncut, noise_cut=1e-16):
        """
        Return flux at 50% completeness of a wavelength slice.  

        Parameters
        ----------
        lambda_low, lambda_high : float
            wavelength slice in Angstrom
            (includes these slices)
        sncut : float
            the detection significance (S/N) cut
            applied to the data
        noise_cut : float (optional)
            remove areas with more noise
            than this. Default: 1e-16 erg/s/cm2
 
        Return
        ------
        f50 : float
            the flux at 50% completeness
            for the given ``sncut`` in this
            wavelength slice

        """

        try:
            if lambda_low < 3000.0 or lambda_low > 6000.0:
                raise WavelengthException("""Odd wavelength value. Are you
                                             sure it's in Angstrom?""")
        except ValueError:
            if any(lambda_low < 3000.0) or any(lambda_low > 6000.0):
                raise WavelengthException("""Odd wavelength value. Are you
                                             sure it's in Angstrom?""")
 
        ix, iy, izlo = self.radecwltoxyz(self.wcs.wcs.crval[0], self.wcs.wcs.crval[1], lambda_low)
        ix, iy, izhigh = self.radecwltoxyz(self.wcs.wcs.crval[0], self.wcs.wcs.crval[1], lambda_high)
        noise = self.sigmas.filled()[izlo:(izhigh + 1), :, :]
        noise = noise[(noise < noise_cut) & (noise > 0)]
  
        f50 = self.f50_from_noise(median(noise), sncut)

        return f50

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

        fits.writeto(filename, self.aper_corr*datascale/self.sigmas.data,
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
    parser.add_argument("--sncut", type=float, default=4.5, 
                        help="S/N cut used")
    parser.add_argument("--fout", type=str, help="Filename to output to",
                        default=None)
    opts = parser.parse_args(args=args)

    coord = SkyCoord(opts.ra, opts.dec)

    scube = SensitivityCube.from_file(opts.filename, [3500.0, 5500.0], [opts.alpha, opts.alpha])
    f50 = scube.get_f50([coord.ra.deg], [coord.dec.deg], [opts.lambda_], opts.sncut)

    fluxes = linspace(0.01*f50, 4*f50, 100)

    compl = scube.return_completeness(array(fluxes), coord.ra.deg*ones(len(fluxes)),
                                      coord.dec.deg*ones(len(fluxes)),
                                      opts.lambda_*ones(len(fluxes)), opts.sncut)
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
    parser.add_argument("--sncut", type=float, default=4.5, 
                        help="S/N cut used") 
    parser.add_argument("--fout", type=str, help="Filename to output to", 
                        default=None)
    opts = parser.parse_args(args=args)

    coord = SkyCoord(opts.ra, opts.dec)
    print("WARNING using fixed alpha=-3.1")
    scube = SensitivityCube.from_file(opts.filename, [3500.0, 5500.0], [-3.1, -3.1])

    wls = linspace(3500, 5490.0, 1000)
    f50 = scube.get_f50(coord.ra.deg*ones(len(wls)), coord.dec.deg*ones(len(wls)), wls, opts.sncut)

    plt.plot(wls, f50/1e-16, "k-", label="Flux at 50% completeness")
    plt.ylabel("Flux $10^{-16}\,$(erg/s/cm$^2$/A)", fontsize=14.0)
    plt.xlabel("Wavelength (A)", fontsize=14.0)
    plt.legend(loc="upper right")

    if opts.fout:
        plt.savefig(opts.fout)
    else:
        plt.show()
