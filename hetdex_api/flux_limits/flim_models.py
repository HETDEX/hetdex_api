"""

This module stores different models 
to convert between the values in the 
sensitivity cubes and the flux at
50% detection completeness. It also
stores the tools to interpolate over
simulation results. 

.. moduleauthor:: Daniel Farrow <dfarrow@mpe.mpg.de>

"""

import matplotlib.pyplot as plt
from glob import glob
from scipy.interpolate import interp1d, interp2d, splrep, griddata, RectBivariateSpline
from numpy import (polyval, mean, median, loadtxt, meshgrid, 
                   array, linspace, tile, ones, array, argmin)

def read_karl_file(fn):
    """ 
    Read in the results from
    Karl Gebhardt's simulation

    Parameters
    ----------
    fn : str
        filename of Karl's
        fcor.use file 

    Returns
    -------
    waves : array
        wavelengths  at
        which simulation has
        results
    f50 : array
        50% completeness at
        each wavelength (in
        10^-17 erg/s/cm2 units)
    compl_curves : 2D array
        2D array of completeness
        curves at every wavelength
    fluxes : array
        flux bins for completeness
         curves
    """
    karl_det_file = loadtxt(fn)
    waves = karl_det_file[0,1:]
    f50 = karl_det_file[1,1:]
    compl_curves = karl_det_file[2:, 1:].T
    fluxes = karl_det_file[2:, 0]

    return waves, f50, compl_curves, fluxes


class SimulationInterpolator(object):
    """
    Interpolate over the results of source
    simulations. Use nearest neighbour 
    interpolation to select the shape of
    the completeness given a S/N cut (i.e.
    which of Karl's files to use)  

    Parameters
    ----------
    fdir : str
        A directory of sn?.?.use files,
        containing Karl's simulation 
        results
    **kwargs : 
        passed to SingleSNSimulationInterpolator
    """

    def __init__(self, fdir, **kwargs):

        snfiles = glob(fdir + "/sn?.?.use")

        sns = []
        self.sninterpolators = []
        for snfile in snfiles:
            self.sninterpolators.append(SingleSNSimulationInterpolator(snfile, **kwargs))
            sns.append(float(snfile[-7:-4]))

        print("Read S/N files for the following cuts: {:s}".format(str(sns)))
        self.sns = array(sns)        


    def __call__(self, flux, f50, wave, sncut, verbose=False):
        """
        Find the nearest S/N versus completeness curve 
        from the simulations and use it to predict
        completeness

        Parameters
        ----------
        flux : array
            fluxes to return completeness
            at
        f50 : array
            50% completeness fluxes
        wave : array
            wavelengths (in A)
        sncut : float
            the S/N cut applied
        """

        diffs = abs(self.sns - sncut)
        iclosest = argmin(diffs)
   
        if verbose:
            print("Using simulation from S/N cut {:f}".format(self.sns[iclosest]))
        
        return self.sninterpolators[iclosest](flux, f50, wave)


class SingleSNSimulationInterpolator(object):
    """
    Interpolate over the results of source
    completeness simulations, of the form
    run by Karl Gebhardt, for one S/N cut only

    Parameters
    ----------
    filename : str
        Karl's fcor.use file

    wl_collapse : bool
        If True collapse
        wavelength bins in 
        the simulation file to
        use one curve for all
        wavelengths (default: False)

    cmax : float
       Maximum completeness to normalize
       the curves to (default: 0.98)

    Attributes
    ----------
    waves : array
        wavelengths at which 
        completeness curves are 
        measured
    f50 : array
        50% completeness values for
        each wavelength
    compl_curves : array
        2D array of completeness 
        in flux and wavelength bins
    fluxes : array
        the fluxes (*1e17 erg/s/cm) at
        which completeness is tabulated
    completeness_model : callable
        model of completeness given 
        flux, wavelength and noise
    f50_interpolator : callable
        given wavelength, return
        50% completeness from the
        simulation files

    """
    def __init__(self, filename, wl_collapse = False, cmax = 1.0):

        self.waves, self.f50, self.compl_curves, self.fluxes =\
            read_karl_file(filename)

        self._wl_collapse = wl_collapse
        self.cmax = cmax 

        # Optionally normalize all curves to cmax
        if cmax:
            for i in range(len(self.waves)):
                 self.compl_curves[i, :] /= max(self.compl_curves[i, :]) 
                 self.compl_curves[i, :] *= cmax

        self.completeness_model = self.interpolated_model()
        self.f50_interpolator = interp1d(self.waves, 1e-17*self.f50) 


    def interpolated_model(self, plot=False):
        """ 
        Generate an interpolated model from
        the completeness curves. Waves have to
        be uniformly spaced
        """

        if plot:
            import matplotlib.pyplot as plt
            import matplotlib as mpl
            mpl.use("tkagg")

        # bins to interpolate all completeness curves to,
        # in these coordinates 50% completeness always 
        # at 1.0 
        fluxes_f50_units = linspace(0, 10, 2000)

        c_all = []

        for tf50, c in zip(self.f50, self.compl_curves):

            if plot:
                plt.plot(self.fluxes/tf50, c, linestyle="--")
 
            # divide so 50% completeness to
            # convert to flux units of per f50
            interpolator = interp1d(self.fluxes/tf50, c, bounds_error = False, 
                                    fill_value=(0.0, c[-1]))

            # interpolate to the coordinates where 50% is at 1.0
            c_all.append(interpolator(fluxes_f50_units))


        c_all = array(c_all)

        # Make a combined model
        if self._wl_collapse:
            cmean = mean(c_all, axis=0)
        
            completeness_model = interp1d(fluxes_f50_units, cmean, 
                                          fill_value=(0.0, cmean[-1]), 
                                          bounds_error=False)

            if plot:
                vals_to_plot = completeness_model(fluxes_f50_units)

        else:
            # waves have to be uniformly spaced for this to work
            interp = RectBivariateSpline(self.waves, fluxes_f50_units, c_all, 
                                         kx=3, ky=3)

            completeness_model = lambda x, y : interp(x, y, grid=False)

            if plot:
                vals_to_plot = completeness_model(self.waves[2]*ones(len(fluxes_f50_units)),
                                                  fluxes_f50_units)


        if plot:
            plt.plot(fluxes_f50_units, vals_to_plot, "k.", lw=2.0)
            plt.xlim(0, 12.0)
            plt.xlabel("Flux/(50% flux) [erg/s/cm^2]")
            plt.ylabel("Normalized Completeness")
            plt.show() 

        return completeness_model


    def __call__(self, flux, f50, wave):
        """
        Return the completeness

        Parameters
        ----------
        flux : array
            fluxes to return completeness
            at
        f50 : array
            50% completeness fluxes
        wave : array
            wavelengths (in A)

        """
        flux = array(flux)
        f50 = array(f50)
        wave = array(wave)

        if self._wl_collapse:
            return self.completeness_model(flux/f50)
        else:
            return self.completeness_model(wave.flatten(), (flux/f50).flatten())

def hdr2pt1pt1_f50_from_noise(noise, sncut):
    """
    Return the 50% completeness
    flux given noise and 
    S/N cut. This uses an empirical
    model calibrated on simulated
    LAE detections inserted into
    the calibration fields. The
    simulations were run by Karl Gebhardt.

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
    try:
        if sncut < 4.8 or sncut > 7.0:
            print("WARNING: model not calibrated for this S/N range")
    except ValueError:
        if any(sncut < 4.5) or any(ncut > 7.5):
            print("WARNING: model not calibrated for this S/N range")

    params = [2.76096687e-03, 2.09732448e-02, 7.21681512e-02, 3.36040017e+00]
    snmult = polyval(params, sncut)

    return noise*snmult


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

    try:
        if sncut < 4.5 or sncut > 7.5:
            print("WARNING: model not calibrated for this S/N range")
    except ValueError:
        if any(sncut < 4.5) or any(ncut > 7.5):
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

