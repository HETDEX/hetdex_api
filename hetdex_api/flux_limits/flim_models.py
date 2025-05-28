"""

This module stores different models 
to convert between the values in the 
sensitivity cubes and the flux at
50% detection completeness. It also
stores the tools to interpolate over
simulation results. 

.. moduleauthor:: Daniel Farrow <dfarrow@mpe.mpg.de>

"""

from glob import glob
from os.path import join
from scipy.interpolate import (interp1d, interp2d, splrep, griddata, RectBivariateSpline,
                               NearestNDInterpolator)
from numpy import (polyval, mean, median, loadtxt, meshgrid, savetxt, 
                   array, linspace, tile, ones, array, argmin, zeros, 
                   sqrt)
from hetdex_api.config import HDRconfig
from hetdex_api.flux_limits import flim_models_old, flim_model_cache

class NoLineWidthModel(Exception):
    pass

class FlimOptionNotSupported(Exception):
    pass

class ModelInfo(object):
    """
    Store information about the flux
    limit model

    Parameters
    ----------
    snfile_version : str
       which version of the sn?.?.use 
       files to use
    snpoly : array
       a polynomial to go from S/N
       cut to noise multiplier
    wavepoly : array or NoneType
       a polynomial to adjust the
       noise as a function of
       wavelength (None to not use)
    """
    def __init__(self, snfile_version, snpoly, wavepoly,
                 interp_sigmas, dont_interp_to_zero, 
                 snlow=4.8, snhigh=7.0, lw_scaling=None,
                 wl_collapse=False, single_snfile=False):

        self.snfile_version = snfile_version
        self.snpoly = snpoly
        self.wavepoly = wavepoly
        self.interp_sigmas = interp_sigmas
        self.snlow = snlow
        self.snhigh = snhigh
        self.dont_interp_to_zero = dont_interp_to_zero
        self.lw_scaling = lw_scaling
        self.wl_collapse = wl_collapse
        self.single_snfile = single_snfile

class NoSNFilesException(Exception):
    pass


def linewidth_f50_scaling_v1(linewidth, sncut):
    """
    Model for how the 50% completeness flux
    changes with linewidth. Based on Karl's
    simulation calibrated model from 
    July 11th 2021, quoting from Karl
    Gebhardt:

    x12=min(x12,12.)
    rfit=sqrt(x12/2.2)*(1.64-(sncut-4.8)/(9.5-4.8))/(1.+0.07*(12.-x12))
    rfit=max(1.,rfit)

    Parameters
    ----------
    linewidth : array or float
        linewidth in A
    sncut : float
        S/N cut

    Returns
    -------
    rfit : array or float
        scaling of 50% flux
        with linewidth
    """
    
    try:
        linewidth = array(linewidth)
        linewidth[linewidth > 12.0] = 12.0
    except TypeError:
        linewidth = max(linewidth, 12.0)

    rfit=sqrt(linewidth/2.2)*(1.64-(sncut-4.8)/(9.5-4.8))/(1.+0.07*(12.-linewidth))

    try:
        rfit[rfit < 1.0] = 1.0
    except TypeError:
        rfit = max(1.0, rfit)

    return rfit

def linewidth_f50_scaling_v2(linewidth, sncut):
    """
    Model for how the 50% completeness flux
    changes with linewidth. Based on section 
    3.1 of Liu et al (2022) 

    Parameters
    ----------
    linewidth : array or float
        linewidth in A
    sncut : float
        S/N cut

    Returns
    -------
    rfit : array or float
        scaling of 50% flux
        with linewidth
    """
    
    try:
        linewidth = array(linewidth)
        linewidth[linewidth > 12.0] = 12.0
    except TypeError:
        linewidth = max(linewidth, 12.0)

    rfit=sqrt(linewidth/2.2)*(15.09 - sncut)
    rfit = rfit*(0.098 - 0.0004 * linewidth + 0.0004 * linewith ** 2)

    try:
        rfit[rfit < 1.0] = 1.0
    except TypeError:
        rfit = max(1.0, rfit)

    # XXX HACK
    return 1.0




def write_karl_file(fn, f50s, fcens, wlcens, 
                    completeness_2d):

    with open(fn, 'w') as fp:
        fp.write("0          ")
        for wl in wlcens:
            fp.write("{:5.1f}     ".format(wl))
        fp.write("\n")

        fp.write("0.500      ")
        for f50 in f50s:
            fp.write("{:5.3f}     ".format(f50))
        fp.write("\n")
 
        out_arr = zeros((completeness_2d.shape[0] + 1, completeness_2d.shape[1]))
        out_arr[0, :] = fcens
        out_arr[1:, :] = completeness_2d
        savetxt(fp, out_arr.T, fmt="%.4f    ")


def write_sn_file(fn, sns, wlcens, completeness_2d):

    out_arr = zeros((completeness_2d.shape[0] + 1, completeness_2d.shape[1] + 1))

    out_arr[0, 0] = 0
    out_arr[1:, 0] = wlcens
    out_arr[0, 1:] = sns
    out_arr[1:, 1:] = completeness_2d

    savetxt(fn, out_arr.T, fmt="%.4f    ")


def read_sn_file(fn):

    compl_file = loadtxt(fn)
    waves = compl_file[0,1:]
    compl_curves = compl_file[1:, 1:].T
    sns = compl_file[1:, 0]

    return waves, sns, compl_curves



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
    single_snfile : bool
        Use the same file of completeness
        curves for all S/N cuts
    **kwargs : 
        passed to SingleSNSimulationInterpolator
    """

    def __init__(self, fdir, dont_interp_to_zero, 
                 single_snfile=False,
                 snmode=False, verbose=False, **kwargs):

        if single_snfile:
            if snmode:
                raise FlimOptionNotSupported("Using single_snfile and ",
                                              "snmode at same time ", 
                                              "not supported")
            if dont_interp_to_zero:
                raise FlimOptionNotSupported("Using dont_interp_to_zero ",
                                              "and single_snfile at same time ",
                                              "not supported")
 

        if single_snfile:
            if verbose:
                print("Using a single S/N file for everything!")
            snfiles = glob(fdir + "/sn_all.use")
        elif snmode:
            snfiles = glob(fdir + "/sn_based_?.?.dat")
        else: 
            snfiles = glob(fdir + "/sn?.?.use")

        self.dont_interp_to_zero = dont_interp_to_zero

        sns = []
        self.sninterpolators = []
        for snfile in snfiles:
            if verbose:
                print(snfile)
            self.sninterpolators.append(SingleSNSimulationInterpolator(snfile, dont_interp_to_zero,
                                                                       snmode=snmode, **kwargs))
            if single_snfile:
                # value doesn't matter
                sns.append(5.5)
            else:
                sns.append(float(snfile[-7:-4]))

        if len(sns) == 0:
            raise NoSNFilesException("Could not find any simulation files! ")

        if verbose:
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
    """
    def __init__(self, filename, dont_interp_to_zero, wl_collapse = False, 
                 cmax = None, snmode = False):

        if not snmode:
            self.waves, self.f50, self.compl_curves, self.fluxes =\
                read_karl_file(filename)
        else:  
            self.waves, self.fluxes, self.compl_curves = \
                read_sn_file(filename)
            # set 50% to one as already in S/N units
            self.f50 = ones(len(self.waves))

        # Check bin spacing uniform
        offs = self.fluxes[1:] - self.fluxes[:-1]
        if max(abs(offs - offs[0])) > 1e-20:
            raise Exception("Bin spacing must be uniform!")

        self._snmode = snmode
        self._wl_collapse = wl_collapse
        self.cmax = cmax 
        self.dont_interp_to_zero = dont_interp_to_zero

        # Optionally normalize all curves to cmax
        if cmax:
            print("Normalizing to {:f}".format(cmax))
            for i in range(len(self.waves)):
                 self.compl_curves[i, :] /= max(self.compl_curves[i, :]) 
                 self.compl_curves[i, :] *= cmax

        self.completeness_model = self.interpolated_model()


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
        # at 1.0, also should in combination will fill_value
        # ensure 0 returns zero completeness in the 
        # RectBivariateSpline 
        fluxes_f50_units = linspace(0, 50, 5000)

        c_all = []

        fgrid_for_mask = []
        wgrid_for_mask = []                
               
 
        # Offset by half a bin brighter, so the mask kicks in
        # at the center location, not 1/2 a bin away. The
        # 0.999 factor is to ensure the actual value itself
        # is not in the mask
        fbsizediv2 = 0.999*(self.fluxes[1] - self.fluxes[0])/2.0

        for twave, tf50, c in zip(self.waves, self.f50, self.compl_curves):

            if plot:
                plt.plot(self.fluxes/tf50, c, linestyle="--")
           
               
            # Shift grid to center of bin, so don't
            # interpolate past value
            fgrid_for_mask.append((self.fluxes + fbsizediv2)/tf50)
            wgrid_for_mask.append(ones(len(self.fluxes))*twave)

            # divide so 50% completeness to
            # convert to flux units of per f50
            interpolator = interp1d(self.fluxes/tf50, c, bounds_error = False, 
                                    fill_value=(0.0, c[-1]))

            # interpolate to the coordinates where 50% is at 1.0
            c_all.append(interpolator(fluxes_f50_units))


        c_all = array(c_all)

        # Make a combined model
        if self._wl_collapse or (len(self.waves) == 1): 

            if len(self.waves) == 1:
                cmean = c_all[0]
            else:
                # Don't use first two wavebins as
                # big scatter
                cmean = mean(c_all[2:], axis=0)
        
            completeness_model = interp1d(fluxes_f50_units, cmean, 
                                          fill_value=(0.0, cmean[-1]), 
                                          bounds_error=False)

            if plot:
                vals_to_plot = completeness_model(fluxes_f50_units)

        else:

            # waves have to be uniformly spaced for this to work (? don't think so?)
            interp = RectBivariateSpline(self.waves, fluxes_f50_units, c_all, 
                                         kx=3, ky=3)

            if self.dont_interp_to_zero:
                # Use this as a mask to not extrapolate toward 0.0
                # if nearest point is zero
                compl_mask = zeros(self.compl_curves.shape)
                compl_mask[self.compl_curves > 0.0] = 1.0
                interp_mask = NearestNDInterpolator(list(zip(array(wgrid_for_mask).ravel(), 
                                                             array(fgrid_for_mask).ravel())), 
                                                             compl_mask.ravel())

                completeness_model = lambda x, y : interp_mask(x, y)*interp(x, y, grid=False)
            else:
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

        if self._snmode:
            return self.completeness_model(f50, flux/f50)
        elif self._wl_collapse or (len(self.waves) == 1):
            return self.completeness_model(flux/f50)
        else:
            return self.completeness_model(wave.flatten(), (flux/f50).flatten())


    # Blue end noise adjustment
    #params_wl = [6.18971170e-07, -5.15522003e-03,  1.17598910e+01] 
    #
    #try:
    #    noise[lambda_ < 4000.0] = noise[lambda_ < 4000]*polyval(params_wl, lambda_[lambda_ < 4000.0])
    #except TypeError:
    #    if lambda_ < 4000.0:
    #        noise = noise*polyval(params_wl, lambda_)



def return_flux_limit_model_old(flim_model):
    """
    Return the noise -> 50% completeness
    scaling and a function for the completeness
    curves
    """

    f50_from_noise = getattr(flim_models_old, "{:s}_f50_from_noise".format(flim_model)) 
   
    return f50_from_noise, None, False


def return_flux_limit_model(flim_model, cache_sim_interp = True, 
                            verbose = False):
    """
    Return the noise -> 50% completeness
    scaling and a function for the 
    completeness curves

    """

    # old models for legacy support
    if flim_model in ["hdr1", "hdr2pt1"]:
        if verbose:
            print("Using flim model: {:s}".format(flim_model))
        return return_flux_limit_model_old(flim_model)    

    models = {
              "one_sigma_nearest_pixel" : ModelInfo("curves_v1", 
                                              [1.0, 0.0], 
                                              None, 
                                              False, False, 
                                              snlow=0.999999, snhigh=1.000001),
              "one_sigma_interpolate" : ModelInfo("curves_v1", 
                                                  [1.0, 0.0], 
                                                  None, 
                                                  True, False, 
                                                  snlow=0.999999, snhigh=1.000001), 
              "hdr2pt1pt1" : ModelInfo("curves_v1", 
                                       [2.76096687e-03, 2.09732448e-02, 7.21681512e-02, 3.36040017e+00], 
                                       None, False, False),
              "hdr2pt1pt3" : ModelInfo("curves_v1",
                                       [6.90111625e-04, 5.99169372e-02, 2.92352510e-01, 1.74348070e+00],
                                       None, False, False),
              "v1" : ModelInfo("curves_v1", 
                               [-8.80650683e-02,  2.03488098e+00, -1.73733048e+01, 
                               6.56038443e+01, -8.84158092e+01], 
                               None, False, True),
              "v1.1" : ModelInfo("curves_v1", 
                               [-8.80650683e-02,  2.03488098e+00, -1.73733048e+01, 
                               6.56038443e+01, -8.84158092e+01], 
                               None, True, True),
              "v2" : ModelInfo("curves_v1",
                               [1.0, 0.0], 
                               [-1.59395767e-14, 3.10388106e-10, 
                                -2.26855051e-06,  7.38507004e-03, 
                                -8.06953973e+00], False, True,
                                lw_scaling=linewidth_f50_scaling_v1),
              "v3" : ModelInfo("curves_v1",
                                [-0.04162407, 0.80167981, -4.11209695, 9.95753597],
                                [1.46963918e-15, -6.68766843e-11,  7.56849155e-07, -3.28661164e-03,
                                 5.95152597e+00],
                                False, True,
                                lw_scaling=linewidth_f50_scaling_v1,
                                wl_collapse=True),
              "v4" : ModelInfo("curves_v2",
                                [1.0, 0.0],
                                [-4.10848216e-11, 6.68571465e-07, 
                                 -3.55515621e-03,  6.99466106e+00],
                                False, False,
                                lw_scaling=linewidth_f50_scaling_v1,
                                single_snfile=True)
              }


    default = "v4" 

    if not flim_model:
        flim_model = default
    
    model =  models[flim_model]
    if verbose:
        print("Using flim model: {:s}".format(flim_model))
 
    if flim_model_cache.cached_model == flim_model and cache_sim_interp:
        sinterp = flim_model_cache.cached_sim_interp
    else:
        conf = HDRconfig()
        fdir = conf.flim_sim_completeness 
        fdir = join(fdir, model.snfile_version)
        sinterp = SimulationInterpolator(fdir, model.dont_interp_to_zero, 
                                         snmode=False, verbose=verbose,
                                         wl_collapse=model.wl_collapse,
                                         single_snfile=model.single_snfile)

    # save model in cache
    if cache_sim_interp:
        flim_model_cache.cached_model = flim_model
        flim_model_cache.cached_sim_interp = sinterp


    def f50_from_noise(noise, lambda_, sncut, linewidth = None):
        """
        Return the 50% completeness
        flux given noise and S/N cut. 
    
        Parameters
        ----------
        noise : float
            the noise from the
            sensitivity cubes
        sncut : float
            the signal to noise
            cut to assume
        linewidth : float
            the linewidth in A,
            only used with supported
            models    

        Returns
        -------
        f50s : array
           the fluxes at 50%
           completeness  
        """
        try:
            if sncut < model.snlow or sncut > model.snhigh:
                print("WARNING: model {:s} not calibrated for this S/N range".format(flim_model))
        except ValueError:
            if any(sncut < 4.5) or any(ncut > 7.5):
                print("WARNING: model {:s} not calibrated for this S/N range".format(flim_model))

        bad = noise > 998

        if model.wavepoly:
            noise = noise*polyval(model.wavepoly, lambda_)
 
        if type(linewidth) != type(None):
            if type(model.lw_scaling) != type(None):
                lw_scale = model.lw_scaling(linewidth, sncut)
            else:
                raise NoLineWidthModel("Linewidth dependence not available for this flim model")

            noise = noise*lw_scale
           
        snmult = polyval(model.snpoly, sncut)
        f50 = snmult*noise

        # keep bad values unscaled
        try:
            f50[bad] = 999
        except TypeError:
            if bad:
                f50 = 999

        return f50

    return f50_from_noise, sinterp, model.interp_sigmas  
