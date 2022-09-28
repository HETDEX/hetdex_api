"""

Derive a scaling between the values in 
from the ShotSensitivity API and the completeness
simulations Karl Gebhardt ran.  

Daniel Farrow (MPE) 2021, 2022

"""
from numpy import (loadtxt, savetxt, transpose, interp, sqrt, exp, 
                   mean, linspace, zeros, array, polyfit, polyval, 
                   abs, std, unique, histogram, histogram2d, interp)
from numpy.random import uniform
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import TABLEAU_COLORS
from astropy.table import Table, vstack
from astropy.stats import biweight_location 

from hetdex_api.extinction import get_2pt1_extinction_fix

from hetdex_api.flux_limits.shot_sensitivity import ShotSensitivity 
from hetdex_api.flux_limits.flim_models import (SimulationInterpolator, read_karl_file, 
                                                return_flux_limit_model)


FS=16.0
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['mathtext.fontset'] = 'dejavuserif'
mpl.rcParams["xtick.labelsize"] = 15.0
mpl.rcParams["ytick.labelsize"] = 15.0
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["ytick.direction"] = "in"
mpl.use('tkagg')

def biweight_one_sigma_from_table(table, wave_bins):
    """
    Compute the biweight midvariance of 
    the noise in wavelength bins in
    an astropy table

    Parameters
    ----------
    table : astropy.table:Table
        a table containing 
        `flux_noise_1sigma` and
        `wave`
    wave_bins : array
        edges of wavelength bins

    """

    waves = 0.5*(wave_bins[1:] + wave_bins[:-1])  

    sigmas = []
    for wl, wh in zip(wave_bins[:-1], wave_bins[1:]):
        noise = table["flux_noise_1sigma"][(table["wave_in"] > wl) & (table["wave_in"] <= wh)]
        sigmas.append(biweight_location(noise[(noise < 100.0) & (noise > 0)], 
                      ignore_nan = True))

    return waves, array(sigmas)


def measure_completeness(table, wave_bins, sn):
    """
    Measure completeness in wavelength
    and flux bins

    """
    fluxes = linspace(0.5, 50, 20)

    flux_cen = 0.5*(fluxes[1:] + fluxes[:-1])
    wave_cen = 0.5*(wave_bins[1:] + wave_bins[:-1])

    hist2d_all = histogram2d(table["flux_in"], table["wave_in"],
                             bins=(fluxes, wave_bins))[0]

    det = (table["sn"] > sn) & (table["chi"] < 2.5)
 
    hist2d_det = histogram2d(table["flux_in"][det], 
                             table["wave_in"][det],
                             bins=(fluxes, wave_bins))[0]
 
    compl = 1.0*hist2d_det/hist2d_all

    f50s = []
    for j in range(compl.shape[1]):
        #compl[:, j] /= compl[-1, j]
        f50 = interp(0.5, compl[:, j]/compl[-1, j], flux_cen)
        f50s.append(f50) 
        print(sn, wave_cen[j], f50)

    return wave_cen, array(f50s), compl, fluxes



def measure_sigma_to_f50(fout, split_files=False, use_karl_file=False):
    """
    Measure the conversion between the
    noise measured in apertures, to
    the flux at 50% completeness from
    the simulation.

    Parameters
    ----------
    fout : str
        filename to output the
        scaling to

    """

    s = ShotSensitivity("20200519v013", flim_model="v4")

    # read in the shots to consider 
    if split_files: 
        with open("datevobs_list_cal", 'r') as fp:
            rshots = fp.readlines()
        shots = [x.strip() for x in rshots]
    else:
        table_full = Table.read("combined_onfly_3.5_3.fits")
        shots = unique(table_full["datevshot"])


    fix = get_2pt1_extinction_fix()

    # get rid of this for > hdr2.1
    fix = lambda x : 1.0
 
    all_tables = []
    wave_bins_fine = linspace(3500, 5500, 200)
    wave_bins = array([3550., 3821., 4093., 4364., 4635., 4907., 5178., 5450.])
    snlist = [4.8, 5.0, 5.5, 6.0, 6.5, 7.0]
 
    for shot in shots:
 
        # For each shot measure the average flux noise 1 sigma
        # This should be a table of API-extracted noise
        # values for the cal shots
        if split_files:
            table = Table.read("{:s}_full_input_0.fits".format(shot))
        else:
            table = table_full[table_full["datevshot"] == shot]


        table = table[table["norm_3_3.5"] > 0.6]
        src_waves = table["wave_in"]

        # above atmosphere units
        src_noise = table["sigma_3_3.5"]
        # Need this later
        table["flux_noise_1sigma"] = src_noise
        #src_noise = table["flux_noise_1sigma"]

        wlav, sigma_av = biweight_one_sigma_from_table(table, wave_bins_fine)
        wlav_binned, sigma_av_binned = biweight_one_sigma_from_table(table, wave_bins)
 
        t = Table()
        t["sn"] = snlist
        t["datevobs"] = shot 
        scales = zeros((len(snlist), len(wlav_binned)))
 
        # Now loop over deriving a scaling 
        for i, sn in enumerate(snlist):
            if use_karl_file:
                waves, f50, compl_curves_base, fluxes_base =\
                    read_karl_file("karl_f50_stuff/headfix{:s}.sn{:2.1f}.fcor".format(shot, sn))
            else:
                waves, f50, compl_curves_base, fluxes_base = measure_completeness(table, wave_bins, sn)


            print(waves, f50, compl_curves_base, fluxes_base)

            f50 = f50/fix(waves)

            scales[i, :] = f50*1e-17/sigma_av_binned

            print(">>>>>")
            print(scales[i, :])
            print(s.f50_from_noise(sigma_av_binned, waves, sn)/sigma_av_binned)
            

            # Correct for the scatter within the shot
            for j in range(scales.shape[1]):
                sel = (src_waves > wave_bins[j]) & (src_waves <= wave_bins[j + 1])
                scales[i, j] = derive_corrections_to_conversion(src_noise[sel], src_waves[sel],
                                                                f50[j]*1e-17, scales[i, j], sn, 
                                                                fluxes_base/fix(waves[j]),
                                                                compl_curves_base[j])
         
        for i, wl in enumerate(wlav_binned):
            t["wave_{:4.0f}".format(wl)] = scales[:, i] 
 
        all_tables.append(t)
 
    table = vstack(all_tables)
    table.write(fout)

    return table

def derive_corrections_to_conversion(flux_noise_1sigma, waves,
                                     f50, guess, sncut,
                                     fluxes_base, compl_curves_base,
                                     maxiter=20, frac=0.02):
    """
    Derive corrections to the 1 sigma to f50 conversion 
    iteratively by evaluating the completeness model
    on a set of sampled positions. The goal is to
    match the results of source simulations.  

    Parameters
    ----------
    flux_noise_1sigma : array 
        the noise in the sampled positions
    waves : array
        the wavelength of the sampled 
        positions
    f50 : float
        the target 50% completeness
        as measured from the source simulations
    guess : float
        initial guess for scaling factor
    sncut : float   
        the S/N threshold
    fluxes_base, compl_curves_base : array 
        the fluxes and completeness curves
        from the source simulation
    maxiter : int
        maximum iterations of evaluating the
        completeness model before giving up
    frac : float
        target fractional accuracy on the
        conversion to f50


    Returns
    ------- 
    guess : float
       the new value for the conversion between
       the source simulation 50% completeness flux
       and the noise in apertures from the API

    """ 

    # datevshot just needed for the completeness interpolator
    s = ShotSensitivity("20200519v013", flim_model="v4")
    flux_bins = linspace(5e-18, 1e-15, 60)
    fbcens = 0.5*(flux_bins[1:] + flux_bins[:-1])

    fluxes = uniform(5e-18, 1e-15, size=len(flux_noise_1sigma)) 
    hist_in = histogram(fluxes, bins=flux_bins)[0]
    difflast = 1e99
    guesslast = 0.
   
    for i in range(maxiter):

        compl = s.return_completeness(fluxes, None, None, waves, sncut,
                                      f50s = guess*flux_noise_1sigma)

 
        hist_out = histogram(fluxes, bins=flux_bins, 
                             weights=compl)[0]
   
        ratios = 1.0*hist_out/hist_in 
        f50_new = interp(0.5, ratios, #/max(ratios),
                         fbcens) 

        #plt.plot(fbcens, ratios/max(ratios))
        #plt.plot(fluxes_base*1e-17, compl_curves_base/max(compl_curves_base), "k:")
        #plt.axvline(f50)
        #plt.axvline(f50_new)
        #plt.show()
 
        diff = (f50_new - f50)/f50_new
        i = i + 1 
        print("Iteration {:d} ".format(i), f50, f50_new, f50/f50_new, diff)

        if abs(diff) < 0.02:
            break

        # If we're getting worse decrease step size
        # and go back
        if abs(diff) > abs(difflast):
            frac = frac/2.0
            guess = guesslast
            diff = difflast
        else:
            difflast = diff
            guesslast = guess

        if diff > 0.0:
            guess = guess - frac*guess
        else:
            guess = guess + frac*guess

    if not abs(diff) < frac:
        print("Couldn't meet spec here: ", min(waves), max(waves), diff)

    return guess

def plot_sigma_to_f50_scaling(table):
    """

    Plot the scaling between the 1 sigma
    value in apertures and the 50%
    completeness flux. Add the current
    model fits as a dashed line.

    Parameters
    ----------
    table : astropy.table:Table
        the table of simulation 
        measured scaling values


    """
    fig = plt.figure(figsize=(9.5,8))
    sns = unique(table["sn"])
    scalecens = []
    scalecens_err = []
    wavetrends = []
    cols = []

 
    f50_from_noise, sinterp, interp_sigmas = \
        return_flux_limit_model("v5")
 
    for color, sn in zip(TABLEAU_COLORS, sns):
        there = table[abs(table["sn"] - sn) < 1e-8]

        waves = []
        mscales = []
        stds = []
        cols.append(color)

        for colname in there.colnames:
           if "wave" in colname:
                waves.append(float(colname.replace("wave_", "")))
                mscales.append(mean(there[colname]))
                stds.append(std(there[colname]))
 
        scalecens.append(mscales[4])
        scalecens_err.append(stds[4])
        wavetrends.append(mscales/mscales[4])
        #plt.axhline(sn, linestyle=":", color=color)

        model_vals = f50_from_noise(0*array(waves) + 1.0, waves, sn)   
        plt.plot(waves, model_vals, "--", color=color)
        plt.errorbar(waves, mscales, yerr=stds, marker="o", 
                     label = "S/N $>$" +  "{:2.1f}".format(sn),
                     color=color, linestyle="None")

    # Parts adapted from the official matplotlib tutorial
    cmap = mpl.colors.ListedColormap(cols)
    bounds = range(len(sns) + 1)
    cens = 0.5*(array(bounds[1:]) + array(bounds[:-1]))
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    #cax = plt.gca().inset_axes([5600, 4.2, 80.0, 3.4],
    #                           transform=plt.gca().transData) 

    cax = plt.gca().inset_axes([5600, 3.98, 80.0, 2.92],
                               transform=plt.gca().transData) 


    cb = fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm),
                      cax=cax, ticks=cens,
                      boundaries=bounds)

    cb.ax.set_yticklabels(sns)
    cb.set_label("S/N cut", fontsize=FS)   
    plt.xlabel("Wavelength (A)", fontsize=FS) 
    plt.ylabel(r"Simulation-derived $f_{\rm 50}/\sigma_{\rm aper}$ ratio", fontsize=FS) 
    plt.tight_layout()

    return sns, waves, scalecens, scalecens_err, wavetrends, mscales, stds

def fit_wavetrend(waves, wavetrend):
    """
    Fit the wavelength trends of the
    1 sigma -> f50 conversion.

    Parameters
    ----------
    wave : array
        the wavelengths
    wavetrend : array
        2D array of the 
        simulation derived f50/1sigma
        ratios for the different 
        calibration shots

    Return
    ------
    vals : array
        the best fitting polynomial
 
    """
    plt.figure(figsize=(9,8))
    for wavetrend in wavetrends:
        plt.plot(waves, wavetrend)
    mean_wavetrend = mean(wavetrends, axis=0)
    vals = polyfit(waves, mean_wavetrend, 4)
    plt.plot(waves, polyval(vals, waves), "k:",
             linewidth=8, label="Best-fit")

    plt.legend(frameon=False, 
               prop={"size" : FS})
    plt.ylabel("Wavelength dependence", 
               fontsize=FS) 
    plt.xlabel("Wavelength (A)", 
               fontsize=FS)
    plt.tight_layout()
 
    return vals

def fit_sntrend(sns, scalecens, scalecens_err):
    """

    Fit the conversion between 1-sigma and
    50% completeness for different S/N
    thresholds.

    Parameters
    ----------
    sns : array
        the S/N theshholds
    scalecens, scalecens_err :
        the ration f50/1 sigma
        for different S/N
        cuts and its error

    Return
    ------
    vals : array
        the best fitting polynomial 

    """
    plt.figure(figsize=(9,8))

    vals = polyfit(sns, scalecens, 3)

    plt.errorbar(sns, scalecens, linestyle="none",
                 yerr=scalecens_err,
                 marker="o", markersize=8)
    plt.plot(sns, polyval(vals, sns), "k:",
             linewidth=4, label="Best-fit")

    plt.legend(frameon=False, 
               prop={"size" : FS})
    plt.ylabel(r"$f_{50}/\sigma_{\rm aper}$", 
               fontsize=FS) 
    plt.xlabel("S/N", 
               fontsize=FS)
    plt.tight_layout()
 
    return vals


if __name__ == "__main__":

    # If refitting make sure to use correct
    # flim model when calling ShotSensitivity

    remeasure = False
    fscale = "scaling.fits"

    if remeasure:
        table = measure_sigma_to_f50(fscale) 
    else:
        table = Table.read(fscale) 


    sns, waves, scalecens, scalecens_err, wavetrends, allscales, allstd = \
        plot_sigma_to_f50_scaling(table)


    wave_vals = fit_wavetrend(waves, wavetrends)
    sn_vals = fit_sntrend(sns, scalecens, scalecens_err)

    print("Wavelength poly: ", wave_vals)
    print("S/N poly:", sn_vals)

    plt.show()
