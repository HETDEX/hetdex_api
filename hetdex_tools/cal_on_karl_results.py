"""

Derive a scaling between the values in the
flux limit cubes and the completeness
simulations Gebhardt ran.

Daniel Farrow (MPE) 2020

"""

from numpy import loadtxt, savetxt, transpose, interp, sqrt, mean, linspace, zeros
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.special import erf
from scipy.optimize import leastsq
from hetdex_api.flux_limits.hdf5_sensitivity_cubes import (return_biweight_one_sigma, return_sensitivity_hdf_path,
                                                           SensitivityCubeHDF5Container)
from hetdex_api.flux_limits.sensitivity_cube import fleming_function
from hetdex_api.flux_limits.flim_models import (hdr2pt1_f50_from_noise, SimulationInterpolator, 
                                                hdr2pt1pt1_f50_from_noise)

mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["ytick.direction"] = "in"
mpl.use('tkagg')

SQRT2 = sqrt(2.0)

def gaussian_det_frac(snr, snr_cut):
    """
    Fraction of sources detected assuming
    a Gaussian noise field

    """ 

    return 0.5*erf((snr - snr_cut)/SQRT2) + 0.5

def func(alpha, fluxes, tf50, data):
    """ Function for scipy leastsq to minimize """
    model = fleming_function(fluxes, tf50, alpha) 
    return model - data

def func_f50s(rescale, f50s, f50_pred):
    """ Function for scipy leastsq to minimize """
    return f50s - rescale*f50_pred

def read_karl_file(fn):
    """ 
    Read in the results from
    Karl's simulation
    """
    karl_det_file = loadtxt("fcor.use")
    waves = karl_det_file[0,1:]
    f50 = karl_det_file[1,1:]
    compl_curves = karl_det_file[2:, 1:].T
    fluxes = karl_det_file[2:, 0]

    return waves, f50, compl_curves, fluxes


def plot_f50(waves, f50, wlcube, sigma_cube):
    """
    Plot the 50% compleness values and the
    sigma values from the flux limit
    cubes
    """
    FS = 12.0
    plt.plot(wlcube, 1e17*hdr2pt1_f50_from_noise(sigma_cube, 5.0), "b:", 
             label="HDR2.1 scaled sigmas (S/N=5.0)")
    plt.plot(wlcube, 5e17*sigma_cube, "k--", 
             label="5$\sigma$ from cubes")

    # fit a new scaling
    sigma_f50_pos = interp(waves, wlcube, sigma_cube)
    best, info = leastsq(func_f50s, [5.0], args=(f50[1:], 1e17*sigma_f50_pos[1:]))
    print(best[0], best[0]/hdr2pt1_f50_from_noise(1.0, 5.0))

    plt.plot(wlcube, 1e17*best[0]*sigma_cube, "g-", 
             label="New scaling ({:3.2f}".format(best[0]) + "$\sigma$)")
    plt.plot(waves, f50, "r*", label="Karl sims", markersize=9.0)

    plt.xlabel("Wavelength [AA]", 
               fontsize=FS)
    plt.ylabel(r"50\% flux [10$^{-17} erg/s/cm$^2$]", 
               fontsize=FS)

def return_completeness_from_shots(shots, fluxes, lambda_low, lambda_high,
                                   sncut, rescale_50 = None):
   """
   Return the completeness over a range
   of shots

   """

   bin_edges = linspace(0, 1e-16, 1000)
   bins = 0.5*(bin_edges[:-1] + bin_edges[1:])
   nbinned = zeros(len(bin_edges) - 1)

   cube_comp_all = []
   for shot in shots:
       fn, fnmask = return_sensitivity_hdf_path(shot, return_mask_fn = True)

       with SensitivityCubeHDF5Container(fn, mask_filename = fnmask) as hdf:
           cube_compl, tnbinned = hdf.return_shot_completeness(fluxes*1e-17, lambda_low, 
                                                               lambda_high, sncut,
                                                               bin_edges = bin_edges)
           # rescale curve so 50% flux at rescale_50
           if rescale_50:
               f50 = interp(0.5, cube_compl/cube_compl[-1], fluxes)
               cube_compl = interp(fluxes, rescale_50*fluxes/f50, cube_compl)  

           nbinned += tnbinned

       cube_comp_all.append(cube_compl)
 
   return mean(cube_comp_all, axis=0), bins, nbinned


def plot_completeness_curves(waves, f50, fluxes, compl_curves, shots):
    
    FS = 13.0
    alpha_fix = -3.5 

    alphas_fit = []
    colors =['r', 'g', 'b', 'k', 'm', 'c', 'y']

    fig1, ax1 = plt.subplots(nrows = 2, ncols = 4, 
                             figsize=(14.0, 7))
    fig2, ax2 = plt.subplots(nrows = 1, ncols = 1)
    ax1 = ax1.flatten()

    for i, (col, wl, tf50, compl)  in enumerate(zip(colors, waves, f50, compl_curves)):

        sinterp = SimulationInterpolator("fcor.use")
        #model = compl[-2]*sinterp(fluxes, tf50) 

        rerun = True
        if rerun:

            wllo = wl - 272
            if wllo < 3500:
                wllo = 3500

            compl_hdfs, bins, noise_hist = return_completeness_from_shots(shots, fluxes, 
                                                                          wllo, wl + 272,
                                                                          5.0,
                                                                          rescale_50 = tf50)
            savetxt("av_cube_completeness_cal_{:4.1f}_rescale.txt".format(wl), 
                    transpose([fluxes, compl_hdfs]))
            #savetxt("noise_hist_{:4.1f}.txt".format(wl), 
            #        transpose([bins, noise_hist])) 
        else:
            fluxes_check, compl_hdfs = loadtxt("av_cube_completeness_cal_{:4.1f}_rescale.txt".format(wl), 
                                               unpack=True)
            bins, noise_hist = loadtxt("noise_hist_{:4.1f}.txt".format(wl),
                                       unpack = True)

            assert all((fluxes_check - fluxes)/fluxes < 1e-6) 

        ax1[i].text(20, 0.2, r"$\lambda=" + "{:4.0f}".format(wl) + r"\,\mathrm{A}$", fontsize=FS)
        ax1[i].plot(fluxes, compl, color=col, label="{:4.0f}".format(wl))
        ax1[i].plot(fluxes, compl_hdfs, "--", color=col)
        ax1[i].axhline(1.0, color="k", linestyle=":")

        if i > 3:
            ax1[i].set_xlabel(r"flux [10$^{-17}$ erg/s/cm$^2$]",
                              fontsize=FS)
 
        ax1[7].plot(fluxes, (compl_hdfs - compl)/compl, color=col, label=None)

        ax2.plot(1e17*bins, noise_hist, color=col, label="{:4.0f}".format(wl))


    ax1[0].set_ylabel("Completeness", fontsize=FS)
    ax1[4].set_ylabel("Completeness", fontsize=FS)

    ax1[7].set_ylabel("(Model - Sim.)/Sim.", fontsize=FS)
    ax1[7].set_xlabel(r"flux [10$^{-17}$ erg/s/cm$^2$]",
                      fontsize=FS)
 
    ax1[7].axhline(0.0, color="k", linestyle=":")
    ax1[7].set_ylim(-1.0, 1.0)
    fig1.subplots_adjust(left=0.07, right=0.95, 
                         wspace=0.35, top=0.98)
    #plt.tight_layout()

    ax2.set_ylabel("N pixels", fontsize=FS) 
    ax2.set_xlabel("Noise [10^{-17} erg/s/cm^{-2}]", fontsize=FS) 
    ax2.legend(frameon=False, prop={"size" : FS})

    return None
    #return alphas_fit

if __name__ == "__main__":

    waves, f50, compl_curves, fluxes = read_karl_file("fcor.use")
    
    with open("calshots", 'r') as fp:
        rshots = fp.readlines()
    
    shots = [x.strip() for x in rshots]

    remeasure = False 
    if remeasure: 
        wlcube, sigma_cube = return_biweight_one_sigma(shots)
        savetxt("average_one_sigma_cal.txt", transpose([wlcube, sigma_cube]))
    else:
        wlcube, sigma_cube = loadtxt("average_one_sigma_cal.txt").T


    #plt.figure(figsize=(11,4))
    #plot_f50(waves, f50, wlcube, sigma_cube)
    #plt.xlim(3600, 5500.)
    #plt.ylim(0, 30)
    #plt.legend(frameon=False, prop={"size" : 12.0})
    #plt.tight_layout()
    
    plot_completeness_curves(waves, f50, fluxes, compl_curves, shots)
    plt.tight_layout()
    plt.show()
    #plt.savefig("curves_with_model.png")    

    #plt.figure()
    #plt.plot(waves, alphas, "k*")
    #plt.xlabel("Wavelength")
    #plt.ylabel("alpha")
    #plt.show()
    
    #plt.show()
