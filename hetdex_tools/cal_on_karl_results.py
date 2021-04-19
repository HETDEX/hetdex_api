"""

Derive a scaling between the values in the
flux limit cubes and the completeness
simulations Karl Gebhardt ran. Also 
produce curves of completeness versus
flux from the flux limit cubes, and
compare to the 

Daniel Farrow (MPE) 2021

"""
FS=16.0
from numpy import (loadtxt, savetxt, transpose, interp, sqrt, 
                   mean, linspace, zeros, array, polyfit, polyval)
import matplotlib as mpl
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['mathtext.fontset'] = 'dejavuserif'
mpl.rcParams["xtick.labelsize"] = 15.0
mpl.rcParams["ytick.labelsize"] = 15.0
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import TABLEAU_COLORS
from scipy.special import erf
from scipy.optimize import leastsq
from scipy.signal import convolve
from hetdex_api.flux_limits.hdf5_sensitivity_cubes import (return_biweight_one_sigma, return_sensitivity_hdf_path,
                                                           SensitivityCubeHDF5Container)
from hetdex_api.flux_limits.sensitivity_cube import fleming_function
from hetdex_api.flux_limits.flim_models import (hdr2pt1_f50_from_noise, SimulationInterpolator, 
                                                hdr2pt1pt1_f50_from_noise, hdr2pt1pt3_f50_from_noise, 
                                                read_karl_file)

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


def plot_f50(waves, f50, wlcube, sigma_cube, sncut, col, marker="*",
             plot_cube = True):
    """
    Plot and fit polynomials to the 50% compleness values 
    and the sigma values from the flux limit cubes

    Parameters
    ----------
    waves, f50 : array
        arrays of wavelength and
        50% completeness from the
        simulations
    wlcube, sigma_cube : array
        average noise and wavelengths
        from the cubes
    sncut : float
        signal to noise cut
    col : str
        colour of markers and lines
    plot_cube : logical (optional)
        plot the noise from the cube,
        scaled by the best fit noise to
        50% completeness value 
    """
    #plt.plot(wlcube, 1e17*hdr2pt1_f50_from_noise(sigma_cube, sncut), "--", 
    #         label="HDR2.1 scaled sigmas", color=col)
    #plt.plot(wlcube, 5e17*sigma_cube, "k--", 
    #         label="5$\sigma$ from cubes")

    # fit a new scaling
    sigma_f50_pos = interp(waves, wlcube, sigma_cube)
    best, info = leastsq(func_f50s, [5.0], args=(f50[2:], 1e17*sigma_f50_pos[2:]))
    print(sncut, best[0], best[0]/hdr2pt1_f50_from_noise(1.0, sncut))

    if plot_cube:
        plt.plot(wlcube, 1e17*best[0]*sigma_cube, "-", color=col,
                 label="New scaling ({:3.2f}".format(best[0]) + "$\sigma$)")

    plt.plot(waves, f50, marker, label="Karl sims", 
             markersize=14.0, color=col)

    plt.xlabel("Wavelength [AA]", 
               fontsize=FS)
    plt.ylabel(r"50% flux [10$^{-17}$ erg/s/cm$^2$]", 
               fontsize=FS)

    return best[0]

def return_completeness_from_shots(shots, fluxes, lambda_low, lambda_high,
                                   sncut, rescale_50 = None):
   """
   Return the completeness over a list
   of shots

   """

   bin_edges = linspace(0, 1e-16, 500)
   bins = 0.5*(bin_edges[:-1] + bin_edges[1:])
   nbinned = zeros(len(bin_edges) - 1)

   cube_comp_all = []
   for shot in shots:
       fn, fnmask = return_sensitivity_hdf_path(shot, return_mask_fn = True)

       with SensitivityCubeHDF5Container(fn, mask_filename = fnmask) as hdf:
           cube_compl, tnbinned = hdf.return_shot_completeness(fluxes*1e-17, lambda_low, 
                                                               lambda_high, sncut,
                                                               bin_edges = bin_edges)

           # rescale to units of 1/f50
           # issues with 
           if rescale_50:
               f50 = interp(0.5, cube_compl/cube_compl[-1], fluxes)
               cube_compl = interp(fluxes/rescale_50, fluxes/f50, cube_compl)  

           nbinned += tnbinned

       cube_comp_all.append(cube_compl)
 
   return mean(cube_comp_all, axis=0), bins, nbinned


def produce_completeness_curves(shots, fluxes, f50, sncut, waves, 
                                write_noise=False):
    """

    Produce curves of the completeness versus flux for a wavelength region,
    across a range of shots. Write them to files.

    """
    wlbsize = waves[2] - waves[1]
    print("Wave bin size: {:f}".format(wlbsize))
    wlbsize = 0.0


    # Loop over wavelength ranges, rescale to match 
    # the simulation 50% completeness
    for tf50, wl in zip(f50, waves):

        wllo = wl - wlbsize/2.0
        if wllo < 3500:
            wllo = 3500

        compl_hdfs, bins, noise_hist = return_completeness_from_shots(shots[:1], fluxes, 
                                                                      wllo, wl + wlbsize/2.0,
                                                                      sncut,
                                                                      rescale_50 = tf50)

        savetxt("av_cube_completeness_cal_{:4.1f}_{:2.1f}_rescale.txt".format(wl, sncut), 
                transpose([fluxes, compl_hdfs]))

        if write_noise:
            savetxt("noise_hist_{:4.1f}.txt".format(wl), 
                    transpose([bins, noise_hist])) 



def plot_completeness_curves(waves, f50, fluxes, compl_curves, shots, sncut):
    
    FS = 13.0
    alpha_fix = -3.5 

    alphas_fit = []
    colors =['r', 'g', 'b', 'k', 'm', 'c', 'y']

    fig1, ax1 = plt.subplots(nrows = 2, ncols = 4, 
                             figsize=(14.0, 7))
    fig2, ax2 = plt.subplots(nrows = 1, ncols = 1, figsize=(8, 4))
    ax1 = ax1.flatten()

    for i, (col, wl, tf50, compl)  in enumerate(zip(colors, waves, f50, compl_curves)):

        #sinterp = SimulationInterpolator("fcor.use")
        #model = compl[-2]*sinterp(fluxes, tf50) 

        fluxes_check, compl_hdfs = loadtxt("av_cube_completeness_cal_{:4.1f}_{:2.1f}_rescale.txt".format(wl, sncut), 
                                           unpack=True)
        bins, noise_hist = loadtxt("noise_hist_{:4.1f}.txt".format(wl),
                                   unpack = True)

        assert all((fluxes_check - fluxes)/fluxes < 1e-6) 

        ax1[i].text(20, 0.2, r"$\lambda=" + "{:4.0f}".format(wl) + r"\,\mathrm{A}$", fontsize=FS)
        ax1[i].plot(fluxes, compl/compl[-1], color=col, label="{:4.0f}".format(wl))
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
    plt.tight_layout()

    ax2.set_ylabel("N pixels", fontsize=FS) 
    ax2.set_xlabel(r"Noise [$10^{-17}$ erg/s/cm$^{-2}$]", fontsize=FS) 
    ax2.legend(frameon=False, prop={"size" : FS})

    return None

def plot_sn_scaling(snlist, scaling):
    """
    Plot the scaling between the cubes and the
    50% flux limit, as derived from sims

    Parameters
    ----------
    snlist : array
        list of S/N cuts
    scaling : array
        scaling factor between
        noise and 50% completeness
        as a function of wavelength
    """

    plt.figure() 
    params = polyfit(snlist, scaling, 3)
    snplot = linspace(snlist[0] - 0.5, snlist[-1] + 0.5, 200) 
    model = polyval(params, snplot)

    model_hdr213 = hdr2pt1pt3_f50_from_noise(1.0, snplot)
    model_hdr211 = hdr2pt1pt1_f50_from_noise(1.0, snplot)

    plt.plot(snlist, scaling, "k*", markersize=12.0,  label="Simulations")
    plt.plot(snplot, model, "r-", label="Polynomial fit")
    plt.plot(snplot, model_hdr213, "k--", label="HDR 2.1.3 Model")

    print(params)

    plt.xlabel("S/N cut", fontsize=FS)
    plt.ylabel("Cube to Simulation 50\% completeness flux scaling",
               fontsize=FS)
    plt.legend(frameon=False, prop={"size" : FS})


def measure_f50_scaling(folder, snlist, wlcube, sigma_cube):
    """
    Loop over the sn?.?.use files deriving a scaling
    between the 50% flux and the noise in the
    flux limit cubes
    
    Parameters
    ----------
    folder : str
        path to folder of sn?.?.use files
    snlist : array
        list of S/N cuts to consider
    wlcube, sigma_cube : array
        arrays of wavelength and
        average noise from the cube - 
        to derive scalings between
        noise and completeness
    """
    plt.figure(figsize=(11*2,4*2.2))

    cols = []
    scaling = []
    for sncut, col in zip(snlist, TABLEAU_COLORS):

        waves, f50, compl_curves, fluxes = read_karl_file(folder + "sn{:2.1f}.use".format(sncut))

        scale = plot_f50(waves, f50, wlcube, sigma_cube, sncut, col)

        scaling.append(scale)
        cols.append(col)

    plt.xlim(3600, 5500.)
    plt.ylim(0, 45)
    plt.legend(frameon=False, prop={"size" : 12.0})
    plt.tight_layout()

    cmap = mpl.colors.ListedColormap(cols)
    bounds = range(len(snlist) + 1)
    cens = 0.5*(array(bounds[1:]) + array(bounds[:-1]))

    cax = plt.gca().inset_axes([4500, 33.5, 820, 1.5],
                               transform=plt.gca().transData) 

    cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
                                    boundaries=bounds,
                                    ticks=cens, 
                                    orientation="horizontal")
    cb.ax.set_xticklabels(snlist)

    cb.set_label("S/N cut", fontsize=FS)   


    legend_elements = [] 
    legend_elements.append(mpl.lines.Line2D([0], [0], color='k',
                                            lw=2.5, linestyle="-", 
                                            label="Scaling to Karl's simulation"))
    legend_elements.append(mpl.lines.Line2D([0], [0], color='k',
                                            linestyle="", marker="*", 
                                            markersize=14.0,
                                            label="Karl's source simulations"))
    plt.gca().legend(handles=legend_elements, frameon=False, loc="lower left",
                     prop={"size" : FS})

    plt.show()
 
    return scaling

if __name__ == "__main__":
   
    folder = "/data/hetdex/u/dfarrow/hetdex_data/hdr2.1/reduction/flim/snfiles/"
    with open("calshots", 'r') as fp:
        rshots = fp.readlines()
    
    shots = [x.strip() for x in rshots]

    snlist = [4.8, 5.0, 5.5, 6.0, 6.5, 7.0]

    remeasure = False 
    if remeasure: 
        wlcube, sigma_cube = return_biweight_one_sigma(shots)
        savetxt("average_one_sigma_cal.txt", transpose([wlcube, sigma_cube]))
    else:
        wlcube, sigma_cube = loadtxt("average_one_sigma_cal.txt").T


    scaling = measure_f50_scaling(folder, snlist, wlcube, sigma_cube)
    plot_sn_scaling(snlist, scaling)

    plt.show()

    #first = True
    #for sncut in snlist:
    #    waves, f50, compl_curves_base, fluxes_base = read_karl_file(folder + "sn{:2.1f}.use".format(sncut)) 
    #  
    #    # Add more flux bins to help with interpolation later
    #    fluxes = linspace(fluxes_base[0], 1.4*fluxes_base[-1], 2000) 
    #    compl_curves = []
    #    for c in compl_curves_base:        
    #        compl_curves.append(interp(fluxes, fluxes_base, c))
    #        

    #    #produce_completeness_curves(shots, fluxes, f50, sncut, waves, write_noise=first)
    #    first = False
    #    plot_completeness_curves(waves, f50, fluxes, compl_curves, shots, sncut)
  
    #plt.tight_layout()
    #plt.show()
  
    #plt.savefig("curves_with_model.png")    

    #plt.figure()
    #plt.plot(waves, alphas, "k*")
    #plt.xlabel("Wavelength")
    #plt.ylabel("alpha")
    #plt.show()    
