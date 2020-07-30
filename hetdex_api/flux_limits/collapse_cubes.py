"""

Module to combine and collapse sensitivity
cubes into single curves of wavelength and flux
versus completeness 


.. moduleauthor:: Daniel Farrow <dfarrow@mpe.mpg.de>

"""
import logging
import argparse
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
from numpy import (ones_like, zeros_like, sum, inf, zeros, linspace, interp, array, 
                   nanmean, nanmax, isfinite, nanmedian, isnan, sqrt)
from scipy.optimize import least_squares
from astropy.table import Table
from astropy.stats import biweight_location, biweight_midvariance
from hetdex_api.flux_limits.sensitivity_cube import fleming_function, SensitivityCube
from hetdex_api.flux_limits.hdf5_sensitivity_cubes import SensitivityCubeHDF5Container

def fleming_diff(alpha, fluxes_in, compl_in, f50):
    """ Compute residuals of a Fleming model fit """    

    pred = fleming_function(fluxes_in, f50, alpha)
   
    return compl_in - pred


def plot_collapsed_cube(f50vals, alphas, lambda_, fluxes, combined_cube, fn_plot):  
    """
    Plot the results of collapsing a cube in the RA, DEC directions
    """

    good_indices = (f50vals < 1.0)

    # Top panel flux at 50%
    plt.subplot(221)
    plt.plot(lambda_[good_indices], f50vals[good_indices]*1e16, 'k.')

    # Plot these as lower limits
    plt.errorbar(lambda_[~good_indices], max(f50vals[good_indices]*1e16)*ones_like(lambda_[~good_indices]), 
                 color='r', linestyle="", yerr=0.1, lolims=True)

    plt.xlabel("Wavelength (Angstrom)")
    plt.ylabel("Flux at 50% Completeness ($10^{-16}$ ergs/s/cm$^2$/A)")

    # Test of recovery versus flux, at a few wavelengths
    plt.subplot(222)
    for color, sample in zip(['r', 'g', 'b'], [100, 500, 800]):
        wl = lambda_[sample]
        compl_here = combined_cube[:, sample]

        # If we hit a bad fit don't plot
        if f50vals[sample] > 1.0:
            continue

        plt.plot(fluxes*1e16, compl_here, marker="*", color=color, label=None)
        model = fleming_function(fluxes, f50vals[sample], alphas[sample])
        plt.plot(fluxes*1e16, model, linestyle="--", color=color, label="{:4.1f} A".format(wl))

    plt.legend(frameon=False, loc="lower right")
    plt.xlabel("Flux ($10^{-16}$  erg/s/cm$^2$/A)")
    plt.ylabel("Detected Fraction")
    plt.xlim(fluxes[0]*1e16, 0.3*fluxes[-1]*1e16)
    plt.ylim(0.0, 1.05)

    # Plot of alpha versus wavelength
    plt.subplot(223)
    plt.plot(lambda_[good_indices], alphas[good_indices], 'k.')
    plt.xlabel("Wavelength (Angstrom)")
    plt.ylabel(r"$\alpha$", fontsize=15.0)

    plt.gcf().set_size_inches(10, 10)
    plt.subplots_adjust(bottom = 0.07, left = 0.08, top = 0.98, right=0.97)
    plt.savefig(fn_plot)

def return_flattened_slice(cube, lambda_):
    """

    Return a flattened slice at wavelength
    lambda_ from sensitivity cube

    Parameters
    ----------
    cube : pyhetdex.selfunc.sensitivity_cube:SensitivityCube
        sensitivity cube to collapse 

    """

    # Look up the relevant iz for lambda
    ra, dec, junk = cube.wcs.all_pix2world(5.0, 5.0, 5.0, 0)
    ix, iy, iz = cube.wcs.all_world2pix(ra, dec, lambda_, 0)

    # Grab a wavelength slice and collapse it
    wavelength_slice = cube.f50vals[int(iz),:,:] 
    slice_flattened = wavelength_slice.flatten()

    return slice_flattened

def return_flattened_wlrange(cube, lambda_1, lambda_2):
    """

    Return a flattened cube between two
    wavelength cuts

    Parameters
    ----------
    cube : pyhetdex.selfunc.sensitivity_cube:SensitivityCube
        sensitivity cube to collapse 
    lambda1, lambda2 : float
        the wavelength range in Angstrom. lambda_2 > lamba_1

    """

    if lambda_2 < lambda_1:
        raise ValueError("The second wavelength must be bigger than the first!")

    # Look up the relevant iz for lambda
    ra, dec, junk = cube.wcs.all_pix2world(5.0, 5.0, 5.0, 0)
    ix, iy, iz1 = cube.wcs.all_world2pix(ra, dec, lambda_1, 0)
    ix, iy, iz2 = cube.wcs.all_world2pix(ra, dec, lambda_2, 0)

    # Grab a wavelength slice and collapse it
    wavelength_slice = cube.f50vals[int(iz1):int(iz2),:,:]
    slice_flattened = wavelength_slice.flatten()

    return slice_flattened



def return_biwt_cmd(args=None):
    """
    Command line tool to return the bi-weight
    of the flux limit at a certain wavelength
    from an HDF5 file

    """

    # Parse the arguments
    parser = argparse.ArgumentParser(description="Compute biweight flux limits from HDF5 file")
    parser.add_argument("--wlrange", nargs=2, default=[4500, 4600], type=int, help="Wavelength range to compute median over")
    parser.add_argument("--hist", action="store_true", help="Plot histograms of flux limit")
    parser.add_argument("--nkeep", type=int, default=10000, help="To remove the edge pixel only the nkeep deepest pixels"  
                                                                 " are considered (default 10000)")
    parser.add_argument("--hist_all", action="store_true", help="Plot histograms flux limit for all inputs")
    parser.add_argument("--fout", default=None, type=str, help="Ascii file to save results to")
    parser.add_argument("--fn-shot-average", default=None, type=str, help="Ascii file to append shot average flim to")
    parser.add_argument("files", help="HDF container(s) of sensitivity cubes", nargs='+')
    opts = parser.parse_args(args=args)
    
    print("Using wavelengths {:f} to {:f} AA".format(*opts.wlrange))
    
    # Loop over the files producing median and percentiles
    biwt_ls = []
    ifu = []
    biwt_vars =[]
    dateshot = []

    for fn in opts.files:

        # Open the HDF5 container
        with SensitivityCubeHDF5Container(fn) as hdfcont:

            # Loop over shots
            shots_groups = hdfcont.h5file.list_nodes(hdfcont.h5file.root) 
            for shot_group in shots_groups:
                str_datevshot = shot_group._v_name
                flims_shot = []

                # Loop over sensitivity cubes
                for ifu_name, scube in hdfcont.itercubes(datevshot=str_datevshot):
                    flims = return_flattened_wlrange(scube, opts.wlrange[0], opts.wlrange[1]) 
                    flims = flims[isfinite(flims)]
                    flims.sort()
                    
                    # Compute statistics and save numbers
                    # from this cube
                    flims_shot.extend(flims[:opts.nkeep])    
                    biwt_ls.append(biweight_location(flims[:opts.nkeep]))
                    biwt_vars.append(biweight_midvariance(flims[:opts.nkeep]))

                    # Save IFU and shot info
                    ifu.append(ifu_name.strip("ifuslot_"))
                    dateshot.append(str_datevshot.strip("/virus_"))

                if opts.fn_shot_average:
                    with open(opts.fn_shot_average, 'a') as fn:
                        fn.write("{:s} {:e} \n".format(str_datevshot.strip("virus_"), biweight_location(flims_shot)))

    table = Table([dateshot, ifu, biwt_ls, sqrt(biwt_vars)], names=["dateshot", "ifu", "biwt_loc", "sqrt_biwt_mdvar"])
        
    # Output or save
    if opts.fout:
        table.write(opts.fout)
    else:
        table.pprint(max_lines=-1)   


def return_spatially_collapsed_cube(cube, fluxes, min_compl=0.99):
    """
    Collapse the x,y axes of a sensitivity cube and return
    the probability of detecting a certain flux
    versus wavelength. Ignore pixels with value=0 in 
    cubes, and also very shallow regions (see min_compl).

    Parameters
    ----------
    cube : pyhetdex.selfunc.sensitivity_cube:SensitivityCube
        The sensitivity cube to spatially collapse
    fluxes : array
        The fluxes to return the completeness versus
        wavelength from

    Returns
    -------
    lambda_ : array
        the wavelengths
    npix : array
        The number of good flux limit pixels as a 
        function of lambda_, also the number of
        pixels averaged over
    compls : array
        2D array of the completeness at lambda_ 
        where the first axes is the fluxes and
        the second the wavelengths
    min_compl : float (Optional)
        The minimum required completeness at
        the brightest flux to be included in the 
        calculation.
    """ 

    nz, ny, nx = cube.f50vals.shape  

    # Generate wavelengths for all the pixels
    junk1, junk2, lambda_ = cube.wcs.all_pix2world(range(nz), range(nz), range(nz), 0)    

    # Generate a cube of alphas for the Fleming function
    alphas = cube.alpha_func(lambda_)
    alpha_cube = zeros_like(cube.f50vals)
    for i, alpha in enumerate(alphas, 0):
        alpha_cube[i, :, :] = alpha

    # The number of good flux limit pixels as a function of lambda_
    masked_indices = (cube.f50vals == inf)
   
    # Remove places so shallow they're not complete at max flux (otherwise curve never reaches 1.0)
    max_flux = fluxes[-1]
    compl_cube = fleming_function(max_flux, cube.f50vals, alpha_cube)

    #for c in compl_cube[~masked_indices].flatten():
    #    print(c, max_flux)

    too_shallow = compl_cube < min_compl

    # Indices to ignore
    bad_indices = masked_indices | too_shallow

    # Number of usable pixels
    npix = sum(~bad_indices, axis=(1,2))

    # Generate completeness for all the pixels
    compls = zeros((len(fluxes), nz))
    for i, flux in enumerate(fluxes, 0):
        compl_cube = fleming_function(flux, cube.f50vals, alpha_cube) 
 
        # Remove empty regions
        compl_cube[bad_indices] = 0.0

        # Compute mean of non-empty pixels
        compls[i, :] = sum(compl_cube, axis=(1, 2))/npix

    return lambda_, npix, compls

def compute_new_fleming_fits(lambda_, fluxes, compls):
    """
    Compute the 50% detection limit of a data cube, collapsed spaitally, as
    a function of wavelength. Do this by interpolating the completeness
    in flux bins. Then use that to fit the alpha value of the Fleming 
    parameterisation

    Parameters
    ----------
    lambda_  : array
        the wavelengths of the flux limits
    fluxes : array
        the fluxes correponding to the
        first axis of compls
    compls : array
        2D array of completeness versus flux
        and wavelength
    cube : pyhetdex.selfunc.sensitivity_cube:SensitivityCube
        cube to collapse and refit
  
    Returns
    -------
    f50vals : array
        the 50% completeness limits 
    alphas : array
        the new alpha values for the 
        Fleming+ (1995) parameterisation on
        completeness
   """
    warn = False
    f50vals = zeros_like(lambda_)
    alphas = zeros_like(lambda_)

    for i, wl in enumerate(lambda_):

        # Find 50% completeness limit
        # First index flux, second lambda
        f50vals[i] = interp(0.5, compls[:, i], fluxes)

        # Catch bad wavelength slices
        if not all(isfinite(compls[:, i])):
            warn = True
            f50vals[i] = 999.9
            alphas[i] = 999.9 
        else:
            # Refit the Fleming function to the combined completness
            opr = least_squares(fleming_diff, [-3.5], args=(fluxes, compls[:, i], f50vals[i]))
            alphas[i] = opr.x[0]

    if warn:
        logging.warning("At least some of the wavelength slices had invalid completeness measurements!")

    return f50vals, alphas



def collapse_datacubes_command(args=None):
    """
    Produce combined flux limit versus 
    wavelength estimates for sensitivity
    cubes passed on command line 
    """
    import argparse 

    # Command line options
    parser = argparse.ArgumentParser(description="""
                                                 Collapse the RA and DEC of a single or 
                                                 set of sensitivity cube(s) to
                                                 produce one file of 50% flux limit versus
                                                 wavelength.
                                                 """)

    parser.add_argument("--plot", type=str, help="Filename for optional plot", default="")

    parser.add_argument("--fmin", help="Minimum flux to consider when interpolating for 50% limit", 
                        type=float, default=1e-17)

    parser.add_argument("--fmax", help="""Maximum flux to consider when interpolating for 50% limit.
                                          Regions not 99% at this flux ignored!""",
                        type=float, default=1e-15)

    parser.add_argument("--nbins", help="Number of flux bin to use when interpolating to measure 50% limit",
                        type=int, default=100)

    parser.add_argument("--alpha", help="The alpha of the Fleming function fit (default=-3.5)",  
                        default=-3.5, type=float)

    parser.add_argument("scubes", nargs='+', help="The sensitivity cube(s) you want to collapse and combine")
    parser.add_argument("output", help="(Ascii) file to output to", type=str)
    opts = parser.parse_args(args=args)

    # Stores a list of the completeness versus flux and lambda for
    # all the cubes
    compl_cube_list = []

    # Flux bins to compute completeness in
    fluxes = linspace(opts.fmin, opts.fmax, opts.nbins)
 

    # A list of the number of good pixels in each cube as 
    # a function of wavelength
    npixes_list = []
    for cube_name in opts.scubes:

        # Currently fix the alpha versus wavelength
        tscube = SensitivityCube.from_file(cube_name, [3500.0, 5500.0], [opts.alpha, opts.alpha])
        lambda_, npix, compls = return_spatially_collapsed_cube(tscube, fluxes)

        compl_cube_list.append(compls)
        npixes_list.append(npix)

    # Produce a combined cube of completness versus flux and lambda, 
    # weighted by the number of good pixels in each cube
    combined_cube = compl_cube_list[0]*npixes_list[0]
    for npix, compl_cube in zip(npixes_list[1:], compl_cube_list[1:]):
        combined_cube += npix*compl_cube
  
    # This should give the correct completeness versus lambda and flux, ignoring empty pixels 
    combined_cube /= sum(array(npixes_list), axis=0)

    # This should give is the 50% flux limits and the new alpha values
    f50vals, alphas = compute_new_fleming_fits(lambda_, fluxes, combined_cube) 

    if opts.plot:
        plot_collapsed_cube(f50vals, alphas, lambda_, fluxes, combined_cube, opts.plot) 

    # Write out
    table = Table([lambda_, f50vals, alphas], names=["wavelength", "f50", "alpha"])
    table.write(opts.output, format="ascii")


