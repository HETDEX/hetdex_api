"""

Module that deals with storing and
extracting sensitivity cubes in HDF5
files

.. moduleauthor:: Daniel Farrow <dfarrow@mpe.mpg.de>

"""

import logging
from os.path import isfile, join
from re import compile
import tables as tb
from numpy.ma import MaskedArray
from numpy import (linspace, mean, arange, array, histogram, zeros, std, 
                   array, sum, pad, ceil, floor, nan)
from astropy.stats import biweight_location
from astropy.table import Table
from hetdex_api.config import HDRconfig
from hetdex_api.flux_limits.sensitivity_cube import SensitivityCube

_logger = logging.getLogger()
_logger.setLevel("WARNING")
hndlr = logging.StreamHandler()
_logger.addHandler(hndlr)


class FileExists(Exception):
    pass


class NoFluxLimsAvailable(Exception):
    pass


def return_sensitivity_hdf_path(datevobs, release="hdr2.1", 
                                return_mask_fn = False):
    """
    Return the full file path
    to a HDF5 container of sensitivity     
    cubes

    Parameters
    ----------
    datevobs : str
        the date of the observation
    release : str (optional)
        the name of the release,
        default: hdr2
    Returns
    -------
    path_to_hdf5 : str
        the absolute path to
        the file

    Raises
    ------
    NoFluxLimsAvailable :
        Returned when flux limits
        not found for that shot

    """
    config = HDRconfig(release)
    flim_dir = config.flim_dir 
    mask_dir = config.flimmask

    path_to_hdf5s = join(flim_dir, "{:s}_sensitivity_cube.h5".format(datevobs))

    if isfile(path_to_hdf5s):
        if return_mask_fn:
            mask_fn = join(mask_dir, "{:s}_mask.h5".format(datevobs))
            return path_to_hdf5s, mask_fn
        else:
            return path_to_hdf5s
    else:
        raise NoFluxLimsAvailable(
            "Cannot find flux limits for that shot!"
            " Tried here {:s}".format(path_to_hdf5s)
        )


class SensitivityCubeHDF5Container(object):
    """
    Handle accessing and writing the sensitivity
    cubes stored together in an HDF5 container.
    Arguments passed to the init are passed on
    to the open_file function of tables.    

    Parameters
    ----------
    filename : string
        the filename of the HDF5
    flim_model : string (optional)
        specifies the flux limit model
        to use in the sensitivity
        cubes. If None assume the latest
        model (default).
    aper_corr : float (optional)
        aperture correction to
        apply to flux limits. Default
        is 1.0, if None then read 
        correction from header.
    mask_filename : str (optional)
        optional mask to apply, path to
        one of the _mask.h5 files

    Attributes
    ----------
    h5file : tables:File 
        the tables File object

    """

    def __init__(self, filename, mode="r", flim_model=None, aper_corr=1.0, 
                 mask_filename=None, verbose=True, **kwargs):

        if (mode == "w") and isfile(filename):
            raise FileExists("Error! Output file {:s} exists!".format(filename))

        self.verbose = verbose

        # Filter for compression, set complevel=4 as higher levels only
        # give a few extra MB per file
        # fletcher32 error checking and zlib
        self.compress_filter = tb.Filters(complevel=4, fletcher32=True, complib="zlib")
        self.h5file = tb.open_file(
            filename, mode=mode, filters=self.compress_filter, **kwargs
        )
        if mask_filename:
            self.h5mask = tb.open_file(mask_filename)
        else:
            self.h5mask = None

        self.filename = filename
        self.flim_model = flim_model
        self.aper_corr = aper_corr

    def add_sensitivity_cube(self, datevshot, ifuslot, scube, flush=False):
        """
        Add a sensitivity cube to the HDF5 file

        Parameters
        ----------
        datevshot : str
            the shot the IFU belongs 
            to
        ifuslot : str
            the IFU identifier
        scube : SensitivityCube
            the sensitivity cube object to
            add
        """

        # Try to get this shot from the database
        # if it doesn't exist add it
        try:
            shot = self.h5file.get_node(self.h5file.root, datevshot)
        except tb.NoSuchNodeError:
            shot = self.h5file.create_group(self.h5file.root, datevshot)

        # Store this in a compressible array
        if type(scube.sigmas) == MaskedArray.__class__:
            array = self.h5file.create_carray(
                               shot, ifuslot, obj=scube.sigmas.data, 
                               title="1 sigma noise")
        else:     
            array = self.h5file.create_carray(
                                shot, ifuslot, obj=scube.sigmas, 
                                title="1 sigma noise")
 
        #  Store what aperture correction has been applied
        array.attrs.aper_corr = scube.aper_corr

        # Store the header as an attribute
        array.attrs.header = scube.header
        array.attrs.wavelengths = scube.wavelengths
        array.attrs.alphas = scube.alphas

        # should always be 1-sigma stored from now on
        array.attrs.nsigma = 1.0

    def list_contents(self):
        """ List the contents of the HDF5 file """
        print(self.h5file)

    def itercubes(self, datevshot=None, cache_sim_interp=True):
        """ 
        Iterate over the IFUs 

        Parameters
        ----------
        datevshot : str (Optional)
            specify a datevshot or
            just take the first shot in
            the HDF5 file

        Yields
        ------
        ifuslot : str
           the IFU slot of the 
           returned IFU
        scube : SensitivityCube
           a sensitivity cube
           object

        """
        # Use first shot if dateshot not specified
        if not datevshot:
            shots = self.h5file.list_nodes(self.h5file.root)
            nshots = len(shots)
            if nshots > 1:
                _logger.warn(
                    """Datevshot not specified but multiple shots in file!
                                Using first in file."""
                )

            shot = shots[0]
        else:
            shot = self.h5file.get_node(self.h5file.root, name=datevshot)

        first = True
        warn = True
        for ifu in shot:

            # Extract the data we need for a sensitivity cube
            header = ifu.attrs.header
            wavelengths = ifu.attrs.wavelengths
            alphas = ifu.attrs.alphas
            sigmas = ifu.read() / ifu.attrs.aper_corr

            if first:
                verbose = self.verbose
                first = False
            else:
                verbose = False

            if self.h5mask:
                mask = self.h5mask.get_node(self.h5mask.root.Mask, 
                                            name=ifu.name).read()
            else:
                mask = None

            try:
                nsigma = ifu.attrs.nsigma
            except AttributeError as e:
                if self.flim_model == "hdr1":
                    nsigma = 6.0
                else:
                    nsigma = 1.0
                if warn:
                    print("No nsigma found, assuming nsigma={:2.1f} (warning will not repeat)".format(nsigma))
                    warn = False

            # XXX HACK HACK HACK to change alpha
            #alphas = [-1.9, -1.9]

            yield ifu.name, SensitivityCube(sigmas, header, wavelengths, alphas, 
                                            flim_model=self.flim_model,
                                            aper_corr=self.aper_corr, 
                                            nsigma=nsigma, mask=mask, 
                                            verbose=verbose, 
                                            cache_sim_interp=cache_sim_interp)

    def extract_ifu_sensitivity_cube(self, ifuslot, datevshot=None, 
                                     cache_sim_interp=True):
        """
        Extract the sensitivity cube
        from IFU (ifuslot). If multiple
        shots are saved in the file specify
        which to use with datevshot

        Parameters
        ----------
        ifuslot : string
            the IFU slot to extract
        datevshot : string (Optional)
            the datevshot if multiple
            shots are stored in the 
            HDF5. If None then test
            that one shot is present 
            and return the IFU
            for that

        Returns
        -------
        scube : hetdex_api.flux_limits.sensitivity_cube:SensitivityCube
            the sensitivity cube
        """

        # Use first shot if dateshot not specified
        if not datevshot:
            shots = self.h5file.list_nodes(self.h5file.root)
            nshots = len(shots)
            if nshots > 1:
                _logger.warn(
                    """Datevshot not specified but multiple shots in file!
                                Using first in file."""
                )

            shot = shots[0]
        else:
            shot = self.h5file.get_node(self.h5file.root, name=datevshot)


        if self.h5mask:
            mask = self.h5mask.get_node(self.h5mask.root.Mask, 
                                        name=ifuslot).read()
        else:
            mask = None

        # Now get the desired IFU
        ifu = self.h5file.get_node(shot, name=ifuslot)

        # Extract the data we need for a sensitivity cube
        header = ifu.attrs.header
        wavelengths = ifu.attrs.wavelengths
        alphas = ifu.attrs.alphas

        # Remove any aperture correction
        sigmas = ifu.read() / ifu.attrs.aper_corr

        try:
            nsigma = ifu.attrs.nsigma
        except AttributeError as e:
            if self.flim_model == "hdr1":
                nsigma = 6.0
            else:
                nsigma = 1.0
            print("No nsigma found, assuming nsigma={:2.1f} ".format(nsigma))

        # Force apcor to be 1.0 here, so we don't double count it
        return SensitivityCube(sigmas, header, wavelengths, alphas, 
                               nsigma=nsigma, flim_model=self.flim_model,
                               aper_corr=self.aper_corr, mask=mask,
                               cache_sim_interp=cache_sim_interp)

    def return_shot_completeness(self, flux, lambda_low, lambda_high, sncut,
                                 bin_edges = None, sigma_clip = False):
        """
        Parameters
        ----------
        flux : array
            flux to return completeness
            for
        lambda_high, lambda_low : float
            wavelength range
        sncut : float
            S/N cut to return completness for
        bin_edges : array (optional)
            if passed, also return the number of
            noise pixels in the supplied 
            bin_edges
        """
        
        means = []
        vals_all = []
        compl_all = []
        noise_hists = []

        for ifuslot, scube in self.itercubes():
            if type(bin_edges) != type(None):
                compl, vals = scube.return_wlslice_completeness(flux,
                                                                lambda_low, lambda_high,
                                                                sncut, return_vals=True)

                if len(compl) > 0:
                    vals_all.extend(vals)
                    means.append(mean(vals))               
                    noise_hists.append(histogram(vals, bins=bin_edges)[0])
            else:
                compl = scube.return_wlslice_completeness(flux,
                                                          lambda_low, lambda_high,
                                                          sncut)

            if len(compl)>0:
                compl_all.append(compl)
            else:
                print("Warning! No good data for {:s}, skipping".format(ifuslot))

        # Remove bad IFUs
        supermean = mean(means)
        superstd = std(means)
        means = array(means)
        compl_all = array(compl_all)
        noise_hists = array(noise_hists)
        print(superstd)
       
        if sigma_clip: 
            compl_all = compl_all[means < supermean + 2.*superstd]
            nbinned = sum(noise_hists[means < supermean + 2.*superstd], axis=0)
        else:
            nbinned = sum(noise_hists, axis=0)

        if type(bin_edges) != type(None):
            return mean(compl_all, axis=0), nbinned
        else:
            return mean(compl_all, axis=0)

 
    def flush(self):
        """ Write all alive leaves to disk """
        self.h5file.flush()

    def close(self):
        """ Close the file and destroy the object """
        self.h5file.close()
        _logger.info("Closed {:s}".format(self.filename))
        if self.h5mask:
            self.h5mask.close()

    def __enter__(self):
        """ Added to support using the `with` statement """
        return self

    def __exit__(self, type_, value, traceback):
        """ Support tidying up after using the `with` statement """
        self.close()

def shot_completeness_plot(args=None):

    import argparse
    import matplotlib as mpl
    mpl.use("Qt4Agg")
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser(
        description="Plot completeness curves in a lambda range",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument('h5file')
    parser.add_argument('--sim-results', help="Filename of sim results to overplot", default=None)
    parser.add_argument('--mask', help="Mask filename, if not passed assume no mask.", default=None)
    parser.add_argument('lambda_low', type=float, help="Lower wavelength limit (A)")
    parser.add_argument('lambda_high', type=float, help="Upper wavelength limit (A)")
    parser.add_argument('sncut', type=float, help="Signal to noise cut")
    parser.add_argument('--flux_low', type=float, default="1e-17")
    parser.add_argument('--flux_high', type=float, default="3e-16")
    opts = parser.parse_args()

    if opts.sim_results:
        table = Table.read(opts.sim_results, format="ascii")
        plt.plot(table["flux"], table["complete"], "r*", label="Sim results")
 
    h5s = SensitivityCubeHDF5Container(opts.h5file, mask_filename = opts.mask) 
 
    fluxes = linspace(opts.flux_low, opts.flux_high, 1000)
    compl = h5s.return_shot_completeness(fluxes, opts.lambda_low, 
                                         opts.lambda_high, opts.sncut)

    plt.plot(fluxes*1e17, compl, "k-", label="Model S/N > {:3.2f}".format(opts.sncut))

    plt.title("{:5.1f}<wavelength (A)<{:5.1f}".format(opts.lambda_low,  opts.lambda_high))
    plt.legend()
    plt.ylabel("Completeness")
    plt.xlabel("Fluxes [10$^{-17}$ erg/s/cm2]")
    plt.show()



def add_sensitivity_cube_to_hdf5(args=None):
    """
    Command line tool to add a sensitivity cube(s) 
    to an HDF5 containter

    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Add sensitivty cubes to HDF5 container",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--regex",
        default=".*(2[0-9]{7}v[0-9]{3})_[0-9]{3}_([0-9]{3})",
        help="""Regex with two capture groups, the first for datevshot the second 
                                for IFU slot""",
    )

    parser.add_argument("--append", action="store_true", help="Append to existing HDF5")

    parser.add_argument(
        "--alpha", type=float, help="Alpha for Fleming function", default=-3.5
    )

    parser.add_argument(
        "fits_files",
        type=str,
        nargs="+",
        help="Files to add, filename must follow convention set out in --regex option",
    )

    parser.add_argument("hdf5_out", type=str, help="HDF5 container to add to")

    opts = parser.parse_args(args=args)

    if opts.append:
        hdfcont = SensitivityCubeHDF5Container(opts.hdf5_out, mode="a")
    else:
        hdfcont = SensitivityCubeHDF5Container(opts.hdf5_out, mode="w")

    # Compile the regex so it's faster
    regex = compile(opts.regex)

    for fn in opts.fits_files:

        m = regex.match(fn)
        datevshot = m.group(1)
        ifuslot = m.group(2)
        _logger.info(
            "Inserting {:s} with dateshot {:s} and IFU slot {:s}".format(
                fn, datevshot, ifuslot
            )
        )
        scube = SensitivityCube.from_file(
            fn, [3500.0, 5500.0], [opts.alpha, opts.alpha]
        )
        hdfcont.add_sensitivity_cube("virus_" + datevshot, "ifuslot_" + ifuslot, scube)

    hdfcont.close()


def extract_sensitivity_cube(args=None):
    """
    Extract a sensitivity cube from an HDF5 file

    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Add sensitivty cubes to HDF5 container",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("h5file", type=str, help="HDF5 container to extract from")
    parser.add_argument(
        "--datevshot",
        type=str,
        help="Datevshot to extract [DATE]v[SHOT]. If None extract first/only shot in file.",
    )
    parser.add_argument("ifuslot", help="IFUSLOT of IFU to extract")
    parser.add_argument("outfile", help="Output filename")
    opts = parser.parse_args(args=args)

    with SensitivityCubeHDF5Container(opts.h5file, mode="r") as hdfcont:
        if opts.datevshot:
            scube = hdfcont.extract_ifu_sensitivity_cube(
                "ifuslot_" + opts.ifuslot, datevshot="virus_" + opts.datevshot
            )
        else:
            scube = hdfcont.extract_ifu_sensitivity_cube("ifuslot_" + opts.ifuslot)

    scube.write(opts.outfile)


def return_biweight_one_sigma(shots, pixlo=5, pixhi=25, wave_bins=None):
    """
    Return the biweight midvariance of the 
    1 sigma values of a series of shots. Trim
    off edge values
   
    Parameters
    ----------
    shots : list
        a list of shots to loop
        over
    pixlo, pixhi : int (optional)
        limits for the 2D ra/dec
        pixel array when computing 
        the values. Default is 5,25.
    """
    sigmas_all = []
    for shot in shots:
        fn, fnmask = return_sensitivity_hdf_path(shot, 
                                                 return_mask_fn = True)

        print(fn)
        with SensitivityCubeHDF5Container(fn, mask_filename = fnmask) as hdf:

            first = True

            for ifuslot, scube in hdf.itercubes():

                 # assume wavelenghth WCS the same for whole HDF
                 if (type(wave_bins) != type(None)) and first:
                     junkx, junky, wlo = scube.wcs.all_world2pix(0., 0., wave_bins[:-1], 0)
                     junkx, junky, whi = scube.wcs.all_world2pix(0., 0., wave_bins[1:], 0)
                     first = False
                     maxsize = int(max(floor(whi) - ceil(wlo))*(pixhi - pixlo)*(pixhi - pixlo))

                 if type(wave_bins) != type(None):
                     sigmas = []
                     for lo, high in zip(wlo, whi):
                         sigmas_slice = scube.sigmas.filled(nan)[int(ceil(lo)):int(floor(high)), 
                                                                 pixlo:pixhi, pixlo:pixhi]
                         sigmas_slice = sigmas_slice.reshape(sigmas_slice.shape[0]*sigmas_slice.shape[1]*sigmas_slice.shape[2])
                         
                         # pad with NaN to keep consistent size
                         if len(sigmas_slice) < maxsize:
                             sigmas_slice = pad(sigmas_slice, (0, maxsize - len(sigmas_slice)), 
                                                'empty')

                         sigmas.append(sigmas_slice)
                 else:
                     sigmas = scube.sigmas.filled(nan)[:, pixlo:pixhi, pixlo:pixhi]
                     sigmas = sigmas.reshape(sigmas.shape[0], sigmas.shape[1]*sigmas.shape[2])

                 sigmas_all.extend(array(sigmas).T)

            # assume wavelenghth WCS the same for whole HDF
            if type(wave_bins) == type(None):
                 pixels = arange(scube.sigmas.shape[0])
                 junkras, junkdecs, waves = scube.wcs.all_pix2world(0*pixels, 0*pixels,  
                                                                    pixels, 0)
            else:
                 waves = 0.5*(wave_bins[1:] + wave_bins[:-1])
 
    biwt = biweight_location(array(sigmas_all), axis=0, ignore_nan=True)

    return waves, biwt 

if __name__ == "__main__":
    extract_sensitivity_cube()
    ## add_sensitivity_cube_to_hdf5()
