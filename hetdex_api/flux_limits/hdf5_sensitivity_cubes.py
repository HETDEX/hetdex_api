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
from hetdex_api.config import HDRconfig
from hetdex_api.flux_limits.sensitivity_cube import SensitivityCube
from numpy.ma import MaskedArray

_logger = logging.getLogger()
_logger.setLevel("INFO")
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
        cubes either hdr1 or hdr2pt1.
        Default is hdr2pt1
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

    def __init__(self, filename, mode="r", flim_model="hdr2pt1", aper_corr=1.0, 
                 mask_filename = None, **kwargs):

        if (mode == "w") and isfile(filename):
            raise FileExists("Error! Output file {:s} exists!".format(filename))

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

    def itercubes(self, datevshot=None):
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

        warn = True
        for ifu in shot:

            # Extract the data we need for a sensitivity cube
            header = ifu.attrs.header
            wavelengths = ifu.attrs.wavelengths
            alphas = ifu.attrs.alphas
            sigmas = ifu.read() / ifu.attrs.aper_corr

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
                                            nsigma=nsigma, mask=mask)

    def extract_ifu_sensitivity_cube(self, ifuslot, datevshot=None):
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
                               aper_corr=self.aper_corr, mask=mask)

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


if __name__ == "__main__":
    extract_sensitivity_cube()
    # add_sensitivity_cube_to_hdf5()
