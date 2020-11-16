# -*- coding: utf-8 -*-
"""

Functions to interact with the shot HDF5 file and the Fibers Class. Requires a
shotid (interger: 20180123009) or datevobs (str: eg. '20180123v009')

author = Erin Mentuch Cooper
"""

import os.path as op
import re
import tables as tb
import numpy as np

import warnings
import sys

from astropy.table import Table
import astropy.units as u
from astropy.coordinates import SkyCoord

from hetdex_api.config import HDRconfig

if not sys.warnoptions:
    warnings.simplefilter("ignore")

    
def open_shot_file(shotid, survey="hdr2.1"):
    """
    Open the H5 file for a shot. This is a global function that allows you
    to open an H5 file based on its shotid and data release. 

    Parameters
    ----------

    shotid
        either in the form of integer shotid (eg. 20180123009) or
        string datevobs (eg. '20180123v009')

    survey : string
            Data release you would like to load, i.e., 'HDR1', 'hdr2', 
            This is case insensitive.

    Example
    -------

    >>> from hetdex_api.shot import open_shot_file
    >>> fileh = open_shot_file(20180123009)

    or

    >>> fileh = open_shot_file('20180123v009')

    """

    global config
    config = HDRconfig(survey=survey.lower())

    if re.search("v", str(shotid)):
        file = op.join(config.data_dir, str(shotid) + ".h5")
    else:
        file = op.join(
            config.data_dir, str(shotid)[0:8] + "v" + str(shotid)[8:11] + ".h5"
        )
    fileh = tb.open_file(file, "r")
    return fileh


class Fibers:
    def __init__(self, shot, survey="hdr2.1"):
        """
        Initialize Fibers Class

        This creates an astropy coordinates array of all
        fiber coordinates. Because of the large number of fibers
        for an individual shot this will take some time, but makes
        positional querying much easier.

        This will also initiate the wave_rect attribute which is
        an array of rectified wavelengths corresponding the the
        'calfib', 'calfibe', 'Amp2Amp', and 'Throughput' datasets

        Parameters
        ----------
        shot
            either in the form of integer shotid (eg. 20180123009) or
            string datevobs (eg. 20180123v009)
        survey
            Data release you would like to load, i.e., 'HDR1', 'hdr2', 'hdr2.1'
            This is case insensitive.

        Attributes
        ----------
        hdfile
            h5 file handle for the specified shotid
        table
            h5 table of Fiber data
        coords
            astropy sky coordinates for each fiber
        wave_rect
            rectified wavelength for interpolated spectral data 'calfib',
            'calfibe'
        """

        self.hdfile = open_shot_file(shot, survey=survey)

        self.table = self.hdfile.root.Data.Fibers

        # Grab attributes from FiberIndex table if survey!='hdr1'

        if survey == "hdr1":
            colnames = self.hdfile.root.Data.Fibers.colnames
            for name in colnames:

                if isinstance(
                    getattr(self.hdfile.root.Data.Fibers.cols, name)[0], np.bytes_
                ):
                    setattr(
                        self,
                        name,
                        getattr(self.hdfile.root.Data.Fibers.cols, name)[:].astype(str),
                    )
                else:
                    setattr(
                        self, name, getattr(self.hdfile.root.Data.Fibers.cols, name)[:]
                    )
        else:
            colnames = self.hdfile.root.Data.FiberIndex.colnames
            for name in colnames:

                if isinstance(
                    getattr(self.hdfile.root.Data.FiberIndex.cols, name)[0], np.bytes_
                ):
                    setattr(
                        self,
                        name,
                        getattr(self.hdfile.root.Data.FiberIndex.cols, name)[:].astype(
                            str
                        ),
                    )
                else:
                    setattr(
                        self,
                        name,
                        getattr(self.hdfile.root.Data.FiberIndex.cols, name)[:],
                    )

        self.coords = SkyCoord(
            self.ra[:] * u.degree, self.dec[:] * u.degree, frame="icrs",
        )
        self.wave_rect = 2.0 * np.arange(1036) + 3470.0

    def query_region(self, coords, radius=3.0 / 3600.0):
        """
        Returns an indexed fiber table for a defined aperture.

        Parameters
        ----------
        self
            Fibers class object
        coords
            astropy coordinate object
        radius
            radius in degrees
        """

        idx = coords.separation(self.coords) < radius * u.degree

        return self.table[idx]

    def query_region_idx(self, coords, radius=3.0):
        """
        Returns an index for a Fibers class object to
        retrieve all fibers in the defined aperture

        self   - Fibers class object
        coords - astropy coordinate object
        radius - astropy quantity object or value in arcsec
        """
        try:
            idx = coords.separation(self.coords) < radius
        except:
            idx = coords.separation(self.coords) < radius * u.arcsec
        return np.where(idx)[0]

    def get_closest_fiber(self, coords, exp=None):
        """
        Returns index to closest fiber to an RA/DEC

        Parameters
        ----------
        coords
            an astropy SkyCoords object
        exp
            exposure number.
        """
        if exp in [1, 2, 3]:
            sel = self.expnum = exp
            fib_idx = coords.match_to_catalog_sky(self.coords[sel])[0]
        else:
            fib_idx = coords.match_to_catalog_sky(self.coords)[0]
        return fib_idx

    def get_image_xy(self, idx, wave_obj):
        """
        Finds the X,Y image position from a row in the
        fibers table

        Note: this returns x,y image coords from 2D
        image arrays produced by Panacea.
        """

        wave_data = self.table[idx]["wavelength"]
        trace_data = self.table[idx]["trace"]

        y = int(np.round(np.interp(wave_obj, wave_data, range(len(wave_data)))))
        x = int(np.round(np.interp(y, range(len(trace_data)), trace_data)))
        return x, y

    def plot_fiber_spectrum(self, idx, type="calfib", xlim=None, ylim=None):
        import matplotlib.pyplot as plt

        if type == "calfib":
            try:
                plt.plot(self.wave_rect, self.table[idx]["calfib"])
                if xlim is not None:
                    plt.xlim(xlim)
                if ylim is not None:
                    plt.ylim(ylim)
                plt.xlabel("wavelength")
                plt.ylabel(type)
            except:
                print("Error plotting calib spectrum")
        else:
            try:
                plt.plot(self.table[idx]["wavelength"], self.table[idx][type])
                if xlim is not None:
                    plt.xlim(xlim)
                if ylim is not None:
                    plt.ylim(ylim)
                plt.xlabel("wavelength")
                plt.ylabel(type)
            except:
                print("Error plotting spectrum")

    def save_fiber_spectrum(self, idx, type="calfib", file="spec.dat"):
        """
        Saves a fiber spectrum at a defined index

        Parameters
        ---------
        self
            Fibers class object
        idx
            index of the fiber in the Fibers Table
        types
            The fiber array type. Options are 'calfib', 'spectrum',
            'sky_spectrum', 'twi_spectrum','error1Dfib', fiber_to_fiber'
        file
            output file name

        Returns
        -------
        An astropy table of the fiber spectra with columns of wavelength,
        spectrum data, error.
        """

        spectab = Table()
        if type == "calfib":
            try:
                print("Saving the flux-calibrated fiber spectrum")
                spectab["wavelength"] = self.wave_rect
                spectab[type] = self.table[idx][type]
                spectab["error"] = self.table[idx]["calfibe"]
            except:
                print("Could not save calibrated fiber spectrum")
        else:
            try:
                spectab["wavelength"] = self.table[idx]["wavelength"]
                spectab[type] = self.table[idx][type]
                spectab["error"] = self.table[idx]["error1Dfib"]
            except:
                print("Could not retrieve Fiber spectrum")
        spectab.write(file, format="ascii", overwrite=True)

    def close(self):
        """ Close the H5 file related to the Fibers call"""
        self.hdfile.close()

    def get_fib_image2D(
        self,
        wave_obj=None,
        fiber_id=None,
        fibnum_obj=None,
        multiframe_obj=None,
        expnum_obj=None,
        width=60,
        height=40,
        imtype="clean_image",
    ):
        """
        Returns an image from the 2D data for a  specific fiber centered at a
        specific wavelength

        Parameters
        ---------
        self
            fibers class object
        wave_obj
            astropy wavelength object
        fiber_id
            fiber_id for a fiber: shotid_exp_multiframe_fibnum
        expnum
            dither exposure number [1,2,3]
        fibnum
            fiber number (This comes from the detections table and relates to
            the fiber index by fibnum = fibidx + 1)
        multiframe
            amp multiframe ID
        imtype
            image option to display
        width
            pixel width to be cutout (image size is 1032 pix)
        height
            pixel height to be cutout (image size is 1032 pix)

        Returns
        -------
        A 2D image of the specified fiber
        """

        if fiber_id:
            idx = np.where(self.fiber_id == fiber_id)[0][0]
            shotid = int(fiber_id[0:11])
            expnum_obj = int(fiber_id[12:13])
            multiframe_obj = fiber_id[14:34]
            fibnum_obj = int(fiber_id[35:39])
        else:
            idx = np.where(
                (self.fibidx == (fibnum_obj - 1))
                * (self.multiframe == multiframe_obj)
                * (self.expnum == expnum_obj)
            )[0][0]
        if np.size(idx) > 1:
            print("Somethings is wrong, found {} fibers".format(np.size(idx)))
            sys.exit()
        elif np.size(idx) == 0:
            print("Could not find a fiber match. Check inputs")
            sys.exit()
        else:
            pass

        x, y = self.get_image_xy(idx, wave_obj)

        im0 = self.hdfile.root.Data.Images.read_where(
            "(multiframe == multiframe_obj) & (expnum == expnum_obj)"
        )

        # create image of forced dims of input width x height

        height = np.minimum(height, 1032)
        width = np.minimum(width, 1032)

        im_base = np.zeros((height, width))

        dx = int(height / 2)
        dy = int(width / 2)

        x1 = np.maximum(0, x - dx)
        x2 = np.minimum(x + dx + (height % 2), 1032)

        y1 = np.maximum(0, y - dy)
        y2 = np.minimum(y + dy + (width % 2), 1032)

        if y1 == 0:
            y1_slice = int(width - (y2 - y1))
            y2_slice = width

        elif y2 == 1032:
            y1_slice = 0
            y2_slice = y2 - y1
        else:
            y1_slice = 0
            y2_slice = width

        if x1 == 0:
            x1_slice = int(height - (x2 - x1))
            x2_slice = height

        elif x2 == 1032:
            x1_slice = 0
            x2_slice = x2 - x1
        else:
            x1_slice = 0
            x2_slice = height

        im_reg = im0[imtype][0][x1:x2, y1:y2]

        im_base[x1_slice:x2_slice, y1_slice:y2_slice] = im_reg

        return im_base

    def return_astropy_table(self):
        """
        Return an astropy table version of the Fibers table
        that can easily be saved

        Returns
        -------
        table : astropy.table:Table
            an astropy table you can save

        """
        table = Table()
        for name in self.hdfile.root.Data.Fibers.colnames:
            table[name] = getattr(self, name)

        return table


def get_fibers_table(
    shot, coords=None, radius=3.0 * u.arcsec, survey="hdr2.1", astropy=True
):
    """
    Returns fiber specta for a given shot.

    Parameters
    ---------
    shot
        either shotid or datevobs
    coords
        astropy coordinate object
    radius
        an astropy quantity object or radius in degrees
    astropy
        flag to make it an astropy table
    survey
        data release you want to access

    Returns
    -------
    A table of fibers within the defined aperture. Will be an astropy table
    object if astropy=True is set

    """

    fileh = open_shot_file(shot, survey=survey.lower())
    fibers = fileh.root.Data.Fibers
    try:
        ra_in = coords.ra.degree
        dec_in = coords.dec.degree
    except:
        print("Coords argument must be an astropy coordinates object")

    try:
        rad_in = radius.to(u.degree)
        rad = radius
    except:
        rad_in = radius / 3600.0
        rad = radius * u.arcsec
        pass

    if survey == "hdr1":
        # search first along ra

        ra_table = fibers.read_where("sqrt((ra - ra_in)**2) < (rad_in + 2./3600)")

        if any(ra_table):
            coords_table = SkyCoord(
                ra_table["ra"] * u.deg, ra_table["dec"] * u.deg, frame="icrs"
            )
            idx = coords.separation(coords_table) < rad
            fibers_table = ra_table[idx]

            fibers_table["calfib"] = fibers_table["calfib"] / 2.0
            fibers_table["calfibe"] = fibers_table["calfibe"] / 2.0

            if astropy:
                fibers_table = Table(fibers_table)

        else:
            fibers_table = None

    else:

        # use FiberIndex table to find fiber_ids
        fiberindex = Fibers(shot, survey=survey)
        fibers_table = fiberindex.query_region(coords, radius=rad_in.value)

        if np.size(fibers_table) > 0:
            if astropy:
                fibers_table = Table(fibers_table)
        else:
            fibers_table = None

    fileh.close()
    return fibers_table


def get_image2D_cutout(
    shot, coords, wave_obj, width=40, height=40, imtype="clean_image", survey="hdr2.1"
):
    """
    Returns an image from the 2D data based on
    ra/dec/wave.

    Parameters
    ----------
    self
        HDF5 file called from open_shot_file('shot')
    coords
        astropy coordinate object
    wave_obj
        astropy wavelength object
    imtype
        image option to display. Options are 'clean_image', 'image', and
        'error'.
    width
        pixel width to be cutout (image size is 1032 pix)
        this is the width X in the wavelength dimension
    height
        pixel height to be cutout (image size is 1032 pix)
        this is height Y in the fiber dimension
    survey
        string identifying data release

    Returns
    -------
    
    """
    fibers = Fibers(shot)

    idx = fibers.get_closest_fiber(coords)
    multiframe_obj = fibers.table.cols.multiframe[idx].astype(str)
    expnum_obj = fibers.table.cols.expnum[idx]
    x, y = fibers.get_image_xy(idx, wave_obj)

    im0 = fibers.hdfile.root.Data.Images.read_where(
        "(multiframe == multiframe_obj) & (expnum == expnum_obj)"
    )

    return im0[imtype][0][
        x - int(np.floor(height / 2)): x + int(np.ceil(height / 2)),
        y - int(np.floor(width / 2)): y + int(np.ceil(width / 2)),
    ]


def get_image2D_amp(
    shot,
    multiframe=None,
    specid=None,
    amp=None,
    ifuslot=None,
    imtype="clean_image",
    expnum=1,
    survey="hdr2.1",
):
    """
    Returns an image from the 2D data based on
    a multiframe or a ifuslot/specid and amp combo

    Parameters
    ---------
    multiframe
        unique amp identifier to display
    specid
        instead of multiframe, you can provide both
        specid and amp or ifuslot and amp
    ifuslot
        provide ifuslot id and amp in place of multiframe
    amp
        inplace of amp you can provide the ifuslot and amp
        or the specid and amp
    imtype
        image option to display
        options are:['spectrum', 'error, clean_image']
    expnum_obj
        integer for the dither/exposure

    Returns
    -------
    a 2D numpy array for the specified amp

    """
    fileh = open_shot_file(shot, survey=survey)

    _expnum = expnum

    if multiframe:
        _multiframe = multiframe
        im0 = fileh.root.Data.Images.read_where(
            "(multiframe == _multiframe) & (expnum == _expnum)"
        )
    elif specid:
        _specid = specid

        if amp:
            _amp = amp
            im0 = fileh.root.Data.Images.read_where(
                "(specid == _specid) & (amp == _amp) & (expnum == _expnum)"
            )
        else:
            print("You must provide both specid and amp")
    elif ifuslot:
        _ifuslot = ifuslot
        if amp:
            _amp = amp
            im0 = fileh.root.Data.Images.read_where(
                "(ifuslot == _ifuslot) & (amp == _amp) & (expnum == _expnum)"
            )
        else:
            print("You must provide both ifuslot and amp")

    else:
        print("You need to provide a multiframe or specid/amp or ifuslot/amp")

    fileh.close()

    return im0[imtype][0]
