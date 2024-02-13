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
from astropy.nddata import bitmask
from astropy.nddata.bitmask import BitFlagNameMap

if not sys.warnoptions:
    warnings.simplefilter("ignore")

try:
    from hetdex_api.config import HDRconfig

    LATEST_HDR_NAME = HDRconfig.LATEST_HDR_NAME
    config = HDRconfig(survey=LATEST_HDR_NAME)

except Exception as e:
    print("Warning! Cannot find or import HDRconfig from hetdex_api!!", e)
    LATEST_HDR_NAME = "hdr4"

    
class CALFIB_DQ(BitFlagNameMap):
    MAIN = 1
    FTF = 2
    CHI2FIB = 4
    BADPIX = 8
    BADAMP = 16
    LARGEGAL = 32
    METEOR = 64
    BADSHOT = 128
    THROUGHPUT = 256
    BADFIB = 512
    SAT = 1024

def open_shot_file(shotid, survey=LATEST_HDR_NAME):
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
    def __init__(
        self,
        shot,
        survey=LATEST_HDR_NAME,
        add_rescor=False,
        add_mask=False,
        mask_version=None,
    ):
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
        add_mask
            boolean option to open mask h5 file handle. Defaults to False for now.
        mask_version
            Force mask version rather than using the default from HDRconfig
        add_rescor
            boolean option to open rescor h5 file. Defaults to False
        survey
            Data release you would like to load, i.e., 'HDR1', 'hdr2', 'hdr2.1'
            This is case insensitive.

        Attributes
        ----------
        hdfile
            h5 file handle for the specified shotid
        maskh5
            h5 file handle for the mask
        rescorh5
            h5 file handle for residual correction file
        table
            h5 table of Fiber data
        coords
            astropy sky coordinates for each fiber
        wave_rect
            rectified wavelength for interpolated spectral data 'calfib',
            'calfibe'
        """

        global config

        self.hdfile = open_shot_file(shot, survey=survey)

        if re.search("v", str(shot)):
            shotid = int(shot.replac("v", ""))
            datevobs = shot
        else:
            shotid = shot
            datevobs = "{}v{}".format(str(shotid)[0:8], str(shotid)[8:11])

        if add_mask:
            try:
                if mask_version is not None:
                    self.mask_version = mask_version
                else:
                    self.mask_version = config.LATEST_MASK_DICT[survey]

                # compose the file handle
                fmask = op.join(
                    config.mask_dir, self.mask_version, "{}.h5".format(datevobs)
                )
                self.maskh5 = tb.open_file(fmask, "r")
            except:
                print(
                    "Could not open mask h5 file for version {}.".format(
                        self.mask_version
                    )
                )
        else:
            self.maskh5 = None
            self.mask_version = None

        if add_rescor:
            filerescor = op.join(config.red_dir, "ffsky_rescor", str(shotid) + ".h5")
            try:
                self.rescorh5 = tb.open_file(filerescor, "r")
            except:
                print(
                    "Could not open {}. Forcing add_rescor to False".format(filerescor)
                )
        else:
            self.rescorh5 = None

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
            self.ra[:] * u.degree,
            self.dec[:] * u.degree,
            frame="icrs",
        )
        self.wave_rect = 2.0 * np.arange(1036) + 3470.0

    def query_region(self, coords, radius=3.5 * u.arcsec):
        """
        Returns an indexed fiber table for a defined aperture.

        Parameters
        ----------
        self
            Fibers class object
        coords
            astropy coordinate object
        radius
            astropy quantity. If no quantity given, assume arcsec
        """
        try:
            idx = coords.separation(self.coords) < radius
        except:
            idx = coords.separation(self.coords) < radius * u.arcsec

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
        """Close the H5 files related to the Fibers call"""
        self.hdfile.close()

        if self.maskh5 is not None:
            self.maskh5.close()

        if self.rescorh5 is not None:
            self.rescorh5.close()

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
            if hasattr(self, name):
                table[name] = getattr(self, name)

        return table


def get_fibers_table(
    shot,
    coords=None,
    ifuslot=None,
    multiframe=None,
    expnum=None,
    radius=3.5 * u.arcsec,
    survey=LATEST_HDR_NAME,
    astropy=True,
    verbose=False,
    rawh5=False,
    F=None,
    fiber_flux_offset=None,
    add_rescor=False,
    add_mask=False,
    mask_options=None,
    mask_version=None,
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
        an astropy quantity object
    astropy
        flag to make it an astropy table. Deprecated. Output is always an astropy table
    survey
        data release you want to access
    rawh5: bool
        if True, this will simply return the fibers from the specified shoth5
        file. If False (the default), any relevent correcctions
        are applied.
    verbose
        print out warnings. Default is False
    F   Fibers class object
        a pre-intiated fibers class object. This is used to limit I/O.
    fiber_flux_offset: 1036 array
        array of values in units of 10**-17 ergs/s/cm2/AA to add
        to each fiber spectrum used in the extraction. Defaults
        to None
    add_rescor
        option to add calfib_ffsky_rescor column generated by Maja Lujan Niemeyer.
        Defaults to False for now
    add_mask
        option to add mask column. Defaults to False for now.
    mask_options
        string array options to select to mask. Default None will select all flags.
        Options are 'MAIN', 'FTF', 'CHI2FIB', 'BADPIX', 'BADAMP', 'LARGEGAL', 'METEOR',
        'BADSHOT', 'THROUGHPUT', 'BADFIB', 'SAT'

    Returns
    -------
    A table of fibers within the defined aperture. Will be an astropy table
    object if astropy=True is set

    """
    if verbose:
        print("Grabbing fibers from {}".format(shot))

    update_F = False  # flag to update Fibers class object if needed

    if F is not None:
        fileh = F.hdfile

        # check if you mask h5 needs to be open

        if add_mask and F.maskh5 is None:
            update_F = True

        if mask_version is None:
            mask_version = F.mask_version
            if verbose:
                print('Using mask_version={}'.format(mask_version))
                                                     
        if F.mask_version != mask_version:
            if verbose:
                print('Updating mask_version to mask_version={}'.format(mask_version))
            update_F = True

        if add_rescor and F.rescorh5 is None:
            update_F = True

        if update_F:

            if verbose:
                print(
                    "Fibers class object is updating for options: add_mask={}, add_rescor={}, mask_version={}".format(
                        add_mask, add_rescor, mask_version
                    )
                )
            
            F = Fibers(
                shot,
                survey=survey.lower(),
                add_rescor=add_rescor,
                add_mask=add_mask,
                mask_version=mask_version,
            )
            
        close_F_at_end = False
    else:
        if verbose:
            print(
                "Opening Fibers class object. If you are running many calls of get_fibers_table consider pre-opening the Fibers class object."
            )
            print(
                    "Fibers class object is opening for options: add_mask={}, add_rescor={}, mask_version={}".format(
                        add_mask, add_rescor, mask_version
                    )
            )

        F = Fibers(
            shot,
            survey=survey.lower(),
            add_rescor=add_rescor,
            add_mask=add_mask,
            mask_version=mask_version,
        )
        fileh = F.hdfile
        close_F_at_end = True

    config = HDRconfig(survey=survey.lower())

    #set index to None for logic later in code
    idx =[]
    
    if coords is not None:
        if verbose:
            print("Gathering fibers in aperture")
        try:
            ra_in = coords.ra.degree
            dec_in = coords.dec.degree
        except:
            print("Coords argument must be an astropy coordinates object")

        idx = F.query_region_idx(coords, radius=radius)

        if len(idx) > 0:
            fibers_table = Table(fileh.root.Data.Fibers.read_coordinates(idx))
        else:
            if verbose:
                print('No fibers found for {} in {}'.format(coords, shot))
            return None
        
    elif multiframe is not None:
        if verbose:
            print("Accessing fibers for {}".format(multiframe))
        multiframe_i = multiframe

        idx = fileh.root.Data.FiberIndex.get_where_list("multiframe == multiframe_i")

        fibers_table = Table(fileh.root.Data.Fibers.read_coordinates(idx))

    elif ifuslot is not None:
        # ensure ifuslot is three digit Unicode string
        ifuslot = ifuslot.zfill(3)
        if isinstance(ifuslot, bytes):
            ifuslot = ifuslot.decode()

        if verbose:
            print("Acessing fibers for ifuslot {}".format(ifuslot))

        multiframe_array = np.unique(F.multiframe)

        ifuslot_array = np.array([x[10:13] for x in multiframe_array])

        idx = []
        for mf in multiframe_array[ifuslot_array == ifuslot]:
            idx.extend(fileh.root.Data.FiberIndex.get_where_list("multiframe == mf"))

        fibers_table = Table(fileh.root.Data.Fibers.read_coordinates(idx))
        
    else:
        if verbose:
            print("Loading full fibers table for shot. This will take some time.")
            print(
                "Consider requesting a single multiframe and expnum or querying by coordinates"
            )
        fibers_table = Table(fileh.root.Data.Fibers.read())

    if add_rescor:
        if verbose:
            print("Adding calfib_ffsky_rescor column")
        ffsky_rescor = F.rescorh5.root.calfib_ffsky_rescor.read()

        if len(idx) > 0:
            fibers_table["calfib_ffsky_rescor"] = ffsky_rescor[idx]
        else:
            fibers_table["calfib_ffsky_rescor"] = ffsky_rescor

    if add_mask:
        if verbose:
            print('Adding mask array using mask_version={}'.format(F.mask_version)) 

        if len(idx) > 0:
            bitmaskDQ = F.maskh5.root.CalfibDQ.read_coordinates( idx, 'calfib_dq')
        else:
            bitmaskDQ = F.maskh5.root.CalfibDQ.read()['calfib_dq']
            
        if mask_options is not None:
            if mask_option=='bitmask':
                #return the mask the full bitmask array
                fiber_table['mask'] = bitmaskDQ
                
            # get mask name dictionary
            mask_names = []
            for i in CALFIB_DQ.__dict__.keys():
                if '_' not in i:
                    mask_names.append(i)

            for mask_option in mask_options:
                if mask_option not in mask_names:
                    print('mask_option={} not in mask_option dictionary for mask_version={}. Mask options are:{}'.format(mask_option, mask_version, mask_names))

        else:
            if verbose:
                print('Masking everything. mask_options set to None')
            bool_mask = bitmask.bitfield_to_boolean_mask(bitmaskDQ, good_mask_value=True)
            fibers_table['mask'] = bool_mask
            
    # down-select on expnum if desired
    if expnum is not None:
        if verbose:
            print("Accessing fibers for expnum {}".format(expnum))
        fibers_table = fibers_table[fibers_table["expnum"] == expnum]
        
    if rawh5:
        intensityunit = 10**-17 * u.erg / (u.cm**2 * u.s * 2 * u.AA)
    else:
        intensityunit = 10**-17 * u.erg / (u.cm**2 * u.s * u.AA)

        if verbose:
            print("Convert to spectral density units in 10^-17 ergs/s/cm^2/AA")

        fibers_table["calfib"] /= 2.0
        fibers_table["calfibe"] /= 2.0

        if survey.lower() == "hdr2.1":
            fibers_table["spec_fullsky_sub"] /= 2.0
        else:
            fibers_table["calfib_ffsky"] /= 2.0

        if survey.lower() == "hdr3":
            if verbose:
                print("Applying spectral correction from WD modelling")
            wd_corr = Table.read(
                config.wdcor, format="ascii.no_header", names=["wave", "corr"]
            )
            fibers_table["calfib"] /= wd_corr["corr"]
            fibers_table["calfib_ffsky"] /= wd_corr["corr"]
            fibers_table["calfibe"] /= wd_corr["corr"]

            if verbose:
                print("Adjusting noise values by 7% where applicable")
            # adjust noise at IFU edges by factor of 1.07
            sel_fib1 = (
                (fibers_table["amp"] == b"RU") | (fibers_table["amp"] == b"LL")
            ) & (fibers_table["fibnum"] <= 12)
            sel_fib2 = (
                (fibers_table["amp"] == b"LU") | (fibers_table["amp"] == b"RL")
            ) & (fibers_table["fibnum"] >= 101)

            sel_fib = sel_fib1 | sel_fib2

            fibers_table["calfibe"][sel_fib] *= 1.07

            if fiber_flux_offset is not None:
                if verbose:
                    print(
                        "Applying supplied fiber_flux_offset: {}".format(
                            fiber_flux_offset
                        )
                    )

                # DD 2023-09-19 replace calfibe 0 spot with nan wo the fiber_flux_offset does not create "data"
                sel_nan = fibers_table["calfibe"] == 0
                fibers_table["calfib"][sel_nan] = np.nan
                fibers_table["calfib_ffsky"][sel_nan] = np.nan

                fibers_table["calfib"] += fiber_flux_offset
                fibers_table["calfib_ffsky"] += fiber_flux_offset
                fibers_table["calfib"][sel_nan] = 0
                fibers_table["calfib_ffsky"][sel_nan] = 0

                if add_rescor:
                    fibers_table["calfib_ffsky_rescor"][sel_nan] = 0

        elif float(survey.lower()[3:]) >= 4.0:  # only HDR4 and up
            if verbose:
                print("Applying spectral correction from WD modelling")
            wd_corr = Table.read(
                config.wdcor, format="ascii.no_header", names=["wave", "corr"]
            )
            fibers_table["calfib"] /= wd_corr["corr"]
            fibers_table["calfib_ffsky"] /= wd_corr["corr"]
            fibers_table["calfibe"] /= wd_corr["corr"]

            if shot < 20210901000:
                early_2019_hdr4 = np.loadtxt(
                    op.join(config.bad_dir, "hdr4_2019.shots"), dtype=int
                )
                if shot not in early_2019_hdr4:
                    # for HDR3 frames adjust noise at IFU by factor of 1.07
                    if verbose:
                        print("Adjusting noise values by 7% where applicable")

                    sel_fib1 = (
                        (fibers_table["amp"] == b"RU") | (fibers_table["amp"] == b"LL")
                    ) & (fibers_table["fibnum"] <= 12)
                    sel_fib2 = (
                        (fibers_table["amp"] == b"LU") | (fibers_table["amp"] == b"RL")
                    ) & (fibers_table["fibnum"] >= 101)

                    sel_fib = sel_fib1 | sel_fib2

                    fibers_table["calfibe"][sel_fib] *= 1.07

            if fiber_flux_offset is not None:
                if verbose:
                    print(
                        "Applying supplied fiber_flux_offset: {}".format(
                            fiber_flux_offset
                        )
                    )

                # DD 2023-09-19 replace calfibe 0 spot with nan wo the fiber_flux_offset does not create "data"
                sel_nan = fibers_table["calfibe"] == 0
                fibers_table["calfib"][sel_nan] = np.nan
                fibers_table["calfib_ffsky"][sel_nan] = np.nan

                fibers_table["calfib"] += fiber_flux_offset
                fibers_table["calfib_ffsky"] += fiber_flux_offset
                fibers_table["calfib"][sel_nan] = 0
                fibers_table["calfib_ffsky"][sel_nan] = 0

    if np.size(fibers_table) > 0:
        if astropy:
            # convert to an astropy table format and add units
            fibers_table = Table(fibers_table)
            # add units
            fibers_table["calfib"].unit = intensityunit
            if survey.lower() == "hdr2.1":
                fibers_table["spec_fullsky_sub"].unit = intensityunit
            else:
                fibers_table["calfib_ffsky"].unit = intensityunit
            if add_rescor:
                fibers_table["calfib_ffsky_rescor"].unit = intensityunit
            fibers_table["calfibe"].unit = intensityunit
            fibers_table["ra"].unit = u.deg
            fibers_table["dec"].unit = u.deg
            fibers_table["wavelength"].unit = u.AA
    else:
        fibers_table = None

    # close h5 file handle if not passed initially
    if close_F_at_end:
        F.close()

    return fibers_table


def get_image2D_cutout(
    shot,
    coords,
    wave_obj,
    width=40,
    height=40,
    imtype="clean_image",
    survey=LATEST_HDR_NAME,
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
        x - int(np.floor(height / 2)) : x + int(np.ceil(height / 2)),
        y - int(np.floor(width / 2)) : y + int(np.ceil(width / 2)),
    ]


def get_image2D_amp(
    shot,
    multiframe=None,
    specid=None,
    amp=None,
    ifuslot=None,
    imtype="clean_image",
    expnum=1,
    survey=LATEST_HDR_NAME,
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
