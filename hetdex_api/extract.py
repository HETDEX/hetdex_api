# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 11:48:55 2019

Note: this uses nway and pandas

@authors: gregz, Erin Mentuch Cooper, Daniel Farrow, Dustin Davis
"""

import numpy as np
from numba import jit

from astropy.table import Table
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.coordinates import SkyCoord
from astropy.modeling.models import Moffat2D, Gaussian2D
from astropy import units as u
from scipy.interpolate import griddata, LinearNDInterpolator

# Need to do this to set backend to avoid crashes when
# importing nway
import matplotlib
matplotlib.use("agg")

NWAY_IMPORTED = False

from hetdex_api.shot import Fibers, open_shot_file, get_fibers_table
from hetdex_api.input_utils import setup_logging

from hetdex_api.extinction import get_2pt1_extinction_fix

try:
    from hetdex_api.config import HDRconfig

    LATEST_HDR_NAME = HDRconfig.LATEST_HDR_NAME
    config = HDRconfig()

except Exception as e:
    print("Warning! Cannot find or import HDRconfig from hetdex_api!!", e)
    LATEST_HDR_NAME = "hdr2.1"
    config = None

@jit(nopython=True)
def sclean(waves, fiber_data, fiber_error, mask,
           npix = 5, error_scale = 4.0):
    """
    A version of `sclean` from 
    `karlspipe/fitrsdecsp.f`. Bad data
    is replaced with the average of a +/-
    `npix` window around the nearest element.
    The error is replaced with `error_scale`
    times the error of the nearest element.
    These value where tuned by K. Gebhardt to
    optimise the detections and minimise false
    detections.
    
    Parameters 
    ----------
    fiber_data, fiber_error : 1D array
        the data and error arrays of
        the fiber
    mask : 1D array
        the mask for this fiber (so
        the bad data can be unmasked
        after replacement)
    npix : int
        the size of the +/- pixel
        range to average to find the 
        replacement value
    error_scale : float
        scale the error by this
        when replacing

    Returns
    -------
    fiber_data_out, fiber_error_out, mask : 1D array
        the input data with bad values 
        replaced and unmasked
    """
    fiber_data_out = np.copy(fiber_data)
    fiber_error_out = np.copy(fiber_error)

    bad = (np.abs(fiber_data) < 1e-40) | (np.abs(fiber_error) < 1e-40) | np.logical_not(mask)
    good = np.logical_not(bad)
    good_indices = np.arange(len(waves))[good]
 
    # If it's all bad, give up!
    if np.all(bad):
        return fiber_data_out, fiber_error_out, mask, bad

    mask[bad] = True

    for idx in np.arange(len(waves))[bad]:
        jnear = good_indices[np.argmin(np.abs(waves[good] - waves[idx]))]
        arr = fiber_data[np.max(np.array([0, jnear - npix])): jnear + npix + 1]
        fiber_data_out[idx] = np.mean(arr[np.abs(arr) > 0])
        #print(np.mean(arr[arr > 0]), idx, any(np.isnan(arr)), len(arr[arr > 0]))
        fiber_error_out[idx] = error_scale*fiber_error[jnear]        

    return fiber_data_out, fiber_error_out, mask, bad


class Extract:
    def __init__(self, wave=None, apply_update=True):
        """
        Initialize Extract class

        Parameters
        ----------
        wave: numpy 1d array
            wavelength of calfib extension for hdf5 files, does not need to be
            set unless needed by development team
        apply_update: bool
            flag to turn on/off automatic API updates. For example in HDR3
            this includes adjustment to noise arrays by 1.07 factor 
            to fibnum=1-12 for LL/RU channel and fibnum=101-112 for LU/RL.
            It also includes spectral calibration adjustment wdcor.txt
            Default is to include it
        
        """
        if wave is not None:
            self.wave = wave
        else:
            self.wave = self.get_wave()
        self.get_ADR()
        self.log = setup_logging("Extract")
        self.fibers = None
        self.set_dither_pattern()
        self.apply_update = apply_update

    def set_dither_pattern(self, dither_pattern=None):
        """ 
        Set dither pattern (default if None given)
        
        Parameters
        ----------
        dither_pattern: numpy array (length of exposures)
            only necessary if the dither pattern isn't the default
        """
        if dither_pattern is None:
            self.dither_pattern = np.array([[0.0, 0.0], [1.27, -0.73], [1.27, 0.73]])
        else:
            self.dither_pattern = dither_pattern

    def get_wave(self):
        """ Return implicit wavelength solution for calfib extension """
        return np.linspace(3470, 5540, 1036)

    def get_ADR(self, angle=0.0):
        """ 
        Use default ADR from Karl Gebhardt (adjusted by Greg Zeimann)

        Parameters
        ----------
        angle: float
            angle=0 along x-direction, appropriate if observing at PA
        dither_pattern: numpy array (length of exposures)
            only necessary if the dither pattern isn't the default
        """
        wADR = [3500.0, 4000.0, 4500.0, 5000.0, 5500.0]
        ADR = [-0.74, -0.4, -0.08, 0.08, 0.20]
        ADR = np.polyval(np.polyfit(wADR, ADR, 3), self.wave)
        self.ADRx = np.cos(np.deg2rad(angle)) * ADR
        self.ADRy = np.sin(np.deg2rad(angle)) * ADR

    def load_shot(
            self, shot_input, survey=LATEST_HDR_NAME, dither_pattern=None, fibers=True, add_mask=True, args=None,
    ):
        """
        Load fiber info from hdf5 for given shot_input
        
        Parameters
        ----------
        shot_input: str
            e.g., 20190208v024 or 20190208024
        """
        self.shot = shot_input
        self.survey = survey

        if fibers:
            self.fibers = Fibers(self.shot, survey=survey, add_mask=add_mask,args=args)
            self.shoth5 = self.fibers.hdfile
        else:
            self.fibers = None
            self.shoth5 = None

        self.set_dither_pattern(dither_pattern=dither_pattern)

    def convert_radec_to_ifux_ifuy(self, ifux, ifuy, ra, dec, rac, decc):
        """
        Input SkyCoord object for sources to extract
        
        Parameters
        ----------
        coord: SkyCoord object
            single ra and dec to turn into ifux and ifuy coordinates
        """
        if len(ifuy) < 2:
            return None, None
        ifu_vect = np.array([ifuy[1] - ifuy[0], ifux[1] - ifux[0]])
        radec_vect = np.array(
            [(ra[1] - ra[0]) * np.cos(np.deg2rad(dec[0])), dec[1] - dec[0]]
        )
        V = np.sqrt(ifu_vect[0] ** 2 + ifu_vect[1] ** 2)
        W = np.sqrt(radec_vect[0] ** 2 + radec_vect[1] ** 2)
        scale_vect = np.array([3600.0, 3600.0])
        v = ifu_vect * np.array([1.0, 1.0])
        w = radec_vect * scale_vect
        W = np.sqrt(np.sum(w ** 2))
        ang1 = np.arctan2(v[1] / V, v[0] / V)
        ang2 = np.arctan2(w[1] / W, w[0] / W)
        ang1 += (ang1 < 0.0) * 2.0 * np.pi
        ang2 += (ang2 < 0.0) * 2.0 * np.pi
        theta = ang1 - ang2
        dra = (rac - ra[0]) * np.cos(np.deg2rad(dec[0])) * 3600.0
        ddec = (decc - dec[0]) * 3600.0
        dx = np.cos(theta) * dra - np.sin(theta) * ddec
        dy = np.sin(theta) * dra + np.cos(theta) * ddec
        yc = dx + ifuy[0]
        xc = dy + ifux[0]
        return xc, yc

    def get_fiberinfo_for_coords(self,
                                 coords,
                                 radius=3.5,
                                 ffsky=False,
                                 return_fiber_info=False,
                                 fiber_lower_limit=3,
                                 verbose=False,
                                 nmax=5000,
                                 fiber_flux_offset=None,
    ):
        """ 
        Grab fibers within a radius and get relevant info,
        optimised for searching for a longer list of 
        coordinates. Uses Nway (Salvato et al 2018)
        
        Parameters
        ----------
        coord: SkyCoord Object
            a SkyCoord object with multiple ra and dec
            positions
        radius:
            radius to extract fibers in arcsec
        ffsky: bool
            Flag to choose local (ffsky=False) or full frame (ffsky=True)
            sky subtraction
        return_fiber_info: bool
            Return full fibers table. This is needed to get additional
            masking and to debug fiberid weights
        fiber_lower_limit : int
            Minimum number of fibers needed in aperture to
            return a result
        nmax : int
            Maximum number of coordinates to consider at once
            during the match. If you put in more it loops
            over them.
        fiber_flux_offset: 1036 array
            array of values in units of 10**-17 ergs/s/cm2/AA to add
            to each fiber spectrum used in the extraction. Defaults
            to None 

        Returns
        -------
        icoord_all : array (length of number of fibers)
            the indices of the input coordinates where
            there are sufficient number of fibers to
            produce an output. To get all the fibers
            matched to a given input coordinate, select
            on this column. 
        seps_all : array (length of number of fibers)
            separation of each fiber from the input
            coordinates
        ifux: numpy array (length of number of fibers)
            ifu x-coordinate accounting for dither_pattern
        ifuy: numpy array (length of number of fibers)
            ifu y-coordinate accounting for dither_pattern
        ra: numpy array (length of number of fibers)
            Right ascension of fibers
        dec: numpy array (length of number of fibers)
            Declination of fibers
        spec: numpy 2d array (number of fibers by wavelength dimension)
            Calibrated spectra for each fiber
        spece: numpy 2d array (number of fibers by wavelength dimension)
            Error for calibrated spectra
        mask: numpy 2d array (number of fibers by wavelength dimension)
            Mask of good values for each fiber and wavelength
        fiberid: numpy array (length of number of fibers)
            array of fiberids for fibers used in the extraction.
            Returned only if return_fiber_info is True
        multiframe_array: numpy array (length of number of fibers)
            array of amp IDs/multiframe values for each fiber.
            Return only if return_fiber_info is True

        """

        global NWAY_IMPORTED

        if not self.fibers:
            raise Exception("Only supported with preloaded Fibers class")
            
        if fiber_lower_limit < 2:
            raise Exception("fiber_lower_limit must be greater than 2")

        import nwaylib
        from nwaylib import _create_match_table
        from nwaylib.logger import NormalLogger, NullOutputLogger
        NWAY_IMPORTED = True

        # remove NaN fibers
        notnan = np.isfinite(self.fibers.coords.ra.value) & np.isfinite(self.fibers.coords.dec.value)

        # Save original indices to use later
        indices_original = np.arange(len(self.fibers.coords.ra.value))[notnan]
        fibers_cat = {"name" : "fibers", "ra" : self.fibers.coords.ra.value[notnan], 
                      "dec" : self.fibers.coords.dec.value[notnan], 
                      "error" : np.ones(sum(notnan))}

        if len(coords) > nmax:
            if verbose:
                print("Number of coords exceeds nmax, splitting into sets.")
            nsets = int(np.ceil(len(coords)/nmax))
        else:
            nsets = 1

        for i in range(nsets):

            # Match the two tables with Nway
            tcoords =  coords[i*nmax:(i+1)*nmax]
            cat1 = {"name" : "sources", "ra" : tcoords.ra.value, "dec" : tcoords.dec.value, 
                    "error" : np.ones(len(tcoords))}

            if verbose:
                ttable, resultstable, separations, errors = _create_match_table([cat1, fibers_cat], radius, 
                                                                               NormalLogger())
            else:
                ttable, resultstable, separations, errors = _create_match_table([cat1, fibers_cat], radius, 
                                                                               NullOutputLogger())        

            # Shift index along 
            ttable["sources"] = ttable["sources"] + i*nmax

            # remove sources with no matched fibers
            ttable = ttable[ttable["fibers"] > -1]

            if i == 0:
                table = ttable
            else:
                table = table.append(ttable)

        # Indices of matched positions and fibers
        icoord_all = table["sources"].to_numpy()
        idx_all = table["fibers"].to_numpy()
        seps_all = table["Separation_sources_fibers"].to_numpy()
           
        # Remember to grab the fibers using their original table indices
        table_here = self.fibers.table.read_coordinates(indices_original[idx_all])
        ifux = table_here["ifux"]
        ifuy = table_here["ifuy"]
        ra = table_here["ra"]
        dec = table_here["dec"]
           
        if ffsky:
            if self.survey == 'hdr2.1':
                spec = table_here["spec_fullsky_sub"] / 2.0
            else:
                spec = table_here["calfib_ffsky"] / 2.0
        else:
            spec = table_here["calfib"] / 2.0

        spece = table_here["calfibe"] / 2.0

        if self.survey == 'hdr2.1':
            if verbose:
                print('Removing extinction of E(B-V)=0.02')
            #apply HDR2.1 E(B-V)=0.02 extinction fix to fiber arrays
            fix = get_2pt1_extinction_fix()
            spec /= fix(self.wave)
            spece /= fix(self.wave)
            
        elif (self.survey == 'hdr3') & self.apply_update:
            if verbose:
                print('Increasing noise by 7% for fib 1-12 on LL/RU and fib 101-112 on LU/RL')

            sel_fib1 = ((table_here['amp'] == b'RU') | (table_here['amp'] == b'LL')) & (table_here['fibnum'] <= 12)
            sel_fib2 = ((table_here['amp'] == b'LU') | (table_here['amp'] == b'RL')) & (table_here['fibnum'] >= 101)
            sel_fib = sel_fib1 | sel_fib2

            spece[sel_fib] *= 1.07
            
            # apply WD correction
            if verbose:
                print("Applying spectral correction from WD modelling")
            wd_corr = Table.read(
                config.wdcor, format="ascii.no_header", names=["wave", "corr"]
            )
            spec /= wd_corr['corr']
            spece /= wd_corr['corr']

            if fiber_flux_offset is not None:
                if verbose:
                    print("Applying supplied fiber_flux_offset: {}".format(fiber_flux_offset))

                # DD 2023-09-19 replace calfibe 0 spot with nan wo the fiber_flux_offset does not create "data"
                sel_nan = spece == 0
                spec[sel_nan] = np.nan
                spec += fiber_flux_offset
                spec[sel_nan] = 0
        elif (float(self.survey[3:]) >= 4.0) and self.apply_update: #HDR4 and up
            #apply the WD correction but not the earlier HDR3 1.07x noise correction to 1st and last 12 fibers

            if verbose:
                print("Applying spectral correction from WD modelling")
            wd_corr = Table.read(
                config.wdcor, format="ascii.no_header", names=["wave", "corr"]
            )
            spec /= wd_corr['corr']
            spece /= wd_corr['corr']

            if fiber_flux_offset is not None:
                if verbose:
                    print("Applying supplied fiber_flux_offset: {}".format(fiber_flux_offset))

                # DD 2023-09-19 replace calfibe 0 spot with nan wo the fiber_flux_offset does not create "data"
                sel_nan = spece == 0
                spec[sel_nan] = np.nan
                spec += fiber_flux_offset
                spec[sel_nan] = 0
                
        ftf = table_here["fiber_to_fiber"]

        if self.survey == "hdr1":
            mask = table_here["Amp2Amp"]
            mask = (mask > 1e-8) * (np.median(ftf, axis=1) > 0.5)[:, np.newaxis]
        else:
            mask = self.fibers.table.read_coordinates(indices_original[idx_all], "calfibe")
            mask = (mask > 1e-8) * (np.median(ftf, axis=1) > 0.5)[:, np.newaxis]

        expn = np.array(
               table_here["expnum"], dtype=int
               )
        mf_array = table_here["multiframe"].astype(str)
        fiber_id_array = table_here["fiber_id"].astype(str)
                
        ifux[:] = ifux + self.dither_pattern[expn - 1, 0]
        ifuy[:] = ifuy + self.dither_pattern[expn - 1, 1]
        
        xc = []
        yc = []
        for ic in np.unique(icoord_all):
            sel = (icoord_all == ic)
            n = len(ifux[sel])
            if n >= fiber_lower_limit:
                txc, tyc = self.convert_radec_to_ifux_ifuy(
                                   ifux[sel], ifuy[sel], ra[sel], dec[sel], 
                                   coords.ra.deg[ic], coords.dec.deg[ic]
                                   )

            else:
                # if not enough fibers set position to 999
                # and flag by setting index to -1
                txc = np.array(999.)
                tyc = np.array(999.)
                icoord_all[sel] = -1
                       
            xc.extend(txc.repeat(n))
            yc.extend(tyc.repeat(n))
                             
        xc = np.array(xc)
        yc = np.array(yc)
                    
        if return_fiber_info:
            return icoord_all, seps_all, ifux, ifuy, xc, yc, ra, dec, spec, spece, mask, fiber_id_array, mf_array
        else:
            return icoord_all, seps_all, ifux, ifuy, xc, yc, ra, dec, spec, spece, mask
    
    
    def get_fiberinfo_for_coord(self,
                                coord,
                                radius=3.5,
                                ffsky=False,
                                ffsky_rescor=False,
                                return_fiber_info=False,
                                fiber_lower_limit=3,
                                verbose=False,
                                fiber_flux_offset=None,
                                add_mask=False,
                                mask_options=None,
    ):
        """ 
        Grab fibers within a radius and get relevant info
    
        Parameters
        ----------
        coord: SkyCoord Object
            a single SkyCoord object for a given ra and dec
        radius:
            radius to extract fibers in arcsec
        ffsky: bool
            Flag to choose local (ffsky=False) or full frame (ffsky=True)
            sky subtraction
        return_fiber_info: bool
            Return full fibers table. This is needed to get additional
            masking and to debug fiberid weights
        fiber_lower_limit : int
            Minimum number of fibers needed in aperture to
            return a result
        fiber_flux_offset: 1036 array
            array of values in units of 10**-17 ergs/s/cm2/AA to add
            to each fiber spectrum used in the extraction. Defaults
            to None 
        ffsky_rescor: bool
            Flag to use updated spectra for ffsky with residual correction applied.
            Corrected fiber arrays generated by Maja Lujan Niemeyer 2023-12-14.
            This is a post-HDR4 correction
        add_mask: bool
            Option to use updated mask model created by EMC 2023-12-15. This is a
            post-HDR4 correction. Will use default model unless mask_options is provided
        mask_options
            string or array of strings as options to select to mask. Default None 
            will select all flags. Set this to 'BITMASK' to return the full bitmask array.
            Options are 'MAIN', 'FTF', 'CHI2FIB', 'BADPIX', 'BADAMP', 'LARGEGAL', 'METEOR',
            'BADSHOT', 'THROUGHPUT', 'BADFIB', 'SAT'
            If BITMASK appears as any element in the list, it overrides all others 
            and returns the full bitmask array.

        Returns
        -------
        ifux: numpy array (length of number of fibers)
            ifu x-coordinate accounting for dither_pattern
        ifuy: numpy array (length of number of fibers)
            ifu y-coordinate accounting for dither_pattern
        ra: numpy array (length of number of fibers)
            Right ascension of fibers
        dec: numpy array (length of number of fibers)
            Declination of fibers
        spec: numpy 2d array (number of fibers by wavelength dimension)
            Calibrated spectra for each fiber
        spece: numpy 2d array (number of fibers by wavelength dimension)
            Error for calibrated spectra
        mask: numpy 2d array (number of fibers by wavelength dimension)
            Mask of good values for each fiber and wavelength
        fiberid: numpy array (length of number of fibers)
            array of fiberids for fibers used in the extraction.
            Returned only if return_fiber_info is True
        multiframe_array: numpy array (length of number of fibers
            array of amp IDs/multiframe values for each fiber.
            Return only if return_fiber_info is True
        """

        if fiber_lower_limit < 2:
            raise Exception("fiber_lower_limit must be greater than 2")
        
        if self.apply_update is False:
            if verbose:
                print('No calfib update applied')
            idx = self.fibers.query_region_idx(coord, radius=radius)

            if len(idx) < fiber_lower_limit:
                return None

            ifux = self.fibers.table.read_coordinates(idx, "ifux")
            ifuy = self.fibers.table.read_coordinates(idx, "ifuy")
            ra = self.fibers.table.read_coordinates(idx, "ra")
            dec = self.fibers.table.read_coordinates(idx, "dec")

            if ffsky:
                if self.survey == 'hdr2.1':
                    spec = self.fibers.table.read_coordinates(idx, "spec_fullsky_sub") / 2.0
                else:
                    spec = self.fibers.table.read_coordinates(idx, "calfib_ffsky") / 2.0
            else:
                spec = self.fibers.table.read_coordinates(idx, "calfib") / 2.0

            spece = self.fibers.table.read_coordinates(idx, "calfibe") / 2.0

            if self.survey.lower() in ['pdr1','pdr2','pdr3']:
                pass
            else:
                ftf = self.fibers.table.read_coordinates(idx, "fiber_to_fiber")

            if self.survey == "hdr1":
                mask = self.fibers.table.read_coordinates(idx, "Amp2Amp")
                mask = (mask > 1e-8) * (np.median(ftf, axis=1) > 0.5)[:, np.newaxis]
            else:
                mask = self.fibers.table.read_coordinates(idx, "calfibe")
                if self.survey.lower() in ['pdr1','pdr2','pdr3']:
                    mask = (mask > 1e-8) * (spec != 0)
                else:
                    mask = (mask > 1e-8) * (np.median(ftf, axis=1) > 0.5)[:, np.newaxis] * (spec != 0)

            expn = np.array(
                self.fibers.table.read_coordinates(idx, "expnum"), dtype=int
            )
            mf_array = self.fibers.table.read_coordinates(
                idx, "multiframe").astype(str)
            fiber_id_array = self.fibers.table.read_coordinates(
                idx, "fiber_id").astype(str)
        else:
            if verbose:
                print('Calfib updates applied')

            if self.fibers is None:
                self.fibers = Fibers(self.shot,
                                     survey=self.survey,
                                     add_rescor=ffsky_rescor,
                                     add_mask=add_mask)

            fib_table = get_fibers_table(
                self.shot, coord,
                survey=self.survey,
                radius=radius,
                verbose=False,
                astropy=True,
                F=self.fibers,
                fiber_flux_offset=fiber_flux_offset,
                add_rescor=ffsky_rescor,
                add_mask=add_mask,
                mask_options=mask_options,
            )

            if np.size(fib_table) < fiber_lower_limit:
                return None

            ifux = fib_table["ifux"]
            ifuy = fib_table["ifuy"]
            ra = fib_table["ra"]
            dec = fib_table["dec"]

            if ffsky_rescor:
                try:
                    spec = fib_table['calfib_ffsky_rescor']
                except:
                    spec = fib_table['calfib_ffsky']
                    
            elif ffsky:
                if self.survey == 'hdr2.1':
                    spec = fib_table["spec_fullsky_sub"]
                else:
                    spec = fib_table["calfib_ffsky"]
            else:
                spec = fib_table["calfib"]
            spece = fib_table["calfibe"]
            
            if self.survey.lower() in ['pdr1','pdr2','pdr3']:
                pass
            else:
                ftf = fib_table["fiber_to_fiber"]
            
            if self.survey == "hdr1":
                mask = fib_table["Amp2Amp"]
                mask = (mask > 1e-8) * (np.median(ftf, axis=1) > 0.5)[:, np.newaxis]
            elif add_mask:
                mask = fib_table['mask']
            else:
                mask = fib_table["calfibe"]
                if self.survey.lower() in ['pdr1','pdr2','pdr3']:
                    mask = (mask > 1e-8) * (spec != 0)
                else:
                    mask = (mask > 1e-8) * (np.median(ftf, axis=1) > 0.5)[:, np.newaxis] * (spec != 0)
                
            expn = np.array(fib_table["expnum"], dtype=int)
            mf_array = fib_table['multiframe'].astype(str)
            try:
                fiber_id_array = fib_table['fiber_id'].astype(str)
            except:
                fiber_id_array = []

        if self.survey == 'hdr2.1':
            if verbose:
                print('Removing extinction of E(B-V)=0.02')
            #apply HDR2.1 E(B-V)=0.02 extinction fix to fiber arrays
            print('applying extinction fix')
            fix = get_2pt1_extinction_fix()
            spec /= fix(self.wave)
            spece /= fix(self.wave)

        ifux[:] = ifux + self.dither_pattern[expn - 1, 0]
        ifuy[:] = ifuy + self.dither_pattern[expn - 1, 1]

        xc, yc = self.convert_radec_to_ifux_ifuy(
            ifux, ifuy, ra, dec, coord.ra.deg, coord.dec.deg
        )
        if return_fiber_info:
            return ifux, ifuy, xc, yc, ra, dec, spec, spece, mask, fiber_id_array, mf_array
        else:
            return ifux, ifuy, xc, yc, ra, dec, spec, spece, mask

    def get_starcatalog_params(self):
        """
        Load Star Catalog coordinates, g' magnitude, and star ID
        
        Returns
        -------
        coords: SkyCoord Object
            SkyCoord object of the ra and dec's of the stars
        gmag: numpy array (float)
            g' magnitude of the stars in the star catalog
        starid: numpy array (int)
            Object ID from the original catalog of the stars (e.g., SDSS)
        """
        if not hasattr(self, "shoth5"):
            pass
        else:
            self.shoth5 = open_shot_file(self.shot, survey=self.survey)
        
            #self.log.warning("Please do load_shot to get star catalog params.")
            #return None

        ras = self.shoth5.root.Astrometry.StarCatalog.cols.ra_cat[:]
        decs = self.shoth5.root.Astrometry.StarCatalog.cols.dec_cat[:]
        gmag = self.shoth5.root.Astrometry.StarCatalog.cols.g[:]
        starid = self.shoth5.root.Astrometry.StarCatalog.cols.star_ID[:]

        coords = SkyCoord(ras * u.deg, decs * u.deg, frame="fk5")
        return coords, gmag, starid

    def intersection_area(self, d, R, r):
        """
        Return the area of intersection of two circles.
        The circles have radii R and r, and their centers are separated by d.
    
        Parameters
        ----------
        d : float
            separation of the two circles
        R : float
            Radius of the first circle
        r : float
            Radius of the second circle
        """

        if d <= abs(R - r):
            # One circle is entirely enclosed in the other.
            return np.pi * min(R, r) ** 2
        if d >= r + R:
            # The circles don't overlap at all.
            return 0

        r2, R2, d2 = r ** 2, R ** 2, d ** 2
        alpha = np.arccos((d2 + r2 - R2) / (2 * d * r))
        beta = np.arccos((d2 + R2 - r2) / (2 * d * R))
        answer = (
            r2 * alpha
            + R2 * beta
            - 0.5 * (r2 * np.sin(2 * alpha) + R2 * np.sin(2 * beta))
        )
        return answer

    def tophat_psf(self, radius, boxsize, scale, fibradius=0.75):
        """
        Tophat PSF profile image 
        (taking fiber overlap with fixed radius into account)
        
        Parameters
        ----------
        radius: float
            Radius of the tophat PSF
        boxsize: float
            Size of image on a side for Moffat profile
        scale: float
            Pixel scale for image
        
        Returns
        -------
        zarray: numpy 3d array
            An array with length 3 for the first axis: PSF image, xgrid, ygrid
        """
        xl, xh = (0.0 - boxsize / 2.0, 0.0 + boxsize / 2.0 + scale)
        yl, yh = (0.0 - boxsize / 2.0, 0.0 + boxsize / 2.0 + scale)
        x, y = (np.arange(xl, xh, scale), np.arange(yl, yh, scale))
        xgrid, ygrid = np.meshgrid(x, y)
        t = np.linspace(0.0, np.pi * 2.0, 360)
        r = np.arange(
            radius - fibradius, radius + fibradius * 7.0 / 6.0, 1.0 / 6.0 * fibradius
        )
        xr, yr, zr = ([], [], [])
        for i, ri in enumerate(r):
            xr.append(np.cos(t) * ri)
            yr.append(np.sin(t) * ri)
            zr.append(
                self.intersection_area(ri, radius, fibradius)
                * np.ones(t.shape)
                / (np.pi * fibradius ** 2)
            )
        xr, yr, zr = [np.hstack(var) for var in [xr, yr, zr]]
        psf = griddata(
            np.array([xr, yr]).swapaxes(0, 1), zr, (xgrid, ygrid), method="linear"
        )
        psf[np.isnan(psf)] = 0.0
        zarray = np.array([psf, xgrid, ygrid])
        zarray[0] /= zarray[0].sum()
        return zarray

    def gaussian_psf(self, xstd, ystd, theta, boxsize, scale):
        """
        Gaussian PSF profile image
        
        Parameters
        ----------
        xstd: float
            Standard deviation of the Gaussian in x before rotating by theta
        ystd: float
            Standard deviation of the Gaussian in y before rotating by theta
        theta: float
            Rotation angle in radians.
            The rotation angle increases counterclockwise
        boxsize: float
            Size of image on a side for Moffat profile
        scale: float
            Pixel scale for image
        
        Returns
        -------
        zarray: numpy 3d array
            An array with length 3 for the first axis: PSF image, xgrid, ygrid
        """
        M = Gaussian2D()
        M.x_stddev.value = xstd
        M.y_stddev.value = ystd
        M.theta.value = theta
        xl, xh = (0.0 - boxsize / 2.0, 0.0 + boxsize / 2.0 + scale)
        yl, yh = (0.0 - boxsize / 2.0, 0.0 + boxsize / 2.0 + scale)
        x, y = (np.arange(xl, xh, scale), np.arange(yl, yh, scale))
        xgrid, ygrid = np.meshgrid(x, y)
        zarray = np.array([M(xgrid, ygrid), xgrid, ygrid])
        zarray[0] /= zarray[0].sum()
        return zarray

    def moffat_psf(self, seeing, boxsize, scale, alpha=3.5):
        """
        Moffat PSF profile image
        
        Parameters
        ----------
        seeing: float
            FWHM of the Moffat profile
        boxsize: float
            Size of image on a side for Moffat profile
        scale: float
            Pixel scale for image
        alpha: float
            Power index in Moffat profile function
        
        Returns
        -------
        zarray: numpy 3d array
            An array with length 3 for the first axis: PSF image, xgrid, ygrid
        """
        M = Moffat2D()
        M.alpha.value = alpha
        M.gamma.value = 0.5 * seeing / np.sqrt(2 ** (1.0 / M.alpha.value) - 1.0)
        xl, xh = (0.0 - boxsize / 2.0, 0.0 + boxsize / 2.0 + scale)
        yl, yh = (0.0 - boxsize / 2.0, 0.0 + boxsize / 2.0 + scale)
        x, y = (np.arange(xl, xh, scale), np.arange(yl, yh, scale))
        xgrid, ygrid = np.meshgrid(x, y)

        Z = self.moffat_psf_integration(xgrid.ravel(), ygrid.ravel(),
                                        seeing,
                                        boxsize=boxsize+1.5,
                                        alpha=alpha)
        Z = np.reshape(Z, xgrid.shape)

        zarray = np.array([Z, xgrid, ygrid])
        #zarray[0] /= zarray[0].sum()
       
        return zarray

        
    def moffat_psf_integration(self, xloc, yloc, seeing, boxsize=14.,
                               scale=0.1, alpha=3.5):
        '''
        Based on Remedy extract.py code from G. Zeimann
        https://github.com/grzeimann/Remedy/blob/master/extract.py
        Moffat PSF profile image
        
        Parameters
        ----------
        seeing: float
            FWHM of the Moffat profile
        boxsize: float
            Size of image on a side for Moffat profile
        scale: float
            Pixel scale for image
        alpha: float
            Power index in Moffat profile function
        
        Returns
        -------
        zarray: numpy 3d array
            An array with length 3 for the first axis: PSF image, xgrid, ygrid
        '''
        M = Moffat2D()
        M.alpha.value = alpha
        M.gamma.value = 0.5 * seeing / np.sqrt(2**(1./ M.alpha.value) - 1.)
        xl, xh = (0. - boxsize / 2., 0. + boxsize / 2. + scale)
        yl, yh = (0. - boxsize / 2., 0. + boxsize / 2. + scale)
        x, y = (np.arange(xl, xh, scale), np.arange(yl, yh, scale))
        xgrid, ygrid = np.meshgrid(x, y)
        Z = M(xgrid, ygrid)
        Z = Z / Z.sum()
        V = xloc * 0.
        cnt = 0
        for xl, yl in zip(xloc, yloc):
            d = np.sqrt((xgrid-xl)**2 + (ygrid-yl)**2)
            sel = d <= 0.75
            adj = np.pi * 0.75**2 / (sel.sum() * scale**2)
            V[cnt] = np.sum(Z[sel]) * adj
            cnt += 1
        return V


    def model_psf(
        self,
        gmag_limit=21.0,
        radius=8.0,
        pixscale=0.25,
        boundary=[-21.0, 21.0, -21.0, 21.0],
        interp_kind="linear",
    ):
        """
        Model the VIRUS on-sky PSF for a set of three exposures
        
        Parameters
        ----------
        gmag_limit: float
            Only stars brighter than this value will be used to model the psf
        radius: float
            Radius in arcseconds used to collect fibers around a given coord
        pixscale: float
            Pixel scale in arcseconds of the psf image
        boundary: list of 4 values
            [x_lower, x_higher, y_lower, y_higher] limits for including a star
            in ifu coordinates
        interp_kind: str
            Kind of interpolation to pixelated grid from fiber intensity

        Returns
        -------
        zarray: numpy 3d array
            An array with length 3 for the first axis: PSF image, xgrid, ygrid
        """
        # Fit a box in a circle (0.95 factor for caution)
        boxsize = int(np.sqrt(2.0) * radius * 0.95 / pixscale) * pixscale

        coords, gmag, starid = self.get_starcatalog_params()

        psf_list = []
        for i, coord in enumerate(coords):
            # Only use stars that are bright enough and not near the edge
            if gmag[i] > gmag_limit:
                self.log.info(
                    "PSF model StarID: {} too faint: {:4.2f}".format(starid[i], gmag[i])
                )
                continue
            result = self.get_fiberinfo_for_coord(coord, radius=radius)
            if result is None:
                continue
            ifux, ifuy, xc, yc, ra, dec, data, datae, mask = result
            in_bounds = (
                (xc > boundary[0])
                * (xc < boundary[1])
                * (yc > boundary[2])
                * (yc < boundary[3])
            )
            if not in_bounds:
                self.log.info(
                    "PSF model StarID: {} on edge: {:4.2f}, {:4.2f}".format(starid[i], xc, yc)
                )
                continue
            psfi = self.make_collapsed_image(
                xc,
                yc,
                ifux,
                ifuy,
                data,
                mask,
                boxsize=boxsize,
                scale=pixscale,
                interp_kind=interp_kind,
            )
            psf_list.append(psfi)

        if len(psf_list) == 0:
            self.log.warning("No suitable stars for PSF")
            self.log.warning('Using default moffat PSF with 1.8" seeing')
            return self.moffat_psf(1.8, boxsize, pixscale)

        self.log.info("%i suitable stars for PSF" % len(psf_list))
        C = np.array(psf_list)
        avg_psf_image = np.median(C[:, 0, :, :], axis=0)
        avg_psf_image[np.isnan(avg_psf_image)] = 0.0
        zarray = np.array([avg_psf_image, C[0, 1], C[0, 2]])
        return zarray

    def make_collapsed_image(
        self,
        xc,
        yc,
        xloc,
        yloc,
        data,
        mask,
        scale=0.25,
        seeing_fac=1.8,
        boxsize=4.0,
        wrange=[3470, 5540],
        nchunks=11,
        convolve_image=False,
        interp_kind="linear",
    ):
        """
        Collapse spectra to make a single image on a rectified grid.  This
        may be done for a wavelength range and using a number of chunks
        of wavelength to take ADR into account.
        
        Parameters
        ----------
        xc: float
            The ifu x-coordinate for the center of the collapse frame
        yc: float
            The ifu y-coordinate for the center of the collapse frame
        xloc: numpy array
            The ifu x-coordinate for each fiber
        yloc: numpy array
            The ifu y-coordinate for each fiber
        data: numpy 2d array
            The calibrated spectra for each fiber
        mask: numpy 2d array
            The good fiber wavelengths to be used in collapsed frame
        scale: float
            Pixel scale for output collapsed image
        seeing_fac: float
            seeing_fac = 2.35 * radius of the Gaussian kernel used 
            if convolving the images to smooth out features. Unit: arcseconds
        boxsize: float
            Length of the side in arcseconds for the convolved image
        wrange: list
            The wavelength range to use for collapsing the frame
        nchunks: int
            Number of chunks used to take ADR into account when collapsing
            the fibers.  Use a larger number for a larger wavelength.
            A small wavelength may only need one chunk
        convolve_image: bool
            If true, the collapsed frame is smoothed at the seeing_fac scale
        interp_kind: str
            Kind of interpolation to pixelated grid from fiber intensity
        
        Returns
        -------
        zarray: numpy 3d array
            An array with length 3 for the first axis: PSF image, xgrid, ygrid
        """
        a, b = data.shape
        N = int(boxsize / scale)
        xl, xh = (xc - boxsize / 2.0, xc + boxsize / 2.0)
        yl, yh = (yc - boxsize / 2.0, yc + boxsize / 2.0)
        x, y = (np.linspace(xl, xh, N), np.linspace(yl, yh, N))
        xgrid, ygrid = np.meshgrid(x, y)
        S = np.zeros((a, 2))
        area = np.pi * 0.75 ** 2
        sel = (self.wave > wrange[0]) * (self.wave <= wrange[1])
        I = np.arange(b)
        ichunk = [np.mean(xi) for xi in np.array_split(I[sel], nchunks)]
        ichunk = np.array(ichunk, dtype=int)
        cnt = 0
        image_list = []
        if convolve_image:
            seeing = seeing_fac / scale
            G = Gaussian2DKernel(seeing / 2.35)
        if interp_kind not in ["linear", "cubic"]:
            self.log.warning('interp_kind must be "linear" or "cubic"')
            self.log.warning('Using "linear" for interp_kind')
            interp_kind = "linear"

        for chunk, mchunk in zip(
            np.array_split(data[:, sel], nchunks, axis=1),
            np.array_split(mask[:, sel], nchunks, axis=1),
        ):
            marray = np.ma.array(chunk, mask=mchunk < 1e-8)
            image = np.ma.median(marray, axis=1)
            image = image / np.ma.sum(image)
            S[:, 0] = xloc - self.ADRx[ichunk[cnt]]
            S[:, 1] = yloc - self.ADRy[ichunk[cnt]]
            cnt += 1
            grid_z = (
                griddata(
                    S[~image.mask],
                    image.data[~image.mask],
                    (xgrid, ygrid),
                    method=interp_kind,
                )
                * scale ** 2
                / area
            )
            if convolve_image:
                grid_z = convolve(grid_z, G)
            image_list.append(grid_z)

        image = np.median(image_list, axis=0)
        image[np.isnan(image)] = 0.0
        zarray = np.array([image, xgrid - xc, ygrid - yc])
        return zarray

    def make_narrowband_image(
        self,
        xc,
        yc,
        xloc,
        yloc,
        data,
        mask,
        error=None,
        scale=0.25,
        seeing_fac=1.8,
        boxsize=4.0,
        wrange=[3470, 5540],
        convolve_image=False,
        interp_kind="linear",
        fill_value=0.0,
        bitmask=False,
    ):
        """
        Sum spectra across a wavelength range or filter to make a single image
        on a rectified grid.  ADR will be calculated at mean wavelength in wrange.
        
        Parameters
        ----------
        xc: float
            The ifu x-coordinate for the center of the collapse frame
        yc: float
            The ifu y-coordinate for the center of the collapse frame
        xloc: numpy array
            The ifu x-coordinate for each fiber
        yloc: numpy array
            The ifu y-coordinate for each fiber
        data: numpy 2d array
            The calibrated spectra for each fiber
        error: numpy 2d array
            Optional 2D array with fiber errors 
        mask: numpy 2d array
            The good fiber wavelengths to be used in collapsed frame
        scale: float
            Pixel scale for output collapsed image
        seeing_fac: float
            seeing_fac = 2.35 * radius of the Gaussian kernel used
            if convolving the images to smooth out features. Unit: arcseconds
        boxsize: float
            Length of the side in arcseconds for the convolved image
        wrange: list
            The wavelength range to use for collapsing the frame
        convolve_image: bool
            If true, the collapsed frame is smoothed at the seeing_fac scale
        interp_kind: str
            Kind of interpolation to pixelated grid from fiber intensity.
            Options are 'linear', 'cubic', 'nearest'. Default is linear.
        fill_value: float, optional
            Value used to fill in for requested points outside of coverage or in a mask
            region. If not provided, then the default is nan.
        bitmask: bool
            Option to include an extension with interpolated bitmask array for slice.
            mask parameter must be in bitmask format.

        Returns
        -------
        zarray: numpy 3d array
            An array with length 3 for the first axis: wavelength summed image
            across wrange in units of 10^-17/ergs/cm^2
        xgrid, ygrid : numpy 3d array
            in relative arcsec from center coordinates
        
        """

        if bitmask:
            if error is None:
                print('You must provide an error and bitmask array for bitmask=True option')
                return None
                
        a, b = data.shape
        N = int(boxsize / scale)
        xl, xh = (xc - boxsize / 2.0, xc + boxsize / 2.0)
        yl, yh = (yc - boxsize / 2.0, yc + boxsize / 2.0)
        x, y = (np.linspace(xl, xh, N), np.linspace(yl, yh, N))
        xgrid, ygrid = np.meshgrid(x, y)
        S = np.zeros((a, 2))
        area = np.pi * 0.75 ** 2
        sel = (self.wave > wrange[0]) * (self.wave <= wrange[1])

        I = np.arange(b)

        if convolve_image:
            seeing = seeing_fac / scale
            G = Gaussian2DKernel(seeing / 2.35)

        if interp_kind not in ["linear", "cubic", "nearest"]:
            self.log.warning('interp_kind must be "linear", "nearest" or "cubic"')
            self.log.warning('Using "linear" for interp_kind')
            interp_kind = "linear"

        if bitmask:
            # only apply native masking
            image_mask = mask[:, sel] == 1
        else:
            image_mask = mask[:, sel] < 1e-8

        marray = np.ma.array(data[:, sel], mask=image_mask)
        image = 2.*np.ma.sum(marray, axis=1)  # multiply for 2AA bins

        if error is not None:
            marray = np.ma.array(error[:, sel], mask=image_mask)
            error_image = np.sqrt(2.*np.ma.sum(marray**2, axis=1))
            # multiply for 2AA bins. Add error in quadrature
        if bitmask:
            bitmask_image = np.max( mask[:, sel], axis=1 )
        
        S[:, 0] = xloc - np.mean(self.ADRx[sel])
        S[:, 1] = yloc - np.mean(self.ADRy[sel])

        if np.sum( image.data[~image.mask]) == 0:
            grid_z = np.zeros_like(xgrid)
        else:
            try:
                grid_z = (
                    griddata(
                        S[~image.mask],
                        image.data[~image.mask],
                        (xgrid, ygrid),
                        method=interp_kind,
                        fill_value=0.0, # 2025-01-21 EMC fix issue where all 0s are returned
                    )
                    * scale ** 2
                    / area
                )
            except Exception:
                grid_z = np.zeros_like(xgrid)

        # convert back to nan
        grid_z[ grid_z==0.0] = np.nan
                
        # propogate mask through to new grid image. Added by EMC 2024-04-16, force to be a bool
        grid_mask = griddata(S, ~image.mask, (xgrid, ygrid), method='nearest').astype(bool)
        
        grid_z[~grid_mask] = np.nan

        if error is not None:

            if np.sum( error_image.data[~error_image.mask]) == 0:
                grid_z_error = np.zeros_like(xgrid)
            else:
                try:
                    grid_z_error = (
                        griddata(
                            S[~error_image.mask],
                            error_image.data[~error_image.mask],
                            (xgrid, ygrid),
                            method=interp_kind,
                            fill_value=0.0, # 2025-01-21 fix issue where all 0s are returned
                        )
                        * scale ** 2
                        / area
                    )
                except:
                    grid_z_error = np.zeros_like(xgrid)
                
            grid_z_error[ grid_z_error<=0.0] = np.nan

            # mask image Added by EMC 2024-04-16
            grid_z_error[~grid_mask] = np.nan
        if bitmask:
            grid_bitmask = (
                griddata(
                    S,
                    bitmask_image,
                    (xgrid, ygrid),
                    method='nearest',
                    fill_value=-1, # no mask value where we have no coverage
                )
            ).astype(int)
                       
        if convolve_image:
            grid_z = convolve(grid_z, G)
            if error is not None:
                grid_z_error = convolve(grid_z_error, G)

        if error is None:
            image = grid_z
            #Added by EMC 2024-04-16 
            image[np.isnan(image)] = fill_value
            zarray = np.array([image, xgrid - xc, ygrid - yc])
        else:
            if bitmask:
                image = grid_z

                image[np.isnan(image)] = fill_value
                image_error = grid_z_error
                image_error[np.isnan(image_error)] = fill_value

                zarray = np.array([image, image_error, grid_bitmask, xgrid - xc, ygrid - yc])
            else:
                
                image = grid_z
                #Added by EMC 2024-04-16 
                image[np.isnan(image)] = fill_value
                image_error = grid_z_error
                image_error[np.isnan(image_error)] = fill_value
                
                zarray = np.array([image, image_error, xgrid - xc, ygrid - yc])
        
        return zarray

    def get_psf_curve_of_growth(self, psf):
        """
        Analyse the curve of growth for an input psf
        
        Parameters
        ----------
        psf: numpy 3d array
            zeroth dimension: psf image, xgrid, ygrid
        
        Returns
        -------
        r: numpy array
            radius in arcseconds
        cog: numpy array
            curve of growth for psf
        """
        r = np.sqrt(psf[1] ** 2 + psf[2] ** 2).ravel()
        inds = np.argsort(r)
        maxr = psf[1].max()
        sel = r[inds] <= maxr
        cog = np.cumsum(psf[0].ravel()[inds][sel]) / np.sum(psf[0].ravel()[inds][sel])
        return r[inds][sel], cog

    
    
    def build_weights(self, xc, yc, ifux, ifuy, psf, I = None, fac = None,
                      return_I_fac = False):
        """
        Build weight matrix for spectral extraction.
        
        Parameters
        ----------
        xc: float
            The ifu x-coordinate for the center of the collapse frame
        yc: float 
            The ifu y-coordinate for the center of the collapse frame
        xloc: numpy array
            The ifu x-coordinate for each fiber
        yloc: numpy array
            The ifu y-coordinate for each fiber
        psf: numpy 3d array
            zeroth dimension: psf image, xgrid, ygrid
        I : callable (Optional)
            a function that returns the PSF function
            for a given x, y. If this is passed the
            psf image is ignored, otherwise it is 
            computed (default: None).
        fac : float (Optional)
            scale factor to account for area
            differences between the PSF image and
            the fibers. If None this is computed 
            from the PSF grid (default: None)
        return_I_fac : bool (Optional)
            return the computed (or input) I and fac 
            values (see above). This saves compute time 
            if you call this function many times with 
            the same PSF, as you can pass these values 
            back in the next time you run. 

        Returns
        -------
        weights: numpy 2d array (len of fibers by wavelength dimension)
            Weights for each fiber as function of wavelength for extraction
        I: PSF function
            returns the PSF function to pass through for next iteration if
            return_I_fac is set to True.
        """
        SX = np.zeros(len(ifux))
        SY = np.zeros(len(ifuy))
        weights = np.zeros((len(ifux), len(self.wave)))

        if not I:
            T = np.array([psf[1].ravel(), psf[2].ravel()]).swapaxes(0, 1)
            I = LinearNDInterpolator(T, psf[0].ravel(), fill_value=0.0)

        if not fac:
            scale = np.abs(psf[1][0, 1] - psf[1][0, 0])
            #area = 0.75 ** 2 * np.pi
            #area = 1.7671458676442586
            #fac = area / scale ** 2
            # scale to fiber size already included in moffat_psf_integration?
            fac = 1.0

        # Avoid using a loop to speed things up, uses more memory though
        SX = np.tile(ifux, len(self.wave)).reshape(len(self.wave), len(ifux)).T
        SY = np.tile(ifuy, len(self.wave)).reshape(len(self.wave), len(ifuy)).T
        ADRx3D = np.repeat(self.ADRx, len(ifux)).reshape(len(self.wave), len(ifux)).T
        ADRy3D = np.repeat(self.ADRy, len(ifuy)).reshape(len(self.wave), len(ifuy)).T

        SX = SX - ADRx3D - xc
        SY = SY - ADRy3D - yc
        weights = fac*I(SX, SY)
  
        if not return_I_fac:
            return weights
        else:
            return weights, I, fac


    def build_weights_old(self, xc, yc, ifux, ifuy, psf):
        """
        Build weight matrix for spectral extraction
        
        Parameters
        ----------
        xc: float
            The ifu x-coordinate for the center of the collapse frame
        yc: float 
            The ifu y-coordinate for the center of the collapse frame
        xloc: numpy array
            The ifu x-coordinate for each fiber
        yloc: numpy array
            The ifu y-coordinate for each fiber
        psf: numpy 3d array
            zeroth dimension: psf image, xgrid, ygrid
        
        Returns
        -------
        weights: numpy 2d array (len of fibers by wavelength dimension)
            Weights for each fiber as function of wavelength for extraction
        I, fac : 
            only returned if return_I_fac = True, see description above.
        """
        S = np.zeros((len(ifux), 2))
        T = np.array([psf[1].ravel(), psf[2].ravel()]).swapaxes(0, 1)
        I = LinearNDInterpolator(T, psf[0].ravel(), fill_value=0.0)
        weights = np.zeros((len(ifux), len(self.wave)))
        scale = np.abs(psf[1][0, 1] - psf[1][0, 0])
        area = 0.75 ** 2 * np.pi
        for i in np.arange(len(self.wave)):
            S[:, 0] = ifux - self.ADRx[i] - xc
            S[:, 1] = ifuy - self.ADRy[i] - yc
            # don't rescale by area factor - already done in moffat_psf_integration?
            weights[:, i] = I(S[:, 0], S[:, 1]) #* area / scale ** 2

        return weights

    def get_spectrum(self, data, error, mask, weights, remove_low_weights = True,
                     sclean_bad = True, return_scleaned_mask = False):
        """
        Weighted spectral extraction
        
        Parameters
        ----------
        data: numpy 2d array (number of fibers by wavelength dimension)
            Flux calibrated spectra for fibers within a given radius of the 
            source.  The radius is defined in get_fiberinfo_for_coord().
        error: numpy 2d array
            Error of the flux calibrated spectra 
        mask: numpy 2d array (bool)
            Mask of good wavelength regions for each fiber
        weights: numpy 2d array (float)
            Normalized extraction model for each fiber
        remove_low_weights : bool 
            remove the contribution from fibers with
            weights less than 10% the median 
            (default: True).
        sclean_bad : bool
            replace bad data in the
            way described by Karl Gebhardt and
            implemented in his `fitradecsp` code
        return_scleaned_mask : bool
            if true return a mask the same shape
            as the data if any of the        

        Returns
        -------
        spectrum: numpy 1d array
            Flux calibrated extracted spectrum
        spectrum_error: numpy 1d array
            Error for the flux calibrated extracted spectrum
        scleaned : numpy 2d array
            Only returned if `return_scleaned_mask` is true.
            A mask which is 1 in elements where sclean replaced
            bad data
        """

        scleaned = np.zeros_like(data)

        # Replace bad data with an average
        if sclean_bad:
            for i in range(data.shape[0]):
                data[i, :], error[i, :], mask[i, :], scleaned[i, :] = sclean(self.wave, data[i,:], error[i,:], mask[i, :])

        spectrum = np.sum(data * mask * weights, axis=0) / np.sum(
            mask * weights ** 2, axis=0
        )
        spectrum_error = np.sqrt(
            np.sum(error ** 2 * mask * weights, axis=0)
            / np.sum(mask * weights ** 2, axis=0)
        )

        # Only use wavelengths with enough weight to avoid large noise spikes
        if remove_low_weights:
            w = np.sum(mask * weights ** 2, axis=0)
            sel = w < 0.05
            spectrum[sel] = np.nan
            spectrum_error[sel] = np.nan

        if return_scleaned_mask:   
            return spectrum, spectrum_error, scleaned, data, error, mask
        else:
            return spectrum, spectrum_error

    def close(self):
        """
        Close open shot h5 file if open
        """
        
        if self.shoth5 is not None:
            self.shoth5.close()
            
        if self.fibers is not None:
            self.fibers.close()

