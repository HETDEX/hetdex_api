# -*- coding: utf-8 -*-
"""

Initiates the Detections class.


Created on 2019/01/28

@author: Erin Mentuch Cooper
"""

from __future__ import print_function
from __future__ import unicode_literals

import sys
import os.path as op
import numpy as np
import tables as tb
import copy

np.warnings.filterwarnings("ignore")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from astropy.table import vstack, Table, Column, join
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
import pickle
import speclite.filters

from hetdex_api.survey import Survey
from hetdex_api.config import HDRconfig
from hetdex_api.mask import *

PYTHON_MAJOR_VERSION = sys.version_info[0]
PYTHON_VERSION = sys.version_info

try:
    from hetdex_api.config import HDRconfig

    LATEST_HDR_NAME = HDRconfig.LATEST_HDR_NAME

except Exception as e:
    print("Warning! Cannot find or import HDRconfig from hetdex_api!!", e)
    LATEST_HDR_NAME = "hdr2.1"


class Detections:
    def __init__(
        self,
        survey=LATEST_HDR_NAME,
        catalog_type="lines",
        curated_version=None,
        loadtable=True,
    ):
        """
        Initialize the detection catalog class for a given data release

        Input
        -----
        survey : string
            Data release you would like to load, i.e., 'hdr2','HDR1'
            This is case insensitive.
        catalog_type : string
            Catalog to laod up. Either 'lines' or 'continuum'. Default is 
            'lines'.
        load_table : bool
           Boolean flag to load all detection table info upon initialization.
           For example, if you just want to grab a spectrum this isn't needed.
        
        """
        survey_options = ["hdr1", "hdr2", "hdr2.1"]
        catalog_type_options = ["lines", "continuum", "broad"]

        if survey.lower() not in survey_options:
            print("survey not in survey options")
            print(survey_options)
            return None

        if catalog_type.lower() not in catalog_type_options:
            print("catalog_type not in catalog_type options")
            print(catalog_type_options)
            return None

        # store to class
        if curated_version is not None:
            self.version = curated_version
            self.loadtable = False
            self.survey = "hdr" + curated_version[0:3]
        else:
            self.version = None
            self.survey = survey
            self.loadtable = loadtable

        global config
        config = HDRconfig(survey=self.survey)

        if catalog_type == "lines":
            self.filename = config.detecth5
        elif catalog_type == "continuum":
            self.filename = config.contsourceh5
        elif catalog_type == "broad":
            try:
                self.filename = config.detectbroadh5
            except:
                print("Could not locate broad line catalog")

        self.hdfile = tb.open_file(self.filename, mode="r")

        # store to class
        if curated_version is not None:
            self.version = curated_version
            self.loadtable = False
            self.survey = "hdr" + curated_version[0:3]
        else:
            self.survey = survey
            self.loadtable = loadtable

        if self.version is not None:

            try:
                catfile = op.join(
                    config.detect_dir, "catalogs", "detect_hdr" + self.version + ".fits"
                )
                det_table = Table.read(catfile)

                for col in det_table.colnames:
                    if isinstance(det_table[col][0], str):
                        setattr(self, col, np.array(det_table[col]).astype(str))
                    else:
                        setattr(self, col, np.array(det_table[col]))

                self.vis_class = -1 * np.ones(np.size(self.detectid))

            except:
                print("Could not open curated catalog version: " + self.version)
                return None

        elif self.loadtable:
            colnames = self.hdfile.root.Detections.colnames
            for name in colnames:
                if isinstance(
                    getattr(self.hdfile.root.Detections.cols, name)[0], np.bytes_
                ):
                    setattr(
                        self,
                        name,
                        getattr(self.hdfile.root.Detections.cols, name)[:].astype(str),
                    )
                else:
                    setattr(
                        self, name, getattr(self.hdfile.root.Detections.cols, name)[:]
                    )
            if self.survey == "hdr2.1":
                # Fix fluxes and continuum values for aperture corrections  
                wave = self.hdfile.root.Detections.cols.wave[:]
                apcor = self.hdfile.root.Spectra.cols.apcor[:]
                wave_spec = self.hdfile.root.Spectra.cols.wave1d[:]

                apcor_array = np.ones_like(wave)
                for idx in np.arange(0, np.size(wave)):
                    sel_apcor = np.where(wave_spec[idx, :] > wave[idx])[0][0]
                    apcor_array[idx]=apcor[idx, sel_apcor]

                self.flux /= apcor_array
                self.flux_err /= apcor_array
                self.continuum /= apcor_array
                self.continuum_err /= apcor_array

            # add in the elixer probabilties and associated info:
            if self.survey == "hdr1" and catalog_type == "lines":

                self.hdfile_elix = tb.open_file(config.elixerh5, mode="r")
                colnames2 = self.hdfile_elix.root.Classifications.colnames
                for name2 in colnames2:
                    if name2 == "detectid":
                        setattr(
                            self,
                            "detectid_elix",
                            self.hdfile_elix.root.Classifications.cols.detectid[:],
                        )
                    else:
                        if isinstance(
                            getattr(self.hdfile_elix.root.Classifications.cols, name2)[
                                0
                            ],
                            np.bytes_,
                        ):
                            setattr(
                                self,
                                name2,
                                getattr(
                                    self.hdfile_elix.root.Classifications.cols, name2
                                )[:].astype(str),
                            )
                        else:
                            setattr(
                                self,
                                name2,
                                getattr(
                                    self.hdfile_elix.root.Classifications.cols, name2
                                )[:],
                            )
            else:

                # add elixer info if node exists
                try:
                    colnames = self.hdfile.root.Elixer.colnames
                    for name in colnames:
                        if name == "detectid":
                            continue
                        if isinstance(
                            getattr(self.hdfile.root.Elixer.cols, name)[0], np.bytes_
                        ):
                            setattr(
                                self,
                                name,
                                getattr(self.hdfile.root.Elixer.cols, name)[:].astype(
                                    str
                                ),
                            )
                        else:
                            setattr(
                                self,
                                name,
                                getattr(self.hdfile.root.Elixer.cols, name)[:],
                            )
                    self.gmag = self.mag_sdss_g
                    self.gmag_err = self.mag_sdss_g
                except:
                    print("No Elixer table found")

            # also assign a field and some QA identifiers
            self.field = np.chararray(np.size(self.detectid), 12, unicode=True)
            self.fwhm = np.zeros(np.size(self.detectid))
            if self.survey == "hdr1":
                self.fluxlimit_4550 = np.zeros(np.size(self.detectid))
            else:
                self.fluxlimit_4540 = np.zeros(np.size(self.detectid))

            self.throughput = np.zeros(np.size(self.detectid))
            self.n_ifu = np.zeros(np.size(self.detectid), dtype=int)

            S = Survey(self.survey)

            for index, shot in enumerate(S.shotid):
                ix = np.where(self.shotid == shot)
                self.field[ix] = S.field[index].astype(str)
                # NOTE: python2 to python3 strings now unicode
                if self.survey == "hdr1":
                    self.fwhm[ix] = S.fwhm_moffat[index]
                    self.fluxlimit_4550[ix] = S.fluxlimit_4550[index]
                else:
                    self.fwhm[ix] = S.fwhm_virus[index]
                try:
                    self.fluxlimit_4540[ix] = S.fluxlimit_4540[index]
                except:
                    pass
                self.throughput[ix] = S.response_4540[index]
                self.n_ifu[ix] = S.n_ifu[index]

                # assign a vis_class field for future classification
                # -2 = ignore (bad detectid, shot)
                # -1 = no assignemnt
                # 0 = artifact
                # 1 = OII emitter
                # 2 = LAE emitter
                # 3 = star
                # 4 = nearby galaxies (HBeta, OIII usually)
                # 5 = other line
            # close the survey HDF5 file
            S.close()

            self.vis_class = -1 * np.ones(np.size(self.detectid))

            if self.survey == "hdr1":
                self.add_hetdex_gmag(loadpickle=True, picklefile=config.gmags)

            if self.survey == "hdr1":
                if PYTHON_MAJOR_VERSION < 3:
                    self.plae_poii_hetdex_gmag = np.array(
                        pickle.load(open(config.plae_poii_hetdex_gmag, "rb"))
                    )
                else:
                    self.plae_poii_hetdex_gmag = np.array(
                        pickle.load(
                            open(config.plae_poii_hetdex_gmag, "rb"), encoding="bytes"
                        )
                    )

        else:
            # just get coordinates and detectid
            self.detectid = self.hdfile.root.Detections.cols.detectid[:]
            self.ra = self.hdfile.root.Detections.cols.ra[:]
            self.dec = self.hdfile.root.Detections.cols.dec[:]

        # set the SkyCoords
        self.coords = SkyCoord(self.ra * u.degree, self.dec * u.degree, frame="icrs")

    def __getitem__(self, indx):
        """ 
        This allows for slicing of the Detections class
        object so that a mask can be applied to
        every attribute automatically by:
        
        detections_sliced = detects[indx]
        
        """

        p = copy.copy(self)
        attrnames = self.__dict__.keys()
        for attrname in attrnames:
            try:
                setattr(p, attrname, getattr(self, attrname)[indx])
            except:
                setattr(p, attrname, getattr(self, attrname))
        return p

    def refine(self, gmagcut=None, remove_large_gal=True, d25scale=3.0):
        """
        Masks out bad and bright detections 
        and returns a refined Detections class
        object

        gmagcut = mag limit to exclude everything
                  brighter, defaults to None
        """

        mask1 = self.remove_bad_amps()
        mask2 = self.remove_bad_detects()
        mask3 = self.remove_bright_stuff(gmagcut)
        mask4 = self.remove_bad_pix()
        mask5 = self.remove_shots()
        mask6 = self.remove_meteors()

        mask = mask1 * mask2 * mask3 * mask4 * mask5 * mask6

        if remove_large_gal:
            galmask = self.remove_large_gal(d25scale=d25scale)
            mask = mask * galmask

        return self[mask]

    def query_by_coords(self, coords, radius):
        """
        Returns mask based on a coordinate search
        
        self = Detections Class object   
        coords - astropy coordinate object
        radius - an astropy Quantity object, or a string 
        that can be parsed into one.  e.g., '1 degree' 
        or 1*u.degree. Will assume arcsec if no units given
        
        """
        sep = self.coords.separation(coords)
        try:
            maskcoords = sep < radius
        except:
            maskcoords = sep.arcmin < radius
        return maskcoords

    def find_match(
        self, coord, radius=5.0 * u.arcsec, wave=None, dwave=5.0, shotid=None
    ):
        """
        Function to cross match another line detection

        Parameters
        ----------
        coord
            an astropy coordinates object
        wave
            central wavelength in AA you want to search. If
            nothing is given, it will search without any
            wavelength contraint
        radius
            search radius. An astropy quantity
        dwave
            delta wavelength to search
        shotid
            optional shotid for a specific observation
        
        Returns
        -------
        match_index
            index of matches
        """

        selmatch = self.query_by_coords(coord, radius)

        if wave is not None:
            selwave = np.abs(self.wave - wave) < dwave
            selmatch = selwave * selmatch

        if shotid is not None:
            selshot = self.shotid == shotid
            selmatch = selshot * selmatch

        return selmatch

    def query_by_dictionary(self, limits):
        """
        Takes a dictionary of query limits
        and reduces the detections database. This 
        can either be taken from the dictionary saved
        by detwidgets.py GUI or can take the following 
        form:


        class Det_limits:
        # limits to be returned to query detections                                    
        def __init__(self):
          self.wave_low = None
          self.wave_high = None
          self.flux_low = None
          self.flux_high = None
          self.linewidth_low = None
          self.linewidth_high = None
          self.sn_low = None
          self.sn_high = None
          self.chi2_low = None
          self.chi2_high = None
          self.cont_low = None
          self.cont_high = None
          self.aperture_flag = False
          self.ra = None
          self.dec = None
          self.rad = None
          self.field = None


        Dictionary Description
        aperture_flat = when True, will query for defined aperture
        ra = right ascension for aperture in degrees
        dec = declination of aperture in degrees
        rad = radius of aperture in arcmin
        field =  ('all', 'dex-spring', 'dex-fall', 'cosmos', 'egs', 'goods-n', 'other')
        
        others should be obvious
        """

        ndets = np.size(self.detectid)

        if limits.wave_low or limits.wave_high:
            maskwave = (self.wave > limits.wave_low) * (self.wave < limits.wave_high)
        else:
            maskwave = np.ones(ndets, dtype=bool)

        if limits.flux_low or limits.flux_high:
            maskflux = (self.flux > limits.flux_low) * (self.flux < limits.flux_high)
        else:
            maskflux = np.ones(ndets, dtype=bool)

        if limits.linewidth_low or limits.linewidth_high:
            masklw = (self.linewidth > limits.linewidth_low) * (
                self.linewidth < limits.linewidth_high
            )
        else:
            masklw = np.ones(ndets, dtype=bool)

        if limits.sn_low or limits.sn_high:
            masksn = (self.sn > limits.sn_low) * (self.sn < limits.sn_high)
        else:
            masksn = np.ones(ndets, dtype=bool)

        if limits.chi2_low or limits.chi2_high:
            maskchi2 = (self.chi2 > limits.chi2_low) * (self.chi2 < limits.chi2_high)
        else:
            maskchi2 = np.ones(ndets, dtype=bool)

        if limits.cont_low or limits.cont_high:
            maskcont = (self.continuum > limits.cont_low) * (
                self.continuum < limits.cont_high
            )
        else:
            maskcont = np.ones(ndets, dtype=bool)

        if limits.aperture_flag:
            coords = SkyCoord(limits.ra * u.degree, limits.dec * u.degree, frame="icrs")
            maskfield = self.query_by_coords(coords, limits.rad)
        else:
            maskfield = np.zeros(ndets, dtype=bool)
            print("Subselecting for field(s):", limits.field)

            for field_index in limits.field:
                if field_index == "all":
                    print("Field = 'all'; not downselecting")
                    maskfield = np.ones(ndets, dtype=bool)
                else:
                    if isinstance(field_index, str):
                        # python2 to python3 issue (pytables also ... as bytes vs unicode)
                        mask_i = self.field.decode() == field_index
                    else:
                        mask_i = self.field == field_index
                    maskfield = np.logical_or(maskfield, mask_i)

        mask = maskwave * masklw * masksn * maskflux * maskchi2 * maskcont * maskfield

        return mask

    def query_by_pickle(self, picklefile):
        """
        
        this function queries the Detections class
        object based on query made with detwidgets.py
        
        Input
        
        self = a detections class object
        picklefile = string filename for a pickle
        created in detwidget.py
        
        """
        # if locally encoded and opened (same version of python)
        if PYTHON_MAJOR_VERSION < 3:
            limits = pickle.load(open(picklefile, "rb"))
        else:
            limits = pickle.load(open(picklefile, "rb"), encoding="bytes")

        mask = self.query_by_dictionary(limits)
        return mask

    def remove_bad_detects(self):
        """
        Reads in the bad detect list from config.py
        and removes those detectids
        Set these to -2 (they are software errors)
        Don't use for machine learning or other
        classifying algorithms tests.
        """
        global config
        # set an empty mask to start
        mask = np.zeros(np.size(self.detectid), dtype=bool)

        baddetects = np.loadtxt(config.baddetect, dtype=int)

        for baddet in baddetects:
            maskdet = self.detectid == baddet
            mask = np.logical_or(mask, maskdet)

        return np.invert(mask)

    def remove_bad_amps(self):
        """
        Reads in the bad amp list from config.py
        and creates a mask to remove those detections.
        It will also assign a 0 to signify an artifact
        in vis_class
        """
        global config
        # set an empty mask to start

        mask = np.zeros(np.size(self.detectid), dtype=bool)

        if self.survey == "hdr1":
            badamps1 = ascii.read(
                config.badamp, names=["ifuslot", "amp", "date_start", "date_end"]
            )

            badamps2 = ascii.read(
                config.badamp, names=["ifuslot", "amp", "date_start", "date_end"]
            )

            badamps = vstack([badamps1, badamps2])

            if self.survey == "hdr2":
                self.date = (self.shotid / 1000).astype(int)

            for row in np.arange(np.size(badamps)):
                if badamps["amp"][row] == "AA":
                    maskamp = (
                        (self.ifuslot == str(badamps["ifuslot"][row]).zfill(3))
                        * (self.date >= badamps["date_start"][row])
                        * (self.date <= badamps["date_end"][row])
                    )
                    mask = np.logical_or(mask, maskamp)
                else:
                    maskamp = (
                        (self.amp == badamps["amp"][row])
                        * (self.ifuslot == str(badamps["ifuslot"][row]).zfill(3))
                        * (self.date >= badamps["date_start"][row])
                        * (self.date <= badamps["date_end"][row])
                    )
                    mask = np.logical_or(mask, maskamp)

            return np.logical_not(mask)
        else:

            # first read in amp_flag.fits file
            badamps = Table.read(config.badamp)

            det_table = self.return_astropy_table()

            join_tab = join(
                det_table, badamps, keys=["shotid", "multiframe"], join_type="left"
            )
            # this is needed to match with detection class object sorting
            join_tab.sort("detectid")

            mask1 = join_tab["flag"] != 0

            del det_table, join_tab

            # add in any newly found badamps that haven't made it into the
            # amp_flag.fits file yet

            mask2 = np.zeros(np.size(self.detectid), dtype=bool)

            badamps2 = Table.read(config.badamp2, format="ascii")

            for row in badamps2:
                selmf = self.multiframe == row["multiframe"]
                seldate = (self.date >= row["date_start"]) * (
                    self.date <= row["date_end"]
                )
                mask2 = np.logical_or(mask2, selmf * seldate)

            mask = mask1 * np.logical_not(mask2)

            return mask

    def remove_bad_pix(self):
        """
        Takes the post-hdr1 list of bad pixels 
        in HETDEX_API/known_issues/hdr1/posthdr1badpix.list
        
        For current development we will use the one in EMC's 
        directory as that will be most up to date:
        
        /work/05350/ecooper/stampede2/HETDEX_API/known_issues/hdr1
        
        and removes any detections within a Detections() class 
        object in the defined regions when
        either the .refine() or .remove_bad_pix() methods are 
        called.
        
        Note: all previously know bad detections are stored in
        
        HETDEX_API/known_issues/hdr1/badpix.list 
        
        **** THESE SHOULD BE REMOVED WHEN USING THE SHOT H5 files
        from HDR1
        
        """

        if True:
            badpixlist = ascii.read(
                config.badpix, names=["multiframe", "x1", "x2", "y1", "y2"]
            )

            mask = np.zeros(np.size(self.detectid), dtype=bool)

            for row in badpixlist:
                maskbadpix = (
                    (self.multiframe == row["multiframe"])
                    * (self.x_raw >= row["x1"])
                    * (self.x_raw <= row["x2"])
                    * (self.y_raw >= row["y1"])
                    * (self.y_raw <= row["y2"])
                )
                mask = np.logical_or(maskbadpix, mask)

            self.vis_class[mask] = 0

        else:  # except:
            mask = np.zeros(np.size(self.detectid), dtype=bool)

        return np.invert(mask)

    def remove_shots(self):
        """
        Takes a list of bad shots and removes them. Assigns -2 
        to detections in these shots so they are not used 
        in any MLing analysis
        """
        global config
        mask = np.zeros(np.size(self.detectid), dtype=bool)
        badshots = np.loadtxt(config.badshot, dtype=int)

        for shot in badshots:
            maskshot = self.shotid == shot
            mask = np.logical_or(maskshot, mask)

        self.vis_class[mask] = -2

        return np.invert(mask)

    def remove_balmerdip_stars(self):
        """
        Applies a cut to the databased to remove all
        stars that show a false emission feature around 3780
        Also assigns a star classification to each detectid

        This is obsolete as gmag cuts get rid of these easily
        """
        mask1 = (self.wave > 3775) * (self.wave < 3785) * (self.continuum > 3)

        mask2 = (self.wave > 4503.0) * (self.wave < 4513.0) * (self.continuum > 10)

        mask = np.logical_or(mask1, mask2)

        self.vis_class[mask] = 3

        return np.invert(mask)

    def remove_bright_stuff(self, gmagcut):
        """
        Applies a cut to remove bright stars based on
        gmag attribute. Assigns a star classification
        to each detectid, although some may be nearby
        galaxies. Will want to improve this if you are
        looking to study nearby galaxies.
        """

        if gmagcut:
            mask = self.gmag < gmagcut
        else:
            mask = np.zeros(np.size(self.detectid), dtype=bool)

        return np.invert(mask)

    def remove_ccd_features(self):
        """
        Remove all objects with very at the
        edges of the detectors and denotes 
        them as artifacts in vis_class
        """

        mask = (self.wave < 3540.0) | (self.wave > 5515.0)
        self.vis_class[mask] = 0

        return np.invert(mask)

    def remove_meteors(self):
        """
        Returns boolean mask with detections landing on meteor
        streaks masked. Use np.invert(mask) to find meteors
        """

        global config

        met_tab = Table.read(config.meteor, format="ascii")

        mask = np.ones_like(self.detectid, dtype=bool)

        for row in met_tab:
            sel_shot = np.where(self.shotid == row["shotid"])[0]
            for idx in sel_shot:
                coord = self.coords[idx]
                mask[idx] = meteor_flag_from_coords(coord, row["shotid"])
        return mask

    def remove_large_gal(self, d25scale=3.0):
        """
        Returns boolean mask with detections landing within
        galaxy defined by d25scale flagged as False.

        Based on check_all_large_gal from hetdex_tools/galmask.py
        written by John Feldmeier
        """

        global config

        galaxy_cat = Table.read(config.rc3cat, format="ascii")

        mask = np.ones_like(self.detectid, dtype=bool)

        # Loop over each galaxy

        for idx, row in enumerate(galaxy_cat):
            gal_coord = SkyCoord(row["Coords"], frame="icrs")
            rlimit = 1.1*d25scale*row["SemiMajorAxis"]*u.arcmin

            if np.any(self.coords.separation(gal_coord) < rlimit):
                galregion = create_gal_ellipse(galaxy_cat,
                                               row_index=idx,
                                               d25scale=d25scale)
                dummy_wcs = create_dummy_wcs(galregion.center,
                                             imsize=2*galregion.height)
                galflag = galregion.contains(self.coords, dummy_wcs)
                mask = mask * np.invert(galflag)

        return mask

    def get_spectrum(self, detectid_i):
        """
        Grabs the 1D spectrum used to measure fitted parameters.
        """
        spectra = self.hdfile.root.Spectra
        spectra_table = spectra.read_where("detectid == detectid_i")
        data = Table()
        intensityunit = u.erg / (u.cm ** 2 * u.s * u.AA)
        data["wave1d"] = Column(spectra_table["wave1d"][0], unit=u.AA)
        data["spec1d"] = Column(
            spectra_table["spec1d"][0], unit=1.0e-17 * intensityunit
        )
        data["spec1d_err"] = Column(
            spectra_table["spec1d_err"][0], unit=1.0e-17 * intensityunit
        )

        # convert from 2AA binning to 1AA binning:
        data["spec1d"] /= 2.0
        data["spec1d_err"] /= 2.0

        return data

    def get_gband_mag(self, detectid_i):
        """
        Calculates the gband magnitude from the 1D spectrum
        """

        spec_table = self.get_spectrum(detectid_i)
        gfilt = speclite.filters.load_filters("sdss2010-g")
        flux, wlen = gfilt.pad_spectrum(
            np.array(1.0e-17 * spec_table["spec1d"]), np.array(spec_table["wave1d"])
        )
        mag = gfilt.get_ab_magnitudes(flux, wlen)

        return mag

    def add_hetdex_gmag(self, loadpickle=True, picklefile="gmags.pickle"):
        """
        Calculates g-band magnitude from spec1d for
        each detectid in the Detections class instance
        If the gmags.pickle file and pickle=True is 
        given then it will just load from previous computation
        """
        if loadpickle:
            # todo: encoding='latin1' is an assumption ... might be better to use bytes?
            if PYTHON_MAJOR_VERSION < 3:
                self.gmag = pickle.load(open(picklefile, "rb"))
            else:
                self.gmag = pickle.load(open(picklefile, "rb"), encoding="bytes")
        else:
            self.gmag = np.zeros(np.size(self.detectid), dtype=float)

            # fastest way is to calculate gmag over whole 1D spec array
            # then populate the detections class

            detectid_spec = self.hdfile.root.Spectra.cols.detectid[:]
            spec1d = self.hdfile.root.Spectra.cols.spec1d[:]
            # convert from ergs/s/cm2 to ergs/s/cm2/AA
            spec1d /= 2.0
            wave_rect = 2.0 * np.arange(1036) + 3470.0
            gfilt = speclite.filters.load_filters("sdss2010-g")
            flux, wlen = gfilt.pad_spectrum(1.0e-17 * spec1d, wave_rect)
            gmags = gfilt.get_ab_magnitudes(flux, wlen)

            for idx, detectid_i in enumerate(self.detectid[:]):
                seldet = np.where(detectid_spec == detectid_i)[0]
                if np.size(seldet) == 1:
                    self.gmag[idx] = gmags[seldet][0][0]
                else:
                    self.gmag[idx] = np.nan

    def get_hetdex_mag(self, detectid_i, filter="sdss2010-g"):
        """ 
        filter = can be any filter used in the speclite
                 package 
                 https://speclite.readthedocs.io/en/latest/api.html

        """
        spec_table = self.get_spectrum(detectid_i)
        filt = speclite.filters.load_filters(filter)
        flux, wlen = filt.pad_spectrum(
            np.array(1.0e-17 * spec_table["spec1d"]), np.array(spec_table["wave1d"])
        )
        mag = filt.get_ab_magnitudes(flux, wlen)[0][0]

        return mag

    def return_astropy_table(self):
        """
        Return an astropy table version of the Detections
        that can easily be saved

        Returns
        -------
        table : astropy.table:Table
            an astropy table you can save

        """
        table = Table()
        for name in self.hdfile.root.Detections.colnames:
            table[name] = getattr(self, name)

        table.add_column(Column(self.fwhm), index=1, name="fwhm")
        table.add_column(Column(self.throughput), index=2, name="throughput")
        table.add_column(Column(self.field), index=4, name="field")
        table.add_column(Column(self.n_ifu), index=5, name="n_ifu")

        if self.survey == "hdr1":
            table.add_column(
                Column(self.fluxlimit_4550), index=3, name="fluxlimit_4550"
            )
            table.add_column(Column(self.gmag), index=6, name="gmag")
            table.add_column(
                Column(self.plae_poii_hetdex_gmag), name="plae_poii_hetdex_gmag"
            )
            table.add_column(Column(self.plae_poii_cat), name="plae_poii_cat")
            table.add_column(Column(self.plae_poii_aperture), name="plae_poii_aperture")

        else:
            try:
                table.add_column(Column(self.gmag), index=6, name="gmag")
                table.add_column(Column(self.gmag_err), index=6, name="gmag_err")
                for name in self.hdfile.root.Elixer.colnames:
                    table[name] = getattr(self, name)
            except:
                print("Could not add elixer columns")
            try:
                table.add_column(
                    Column(self.fluxlimit_4540), index=3, name="fluxlimit_4540"
                )
            except:
                print("Could not add average flux limit")

        return table

    def save_spectrum(self, detectid_i, outfile=None):
        spec_data = self.get_spectrum(detectid_i)
        if outfile:
            ascii.write(spec_data, outfile, overwrite=True)
        else:
            ascii.write(spec_data, "spec_" + str(detectid_i) + ".dat", overwrite=True)

    def plot_spectrum(self, detectid_i, xlim=None, ylim=None):
        spec_data = self.get_spectrum(detectid_i)
        plt.figure(figsize=(8, 6))
        plt.errorbar(
            spec_data["wave1d"], spec_data["spec1d"], yerr=spec_data["spec1d_err"]
        )
        plt.title("DetectID " + str(detectid_i))
        plt.xlabel("wave (AA)")
        plt.ylabel("flux (1e-17 erg/s/cm^2/AA)")
        if xlim is not None:
            plt.xlim(xlim)
            if ylim is not None:
                plt.ylim(ylim)
        plt.show()

    def close(self):
        self.hdfile.close()
