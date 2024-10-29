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
from scipy import interpolate

from astropy.table import vstack, Table, Column, join
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
import pickle
import speclite.filters

from hetdex_api.shot import Fibers, get_fibers_table
from hetdex_api.survey import Survey
from hetdex_api.config import HDRconfig
from hetdex_api.mask import *
from hetdex_api.extinction import *
import extinction

import healpy as hp


PYTHON_MAJOR_VERSION = sys.version_info[0]
PYTHON_VERSION = sys.version_info

try:
    LATEST_HDR_NAME = HDRconfig.LATEST_HDR_NAME
    from dustmaps.config import config as dustmaps_config

    if dustmaps_config["data_dir"] is None:
        print("Populating dustmaps config with {}".format(config.dustmaps))
        dustmaps_config["data_dir"] = config.dustmaps
    from dustmaps.sfd import SFDQuery

except Exception as e:
    print("Warning! Cannot find or import HDRconfig from hetdex_api!!", e)
    LATEST_HDR_NAME = "hdr2.1"


class Detections:
    def __init__(
        self,
        survey=LATEST_HDR_NAME,
        catalog_type="lines",
        curated_version=None,
        loadtable=False,
        verbose=False,
        searchable=False,
        recentbadamp=False,
        pdr=False,
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
        curated_version: string
            upload the curated detection catalog. e.g. '2.1.3'. This is
            much faster to load.
        load_table : bool
           Boolean flag to load all detection table info upon initialization.
           For example, if you just want to grab a spectrum this isn't needed.
        searchable: bool
            Boolean flag, related to loadtable. If loadtable is True, searchable
            is ignored. If loadtable is False, searchable toggles whether to load
            all detectids, ra, dec, waves from the h5 file to make the detections
            searchable. This can be time consuming and should be False if you know
            the detecids and only want to retrieve individual detections.
        recentbadamp: bool
            add in any recent bad amps that have not yet been added to amp_flag.fits.
            This is quite slow on the full database so is only recommended when using
            curated_version option.
        pdr: bool
            if TRUE the curated version string applies to the PDR (Public Data Release) instead
            of the internal HDR (HETDEX Data Release)
        """
        survey_options = ["hdr1", "hdr2", "hdr2.1", "hdr3", "hdr4", "hdr5", "pdr1"]
        catalog_type_options = ["lines", "continuum", "broad", "index"]

        if survey.lower() not in survey_options:
            print("survey not in survey options")
            print(survey_options)
            return None

        if catalog_type.lower() not in catalog_type_options:
            print("catalog_type not in catalog_type options")
            print(catalog_type_options)
            return None
        # store flag to update recent badamp list
        self.recentbadamp = recentbadamp
        self.catalog_type = catalog_type.lower()
        # store to class
        if curated_version is not None:
            self.version = curated_version
            self.loadtable = False

            if pdr:
                self.survey = "pdr" + curated_version  # [0:3]
            else:
                if curated_version[0] == "3":
                    # for now I'm assuming we don't have to do hdr3.0
                    self.survey = "hdr3"
                elif curated_version[0] == "4":
                    self.survey = "hdr4"
                elif curated_version[0] == "5":
                    self.survey = "hdr5"
                else:
                    self.survey = "hdr" + curated_version  # [0:3]
        else:
            self.version = None
            self.survey = survey
            self.loadtable = loadtable
            self.searchable = searchable
        # create attributes to store masking catalogs for get_detection_flags later
        self.badamps = None
        self.badpix = None
        self.galaxy_cat = None
        self.badfib = None
        self.badfib_transient = None
        self.cal5200_tab = None
        self.cal5460_tab = None

        self.config = HDRconfig(survey=self.survey)

        if float(self.survey[3:]) >= 3.0:  # == "hdr3":
            # open wd correction curve and create fit for quick assignment
            try:
                self.wd_corr = Table.read(
                    self.config.wdcor, format="ascii.no_header", names=["wave", "corr"]
                )
                self.wd_corr_f = interpolate.interp1d(
                    self.wd_corr["wave"], self.wd_corr["corr"]
                )
            except:  # right now this is pretty fatal ... will fail later
                print(f"Could not locate WD correction file: {self.config.wdcor}")

        if catalog_type == "lines":
            self.filename = self.config.detecth5
        elif catalog_type == "continuum":
            self.filename = self.config.contsourceh5
        elif catalog_type == "broad":
            try:
                self.filename = self.config.detectbroadh5
            except:
                print("Could not locate broad line catalog")
        elif catalog_type == "index":
            if self.survey not in ["hdr4", "hdr5", "pdr1"]:
                print("The Detection Index file only exists beyond hdr3")
            else:
                self.filename = self.config.detectindexh5

        try:
            self.hdfile = tb.open_file(self.filename, mode="r")
        except (
            Exception
        ) as e:  # this is essentially fatal, set the file to None, but might still be able to get the survey
            print(f"Could not locate catalog file: {self.filename}")
            self.hdfile = None
            raise e

        self.surveyh5 = tb.open_file(self.config.surveyh5, "r")

        if self.version is not None:
            try:
                if self.survey == "hdr2.1":
                    catfile = op.join(
                        self.config.detect_dir,
                        "catalogs",
                        "detect_hdr{}.fits".format(self.version),
                    )
                else:
                    catfile = op.join(
                        self.config.hdr_dir[self.survey],
                        "catalogs",
                        "detect_hdr{}.fits".format(self.version),
                    )

                det_table = Table.read(catfile)

                for col in det_table.colnames:
                    if isinstance(det_table[col][0], str):
                        setattr(self, col, np.array(det_table[col]).astype(str))
                    else:
                        setattr(self, col, np.array(det_table[col]))

                self.vis_class = -1 * np.ones(np.size(self.detectid))
            except Exception:
                print("Could not open curated catalog version: " + self.version)
                return None

        elif self.loadtable:
            if catalog_type == "index":
                colnames = self.hdfile.root.DetectIndex.colnames

                for name in colnames:
                    if isinstance(
                        getattr(self.hdfile.root.DetectIndex.cols, name)[0], np.bytes_
                    ):
                        if name == "survey":
                            outname = "hdr"
                        else:
                            outname = name
                        setattr(
                            self,
                            outname,
                            getattr(self.hdfile.root.DetectIndex.cols, name)[:].astype(
                                str
                            ),
                        )
                    else:
                        setattr(
                            self,
                            name,
                            getattr(self.hdfile.root.DetectIndex.cols, name)[:],
                        )
            else:
                colnames = self.hdfile.root.Detections.colnames
                for name in colnames:
                    if isinstance(
                        getattr(self.hdfile.root.Detections.cols, name)[0], np.bytes_
                    ):
                        setattr(
                            self,
                            name,
                            getattr(self.hdfile.root.Detections.cols, name)[:].astype(
                                str
                            ),
                        )
                    else:
                        setattr(
                            self,
                            name,
                            getattr(self.hdfile.root.Detections.cols, name)[:],
                        )

            if self.survey == "hdr3":
                if catalog_type == "lines":
                    if verbose:
                        print("Adjusting noise values by 7% where applicable")
                    # adjust noise at IFU edges by factor of 1.07. This will affect the
                    # sn measures and the flux_noise_1sigma values in the lines catalog only
                    sel_fib1 = ((self.amp == "RU") | (self.amp == "LL")) & (
                        self.fibnum <= 12
                    )
                    sel_fib2 = ((self.amp == "LU") | (self.amp == "RL")) & (
                        self.fibnum >= 101
                    )
                    sel_fib = sel_fib1 | sel_fib2

                    self.sn[sel_fib] /= 1.07
                    self.sn_3fib[sel_fib] /= 1.07
                    self.sn_3fib_cen[sel_fib] /= 1.07
                    self.sn_cen[sel_fib] /= 1.07
                    self.flux_noise_1sigma[sel_fib] *= 1.07

                    if verbose:
                        print("Applying HDR3 flux corrections")

                    self.flux /= self.wd_corr_f(self.wave)
                    self.flux_err /= self.wd_corr_f(self.wave)
                    self.continuum /= self.wd_corr_f(self.wave)
                    self.continuum_err /= self.wd_corr_f(self.wave)
                    self.flux_noise_1sigma /= self.wd_corr_f(self.wave)
            elif float(self.survey[3:]) >= 4.0:  # == "hdr3":
                if verbose:
                    print("Applying HDR4 and up WD flux corrections")

                    self.flux /= self.wd_corr_f(self.wave)
                    self.flux_err /= self.wd_corr_f(self.wave)
                    self.continuum /= self.wd_corr_f(self.wave)
                    self.continuum_err /= self.wd_corr_f(self.wave)
                    self.flux_noise_1sigma /= self.wd_corr_f(self.wave)

            elif self.survey == "hdr2.1":
                # Fix fluxes and continuum values for aperture corrections
                wave = self.hdfile.root.Detections.cols.wave[:]
                apcor = self.hdfile.root.Spectra.cols.apcor[:]
                wave_spec = self.hdfile.root.Spectra.cols.wave1d[:]

                apcor_array = np.ones_like(wave)
                for idx in np.arange(0, np.size(wave)):
                    sel_apcor = np.where(wave_spec[idx, :] > wave[idx])[0][0]
                    apcor_array[idx] = apcor[idx, sel_apcor]
                self.apcor = apcor_array

                if catalog_type == "lines":
                    # remove E(B-V)=0.02 screen extinction
                    fix = get_2pt1_extinction_fix()

                    self.flux /= fix(wave)
                    self.flux_err /= fix(wave)
                    self.continuum /= fix(wave)
                    self.continuum_err /= fix(wave)

                    # store observed flux values in new columns
                    self.flux_obs = self.flux.copy()
                    self.flux_err_obs = self.flux_err.copy()
                    self.continuum_obs = self.continuum.copy()
                    self.continuum_err_obs = self.continuum_err.copy()

                    # apply extinction to observed values to get
                    # dust corrected values
                    # Apply S&F 2011 Extinction correction from SFD Map
                    # https://iopscience.iop.org/article/10.1088/0004-637X/737/2/103#apj398709t6

                    self.coords = SkyCoord(
                        self.ra * u.degree, self.dec * u.degree, frame="icrs"
                    )

                    sfd = SFDQuery()
                    self.ebv = sfd(self.coords)
                    Rv = 3.1
                    corr_SF2011 = 2.742  # Landolt V
                    ext = []

                    self.Av = corr_SF2011 * self.ebv

                    for index in np.arange(np.size(self.detectid)):
                        src_wave = np.array([np.double(self.wave[index])])
                        ext_i = extinction.fitzpatrick99(src_wave, self.Av[index], Rv)[
                            0
                        ]
                        ext.append(ext_i)

                    deredden = 10 ** (0.4 * np.array(ext))

                    self.flux *= deredden
                    self.flux_err *= deredden
                    self.continuum *= deredden
                    self.continuum_err *= deredden

            # add in the elixer probabilties and associated info:
            if self.survey == "hdr1" and catalog_type == "lines":
                self.hdfile_elix = tb.open_file(self.config.elixerh5, mode="r")
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
                    pass
                    # print("No Elixer table found")

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
                self.add_hetdex_gmag(loadpickle=True, picklefile=self.config.gmags)

            if self.survey == "hdr1":
                if PYTHON_MAJOR_VERSION < 3:
                    self.plae_poii_hetdex_gmag = np.array(
                        pickle.load(open(self.config.plae_poii_hetdex_gmag, "rb"))
                    )
                else:
                    self.plae_poii_hetdex_gmag = np.array(
                        pickle.load(
                            open(self.config.plae_poii_hetdex_gmag, "rb"),
                            encoding="bytes",
                        )
                    )

        elif searchable:
            if catalog_type == "index":
                # just get coordinates, wavelength and detectid
                self.detectid = self.hdfile.root.DetectIndex.cols.detectid[:]
                self.ra = self.hdfile.root.DetectIndex.cols.ra[:]
                self.dec = self.hdfile.root.DetectIndex.cols.dec[:]
                self.wave = self.hdfile.root.DetectIndex.cols.wave[:]
            else:
                # just get coordinates, wavelength and detectid
                self.detectid = self.hdfile.root.Detections.cols.detectid[:]
                self.ra = self.hdfile.root.Detections.cols.ra[:]
                self.dec = self.hdfile.root.Detections.cols.dec[:]
                self.wave = self.hdfile.root.Detections.cols.wave[:]

        # set the SkyCoords
        if searchable or loadtable:
            self.coords = SkyCoord(
                self.ra * u.degree, self.dec * u.degree, frame="icrs"
            )

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

    def refine(self, gmagcut=None, remove_large_gal=True, d25scale=1.5):
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

    def query_by_coord(self, coord, radius=None, astropy=True):
        """
        Returns list of detecitds or astropy table object
        based on a coordinate search

        Paramters
        ---------

        self
            Detections Class object. If catalog_type=='index', healpix
            searching is done and is most efficient. Can also be
            applied to other catalog_types with the searchable==True
            or loadtable==True option set
        coord - astropy coordinate object
            only works for a single coordinate objects
        radius
           an astropy Quantity object
        astropy bool
           if True will return an astropy table of detection info
           if False will return detectid list

        Returns
        -------
        Detectid or Table
        """
        if radius is None:
            radius = 3.0 * u.arcsec

        if self.catalog_type in ["lines", "broad", "continuum"]:
            if self.loadtable:
                sep = self.coords.separation(coord)

                maskcoords = sep < radius

                if astropy:
                    return self[maskcoords].return_astropy_table()
                else:
                    return self.detectid[maskcoords]
            else:
                print(
                    "You must use loadtable==True for catolog_type=={}".format(
                        self.catalog_type
                    )
                )

        elif self.catalog_type == "index":
            Nside = 2**15

            ra_obj = coord.ra.deg
            dec_obj = coord.dec.deg

            ra_sep = radius.to(u.degree).value + 6.0 / 3600.0

            vec = hp.ang2vec(ra_obj, dec_obj, lonlat=True)

            pix_region = hp.query_disc(Nside, vec, (ra_sep * np.pi / 180))

            seltab = Table()

            for hpix in pix_region:
                tab = Table(self.hdfile.root.DetectIndex.read_where("healpix == hpix"))
                seltab = vstack([seltab, tab])

            # query sub-table based on coordinates

            tab_coords = SkyCoord(ra=seltab["ra"] * u.deg, dec=seltab["dec"] * u.deg)
            sep = tab_coords.separation(coord)
            seltab["separation"] = sep.arcsec
            maskcoords = sep < radius

            tab = seltab[maskcoords]
            tab.sort("detectid")

            if len(tab) > 0:
                if astropy:
                    return tab
                else:
                    return list(tab["detectid"])
            else:
                return None

    def find_match(self, coord, radius=None, wave=None, dwave=5.0, shotid=None):
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

        if radius is None:
            radius = 5.0 * u.arcsec
        selmatch = self.query_by_coord(coord, radius)

        if wave is not None:
            selwave = np.abs(self.wave - wave) < dwave
            selmatch = selwave * selmatch

        if shotid is not None:
            selshot = self.shotid == shotid
            selmatch = selshot * selmatch

        return selmatch

    def query_by_dictionary(self, limits):
        """
        Note: This function is pretty obsolete. You should
        use the source catalog

        Takes a dictionary of query limits
        and reduces the detections database. This
        can either be taken from the dictionary saved
        by detwidgets.py GUI or can take the following
        form:

        # limits to be returned to query detections

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
            maskfield = self.query_by_coord(coords, limits.rad)
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

        # set an empty mask to start
        mask = np.zeros(np.size(self.detectid), dtype=bool)

        baddetects = np.loadtxt(self.config.baddetect, dtype=int)

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

        # set an empty mask to start

        mask = np.zeros(np.size(self.detectid), dtype=bool)

        if self.survey == "hdr1":
            badamps1 = ascii.read(
                self.config.badamp, names=["ifuslot", "amp", "date_start", "date_end"]
            )

            badamps2 = ascii.read(
                self.config.badamp, names=["ifuslot", "amp", "date_start", "date_end"]
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
            badamps = Table.read(self.config.badamp)

            det_table = self.return_astropy_table()

            join_tab = join(
                det_table, badamps, keys=["shotid", "multiframe"], join_type="left"
            )
            # this is needed to match with detection class object sorting
            join_tab.sort("detectid")

            mask1 = join_tab["flag"] != 0

            del det_table, join_tab

            if self.recentbadamp:  # This is really slow on the full catalog.
                print("Adding in newly found badamps")
                # add in any newly found badamps that haven't made it into the
                # amp_flag.fits file yet

                mask2 = np.zeros(np.size(self.detectid), dtype=bool)

                badamps2 = Table.read(self.config.badamp2, format="ascii")

                for row in badamps2:
                    selmf = self.multiframe == row["multiframe"]
                    seldate = (self.date >= row["date_start"]) * (
                        self.date <= row["date_end"]
                    )
                    mask2 = np.logical_or(mask2, selmf * seldate)

                badamps_single = Table.read(config.badamp_single, format="ascii")

                mask3 = np.zeros(np.size(self.detectid), dtype=bool)

                for row in badamps_single:
                    selmf = self.multiframe == row["multiframe"]
                    selshotid = self.shotid == row["shotid"]

                    mask3 = np.logical_or(mask3, selmf * selshotid)

                mask = mask1 * np.logical_not(mask2) * np.logical_not(mask3)

            else:
                mask = mask1

            return mask

    def remove_bad_pix(self):
        """
        Refines catalog based on badpixel lists found for
        the release in hetdex-api:
        https://github.com/HETDEX/hetdex_api/tree/master/known_issues/hdr4/badpix.list
        """

        if True:
            badpixlist = ascii.read(
                self.config.badpix, names=["multiframe", "x1", "x2", "y1", "y2"]
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

        mask = np.zeros(np.size(self.detectid), dtype=bool)
        badshots = np.loadtxt(self.config.badshot, dtype=int)

        for shot in badshots:
            maskshot = self.shotid == shot
            mask = np.logical_or(maskshot, mask)

        self.vis_class[mask] = -2

        if self.survey == "hdr4":
            # remove shots that were incorrectly included
            hdr4_shots_to_remove = [
                20190802013,
                20190802026,
                20190802027,
                20190802028,
                20190802029,
            ]
            for shot in hdr4_shots_to_remove:
                maskshot = self.shotid == shot
                mask = np.logical_or(maskshot, mask)

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

        met_tab = Table.read(self.config.meteor, format="ascii")

        mask = np.ones_like(self.detectid, dtype=bool)

        for row in met_tab:
            sel_shot = np.where(self.shotid == row["shotid"])[0]
            for idx in sel_shot:
                coord = self.coords[idx]
                mask[idx] = meteor_flag_from_coords(coord, row["shotid"])
        return mask

    def remove_large_gal(self, d25scale=1.5):
        """
        Returns boolean mask with detections landing within
        galaxy defined by d25scale flagged as False.

        Based on check_all_large_gal from hetdex_tools/galmask.py
        written by John Feldmeier
        """

        galaxy_cat = Table.read(self.config.rc3cat, format="ascii")

        mask = np.ones_like(self.detectid, dtype=bool)

        # Loop over each galaxy

        for idx, row in enumerate(galaxy_cat):
            gal_coord = SkyCoord(row["Coords"], frame="icrs")
            rlimit = 1.1 * d25scale * row["SemiMajorAxis"] * u.arcmin

            if np.any(self.coords.separation(gal_coord) < rlimit):
                galregion = create_gal_ellipse(
                    galaxy_cat, row_index=idx, d25scale=d25scale
                )
                dummy_wcs = create_dummy_wcs(
                    galregion.center, imsize=2 * galregion.height
                )
                galflag = galregion.contains(self.coords, dummy_wcs)
                mask = mask * np.invert(galflag)

        return mask

    def get_detection_info(self, detectid_i, rawh5=False, verbose=False):
        """
        Returns Detections table information from H5 file
        Applies relevent corrections such as noise model fix
        and spectral adjustment if rawh5 is False

        Parameters
        ----------
        detectid: int
            detectid (integer ID) of the detection you want to
            get information for
        rawh5: bool
            if True, this will simply return the row from the detecth5
            file. If False (the default), any relevent correcctions
            are applied.
        verbose: bool
            provide info statements if set to True. Default is False

        Returns
        -------
        det_info: ndarray object
            returns detection info
            column info given by det_info.dtype attribute
        """

        # notice ... this is unhandled so caller gets the exception
        if self.catalog_type == "index":
            det_row = self.hdfile.root.DetectIndex.read_where("detectid == detectid_i")
            return det_row
        else:
            det_row = self.hdfile.root.Detections.read_where("detectid == detectid_i")

        if rawh5:
            if verbose:
                print("Returning raw H5 Detections table row")
        else:
            if verbose:
                print("Returning updated Detections table row")
            if (self.catalog_type == "lines") & (self.survey == "hdr3"):
                if verbose:
                    print("Adjusting noise values by 7% where applicable")
                # adjust noise at IFU edges by factor of 1.07. This will affect the
                # sn measures and the flux_noise_1sigma values in the lines catalog only
                sel_fib1 = ((det_row["amp"] == b"RU") | (det_row["amp"] == b"LL")) & (
                    det_row["fibnum"] <= 12
                )
                sel_fib2 = ((det_row["amp"] == b"LU") | (det_row["amp"] == b"RL")) & (
                    det_row["fibnum"] >= 101
                )

                if sel_fib1 | sel_fib2:
                    if verbose:
                        print("Noise model update is required. Update performed.")
                    det_row["sn"] /= 1.07
                    det_row["sn_3fib"] /= 1.07
                    det_row["sn_3fib_cen"] /= 1.07
                    det_row["sn_cen"] /= 1.07
                    det_row["flux_noise_1sigma"] *= 1.07

                if verbose:
                    print("Applying HDR3 flux corrections")
                det_row["flux"] /= self.wd_corr_f(det_row["wave"])
                det_row["flux_err"] /= self.wd_corr_f(det_row["wave"])
                det_row["continuum"] /= self.wd_corr_f(det_row["wave"])
                det_row["continuum_err"] /= self.wd_corr_f(det_row["wave"])
                det_row["flux_noise_1sigma"] /= self.wd_corr_f(det_row["wave"])
            elif (self.catalog_type == "lines") and float(self.survey[3:]) >= 4.0:
                # specifically HDR3 handled just above, but only the lines catalog
                # since these are not the spectra, just the emission line
                if verbose:
                    print("Applying HDR4 and up WD flux corrections")
                det_row["flux"] /= self.wd_corr_f(det_row["wave"])
                det_row["flux_err"] /= self.wd_corr_f(det_row["wave"])
                det_row["continuum"] /= self.wd_corr_f(det_row["wave"])
                det_row["continuum_err"] /= self.wd_corr_f(det_row["wave"])
                det_row["flux_noise_1sigma"] /= self.wd_corr_f(det_row["wave"])

        return det_row

    def get_survey_info(self, detectid_i):
        """
        Return Survey class info for detectid

        Parameters
        ----------
        detectid: int
            detectid (integer ID) of the detection you want to
            grab the spectrum for

        Returns
        -------
        survey_info: ndarray object
            row of information from the survey_hdrX.h5 table
            column info given by survey_info.dtype attribute
        """

        det_row = self.get_detection_info(detectid_i)

        shotid_i = det_row["shotid"]

        survey_row = self.surveyh5.root.Survey.read_where("shotid == shotid_i")

        return survey_row

    def get_coord(self, detectid_i):
        """
        Return SkyCoord object for a single detection

        Parameters
        ----------
        detectid: int
            detectid (integer ID) of the detection you want to
            get SkyCoord object for

        Return
        ------
        coords: SkyCoord Class object

        """

        det_row = self.get_detection_info(detectid_i)

        coord = SkyCoord(ra=det_row["ra"][0] * u.deg, dec=det_row["dec"][0] * u.deg)

        return coord

    def get_fiber_info(self, detectid_i):
        """
        Return Fibers table for a single detectid

        Parameters
        ----------
        detectid: int
            detectid (integer ID) of the detection you want to
            get the fiber information

        Returns
        -------
        fiber_info: ndarray object
            table of fiber information for the specific
            detection. Use fiber_info.dtype to get column info

        """

        fib_info = self.hdfile.root.Fibers.read_where("detectid == detectid_i")

        return fib_info

    def get_spectrum(
        self,
        detectid_i,
        deredden=False,
        apply_extinction_fix=True,
        add_apcor=False,
        rawh5=False,
        verbose=False,
        ffsky=False,
    ):
        """
        Grabs the 1D spectrum used to measure fitted parameters.

        Parameters
        ----------
        detectid: int
            detectid (integer ID) of the detection you want to
            grab the spectrum for
        deredden: bool
            flag to return dust corrected spectrum
        apply_extinction_fix: bool
            HDR2 correction. Does not apply to any other release
        add_apcor: bool
            add a column with the applied aperture correction
        rawh5: bool
            if False, this will convert to 1AA binning and apply
            any relevent spectral correcctions. This is the default
            if True, it will only pull the spetral data from the detections
            h5 file
        verbose: bool
            provide info statements if set to True. Default is False
        ffsky: bool
            option to use full frame sky calibration. Default is False

        Returns
        -------
        spec_table: astropy table
            2 or 3 column astropy table with spectral data

        """
        spectra = self.hdfile.root.Spectra
        spectra_table = spectra.read_where("detectid == detectid_i")

        data = Table()

        if rawh5:
            intensityunit = u.erg / (u.cm**2 * u.s * 2 * u.AA)
        else:
            intensityunit = u.erg / (u.cm**2 * u.s * u.AA)

        data["wave1d"] = Column(spectra_table["wave1d"][0], unit=u.AA)

        if ffsky:
            data["spec1d"] = Column(
                spectra_table["spec1d_ffsky"][0], unit=1.0e-17 * intensityunit
            )
        else:
            data["spec1d"] = Column(
                spectra_table["spec1d"][0], unit=1.0e-17 * intensityunit
            )
        data["spec1d_err"] = Column(
            spectra_table["spec1d_err"][0], unit=1.0e-17 * intensityunit
        )

        if add_apcor:
            data["apcor"] = Column(spectra_table["apcor"][0])

        if rawh5 is False:
            if verbose:
                print("Converting from 2AA to 1AA binning")
            # convert from 2AA binning to 1AA binning:
            data["spec1d"] /= 2.0
            data["spec1d_err"] /= 2.0
        else:
            if verbose:
                print("Units are in 2AA binning.")

        if self.survey == "hdr2.1" and apply_extinction_fix:
            # remove E(B-V)=0.02 screen extinction
            fix = get_2pt1_extinction_fix()

            flux_corr = fix(data["wave1d"])
            data["spec1d"] /= flux_corr
            data["spec1d_err"] /= flux_corr

        elif self.survey == "hdr3":
            if rawh5 is False:
                if verbose:
                    print("Applying spectral correction")
                try:
                    data["spec1d"] /= self.wd_corr["corr"]
                    data["spec1d_err"] /= self.wd_corr["corr"]
                except:
                    print("Warning! Could not apply wd_corr.")



                # Apply HDR3 noise model correction
                if verbose:
                    print("Applying HDR3 noise model update")
                det_row = self.get_detection_info(detectid_i)

                # adjust noise at IFU edges by factor of 1.07.

                sel_fib1 = ((det_row["amp"] == "RU") | (det_row["amp"] == "LL")) & (
                    det_row["fibnum"] <= 12
                )
                sel_fib2 = ((det_row["amp"] == "LU") | (det_row["amp"] == "RL")) & (
                    det_row["fibnum"] >= 101
                )

                if sel_fib1 | sel_fib2:
                    if verbose:
                        print("Noise model is required. Update performed.")
                    data["spec1d_err"] *= 1.07
                else:
                    if verbose:
                        print("Noise model adjustment not required")

        elif (
            float(self.survey[3:]) >= 4.0 or (self.survey[0:3]).lower() == "pdr"
        ):  # == "hdr4":
            if rawh5 is False:
                if verbose:
                    print("Applying spectral correction")
                try:
                    data["spec1d"] /= self.wd_corr["corr"]
                    data["spec1d_err"] /= self.wd_corr["corr"]
                except:
                    print("Warning! Could not apply wd_corr.")

        if deredden:
            if verbose:
                print("Applying dust correction to spectrum")

            det_row = self.hdfile.root.Detections.read_where("detectid == detectid_i")
            coords = SkyCoord(ra=det_row["ra"][0], dec=det_row["dec"][0], unit="deg")
            deredden_corr = deredden_spectra(data["wave1d"], coords)
            data["spec1d"] *= deredden_corr
            data["spec1d_err"] *= deredden_corr

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
        try:
            if self.survey == "hdr2.1":
                table.add_column(Column(self.apcor), name="apcor")
        except AttributeError:
            pass

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
                pass
            try:
                table.add_column(
                    Column(self.fluxlimit_4540), index=3, name="fluxlimit_4540"
                )
            except:
                print("Could not add average flux limit")

        try:
            table.add_column(self.fiber_ratio, name="fiber_ratio")
        except:
            pass

        # add new elixer info and fiber_ratio column
        try:
            table.add_column(self.best_z, name="best_z")
            table.add_column(self.best_pz, name="best_pz")
            table.add_column(self.flags_elixer, name="flags_elixer")
            table.add_column(self.multiline_name, name="multiline_name")
            table.add_column(self.classification_labels, name="classification_labels")
            table.add_column(self.counterpart_mag, name="counterpart_mag")
            table.add_column(self.counterpart_mag_err, name="counterpart_mag_err")
            table.add_column(self.counterpart_dist, name="counterpart_dist")
            table.add_column(
                self.counterpart_catalog_name, name="counterpart_catalog_name"
            )
            table.add_column(
                self.counterpart_filter_name, name="counterpart_filter_name"
            )
            table.add_column(self.forced_mag, name="forced_mag")
            table.add_column(self.forced_mag_err, name="forced_mag_err")
            table.add_column(self.forced_catalog_name, name="forced_catalog_name")
            table.add_column(self.forced_filter_name, name="forced_filter_name")
            table.add_column(self.forced_radius, name="forced_radius")

        except:
            pass

        if self.survey == "hdr2.1":
            try:
                table.add_column(self.Av, name="Av")
                table.add_column(self.ebv, name="ebv")
                table.add_column(self.flux_obs, name="flux_obs")
                table.add_column(self.flux_err_obs, name="flux_obs_err")
                table.add_column(self.continuum_obs, name="continuum_obs")
                table.add_column(self.continuum_err_obs, name="continuum_obs_err")
            except Exception:
                pass

        return table

    def save_spectrum(self, detectid_i, outfile=None):
        spec_data = self.get_spectrum(detectid_i)
        if outfile:
            ascii.write(spec_data, outfile, overwrite=True)
        else:
            ascii.write(spec_data, "spec_" + str(detectid_i) + ".dat", overwrite=True)

    def plot_spectrum(self, detectid_i, xlim=None, ylim=None):
        print("We have removed this code bit from the API for efficiency")

    #        spec_data = self.get_spectrum(detectid_i)
    #        plt.figure(figsize=(8, 6))
    #        plt.errorbar(
    #            spec_data["wave1d"], spec_data["spec1d"], yerr=spec_data["spec1d_err"]
    #        )
    #        plt.title("DetectID " + str(detectid_i))
    #        plt.xlabel("wave (AA)")
    #        plt.ylabel("flux (1e-17 erg/s/cm^2/AA)")
    #        if xlim is not None:
    #            plt.xlim(xlim)
    #            if ylim is not None:
    #                plt.ylim(ylim)
    #        plt.show()

    def __len__(self):
        try:
            # if instatiated w/loadtable=False, there is no length (no records are loaded, self.ra does not exist)
            # NOTICE: other calls will also fail, like remove_bad_amps(), etc
            return len(self.ra)
        except:
            # could alternately use None or np.nan, but 0 may be the preferred use case.
            # treat as empty rather than an undefined length
            return 0

    def close(self):
        self.hdfile.close()

    def get_detection_flags(self, detectid, F=None, FI=None):
        """
        Parameters
        ----------
        detectid int

        F: Fibers class object
        Optional to pass through the Fibers class so that it does not have to be reopen for each detectid
        FI : FiberIndex class object
        Optional to use for flagging 

        Returns
        -------
        flag diction
        keys are flag_pixmask, flag_badamp, flag_badpix, flag_badfib, flag_meteor, flag_largegal, flag_chi2fib
        1 if good, 0 if flagged bad
        """

        # close Fiber class object at end if not passed in
        closeF = False

        if F is None:
            closeF = True

        flag_dict = {
            "flag_badamp": 1,
            "flag_badpix": 1,
            "flag_badfib": 1,
            "flag_pixmask": 1,
            "flag_meteor": 1,
            "flag_largegal": 1,
            "flag_chi2fib": 1,
            "flag_satellite": 1,
            "flag_badcal": 1,
            'flag_throughput': 1,
            'flag_shot': 1,
        }

        det_info = self.get_detection_info(detectid)[0]
        det_fib_info = self.get_fiber_info(detectid)
        det_coords = self.get_coord(detectid)
        shotid = det_info["shotid"]
        expnum = det_info["expnum"]

        # find detection wavelength index
        wave = det_info["wave"]
        # wave_i = np.argmin(np.abs(wave_rect - wave))

        if FI is None:
            FI = FiberIndex(self.survey)

        flag_dict_FI = FI.get_fiber_flags(det_coords, shotid=shotid)

        for key in ['flag_meteor', 'flag_satellite', 'flag_largegal','flag_throughput', 'flag_shot']:
            flag_dict[key] = flag_dict_FI[key]


        if F is None:
            F = Fibers(shotid, survey=self.survey)

        # grab fibers table to get calfibe info
        fibers_table = get_fibers_table(shotid, coords=det_coords, F=F)

        ifiber = det_fib_info["weight"].argsort()[::-1]

        # store list of 3 highest weight fibers and their fiber info
        fiberids = det_fib_info[ifiber]["fiber_id"][0:3].astype(str)
        mfs = det_fib_info[ifiber]["multiframe"][0:3].astype(str)
        xs = det_fib_info[ifiber]["x_raw"][0:3]
        ys = det_fib_info[ifiber]["y_raw"][0:3]
        fibnums = det_fib_info[ifiber]["fibnum"][0:3]

        # check for flag_badamp
        if self.badamps is None:
            self.badamps = Table.read(self.config.badamp)

        for mf in np.unique(mfs):
            sel_row = (self.badamps["multiframe"] == mf) & (
                self.badamps["shotid"] == shotid
            )

            if np.any(self.badamps["flag"][sel_row] == 0):
                flag_dict["flag_badamp"] = 0
            elif np.sum(sel_row) == 0:
                flag_dict["flag_badamp"] = 0

        # check for flag_badpix and flag_badfib
        if self.badpix is None:
            self.badpix = Table.read(
                self.config.badpix,
                format="ascii",
                names=["multiframe", "x1", "x2", "y1", "y2"],
            )

        if self.badfib is None:
            self.badfib = Table.read(self.config.badfib, format="ascii")

        if self.badfib_transient is None:
            self.badfib_transient = Table.read(
                self.config.badfib_transient,
                format="ascii",
                names=["multiframe", "fibnum", "date_start", "date_end"],
            )

        for i in np.arange(0, 2):
            # only check for bad pix for line sources
            if self.catalog_type == "lines":
                sel_mf = self.badpix["multiframe"] == mfs[i]
                sel_x = (xs[i] >= self.badpix["x1"]) & (xs[i] <= self.badpix["x2"])
                sel_y = (ys[i] >= self.badpix["y1"]) & (ys[i] <= self.badpix["y2"])
                sel_badpix = sel_mf & sel_x & sel_y

                if np.sum(sel_badpix) > 0:
                    flag_dict["flag_badpix"] = 0

            sel_fib = (self.badfib["multiframe"] == mfs[i]) & (
                self.badfib["fibnum"] == fibnums[i]
            )

            if np.sum(sel_fib) > 0:
                flag_dict["flag_badfib"] = 0

            sel_fib2 = (self.badfib_transient["multiframe"] == mfs[i]) & (
                self.badfib_transient["fibnum"] == fibnums[i]
            )
            sel_date = (det_info["date"] >= self.badfib_transient["date_start"]) * (
                det_info["date"] <= self.badfib_transient["date_end"]
            )

            if np.sum(sel_fib2 * sel_date) > 0:
                flag_dict["flag_badfib"] = 0

        # flag detection if any 3 spectral resolution elements are masked
        if self.catalog_type == "lines":  # not needed for continuum sources
            chi2fib = np.ndarray((3, 3))

            for i, fib_i in enumerate(fiberids):
                sel_fib = np.where(fibers_table["fiber_id"] == fib_i)[0]

                if len(sel_fib) == 0:
                    continue

                # access native wavelength and spectrum
                wavelength = fibers_table["wavelength"][sel_fib][0]
                spectrum = fibers_table["spectrum"][sel_fib][0]
                chi2 = fibers_table["chi2"][sel_fib][0]
                # error1D = fibers_table['error1D'][sel_fib][0]

                wave_io = np.argmin(np.abs(wavelength - wave))

                if np.any(spectrum[wave_io - 1 : wave_io + 2] == 0.0):
                    flag_dict["flag_pixmask"] = 0

                selwave = np.abs(wavelength - wave) < 5
                if np.sum(selwave) > 0:
                    chi2fib[i] = np.max(chi2[np.where(selwave)])

            flag_dict["chi2fib"] = chi2fib

            if np.max(chi2fib) >= 4.5:
                flag_dict["flag_chi2fib"] = 0

        # add calflagging
        if self.catalog_type == "lines":
            if self.cal5200_tab is None:
                # open table if not yet opened
                self.cal5200_tab = Table.read(
                    self.config.cal5200,
                    format="ascii",
                    names=["shotid", "multiframe", "expnum"],
                )

            if self.cal5460_tab is None:
                # open table if not yet opened
                self.cal5460_tab = Table.read(
                    self.config.cal5460,
                    format="ascii",
                    names=["shotid", "multiframe", "expnum"],
                )

            for mf in np.unique(mfs):
                flag_cal = cal_flag_for_amp_wave(
                    wave,
                    mf,
                    shotid,
                    cal5460_tab=self.cal5460_tab,
                    cal5200_tab=self.cal5200_tab,
                )

                if flag_cal is False:
                    flag_dict["flag_badcal"] = 0

        if closeF:
            F.close()

        return flag_dict
