# -*- coding: utf-8 -*-
"""

Initiates survey class and provides an API to query the survey class based
on astropy coordinates

Added 2020/03/29 = FiberIndex class to query all fibers in the survey

Created on Tue Jan 22 11:02:53 2019
@author: gregz/Erin Mentuch Cooper
"""
from __future__ import print_function

import numpy as np
import tables as tb
import numpy
import copy
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack

import healpy as hp
from hetdex_api.config import HDRconfig

try:
    LATEST_HDR_NAME = HDRconfig.LATEST_HDR_NAME
except Exception as e:
    print("Warning! Cannot find or import HDRconfig from hetdex_api!!", e)
    LATEST_HDR_NAME = "hdr2.1"


class Survey:
    def __init__(self, survey=LATEST_HDR_NAME):
        """
        Initialize the Survey class for a given data release

        Parameters
        ----------
        survey : string
            Data release you would like to load, i.e., 'hdr1','HDR2'
            This is case insensitive.

        Returns
        -------
        Survey class object with the following attributes
        """

        global config
        config = HDRconfig(survey=survey.lower())

        self.filename = config.surveyh5
        self.hdfile = tb.open_file(self.filename, mode="r")
        colnames = self.hdfile.root.Survey.colnames

        for name in colnames:
            if name == "ra_flag":
                setattr(
                    self,
                    name,
                    getattr(self.hdfile.root.Survey.cols, name)[:].astype(str),
                )
            elif isinstance(getattr(self.hdfile.root.Survey.cols, name)[0], np.bytes_):
                setattr(
                    self,
                    name,
                    getattr(self.hdfile.root.Survey.cols, name)[:].astype(str),
                )
            else:
                setattr(self, name, getattr(self.hdfile.root.Survey.cols, name)[:])

        # set the SkyCoords
        self.coords = SkyCoord(self.ra * u.degree, self.dec * u.degree, frame="icrs")

        # append flux limits
        if survey == "hdr2.1":

            flim = Table.read(
                config.flim_avg,
                format="ascii",
                names=["datevobs", "col2", "fluxlimit_4540"],
            )
            fluxlimit = []
            for datevobs in self.datevobs:
                sel = flim["datevobs"] == datevobs
                if np.sum(sel) == 1:
                    fluxlimit.extend(flim["fluxlimit_4540"][sel])
                elif np.sum(sel) > 1:
                    print("Found two fluxlimits for ", datevobs)
                    fluxlimit.extend(flim["fluxlimit_4540"][sel][0])
                else:
                    fluxlimit.append(np.nan)

            self.fluxlimit_4540 = np.array(fluxlimit)

    def __getitem__(self, indx):
        """ 
        This allows for slicing of the survey class
        object so that a mask can be applied to
        every attribute automatically by:

        Example
        -------
        survey_sliced = survey[indx]
        """

        p = copy.copy(self)
        attrnames = self.__dict__.keys()
        for attrname in attrnames:
            try:
                setattr(p, attrname, getattr(self, attrname)[indx])
            except:
                setattr(p, attrname, getattr(self, attrname))
        return p

    def slice(self):
        """
        This will slice the survey class upon
        initialization to remove
        any bad shots from the survey so that
        only one class variable is used.

        Example
        -------
        survey = Survey('hdr1').slice()

        """

        maskshots = self.remove_shots()
        return self[maskshots]

    def remove_shots(self):
        """
        Function to remove shots that should be exluded from most analysis.

        Parameters
        ----------
        self
           Survey class object

        Returns
        -------
        ind_good
           Boolean mask of good shots

        Examples
        --------
        S = Survey('hdr2')
        ind_good_shots = S.remove_shots()
        S_good = S[ind_good_shots]
        
        """
        global config

        mask = np.zeros(np.size(self.shotid), dtype=bool)
        badshots = np.loadtxt(config.badshot, dtype=int)
        for shot in badshots:
            maskshot = self.shotid == shot
            mask = mask | maskshot

        notvalid = self.shotid < 20170000
        mask = mask | notvalid

        return np.invert(mask)

    def get_shotlist(self, coords, radius=None, width=None, height=None):
        """
        Returns a list of shots for defined region of sky.

        Parameters
        ----------
        self
            Survey Class object
        coords
            astropy coordinate object
        radius
            an astropy Quantity object, or a string that can be parsed into
            one.  e.g., '1 degree' or 1*u.degree.  If radius is specified,
            the shape is assumed to be a circle
        width
            in degress.  Specifies the edge length of a square box. If a radius
            is not specified a width and height of a region must be given.
        height
            in degrees.  Specifies the height of a rectangular box.  Must be
            passed with width.

        Examples
        --------
        S = Survey('hdr1')
        coords = SkyCoord(150.025513 * u.deg, 2.087767 * u.deg, frame='icrs')
        shotlist = S.get_shotlist(coords, radius=0.5*u.deg)

        or for a rectangular region around the point defined in coords

        shotlist = S.get_shotlist(coords, width=0.5, height=0.1)
        """

        try:
            self.coords
        except:
            self.coords = SkyCoord(
                self.ra * u.degree, self.dec * u.degree, frame="icrs"
            )

        if radius is not None:
            try:
                idx = self.coords.separation(coords) < radius
            except:
                print("Assuming radius in degrees")
                idx = self.coords.separation(coords) < radius * u.degree
        else:
            try:
                idx1 = abs(self.ra - coords.ra.value) < width / 2.0
                idx2 = abs(self.dec - coords.dec.value) < height / 2.0
                idx = idx1 * idx2
            except:
                print("Provide both width and height of sky region in degrees.")

        return self.shotid[idx]

    def return_astropy_table(self, return_good=True):
        """
        Function to return an astropy table that is machine readable

        Parameters
        ----------
        self
            Survey Class object 
        return_good
            Boolean flag to only exclude shots in config.badshot

        Returns
        -------
        survey_table
            astropy table of the Survey class object
        """

        survey_table = Table(self.hdfile.root.Survey.read())

        good_shots = self.remove_shots()

        survey_table["shot_flag"] = good_shots

        survey_table["mjd"] = self.mjd[:, 0]
        survey_table["exptime"] = np.mean(self.exptime, axis=1)
        survey_table["fluxlimit_4540"] = self.fluxlimit_4540

        for col in survey_table.colnames:
            try:
                if np.shape(survey_table[col])[1] == 3:
                    survey_table.remove_column(col)
            except:
                pass

        if return_good:
            return survey_table[good_shots]
        else:
            return survey_table

    def close(self):
        """
        Close the HDF5 file when you are done using
        it to release anything that might be in memory

        Example
        -------

        S.close()
        """

        self.hdfile.close()


class FiberIndex:
    def __init__(self, survey=LATEST_HDR_NAME, loadall=False):
        """
        Initialize the Fiber class for a given data release
        
        Parameters
        ----------
        survey : string
            Data release you would like to load, i.e., 'hdr1','HDR2'
            This is case insensitive.

        Returns
        -------
        FiberIndex class object
        """
        self.survey = survey

        if self.survey == "hdr1":
            print("Sorry there is no FiberIndex for hdr1")
            return None

        global config
        config = HDRconfig(survey=survey.lower())

        self.filename = config.fiberindexh5
        self.hdfile = tb.open_file(self.filename, mode="r")

        if loadall:
            colnames = self.hdfile.root.FiberIndex.colnames

            for name in colnames:

                if isinstance(
                    getattr(self.hdfile.root.FiberIndex.cols, name)[0], np.bytes_
                ):
                    setattr(
                        self,
                        name,
                        getattr(self.hdfile.root.FiberIndex.cols, name)[:].astype(str),
                    )
                else:
                    setattr(
                        self, name, getattr(self.hdfile.root.FiberIndex.cols, name)[:]
                    )

            self.coords = SkyCoord(
                self.ra[:] * u.degree, self.dec[:] * u.degree, frame="icrs"
            )

    def query_region(self, coords, radius=3.0 * u.arcsec, shotid=None):
        """
        Function to retrieve the indexes of the FiberIndex table
        for a specific region

        Parameters
        ----------
        self
            the FiberIndex class for a specific survey
        coords
            center coordinate you want to search. This should
            be an astropy SkyCoord object
        radius
            radius you want to search. An astropy quantity object
        shotid
            Specific shotid (dtype=int) you want

        Returns
        -------
        An astropy table of Fiber infomation in queried aperture
        """

        Nside = 2 ** 15

        ra_obj = coords.ra.deg
        dec_obj = coords.dec.deg

        ra_sep = radius.to(u.degree).value + 3.0 / 3600.0

        vec = hp.ang2vec(ra_obj, dec_obj, lonlat=True)

        pix_region = hp.query_disc(Nside, vec, (ra_sep * np.pi / 180))

        seltab = Table()
        for hpix in pix_region:
            if shotid:
                h_tab = self.get_fib_from_hp(hpix, shotid=shotid)
            else:
                h_tab = self.get_fib_from_hp(hpix)
            seltab = vstack([seltab, h_tab])

        fibcoords = SkyCoord(
            seltab["ra"] * u.degree, seltab["dec"] * u.degree, frame="icrs"
        )

        idx = coords.separation(fibcoords) < radius
        
        return seltab[idx]

        
    def get_fib_from_hp(self, hp, shotid=None, astropy=True):

        if astropy:

            tab= Table(self.hdfile.root.FiberIndex.read_where("healpix == hp"))
            
            if shotid:
                sel_shot = tab['shotid'] == shotid
                return tab[sel_shot]
            else:
                return tab
        else:
            if shotid:
                sid = shotid
                return self.hdfile.root.FiberIndex.read_where("(healpix == hp) & (shoti d== sid)")
            else:
                return self.hdfile.root.FiberIndex.read_where("(healpix == hp)")
                
    def get_closest_fiberid(self, coords, shotid=None, maxdistance=8.*u.arcsec):
        """
        Function to retrieve the closest fiberid in a shot
        
        Parameters
        ----------
        self
            the FiberIndex class for a specific survey
        coords
            center coordinate you want to search. This should
            be an astropy SkyCoord object
        shotid
            Specific shotid (dtype=int) you want
        maxdistance
            The max distance you want to search for a nearby fiber.
            Default is 8.*u.arcsec 
        
        Returns
        -------
        fiberid
            unique fiber identifier for the shot
        
        """

        #start searching at small radius to search more efficiently
        search = 2.*u.arcsec

        while search <= maxdistance:
            fiber_table = self.query_region(coords, radius=search, shotid=shotid)
            search = search + 4.*u.arcsec
            if np.size(fiber_table) > 0:
                break

        if np.size(fiber_table) > 0:
            fibcoords = SkyCoord(
                fiber_table["ra"] * u.degree, fiber_table["dec"] * u.degree, frame="icrs"
            )
            
            idx = np.argmin(coords.separation(fibcoords))
            fiberid = fiber_table['fiber_id'][idx]
        
            return fiber_table['fiber_id'][idx]
        else:
            return None

    def close(self):
        """
        Close the hdfile when done
        """
        self.hdfile.close()
