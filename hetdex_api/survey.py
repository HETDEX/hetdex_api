# -*- coding: utf-8 -*-
"""

Initiates survey class and provides an API to query the survey class based
on astropy coordinates

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

from hetdex_api.config import HDRconfig


class Survey:
    def __init__(self, survey='hdr1'):
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
        function to remove shots that should be exluded from most analysis.
        They are engineering shots that were accidentally included in the
        release.
        """
        global config

        mask = np.zeros(np.size(self.shotid), dtype=bool)
        badshots = np.loadtxt(config.badshot, dtype=int)
        for shot in badshots:
            maskshot = self.shotid == shot
            mask = mask | maskshot

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

        if radius:
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

    def close(self):
        """
        Close the HDF5 file when you are done using
        it to release anything that might be in memory

        Example
        -------

        S.close()
        """

        self.hdfile.close()
