# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 11:02:53 2019

@author: gregz/Erin Mentuch Cooper
"""

import tables as tb
import numpy
import astropy.units as u
from astropy.coordinates import SkyCoord

import config

class Survey:
    def __init__(self, survey):
        '''
        Initialize the Survey class for a given data release

        Input
        -----
        survey : string
            Data release you would like to load, i.e., 'DR1' or 'Parallel'.
            This is case insensitive.
        '''
        survey_options = {'hdr1': config.surveyh5,
                          'parallel': 'PATHNAME'}
        if survey.lower() not in survey_options:
            print('survey not in survey options')
            print(survey_options)
            return None
        self.filename = survey_options[survey]
        self.hdfile = tb.open_file(self.filename, mode='r')
        colnames = self.hdfile.root.Survey.colnames
        for name in colnames:
            setattr(self, name,
                    getattr(self.hdfile.root.Survey.cols, name)[:])
            
        # set the SkyCoords
        self.coords = SkyCoord(self.ra * u.degree, self.dec * u.degree, frame='icrs')


    def get_shotlist(self, coords, radius=None, width=None, height=None):
        """
        
        Returns a list of shots for a 
        defined region of sky.
        
        self = Survey Class object        
        coords - astropy coordinate object
        radius - an astropy Quantity object, or a string 
                 that can be parsed into one.  e.g., '1 degree' 
                 or 1*u.degree.  If radius is specified, the 
                 shape is assumed to be a circle
        width -  in degress.  Specifies the edge length of a square box
        height - in degrees.  Specifies the height of a
                 rectangular box.  Must be passed with width.

        
        Example:
        
        S = Survey('hdr1')
        coords = SkyCoord(150.025513 * u.deg, 2.087767 * u.deg, frame='icrs')
        shotlist = S.get_shotlist(coords, radius=0.5*u.deg)
        
        or for a rectangular region around the point defined in coords

        shotlist = S.get_shotlist(coords, width=0.5, height=0.1)
        
        
        """
        
        try:
            self.coords
        except:
            self.coords = SkyCoord(self.ra * u.degree, self.dec * u.degree, frame='icrs')
    
        if radius:
            try:
                idx = self.coords.separation(coords) < radius
            except:
                print "Assuming radius in degrees"
                idx = self.coords.separation(coords) < radius * u.degree
        else:
            try:
                idx1 = (abs(self.ra  - coords.ra.value) < width/2.) 
                idx2 = (abs(self.dec  - coords.dec.value) < height/2.)
                idx = idx1 * idx2
            except:
                print "Provide both width and height of sky region in degrees."

        return self.shotid[idx]


    def close():
        '''
        Be sure to close the HDF5 file when you are done using
        it to release anything that might be in memory

        Example:

        S.close()
        '''
        
        self.hdfile.close()
