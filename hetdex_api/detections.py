# -*- coding: utf-8 -*-
"""

Initiates the Detections class.
An API may or may not be developed in the future.
We recommend just using pytables since this is only one file.

Created on 2019/01/28

@author: Erin Mentuch Cooper
"""

import numpy as np
import tables as tb
from astropy.table import Table, Column
import astropy.units as u
from astropy.coordinates import SkyCoord
import pickle

from hetdex_api.survey import Survey
from hetdex_api import config


class Detections:
    def __init__(self, survey):
        '''
        Initialize the detection catalog class for a given data release

        Input
        -----
        survey : string
            Data release you would like to load, i.e., 'DR1' or 'cont_sources'.
            This is case insensitive.
        '''
        survey_options = {'hdr1': config.detecth5,
                          'cont_sources': config.contsourceh5}
        if survey.lower() not in survey_options:
            print('survey not in survey options')
            print(survey_options)
            return None

        self.filename = survey_options[survey]
        self.hdfile = tb.open_file(self.filename, mode='r')
        colnames = self.hdfile.root.Detections.colnames
        for name in colnames:
            setattr(self, name,
                    getattr(self.hdfile.root.Detections.cols, name)[:])
        
        # set the SkyCoords
        self.coords = SkyCoord(self.ra * u.degree, self.dec * u.degree, frame='icrs')

        # add in the elixer probabilties and associated info:
        self.hdfile_elix = tb.open_file(config.elixerh5, mode='r')
        colnames2 = self.hdfile_elix.root.Classifications.colnames
        for name2 in colnames2:
            if name2 == 'detectid':
                setattr(self, 'detectid_elix', self.hdfile_elix.root.Classifications.cols.detectid[:])
            else:
                setattr(self, name2,
                        getattr(self.hdfile_elix.root.Classifications.cols, name2)[:])
        
        # also assign a field identifier
        self.field = np.chararray(np.size(self.detectid),12)
        S = Survey(survey)
        for index, shot in enumerate(S.shotid):
            ix = np.where(self.shotid == shot)
            self.field[ix] = S.field[index]


    def query_by_coords(self, coords, radius):
        '''
        Returns mask based on a coordinate search
        
        self = Detections Class object   
        coords - astropy coordinate object
        radius - an astropy Quantity object, or a string 
        that can be parsed into one.  e.g., '1 degree' 
        or 1*u.degree. Will assume arcsec if no units given
        
        '''
        sep = self.coords.separation(coords) 
        try:
            maskcoords = sep < radius
        except:
            maskcoords = sep.arcmin < radius
        return maskcoords


    def query_by_dictionary(self, limits):
        
        '''
        Takes a dictionary of query limits
        and reduces the detections database
        
        '''
        
        ndets = np.size(self.detectid)
        
        if limits.wave_low or limits.wave_high: 
            maskwave = (self.wave > limits.wave_low) * (self.wave < limits.wave_high)
        else: 
            maskwave = np.ones(np.size(ndets), dtype=bool)
            
        if limits.flux_low or limits.flux_high:
            maskflux = (self.flux > limits.flux_low) * (self.flux < limits.flux_high)
        else:
            maskflux = np.ones(np.size(ndets), dtype=bool)
            
        if limits.linewidth_low or limits.linewidth_high:
            masklw = (self.linewidth > limits.linewidth_low) * (self.linewidth < limits.linewidth_high)
        else:
            masklw = np.ones(np.size(ndets), dtype=bool)
            
        if limits.sn_low or limits.sn_high:
            masksn = (self.sn > limits.sn_low) * (self.sn < limits.sn_high)
        else:
            masksn = np.ones(np.size(ndets), dtype=bool)
            
        if limits.chi2_low or limits.chi2_high:
            maskchi2 = (self.chi2 > limits.chi2_low) * (self.chi2 < limits.chi2_high)
        else:
            maskchi2 = np.ones(np.size(ndets), dtype=bool)
            
            
        if limits.aperture_flag:
            coords = SkyCoord(limits.ra * u.degree, limits.dec * u.degree, frame='icrs')
            maskfield = self.query_by_coords(coords, limits.rad)
        else:
            try:
                maskfield = np.ones(np.size(ndets), dtype=bool)
                
                if any(limits.field) != 'all':
                    maskfield = np.zeros(np.size(ndets), dtype=bool)
                    print 'Subselecting for field(s):', limits.field
                    for field_index in limits.field:
                        mask_i = (self.field == field_index)
                        maskfield = maskfield | mask_i
            except:
                print('Error downselecting based on field selection')

        mask = maskwave * masklw * masksn * maskflux * maskchi2 * maskfield
        
        return mask
                
    def query_by_pickle(self, picklefile):
        '''
        
        this function queries the Detections class
        object based on query made with detwidgets.py
        
        Input
        
        self = a detections class object
        picklefile = string filename for a pickle
        created in detwidget.py
        
        '''
        limits = pickle.load( open( picklefile, "rb" ) )
        mask = self.query_by_dictionary(limits)
        return mask

        
    def remove_bad_detects(self):
        '''
        Reads in the bad detect list from config.py
        and removes those detectids
        '''


    def remove_shots(self, shotlist):
        '''
        Takes a list of shots and removes them
        ''

