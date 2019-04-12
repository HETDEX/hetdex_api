# -*- coding: utf-8 -*-
"""

Initiates the Detections class.
An API may or may not be developed in the future.
We recommend just using pytables since this is only one file.

Created on 2019/01/28

@author: Erin Mentuch Cooper
"""

import sys
import os
import os.path as op
import numpy as np
import tables as tb
from astropy.table import Table, Column
import astropy.units as u
from astropy.coordinates import SkyCoord
import pickle
import tarfile

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
        self.fwhm = np.zeros(np.size(self.detectid))
        self.flux_limit = np.zeros(np.size(self.detectid))
        self.throughput = np.zeros(np.size(self.detectid))

        S = Survey(survey)
        for index, shot in enumerate(S.shotid):
            ix = np.where(self.shotid == shot)
            self.field[ix] = S.field[index]
            self.fwhm[ix] = S.fwhm_moffat[index]
            self.flux_limit[ix] = S.fluxlimit_4550[index]
            self.throughput[ix] = S.response_4540[index]
    

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
        '''
        
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
            masklw = (self.linewidth > limits.linewidth_low) * (self.linewidth < limits.linewidth_high)
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
            
            
        if limits.aperture_flag:
            coords = SkyCoord(limits.ra * u.degree, limits.dec * u.degree, frame='icrs')
            maskfield = self.query_by_coords(coords, limits.rad)
        else:
            maskfield = np.zeros(ndets, dtype=bool)
            print 'Subselecting for field(s):', limits.field
        
            for field_index in limits.field:
                if field_index == 'all':
                    print "Field = 'all'; not downselecting"
                    maskfield = np.ones(ndets, dtype=bool)
                else:
                    mask_i = (self.field == field_index)
                    maskfield = maskfield | mask_i
        
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
        '''

    def remove_balmerdip_stars(self):
        '''
        Applies a cut to the databased to remove all
        stars that show a false emission feature around 3780
        '''
        maskwave = (self.wavelength > 3775) * (self.wavelength < 3785) * self.continuum > 3)
        
        return np.invert(maskwave)
        
def show_elixer(detectid):
    '''
    Takes a detectid and pulls out the elixer PDF from the
    elixer tar files on hdr1 and shows it in matplotlib
    '''

    tarname = op.join(config.elix_dir,"erg_%d.tar" %(detectid//100000))
    if detectid < 1000690662:
        file_pdf = op.join("egs_%d" %(detectid//100000), str(detectid) + '.pdf')
    elif detectid >= 1000690662 and detectid <= 1000690798:
        file_pdf = op.join(str(detectid) + '.pdf')

    os.system('tar -xvf '+ tarname + ' ' + file_pdf) 
    os.system('pdftoppm ' + file_pdf + ' ' + str(detectid) + ' -png -singlefile')
    file_png = str(detectid) + '.png'
    elixim = plt.imread(file_png)
     
