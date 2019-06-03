# -*- coding: utf-8 -*-
"""

Initiates the Detections class.


Created on 2019/01/28

@author: Erin Mentuch Cooper
"""

from __future__ import print_function
from __future__ import unicode_literals

import sys
import os
import os.path as op
import numpy as np
import tables as tb
import copy
import matplotlib.pyplot as plt

from astropy.table import Table, Column
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
import pickle
import speclite.filters

from hetdex_api.survey import Survey
from hetdex_api import config

np.warnings.filterwarnings('ignore')

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
        if survey == 'hdr1':
            self.hdfile_elix = tb.open_file(config.elixerh5, mode='r')
            colnames2 = self.hdfile_elix.root.Classifications.colnames
            for name2 in colnames2:
                if name2 == 'detectid':
                    setattr(self, 'detectid_elix', self.hdfile_elix.root.Classifications.cols.detectid[:])
                else:
                    setattr(self, name2,
                            getattr(self.hdfile_elix.root.Classifications.cols, name2)[:])
        
        # also assign a field and some QA identifiers
        self.field = np.chararray(np.size(self.detectid),12)
        self.fwhm = np.zeros(np.size(self.detectid))
        self.flux_limit = np.zeros(np.size(self.detectid))
        self.throughput = np.zeros(np.size(self.detectid))
        self.n_ifu = np.zeros(np.size(self.detectid), dtype=int)

        S = Survey('hdr1')
        
        for index, shot in enumerate(S.shotid):
            ix = np.where(self.shotid == shot)
            self.field[ix] = S.field[index] #NOTE: python2 to python3 strings now unicode
            self.fwhm[ix] = S.fwhm_moffat[index]
            self.flux_limit[ix] = S.fluxlimit_4550[index]
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

        self.vis_class = -1 * np.ones(np.size(self.detectid))


        if survey == 'hdr1':
            self.add_hetdex_gmag(loadpickle=True, 
                                picklefile=config.gmags)
        elif survey == 'cont_sources':
            self.add_hetdex_gmag(loadpickle=True, 
                                 picklefile=config.gmags_cont)

    def __getitem__(self, indx):
        ''' 
        This allows for slicing of the Detections class
        object so that a mask can be applied to
        every attribute automatically by:
        
        detections_sliced = detects[indx]
        
        '''
        
        p = copy.copy(self)
        attrnames = self.__dict__.keys()
        for attrname in attrnames:
            try:
                setattr(p, attrname, getattr(self, attrname)[indx])
            except:
                setattr(p, attrname, getattr(self, attrname))
        return p

    def refine(self, gmagcut=None, removebalmerstars=False):
        '''
        Masks out bad and bright detections 
        and returns a refined Detections class
        object

        gmagcut = mag limit to exclude everything
                  brighter, defaults to None
        '''

        mask1 = self.remove_bad_amps() 
        mask2 = self.remove_bright_stuff(gmagcut)
        mask3 = self.remove_ccd_features()
        if removebalmerstars:
            mask4 = self.remove_balmerdip_stars()
        else:
            mask4 = np.ones(np.size(self.detectid), dtype=bool)
  
        mask5 = self.remove_bad_detects()
        mask6 = self.remove_shots()
        mask7 = self.remove_bad_pix()

        mask = mask1 * mask2 * mask3 * mask4 * mask5 * mask6 * mask7

        return self[mask]

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
            
        if limits.cont_low or limits.cont_high:
            maskcont = (self.continuum > limits.cont_low) * (self.continuum < limits.cont_high)
        else:
            maskcont = np.ones(ndets, dtype=bool)
            
        if limits.aperture_flag:
            coords = SkyCoord(limits.ra * u.degree, limits.dec * u.degree, frame='icrs')
            maskfield = self.query_by_coords(coords, limits.rad)
        else:
            maskfield = np.zeros(ndets, dtype=bool)
            print ('Subselecting for field(s):', limits.field)
        
            for field_index in limits.field:
                if field_index == 'all':
                    print ("Field = 'all'; not downselecting")
                    maskfield = np.ones(ndets, dtype=bool)
                else:
                    if isinstance(field_index, str):
                        #python2 to python3 issue (pytables also ... as bytes vs unicode)
                        mask_i = (self.field.decode() == field_index)
                    else:
                        mask_i = (self.field == field_index)
                    maskfield = maskfield | mask_i
        
        mask = maskwave * masklw * masksn * maskflux * maskchi2 * maskcont * maskfield
                    
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
        #todo: encoding='latin1' is an assumption ... might be better to use bytes?
        try:
            #if locally encoded and opened (same version of python)
            limits = pickle.load( open( picklefile, "rb" ))
        except:
            #if this is a python2 to python3 issue
            limits = pickle.load(open(picklefile, "rb"), encoding='bytes')
        mask = self.query_by_dictionary(limits)
        return mask

        
    def remove_bad_detects(self):
        '''
        Reads in the bad detect list from config.py
        and removes those detectids
        Set these to -2 (they are software errors)
        Don't use for machine learning or other
        classifying algorithms tests.
        '''
        
        # set an empty mask to start
        mask = np.zeros(np.size(self.detectid), dtype=bool)
        baddetects = np.loadtxt(config.baddetect, dtype=int)

        for baddet in baddetects:
            maskdet = self.detectid == baddet
            mask = mask | maskdet

        self.vis_class[mask] = -2

        return np.invert(mask)


    def remove_bad_amps(self):
        '''
        Reads in the bad amp list from config.py
        and creates a mask to remove those detections.
        It will also assign a 0 to signify an artifact
        in vis_class
        '''

        # set an empty mask to start

        mask = np.zeros(np.size(self.detectid), dtype=bool)

        badamps = ascii.read(config.badamp, 
                             names=['ifuslot', 'amp','date_start', 'date_end'])

        for row in np.arange(np.size(badamps)):
            if badamps['amp'][row] == 'AA':
                 maskamp = (self.ifuslot == str(badamps['ifuslot'][row]).zfill(3)) * (self.date >= badamps['date_start'][row]) * (self.date <= badamps['date_end'][row])
                 mask = maskamp | mask
            else:
                maskamp = (self.amp == badamps['amp'][row]) * (self.ifuslot == str(badamps['ifuslot'][row]).zfill(3)) * (self.date >= badamps['date_start'][row]) * (self.date <= badamps['date_end'][row])
                mask = maskamp | mask

        self.vis_class[mask] = 0

        return np.invert(mask)

    def remove_bad_pix(self):
        '''
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
        
        '''
        
        badpixlist = ascii.read(config.badpix,
                                names=['multiframe','x1','x2','y1','y2'])
    
        mask = np.zeros(np.size(self.detectid), dtype=bool)
        
        for row in badpixlist:
            maskbadpix = (self.multiframe == row['multiframe']) * (self.x_raw > row['x1']) * (self.x_raw < row['x2']) * (self.y_raw > row['y1']) * (self.y_raw < row['y2'])
            mask = maskbadpix | mask
            
        self.vis_class[mask] = 0
        
        return np.invert(mask)


    def remove_shots(self):
        '''
        Takes a list of bad shots and removes them. Assigns -2 
        to detections in these shots so they are not used 
        in any MLing analysis
        '''

        mask = np.zeros(np.size(self.detectid), dtype=bool)
        badshots = np.loadtxt(config.badshot, dtype=int)

        for shot in badshots:
            maskshot = (self.shotid == shot)
            mask = mask | maskshot
        
        self.vis_class[mask] = -2

        return np.invert(mask)


    def remove_balmerdip_stars(self):
        '''
        Applies a cut to the databased to remove all
        stars that show a false emission feature around 3780
        Also assigns a star classification to each detectid

        This is obsolete as gmag cuts get rid of these easily
        '''
        mask1 = (self.wave > 3775) * (self.wave < 3785) * (self.continuum > 3)

        mask2 = (self.wave > 4503.) * (self.wave <4513.) * (self.continuum > 10)
 
        mask = mask1 | mask2 

        self.vis_class[mask] = 3

        return np.invert(mask)

        
    def remove_bright_stuff(self, gmagcut):
        '''
        Applies a cut to remove bright stars based on
        gmag attribute. Assigns a star classification
        to each detectid, although some may be nearby
        galaxies. Will want to improve this if you are
        looking to study nearby galaxies.
        '''

        if gmagcut:
            mask = (self.gmag < gmagcut) 
        else:
            mask = np.zeros(np.size(self.detectid), dtype=bool)
        
        return np.invert(mask)
    
    def remove_ccd_features(self):
        '''
        Remove all objects with very at the
        edges of the detectors and denotes 
        them as artifacts in vis_class
        '''

        mask = (self.wave < 3540.) | (self.wave > 5515.) 
        self.vis_class[mask] = 0
        
        return np.invert(mask)

    def get_spectrum(self, detectid_i):
        spectra = self.hdfile.root.Spectra
        spectra_table = spectra.read_where("detectid == detectid_i")
        data = Table()
        intensityunit = u.erg / (u.cm ** 2 * u.s * u.AA)
        data['wave1d'] = Column( spectra_table['wave1d'][0], unit= u.AA)
        data['spec1d'] = Column( spectra_table['spec1d'][0], unit= 1.e-17 * intensityunit)
        data['spec1d_err'] = Column( spectra_table['spec1d_err'][0], unit= 1.e-17 * intensityunit)

        # convert from 2AA binning to 1AA binning:
        data['spec1d'] /= 2.
        data['spec1d_err'] /= 2.

        return data        

    def get_gband_mag(self, detectid_i):
        '''
        Calculates the gband magnitude from the 1D spectrum
        '''

        spec_table = self.get_spectrum(detectid_i)
        gfilt = speclite.filters.load_filters('sdss2010-g')
        flux, wlen = gfilt.pad_spectrum(np.array(1.e-17*spec_table['spec1d']), 
                                        np.array(spec_table['wave1d']))
        mag = gfilt.get_ab_magnitudes(flux, wlen)
        
        return mag

    def add_hetdex_gmag(self, loadpickle=True, picklefile='gmags.pickle'):
        '''
        Calculates g-band magnitude from spec1d for
        each detectid in the Detections class instance
        If the gmags.pickle file and pickle=True is 
        given then it will just load from previous computation
        '''
        if loadpickle:
            # todo: encoding='latin1' is an assumption ... might be better to use bytes?
            self.gmag = pickle.load(open(picklefile, 'rb'),encoding="bytes")
        else:
            self.gmag = np.zeros(np.size(self.detectid), dtype=float)
        
            # fastest way is to calculate gmag over whole 1D spec array
            # then populate the detections class

            detectid_spec = self.hdfile.root.Spectra.cols.detectid[:]
            spec1d = self.hdfile.root.Spectra.cols.spec1d[:]
            # convert from ergs/s/cm2 to ergs/s/cm2/AA
            spec1d /= 2.
            wave_rect =  2.0 * np.arange(1036) + 3470.
            gfilt = speclite.filters.load_filters('sdss2010-g')
            flux, wlen = gfilt.pad_spectrum(1.e-17*spec1d, wave_rect)
            gmags = gfilt.get_ab_magnitudes(flux, wlen)

            for idx, detectid_i in enumerate(self.detectid[:]):
                seldet = np.where(detectid_spec == detectid_i)[0]
                if np.size(seldet) == 1:
                    self.gmag[idx] = gmags[seldet][0][0]
                else:
                    self.gmag[idx] = np.nan


    def get_hetdex_mag(self, detectid_i, filter='sdss2010-g'):
        ''' 
        filter = can be any filter used in the speclite
                 package 
                 https://speclite.readthedocs.io/en/latest/api.html

        '''
        spec_table = self.get_spectrum(detectid_i)
        filt = speclite.filters.load_filters(filter)
        flux, wlen = filt.pad_spectrum(np.array(1.e-17*spec_table['spec1d']),
                                        np.array(spec_table['wave1d']))
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

        return table

    def save_spectrum(self, detectid_i, outfile=None):
        spec_data = self.get_spectrum(detectid_i)
        if outfile:
            ascii.write(spec_data, outfile, overwrite=True)
        else:
            ascii.write(spec_data, 'spec_'+str(detectid_i)+'.dat', overwrite=True)


    def plot_spectrum(self, detectid_i, xlim=None, ylim=None):
        spec_data = self.get_spectrum(detectid_i)
        plt.figure(figsize=(8, 6))
        plt.errorbar(spec_data['wave1d'], spec_data['spec1d'], yerr=spec_data['spec1d_err'])
        plt.title("DetectID "+str(detectid_i))
        plt.xlabel('wave (AA)')
        plt.ylabel('flux (1e-17 erg/s/cm^2/AA)')
        if xlim is not None:
            plt.xlim(xlim)
            if ylim is not None:
                plt.ylim(ylim)
        plt.show()


def show_elixer(detectid):
    '''
    Takes a detectid and pulls out the elixer PDF from the
    elixer tar files on hdr1 and shows it in matplotlib
    '''
    elix_dir =  '/work/05350/ecooper/stampede2/elixer/jpgs/'
    file_jpg = op.join(elix_dir, "egs_%d" %(detectid//100000), str(detectid) + '.jpg')
    plt.figure(figsize=(10,8))
    im = plt.imread(file_jpg)
    plt.imshow(im)
