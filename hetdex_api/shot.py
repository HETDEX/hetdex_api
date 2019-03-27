# -*- coding: utf-8 -*-
"""

Functions to interact with the shot HDF5 file
and the Fibers Class. 
Requires a shotid or datevobs

author = Erin Mentuch Cooper
"""

import os.path as op
import glob
import re
import tables as tb
import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table, Column
import astropy.units as u
from astropy.coordinates import SkyCoord

import config

path_data = config.data_dir

def open_shot_file(shot):
    """
    get the file for a shot

    Input

    shot: either in the form of integer
          shotid (eg. 20180123009) or                                               
          string datevobs (eg. 20180123v009)

    Example:

    fileh = open_shot_file(20180123009)
    fileh = open_shot_file('20180123v009')
   
    """

    if re.search('v', str(shot)):
        file = op.join(path_data, str(shot)+'.h5')
    else:
        file = op.join(path_data, str(shot)[0:8]
                            + 'v' + str(shot)[8:11] + '.h5')
    fileh = tb.open_file(file, 'r')
    return fileh


class Fibers:
    def __init__(self, shot):
        '''
        Initialize Fibers Class

        This creates an astropy coordinates array of all
        fiber coordinates. Because of the large number of fibers
        for an individual shot this will take some time, but makes
        positional querying much easier.

        This will also initiate the wave_rect attribute which is
        an array of rectified wavelengths corresponding the the
        'calfib' and 'calfibe' datasets

        '''
        
        self.hdfile = open_shot_file(shot)
        self.table  = self.hdfile.root.Data.Fibers
        self.coords = SkyCoord(self.table.cols.ra[:] * u.degree,
                               self.table.cols.dec[:] * u.degree,
                               frame='icrs')
        self.wave_rect = 2.0 * np.arange(1036) + 3470.


    def query_region(self, coords, radius=3./3600.):
        """
        returns an indexed fiber table
        for a defined aperture
        
        self = Fibers class object
        coords = astropy coordinate object
        radius = astropy quantity object
        """
        
        idx = coords.separation(self.coords) < radius * u.degree
        
        return self.table[idx]

    def query_region_idx(self, coords, radius=3./3600.):
        """
        Returns an index for a Fibers class object to
        retrieve all fibers in the defined aperture

        self   - Fibers class object
        coords - astropy coordinate object
        radius - astropy quantity object                                                  
        """
        idx = coords.separation(self.coords) < radius * u.degree
        return np.where(idx)[0]

    def plot_fibertable_spectra(self, xlim=None, ylim=None):
        """
        Plots up series of spectra in a fibertable
        This could be messy if there are a lot of fibers!

        Inputs
        -----------
        self - a fibers table object, either 
               called by Fibers().table
               or fileh.root.Data.Fibers
        """
       
#        for row in self.table:

    def get_closest_fiber(self, coords):
        """
        Returns index to closest fiber to an RA/DEC
        """
        return coords.match_to_catalog_sky(self.coords)[0]
    
    
    def get_image_xy(self, idx, wave_obj):
        """
        Finds the X,Y image position from a row in the
        fibers table
        
        Note: this returns x,y image coords from 2D
        image arrays produced by Panacea. 
        """
            
        wave_data = self.table[idx]['wavelength']
        trace_data = self.table[idx]['trace']
        
        y = int(round(np.interp(wave_obj,wave_data,range(len(wave_data)))))
        x = int(round(np.interp(y,range(len(trace_data)),trace_data)))
        return x, y


    def plot_fiber_spectrum(self, idx, type='calfib', xlim=None, ylim=None):
        if type == 'calfib':
            try:
                plt.plot(self.wave_rect, self.table[idx]['calfib'])
                if xlim is not None:
                    plt.xlim(xlim)
                if ylim is not None:
                    plt.ylim(ylim)
                plt.xlabel('wavelength')
                plt.ylabel(type)
            except:
                print("Error plotting calib spectrum")
        else:
            try:
                plt.plot(self.table[idx]['wavelength'], self.table[idx][type])
                if xlim is not None:
                    plt.xlim(xlim)
                if ylim is not None:
                    plt.ylim(ylim)
                plt.xlabel('wavelength')
                plt.ylabel(type)
            except:
                print("Error plotting spectrum")
             
    def save_fiber_spectrum(self, idx, type='calfib', file='spec.dat'):
        
        spectab = Table()
        if type == 'calfib':
            try:
                print("Saving the flux-calibrated fiber spectrum")
                spectab['wavelength'] = self.wave_rect
                spectab[type] = self.table[idx][type]
                spectab['error'] = self.table[idx]['calfibe']
            except:
                print("Could not save calibrated fiber spectrum")
        else:
            try:
                spectab['wavelength'] = self.table[idx]['wavelength']
                spectab[type] = self.table[idx][type] 
                spectab['error'] = self.table[idx]['error1Dfib']
            except:
                print("Could not retrieve Fiber spectrum")
        spectab.write(file, format='ascii', overwrite=True)

    def close(self):
        self.hdfile.close()


def get_fibers_table(shot, coords, radius):
    """
    Returns fiber specta for defined aperture
    
    shot - either shotid or datevobs
    coords - astropy coordinate object
    radius - an astropy quantity object
    or radius in degrees
    
    """
    fileh = open_shot_file(shot)
    fibers = fileh.root.Data.Fibers
    ra_in = coords.ra.deg
    dec_in = coords.dec.deg
    rad = radius.degree

    fibers_table = fibers.read_where("sqrt((ra - ra_in)**2 + (dec - dec_in)**2) < rad")
    return fibers_table


def get_image2D_cutout(shot, coords, wave_obj, width=40, height=40, imtype='clean_image'):    
    """
    Returns an image from the 2D data based on
    ra/dec/wave.

    self - HDF5 file called from open_shot_file('shot')
    coords - astropy coordinate object
    wave_obj - astropy wavelength object
    imtype - image option to display
             options are 'clean_image', 'image', 'error':
    width - pixel width to be cutout (image size is 1032 pix)
    height - pixel height to be cutout (image size is 1032 pix)

    """
    fibers = Fibers(shot)

    idx = fibers.get_closest_fiber(coords)
    multiframe_obj = fibers.table.cols.multiframe[idx]
    expnum_obj = fibers.table.cols.expnum[idx]
    x, y = fibers.get_image_xy(idx, wave_obj)

    im0 = fibers.hdfile.root.Data.Images.read_where("(multiframe == multiframe_obj) & (expnum == expnum_obj)")
    
    return im0[imtype][0][x-int(width/2):x+int(width/2), y-int(height/2):y+int(height/2)]



def get_image2D_amp(shot, multiframe_obj, imtype='clean_image', expnum_obj=1):
    """
    Returns an image from the 2D data based on 
    an multiframe or a specid/amp combo
    
    multiframe - unique amp identifier to display
    imtype - image option to display
             options are:['spectrum', 'wavelength', 'fiber_to_fiber', 'twi_spectrum',
                         'sky_subtracted', 'trace', 'error1Dfib', 'calfib', 'calfibe',
                         'Amp2Amp', 'Throughput']
    expnum_obj - integer for which dither/exposure

    """
    fileh = open_shot_file(shot)
    im0 = fileh.root.Data.Images.read_where("(multiframe == multiframe_obj) & (expnum == expnum_obj)")
    
    return im0[imtype][0]
