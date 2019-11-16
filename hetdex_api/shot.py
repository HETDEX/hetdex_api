# -*- coding: utf-8 -*-
"""

Functions to interact with the shot HDF5 file
and the Fibers Class. 
Requires a shotid or datevobs

author = Erin Mentuch Cooper
"""

import os.path as op
import re
import tables as tb
import numpy as np

import warnings
import sys
if not sys.warnoptions:
    warnings.simplefilter("ignore")

from astropy.table import Table, Column
import astropy.units as u
from astropy.coordinates import SkyCoord

from hetdex_api import config

path_data = config.data_dir

def open_shot_file(shot):
    """
    get the file for a shot. This is a
    global function. It basically just
    replaces dealing with a path name.

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
        'calfib', 'calfibe', 'Amp2Amp', and 'Throughput' datasets

        '''

        self.hdfile = open_shot_file(shot)
        self.table  = self.hdfile.root.Data.Fibers
        self.coords = SkyCoord(self.table.cols.ra[:] * u.degree,
                               self.table.cols.dec[:] * u.degree,
                               frame='icrs')
        self.wave_rect = 2.0 * np.arange(1036) + 3470.

        colnames = self.hdfile.root.Data.Fibers.colnames
        
        for name in colnames:
            if isinstance(getattr(self.hdfile.root.Data.Fibers.cols, name)[0], np.bytes_):
                setattr(self, name,
                        getattr(self.hdfile.root.Data.Fibers.cols, name)[:].astype(str))
            else:
                setattr(self, name,
                        getattr(self.hdfile.root.Data.Fibers.cols, name)[:])
        

    def query_region(self, coords, radius=3./3600.):
        """
        returns an indexed fiber table
        for a defined aperture.

        self = Fibers class object
        coords = astropy coordinate object
        radius = radius in degrees
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


    def get_closest_fiber(self, coords, exp=None):
        """
        Returns index to closest fiber to an RA/DEC
        exp=dither number (used primarily for building 2D image
        """
        if exp in [1,2,3]:
            sel = self.expnum=exp
            fib_idx = coords.match_to_catalog_sky(self.coords[sel])[0] 
        else:
            fib_idx = coords.match_to_catalog_sky(self.coords)[0]
        return fib_idx


    def get_image_xy(self, idx, wave_obj):
        """
        Finds the X,Y image position from a row in the
        fibers table

        Note: this returns x,y image coords from 2D
        image arrays produced by Panacea.
        """

        wave_data = self.wavelength[idx]
        trace_data = self.trace[idx]

        y = int(np.round(np.interp(wave_obj,wave_data,range(len(wave_data)))))
        x = int(np.round(np.interp(y,range(len(trace_data)),trace_data)))
        return x, y


    def plot_fiber_spectrum(self, idx, type='calfib', xlim=None, ylim=None):
        import matplotlib.pyplot as plt

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
        """
        Saves a fiber spectrum

        self = Fibers class object
        idx = index of the fiber in the Fibers Table
        types = ['calfib', 'spectrum', 'sky_spectrum', 'twi_spectrum',
                 'error1Dfib', fiber_to_fiber']
        file = output file name

        """

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

    
    def get_fib_image2D(self, wave_obj, fibnum_obj, multiframe_obj, expnum_obj, width=60, height=40, imtype='clean_image'):
        """
        Returns an image from the 2D data for 
        a specific fiber centered at a specifc wavelength

        self - fibers class object
        wave_obj - astropy wavelength object
        expnum - dither exposure number [1,2,3]
        fibnum - fiber number
        multiframe - amp multiframe ID

        imtype - image option to display
        width - pixel width to be cutout (image size is 1032 pix)
        height - pixel height to be cutout (image size is 1032 pix)

        """
        
        idx = np.where((self.fibidx == (fibnum_obj - 1) ) * (self.multiframe == multiframe_obj) * (self.expnum == expnum_obj))[0][0]
        
        x, y = self.get_image_xy(idx, wave_obj)
            
        im0 = self.hdfile.root.Data.Images.read_where("(multiframe == multiframe_obj) & (expnum == expnum_obj)")

        #create image of forced dims of input width x height

        height= np.minimum(height, 1032)
        width = np.minimum(width, 1032)
        
        im_base = np.zeros((height, width))

        dx = int(height/2)
        dy = int(width/2)

        x1 = np.maximum(0, x-dx)
        x2 = np.minimum(x+ dx + (height % 2), 1032)

        y1 = np.maximum(0, y-dy)
        y2 = np.minimum(y+ dy + (width % 2), 1032)
        
        x1_slice = np.minimum(0, height - (x2 - x1))
        x2_slice = x2-x1
        y1_slice = np.minimum(0, width - (y2 - y1))
        y2_slice = y2-y1
        
        im_reg = im0[imtype][0][x1:x2,y1:y2]

        im_base[ x1_slice:x2_slice, y1_slice:y2_slice] = im_reg
        
        return im_base

    def return_astropy_table(self):
        """
        Return an astropy table version of the Fibers table
        that can easily be saved

        Returns
        -------
        table : astropy.table:Table
            an astropy table you can save

        """
        table = Table()
        for name in self.hdfile.root.Data.Fibers.colnames:
            table[name] = getattr(self, name)
        
        return table

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
    try:
        ra_in = coords.ra.degree
        dec_in = coords.dec.degree
    except:
        print("Coords argument must be an astropy coordinates object")
        
    try:
        rad_in = radius.to(u.degree)
        rad = radius
    except:
        print('Assuming radius in arcsec')
        rad_in = radius/3600.
        rad  = radius * u.arcsec
        pass
        
    #search first along ra 

    ra_table = fibers.read_where("sqrt((ra - ra_in)**2) < (rad_in + 2./3600)")

    if any(ra_table):
        coords_table = SkyCoord(ra_table['ra']*u.deg, ra_table['dec']*u.deg, frame='icrs')
        idx = coords.separation(coords_table) < rad
        print(idx)
        fibers_table = ra_table[idx]
    else:
        fibers_table = None

    fileh.close()
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
    multiframe_obj = fibers.table.cols.multiframe[idx].astype(str)
    expnum_obj = fibers.table.cols.expnum[idx]
    x, y = fibers.get_image_xy(idx, wave_obj)

    im0 = fibers.hdfile.root.Data.Images.read_where("(multiframe == multiframe_obj) & (expnum == expnum_obj)")

    return im0[imtype][0][x-int(width/2):x+int(width/2), y-int(height/2):y+int(height/2)]



def get_image2D_amp(shot, multiframe_obj, imtype='clean_image', expnum_obj=1):
    """
    Returns an image from the 2D data based on
    an multiframe or a specid/amp/expnum combo

    multiframe - unique amp identifier to display
    imtype - image option to display
             options are:['spectrum', 'wavelength', 'fiber_to_fiber', 'twi_spectrum',
                         'sky_subtracted', 'trace', 'error1Dfib', 'calfib', 'calfibe',
                         'Amp2Amp', 'Throughput']
    expnum_obj - integer for the dither/exposure

    """
    fileh = open_shot_file(shot)
    im0 = fileh.root.Data.Images.read_where("(multiframe == multiframe_obj) & (expnum == expnum_obj)")
    fileh.close()

    return im0[imtype][0]
