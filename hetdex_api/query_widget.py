from __future__ import print_function

"""

Widget to query HETDEX spectra via elixer catalog API and 
HETDEX API tools

Authors: Erin Mentuch Cooper

Date: November 9, 2019

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import os.path as op

from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import NDData

from astrowidgets import ImageWidget
import ipywidgets as widgets
from ipywidgets import interact, Layout 

from hetdex_api.shot import *
from get_spec import get_spectra

from astroquery.sdss import SDSS

sys.path.append('/work/03946/hetdex/hdr1/software/elixer')
import catalogs

class QueryWidget():

    def __init__(self, coords=None, survey='hdr1', aperture=3.*u.arcsec, cutout_size=5.*u.arcmin):

        self.survey = survey.lower()

        self.aperture = aperture

        self.cutout_size = cutout_size

        self.catlib = catalogs.CatalogLibrary()

        if coords:
            self.coords = coords
        else:
            self.coords = SkyCoord(210.087982 * u.deg, 51.701839 * u.deg, frame='icrs')

        #initialize the image widget from astrowidgets
        self.imw = ImageWidget(image_width=400, image_height=400)

        self.survey_widget = widgets.Dropdown(options=['HDR1', 'HDR2'], value=self.survey.upper(), layout=Layout(width='10%'))        
        self.im_ra = widgets.FloatText(value=self.coords.ra.value, description='RA (deg):', layout=Layout(width='20%'))
        self.im_dec = widgets.FloatText(value=self.coords.dec.value, description='DEC (deg):', layout=Layout(width='20%'))
        self.pan_to_coords = widgets.Button(description="Pan to coords", disabled=False, button_style='success')
        self.marking_button = widgets.Button(description='Mark Sources', button_style='success') 
        self.marker_table 
        
        self.textimpath = widgets.Text(description='Source: ', value='', layout=Layout(width='90%'))

        self.load_image()

        self.leftbox = widgets.VBox([self.imw, self.textimpath])

        self.rightbox = widgets.VBox([self.marking_button, self.marker_table])

        display(widgets.HBox([self.survey_widget, self.im_ra, self.im_dec, self.pan_to_coords]))
        display(widgets.HBox([self.leftbox, self.rightbox]))
                
        self.pan_to_coords.on_click(self.pan_to_coords_click)

        self.marking_button.on_click(self.marking_on_click)


    def update_coords(self):
        self.coords = SkyCoord(self.im_ra.value * u.deg, self.im_dec.value * u.deg, frame='icrs')

    def pan_to_coords_click(self, b):
        self.update_coords()
        if self.coords.separation(self.orig_coords) < self.cutout_size:
            self.imw.center_on(self.coords)
        else:
            self.load_image()
            
    def load_image(self):

        self.imw.zoom = 3
        im_size = self.cutout_size.to(u.arcsec).value
        mag_aperture = self.aperture.to(u.arcsec).value

        self.cutouts = self.catlib.get_cutouts(position=self.coords, radius=im_size, aperture=mag_aperture, dynamic=False)
        
        # keep original coords of image for bounds checking later
        self.orig_coords = self.coords

        cutout_index = -1

        for index in np.arange(0, np.size(self.cutouts)):
            if self.cutouts[index]['filter'] == 'g':
                cutout_index = index
                break
            else:
                cutout_index = 0

        if cutout_index >= 0:    

            im = NDData( self.cutouts[cutout_index]['cutout'].data, wcs=self.cutouts[cutout_index]['cutout'].wcs)
            self.im_path = self.cutouts[cutout_index]['path']
            self.imw.load_nddata(im)

        else:
            try:
                sdss_im = SDSS.get_images(coordinates=self.coords, band='g')
                im = sdss_im[0][0]
            except:
                sdss_im = SDSS.get_images(coordinates=self.coords, band='g', radius=30.*u.arcsec)
                im = sdss_im[0][0]
            
            self.im_path = "SDSS Astroquery result"
            self.imw.load_fits(im)

        self.textimpath.value= self.im_path
            

    def marking_on_click(self, b):

        if self.marking_button.button_style=='success':

            self.imw.start_marking(marker={'color': 'red', 'radius': 3, 'type': 'circle'},
                                   marker_name='clicked markers'
                               )
            self.marking_button.description='Stop Marking'
            self.marking_button.button_style='danger'
        
        else:
            self.imw.stop_marking()
            self.marking_button.description='Mark Sources'
            self.marking_button.button_style='success'

#        self.marker_tab 
