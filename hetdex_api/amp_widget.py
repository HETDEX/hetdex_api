#!/usr/bin/env python
# coding: utf-8
#
# AmpWidget Widget Object
# Author: Erin Mentuch Cooper
# Date: 2020-11-11

import numpy as np

import tables as tb

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astrowidgets import ImageWidget
import ipywidgets as widgets
from ipywidgets import Layout

from hetdex_api.config import HDRconfig
from hetdex_api.shot import Fibers, get_image2D_amp, open_shot_file
from hetdex_api.survey import Survey, FiberIndex


try:  # using HDRconfig
    LATEST_HDR_NAME = HDRconfig.LATEST_HDR_NAME
    CONFIG_HDR2 = HDRconfig('hdr2.1')
    CONFIG_HDR3 = HDRconfig('hdr3')
    CONFIG_HDR4 = HDRconfig('hdr4')
    CONFIG_HDR5 = HDRconfig('hdr5')
    
except Exception as e:
    print("Warning! Cannot find or import HDRconfig from hetdex_api!!", e)
    LATEST_HDR_NAME = "hdr2.1"

OPEN_DET_FILE = None
DET_HANDLE = None

class AmpWidget:
    def __init__(
        self,
        survey=LATEST_HDR_NAME,
        coords=None,
        radius=3.0,
        detectid=None,
        wave=None,
        shotid=None,
        multiframe=None,
        imtype="clean_image",
        expnum=1,
    ):

        self.survey = survey.lower()
        
        self.detectid = detectid
        self.shotid = shotid
        self.multiframe = multiframe

        if self.shotid is None:
            if self.multiframe is None:
                if self.detectid is None:
                    self.detectid = 4013993285 #random default LAE
            else:
                print('You must provide a shotid/multiframe combo or a detectid')
                
        self.imtype = imtype
        self.expnum = expnum
        self.coords = coords
        self.radius = radius
        self.wave = wave
        self.detfile = None
        
        self.survey_widget = widgets.Dropdown(
                        options=["HDR1","HDR2.1", "HDR3", "HDR4", "HDR5"],
                        value=self.survey.upper(),
                        layout=Layout(width="10%"),
                    )

        # populate hetdex_api config with path files
        self.update_hdr_config()

        # initialize the image widget from astrowidgets
        self.imw = ImageWidget()  # image_width=600, image_height=600)

        self.detectbox = widgets.BoundedIntText(
            value=self.detectid,
            min=2100000000,
            max=6000000000,
            step=1,
            description="DetectID:",
            disabled=False,
            # layout=Layout(width='15%')
        )

        self.det_button = widgets.Button(
            description="Go",
            toolkit="This will bring up the amp for the highest weight fiber",
            button_style="success",
            disabled=False,
        )

        if self.coords is not None:
            self.im_ra = widgets.FloatText(
                value=self.coords.ra.value, description="RA (deg):"
            )
            self.im_dec = widgets.FloatText(
                value=self.coords.dec.value, description="DEC (deg):"
            )
        else:
            self.im_ra = widgets.FloatText(value=0.0, description="RA (deg):")
            self.im_dec = widgets.FloatText(value=0.0, description="DEC (deg):")

        self.select_coords = widgets.Button(
            description="Query coords",
            disabled=False,
            button_style="success",
            toolkit="This will narrow down amps to those within the defined region",
        )

        if self.wave is None:
            wave0 = 3500
        else:
            wave0 = self.wave

        self.wave_widget = widgets.FloatText(
            value=wave0, description="Wave (A)", min=3500.0, max=5500.0, step=1.0
        )

        self.show_button = widgets.Button(
            description="Show Region",
            toolkit="This will highlight the RA/DEC/WAVE input",
            button_style="success",
            disabled=False,
        )

        self.shotid_widget = widgets.Dropdown(
            description="ShotID", options=self.survey_class.shotid, value=self.shotid,
        )

        if self.shotid is not None:
#            self.shoth5 = open_shot_file(self.shotid, survey=self.survey)
            sel_shot = self.ampflag_table["shotid"] == self.shotid
            mflist = np.unique(self.ampflag_table["multiframe"][sel_shot])
            self.multiframe_widget = widgets.Dropdown(
                description="MultiframeID", options=mflist, value=mflist[0]
            )
        else:
            mflist = np.unique(self.ampflag_table["multiframe"])
            self.multiframe_widget = widgets.Dropdown(
                description="MultiframeID", options=mflist, value=mflist[0]
            )

        if self.multiframe is not None:
            self.multiframe_widget.value = self.multiframe
        
        self.expnum_widget = widgets.Dropdown(
            description="ExpNum", options=[1, 2, 3], value=self.expnum,
        )

        self.imtype_widget = widgets.Dropdown(
            description="ImType",
            options=["clean_image", "image", "error"],
            value=self.imtype,
        )

        self.topbox = widgets.HBox(
            [
                self.survey_widget,
                self.detectbox,
                self.det_button,
                #                self.im_ra,
                #                self.im_dec,
                #                self.select_coords,
                #                self.wave_widget,
            ]
        )
        box_layout = Layout()

        label1 = widgets.Label(value="Enter a coordinate to down-select shots")
        label2 = widgets.Label(value="Find a coordinate/wavelength region")

        self.get_cursor_button = widgets.Button(
                        description="Get RA/DEC/WAVE",
                        toolkit="This will map XY coordinate to RA/DEC/WAVE above",
                        button_style="success",
                        disabled=False,
                    )
        
        
        self.boxside = widgets.VBox(
            [
                label1,
                self.im_ra,
                self.im_dec,
                self.select_coords,
                label2,
                self.wave_widget,
                self.show_button,
                self.shotid_widget,
                self.multiframe_widget,
                self.expnum_widget,
                self.imtype_widget,
                self.get_cursor_button,
            ]
        )
        self.midbox = widgets.HBox([self.imw, self.boxside], layout=box_layout)

        self.bottombox = widgets.Output(layout={"border": "1px solid black"})
        
        if self.coords is not None:
            if self.detectid is None:
                
                self.shotid = None
                
                fiber_table_region = self.FibIndex.query_region(
                    self.coords, radius=self.radius * u.arcsec, shotid=self.shotid
                )
                shotlist = np.unique(fiber_table_region["shotid"])
                
                self.shotid_widget = widgets.Dropdown(
                    options=shotlist, value=shotlist[0]
                )
            else:
                pass
        elif self.detectid is not None:
            # open up detections h5 file to get info for
            # highest weight fiber
            self.get_amp_info_from_det()

        self.im = get_image2D_amp(
            self.shotid_widget.value,
            multiframe=self.multiframe,
            imtype=self.imtype,
            expnum=self.expnum,
            survey=self.survey,
        )
        
        self.imw.load_array(self.im)

        display( widgets.VBox([self.topbox, self.midbox, self.bottombox]))

        # action calls
        self.survey_widget.observe(self.survey_widget_change)
        self.multiframe_widget.observe(self.im_widget_change)
        self.expnum_widget.observe(self.im_widget_change)
        self.imtype_widget.observe(self.im_widget_change)
        self.shotid_widget.observe(self.shotid_widget_change)
        self.select_coords.on_click(self.coord_change)
        self.show_button.on_click(self.show_region)
        self.det_button.on_click(self.on_det_go)
        self.get_cursor_button.on_click(self.get_ra_dec_wave)

    def survey_widget_change(self, b):
        self.update_hdr_config()
        self.shotid_widget.options = self.survey_class.shotid
        self.shotid_widget.value = self.shotid
        self.multiframe_widget.value = self.multiframe
        self.update_amp_image()
        
    def im_widget_change(self, b):
        self.update_amp_image()

    def shotid_widget_change(self, b):
        self.bottombox.clear_output()
        #self.shoth5.close()
        self.coords = None
        self.detectid = None
        self.im_ra.value = 0.0
        self.im_dec.value = 0.0
        
        self.shotid = self.shotid_widget.value
        self.multiframe = self.multiframe_widget.value
        
        #self.shoth5 = open_shot_file(self.shotid_widget.value,
        #                             survey=self.survey)
        
        sel_shot = self.ampflag_table["shotid"] == self.shotid
        mflist = np.unique(self.ampflag_table["multiframe"][sel_shot])

        self.multiframe_widget.options = mflist
        
        try:
            sel = (self.ampflag_table["shotid"] == self.shotid) * (
                self.ampflag_table["multiframe"] == self.multiframe
            )
            flag = self.ampflag_table["flag"][sel][0]
        except:
            self.multiframe = mflist[0]
            self.multiframe_widget.value = mflist[0]
            sel = (self.ampflag_table["shotid"] == self.shotid) * (
                self.ampflag_table["multiframe"] == self.multiframe
            )
            flag = self.ampflag_table["flag"][sel][0]

        if flag:
            self.midbox.layout = Layout()
        else:
            box_layout = Layout(
                display="flex",
                # flex_flow='column',
                # align_items='stretch',
                border="5px solid red",
            )
            self.midbox.layout = box_layout

        self.update_amp_image()

    def update_amp_image(self):
        self.bottombox.clear_output()
        # add in bad amp check
        try:
            self.imw.remove_markers()
        except Exception:
            pass

        self.shotid = self.shotid_widget.value
        self.multiframe = self.multiframe_widget.value
        self.expnum = self.expnum_widget.value
        try:
            sel = (self.ampflag_table["shotid"] == self.shotid) * (
                self.ampflag_table["multiframe"] == self.multiframe
            )
            flag = self.ampflag_table["flag"][sel][0]
        except Exception:
            with self.bottombox:
                print("Could not find amp in amp flag table")
            flag = True

        if flag:
            box_layout = Layout()
            self.midbox.layout = box_layout
        else:
            box_layout = Layout(
                display="flex",
                # flex_flow='column',
                # align_items='stretch',
                border="5px solid red",
            )
            self.midbox.layout = box_layout

        try:
            self.im = get_image2D_amp(
                self.shotid_widget.value,
                multiframe=self.multiframe_widget.value,
                imtype=self.imtype_widget.value,
                expnum=self.expnum_widget.value,
                survey=self.survey,
            )
        
            self.imw.load_array(self.im)
        except:
            with self.bottombox:
                print('Could not open amp image')
            
        #if self.coords is not None:
            #self.imw.center_on(self.coords)
            #self.imw.zoom_level = 4

    def on_det_go(self, b):
        self.bottombox.clear_output()
        self.detectid = self.detectbox.value
        self.get_amp_info_from_det()
        
    def get_amp_info_from_det(self):

        global CONFIG_HDR2, CONFIG_HDR3, OPEN_DET_FILE, DET_HANDLE

        if (self.detectid >= 2100000000) * (self.detectid < 2190000000):
            self.det_file = CONFIG_HDR2.detecth5
        elif (self.detectid >= 2100000000) * (self.detectid < 2190000000):
            self.det_file = CONFIG_HDR2.contsourceh5
        elif (self.detectid >= 3000000000) * (self.detectid < 3090000000):
            self.det_file = CONFIG_HDR3.detecth5
        elif (self.detectid >= 3090000000) * (self.detectid < 3100000000):
            self.det_file = CONFIG_HDR3.contsourceh5
        elif (self.detectid >= 4000000000) * (self.detectid < 4090000000):
            self.det_file = CONFIG_HDR4.detecth5
        elif (self.detectid >= 4090000000) * (self.detectid < 4100000000):
            self.det_file = CONFIG_HDR4.contsourceh5
        elif (self.detectid >= 5000000000) * (self.detectid < 5090000000):
            self.det_file = CONFIG_HDR5.detecth5
        elif (self.detectid >= 5090000000) * (self.detectid < 5100000000):
            self.det_file = CONFIG_HDR5.contsourceh5

        if OPEN_DET_FILE is None:

            OPEN_DET_FILE = self.det_file
            DET_HANDLE = tb.open_file(self.det_file, 'r')

        else:
            if self.det_file == OPEN_DET_FILE:
                pass
            else:
                DET_HANDLE.close()
                OPEN_DET_FILE = self.det_file
                try:
                    DET_HANDLE = tb.open_file(self.det_file, 'r')
                except Exception:
                    with self.bottombox:
                        print("Could not open {}".format(self.det_file))

        detectid_obj = self.detectid

        try:
            det_row = DET_HANDLE.root.Detections.read_where("detectid == detectid_obj")[0]
            self.im_ra.value = det_row["ra"]
            self.im_dec.value = det_row["dec"]
            
            self.wave = det_row["wave"]
            self.wave_widget.value = self.wave
            self.coords = SkyCoord(
                det_row["ra"] * u.deg, det_row["dec"] * u.deg, frame="icrs"
            )
            
            if self.shotid != det_row["shotid"]:
                self.shotid = det_row['shotid']
                
            self.shotid_widget.value = self.shotid
            
            # get MF array for shot

            sel_shot = self.ampflag_table["shotid"] == self.shotid
            mflist = np.unique(self.ampflag_table["multiframe"][sel_shot])
            self.multiframe_widget.options = mflist
            self.multiframe = det_row["multiframe"].decode()
            self.multiframe_widget.value = self.multiframe
            self.expnum = det_row["expnum"]
            self.expnum_widget.value = self.expnum
            
            # update amp image
            self.update_amp_image()
            
            x = det_row["x_raw"]
            y = det_row["y_raw"]
            
            self.imw.marker = {"color": "red", "radius": 10, "type": "circle"}
            self.imw.add_markers(Table([[x - 1], [y - 1]], names=["x", "y"]))
        except IndexError:
            with self.bottombox:
                print('Detectid:{} is not found in database'.format(detectid_obj))
                
    def coord_change(self, b):
        self.shotid = None

        self.coords = SkyCoord(
            ra=self.im_ra.value * u.deg, dec=self.im_dec.value * u.deg
        )
        
        fiber_table_region = FIBINDEX.query_region(
            self.coords, radius=self.radius * u.arcsec, shotid=self.shotid
        )
        shotlist = np.unique(fiber_table_region["shotid"])
        if np.size(shotlist) > 0:
            self.shotid_widget = widgets.Dropdown(options=shotlist, value=shotlist[0])
        else:
            with self.bottombox:
                print('No observations near {}'.format(self.coords))
            
    def show_region(self, b):
        # get closest fiber:

        self.wave = self.wave_widget.value

        fibers = Fibers(self.shotid_widget.value)

        self.coords = SkyCoord(
            ra=self.im_ra.value * u.deg, dec=self.im_dec.value * u.deg
        )
        idx = fibers.get_closest_fiber(self.coords)
        multiframe_obj = fibers.table.cols.multiframe[idx].astype(str)
        self.multiframe = multiframe_obj
        self.multiframe_widget.value = self.multiframe

        expnum_obj = fibers.table.cols.expnum[idx]
        x, y = fibers.get_image_xy(idx, self.wave)

        self.imw.marker = {"color": "green", "radius": 10, "type": "circle"}
        self.imw.add_markers(Table([[x - 1], [y - 1]], names=["x", "y"]))

    def update_hdr_config(self):
        self.survey = self.survey_widget.value.lower()
        self.hetdex_api_config = HDRconfig( survey=self.survey)
        
        self.FibIndex = FiberIndex(self.survey)
        self.ampflag_table = Table.read(self.hetdex_api_config.badamp)
        # update survey class and shot list
   
        self.survey_class = Survey(self.survey)

            
    def get_ra_dec_wave(self, b):
        with self.bottombox:
            print('This is not functioning yet')
#        self.reset_markers()
#        self.imw.start_marking(marker={'color': 'red',
#                                       'radius': 5,
#                                       'type': 'cross'},
#                               marker_name='clicked markers',
#        )
#        while is_marking:
#            
#        self.imw.stop_marking()
