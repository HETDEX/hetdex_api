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
    HETDEX_API_CONFIG = HDRconfig(survey=LATEST_HDR_NAME)
    HDR_BASEPATH = HETDEX_API_CONFIG.hdr_dir[LATEST_HDR_NAME]
    HETDEX_DETECT_HDF5_FN = HETDEX_API_CONFIG.detecth5
    HETDEX_DETECT_HDF5_HANDLE = None
    CONT_H5_FN = HETDEX_API_CONFIG.contsourceh5
    CONT_H5_HANDLE = None
    FIBINDEX = FiberIndex(LATEST_HDR_NAME)
    AMPFLAG_TABLE = Table.read(HETDEX_API_CONFIG.badamp)

except Exception as e:
    print("Warning! Cannot find or import HDRconfig from hetdex_api!!", e)
    LATEST_HDR_NAME = "hdr2.1"


class AmpWidget:
    def __init__(
        self,
        survey=LATEST_HDR_NAME,
        coords=None,
        radius=3.0,
        detectid=
        wave=None,
        shotid=None,
        multiframe=None,
        imtype="clean_image",
        expnum=1,
    ):

        global LATEST_HDR_NAME, HETDEX_DETECT_HDF5_FN, HETDEX_DETECT_HDF5_HANDLE
        global AMPFLAG_TABLE

        self.survey = survey.lower()
        self.detectid = detectid
        self.shotid = shotid
        self.multiframe = multiframe
        self.imtype = imtype
        self.expnum = expnum
        self.coords = coords
        self.radius = radius
        self.wave = wave

        s = Survey(self.survey)
        gs = s.remove_shots()

        self.survey_class = s[gs]

        # initialize the image widget from astrowidgets
        self.imw = ImageWidget()  # image_width=600, image_height=600)

        self.survey_widget = widgets.Dropdown(
            options=["HDR1", "HDR2.1"],
            value=self.survey.upper(),
            layout=Layout(width="10%"),
        )

        #        monthlist = np.unique(np.array(s[gs].date/100).astype(int))
        shotlist = np.unique(np.array(s[gs].shotid))

        #        self.month_widget = widgets.Dropdown(
        #            options=monthlist,
        #            value=202001,
        #            description='Month',
        #            disabled=False,
        #            layout=Layout(width='20%'),
        #        )

        self.detectbox = widgets.BoundedIntText(
            value=self.detectid,
            min=2100000000,
            max=3000000000,
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
            description="ShotID", options=shotlist, value=self.shotid,
        )

        if self.shotid is not None:
            self.shoth5 = open_shot_file(self.shotid, survey=self.survey)
            sel_shot = AMPFLAG_TABLE["shotid"] == self.shotid
            mflist = np.unique(AMPFLAG_TABLE["multiframe"][sel_shot])
            self.multiframe_widget = widgets.Dropdown(
                description="MultiframeID", options=mflist, value=mflist[0]
            )
        else:
            mflist = np.unique(AMPFLAG_TABLE["multiframe"])
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
            ]
        )
        self.midbox = widgets.HBox([self.imw, self.boxside], layout=box_layout)

        self.bottombox = widgets.Output(layout={"border": "1px solid black"})
        
        if self.coords is not None:
            if self.detectid is not None:
                
                self.shotid = None
                
                fiber_table_region = FIBINDEX.query_region(
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
        )
        
        self.imw.load_array(self.im)

        display(widgets.VBox([self.topbox, self.midbox, self.bottombox]))

        # plot region for detection
        self.multiframe_widget.observe(self.im_widget_change)
        self.expnum_widget.observe(self.im_widget_change)
        self.imtype_widget.observe(self.im_widget_change)
        self.shotid_widget.observe(self.shotid_widget_change)
        self.select_coords.on_click(self.coord_change)
        self.show_button.on_click(self.show_region)
        self.det_button.on_click(self.on_det_go)

    def im_widget_change(self, b):
        self.update_amp_image()

    def shotid_widget_change(self, b):
        self.bottombox.clear_output()
        self.shoth5.close()
        self.coords = None
        self.detectid = None
        self.im_ra.value = 0.0
        self.im_dec.value = 0.0
        
        self.shotid = self.shotid_widget.value

        self.shoth5 = open_shot_file(self.shotid_widget.value, survey=self.survey)
        sel_shot = AMPFLAG_TABLE["shotid"] == self.shotid
        mflist = np.unique(AMPFLAG_TABLE["multiframe"][sel_shot])

        self.multiframe_widget.options = mflist
        self.multiframe_widget.value = self.multiframe

        try:
            sel = (AMPFLAG_TABLE["shotid"] == self.shotid) * (
                AMPFLAG_TABLE["multiframe"] == self.multiframe
            )
            flag = AMPFLAG_TABLE["flag"][sel][0]
        except:
            self.multiframe = mflist[0]
            sel = (AMPFLAG_TABLE["shotid"] == self.shotid) * (
                AMPFLAG_TABLE["multiframe"] == self.multiframe
            )
            flag = AMPFLAG_TABLE["flag"][sel][0]

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

        try:
            sel = (AMPFLAG_TABLE["shotid"] == self.shotid) * (
                AMPFLAG_TABLE["multiframe"] == self.multiframe
            )
            flag = AMPFLAG_TABLE["flag"][sel][0]
        except Exception:
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

        self.im = get_image2D_amp(
            self.shotid_widget.value,
            multiframe=self.multiframe_widget.value,
            imtype=self.imtype_widget.value,
            expnum=self.expnum_widget.value,
        )
        self.imw.load_array(self.im)

        if self.coords is not None:
            self.imw.center_on(self.coords)
            self.imw.zoom_level = 4

    def on_det_go(self, b):
        self.bottombox.clear_output()
        self.detectid = self.detectbox.value
        self.get_amp_info_from_det()
        
    def get_amp_info_from_det(self):

        global CONT_H5_HANDLE, HETDEX_DETECT_HDF5_HANDLE
        global CONT_H5_FN, HETDEX_DETECT_HDF5_FN

        if self.detectid >= 2190000000:
            if CONT_H5_HANDLE is None:
                try:
                    CONT_H5_HANDLE = tb.open_file(CONT_H5_FN, "r")
                except Exception:
                    print("Could not open continuum database")
            det_handle = CONT_H5_HANDLE

        elif self.detectid >= 2100000000:
            if HETDEX_DETECT_HDF5_HANDLE is None:
                try:
                    HETDEX_DETECT_HDF5_HANDLE = tb.open_file(
                        HETDEX_DETECT_HDF5_FN, "r")
                except Exception:
                    print("Could not open detections database")
            det_handle = HETDEX_DETECT_HDF5_HANDLE

        detectid_obj = self.detectid

        if True:
            det_row = det_handle.root.Detections.read_where("detectid == detectid_obj")[0]
            self.im_ra.value = det_row["ra"]
            self.im_dec.value = det_row["dec"]
            self.wave = det_row["wave"]
            self.wave_widget.value = self.wave
            self.coords = SkyCoord(
                det_row["ra"] * u.deg, det_row["dec"] * u.deg, frame="icrs"
            )
            
            if self.shotid != det_row["shotid"]:
                self.shoth5.close()
                self.shotid = det_row['shotid']
                
            self.shotid_widget.value = self.shotid
            
            # get MF array for shot
            self.shoth5 = open_shot_file(self.shotid_widget.value, survey=self.survey)
            sel_shot = AMPFLAG_TABLE["shotid"] == self.shotid
            mflist = np.unique(AMPFLAG_TABLE["multiframe"][sel_shot])

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
        else:#except IndexError:
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

        self.shotid_widget = widgets.Dropdown(options=shotlist, value=shotlist[0])

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
