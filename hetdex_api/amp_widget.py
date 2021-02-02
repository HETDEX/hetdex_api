A#!/usr/bin/env python
# coding: utf-8
#
# AmpWidget Widget Object
# Author: Erin Mentuch Cooper
# Date: 2020-11-11

import numpy as np

from astropy.table import Table
from astrowidgets import ImageWidget
import ipywidgets as widgets
from ipywidgets import Layout

from hetdex_api.config import HDRconfig
from hetdex_api.shot import get_image2D_amp, open_shot_file
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
        radius=None,
        detectid=None,
        wave=None,
        shotid=20200422015,
        multiframe="multi_406_071_081_LU",
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

        if self.detectid is not None:
            # open up detections h5 file to get info for
            # highest weight fiber
            self.get_amp_info_from_det()

    
        self.detectbox = widgets.BoundedIntText(
            value=self.detectid,
            min=2100000000,
            max=3000000000,
            step=1,
            description="DetectID:",
            disabled=False,
            # layout=Layout(width='15%')
        )
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

        self.shotid_widget = widgets.Dropdown(
            description="ShotID", options=shotlist, value=self.shotid,
        )

        self.shoth5 = open_shot_file(self.shotid, survey=self.survey)
        sel_shot = AMPFLAG_TABLE["shotid"] == self.shotid
        mflist = np.unique(AMPFLAG_TABLE["multiframe"][sel_shot])

        self.multiframe_widget = widgets.Dropdown(
            description="MultiframeID", options=mflist, value=mflist[0]
        )

        self.expnum_widget = widgets.Dropdown(
            description="ExpNum", options=[1, 2, 3], value=1,
        )

        self.imtype_widget = widgets.Dropdown(
            description="ImType",
            options=["clean_image", "image", "error"],
            value=self.imtype,
        )

        if self.coords is not None:
            self.im_ra = widgets.FloatText(
                value=self.coords.ra.value,
                description="RA (deg):",
                layout=Layout(width="20%"),
            )
            self.im_dec = widgets.FloatText(
                value=self.coords.dec.value,
                description="DEC (deg):",
                layout=Layout(width="20%"),
            )
        else:
            self.im_ra = widgets.FloatText(
                value=0.0, description="RA (deg):", layout=Layout(width="20%")
            )
            self.im_dec = widgets.FloatText(
                value=0.0, description="DEC (deg):", layout=Layout(width="20%")
            )

        self.select_coords = widgets.Button(
            description="Select coords", disabled=False, button_style="success"
        )

        if self.wave is None:
            wave0 = 3500
        else:
            wave0 = self.wave

        self.wave_widget = widgets.FloatText(
            value=wave0,
            description="Wavelength (A)",
            min=3500.0,
            max=5500.0,
            step=1.0)
            
        self.im = get_image2D_amp(
            self.shotid_widget.value,
            multiframe=self.multiframe,
            imtype=self.imtype,
            expnum=self.expnum,
        )

        self.imw.load_array(self.im)

        self.topbox = widgets.HBox(
            [
                self.survey_widget,
                self.detectbox,
                self.im_ra,
                self.im_dec,
                self.select_coords,
                self.wave_widget,
            ]
        )
        box_layout = Layout()

        self.boxside = widgets.VBox(
            [  # self.month_widget,
                self.shotid_widget,
                self.multiframe_widget,
                self.expnum_widget,
                self.imtype_widget,
            ]
        )
        self.midbox = widgets.HBox([self.imw, self.boxside], layout=box_layout)
        self.bottombox = widgets.Output(layout={"border": "1px solid black"})
        display(widgets.VBox([self.topbox, self.midbox, self.bottombox]))

        # plot region for detection
        self.multiframe_widget.observe(self.im_widget_change)
        self.expnum_widget.observe(self.im_widget_change)
        self.imtype_widget.observe(self.im_widget_change)
        self.shotid_widget.observe(self.shotid_widget_change)

    def im_widget_change(self, b):
        self.update_amp_image()

    def shotid_widget_change(self, b):
        self.shoth5.close()
        self.coords = None
        self.detectid = None

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
        except:
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
            
    def on_det_change(self, b):            
        self.bottombox.clear_output()
        self.detectid = self.detectbox.value
        self.update_det_coords()
        # get info from highest weight fiber
        
        self.im_ra.value = self.coords.ra.value
        self.im_dec.value = self.coords.dec.value

    def update_det_coords(self):
        detectid_i = self.detectid
        det_row = self.fileh5dets.root.Detections.read_where("detectid == detectid_i")
        if np.size(det_row) > 0:
            self.coords = SkyCoord(det_row["ra"][0] * u.deg, det_row["dec"][0] * u.deg)
        else:
            with self.bottombox:
                print(
                    "{} is not in the {} detect database".format(
                        detectid_i, self.survey
                    )
                )
                
    def get_amp_info_from_det():
            
        if self.detectid >= 2190000000:
            if CONT_H5_HANDLE is None:
                try:
                    CONT_H5_HANDLE = tb.open_file(CONT_H5_FN, "r")
                except:
                    print("Could not open continuum database")
            det_handle = CONT_H5_HANDLE
                    
        elif self.detectid >= 2100000000:
            if HETDEX_DETECT_HDF5_HANDLE is None:
                try:
                    HETDEX_DETECT_HDF5_HANDLE = tb.open_file(HETDEX_DETECT_HDF5_FN, "r")
                except:
                    print("Could not open detections database")
            det_handle = HETDEX_DETECT_HDF5_HANDLE

        detectid_obj = self.detecitd
        det_row = det_handle.root.Detections.read_where(
            "detectid == detectid_obj"
        )
        self.im_ra.value = det_row['ra']
        self.im_dec.value = det_row['dec']
        self.wave = det_row['wave']
        self.wave_widget.value = self.wave
        self.coords = SkyCoord(det_row["ra"], det_row["dec"], frame="icrs")
        self.multiframe = det_row["multiframe"].decode()
        self.multiframe_widget.value = self.multiframe
        self.shotid = det_row["shotid"]
        self.shotid_widget.value = self.shotid
        self.expnum = det_row['expnum']
        self.expnum_widget.value = self.shotid

        #update amp image
        self.update_amp_image()

        x = det_row['x_row']
        y = det_row['y_row']

        self.imw.marker = {'color': 'red', 'radius': 5, 'type': 'circle'}
        self.imw_add_markers(Table([x-1,y-1], names=['x','y']))
