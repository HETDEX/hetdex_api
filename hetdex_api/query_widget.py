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
from astropy.wcs.utils import skycoord_to_pixel
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import NDData
from astropy.table import Table

from astrowidgets import ImageWidget
import ipywidgets as widgets
from ipywidgets import interact, Layout, AppLayout

import plotly.graph_objects as go

from hetdex_api.shot import *
from hetdex_api.config import HDRconfig
from hetdex_tools.get_spec import get_spectra

from astroquery.sdss import SDSS
from elixer import catalogs

try:  # using HDRconfig
    LATEST_HDR_NAME = HDRconfig.LATEST_HDR_NAME
    CONFIG_HDR2 = HDRconfig('hdr2.1')
    CONFIG_HDR3 = HDRconfig('hdr3')
    
except Exception as e:
    print("Warning! Cannot find or import HDRconfig from hetdex_api!!", e)
    LATEST_HDR_NAME = "hdr2.1"

OPEN_DET_FILE = None
DET_HANDLE = None

class QueryWidget:
    def __init__(
        self,
        coords=None,
        detectid=None,
        survey=LATEST_HDR_NAME,
        aperture=3.0 * u.arcsec,
        cutout_size=5.0 * u.arcmin,
        zoom=3,
    ):

        self.survey = survey.lower()

        self.detectid = detectid
        self.aperture = aperture
        self.cutout_size = cutout_size
        self.zoom = zoom

        config = HDRconfig(survey=survey)

        self.catlib = catalogs.CatalogLibrary()

        if coords:
            self.coords = coords
            self.detectid = 1000000000
        elif detectid:
            self.detectid = detectid
            self.update_det_coords()
        else:
            self.coords = SkyCoord(191.663132 * u.deg, 50.712696 * u.deg, frame="icrs")
            self.detectid = 2101848640

        # initialize the image widget from astrowidgets
        self.imw = ImageWidget(image_width=600, image_height=600)

        self.survey_widget = widgets.Dropdown(
            options=["HDR1", "HDR2", "HDR2.1", "HDR3"],
            value=self.survey.upper(),
            layout=Layout(width="10%"),
        )

        self.detectbox = widgets.BoundedIntText(
            value=self.detectid,
            min=1000000000,
            max=4000000000,
            step=1,
            description="DetectID:",
            disabled=False,
        )
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

        self.pan_to_coords = widgets.Button(
            description="Pan to coords", disabled=False, button_style="success"
        )
        self.marking_button = widgets.Button(
            description="Mark Sources", button_style="success"
        )
        self.reset_marking_button = widgets.Button(
            description="Reset", button_style="success"
        )
        self.extract_button = widgets.Button(
            description="Extract Object", button_style="success"
        )

        self.marker_table_output = widgets.Output(layout={"border": "1px solid black"})
        #        self.spec_output = widgets.Output(layout={'border': '1px solid black'})

        self.spec_output = widgets.Tab(
            description="Extracted Spectra:", layout={"border": "1px solid black"}
        )
        self.textimpath = widgets.Text(
            description="Source: ", value="", layout=Layout(width="90%")
        )

        self.topbox = widgets.HBox(
            [
                self.survey_widget,
                self.detectbox,
                self.im_ra,
                self.im_dec,
                self.pan_to_coords,
            ]
        )
        self.leftbox = widgets.VBox(
            [self.imw, self.textimpath], layout=Layout(width="800px")
        )
        self.rightbox = widgets.VBox(
            [
                widgets.HBox(
                    [
                        self.marking_button,
                        self.reset_marking_button,
                        self.extract_button,
                    ]
                ),
                self.marker_table_output,
                self.spec_output,
            ],
            layout=Layout(width="800px"),
        )

        self.bottombox = widgets.Output(layout={"border": "1px solid black"})

        self.load_image()

        self.all_box = widgets.VBox([self.topbox,
                                     widgets.HBox([self.leftbox, self.rightbox]),
                                     #self.spec_output,
                                     self.bottombox])
        display(self.all_box)
        self.detectbox.observe(self.on_det_change)
        self.pan_to_coords.on_click(self.pan_to_coords_click)
        self.marking_button.on_click(self.marking_on_click)
        self.reset_marking_button.on_click(self.reset_marking_on_click)
        self.extract_button.on_click(self.extract_on_click)
        self.survey_widget.observe(self.on_survey_change)

    def on_survey_change(self, b):
        self.survey = self.survey_widget.value.lower()

    def update_coords(self):
        self.coords = SkyCoord(
            self.im_ra.value * u.deg, self.im_dec.value * u.deg, frame="icrs"
        )

    def on_det_change(self, b):
        self.bottombox.clear_output()
        self.detectid = self.detectbox.value
        self.update_det_coords()
        self.im_ra.value = self.coords.ra.value
        self.im_dec.value = self.coords.dec.value

    def update_det_coords(self):
        detectid_i = self.detectid
        
        global CONFIG_HDR2, CONFIG_HDR3, OPEN_DET_FILE, DET_HANDLE

        if (self.detectid >= 2100000000) * (self.detectid < 2190000000):
            self.det_file = CONFIG_HDR2.detecth5
        elif (self.detectid >= 2100000000) * (self.detectid < 2190000000):
            self.det_file = CONFIG_HDR2.contsourceh5
        elif (self.detectid >= 3000000000) * (self.detectid < 3090000000):
            self.det_file = CONFIG_HDR3.detecth5
        elif (self.detectid >= 3090000000) * (self.detectid < 3100000000):
            self.det_file = CONFIG_HDR3.contsourceh5

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
        try:
            det_row = DET_HANDLE.root.Detections.read_where("detectid == detectid_i")

            if np.size(det_row) > 0:
                self.coords = SkyCoord(det_row["ra"][0] * u.deg, det_row["dec"][0] * u.deg)
            else:
                with self.bottombox:
                    print(
                        "{} is not in the {} detect database".format(
                            detectid_i, self.survey
                        )
                    )
        except:
            with self.bottombox:
                print(
                    "{} is not in the {} detect database".format(
                        detectid_i, self.survey
                    )
                )
                
    def pan_to_coords_click(self, b):
        self.bottombox.clear_output()
        self.update_coords()

        x, y = skycoord_to_pixel(self.coords, self.cutout["cutout"].wcs)
        if (x > 0) and (y > 0):
            try:
                value = self.cutout["cutout"].data[int(x), int(y)]
                self.imw.center_on(self.coords)
            except:
                self.load_image()
        else:
            self.load_image()

    def load_image(self):

        im_size = self.cutout_size.to(u.arcsec).value
        mag_aperture = self.aperture.to(u.arcsec).value

        # keep original coords of image for bounds checking later
        self.orig_coords = self.coords

        with self.bottombox:
            try:
                self.cutout = self.catlib.get_cutouts(
                    position=self.coords,
                    side=im_size,
                    aperture=mag_aperture,
                    dynamic=False,
                    filter=["r", "g", "f606W"],
                    first=True,
                )[0]

                im = NDData(self.cutout["cutout"].data, wcs=self.cutout["cutout"].wcs)
                self.im_path = self.cutout["path"]
                self.imw.load_nddata(im)

            except:
                try:
                    sdss_im = SDSS.get_images(coordinates=self.coords, band="g")
                    im = sdss_im[0][0]
                except:
                    sdss_im = SDSS.get_images(
                        coordinates=self.coords, band="g", radius=30.0 * u.arcsec
                    )
                    im = sdss_im[0][0]

                self.im_path = "SDSS Astroquery result"
                self.imw.load_fits(im)

            self.imw.center_on(self.coords)
            self.imw.zoom_level = self.zoom
            self.textimpath.value = self.im_path

    def marking_on_click(self, b):

        if self.marking_button.button_style == "success":
            self.marker_table_output.clear_output()
            self.imw.reset_markers()
            self.imw.start_marking(
                marker={"color": "red", "radius": 3, "type": "circle"},
                marker_name="clicked markers",
            )
            self.marking_button.description = "Stop Marking"
            self.marking_button.button_style = "danger"

        else:
            self.imw.stop_marking()
            self.marking_button.description = "Mark Sources"
            self.marking_button.button_style = "success"

            self.marker_tab = self.imw.get_markers(marker_name="clicked markers")

            with self.marker_table_output:
                print(
                    "{:^5s} {:^8s} {:^8s} {:^28s}".format(
                        "OBJECT", "X", "Y", "Coordinates"
                    )
                )

                for index, row in enumerate(self.marker_tab):
                    c = row["coord"].to_string("hmsdms")

                    print(
                        "{:^5s} {:8.2f} {:8.2f} {}".format(
                            str(index + 1), row["x"], row["y"], c
                        )
                    )

    def reset_marking_on_click(self, b):

        self.marking_button.button_style = "success"
        self.marking_button.description = "Mark Sources"
        self.marker_table_output.clear_output()
        self.imw.reset_markers()
        self.spec_output.children = []
        self.bottombox.clear_output()

    def extract_on_click(self, b):
        self.bottombox.clear_output()

        with self.bottombox:
            self.spec_table = get_spectra(self.marker_tab["coord"], survey=self.survey)

        # set up tabs for plotting
        ID_list = np.unique(self.spec_table["ID"])

        self.out_widgets = []

        for i, obj in enumerate(ID_list):
            self.out_widgets.append(widgets.Output(description=obj))

        self.spec_output.children = self.out_widgets

        for i, obj in enumerate(ID_list):
            self.spec_output.set_title(i, "Obj=" + str(obj))
            with self.out_widgets[i]:
                self.plot_spec(obj)

    def plot_spec(self, obj):

        selobj = self.spec_table["ID"] == obj
        fig = go.FigureWidget()
        for row in self.spec_table[selobj]:
            x = row["wavelength"]
            y = row["spec"]
            if row['flag'] == 1:
                fig.add_trace(
                    go.Scatter(x=x, y=y, mode="lines",
                               name=str(row["shotid"]))
                )
            else:
                fig.add_trace(
                    go.Scatter(x=x, y=y, mode="lines",
                            opacity=0.25,
                            name=str(row["shotid"]) + ' (poor data)'))

        fig.update_traces(hoverinfo="text+name", mode="lines")
        fig.update_layout(
#            title="Object {}".format(row["ID"]),
            xaxis_title="wavelength (A)",
            yaxis_title="f_lambda (1e-17 ergs/s/cm^2/A)",
        )
        fig.update_layout(legend=dict(x=0.95, y=0.99, xanchor="right", yanchor="top"),
                          margin=dict(l=5, r=5, t=30, b=5),
                          autosize=True,
                          height=200,
                          width=500,
                          font_size=8, xaxis_range=[3540, 5510])
        fig.show()
        # fig, ax = plt.subplots(figsize=(8,2))
        # ax.plot(row['wavelength'], row['spec'])
        # ax.set_title('Object ' + str(row['ID']) + '       SHOTID = ' + str(row['shotid']))
        # ax.set_xlabel('wavelength (A)')
        # ax.set_ylabel('spec (1e-17 ergs/s/cm^2/A)')
        # plt.show(fig)
