#!/usr/bin/env python
# coding: utf-8
"""

Widget to view data cubes. This is specifically designed
To work with data cubes created by hetdex_toos.interpolate.make_data_cube

Author: Erin Mentuch Cooper

Date: November 23, 2020


Examples
--------
>>> hdu = make_data_cube(detectid=2100000335, subcont=False, dwave=10, imsize=40.*u.arcsec, pixscale=0.5*u.arcsec)#, dwave=10 )
>>> w = CubeWidget(hdu=hdu)
>>> help(make_data_cube)
"""


import numpy as np
from astropy import units as u
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import NDData
from astropy.visualization import ZScaleInterval
from astropy.coordinates import SkyCoord

from matplotlib import pyplot as plt

plt.style.use("fivethirtyeight")

from ginga.AstroImage import AstroImage
from astrowidgets import ImageWidget
import ipywidgets as widgets
from ipywidgets import interact, Layout, AppLayout
from IPython.display import display, clear_output

from hetdex_tools.interpolate import make_data_cube


class CubeWidget(ImageWidget):
    def __init__(self, hdu=None, im=None, wcs=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # self._4d_idx = 0  # Lock 4th dim to this for now

        if hdu is not None:
            self.im = hdu.data
            self.wcs = WCS(hdu.header)
        elif im is not None:
            self.im = im
            self.wcs = wcs
        else:
            print("Provide a 3D HDU or image and wcs object")

        self.nddata = NDData(self.im, wcs=self.wcs)
        self.load_nddata(self.nddata, n=0)

        # get wave info:

        self.dwave = self.wcs.wcs.cdelt[2]
        self.wave_start = self.wcs.wcs.crval[2]
        self.nwave = np.shape(self.im)[0]
        self.wave_end = self.wave_start + self.nwave * self.dwave

        # self.set_colormap('Greys')

        zscale = ZScaleInterval(contrast=0.3, krej=2.5)

        vmin, vmax = zscale.get_limits(values=self.im)

        self.cuts = (vmin, vmax)

        self.wave_widget = widgets.IntSlider(
            min=self.wave_start,
            max=self.wave_end,
            step=self.dwave,
            value=self.wave_start,
            continuous_update=False,
        )

        self.slider = widgets.interactive(self.show_slice, wave=self.wave_widget)

        self.animate_button = widgets.Button(
            description="Scan Cube",
            disabled=False,
            button_style="success",
            tooltip="Click this to scan in wavelength dimension",
        )

        # For line profile plot
        self._cur_islice = None
        self._cur_ix = None
        self._cur_iy = None
        self.line_out = widgets.Output()
        self.line_plot = None
        self.plot_xlabel = "Wavelength (A)"
        self.plot_ylabel = "Spec (ergs/s/cm^2 -per spaxel)"

        # If plot shows, rerun cell to hide it.
        ax = plt.gca()
        self.line_plot = ax

        self.scan = widgets.Play(
            value=self.wave_start,
            min=self.wave_start,
            max=self.wave_end,
            step=self.dwave,
            #            interval=500,
            description="Scan Cube",
            disabled=False,
        )

        widgets.jslink((self.scan, "value"), (self.wave_widget, "value"))

        left_panel = widgets.VBox([widgets.HBox([self.wave_widget, self.scan]), self])

        display(widgets.HBox([left_panel, self.line_out]))

    def load_nddata(self, nddata, n=0):  # update this for wavelength later

        self.image = AstroImage()
        self.image.load_nddata(nddata, naxispath=[n])
        self._viewer.set_image(self.image)

    def _mouse_click_cb(self, viewer, event, data_x, data_y):

        self._cur_ix = int(round(data_x))
        self._cur_iy = int(round(data_y))
        self.plot_line_profile()
        # Ensure only active marker is shown
        self.reset_markers()

        if self._cur_ix is not None:
            mrk_tab = Table(names=["x", "y"])
            mrk_tab.add_row([self._cur_ix, self._cur_iy])
            self.marker = {"color": "red", "radius": 1, "type": "circle"}
            self.add_markers(mrk_tab)
        # self.reset_markers()

        super()._mouse_click_cb(viewer, event, data_x, data_y)

    def plot_line_profile(self):
        if self.line_plot is None or self._cur_ix is None or self._cur_iy is None:
            return

        #        image = self._viewer.get_image()
        if self.image is None:
            return

        with self.line_out:
            mddata = self.image.get_mddata()
            self.line_plot.clear()

            self.line_plot.plot(
                (self.wave_start + self.dwave * np.arange(self.nwave)),
                mddata[:, self._cur_iy, self._cur_ix],
                "b-",
                linewidth=1.2,
            )

            if self._cur_islice is not None:
                y = mddata[self._cur_islice, self._cur_iy, self._cur_ix]
                x = self.wave_start + self.dwave * self._cur_islice
                self.line_plot.axvline(x=x, color="r", linewidth=1)

            # self.line_plot.set_title(f'X={self._cur_ix + 1} Y={self._cur_iy + 1}')
            self.line_plot.set_xlabel(self.plot_xlabel)
            self.line_plot.set_ylabel(self.plot_ylabel)
            self.line_plot.set_xlim(self.wave_start, self.wave_end)

            clear_output(wait=True)
            display(self.line_plot.figure)

    def image_show_slice(self, n):
        # image = self._viewer.get_image()
        self.image.set_naxispath([n])
        self._viewer.redraw(whence=0)
        self._cur_islice = n

    def show_slice(self, wave):

        n = int((wave - self.wave_start) / self.dwave)

        self.image_show_slice(n - 1)
        self.plot_line_profile()
