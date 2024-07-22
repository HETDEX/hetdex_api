from __future__ import print_function
from IPython.display import display
import ipywidgets as widgets
import pickle


class Det_limits:
    # Dictionary of limits to be returned to query detections
    def __init__(self):
        """
        Initialize the Det_limits class
        """
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


class DetWidget:
    def __init__(self):

        self.limits = Det_limits()

        # set up button definitions and layout
        but_layout = widgets.Layout(width="80%", height="40px")
        text_layout = widgets.Layout(width="100%", height="35px")
        self.start = widgets.Label(value="low", layout=text_layout)
        self.end = widgets.Label(value="high", layout=text_layout)
        blank = widgets.Label(value="", layout=text_layout)
        self.wave_low = widgets.FloatText(value=None, layout=but_layout)
        self.wave_high = widgets.FloatText(value=None, layout=but_layout)
        self.flux_low = widgets.FloatText(value=None, layout=but_layout)
        self.flux_high = widgets.FloatText(value=None, layout=but_layout)
        self.linewidth_low = widgets.FloatText(value=None, layout=but_layout)
        self.linewidth_high = widgets.FloatText(value=None, layout=but_layout)
        self.sn_low = widgets.FloatText(value=None, layout=but_layout)
        self.sn_high = widgets.FloatText(value=None, layout=but_layout)
        self.chi2_low = widgets.FloatText(value=None, layout=but_layout)
        self.chi2_high = widgets.FloatText(value=None, layout=but_layout)
        self.cont_low = widgets.FloatText(value=None, layout=but_layout)
        self.cont_high = widgets.FloatText(value=None, layout=but_layout)

        # set up labels
        wave_label = widgets.Label(value="wavelength (AA)")
        flux_label = widgets.Label(
            value="lineflux (" + r"\(10^{-17}\) ergs/s/cm" + r"\(^2\))"
        )
        linewidth_label = widgets.Label(value="linewidth (AA)")
        sn_label = widgets.Label(value="S/N")
        chi2_label = widgets.Label(value=r"\(\chi^2\)")
        cont_label = widgets.Label(value="Continuum")

        # Set up GUI columns
        textcol = widgets.VBox([blank, self.start, self.end])
        col1 = widgets.VBox([wave_label, self.wave_low, self.wave_high])
        col2 = widgets.VBox([flux_label, self.flux_low, self.flux_high])
        col3 = widgets.VBox([linewidth_label,
                             self.linewidth_low,
                             self.linewidth_high])
        col4 = widgets.VBox([sn_label, self.sn_low, self.sn_high])
        col5 = widgets.VBox([chi2_label, self.chi2_low, self.chi2_high])
        col6 = widgets.VBox([cont_label, self.cont_low, self.cont_high])

        self.toggles = widgets.HBox([textcol, col1, col2,
                                     col3, col4, col5, col6])
        self.button = widgets.Button(
            description="Select Detections", button_style="success"
        )

        # Set up Field Selection Toggles
        self.field = widgets.SelectMultiple(
            options=[
                "all",
                "dex-spring",
                "dex-fall",
                "cosmos",
                "egs",
                "goods-n",
                "other",
            ],
            value=["all"],
            rows=7,
            description="Field:",
            disabled=False,
        )
        # set up coord vals
        ra_label = widgets.Label(value="RA (degrees)", layout=text_layout)
        dec_label = widgets.Label(value="Dec (degrees)", layout=text_layout)
        rad_label = widgets.Label(value="Radius (arcmin)", layout=text_layout)
        self.ra = widgets.FloatText(value=None, layout=but_layout)
        self.dec = widgets.FloatText(value=None, layout=but_layout)
        self.rad = widgets.FloatText(value=None, layout=but_layout)

        print("Either select a specific field or multiple fields (command-click to select multiple fields):\n")
        display(self.field)

        print("Or manually enter an RA, DEC and aperture:\n")
        labelcol = widgets.VBox([ra_label, dec_label, rad_label])
        col = widgets.VBox([self.ra, self.dec, self.rad])
        aperture = widgets.HBox([labelcol, col])
        display(aperture)

        self.check = widgets.Checkbox(
            value=False,
            description="Check to Select from aperture",
            disabled=False
        )
        display(self.check)

        print("You may also wish to consider some down selections:")
        display(self.toggles, self.button)
        self.button.on_click(self.on_button_clicked)

    def on_button_clicked(self, b):
        print("Query limits saved to query.pickle")
        self.limits.wave_low = self.wave_low.value
        self.limits.wave_high = self.wave_high.value
        self.limits.flux_low = self.flux_low.value
        self.limits.flux_high = self.flux_high.value
        self.limits.linewidth_low = self.linewidth_low.value
        self.limits.linewidth_high = self.linewidth_high.value
        self.limits.sn_low = self.sn_low.value
        self.limits.sn_high = self.sn_high.value
        self.limits.chi2_low = self.chi2_low.value
        self.limits.chi2_high = self.chi2_high.value
        self.limits.cont_low = self.cont_low.value
        self.limits.cont_high = self.cont_high.value
        self.limits.aperture_flag = self.check.value
        self.limits.ra = self.ra.value
        self.limits.dec = self.dec.value
        self.limits.rad = self.rad.value
        self.limits.field = self.field.value
        outfile = open("query.pickle", "wb")
        pickle.dump(self.limits, outfile)
        outfile.close()
