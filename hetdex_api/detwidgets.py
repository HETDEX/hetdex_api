from IPython.display import display
import ipywidgets as widgets
import numpy as np
import pickle


def on_button_clicked(b):
    print("Query limits saved to query.pickle")
    limits.wave_low = wave_low.value
    limits.wave_high = wave_high.value
    limits.flux_low = flux_low.value
    limits.flux_high = flux_high.value
    limits.linewidth_low = linewidth_low.value
    limits.linewidth_high = linewidth_high.value
    limits.sn_low = sn_low.value
    limits.sn_high = sn_high.value
    limits.chi2_low = chi2_low.value
    limits.chi2_high = chi2_high.value
    limits.cont_low = cont_low.value
    limits.cont_high = cont_high.value
    limits.aperture_flag = check.value
    limits.ra = ra.value
    limits.dec = dec.value
    limits.rad = rad.value
    limits.field = field.value
    outfile = open('query.pickle', 'wb')
    pickle.dump(limits, outfile)
    outfile.close()


class Det_limits:
    # limits to be returned to query detections
    def __init__(self):
        '''
        Initialize the Det_limits class
        '''
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

limits = Det_limits()

# set up button definitions and layout
but_layout = widgets.Layout(width='80%', height='40px')
text_layout = widgets.Layout(width='100%', height='35px')
start = widgets.Label(value='start', layout=text_layout)
end = widgets.Label(value='end', layout=text_layout)
blank = widgets.Label(value='',layout=text_layout)
wave_low = widgets.FloatText(value=None, layout=but_layout)
wave_high = widgets.FloatText(value=None, layout=but_layout)
flux_low = widgets.FloatText(value=None, layout=but_layout)
flux_high = widgets.FloatText(value=None, layout=but_layout)
linewidth_low = widgets.FloatText(value=None, layout=but_layout)
linewidth_high = widgets.FloatText(value=None, layout=but_layout)
sn_low = widgets.FloatText(value=None, layout=but_layout)
sn_high = widgets.FloatText(value=None, layout=but_layout)
chi2_low = widgets.FloatText(value=None, layout=but_layout)
chi2_high = widgets.FloatText(value=None, layout=but_layout)
cont_low = widgets.FloatText(value=None, layout=but_layout)
cont_high = widgets.FloatText(value=None, layout=but_layout)

# set up labels
wave_label = widgets.Label(value='wavelength (AA)')
flux_label = widgets.Label(value='lineflux ('+ r'\(10^{-17}\) ergs/s/cm'+r'\(^2\))')
linewidth_label = widgets.Label(value='linewidth (AA)')
sn_label = widgets.Label(value='S/N')
chi2_label = widgets.Label(value=r'\(\chi^2\)')
cont_label = widgets.Label(value='Continuum')        

# Set up GUI columns
textcol = widgets.VBox([blank, start, end])
col1 = widgets.VBox([wave_label, wave_low, wave_high])
col2 = widgets.VBox([flux_label, flux_low, flux_high])
col3 = widgets.VBox([linewidth_label, linewidth_low, linewidth_high])
col4 = widgets.VBox([sn_label, sn_low, sn_high])
col5 = widgets.VBox([chi2_label, chi2_low, chi2_high])
col6 = widgets.VBox([cont_label, cont_low, cont_high])

toggles = widgets.HBox([textcol,col1, col2, col3, col4, col5, col6])
button = widgets.Button(description="Select Detections", button_style='success')

# Set up Field Selection Toggles
field = widgets.SelectMultiple(
    options=['all','dex-spring','dex-fall','cosmos','egs','goods-n','other'],
    value=['all'],
    rows=7,
    description='Field:',
    disabled=False
    )
# set up coord vals
ra_label = widgets.Label(value='RA (degrees)', layout=text_layout)
dec_label = widgets.Label(value='Dec (degrees)', layout=text_layout)
rad_label = widgets.Label(value='Radius (arcmin)', layout=text_layout)
ra = widgets.FloatText(value=None, layout=but_layout)
dec = widgets.FloatText(value=None, layout=but_layout)
rad = widgets.FloatText(value=None,layout=but_layout)

print('Either select a specific field or multiple fields (command-click to select multiple fields):\n')
display(field)

print("Or manually enter an RA, DEC and aperture:\n")
labelcol = widgets.VBox([ra_label, dec_label, rad_label])
col = widgets.VBox([ra, dec, rad])
aperture = widgets.HBox([labelcol, col])
display(aperture)

check = widgets.Checkbox(
    value=False,
    description='Check to Select from aperture',
    disabled=False
    )
display(check)

print("You may also wish to consider some down selections:")
display(toggles, button)
button.on_click(on_button_clicked)


