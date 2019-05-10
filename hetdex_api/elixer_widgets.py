import sys
import numpy as np
import matplotlib.pyplot as plt
import os.path as op

from astropy.io import ascii
from astropy.table import Table, Column

import ipywidgets as widgets
from IPython.display import Image
from ipywidgets import interact, interactive
from hetdex_api.detections import *

# set up classification dictionary and associated widget 
# the widget takes an optional detection list as input either
# as an array of detectids or a text file that can be loaded in
# You may also initiate with no variables to just open any ELiXeR 
# on demand

class ElixerWidget():
    '''
    intialize the widget class

    Optional Inputs

    detectlist = numpy array of DetectIDs, should be Integer type
                 Often this will be some subselection for objects
                 from the Detections() class object.
    
    detectfile = path to a file of detectIDs. Should have been 
                 saved as np.savetxt(detects.detectid[sel], detectfile)
    
    None       = if no list is given it will load in all detectIDs 
                 from HDR1. There will be no finessing so beward
                 many ELiXers will be junky.

    savedfile  = astropy table saved to ascii file containing columns
                 of detectid, vis_class and comments

    outfile    = name of file to store classifcations, will use elixer_classifications.dat
                 if not provided. Not needed if savedfile is provided.
    
    For now we are using the following integer flags for vis_class. We also store
    the classification type as a string in vis_type

        # assign a vis_class field for future classification
        # -2 = ignore (bad detectid, software issue, would not want this for testing)
        # -1 = no assignemnt
        # 0 = artifact (pixel artifact/bad amp/cosmic)
        # 1 = OII emitter
        # 2 = LAE emitter
        # 3 = star
        # 4 = nearby galaxies (HBeta, OIII usually)
        # 5 = other line (CIV in a LAE for example)
        # 6 = noise - low sn detection, prob not real
        # 7 = unknown - weird object/hard to tell/flag to follow up
    
    '''

    def __init__(self, detectfile=None, detectlist=None, savedfile=None, outfile=None, resume=False):

        if detectfile:
            self.detectid = np.loadtxt(detectfile, dtype=np.int32)
            self.vis_type = np.zeros(np.size(self.detectid), dtype='|S15')
            self.vis_class = -1*np.ones(np.size(self.detectid), dtype=int)
            self.comment = np.zeros(np.size(self.detectid), dtype='|S30')
        elif savedfile:
            try: 
                saved_data = ascii.read(savedfile)
                self.detectid = np.array(saved_data['detectid'], dtype=int)
                self.vis_type = np.array(saved_data['vis_type'], dtype='|S15')
                self.vis_class = np.array(saved_data['vis_class'], dtype=int)
                self.comment = np.array(saved_data['comments'], dtype='|S30')
            except:
                print("Could not open and read in savedfile. Are you sure its in astropy table format")
        elif type(detectlist) is np.ndarray:
            self.detectid = detectlist
            self.vis_type = np.zeros(np.size(self.detectid), dtype='|S15')
            self.vis_class = -1*np.ones(np.size(self.detectid), dtype=int)
            self.comment = np.zeros(np.size(self.detectid), dtype='|S30')
        else:
            self.detectid = np.arange(1000000000,1000690799,1)
            self.vis_type = np.zeros(np.size(self.detectid), dtype='|S15')
            self.vis_class = -1*np.ones(np.size(self.detectid), dtype=int)
            self.comment = np.zeros(np.size(self.detectid), dtype='|S30')

        # store outfile name if given
        if outfile:
            self.outfilename = outfile
        elif savedfile:
            self.outfilename = savedfile
        else:
            self.outfilename = 'elixer_classifications.dat'

        self.resume=resume

        self.setup_widget()

        interact(self.main_display, x=self.detectbox)

    def add(self, detectid_i, vis_type_i, comment_i):
        ix = np.where(self.detectid == detectid_i)
        self.vis_type[ix] = vis_type_i
        self.comment[ix] = comment_i
        vis_dict={'OII Galaxy' : 1,
                  'LAE Galaxy' : 2, 
                  'Star' :3, 
                  'Nearby Galaxy': 4,
                  'Other Line': 5,
                  'Artifact' : 0,
                  'Noise' : 6, 
                  'Unknown': 7}

        self.vis_class[ix] = vis_dict.get(vis_type_i)

    def main_display(self, x):
        detectid = x
        elix_dir = '/work/05350/ecooper/stampede2/elixer/jpgs/'
        file_jpg = op.join(elix_dir, "egs_%d" %(detectid//100000), str(detectid) + '.jpg')
        display(widgets.HBox([self.previousbutton,self.nextbutton]))
        display(Image(file_jpg))
        display(self.classification)
        display(self.comments)
        self.previousbutton.on_click(self.on_previous_click)
        self.nextbutton.on_click(self.on_next_click)
        display(widgets.HBox([self.submitbutton, self.savebutton]))
        self.submitbutton.on_click(self.on_button_click)
        self.savebutton.on_click(self.on_save_click)
        

    def setup_widget(self):
        if self.resume:
            detectstart = self.detectid[np.max(np.where(self.vis_class != -1)) + 1] 
        else:
            detectstart = np.min(self.detectid)
        self.detectbox = widgets.BoundedIntText(
            value=detectstart,
            min=1000000000,
            max=1000690799,
            step=1,
            description='DetectID:',
            disabled=False
        )
        self.previousbutton = widgets.Button(description='Previous DetectID', button_style='success')
        self.nextbutton = widgets.Button(description='Next DetectID', button_style='success')
        self.detectwidget = widgets.HBox([self.detectbox, self.nextbutton])
        self.classification = widgets.ToggleButtons(
            options=['OII Galaxy', 'LAE Galaxy', 'Star', 'Nearby Galaxy', 'Other Line', 'Artifact', 'Noise', 'Unknown'],
            description='Type:',
            disabled=False,
            button_style='', # 'success', 'info', 'warning', 'danger' or ''
            tooltips=['Low-z OII[3727] emitters', 
                      'Distant high-z LAE[1216] emitters',
                      'Object is a star, spectrum has no emission',
                      'Nearby galaxy, likely Hbeta[4861] or OIII[5007] emitter',
                      'Choose this if you are sure its CIV or some other emission',
                      'Detector Artifact/cosmic ray',
                      'Junky - low SN, very hard to tell, software issue, not suitable for machine learning']
        )
        self.comments = widgets.Text(
            value='',
            placeholder='Enter any comments here',
            description='Comments:',
            disabled=False
        )
        self.submitbutton = widgets.Button(description="Submit Classification", button_style='success')
        self.savebutton = widgets.Button(description="Save Progress", button_style='success')


    def goto_previous_detect(self):
        ix = np.max(np.where(self.detectid <= self.detectbox.value))
        if ix-1 >= 0:
            self.detectbox.value = self.detectid[ix-1]
        else:
            print("At the beginning of the DetectID List")
        # clear the comment box                                                                 
        self.comments.placeholder='Enter any comments here'
        self.comments.value = ''

    def goto_next_detect(self):
        ix = np.where(self.detectid >= self.detectbox.value)[0][0]
        if ix+1 < np.size(self.detectid):
            self.detectbox.value = self.detectid[ix+1]
        else:
            print("At the end of the DetectID List")
        # clear the comment box
        self.comments.placeholder='Enter any comments here'
        self.comments.value = ''
        
    def on_button_click(self, b):
        self.add(self.detectbox.value, self.classification.value, self.comments.value)
        self.goto_next_detect()

    def on_previous_click(self, b):
        self.goto_previous_detect()

    def on_next_click(self, b):
        self.goto_next_detect()

    def on_save_click(self, b):
        self.output = Table()
        self.output.add_column(Column(self.detectid, name='detectid', dtype=int))
        self.output.add_column(Column(self.vis_type, name='vis_type', dtype='|S15'))
        self.output.add_column(Column(self.vis_class, name='vis_class', dtype=int))
        self.output.add_column(Column(self.comment, name='comments', dtype='|S30'))
        ascii.write(self.output, self.outfilename, overwrite=True)
