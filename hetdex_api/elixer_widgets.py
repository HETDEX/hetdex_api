import sys
import numpy as np
import matplotlib.pyplot as plt
import os.path as op
import pickle

import ipywidgets as widgets
from IPython.display import Image
from ipywidgets import interact, interactive
from hetdex_api.detections import *

# set up classification dictionary and associated widget 
# the widget takes an optional detection list as input either
# as an array of detectids or a pickle file that can be loaded in
# You may also initiate with no variables to just open any ELiXeR 
# on demand

class ElixerWidget():
    '''
    intialize the widget class

    Optional Inputs
    
    picklefile = path to a pickle file of detectIDs (since these
                 binary files perserve python datatypes they load
                 and save very quickly. We recommend using this 
                 method
    
    detectlist = numpy array of DetectIDs, should be Integer type

    None       = if no list is given it will load in all detectIDs 
                 from HDR1. There will be no finessing so beward
                 many ELiXers will be junky

    savedpickle = pickle file of previously saved work, will 
                 contain arrays of both detectids and vis_class

    '''

    def __init__(self, picklefile=None, detectlist=None, savedpickle=None):

        if picklefile:
            self.detectid = pickle.load( open(picklefile, "rb"))
            self.vis_class = np.zeros(np.size(self.detectid), dtype='|S15')
        elif detectlist:
            self.detectid = detectlist
            self.vis_class = np.zeros(np.size(self.detectid), dtype='|S15')
        elif savedpickle:
            saved_data = open(savedpickle, "rb")
            data = pickle.load(saved_data)
            self.detectid = data[0]
            self.vis_class = data[1]
            saved_data.close()
        else:
            self.detectid = np.arange(1000000000,1000690799,1)
            self.vis_class = np.chararray(np.size(self.detectid), 15)

        self.comment = np.zeros(np.size(self.detectid), dtype='|S30') 

        self.setup_widget()

        interact(self.main_display, x=self.detectbox)

    def add(self, detectid_i, vis_class_i, comment_i):
        ix = np.where(self.detectid == detectid_i)
        self.vis_class[ix] = vis_class_i
        self.comment[ix] = comment_i

    def main_display(self, x):
        detectid = x
        elix_dir = '/scratch/03261/polonius/jpgs'
        file_jpg = op.join(elix_dir, "egs_%d" %(detectid//100000), str(detectid) + '.jpg')
        display(self.nextbutton)
        display(Image(file_jpg))
        display(self.classification)
        display(self.comments)
        self.nextbutton.on_click(self.on_next_click)
        display(widgets.HBox([self.submitbutton, self.savebutton]))
        self.submitbutton.on_click(self.on_button_click)
        self.savebutton.on_click(self.on_save_click)

    def setup_widget(self):

        self.detectbox = widgets.BoundedIntText(
            value=np.min(self.detectid),
            min=np.min(self.detectid),
            max=np.max(self.detectid),
            step=1,
            description='DetectID:',
            disabled=False
        )
        self.nextbutton = widgets.Button(description='Next DetectID', button_style='success')
        self.detectwidget = widgets.HBox([self.detectbox, self.nextbutton])
        self.classification = widgets.ToggleButtons(
            options=['OII Galaxy', 'LAE Galaxy', 'Star', 'Nearby Galaxy', 'Other Line', 'Artifact','Junk'],
            description='Type:',
            disabled=False,
            button_style='', # 'success', 'info', 'warning', 'danger' or ''
            tooltips=['Low-z OII[3727] emitters', 
                      'Distant high-z LAE[1216] emitters',
                      'Object is a star, spectrum has no emission',
                      'Nearby galaxy, likely Hbeta[4861] or OIII[5007] emitter',
                      'Choose this if you are sure its CIV or some other emission',
                      'Detector Artifact/cosmic ray',
                      'Junky - low SN, very hard to tell']
        )
        self.comments = widgets.Text(
            value='',
            placeholder='Enter any comments here',
            description='Comments:',
            disabled=False
        )
        self.submitbutton = widgets.Button(description="Submit Classification", button_style='success')
        self.savebutton = widgets.Button(description="Save Progress", button_style='success')



    def goto_next_detect(self):
        ix = np.where(self.detectid >= self.detectbox.value)[0][0]
        if ix+1 < np.size(self.detectid):
            self.detectbox.value = self.detectid[ix+1]
        else:
            print("At the end of the self.detectid")
        # clear the comment box
        self.comments.placeholder='Enter any comments here'
        self.comments.value = ''
        
    def on_button_click(self, b):
        self.add(self.detectbox.value, self.classification.value, self.comments.value)
        self.goto_next_detect()

    def on_next_click(self, b):
        self.goto_next_detect()

    def on_save_click(self, b):
          outfile = open('elixer_widget.pickle', 'wb')
          pickle.dump([self.detectid, self.vis_class], outfile)
          outfile.close()
