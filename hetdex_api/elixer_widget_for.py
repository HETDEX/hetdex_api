from __future__ import print_function
"""
Based on Dr. Erin Mentuch Cooper's original elixer_widgets.py

This is a simplified version that only presents a scaled binary selection for the user to mark his/her confidence
that the presented emission line detection is fake or real. Fake refers to some data artifact that does not correspond
to real (astrophysical) photons as an EMISSION LINE.  A real detection typically belongs to a galaxy (but could be an 
emission nebula or even a meteor).

Fake: 
* random noise in the spectrum
* sky line
* interference pattern
* cosmic ray strike
* hot CCD pixel
* adjacent (real) absorpotion features that leave a "peak" going back to the continuum level misterpreted as emission

Real:
* any emission line
* emission line on top of continuum (usually a nearby galaxy)
* emission line from a transient (meteor, etc)
"""

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

elix_dir = '/work/05350/ecooper/stampede2/elixer/jpgs/'
# set up classification dictionary and associated widget
# the widget takes an optional detection list as input either
# as an array of detectids or a text file that can be loaded in
# You may also initiate with no variables to just open any ELiXeR
# on demand

class ElixerWidget():


    def __init__(self, detectfile=None, detectlist=None, savedfile=None, outfile=None, resume=False):

        if detectfile:
            self.detectid = np.loadtxt(detectfile, dtype=np.int32)
            self.vis_class = np.zeros(np.size(self.detectid), dtype=int)

        elif savedfile:
            try:
                saved_data = ascii.read(savedfile)
                self.detectid = np.array(saved_data['detectid'], dtype=int)
                self.vis_class = np.array(saved_data['vis_class'], dtype=int)

            except:
                print("Could not open and read in savedfile. Are you sure its in astropy table format")
        elif type(detectlist) is np.ndarray:
            self.detectid = detectlist
            self.vis_class = np.zeros(np.size(self.detectid), dtype=int)

        else:
            self.detectid = np.arange(1000000000, 1000690799, 1)
            self.vis_class = np.zeros(np.size(self.detectid), dtype=int)


        # store outfile name if given
        if outfile:
            print(
                "Careful with this option, it likely won't work properly. You are better off using the savedfile option")
            self.outfilename = outfile
        elif savedfile:
            self.outfilename = savedfile
        else:
            self.outfilename = 'elixer_for.dat'

        self.resume = resume
        self.setup_widget()

        interact(self.main_display, x=self.detectbox)

    def add(self, detectid_i, classification):
        ix = np.where(self.detectid == detectid_i)
        vis_dict = {
            'Fake -5':-5,
            '-4':-4,
            '-3':-3,
            '-2':-2,
            '-1':-1,
            '0': 0,
            '1':1,
            '2':2,
            '3':3,
            '4':4,
            '5 Real':5
          }
        self.vis_class[ix] = vis_dict.get(classification)


    def main_display(self, x):
        detectid = x

        try:
            objnum = np.where(self.detectid == detectid)[0][0]
            print('On ELiXer Report '+ str(objnum+1) + '/' + str(np.size(self.detectid)))
        except: 
            print('Current object not in original list. Go to Next or Previous DetectID to return to input Detectlist')

        self.rest_widget_values(objnum)

        file_jpg = op.join(elix_dir, "egs_%d" %(detectid//100000), str(detectid) + '.jpg')


        display(widgets.HBox([self.previousbutton, self.nextbutton]))
        display(Image(file_jpg))
        display(widgets.HBox([self.classification,]))


        self.previousbutton.on_click(self.on_previous_click)
        self.nextbutton.on_click(self.on_next_click)
        display(widgets.HBox([self.submitbutton, self.savebutton]))
        self.submitbutton.on_click(self.on_button_click)
        self.savebutton.on_click(self.on_save_click)

    def setup_widget(self):
        if self.resume:
            i_start = np.max(np.where(self.vis_class != 0)) + 1
            if i_start < np.size(self.detectid):
                detectstart = self.detectid[i_start]
            else:
                detectstart = self.detectid[0] #np.min(self.detectid)
        else:
            i_start = 0
            detectstart = self.detectid[0]  # np.min(self.detectid)

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

        vis_dict = {
            -5: 'Fake -5',
            -4: '-4',
            -3: '-3',
            -2: '-2',
            -1: '-1',
            0: '0',
            1: '1',
            2: '2',
            3: '3',
            4: '4',
            5: '5 Real'
        }

        self.classification = widgets.ToggleButtons(
                options=['Fake -5','-4','-3','-2','-1','0','1','2','3','4','5 Real'],
                description='',
                value=vis_dict[self.vis_class[i_start]],
                disabled=False,
                button_style='',
                tooltips=['']
            )
        self.classification.style.button_width = 'initial'

        self.submitbutton = widgets.Button(description="Submit Classification", button_style='success')
        self.savebutton = widgets.Button(description="Save Progress", button_style='success')

    def goto_previous_detect(self):
        ix = np.where(self.detectid == self.detectbox.value)[0][0]

        if ix - 1 >= 0:
            ix -= 1
        else:
            print("At the beginning of the DetectID List")
            return

        self.rest_widget_values(idx=ix)

        #causes dirty flag
        self.detectbox.value = self.detectid[ix]

    def goto_next_detect(self):

        ix = np.where(self.detectid == self.detectbox.value)[0][0] #current position
        if ix+1 < np.size(self.detectid):
            ix += 1

        else:
            print("At the end of the DetectID List")
            return

        self.rest_widget_values(idx=ix)
        self.detectbox.value = self.detectid[ix]




    def rest_widget_values(self,idx=0):
        # very sloppy ... reverse of vis_dict later
        vis_dict = {
            -5:'Fake -5',
            -4:'-4',
            -3:'-3',
            -2:'-2',
            -1:'-1',
            0:'0',
            1:'1',
            2:'2',
            3:'3',
            4:'4',
            5:'5 Real'
        }

        self.classification.value = vis_dict[self.vis_class[idx]]



    def on_button_click(self, b):
        self.add(self.detectbox.value, self.classification.value)
        self.goto_next_detect()

    def on_previous_click(self, b):
        self.goto_previous_detect()

    def on_next_click(self, b):
        self.goto_next_detect()

    def on_save_click(self, b):
        self.output = Table()
        self.output.add_column(Column(self.detectid, name='detectid', dtype=int))
        self.output.add_column(Column(self.vis_class, name='vis_class', dtype=int))

        ascii.write(self.output, self.outfilename, overwrite=True)
