from __future__ import print_function
"""
ELiXer Widget Classify Supriod Detections

Based on Dr. Erin Mentuch Cooper's original elixer_widgets.py

Classify the type of spurious detection


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
#from IPython.display import clear_output
from hetdex_api.detections import *

elix_dir = '/work/05350/ecooper/stampede2/elixer/jpgs/'
# set up classification dictionary and associated widget
# the widget takes an optional detection list as input either
# as an array of detectids or a text file that can be loaded in
# You may also initiate with no variables to just open any ELiXeR
# on demand


class ElixerWidget():


    def __init__(self, detectfile=None, detectlist=None, savedfile=None, outfile=None, resume=False):

        self.current_idx = 0

        if detectfile:
            self.detectid = np.loadtxt(detectfile, dtype=np.int32,ndmin=1)
            self.vis_class = np.zeros(np.size(self.detectid), dtype=int)
            self.flag = np.zeros(np.size(self.detectid),dtype=int)
                #hidden flag, distinguish vis_class 0 as unset vs reviewed & fake
                #and possible future use as followup

        elif savedfile:
            try:
                saved_data = ascii.read(savedfile)
                self.detectid = np.array(saved_data['detectid'], dtype=int)
                self.vis_class = np.array(saved_data['vis_class'], dtype=int)

                #could have a flag
                try:
                    self.flag = np.array(saved_data['flag'],dtype=int)
                except:
                    self.flag = np.zeros(np.size(self.detectid), dtype=int)

            except:
                print("Could not open and read in savedfile. Are you sure its in astropy table format")


       #now, just a list of detections
        elif type(detectlist) is np.ndarray:
            self.detectid = detectlist
            self.vis_class = np.zeros(np.size(self.detectid), dtype=int)
            self.flag = np.zeros(np.size(self.detectid), dtype=int)

        else:
            self.detectid = np.arange(1000000000, 1000690799, 1)
            self.vis_class = np.zeros(np.size(self.detectid), dtype=int)
            self.flag = np.zeros(np.size(self.detectid), dtype=int)


        # store outfile name if given
        if outfile:
            print(
                "Careful with this option, it likely won't work properly. You are better off using the savedfile option")
            self.outfilename = outfile
        elif savedfile:
            self.outfilename = savedfile
        else:
            self.outfilename = 'elixer_spur.dat'

        self.resume = resume
        self.setup_widget()

        interact(self.main_display, x=self.detectbox)

    def main_display(self, x):

        detectid = x
        show_selection_buttons = True

        try:
            objnum = np.where(self.detectid == detectid)[0][0]
            print('On ELiXer Report ' + str(objnum + 1) + '/' + str(np.size(self.detectid)))
        except:
            print('Current object not in original list. Go to Next or Previous DetectID to return to input Detectlist')
            show_selection_buttons = False


        display(widgets.HBox([self.previousbutton, self.nextbutton]))
        self.previousbutton.on_click(self.on_previous_click)
        self.nextbutton.on_click(self.on_next_click)

        if show_selection_buttons:
            # clear_output()
            self.rest_widget_values(objnum)
            display(widgets.HBox([self.s0_button, self.s1_button, self.s2_button, self.s3_button,
                                  self.s4_button, self.s5_button,self.s6_button]))
            self.s0_button.on_click(self.s0_button_click)
            self.s1_button.on_click(self.s1_button_click)
            self.s2_button.on_click(self.s2_button_click)
            self.s3_button.on_click(self.s3_button_click)
            self.s4_button.on_click(self.s4_button_click)
            self.s5_button.on_click(self.s5_button_click)
            self.s6_button.on_click(self.s6_button_click)

        try:
            fname = op.join(elix_dir, "egs_%d" % (detectid // 100000), str(detectid) + '.jpg')

            if op.exists(fname):
                display(Image(fname))
            else:
                print("Cannot load ELiXer Report image: ", fname)
        except:
            print("Cannot load ELiXer Report image: ", fname)


    def setup_widget(self):
        if self.resume:
            try:
                i_start = np.max(np.where(self.vis_class != 0)) + 1

                if i_start is None:
                    i_start = 0
                    detectstart = self.detectid[i_start]
                elif i_start < np.size(self.detectid):
                    detectstart = self.detectid[i_start]
                else:
                    i_start = 0
                    detectstart = self.detectid[i_start]  # np.min(self.detectid)
            except:
                i_start = 0
                detectstart = self.detectid[i_start]

        else:
            i_start = 0
            detectstart = self.detectid[i_start]  # np.min(self.detectid)

        self.current_idx = i_start

        self.detectbox = widgets.BoundedIntText(
            value=detectstart,
            #min=1,
            min=1000000000,
            max=1000690799,
            step=1,
            description='DetectID:',
            disabled=False
        )
        self.previousbutton = widgets.Button(description='Previous DetectID', button_style='success')
        self.nextbutton = widgets.Button(description='Next DetectID', button_style='success')
        self.detectwidget = widgets.HBox([self.detectbox, self.nextbutton])


        #buttons as classification selection
        self.s0_button = widgets.Button(description=' Sky Sub ', button_style='success')
        self.s1_button = widgets.Button(description=' Pixel Flat ', button_style='success')
        self.s2_button = widgets.Button(description=' Bad Amp ', button_style='success')
        self.s3_button = widgets.Button(description=' Cosmic ', button_style='success')
        self.s4_button = widgets.Button(description=' Noise ', button_style='success')
        self.s5_button = widgets.Button(description=' Other ', button_style='success')
        self.s6_button = widgets.Button(description=' REAL ', button_style='success')


        #self.submitbutton = widgets.Button(description="Submit Classification", button_style='success')
        #self.savebutton = widgets.Button(description="Save Progress", button_style='success')
    def goto_previous_detect(self):

        try:
            ix = np.where(self.detectid == self.detectbox.value)[0][0]

            if ix - 1 >= 0:
                ix -= 1
            else:
                print("At the beginning of the DetectID List")
                return

        except:
            #invalid index ... the report displayed is not in the operating list
            #so use the last good index:
            ix = self.current_idx

        self.rest_widget_values(idx=ix)

        #causes dirty flag
        self.current_idx = ix
        self.detectbox.value = self.detectid[ix]


    def goto_next_detect(self):
        try:
            ix = np.where(self.detectid == self.detectbox.value)[0][0]
            if ix+1 < np.size(self.detectid):
                ix += 1
            else:
                print("At the end of the DetectID List")
                return
        except:
            #invalid index ... the report displayed is not in the operating list
            #so use the last good index:
            ix = self.current_idx

        self.rest_widget_values(idx=ix)
        self.current_idx = ix
        self.detectbox.value = self.detectid[ix]



    def set_classification(self,value=0):
        self.current_idx = np.where(self.detectid == self.detectbox.value)[0][0]  # current position

        self.vis_class[self.current_idx] = value
        self.flag[self.current_idx] = 1

        self.on_save_click(None)
        self.goto_next_detect()


    def rest_widget_values(self,idx=0):

        if self.detectbox.value < 1000000000: #assume an indedx
            self.detectbox.value = self.detectid[idx]
            return

        # reset all to base
        self.s0_button.icon = ''
        self.s1_button.icon = ''
        self.s2_button.icon = ''
        self.s3_button.icon = ''
        self.s4_button.icon = ''
        self.s5_button.icon = ''
        self.s6_button.icon = ''

        # mark the label on the button
        if self.flag[idx] != 0:
            if self.vis_class[idx] == 0:
                self.s0_button.icon = 'check'
            elif self.vis_class[idx] == 1:
                self.s1_button.icon = 'check'
            elif self.vis_class[idx] == 2:
                self.s2_button.icon = 'check'
            elif self.vis_class[idx] == 3:
                self.s3_button.icon = 'check'
            elif self.vis_class[idx] == 4:
                self.s4_button.icon = 'check'
            elif self.vis_class[idx] == 5:
                self.s5_button.icon = 'check'
            elif self.vis_class[idx] == 6:
                self.s6_button.icon = 'check'


    def on_previous_click(self, b):
        self.goto_previous_detect()

    def on_next_click(self, b):
        self.goto_next_detect()


    def s0_button_click(self, b):
        self.set_classification(0)

    def s1_button_click(self, b):
        self.set_classification(1)

    def s2_button_click(self, b):
        self.set_classification(2)

    def s3_button_click(self, b):
        self.set_classification(3)

    def s4_button_click(self, b):
        self.set_classification(4)

    def s5_button_click(self, b):
        self.set_classification(5)

    def s6_button_click(self, b):
        self.set_classification(6)


    def on_save_click(self, b):
        self.output = Table()
        self.output.add_column(Column(self.detectid, name='detectid', dtype=int))
        self.output.add_column(Column(self.vis_class, name='vis_class', dtype=int))
        self.output.add_column(Column(self.flag, name='flag', dtype=int))

        ascii.write(self.output, self.outfilename, overwrite=True)
