from __future__ import print_function



"""
ELiXer Widget Classify Supriod Detections

Based on Dr. Erin Mentuch Cooper's original elixer_widgets.py


"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import os.path as op

from astropy.io import ascii
from astropy.table import Table, Column

import ipywidgets as widgets
#from IPython import display#, HTML
from IPython.display import Image
from ipywidgets import interact, Layout #Style #, interactive
#from IPython.display import clear_output
from hetdex_api.detections import *
import tables

#needed only if detection observered wavelength is not supplied
HETDEX_DETECT_HDF5_FN = "/work/03946/hetdex/hdr1/detect/detect_hdr1.h5"
HETDEX_DETECT_HDF5_HANDLE = None

elix_dir_archive = '/work/05350/ecooper/stampede2/elixer/jpgs/'
elix_dir = '/work/03261/polonius/hdr1_classify/all_pngs/'
# set up classification dictionary and associated widget
# the widget takes an optional detection list as input either
# as an array of detectids or a text file that can be loaded in
# You may also initiate with no variables to just open any ELiXeR
# on demand



current_wavelength = -1.0

line_id_dict_lae = "1216 LyA"
line_id_dict_lae_wave = 1216.0
line_id_dict_sep = "--------------"
line_id_dict_other = "Other (xxxx)"

line_id_dict_default = "Unk (????)"
line_id_dict = {line_id_dict_default:-1.0,
                line_id_dict_lae:1216.0,
                "3727 OII":3727.0,
                "4863 H-beta": 4862.68,
                "4960 OIII":4960.295,
                "5008 OIII":5008.240,
                line_id_dict_sep:-1.0,
                "1241 NV":1240.81,
                "1260 SiII": 1260.0,
                "1549 CIV":1549.48,
                "1640 HeII": 1640.4,
                "1909 CII":1908.734,
                "2799 MgII":2799.117,
                "4102 H-delta": 4102.0,
                "4342 H-gamma": 4341.68,
                line_id_dict_other:-1.0
                }

class ElixerWidget():


    def __init__(self, detectfile=None, detectlist=None, savedfile=None, outfile=None, resume=False, img_dir=None):

        global elix_dir

        self.current_idx = 0

        if img_dir is not None:
            if op.exists(img_dir):
                elix_dir = img_dir

        if detectfile:
            self.detectid = np.loadtxt(detectfile, dtype=np.int32,ndmin=1)
            self.vis_class = np.zeros(np.size(self.detectid), dtype=int)
            self.flag = np.zeros(np.size(self.detectid),dtype=int)
            self.z = np.full(np.size(self.detectid),-1.0)
            self.comment = np.full(np.size(self.detectid), '?',dtype='|S80').astype(str)
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

                #could have z
                try:
                    self.z = np.array(saved_data['z'],dtype=float)
                except:
                    self.z = np.full(np.size(self.detectid), -1.0)

                #could hve comment
                try:
                    self.comment = np.array(saved_data['comments'], dtype='|S80').astype(str)
                except:
                    self.comment = np.full(np.size(self.detectid), '?',dtype='|S80').astype(str)

            except:
                print("Could not open and read in savedfile. Are you sure its in astropy table format")


       #now, just a list of detections
        elif type(detectlist) is np.ndarray:
            self.detectid = detectlist
            self.vis_class = np.zeros(np.size(self.detectid), dtype=int)
            self.flag = np.zeros(np.size(self.detectid), dtype=int)
            self.z = np.full(np.size(self.detectid), -1.0)
            self.comment = np.full(np.size(self.detectid), '?',dtype='|S80').astype(str)

        else:
            self.detectid = np.arange(1000000000, 1000690799, 1)
            self.vis_class = np.zeros(np.size(self.detectid), dtype=int)
            self.flag = np.zeros(np.size(self.detectid), dtype=int)
            self.z = np.full(np.size(self.detectid), -1.0)
            self.comment = np.full(np.size(self.detectid), '?', dtype='|S80').astype(str)


        # store outfile name if given
        if outfile:
            print(
                "Careful with this option, it likely won't work properly. You are better off using the savedfile option")
            self.outfilename = outfile
        elif savedfile:
            self.outfilename = savedfile
        else:
            self.outfilename = 'elixer_lae.dat'

        self.resume = resume
        self.setup_widget()

        interact(self.main_display, x=self.detectbox)


    def main_display(self, x):

        detectid = x
        show_selection_buttons = True

        try:
            objnum = np.where(self.detectid == detectid)[0][0]
            print('On ELiXer Report '+ str(objnum+1) + '/' + str(np.size(self.detectid)))
        except:
            print('Current object not in original list. Go to Next or Previous DetectID to return to input Detectlist')
            show_selection_buttons = False

        self.previousbutton.on_click(self.on_previous_click)
        self.nextbutton.on_click(self.on_next_click)
        self.elixerNeighborhood.on_click(self.on_elixer_neighborhood)

        if show_selection_buttons:
            # clear_output()
            self.rest_widget_values(objnum)
            display(widgets.HBox([self.previousbutton, self.nextbutton, self.elixerNeighborhood,
                                  self.line_id_drop, self.wave_box, self.z_box,self.comment_box]))

            display(widgets.HBox([self.sm1_button,self.s0_button,self.s1_button,self.s2_button,self.s3_button,
                                  self.s4_button,self.s5_button]))

            self.sm1_button.on_click(self.sm1_button_click)
            self.s0_button.on_click(self.s0_button_click)
            self.s1_button.on_click(self.s1_button_click)
            self.s2_button.on_click(self.s2_button_click)
            self.s3_button.on_click(self.s3_button_click)
            self.s4_button.on_click(self.s4_button_click)
            self.s5_button.on_click(self.s5_button_click)
        else:
            display(widgets.HBox([self.previousbutton, self.nextbutton, self.elixerNeighborhood]))

        try:
            fname = op.join(elix_dir, "%d.png" % (detectid))

            if op.exists(fname):
                display(Image(fname))
            else: #try the archive location
                print("Cannot load ELiXer Report image: ", fname)
                print("Trying archive location...")
                fname = op.join(elix_dir_archive, "egs_%d" % (detectid // 100000), str(detectid) + '.jpg')
                if op.exists(fname):
                    display(Image(fname))
                else:
                    print("Cannot load ELiXer Report image: ", fname)
        except:
            print("Cannot load ELiXer Report image: ", fname)

    def setup_widget(self):
        if self.resume:
            try:
                i_start = np.max(np.where(self.flag != 0)) + 1

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
        self.previousbutton = widgets.Button(layout=Layout(width='5%'))#description='Previous DetectID')
        self.nextbutton = widgets.Button(layout=Layout(width='5%'))#description='Next DetectID')

        #see https://fontawesome.com/icons?d=gallery
        self.previousbutton.style.button_color = 'darkgray'
        self.previousbutton.icon = 'arrow-circle-left'
        self.nextbutton.style.button_color = 'darkgray'
        self.nextbutton.icon = 'arrow-circle-right'
        #self.nextbutton.layout = Layout()

        self.elixerNeighborhood = widgets.Button(description='Neighbors',layoout=Layout(width="10%"))#, button_style='info')
        self.elixerNeighborhood.style.button_color='darkgray'
        #self.detectwidget = widgets.HBox([self.detectbox, self.nextbutton])

        self.line_id_drop = widgets.Dropdown(options=line_id_dict.keys(),
                                             value=line_id_dict_default,description='Line',
                                             layout=Layout(width="20%"),
                                             disabled=False)
        self.line_id_drop.observe(self._handle_line_id_selection, names='value')

        self.wave_box = widgets.FloatText(value=-1.0,step=0.00001,description=r"$\lambda$ rest",
                                          layout=Layout(width="17%"),disabled=False,indent=False)
        # self.wave_box = widgets.Text(
        #     value='-1.0',
        #     description=r"$\lambda$ rest",
        #     layout=Layout(width="20%")
        #     disabled=False),


        self.wave_box.observe(self._handle_wave_box_selection,names='value')
        self.z_box = widgets.FloatText(value=-1.0,step=0.00001,description="z",
                                       layout=Layout(width="17%"),disabled=False,indent=False)
        self.z_box.add_class("left-spacing-class")
        #display(HTML("<style>.left-spacing-class {margin-left: 10px;}</style>"))

        self.comment_box = widgets.Text(
            value='',
            placeholder='Enter any comments here',
            description='Comments:',
            disabled=False)#,layout=Layout(width='20%'))

        #buttons as classification selection
        # self.s0_button = widgets.Button(description=' No Imaging ', button_style='success')
        # self.s1_button = widgets.Button(description=' Spurious ', button_style='success')
        # self.s2_button = widgets.Button(description=' Not LAE ', button_style='success')
        # self.s3_button = widgets.Button(description=' Unknown ', button_style='success')
        # self.s4_button = widgets.Button(description=' Maybe LAE ', button_style='success')
        # self.s5_button = widgets.Button(description=' Definite LAE ', button_style='success')

        self.sm1_button = widgets.Button(description='Spurious',
                                         button_style='danger',
                                         layout=Layout(width='10%'))

        self.s0_button = widgets.Button(description='  Not LAE (0) ', button_style='success')
        self.s1_button = widgets.Button(description='          (1) ', button_style='success')
        self.s2_button = widgets.Button(description='          (2) ', button_style='success')
        self.s3_button = widgets.Button(description='          (3) ', button_style='success')
        self.s4_button = widgets.Button(description='          (4) ', button_style='success')
        self.s5_button = widgets.Button(description=' YES! LAE (5) ', button_style='success')



        #self.submitbutton = widgets.Button(description="Submit Classification", button_style='success')
        #self.savebutton = widgets.Button(description="Save Progress", button_style='success')


    def get_observed_wavelength(self):
        global HETDEX_DETECT_HDF5_HANDLE,HETDEX_DETECT_HDF5_FN,current_wavelength

        current_wavelength = -1.0

        if HETDEX_DETECT_HDF5_HANDLE is None:
            HETDEX_DETECT_HDF5_HANDLE = tables.open_file(HETDEX_DETECT_HDF5_FN)

        dtb = HETDEX_DETECT_HDF5_HANDLE.root.Detections
        q_detectid = self.detectbox.value
        rows = dtb.read_where('detectid==q_detectid')
        if (rows is not None) and (rows.size == 1):
            current_wavelength = rows[0]['wave']

        return current_wavelength

    def _handle_wave_box_selection(self,event):

        global current_wavelength, line_id_dict_lae_wave

        _old_wave = event['old']
        _new_wave = event['new']

        if _old_wave == _new_wave:
            return

        if _new_wave < 0:  # i.e. a -1.00
            # do not change the z and lambda_rest values
            self.z_box.value = -1.0
            self.wave_box.value = -1.0
        else:
            if current_wavelength < 0:  # need to find it
                self.get_observed_wavelength()

            self.z_box.value = np.round(current_wavelength / _new_wave - 1.0, 5)
            self.wave_box.value = np.round(_new_wave, 5)

            #do NOT reset the line id label ... can cause a loop
            #self.get_line_match(self.z_box.value, current_wavelength )



    def _handle_line_id_selection(self,event):
        """
        event': 'value',
        'old': 'Unk (????)',
        'new': '4960 OIII',
        'owner': Dropdown(description='Line ID', index=13, options=('Unk (????)', 'Other (xxxx)', '1216 LyA', '1241 NV', '1260 SiII', '1549 CIV', '1640 HeII', '1909 CII', '2799 MgII', '3727 OII', '4102 H-delta', '4342 H-gamma', '4863 H-beta', '4960 OIII', '5008 OIII'), value='4960 OIII'),
        'type': 'change'
        :param event:
        :return:
        """
        global current_wavelength,line_id_dict_lae_wave

        _old = event['old']
        _old_wave = line_id_dict[_old]
        _new = event['new']
        _new_wave = line_id_dict[_new]

        if _old_wave == _new_wave:
            return

        if _new_wave < 0: #i.e. a -1.00
            #do not change the z and lambda_rest values
            if _new == line_id_dict_default:
                self.z_box.value = -1.0
                self.wave_box.value = -1.0
        else:
            if current_wavelength < 0: #need to find it
                self.get_observed_wavelength()

            self.z_box.value = np.round(current_wavelength / _new_wave - 1.0, 5)
            self.wave_box.value = np.round(_new_wave,5)

            # if (self.z_box.value < 0.0) and (self.line_id_drop.value != line_id_dict_default):
            #     self.z_box.style = 'warning'


            #what button value to hit?
            #LAE lines would be LyA, CIV, HeII
            #set these to '5' or '4' or what?
            #This causes some syncrhonization problems ... probably should skip it
            # if abs(self.wave_box.value -line_id_dict_lae_wave) < 0.1:
            #     if self.s5_button.icon == '': #only trip the button if not already at 5
            #         self.s5_button_click(None)





        #causes dirty flag
        # self.current_idx = ix
        # self.detectbox.value = self.detectid[ix]

        # print("Drop down event",event)
        # print(event['old'],event['new'],line_id_dict[event['new']])


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

        #print("Prev detect idx", ix)

        self.rest_widget_values(idx=ix)

        #print("Prev detect post reset idx", )

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

        #make sure the values match
        if self.line_id_drop.value == line_id_dict_default:
            self.wave_box.value = -1.0
            self.z_box.value = -1.0

        self.z[self.current_idx] =  self.z_box.value
        self.comment[self.current_idx] = self.comment_box.value

        self.on_save_click(None)
        self.goto_next_detect()


    def get_line_match(self,z,w):
        global current_wavelength

        if z > 0:
            self.line_id_drop.value = line_id_dict_other
            self.wave_box.value = np.round(w / (1. + z),5)
        else:
            self.line_id_drop.value = line_id_dict_default
            self.wave_box.value = -1.0

        if z is None or z < 0 or w is None or w < 0:
            return  self.line_id_drop.value, self.wave_box.value

        w_rest = np.round(w / (1. + z),5)
        #find the match in the line_id_dict
        for k in line_id_dict.keys():
            if abs(line_id_dict[k] - w_rest) < 1.0:
                self.line_id_drop.value = k
                self.wave_box.value = w_rest
                break


        return self.line_id_drop.value, self.wave_box.value



    def rest_widget_values(self,idx=0):

        global current_wavelength

        #print("Reset idx",idx,"Current w",current_wavelength)


        if self.detectbox.value < 1000000000: #assume an index
            self.detectbox.value = self.detectid[idx]
            return


        self.z_box.value = self.z[idx] #-1.0
        current_wavelength = self.get_observed_wavelength()
        self.line_id_drop.value,self.wave_box.value = self.get_line_match(self.z_box.value,current_wavelength)

        self.comment_box.value = self.comment[idx]

        #print("Updated Reset idx", idx, "Current w", current_wavelength)


        # reset all to base
        self.sm1_button.icon = ''
        self.s0_button.icon = ''
        self.s1_button.icon = ''
        self.s2_button.icon = ''
        self.s3_button.icon = ''
        self.s4_button.icon = ''
        self.s5_button.icon = ''

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
            elif self.vis_class[idx] == -1:
                self.sm1_button.icon = 'check'



    def on_previous_click(self, b):
        self.goto_previous_detect()

    def on_next_click(self, b):
        self.goto_next_detect()

    def sm1_button_click(self, b):
        global line_id_dict_lae
        self.line_id_drop.value == line_id_dict_default
        self.wave_box.value = -1.0
        self.z_box.value = -1.0

        self.set_classification(-1)
    
    def s0_button_click(self, b):
        global line_id_dict_lae
        if self.is_consistent_with_lae():
            #NOT okay, if you say is not LAE
            self.line_id_drop.value == line_id_dict_default
            self.wave_box.value = -1.0
            self.z_box.value = -1.0
        else:
            pass #okay, already NOT consistent with LAE

        self.set_classification(0)

    def s1_button_click(self, b):
        global line_id_dict_lae
        if self.is_consistent_with_lae():
            #NOT okay, if you say is not LAE
            self.line_id_drop.value == line_id_dict_default
            self.wave_box.value = -1.0
            self.z_box.value = -1.0
        else:
            pass #okay, already NOT consistent with LAE

        self.set_classification(1)

    def s2_button_click(self, b):
        global line_id_dict_lae

        if self.is_consistent_with_lae():
            #NOT okay, if you say is not LAE
            self.line_id_drop.value == line_id_dict_default
            self.wave_box.value = -1.0
            self.z_box.value = -1.0
        else:
            pass #okay, already NOT consistent with LAE

        self.set_classification(2)

    def s3_button_click(self, b):
        #leave the drop box where ever it is
        self.set_classification(3)

    def s4_button_click(self, b):
        global line_id_dict_lae
        if self.is_consistent_with_lae():
            pass  # already okay, could be CIV, etc
        else: #not consistent with LAE, so reset
            self.line_id_drop.value == line_id_dict_default
            self.wave_box.value = -1.0
            self.z_box.value = -1.0

        self.set_classification(4)

    def s5_button_click(self, b):
        global line_id_dict_lae

        if self.is_consistent_with_lae():
            pass #already okay, could be CIV, etc
        else: #not consistent with LAE, so reset
            self.line_id_drop.value == line_id_dict_default
            self.wave_box.value = -1.0
            self.z_box.value = -1.0

        self.set_classification(5)


    def is_consistent_with_lae(self):
        #is the emission line consistent with this being an LAE in HETDEX
        #basically, a cheat ... is  1.9 < z < 3.5
        if 1.8 < self.z_box.value < 3.6:
            return True
        else:
            return False


    def on_save_click(self, b):
        self.output = Table()
        self.output.add_column(Column(self.detectid, name='detectid', dtype=int))
        self.output.add_column(Column(self.vis_class, name='vis_class', dtype=int))
        self.output.add_column(Column(self.flag, name='flag', dtype=int))
        self.output.add_column(Column(self.z,name='z',dtype=float))
        self.output.add_column(Column(self.comment,name='comments',dtype='|S80'))

        ascii.write(self.output, self.outfilename, overwrite=True)

    def on_elixer_neighborhood(self,b):
        detectid = self.detectbox.value
        path = os.path.join(elix_dir, "%dnei.png" % (detectid))

        if not os.path.isfile(path):
            print("%s not found" % path)
        else:
            display(Image(path))
