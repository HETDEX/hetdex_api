import numpy as np
import matplotlib.pyplot as plt
import os.path as op
import ipywidgets as widgets
from IPython.display import Image
from ipywidgets import interact, interactive


# set up classificaiton dictionary
class class_dict:
    def __init__(self, detectlist):
        self.detectid = detectlist
#        self.vis_class = -1*np.ones(np.size(detectlist))
        self.vis_class = np.chararray(np.size(detectlist), 15)
        self.comment = np.chararray(np.size(detectlist), 60)
    def add(self, detectid_i, vis_class_i, comment_i):
        ix = np.where(self.detectid == detectid_i)
        self.vis_class[ix] = vis_class_i
        self.comment[ix] = comment_i



def main(x):
    detectid = detectbox.value
    elix_dir = '/scratch/03261/polonius/jpgs'
    file_jpg = op.join(elix_dir, "egs_%d" %(detectid//100000), str(detectid) + '.jpg')
    display(Image(file_jpg))
    display(classification)
    display(notes)
    display(submitbutton)
    submitbutton.on_click(on_button_click)

def on_button_click(b):
    class_dict.add(detectbox.value, classification.value, notes.value)
    detectbox.value +=1


detectbox = widgets.BoundedIntText(
    value=7,
    min=1000000000,
    max=1000690798,
    step=1,
    description='DetectID:',
    disabled=False
)

classification = widgets.ToggleButtons(
    options=['OII Galaxy', 'LAE Galaxy', 'Star', 'Nearby Galaxy', 'Other Line', 'Artifact'],
    description='Type:',
    disabled=False,
    button_style='', # 'success', 'info', 'warning', 'danger' or ''
    tooltips=['Description of slow', 'Description of regular', 'Description of fast'],
#     icons=['check'] * 3
)
notes = widgets.Text(
    value='',
    placeholder='Enter any comments here',
    description='Comments:',
    disabled=False
)
submitbutton = widgets.Button(description="Select Detections", button_style='success')

interact(main, x=detectbox)
