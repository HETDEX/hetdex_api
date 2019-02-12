# -*- coding: utf-8 -*-
"""
Created: 2019/01/23

@author: Erin Mentuch Cooper

"""

import glob
import sys
import os
import os.path as op
import subprocess
import numpy as np
import tables as tb
import h5py

from astropy.io import ascii
from astropy.io import fits

path_amp2amp = "/work/00115/gebhardt/maverick/getampnorm/all/"
path_throughput = "/work/00115/gebhardt/maverick/detect/tp/"

class AmpToAmp(tb.IsDescription):
    ampid = tb.StringCol(20)
    ifuslot = tb.StringCol(3)
    ifuid = tb.StringCol(3)
    specid = tb.StringCol(3)
    amp = tb.StringCol(2)
    normdata = tb.Float32Col((2, 100))


fileh = tb.open_file("cal_test.h5", mode="w", title="Calibration test file ")

group = fileh.create_group(fileh.root, 'Calibration', 'HETDEX Calibration Info')
groupAmpToAmp = fileh.create_group(group, 'AmpToAmp', 'Amp to amp Fiber Normalization')
groupThroughput = fileh.create_group(group, 'Throughput', 'Throughput Curves')

# populate the AmpToAmp group with normalization curves and metadata

amptoampfiles = glob.glob(op.join(path_amp2amp,'multi*.norm'))

for file in amptoampfiles:
    
    idx = file.find('multi_')
    ampid = file[idx:idx+20]
    data = ascii.read(file, names=['wave','norm'])
    fileh.create_table(groupAmpToAmp, str(ampid), data.as_array())

# populate Throughput group with throughput curves

throughputfiles = glob.glob(op.join(path_throughput, "*sedtp_f.dat"))

for file in throughputfiles:
    idx = file.find('/20')
    datevshot = file[idx+1:idx+13]
    tp_data = ascii.read(file)
    data = tp_data['col1','col2','col3']
    data.names = ['waves', 'response', 'response_err']
    fileh.create_table(groupThroughput, datevshot, data.as_array())


fileh.close()
