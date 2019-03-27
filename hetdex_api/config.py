#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 12:21:35 2019

@author: gregz
"""
import os.path as op

hdr_dir = '/work/03946/hetdex/hdr1'
software_dir = op.join(hdr_dir, 'software')
red_dir = op.join(hdr_dir, 'reduction')
data_dir = op.join(red_dir, 'data')
tp_dir = op.join(red_dir, 'throughput')
calib_dir = op.join(hdr_dir, 'calib')
raw_dir = op.join(hdr_dir, 'raw')

path_gpinfo = op.join(calib_dir,'MasterFWHM.txt')
path_acc_flags = op.join(red_dir, 'status_summary_hdr1.txt')
path_radec = op.join(calib_dir, 'radec.all')

survey_list = op.join(red_dir, 'hdr1.scilist')
cal_list = op.join(red_dir, 'hdr1.callist')

surveyh5 = op.join(hdr_dir,'survey','survey_hdr1.h5')
detecth5 = op.join(hdr_dir,'detects','detect_cosmos.h5')
