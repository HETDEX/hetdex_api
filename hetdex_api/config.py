"""
Config file for HETDEX data release paths
"""

import os.path as op

hdr_dir = '/work/03946/hetdex/hdr1'
software_dir = op.join(hdr_dir, 'software')
red_dir = op.join(hdr_dir, 'reduction')
data_dir = op.join(red_dir, 'data')
tp_dir = op.join(red_dir, 'throughput')
calib_dir = op.join(hdr_dir, 'calib')
raw_dir = op.join(hdr_dir, 'raw')
flim_dir = op.join(red_dir, 'flim')
elix_dir = op.join(hdr_dir, 'detect', 'ergfiles')

path_gpinfo = op.join(calib_dir,'DR1FWHM.txt')
path_acc_flags = op.join(red_dir, 'status_summary_hdr1.txt')
path_radec = op.join(calib_dir, 'radec.all')

survey_list = op.join(red_dir, 'hdr1.scilist')
cal_list = op.join(red_dir, 'hdr1.callist')

surveyh5 = op.join(hdr_dir,'survey','survey_hdr1.h5')
detecth5 = op.join(hdr_dir,'detect','detect_hdr1.h5')
elixerh5 = op.join(hdr_dir,'detect','elixer.h5')
contsourceh5 = op.join(hdr_dir,'detect','continuum_sources.h5')

# here are files that are changing since HDR1 release
bad_dir = '/work/05350/ecooper/stampede2/HETDEX_API/known_issues/hdr1'
baddetect = op.join(bad_dir, 'baddetects.list')
badshot = op.join(bad_dir, 'badshots.list')
badamp = op.join(bad_dir, 'badamps.list')
badpix = op.join(bad_dir, 'posthdr1badpix.list')
gmags = op.join(bad_dir,'gmags.pickle')
gmags_cont = op.join(bad_dir, 'gmags_cont.pickle')
plae_poii_hetdex_gmag = op.join(bad_dir,'plae_poii_hetdex_gmag.pickle')
