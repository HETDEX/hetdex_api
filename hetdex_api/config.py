"""
Config file for HETDEX data release paths
"""

import os.path as op


class HDRconfig:
    def __init__(self, survey='hdr1'):
        self.hdr_dir = {'hdr1': '/work/03946/hetdex/hdr1',
                        'hdr2': '/scratch/05350/ecooper/hdr2'}
        self.software_dir = op.join(self.hdr_dir[survey], 'software')
        self.red_dir = op.join(self.hdr_dir[survey], 'reduction')
        self.data_dir = op.join(self.red_dir, 'data')
        self.tp_dir = op.join(self.red_dir, 'throughput')
        self.calib_dir = op.join(self.hdr_dir[survey], 'calib')
        self.raw_dir = op.join(self.hdr_dir[survey], 'raw')
        self.flim_dir = op.join(self.red_dir, 'flim')
        self.elix_dir = op.join(self.hdr_dir[survey], 'detect', 'ergfiles')
        self.path_gpinfo = op.join(self.calib_dir, 'DR1FWHM.txt')
        self.path_acc_flags = op.join(self.red_dir, 'status_summary_hdr1.txt')
        self.path_radec = op.join(self.calib_dir, 'radec.all')
        self.survey_list = op.join(self.red_dir, 'hdr1.scilist')
        self.cal_list = op.join(self.red_dir, 'hdr1.callist')
        self.surveyh5 = op.join(self.hdr_dir[survey],'survey','survey_' + survey + '.h5')
        self.detecth5 = op.join(self.hdr_dir[survey],'detect','detect_' + survey + '.h5')
        self.elixerh5 = op.join(self.hdr_dir[survey],'detect','elixer.h5')
        self.contsourceh5 = op.join(self.hdr_dir[survey],'detect','continuum_sources.h5')

        if (survey=='hdr1'):
            # here are files that are changing since HDR1 release
            self.bad_dir = '/work/05350/ecooper/hdr1/HETDEX_API/known_issues/hdr1'
            self.baddetect = op.join(self.bad_dir, 'baddetects.list')
            self.badshot = op.join(self.bad_dir, 'badshots.list')
            self.badamp = op.join(self.bad_dir, 'badamps.list')
            self.badpix = op.join(self.bad_dir, 'posthdr1badpix.list')
            self.gmags = op.join(self.bad_dir, 'gmags.pickle')
            self.gmags_cont = op.join(self.bad_dir, 'gmags_cont.pickle')
            self.plae_poii_hetdex_gmag = op.join(self.bad_dir, 'plae_poii_hetdex_gmag.pickle')
