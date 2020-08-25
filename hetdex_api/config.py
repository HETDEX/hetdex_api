"""
Config file for HETDEX data release paths
"""

import os.path as op
#import re
#import socket
import sys



class HDRconfig:

    LATEST_HDR_NAME = "hdr2.1"

    def __init__(self, survey=LATEST_HDR_NAME):
        # find out which cpu cluster is in use
        # hostname = socket.gethostname()
        # if re.search("wrangler", str(hostname)):
        #     self.host_dir = "/data/05350/ecooper"
        # elif re.search("stampede2", str(hostname)):
        #     self.host_dir = "/scratch/03946/hetdex"
        # else:
        #     sys.exit('Edit hetdex_api/config.py for your local dir')

        if op.exists("/data/05350/ecooper"):
            self.host_dir = "/data/05350/ecooper"
        elif op.exists("/scratch/03946/hetdex"):
            self.host_dir = "/scratch/03946/hetdex"
        else:
            sys.exit('Edit hetdex_api/config.py for your local dir')

        self.hdr_dir = {
            "hdr1": "/work/03946/hetdex/hdr1",
            "hdr2": "/data/05350/ecooper/hdr2",
            "hdr2.1": op.join(self.host_dir, "hdr2.1")}
        self.software_dir = op.join(self.hdr_dir[survey], "software")
        self.red_dir = op.join(self.hdr_dir[survey], "reduction")
        self.data_dir = op.join(self.red_dir, "data")
        self.tp_dir = op.join(self.red_dir, "throughput")
        self.calib_dir = op.join(self.hdr_dir[survey], "calib")
        self.pixflat_dir = op.join(self.hdr_dir[survey], "calib/lib_pflat")
        self.raw_dir = op.join(self.hdr_dir[survey], "raw")
        self.flim_dir = op.join(self.red_dir, "flim")
        self.elix_dir = op.join(self.hdr_dir[survey], "detect", "ergfiles")
        self.detect_dir = op.join(self.hdr_dir[survey], "detect")
        self.path_gpinfo = op.join(self.calib_dir, "DR1FWHM.txt")
        self.path_acc_flags = op.join(self.red_dir, "status_summary_hdr1.txt")
        self.path_radec = op.join(self.calib_dir, "radec.all")
        self.survey_list = op.join(self.red_dir, "hdr1.scilist")
        self.cal_list = op.join(self.red_dir, "hdr1.callist")
        self.surveyh5 = op.join(
            self.hdr_dir[survey], "survey", "survey_" + survey + ".h5"
        )
        self.detecth5 = op.join(
            self.hdr_dir[survey], "detect", "detect_" + survey + ".h5"
        )
        try:
            self.detectbroadh5 = op.join(
                self.hdr_dir[survey], "detect", "detect_broad_" + survey + ".h5"
            )
        except:
            pass
        self.elixerh5 = op.join(self.hdr_dir[survey], "detect", "elixer.h5")
        self.imaging_dir = op.join(self.hdr_dir[survey], "imaging")
        self.contsourceh5 = op.join(
            self.hdr_dir[survey], "detect", "continuum_sources.h5"
        )
        self.fiberindexh5 = op.join(
            self.hdr_dir[survey], "survey", "fiber_index_" + survey + ".h5"
        )
        self.detectml = op.join(self.hdr_dir[survey], "detect", "detect_ml_" + survey + ".h5" )
        self.elix_dir = op.join(self.hdr_dir[survey], "detect", "image_db")
               
        if survey == "hdr1":
            # here are files that are changing since HDR1 release
            self.bad_dir = "/work/05350/ecooper/hdr1/HETDEX_API/known_issues/hdr1"
            self.baddetect = op.join(self.bad_dir, "baddetects.list")
            self.badshot = op.join(self.bad_dir, "badshots.list")
            self.badamp = op.join(self.bad_dir, "badamps.list")
            self.badpix = op.join(self.bad_dir, "posthdr1badpix.list")
            self.gmags = op.join(self.bad_dir, "gmags.pickle")
            self.gmags_cont = op.join(self.bad_dir, "gmags_cont.pickle")
            self.plae_poii_hetdex_gmag = op.join(
                self.bad_dir, "plae_poii_hetdex_gmag.pickle"
            )

        if survey == 'hdr2':
            self.bad_dir = "/work/05350/ecooper/wrangler/hetdex_api/known_issues/hdr2"
            self.baddetect = op.join(self.bad_dir, "baddetects.list")
            self.badshot = op.join(self.bad_dir, "badshots.list")
            self.badamp = op.join(self.bad_dir, "badamps.list")
            self.badpix = op.join(self.bad_dir, "badpix.list")
            #self.elixerh5 = "/data/03261/polonius/hdr2/detect/elixer.h5"
            #self.elix_dir = "/data/03261/polonius/hdr2/detect/image_db"
            #self.imaging_dir = "/data/03261/polonius/hdr2/imaging"

        if survey == 'hdr2.1':
            self.bad_dir = "/work/05350/ecooper/wrangler/hetdex_api/known_issues/hdr2.1"
            self.baddetect = op.join(self.bad_dir, "baddetects.list")
            self.badshot = op.join(self.bad_dir, "badshots.list")
            self.badamp = op.join(self.hdr_dir[survey], "survey", "amp_flag.fits")
            self.badpix = op.join(self.bad_dir, "badpix.list")
            #self.elixerh5 = "/data/03261/polonius/hdr2.1.run/detect/elixer.h5"
            #self.elix_dir = "/data/03261/polonius/hdr2.1.run/detect/image_db"
            #self.imaging_dir = "/data/03261/polonius/hdr2/imaging"
            self.baddetectmask = op.join(self.hdr_dir[survey], "detect", "baddets_hdr2.1.0.p")
            self.flim_avg = op.join(self.hdr_dir[survey], "survey", "flux_limits_all.txt")
