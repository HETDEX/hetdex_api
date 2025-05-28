"""
Config file for HETDEX data release paths
"""

import os
import os.path as op


class HDRconfig:

    LATEST_HDR_NAME = 'hdr3'

    def __init__(self, survey=LATEST_HDR_NAME):
        # Check stampede2 first
        if op.exists("/scratch/03946/hetdex"):
            self.host_dir = "/scratch/03946/hetdex"
        elif op.exists("/corral-repl/utexas/Hobby-Eberly-Teelsco"):
            self.host_dir = "/corral-repl/utexas/Hobby-Eberly-Teelsco"
        elif op.exists("/home/jovyan/Hobby-Eberly-Telesco"):
            self.host_dir = "/home/jovyan/Hobby-Eberly-Telesco"
        elif op.exists("/home/idies/workspace/HETDEX"):
            self.host_dir = "/home/idies/workspace/HETDEX"
        elif op.exists("/mnt/hetdex_release/"):
            self.host_dir = "/mnt/hetdex_release/"
        elif op.exists("/data/hetdex/u/dfarrow/hetdex_data"):
            self.host_dir = "/data/hetdex/u/dfarrow/hetdex_data"
        else:
            self.host_dir = os.getcwd()
            print('Could not find {} directory path'.format(survey))

        self.hdr_dir = {
            "hdr1": op.join(self.host_dir, "hdr1"),
            "hdr2": op.join(self.host_dir, "hdr2"),
            "hdr2.1": op.join(self.host_dir, "hdr2.1"),
            "hdr3": op.join(self.host_dir, "hdr3"),
        }

        self.survey_dir = op.join(self.hdr_dir[survey], "survey")
        self.software_dir = op.join(self.hdr_dir[survey], "software")
        self.red_dir = op.join(self.hdr_dir[survey], "reduction")
        self.data_dir = op.join(self.red_dir, "data")
        self.tp_dir = op.join(self.red_dir, "throughput")
        self.calib_dir = op.join(self.hdr_dir[survey], "calib")
        self.dustmaps = op.join(self.calib_dir, 'dustmaps')
        self.pixflat_dir = op.join(self.hdr_dir[survey], "calib/lib_pflat")
        self.raw_dir = op.join(self.hdr_dir[survey], "raw")
        self.flim_dir = op.join(self.red_dir, "flim")
        self.flim_sim_completeness = op.join(self.flim_dir, "snfiles")
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
        if survey == 'hdr2.1':
            self.elixerh5 = op.join(self.hdr_dir[survey], "detect", "elixer.h5")
        elif survey == 'hdr3':
            self.elixerh5 = op.join(self.hdr_dir[survey], "detect", "elixer_hdr3_emis_cat.h5")
        self.imaging_dir = op.join(self.host_dir, "imaging")
        self.contsourceh5 = op.join(
            self.hdr_dir[survey], "detect", "continuum_sources.h5"
        )
        self.fiberindexh5 = op.join(
            self.hdr_dir[survey], "survey", "fiber_index_" + survey + ".h5"
        )
        self.detectml = op.join(
            self.hdr_dir[survey], "detect", "detect_ml_" + survey + ".h5"
        )
        self.elix_dir = op.join(self.hdr_dir[survey], "detect", "image_db")

        if survey == "hdr1":
            if op.exists("/home/jovyan/software/hetdex_api"):
                self.bad_dir = "/home/jovyan/software/hetdex_api/known_issues/hdr1"
            else:
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

        if survey == "hdr2":
            if op.exists("/home/jovyan/software/hetdex_api"):
                self.bad_dir = "/home/jovyan/software/hetdex_api/known_issues/hdr2"
            else:
                self.bad_dir = (
                    "/work/05350/ecooper/wrangler/hetdex_api/known_issues/hdr2"
                )
            self.baddetect = op.join(self.bad_dir, "baddetects.list")
            self.badshot = op.join(self.bad_dir, "badshots.list")
            self.badamp = op.join(self.bad_dir, "badamps.list")
            self.badpix = op.join(self.bad_dir, "badpix.list")

        if survey in ["hdr2.1",'hdr3']:
            if op.exists("/home/jovyan/software/hetdex_api"):
                self.bad_dir = "/home/jovyan/software/hetdex_api/known_issues/{}".format(survey)
            elif op.exists("/afs/ipp-garching.mpg.de/home/d/dfarrow/hetdex/code/HETDEX_API/known_issues/{}".format(survey)):
                self.bad_dir = "/afs/ipp-garching.mpg.de/home/d/dfarrow/hetdex/code/HETDEX_API/known_issues/{}".format(survey)
            elif op.exists("/home/idies/workspace/HETDEX/hetdex_api/known_issues/{}".format(survey)):
                self.bad_dir = "/home/idies/workspace/HETDEX/hetdex_api/known_issues/{}".format(survey)
            else:
                self.bad_dir = (
                    "/work/05350/ecooper/stampede2/hetdex_api/known_issues/{}".format(survey)
                )

            self.baddetect = op.join(self.bad_dir, "baddetects.list")
            self.badshot = op.join(self.bad_dir, "badshots.list")
            self.badamp = op.join(self.hdr_dir[survey], "survey", "amp_flag.fits")
            self.badamp2 = op.join(self.bad_dir, "badamps.list")
            self.badamp_calfib = op.join(self.bad_dir, "badamps_missingcalfib.list")
            self.badamp_single = op.join(self.bad_dir, "badamps_single.list")
            self.badpix = op.join(self.bad_dir, "badpix.list")
            self.baddetectmask = op.join(
                self.hdr_dir[survey], "detect", "baddets_hdr2.1.0.p"
            )
            self.flim_avg = op.join(
                self.hdr_dir[survey], "survey", "flux_limits_all.txt"
            )
            self.meteor = op.join(self.bad_dir, "meteor.txt")
            self.flimmask = op.join(self.flim_dir, "masks")
            self.lowtpshots = op.join(self.bad_dir, "survey_shots_low_response.txt")
            self.rc3cat = op.join(self.bad_dir, "rc3.ugc.hetdex.both.v8.csv")
            self.agncat = op.join(self.detect_dir, "catalogs", "AGN_v5.2.dat")
            self.gaiacat = op.join(
                self.host_dir,
                "gaia_hetdex_value_added_catalog",
                "HDR2.1_Gaia_final_table.fits",
            )
            self.galaxylabels = op.join(self.bad_dir, 'galaxies.txt')
            self.starlabels = op.join(self.bad_dir, 'stars.txt')
            self.sdsscat = op.join(self.imaging_dir, 'sdss', 'specObj-dr16-trim.fits')
            self.extinction_fix = op.join(self.bad_dir, 'extinction')
            self.fibermaskh5 = op.join(self.survey_dir, 'fiber_mask_3.0.1.h5')

            if survey == 'hdr3':
                self.wdcor = op.join(self.bad_dir, 'wdcor.txt')
