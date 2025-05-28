"""
Config file for HETDEX data release paths
"""

import os
import os.path as op




class HDRconfig:

    LATEST_HDR_NAME = 'hdr5'
    LAST_GOOD_HDR_NAME = 'hdr5' #add for use when LATEST_HDR_NAME is overridden

# Will just use current for now but I'll leave this here in case we go back 2024-12-06 - EMC
#    LATEST_MASK_DICT = {'hdr3': '1.0',
#                        'hdr4': '1.0',
#                        'hdr5': '1.0',
#                        'pdr1': '1.0'}


    def find_host_directory(self):
        # in descending priority
        possible_directories = [
            "/scratch/projects/hetdex",
            #remove from possible options EMC 2024-09-03
#            "/scratch/03946/hetdex",
            "/corral-repl/utexas/Hobby-Eberly-Telesco",
            "/home/jovyan/Hobby-Eberly-Telesco",
            "/home/idies/workspace/HETDEX",
            "/mnt/hetdex_release/"
        ]


        def check_directory(directory):
            return op.exists(directory) and os.access(directory, os.R_OK)
    
        for directory in possible_directories:
            if check_directory(directory):
                return directory
    
        print(f'Could not find a suitable directory path with read permission.')
        return os.getcwd()

    
    def __init__(self, survey=LATEST_HDR_NAME):
        
        # a number of downstream operations use LATEST_HDR_NAME explicitly, regardless of
        # what was pssed in for survey ...
        # overwriting here should handle that _BUT_ might introduce other problems if we do actually want to
        # use the pre-defined LATEST_HDR_NAME even if a different survey was requested.
        if survey is not None:
            self.LATEST_HDR_NAME=survey

        self.host_dir = self.find_host_directory()

        self.hdr_dir = {
            "hdr1": op.join(self.host_dir, "hdr1"),
            "hdr2": op.join(self.host_dir, "hdr2"),
            "hdr2.1": op.join(self.host_dir, "hdr2.1"),
            "hdr3": op.join(self.host_dir, "hdr3"),
            "hdr4": op.join(self.host_dir, "hdr4"),
            "hdr5": op.join(self.host_dir, "hdr5"),
            "pdr1": op.join(self.host_dir, "pdr1"),
        }

        self.survey_dir = op.join(self.hdr_dir[survey], "survey")
        self.software_dir = op.join(self.hdr_dir[survey], "software")
        self.red_dir = op.join(self.hdr_dir[survey], "reduction")
        self.data_dir = op.join(self.red_dir, "data")
        self.mask_dir = op.join(self.host_dir, 'hdr5/reduction/mask')
        self.tp_dir = op.join(self.red_dir, "throughput")
        #self.calib_dir = op.join(self.hdr_dir[survey], "calib") #DD 2023-08-05
        self.calib_dir = op.join(self.host_dir, "lib_calib")
        #self.dustmaps = op.join(self.calib_dir, 'dustmaps')
        self.dustmaps = op.join(self.host_dir, 'hdr3/calib/dustmaps') #DD 2023-08-05
        #self.pixflat_dir = op.join(self.hdr_dir[survey], "calib/lib_pflat")
        self.pixflat_dir = op.join(self.calib_dir, "lib_pflat") #DD 2023-08-05
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
        self.detectindexh5 = op.join(
            self.hdr_dir[survey], "survey", "detect_index_" + survey + ".h5"
        )
        self.elixerh5 = op.join(self.hdr_dir[survey], "detect", "elixer.h5")
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

        if survey in ["hdr2.1",'hdr3', 'hdr4', 'hdr5', 'pdr1']:
            if op.exists("/home/jovyan/software/hetdex_api"):
                self.bad_dir = "/home/jovyan/software/hetdex_api/known_issues/{}".format(survey)
            elif op.exists('/home1/05350/ecooper/hetdex_api/known_issues/{}'.format(survey)):
                self.bad_dir = '/home1/05350/ecooper/hetdex_api/known_issues/{}'.format(survey)
            elif op.exists("/data/hetdex/u/dfarrow/hetdex_data/{}/{}_issues".format(survey, survey)):
                self.bad_dir = "/data/hetdex/u/dfarrow/hetdex_data/{}/{}_issues".format(survey, survey)
            elif op.exists("/home/idies/workspace/HETDEX/hetdex_api/known_issues/{}".format(survey)):
                self.bad_dir = "/home/idies/workspace/HETDEX/hetdex_api/known_issues/{}".format(survey)
            else:
                self.bad_dir = (
                    "/work/05350/ecooper/stampede2/hetdex_api/known_issues/{}".format(survey)
                )

            self.badfib = op.join(self.bad_dir, "badfib.tab")
            self.badfib_transient = op.join(self.bad_dir, 'badfib_bydate.tab')
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
                self.hdr_dir[survey], "survey", "flim4600_averages.txt"
            )

            self.meteor = op.join(self.bad_dir, "meteor.txt")
            self.flimmask = op.join(self.flim_dir, "masks")
            self.lowtpshots = op.join(self.bad_dir, "survey_shots_low_response.txt")
            self.rc3cat = op.join(self.bad_dir, "rc3.ugc.hetdex.both.v8.csv")
            #self.agncat = op.join(self.hdr_dir['hdr4'], "catalogs", "hdr234_agn_v4.5.1.cat")
            #this now needs to be softlinks if there is not a unique file for the release
            #updated 2025-05-06 by EMC
            self.agncat = op.join(self.hdr_dir[survey], "catalogs", "hdr345_agn_v5.1.cat")
            #self.agncat = op.join(self.hdr_dir[survey], "catalogs", "hdr234_agn_v4.5.4.cat")
            self.gaiacat = op.join(
                self.host_dir,
                "gaia_hetdex_value_added_catalog",
                "HDR2.1_Gaia_final_table.fits",
            )
            self.galaxylabels = op.join(self.bad_dir, 'galaxies.txt')
            self.starlabels = op.join(self.bad_dir, 'stars.txt')
            self.sdsscat = op.join(self.imaging_dir, 'sdss', 'specObj-dr16-trim.fits')
            self.extinction_fix = op.join(self.bad_dir, 'extinction')
            self.fibermaskh5 = op.join(self.survey_dir, 'fiber_mask.h5')

            #if survey == 'hdr3':
            try:
                if float(survey[3:]) >= 3.0:
                    self.wdcor = op.join(self.bad_dir, 'wdcor.txt')
            except:
                self.wdcor = None

            #added 20241004: average fit to p-cygni like feature from bad sky calibration at 5200 A
            self.cal_issue = op.join(self.bad_dir, 'spec5200fit.npy')
            self.cal5200 = op.join(self.bad_dir, 'flag5200.txt')
            self.cal5460 = op.join(self.bad_dir, 'flag5460.txt')
            self.satellite = op.join(self.bad_dir, 'satellite_tracks.txt')
